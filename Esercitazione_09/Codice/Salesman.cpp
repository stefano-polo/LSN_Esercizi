#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "TLC.h"
#include "Salesman.h"


using namespace std;

double signum (double x){
    if (x>=0){return 1.;}else{return -1.;}
}

Salesman::Salesman(int NumberOfCity, int NumberOfParents, int n_generation): _FitnessDistances(NumberOfParents){

  Start(&_rnd);
   _n_generation= n_generation;
  _NumberOfCity = NumberOfCity;
  _NumberOfParents = NumberOfParents;

    _best_index=0;

  srand(10);
  vector<int> cities(_NumberOfCity);
  for (int i=0; i<_NumberOfCity; i++){
    cities.at(i) = i;
  }

  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());  //faccio una permutazione differenze per ogni stringa dei numeri generati
    _StringsOfCities.push_back(cities);}  //genero i genitori

}


Salesman::~Salesman(){
  _rnd.SaveSeed();
}


void Salesman::Print(void){

  for(int j=0; j<_NumberOfParents; j++){
		for(int i=0; i<_NumberOfCity;i++){
			cout << _StringsOfCities[j][i] << " ";
		}
    cout << endl;
	}
  return;
}

void Salesman::Square(void){ //genero randomicamente delle posizioni in un quadrato di lato 1
  double L=1.; //lato del quadrato
  vector<double> tram(3);
  for (int i=0; i<_NumberOfCity;i++){
    double x=_rnd.Rannyu(0,L);
    double y=_rnd.Rannyu(0,L);
    tram[0]=i; tram[1]=x; tram[2]=y;
    _Positions.push_back(tram);
  }
  return;
}



void Salesman::Circle(void){  //genero randomicamente delle posizioni su una circonferenza di raggio 1
  double R=1.; //raggio circonferenza
  for (int i=0; i<_NumberOfCity;i++){
    double x=_rnd.Rannyu(-R,R);
    double y= sqrt(R*R -x*x)*signum(_rnd.Rannyu(-1,1));
    vector<double> tram(3);
    tram[0]=i; tram[1]=x; tram[2]=y;
    _Positions.push_back(tram);
  }

  return;
}


void Salesman::Fitness(void){  //calcolo delle distanze totali di ciascuna stringa
  for(int j=0; j<_NumberOfParents; j++){
    double L=0;
    for (int i=0; i<_NumberOfCity-1;i++){
      int index = _StringsOfCities[j][i];
      int index_plus = _StringsOfCities[j][i+1];
      L+=pow(_Positions[index][1]-_Positions[index_plus][1],2) + pow(_Positions[index][2]-_Positions[index_plus][2],2);
    }
    L+=pow(_Positions[_StringsOfCities[j][0]][1]-_Positions[_StringsOfCities[j][_NumberOfCity-1]][1],2)+pow(_Positions[_StringsOfCities[j][0]][2]-_Positions[_StringsOfCities[j][_NumberOfCity-1]][2],2);  //in questo punto faccio in modo che la fine della stringa si ricolleghi con l'inizio
    _FitnessDistances[j] = L;

  }
  return;
}

void Salesman::Crossover(int parent1, int parent2){ //metodo di crossover che fa generare un solo figlio a partire da due stringhe genitrici scelte dalla matrice con indici di riga (parent1 e parent2)

  unsigned int index_start = int(_rnd.Rannyu(1,_NumberOfCity));
  unsigned int index_end = int(_rnd.Rannyu(index_start,_NumberOfCity));


  vector<int> son(_NumberOfCity, 0);
  for(unsigned int i=index_start; i<index_end; i++){
    son.at(i) = _StringsOfCities[parent1][i];

  }
//Prima parte del figlio
  int actual_position=0;
  for(unsigned int i=0; i<index_start; i++){

    for(int j=actual_position; j<_NumberOfCity; j++){//cambiamento

      int counter = 0;
      for(unsigned int k=index_start; k<index_end; k++){

          if(_StringsOfCities[parent2][j] == son.at(k)){

              counter++;}

      }

      if(counter == 0){

          son.at(i) = _StringsOfCities[parent2][j];
          actual_position = j+1;//crucial
          break;
      }
    }
  }
//Seconda parte del figlio
  for(int i=index_end; i<_NumberOfCity; i++){

    for(int j=actual_position; j<_NumberOfCity; j++){
        int counter = 0;
      for(unsigned int k=index_start; k<index_end; k++){
        if(_StringsOfCities[parent2][j] == son.at(k)) counter++;
      }
        if(counter == 0){ son.at(i) = _StringsOfCities[parent2][j];actual_position = j+1; break;}//crucial
    }
  }


  _NewGeneration.push_back(son);
  return;

}


void Salesman::Mutation(double p ){  //mutazione nella generazione figlia con probabilità p
    double random=_rnd.Rannyu(0.,1.);

    if(random<p){

    for (unsigned int i=0; i <_NewGeneration.size();i++){
        unsigned int index_start=int(_rnd.Rannyu(1,_NumberOfCity));
        unsigned int index_end = int(_rnd.Rannyu(1,_NumberOfCity));
        int tram=_NewGeneration[i][index_start];
        _NewGeneration[i][index_start]=_NewGeneration[i][index_end];
        _NewGeneration[i][index_end]=tram;

    }}
    return;
}

void Salesman::Breeders_selection(void){
    //metodo che mi seleziona la classe fertile (ossia che sono in grado di fare crossover) all'interno della generazione genitrice
    vector <int> breeders; //indici delle stringhe genitrici fertili
    vector <double> ordered_distances=_FitnessDistances; //vettore delle distanze in ordine crescente
    sort(ordered_distances.begin(),ordered_distances.end());


    for(int i=0;i<int(ordered_distances.size()/4);i++){ //vado a prendere 1/4 delle migliori distanze
        double key=ordered_distances.at(i);
        vector<double>::iterator walker = find(_FitnessDistances.begin(), _FitnessDistances.end(), key);
        breeders.push_back(distance(_FitnessDistances.begin(), walker));
    }

    for(unsigned int i=ordered_distances.size()/5;i<ordered_distances.size()/2;i++){ //vado a pescare randomicamente tra le stringhe che non ho preso fino alla metà delle distanze più brevi
        double random=int(_rnd.Rannyu(ordered_distances.size()/5,ordered_distances.size()/2));
        double key=ordered_distances.at(random);
        vector<double>::iterator walker = find(_FitnessDistances.begin(), _FitnessDistances.end(), key);
        breeders.push_back(distance(_FitnessDistances.begin(), walker));
    }

    _breeders=breeders;
    return;
}

int Salesman::Selection(void){
    double p=int(_rnd.Rannyu(0,_breeders.size()));  //prendo random tra gli individui fertili
    return _breeders.at(p);
}


void Salesman::Run(double probability_mutation){

    for(int i=0; i <_n_generation;i++){

        Best_string();

        Update_Distance50(i);
        _remember_best_distance.push_back(_FitnessDistances.at(_best_index));
        for(int j=0; j<_NumberOfParents; j++){
            Breeders_selection();
            int dad=Selection(); //seleziono gli indici dei genitori che si accoppiano
            int mom=Selection();
            Crossover( dad, mom);  //creo la nuova generazione figlia

        }
        Mutation(probability_mutation);
        _StringsOfCities=_NewGeneration;  //sostituisco la vecchia generazione con la nuova
        _NewGeneration.erase(_NewGeneration.begin(),_NewGeneration.end()); //resetto la matrice della nuova generazione

    }

    Best_string();  //aggiorno la stringa con la distanza minore percorsa per il grafico del path
    return;
}

void Salesman::Update_Distance50(int generazione) {
  vector <double> ordered_distances=_FitnessDistances; //vettore delle distanze in ordine crescente
  sort(ordered_distances.begin(),ordered_distances.end());
  vector<double> distanze_50(int(_NumberOfParents/2));
  for(int i = 0; i< int(_NumberOfParents/2); i++){
    distanze_50[i] = ordered_distances[i];
  }
  vector<double> appoggio(3);
  appoggio[0] = generazione;
  appoggio[1] = Mean(distanze_50); //media delle distanze al quadrato delle prime 50 migliori configurazioni
  appoggio[2] = sqrt(Error(distanze_50));
  _best_distances50.push_back(appoggio);
  return;
}

void Salesman::Best_string(void){
    Fitness();
    _best_index=min_element(_FitnessDistances.begin(),_FitnessDistances.end()) - _FitnessDistances.begin();

}

void Salesman::Print_Best_distance50(const char * file_name) {
  ofstream results;
  results.open(file_name);

  for(unsigned int i=0; i<_best_distances50.size();i++){
    results<<_best_distances50[i][0]<<" "<<_best_distances50[i][1]<<" "<<_best_distances50[i][2]<<endl;

  }



  results.close();
  return;

}

void Salesman::Best_path(const char * file_name){

    ofstream results;
    results.open(file_name);

    Best_string();


    //stampo  x, y delle città in ordine di passaggio dalla prima all'ultima
    for (int i=0; i <_NumberOfCity;i++){
        results<<_Positions[_StringsOfCities[_best_index][i]][1]<<" "<<_Positions[_StringsOfCities[_best_index][i]][2]<<endl;
    }
    //stampo come posizione finale la hometown
    results<<_Positions[_StringsOfCities[_best_index][0]][1]<<" "<<_Positions[_StringsOfCities[_best_index][0]][2]<<endl;


    results.close();



    return;
}


void Salesman::Best_distance(const char * file_name){
    ofstream results;
    results.open(file_name);

    for(int i=0; i<_n_generation;i++){
        results<<_remember_best_distance[i]<<endl;
    }



    results.close();
    return;

}

void Salesman::Restart(void){
  _StringsOfCities.erase(_StringsOfCities.begin(),_StringsOfCities.end());
  _NewGeneration.erase(_NewGeneration.begin(),_NewGeneration.end());
  _Positions.erase(_Positions.begin(),_Positions.end());
  _best_index=0;
  _remember_best_distance.erase(_remember_best_distance.begin(),_remember_best_distance.end());
  srand(10);
  vector<int> cities(_NumberOfCity);
  for (int i=0; i<_NumberOfCity; i++){
    cities.at(i) = i;
  }

  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());
    _StringsOfCities.push_back(cities);}


}
