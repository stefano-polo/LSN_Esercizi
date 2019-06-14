#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "Salesman.h"


using namespace std;

double signum (double x){
    if (x>=0){return 1.;}else{return -1.;}
}

Salesman::Salesman(int NumberOfCity, int NumberOfParents, int n_generation,int seed): _FitnessDistances(NumberOfParents){

  Start(&_rnd);
   _n_generation= n_generation;
  _NumberOfCity = NumberOfCity;
  _NumberOfParents = NumberOfParents;

    _best_index=0;

  srand(seed);
  vector<int> cities(_NumberOfCity);
  for (int i=0; i<_NumberOfCity; i++){
    cities.at(i) = i;
  }

  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());  //faccio una permutazione differenze per ogni stringa dei numeri generati
    _StringsOfCities.push_back(cities);}  //genero i genitori

  _NewGeneration = _StringsOfCities;
  _accepted = 0;
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

void Salesman::Print_Positions(void){
  cout<<"x  y"<<endl;
  for(int j=0; j<_NumberOfCity;j++){
    for (int i = 1; i < 3; i++) {
      cout<<_Positions[j][i]<<"  ";
    }
    cout<<endl;
  }


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


void Salesman::Fitness(void){  //calcolo delle distanze di ciascuna stringa
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

    for(unsigned int i=0; i<_remember_best_distance.size();i++){
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

void Salesman::Metropolis(double beta, int number_pair_mutation){
   Best_string();
   _remember_best_distance.push_back(_FitnessDistances.at(_best_index));

  vector<double> old_fitness = Fitness_Annealing(_StringsOfCities);

  Mutation_Annealing(number_pair_mutation);
  vector<double> new_fitness = Fitness_Annealing(_NewGeneration);

  for(unsigned int i = 0; i < _FitnessDistances.size(); i++){
   double boltzmann = exp(-beta*(new_fitness[i]-old_fitness[i]));
   if(boltzmann>=1){
	_StringsOfCities[i] = _NewGeneration[i];
	_accepted++;
   }
   else {
	    double number_trial = _rnd.Rannyu();
	    if(number_trial<=boltzmann){
		_StringsOfCities[i] = _NewGeneration[i];
		_accepted++;
	}

  }
  }
  _NewGeneration = _StringsOfCities;

  return;
}

void Salesman::Run_Metropolis(double beta, int number_pair_mutation) {

	for(int i = 0; i < _n_generation; i++) {
  	 Metropolis(beta, number_pair_mutation);
     _Best_People.push_back(_StringsOfCities[_best_index]);
  }
	Best_string();

}

void Salesman::Equilibration(int numero_equilibrazioni,double beta, int number_pair_mutation) {

	for(int i = 0; i < numero_equilibrazioni; i++) {
		Metropolis( beta,  number_pair_mutation);

	}
	_remember_best_distance.erase(_remember_best_distance.begin(),_remember_best_distance.end());
	_accepted = 0;
}



void Salesman::Best_string(void){
    Fitness();
    _best_index=min_element(_FitnessDistances.begin(),_FitnessDistances.end()) - _FitnessDistances.begin();

}

vector<double> Salesman::Fitness_Annealing(vector<vector<int>> v) {
  vector<double> Distance(_NumberOfParents);
  for(int j=0; j<_NumberOfParents; j++){
    double L=0;
    for (int i=0; i<_NumberOfCity-1;i++){
      int index = v[j][i];
      int index_plus = v[j][i+1];
      L+=pow(_Positions[index][1]-_Positions[index_plus][1],2) + pow(_Positions[index][2]-_Positions[index_plus][2],2);
    }
    L+=pow(_Positions[v[j][0]][1]-_Positions[v[j][_NumberOfCity-1]][1],2)+pow(_Positions[v[j][0]][2]-_Positions[v[j][_NumberOfCity-1]][2],2);  //in questo punto faccio in modo che la fine della stringa si ricolleghi con l'inizio
    Distance[j] = L;

  }
  return Distance;
}

void Salesman::Mutation_Annealing(int number_mutation){
  for(int j = 0; j<number_mutation; j++){
  for (unsigned int i=0; i <_NewGeneration.size();i++){
        unsigned int index_start=int(_rnd.Rannyu(1,_NumberOfCity));
        unsigned int index_end = int(_rnd.Rannyu(1,_NumberOfCity));
        int tram=_NewGeneration[i][index_start];
        _NewGeneration[i][index_start]=_NewGeneration[i][index_end];
        _NewGeneration[i][index_end]=tram;

    }
  }
    return;
}

void Salesman::PrintAccepted(){
	cout<<"Acceptance rate: "<<((((double)_accepted/(double)_n_generation))/(double)_NumberOfParents)*100.<<"%"<<endl;
        _accepted = 0;
}

void Salesman::Random_Search(void){
	 for(int i=0; i<_NumberOfParents; i++){
   	 random_shuffle(_StringsOfCities[i].begin()+1, _StringsOfCities[i].end());
   	 }

	return;
}

double Salesman::Minimum_distance(void) {
	int index_min = min_element(_remember_best_distance.begin(),_remember_best_distance.end()) - _remember_best_distance.begin();
	double minimum = _remember_best_distance[index_min];

	return minimum;
}

int* Salesman::Best_Cammino(void) {
  int * best = new int[_NumberOfCity];
  int index = min_element(_remember_best_distance.begin(),_remember_best_distance.end()) - _remember_best_distance.begin();
  for(int i = 0; i<_NumberOfCity; i++){
    best[i] = _Best_People[index][i];
  }
  return best;

}

vector<vector<double>> Salesman::Get_Positions(void){
  return _Positions;
}
