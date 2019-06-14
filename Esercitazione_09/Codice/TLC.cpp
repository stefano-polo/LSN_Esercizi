#include "TLC.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void Start(Random *rnd){

	int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()){
     Primes >> p1 >> p2 ;
  } else cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()){
     while ( !input.eof() ){
        input >> property;
        if( property == "RANDOMSEED" ){
           input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
           rnd->SetRandom(seed,p1,p2);
        }
     }
     input.close();
  } else cerr << "PROBLEM: Unable to open seed.in" << endl;

}

double Sum(vector<double> v, int start, int end){

	double sum = 0;
	for(int i=start; i<end; i++) sum += v.at(i);
	return sum;
}

double Mean(vector<double> v) {
	double sum = 0;
	for(unsigned int i=0; i<v.size(); i++) sum += v.at(i);
	return sum/v.size();
}

double Error(vector<double> v) {
	double sum = 0;
	for(unsigned int i=0; i<v.size(); i++) sum += pow(Mean(v)-v.at(i), 2.);
	return sum/(v.size()-1);
}

void TLC( int N, vector<double> v, const char*namefile){         // Mi permette di trovare le medie e le deviazioni standard della media attrsaverso il metodo dei blocchi
	unsigned int M = v.size();
  int L = int(M/N);   //Numero di esperimenti per blocco
  vector<double> vmean(N,0);   //Vettore che mi contiene le medie dei singoli blocchi
  vector<double> vmean2(N,0);   //Vettore che contiene i quadrati delle medie dei singoli blocchi
  vector<double> mean_cum(N,0);  //Vettore che mi contiene la media cumulativa
  vector<double> mean2_cum(N,0);   //Vettore che contiene i quadrati delle medie cumulativi
  vector<double> err_cum(N,0);      //vettore che mi contiene  gli errori

  for(int i=0; i<N; i++){
    double temp = Sum(v,i*L,(i+1)*L)/L;
    vmean.at(i) = temp;
    vmean2.at(i) = pow(temp , 2);
	  }

  for(int i=0; i<N; i++){
    mean_cum.at(i) = Sum(vmean, 0, i+1)/(i+1);
    mean2_cum.at(i) = Sum(vmean2, 0, i+1)/(i+1);
    err_cum.at(i) = sqrt( ( mean2_cum.at(i) - pow(mean_cum.at(i),2) ) / (i+1)  );
    }

  ofstream results;
  results.open(namefile);

  for(int i=0; i<N; i++){

   if(results.is_open()) results << i*L << " " << mean_cum.at(i) << " " << err_cum.at(i) << endl;
   else cerr << "PROBLEM: Unable to open result.txt" << endl;
  }

  results.close();

}

void AnalysisMatrix(int N, vector<vector <double>> v, double scale_factor, double index_traslation, const char* namefile)
{

	//salvataggio dati su file
	ofstream results;
	results.open(namefile, ios::app);

	//parametri per il datta blocking
	int nrighe=v.size(); // numero di vettori su cui si fa il data blocking
	int M= v[0].size(); // lunghezza di ogni vettore su cui si effettua data blocking
	int L = int(M/N); // lnghezza di ogni blocco

	for(int j=0; j<nrighe; j++){ //data blocking su ogni riga della matriche

		vector<double> ai(N);
		vector<double> ai_2(N);
	 	double a; // variabile di appoggio

	  	for(int i=0; i<N;i++){
			a= Sum(v[j],i*L,(i+1)*L)/(L);
			ai.at(i) = a;
			ai_2.at(i) = pow(a,2);
		}

	//calcolo la media progressiva del vettore mean e le deviazione standard della media progressiva
	 	vector<double> mean_prog(N);
		vector<double> mean_2_prog(N);
		vector<double> mean_std(N);

		for(int i=0;i<N;i++){
			mean_prog.at(i) = (Sum(ai,0,i+1)/(i+1));
			mean_2_prog.at(i) = Sum(ai_2,0,i+1)/(i+1) ;
			mean_std.at(i) = sqrt( ( mean_2_prog.at(i) - pow(mean_prog.at(i),2) ) / (i+1)  );
		}

	// Stampo su file i risultati nel seguente formato:indice, media progressiva, deviazione standard della media progressiva con il massimo numero di blocchi

		if (results.is_open()){
	   		results << j*scale_factor - index_traslation << " "<< mean_prog.at(N-1) << " " <<mean_std.at(N-1)<< endl;
	  	 }
		else cerr << "PROBLEM: Unable to open the file.out" << endl;


	}

   results.close();
}
