#include "random.h"
#include <iostream>
#include <cmath>
#include "Funzioni.h"
#include <cfloat>
#include <stdio.h>   //per il confronto delle stringhe per la funzione Run()
#include <string.h>
#include <fstream>
#include "Metropolis.h"

using namespace std;

Metropolis::Metropolis(vector<double> posizioni0, int N, double passo, int N_dimension, FunzioneBase * f){
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
            	_rnd.SetRandom(seed,p1,p2);
        }
        }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


	_f = f;
	_passo = passo;
	_N = N;
	_posizioni.push_back(posizioni0);  //attacco la posizione di partenza alla matrice
	_Ndim = N_dimension;
	_accepted = 0;
}

Metropolis::~Metropolis(){
	_rnd.SaveSeed();
}

void Metropolis::Run(const char * metodo){   //metodo può essere uguale a "gauss" o "uniform" a seconda che io voglia un passo gaussiano o uniforme
	double incremento = 0;	  // se non lo inizializzo mi dà un avviso il compilatore (tanto non cambia nulla ai fini dell'algoritmo)
	for(int i = 0; i< _N; i++){
        vector<double> possible;
		for(int j = 0; j< _Ndim; j++){

			if(strcmp(metodo,"uniform") == 0){    //metodo per confrontare due stringhe se sono uguali
				incremento = _rnd.Rannyu(-_passo,_passo);  //incremento uniforme
			}
			else if(strcmp(metodo,"gauss") == 0){
				incremento = _rnd.Gauss(0,_passo);   //incremento gaussiano
			}

			possible.push_back(incremento+_posizioni[i][j]);   //possibile posizione futura
		}
		double prob_ratio = _f->Eval(possible)/_f->Eval(_posizioni[i]);
		if(prob_ratio>=1){
			_posizioni.push_back(possible);  //mi muovo
			_accepted++;

		}
		else{
			double estrazione =  _rnd.Rannyu(0.,1);
			if(estrazione<=prob_ratio){
				_posizioni.push_back(possible);	//mi muovo
				_accepted++;
			}
			else{
				_posizioni.push_back(_posizioni[i]);  //rimango fermo
			}

		}

	}
}

void Metropolis::PrintAccepted(){
	cout<<"Acceptance rate: "<<((double)_accepted/(double)_N)*100.<<"%"<<endl;
	_accepted = 0;
}

void Metropolis::Equilibrate(int nstep, const char * metodo){ //equilibra l'algoritmo in nstep
	int passi = _N;  //così non perdo il numero di passi con cui voglio eseguire la simulazione (numero che ho impostato nel costruttore)
	_N = nstep;
	Run(metodo);
	Restart();
	_N = passi;
}



void Metropolis::Restart(){
	_posizioni.erase(_posizioni.begin(),_posizioni.end()-1);  //elimina tutti i passi fatti tranne la posizione finale ottenuta
	_accepted = 0;
}

void Metropolis::Print(const char * namefile) {
	ofstream Risultati;
	Risultati.open(namefile);

	for(unsigned int i=0;i<_posizioni.size();i++){
		for(int j = 0; j<_Ndim;j++){
		if (Risultati.is_open()){
   			Risultati <<_posizioni[i][j]<<" ";
   		}
	}
	Risultati<<endl;


	}

Risultati.close();

}


vector<double> Metropolis::Raggi(){
	vector<double> r;
	double appo;
	for(unsigned int i=0;i<_posizioni.size();i++){
		appo =0;
		for(int j = 0; j<_Ndim;j++){
			appo+= pow(_posizioni[i][j],2);
		}
		r.push_back(sqrt(appo));
	}
	return r;
}
