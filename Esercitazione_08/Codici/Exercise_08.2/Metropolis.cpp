#include "random.h"
#include <iostream>
#include <cmath>
#include "Funzioni.h"
#include <cfloat>
#include <stdio.h>   //per il confronto delle stringhe per la funzione Run()
#include <string.h>
#include <fstream>
#include "Metropolis.h"
#include "TLC.h"

using namespace std;

Metropolis::Metropolis(vector<double> posizioni0, int N, double passo, const char * metodo, int N_dimension, FunzioneBase * f, double a, double b): _walker(_Nbins) {
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
	_metodo = metodo;

	//parametri per l'istogramma

	_a = a;
	_b = b;

	//inizializzo a zero tutti i bin dell'istogramma
	for(int j = 0; j<_Nbins; j++){
			_walker[j] = 0;
		}

}

Metropolis::~Metropolis(){
	_rnd.SaveSeed();

}

void Metropolis::Run(){   //metodo può essere uguale a "gauss" o "uniform" a seconda che io voglia un passo gaussiano o uniforme

	double incremento = 0;	  // se non lo inizializzo mi dà un avviso il compilatore (tanto non cambia nulla ai fini dell'algoritmo)
	for(int i = 0; i< _N; i++){
        vector<double> possible;
		for(int j = 0; j< _Ndim; j++){

			if(strcmp(_metodo,"uniform") == 0){    //metodo per confrontare due stringhe se sono uguali
				incremento = _rnd.Rannyu(-_passo,_passo);  //incremento uniforme
			}
			else if(strcmp(_metodo,"gauss") == 0){
				incremento = _rnd.Gauss(0,_passo);   //incremento gaussiano
			}

			possible.push_back(incremento+_posizioni[i][j]);   //possibile posizione futura
		}
		double prob_ratio =pow(_f->Eval(possible),2)/pow(_f->Eval(_posizioni[i]),2);
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


}

void Metropolis::Restart(){

	_posizioni.erase(_posizioni.begin(),_posizioni.end()-1);  //elimina tutti i passi fatti tranne la posizione finale ottenuta
	_accepted = 0.;
	for(int j = 0; j<_Nbins; j++){
			_walker[j] = 0;
		}

}

void Metropolis::Equilibrate(int nstep){ //equilibra l'algoritmo in nstep
	int passi = _N;  //così non perdo il numero di passi con cui voglio eseguire la simulazione (numero che ho impostato nel costruttore)
	_N = nstep;
	Run();
	Restart();
	_N = passi;
}

double Metropolis::Integrate(FunzioneBase *hamiltonian) {
	double sum = 0;
	for(int i = 0; i<_N; i++){
		sum += hamiltonian->Eval(_posizioni[i]);  //valuta la funzione hamiltoniana con i punti campionati
	}
	return sum/(double)_N;
}


void Metropolis::Print(const char * namefile) {
	ofstream Risultati;
	Risultati.open(namefile,ios::app);

	for(unsigned int i=0;i<_posizioni.size();i++){
		for(int j = 0; j<_Ndim; j++) {
			Risultati <<_posizioni[i][j]<<endl;

		}

	}
	Risultati.close();

}

void Metropolis::UpdateHistogram() {
	double binsize = abs(_b-_a)/(double)_Nbins; //dimensione di ciascun bin
	int index;  //mi indica in quale posizione del bin mi trovo
	double sum = 0.;
	for(int i =0; i<_N; i++){
		if(abs(_posizioni[i][0])<_b){
			index = int((_b +_posizioni[i][0])/binsize);
			_walker[index] += 1;
			sum += 1.;
		}
	}
	Normalization(sum);
	_hist.push_back(_walker);


}

void Metropolis::Normalization(double sum) {
	double binsize = abs(_b-_a)/(double)_Nbins;
	for(int i=0; i<_Nbins; i++) {
		_walker[i] /=(sum*binsize);
	}

}

void Metropolis::PrintHistogram(){
	for(unsigned int i=0; i<_hist.size(); i++){
		for(unsigned int j = 0; j<_walker.size(); j++){
			cout<<_hist[i][j]<<" ";
		}
		cout<<endl;

	}

}

void Metropolis::BlockingHistogram(int Nblock, const char * namefile) {
	double binsize = abs(_b-_a)/(double)_Nbins;
	AnalisiIstogramma(Nblock, _Nbins, binsize, _a, _hist, namefile);

}
