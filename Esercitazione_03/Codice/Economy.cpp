#include <algorithm>
#include "random.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "Economy.h"
using namespace std;

Black_Scholes::Black_Scholes(double S0, double K, double r, double sigma):Random(){
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


	_sigma = sigma;
	_K = K;
	_S0 = S0;
	_St = S0;
	_r = r;

}

Black_Scholes::~Black_Scholes(){
	_rnd.SaveSeed();
}

void Black_Scholes::Price(double ti, double tf) {   //calcola il spot price al tempo tf a partire da quello in ti
	_St= _St*exp( (_r-0.5*pow(_sigma,2))*(tf-ti)+ _sigma*_rnd.Gauss(0.,1)*sqrt(tf-ti) );
}

void Black_Scholes::Restart_Price() {   //funzione che mi restarta il prezzo al tempo zero
	_St = _S0;
}

vector<double> Black_Scholes::Call_Put_Option(double ti, double tf, int passi){  //metodo che mi calcola la call e put option
	vector<double> v(2);
	double L = (tf-ti)/(double(passi));  //dimensione degli intervallini temporali
	for(int j=0;j<passi;j++){
		Price(ti+j*L,ti+(j+1)*L);
	}
	v[0] = exp(-_r*(tf-ti))*max(0.,_St-_K);  //call option europea
	v[1] = exp(-_r*(tf-ti))*max(0.,_K-_St);  //put option europea
	Restart_Price();  //riparto dal prezzo al tempo zero
	return v;

}
