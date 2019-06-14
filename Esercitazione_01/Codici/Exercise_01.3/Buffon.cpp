#include <algorithm>
#include "random.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "Buffon.h"
using namespace std;

Buffon::Buffon(double L, double d):Random(){
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


	_L = L; //lunghezza ago
	_d = d; //distanza tra le righe


}

Buffon::~Buffon(){
	_rnd.SaveSeed();

}

double Buffon::Sin_Angle(){    //funzione per valutare il seno dell'angolo dell'ago che forma con la linea
	double x;
	double y;


	do{
		x = _rnd.Rannyu(-1,1);
		y = _rnd.Rannyu(0,1);

	}
	while((pow(x,2)+pow(y,2)) >= 1);  //esce dal ciclo quando la condizione di controllo è falsa ossia quando è dentro la circonferenza

	return y/sqrt(pow(x,2)+pow(y,2) );
}

double Buffon::Experiment(int N_throws){
	double center;
	double sen;
	int N_hit=0;
	double y1,y2;
	for(int i=0; i<N_throws; i++){
		center = _rnd.Rannyu(-_d*0.5,_d*0.5);   //distanza del centro dalla linea
		sen = Sin_Angle();
		y1 = center+sen*_L*0.5;  //estremi dell'ago
		y2 = center-sen*_L*0.5;
		if((y1>=0 && y2<=0) || (y1<=0 && y2>=0)){ //controllo se
			N_hit++;
		}
	}

	return (2.*_L*N_throws)/(N_hit*_d);
}
