#include <algorithm>
#include "random.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>
#include "emcee.h"
using namespace std;

Walker::Walker(double a):Random(){
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
	

	_x=0;
	_y=0;
	_z=0;
	_a=a;
	
	
}

Walker::~Walker(){
	_rnd.SaveSeed();
}

void Walker::Passo(){
	int r = _rnd.Rannyu(0,6);
	if(r==0){_x++;}
	if(r==1){_y++;}
	if(r==2){_z++;}
	if(r==3){_x--;}
	if(r==4){_y--;}
	if(r==5){_z--;}
	
}

double Walker::RndWalk(int N){   //passo discreto
	for(int i=0;i<N;i++){
	Passo();
	}
	double distanza =_a*sqrt(pow(_x,2)+pow(_y,2)+pow(_z,2));
	Restart();  //riparte dalla posizione iniziale
	return distanza;
}

void Walker::Isos(){   //passo continuo
	double theta = _rnd.Rannyu(0,M_PI);
	double phi = _rnd.Rannyu(0,2.*M_PI);
	_x+=sin(theta)*cos(phi);
	_y+=sin(theta)*sin(phi);
	_z+=cos(phi);
}

double Walker::Diffusione(int N){
	for(int i=0;i<N;i++){
	Isos();
	}
	double distanza = _a*sqrt(pow(_x,2)+pow(_y,2)+pow(_z,2));
	Restart();
	return distanza;
}
