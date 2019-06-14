#include "Integrazione.h"
#include <algorithm>
#include "random.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>
using namespace std;



Integral::Integral(double a, double b, FunzioneBase * function): Random() {
	
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

	_integrand = function;
	_a = min(a,b);
	_b = max(a,b);
		
	if (a > b) _sign = -1;
	else _sign = 1;
	
	
}

Integral::~Integral () { 
	_rnd.SaveSeed();
} 


double Integral::mediaUn(int N_estrazioni) {  //integrazione con il metodo della media
	_sum = 0;
	double media = 0;
	double intervallo = _b-_a;
	
	
	for(int i = 1; i<= N_estrazioni; i++) {
	double x = _rnd.Rannyu(_a, _b);
	double f = _integrand->Eval(x);
	_sum+=f;
	}
	media = (_sum)/(1.*N_estrazioni);

	_integral = media*intervallo; 

	return _sign*_integral;

}

double Integral::mediaRe(int N_estrazioni) {  //integrazione usando l'importance sampling
	_sum = 0;
	double media = 0;
	double intervallo = _b-_a;
	
	
	for(int i = 1; i<= N_estrazioni; i++) {
	double x = _rnd.Retta();	//estraggo numeri distribuiti secondo la distribuzione scelta per l'importance sampling
	double f = _integrand->Eval(x);
	_sum+=f;
	}
	media = (_sum)/(1.*N_estrazioni);

	_integral = media*intervallo; 

	return _sign*_integral;

}
