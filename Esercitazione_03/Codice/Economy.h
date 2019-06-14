#ifndef _Economy_h_
#define _Economy_h_


#include "random.h"


using namespace std;

class Black_Scholes:public Random {
private:
	Random _rnd;
	double _sigma;
	double _K;
	double _S0;
	double _r;
	double _St;




	public:

	Black_Scholes(double S0, double K, double r, double sigma); //preparazione dell'esperimento
	~Black_Scholes();
	void Restart_Price(); //faccio ripartire il prezzo ad _S0
	void Price(double ti, double tf); //calcolo il prezzo dato un intervallo di tempo
	vector<double> Call_Put_Option(double ti, double tf, int passi); //calcolo call option, e put option in un intervallo di tempo dato un numero di passi
};


#endif
