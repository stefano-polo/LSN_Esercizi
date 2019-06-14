#ifndef _Integrator_h_
#define _Integrator_h_

#include "Funzioni.h"
#include "random.h"

class Integral:public Random {
private:
	Random _rnd;
	double _a, _b;
	double _sum;
	int _sign;
	double _integral;
	FunzioneBase * _integrand;



	public:
	Integral(double a, double b, FunzioneBase * function);
	~Integral ();
	double mediaUn(int N_estrazioni);  //calcolo integrale con variabili di sampling ottenute da uniforme
	double mediaRe(int N_estrazioni); //calcolo integrale con variabili di sampling ottenute dalla retta

};




#endif
