#include <iostream>
#include <cmath>
#include "Funzioni.h"
#include <cfloat>
#include <fstream>

using namespace std;

double Coseno::Eval (double x) const {
	return (0.5*M_PI)*cos(0.5*x*M_PI);
}

double G::Eval(double x) const{  //nuova funzione per il calcolo dell'integrale con importance sampling
	return ((0.5*M_PI)*cos(0.5*x*M_PI))/(2*(1-x));
}
