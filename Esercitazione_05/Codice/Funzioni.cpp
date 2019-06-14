#include <iostream>
#include <cmath>
#include "Funzioni.h"
#include <cfloat>
#include <fstream>
#include <vector>


#define a0 1//0.0529e-9


using namespace std;

double S_uno::Eval (vector<double> v) const {     //v[0]=x v[1]=y v[2]=z
	double r = sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2) );	
	return pow(((pow(a0,-3./2.)/sqrt(M_PI)))*exp(-r/a0),2);
} //ricordati di fare il modulo quadro

double P_due::Eval(vector<double> v) const{
	double r = sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2) );	
	return pow((pow(a0,-5./2.)/8.)*sqrt(2./M_PI)*v[2]*exp(-0.5*(r/a0)),2); 	//cos = z/r;
}




