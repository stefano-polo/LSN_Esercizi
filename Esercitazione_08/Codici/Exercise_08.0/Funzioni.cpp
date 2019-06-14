#include <iostream>
#include <cmath>
#include "Funzioni.h"
#include <cfloat>
#include <fstream>
#include <vector>

using namespace std;

Psi::Psi(double mu, double sigma){
	_mu = mu;
	_sigma = sigma;

}

Psi::~Psi()  {		//distruttore

} 

double Psi::Eval (vector<double> v) const {     //v[0]=x 
	double pezzo1 = exp(-pow(v[0]-_mu,2)/(2.*pow(_sigma,2)));
	double pezzo2 = exp(-pow(v[0]+_mu,2)/(2.*pow(_sigma,2)));
	return pezzo1+pezzo2;
} //Il modulo quadro lo eseguo nel Metropolis


Hamiltoniana::~Hamiltoniana(){

}

double Hamiltoniana::Eval (vector<double> v) const {  //metodo che calcola (H*psi)/psi
	double mu = _psi->GetMu();
	double sigma = pow(_psi->GetSigma(),2);
	double Potential = pow(v[0],4) - (5./2.)*pow(v[0],2);
	double e1 = pow(v[0]-mu, 2)/sigma;
	double e2 = pow(v[0]+mu, 2)/sigma;
	
	double psi_sec = ( -_psi->Eval(v) + e1 * exp(-e1/2.) + e2 * exp(-e2/2.) )/sigma;
	
	return -0.5 * psi_sec/_psi->Eval(v) + Potential;

	
}





