#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "random.h"
#include "Funzioni.h"
#include "TLC.h"
#include "Metropolis.h"

using namespace std;

int main (int argc, char *argv[]){

	int N_dimensioni = 1;//campiono solo la x
	double passo = 2.5; //passo metropolis
	int N_passi = 1e3;  //numero di passi con cui eseguire il campionamento di psi con il Metropolis
	int N_equilibrio = 1000; //numero di passi con cui equilibrare il Metropolis
	int N_integrali = 1e4;  //numero di integrali da calcolare
	vector<double> posizione0(1); //posizione iniziale da cui parte l'algoritmo
	posizione0[0] = 0.;
	double mu = 0.;
	double sigma = 0.;
	Psi * psi = new Psi(mu,sigma);
	FunzioneBase * h = new Hamiltoniana(psi);
	Metropolis walker(posizione0, N_passi, passo, "uniform", N_dimensioni, psi);
	vector<double> integral_result(N_integrali); //vettore risultati integrali

	walker.Equilibrate(N_equilibrio);  //equilibrio l'algoritmo

	ofstream Parameters;
	Parameters.open("../../Risultati/ottimizzazione_parametri.out");
	for(double mi=0.795; mi<=0.81; mi+=0.001){
		psi->SetMu(mi);
		cout<<mi<<endl;
		for(double sigma=0.61; sigma<=0.62; sigma+=0.001){

			psi->SetSigma(sigma);
			for(int i=0; i<N_integrali; i++){
				walker.Run();      //campionamento dei punti
				integral_result.at(i) = walker.Integrate(h); //calcolo dell'integrale
				//if(i%N_integrali==0) {walker.PrintAccepted();}
				walker.Restart();
			}

			Parameters << mi << setw(10) << sigma << setw(10) << Media(integral_result) << endl;
		}
	}

	Parameters.close();
	return 0;
}
