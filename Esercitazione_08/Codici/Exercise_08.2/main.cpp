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

	//Parametri per algoritmo metropolis e calcolo dell'integrale
	int N_block = 100;  //numero di blocchi per il data blocking
	int N_dimensioni = 1; //ho da campionare solo una variabile x col Metropolis
	double passo = 2.5;  //Passo del metropolis per 50% di accettazione
	int N_passi = 1e4;  //numero di passi con cui eseguire il campionamento di psi con il Metropolis
	int N_equilibrio = 1000; //numero di passi con cui equilibrare il Metropolis
	int M_integrali = 1e4;  //numero di integrali da calcolare
	vector<double> posizione0(1); //posizione iniziale da cui parte l'algoritmo
	posizione0[0] = 0.;  //posizione iniziale algoritmo
	double mu = 0.797;
	double sigma =0.614;
	Psi * psi = new Psi(mu,sigma);  //distribuzione di probailità da campionare
	FunzioneBase * h = new Hamiltoniana(psi); //funzione per il calcolo dell'integrale

	//Parametri per istogramma della |psi|^2

	double a = -3.; //estremo sinistro dell'intervallo dell'istogramma
	double b = 3.; //estremo dx
  //il numero di bin si fissa da Metropolis.h
	Metropolis walker(posizione0, N_passi, passo, "uniform", N_dimensioni, psi,a,b);
	vector<double> integral_result(M_integrali); //vettore risultati integrali

	walker.Equilibrate(N_equilibrio);  //equilibrare l'algoritmo

	for(int i=0; i<M_integrali; i++){
		walker.Run();      //campionamento dei punti
		integral_result.at(i) = walker.Integrate(h); //calcolo dell'integrale
		if(i%1000==0) {
			walker.PrintAccepted();	  //ogni 1000 calcoli dell'integrale stampo l'acceptance ratio
			walker.Print("../Risultati/Psi.out");
		}
		walker.UpdateHistogram();  //faccio un istogramma dei punti campionati
		walker.Restart();  //resetto i punti campionati e ricalcolo l'integrale
	}

	TLC( N_block, integral_result, "../../Risultati/Integral.out");  //datablocking su integrali
	walker.BlockingHistogram(N_block, "../../Risultati/Histogram_psi.out"); //datablocking su istrogramma
	return 0;
}
