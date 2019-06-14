#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include "Economy.h"
#include <cmath>
#include "TLC.h"

using namespace std;

int main (int argc, char *argv[]){


	double S0 = 100.;  //prezzo iniziale
	double K = 100.;
	double T = 1.;  //tempo finale
	double t = 0.;  //tempo iniziale
	double r = 0.1;
	double sigma = 0.25;
	int N = 100;   //numero di blocchi
	int M = 1e4;   //numero di prezzi (numero di lanci della simulazione per entrambi i casi)

	Black_Scholes DB(S0, K, r, sigma);  //preparo l'esperimento con i parametri fissati
	vector<double> dati(2);  //vettore di appoggio
	vector<double> C_Dr(M); //call option diretta
	vector<double> C_Ds(M); //call option discreta
	vector<double> P_Dr(M);	//pull option diretta
	vector<double> P_Ds(M); //pull option discreta


	for(int i=0; i<M;i++){
		dati = DB.Call_Put_Option(t, T,1);  //eseguo il calcolo per il caso diretto (mi restituisce un vettore a due componenti con call in 0 e put in 1) il caso diretto equivale a fissare il numero di passi a 1
		C_Dr[i] = dati[0];
		P_Dr[i] =dati[1];
		dati = DB.Call_Put_Option(t, T,100); //eseguo il caso discreto fissando il numero di passi a 100
		C_Ds[i] = dati[0];
		P_Ds[i] =dati[1];

	}

	TLC(N, M, C_Dr,"../Risultati/Risultati_diretto_Call.out");  //data blocking per ogni caso----> la funzione si trova in TLC.cpp
	TLC(N, M, P_Dr,"../Risultati/Risultati_diretto_Put.out");
	TLC(N, M, C_Ds,"../Risultati/Risultati_discreto_Call.out");
	TLC(N, M, P_Ds,"../Risultati/Risultati_discreto_Put.out");






}
