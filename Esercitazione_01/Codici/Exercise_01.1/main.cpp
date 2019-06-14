#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include "TLC.h"

using namespace std;

int main (int argc, char *argv[]){

  Random rnd;
  SetRandomGenerator(&rnd);  //inzializzo il seme del generatore di numeri casuali


  int M=1e4;  //numero di dati
  int N=100;  //numero di blocchi
  vector<double> s(M);

  for(int j=0; j<M; j++){
	    s[j]=rnd.Rannyu();
  }


  TLC(N,s,"../../Risultati/Risultati01.1.1.out");// punto 1 dell'esercizio
  TLC_modified(N,s,"../../Risultati/Risultati01.1.2.out"); //punto 2 dell'esercizio


  //punto 3 dell'esercizio: Test del chi quadro

  M = 100;  //numero di intervalli in cui divido l'intervallo [0,1]
  int dim = 1e6;  //numero di numeri casuali da generare
  vector<double> chi2(M);
  vector<double> dati(dim);

  for(int i=0; i<dim; i++){
	   dati[i] = rnd.Rannyu();  //carico il vettore con 1e6 numeri generati uniformemente tra [0,1)
  }
  Chi2(M, dim, dati,"../../Risultati/Risultati01.1.3.out");  //eseguo il test del chi2
  rnd.SaveSeed();
  return 0;
}
