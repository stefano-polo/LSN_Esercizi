#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "emcee.h"
#include "TLC.h"
#include "random.h"


using namespace std;

int main (int argc, char *argv[]){

	int N=100;  //numero di blocchi == numero di passi del singolo random walk
	int M=1e4;  //numero dimensione vettori (1e4 perche voglio 1e4 random walk per ciascun passo)

	vector<double> d_2_disc(M); //distanza per il rndwalk discreto
	vector<double> d_2_cont(M); //distanza per il rndwalk continuo

	Walker w(1); //1 = passo del reticolo
	vector<vector<double>> Risultati_discreto;
  vector<vector<double>> Risultati_continuo;


	for(int i=0; i< N;i++){  //ciclo sul numero di passi del singolo random walk
      for(int j=0; j<M; j++){  //ciclo fino a M
			    d_2_disc[j]=pow(w.RndWalk(i),2); //valuto distanza quadratica rndwalk
					d_2_cont[j]=pow(w.Diffusione(i),2);
      }
			if(i%20==0)cout<<"Passo "<<i<<endl;
			Risultati_discreto.push_back(TLC(N,d_2_disc)); //eseguo datablocking con 100 blocchi
			Risultati_continuo.push_back(TLC(N,d_2_cont));
	}

	Print(Risultati_discreto,"../../Risultati/Exercise_02.2/Risultati02.2_discreto.out");
	Print(Risultati_continuo,"../../Risultati/Exercise_02.2/Risultati02.2_continuo.out");



	return 0;
}
