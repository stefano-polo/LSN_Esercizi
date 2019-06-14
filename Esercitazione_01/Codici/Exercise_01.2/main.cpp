#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <vector>
#include <cmath>
#include "Statistica.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   SetRandomGenerator(&rnd);
   ofstream Risultati2;
   Risultati2.open("../../Risultati/Risultati01.2.out");

   int M=1e4;  //numero di lanci
   int F[4] = {1,2,10,100}; //le varie configurazioni che voglio simulare: 1 dado, 2 dadi, etc..
   vector<double> dado_mean(M); //lanci uniformi
   vector<double> exp_mean(M);  //lanci esponenziali
   vector<double> lor_mean(M);  //lanci lorentziani


   for(int m=0;m<4;m++){  //ciclo sulla dimensione del vettore F[4]--> ciclo su quanti dadi lancio
     vector<double> lor(F[m]);  //preparo vettori di appoggio con la dimensione corrispondente al numero di dadi da lanciare
     vector<double> exp(F[m]);
     vector<double> s(F[m]);
     for(int i=0;i<M;i++){  //ciclo sul numero di lanci
   	    for(int j=0; j<F[m]; j++){
            s[j]=int(rnd.Rannyu(1,7));  //prendo numeri interi uniformi tra [1,6]---> il 6 in questo caso è incluso
            exp[j]=rnd.Exponential(1.);  //numeri esponenziali
		        lor[j]=rnd.Lorentz(1.,0.);  //numeri gaussiani
   	     }
      dado_mean[i]=Media(s);  //faccio la media dei lanci dei F[m] dadi
      exp_mean[i]=Media(exp);
      lor_mean[i]=Media(lor);

      if (Risultati2.is_open()){
   		   Risultati2 <<dado_mean[i]<< " "<< exp_mean[i] << " " <<lor_mean[i]<< endl;
       }
	     else cerr << "PROBLEM: Unable to open random.out" << endl;

  }

  }

    Risultati2.close();

   rnd.SaveSeed();
   return 0;
}
