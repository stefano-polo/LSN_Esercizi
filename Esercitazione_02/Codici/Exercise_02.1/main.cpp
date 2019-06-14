#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "Funzioni.h"
#include "Integrazione.h"
#include "TLC.h"
#include <cmath>
#include <vector>

using namespace std;

int main (int argc, char *argv[]){

   FunzioneBase *f = new Coseno();
   Integral * integrale = new Integral(0.,1.,f); //Integral(estremo sinistro, estremo destro, funzione da valutare)

   int M=1e4;     //numero di integrali eseguiti
   int N=100;    //numero di blocchi
   
   vector<double> s(M);  //vettore su cui calcolo integrale
   for(int j=0; j<M; j++){
	    s[j]=integrale->mediaUn(M); //calcolo M integrali usando M variabili di sampling
   }

   TLC(N,s,"../../Risultati/Exercise_02.1/Risultati02.1.out");

   //punto due dell'esercizio 1--->IMPORTANCE SAMPLING
   FunzioneBase *g = new G();  //nuova funzione per importance sampling
   Integral * calcolo = new Integral(0.,1.,g);
   for(int j=0; j<M; j++){
	    s[j]=calcolo->mediaRe(M);
   }
   TLC(N,s,"../../Risultati/Exercise_02.1/Risultati02.2.out");
   return 0;
}
