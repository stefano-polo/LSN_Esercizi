#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "Buffon.h"
#include <cmath>
#include <vector>
#include  "TLC.h"
using namespace std;

int main (int argc, char *argv[]){
	int N_blocks = 100;    //numero dei blocchi
	int N_throws = 1e4;  //numero di esperimenti

	double d = 6.;  //lunghezza dell'ago
	double L = 4.; //distanza tra le righe 
	vector<double> pi(N_throws);


	Buffon b(L,d); //preparo esperimento

	for(int i=0; i<N_throws; i++){
		pi[i] = b.Experiment(3e4);  //ciascun esperimento è eseguito con 5e4 lanci

	}

	TLC(N_blocks,pi,"../../Risultati/Risultati01.3.out");


	return 0;
}
