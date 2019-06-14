#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Salesman.h"
#include <algorithm>

using namespace std;

int main (){

	Salesman Walker(30,1,1000);  //1000 campionamenti metropolis, 30 città, 1 stringa
	Walker.Circle(); //genero su cerchi
	int number_pair_mutation = 1;  //numero di mutazioni per il metropolis
	cout<<"Città su circonferenza di raggio 1"<<endl;
	for(double beta = 0; beta<60; beta += 0.001){
		Walker.Run_Metropolis(beta, number_pair_mutation);
	}
	Walker.Best_path("../../Risultati/Exercise_10.1/path_circ.dat"); //stampo su file il miglior percorso
  Walker.Best_distance("../../Risultati/Exercise_10.1/distance_circ.dat"); //stampo su file le distanze in funzione di beta

	Salesman Walker2(30,1,1000);
	Walker2.Square();
	number_pair_mutation = 1;
	cout<<"Città dentro quadrato di lato 1"<<endl;
	for(double beta = 0; beta<60; beta += 0.005){
		Walker2.Run_Metropolis(beta, number_pair_mutation);
	}
	Walker2.Best_path("../../Risultati/Exercise_10.1/path_sq.dat");
  Walker2.Best_distance("../../Risultati/Exercise_10.1/distance_sq.dat");
 	return 0;

}
