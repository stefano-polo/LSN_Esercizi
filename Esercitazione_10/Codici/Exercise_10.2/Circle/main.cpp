#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Salesman.h"
#include <algorithm>
#include "mpi.h"

using namespace std;

int main (){

//inizio della programmazione parallela
MPI::Init();
int size = MPI::COMM_WORLD.Get_size();
int rank = MPI::COMM_WORLD.Get_rank();

if(size>4){cout<<"Hai scelto troppi nodi"<<endl;
return 1;}


//preparo parametri dell'esperimento in modo uguale per tutti i nodi
int NumberOfCity = 30; //numero di città
int NumberOfParents = 1; //numero di stringhe (cromosomi)
int n_generation = 150; //numero di passi dell'algoritmo Metropolis del simulated Annealing
int number_pair_mutation = 1;
//vettori di raccolta risultati di ciascun nodo
int path_irecv[size][NumberOfCity]; //matrice che raccoglie i migliori cammini di ciascun nodo
vector<double> distance_irecv(size,0);//vettore che raccoglie le migliori distanze di ciascun nodo
//preparo la configurazione di città uguale per ciascun nodo
Salesman walker(NumberOfCity, NumberOfParents, n_generation,rank);//(rank inizializza in modo diverso il seme della permutazione di ciascun nodo)
walker.Circle(); //importante che in ciascun core le città abbiamo la medesima disposizione
vector<vector<double>> Posizioni_citta = walker.Get_Positions();

//eseguo l'esperimento
for(double beta = 0.; beta<60; beta += 0.005){ //1.,20,0.0005
	walker.Run_Metropolis(beta, number_pair_mutation);
}
//prendo le distanze minime di ciascun nodo e anche il corrispondente cammino
double distanza_minima = walker.Minimum_distance();
int * best_path = new int[NumberOfCity];
best_path = walker.Best_Cammino();
//passo i dati al nodo zero

MPI_Gather(best_path,NumberOfCity,MPI_INT,path_irecv[rank],NumberOfCity,MPI_INT,0,MPI::COMM_WORLD);
MPI_Gather(&distanza_minima,1,MPI_DOUBLE,&distance_irecv[rank],1,MPI_DOUBLE,0,MPI::COMM_WORLD);


if(rank==0){
	int index_min = min_element(distance_irecv.begin(),distance_irecv.end()) - distance_irecv.begin();
  cout<<"Best distance: "<<distance_irecv[index_min]<<endl;
	//cout<<distance_irecv[0]<<endl<<distance_irecv[1]<<endl<<distance_irecv[2]<<endl<<distance_irecv[3]<<endl; // SCOMMENTARE SE SI DESIDERA VEDERE LE DISTANZE OTTENUTE DA CIASCUN NODO
	ofstream results;
    results.open("../../../Risultati/Exercise_10.2/path_circle.out");
		results<<distance_irecv[index_min]<< " "<<distance_irecv[index_min]<<endl;  //stampo come prima cosa il valore della distanza minore percorsa
	//stampo  x, y delle città in ordine di passaggio dalla prima all'ultima
    for (int i=0; i <NumberOfCity;i++){
        results<<Posizioni_citta[path_irecv[index_min][i]][1]<<" "<<Posizioni_citta[path_irecv[index_min][i]][2]<<endl;
    }
    //stampo come posizione finale la hometown
    results<<Posizioni_citta[path_irecv[index_min][0]][1]<<" "<<Posizioni_citta[path_irecv[index_min][0]][2]<<endl;


    results.close();
}

MPI::Finalize();
//fine della programmazione parallela
	return 0;

}
