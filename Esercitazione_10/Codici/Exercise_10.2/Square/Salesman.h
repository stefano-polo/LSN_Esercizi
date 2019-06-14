#ifndef _Salesman_h
#define _Salesman_h

#include "random.h"
#include <vector>


using namespace std;

class Salesman{

private:

  Random _rnd;
  int _NumberOfCity;//numero di città
  int _NumberOfParents; //numero di cromosomi per generazione(iterazione del metropolis)
  int _seed;//seme per il random search

  vector<vector<int>> _StringsOfCities;  //matrice genitori
  vector<vector<double>> _Positions;  //matrice delle posizioni di ciascuna città (x,y)
  vector<double> _FitnessDistances; //calcolo della distanza percorsa per ciascuna stringa
  vector<vector<int>> _NewGeneration; //nuove stringhe generate dal metropolis
  int _best_index;//indice della stringa con miglior distanza per ogni passo del metropolis
  int _n_generation;   //numero di passi campionamento del metropolis
  vector<double> _remember_best_distance; //vettore contenenti le migliori distanze
  int _accepted; //variabile per accpetance ratio del Metropolis
  vector<vector<int>> _Best_People;  //Matrice delle migliori stringhe per ogni passo metropolis


public:

  Salesman(int NumberOfCity, int NumberOfParents, int n_generation,int seed); //Numero di città, numero di stringhe, numero di passi di campionamento metropolis
  ~Salesman();
  void Square(void);  //genero posizioni città in un quadrato
  void Circle(void); //genero posizioni città su circonferenza
  void Fitness(void); //calcola la distanza su matrice _StringsOfCities
  void Print(void); //stampa a video la matrice _StringsOfCities
  void Print_Positions(void); //stampo a video le posizioni
  void Best_path(const char * file_name); //stampo su file il miglior cammino
  void Best_string(void); //valuta la stringa con la distanza migliore (riempie l'indice _best_index)
  void Best_distance(const char * file_name); //stampa su file le migliori distanze in funzione del passo metropolis
  void Restart(void);//restart dell'algoritmo
  vector<double> Fitness_Annealing(vector<vector<int>>);//è come Fitness solo che la calcola per qualsiasi matrice e in più restituisce un vettore di distanze
  void Mutation_Annealing(int number_mutation);//mutazioni per il metropolis (fissato numero di mutazioni)
  void PrintAccepted(); //stampo acceptance ratio
  void Metropolis(double beta, int number_pair_mutation); //metodo del simulated annealing
  void Run_Metropolis(double beta, int number_pair_mutation);//run del simulated annealing
  void Equilibration(int numero_equilibrazioni,double beta, int number_pair_mutation);
  void Random_Search(void); //permutazione delle stringhe
  double Minimum_distance(void); //restituisce la distanza minima registrata in una run
  vector<vector<double>> Get_Positions(void); //restituisce il vettore delle posizioni delle città
  int* Best_Cammino(void); //restituisce il cromosoma con il percorso migliore di una run

};

#endif
