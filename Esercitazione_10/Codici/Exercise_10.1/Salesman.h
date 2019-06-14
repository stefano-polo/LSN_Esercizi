#ifndef _Salesman_h
#define _Salesman_h

#include "random.h"
#include <vector>


using namespace std;

class Salesman{

private:

  Random _rnd;
  int _NumberOfCity;  //numero di città
  int _NumberOfParents; //numero di cromosomi per generazione(iterazione del metropolis)

  vector<vector<int>> _StringsOfCities;  //matrice genitori
  vector<vector<double>> _Positions;  //matrice delle posizioni (x,y) di ciascuna città
  vector<double> _FitnessDistances; //calcolo della distanza percorsa per ciascuna stringa
  vector<vector<int>> _NewGeneration;  //nuove stringhe generate dal metropolis
  int _best_index;  //indice della stringa con miglior distanza per ogni passo del metropolis
  int _n_generation;  //numero di passi campionamento del metropolis
  vector<double> _remember_best_distance; //vettore contenenti le migliori distanze
  int _accepted;//variabile per accpetance ratio del Metropolis


public:

  Salesman(int NumberOfCity, int NumberOfParents, int n_generation); //Numero di città, numero di stringhe, numero di passi di campionamento metropolis
  ~Salesman();
  void Square(void);  //genero posizioni città in un quadrato
  void Circle(void); //genero posizioni città su circonferenza
  void Fitness(void); //calcola la distanza su matrice _StringsOfCities
  void Print(void); //stampa a video la matrice _StringsOfCities
  void Best_path(const char * file_name); //scrive su file il miglior percorso (posizioni x e y)
  void Best_string(void); //valuta la stringa con la distanza migliore (riempie l'indice _best_index)
  void Best_distance(const char * file_name); //stampa su file le migliori distanze in funzione di beta
  vector<double> Fitness_Annealing(vector<vector<int>>); //è come Fitness solo che la calcola per qualsiasi matrice e in più restituisce un vettore di distanze
  void Mutation_Annealing(int number_mutation);  //mutazioni per il metropolis (fissato numero di mutazioni)
  void PrintAccepted(); //stampo acceptance ratio
  void Metropolis(double beta, int number_pair_mutation);  //metodo del simulated annealing
  void Run_Metropolis(double beta, int number_pair_mutation); //run del simulated annealing
  void Equilibration(int numero_equilibrazioni,double beta, int number_pair_mutation); //equilibrazione metropolis

};

#endif
