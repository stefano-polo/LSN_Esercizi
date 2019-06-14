#ifndef _Salesman_h
#define _Salesman_h

#include "random.h"
#include <vector>


using namespace std;

class Salesman{

private:

  Random _rnd;
  int _NumberOfCity;  //numero di città nel sistema == numero di geni
  int _NumberOfParents; //numero di cromosomi per generazione

  vector<vector<int>> _StringsOfCities;  //matrice genitori della data generazione (riga = cromosoma; colonna = città visitata a un certo tempo)
  vector<vector<double>> _Positions;  //matrice delle posizioni di ciascuna città (x,y). la riga 0 == posizione città 0; la riga 1 == posizione città 1....
  vector<double> _FitnessDistances; //vettore della distanza percorsa per ciascuna stringa (cromosoma) contenuto nella matrice _StringsOfCities
  vector<vector<int>> _NewGeneration;  //matrice figlia a una data generazione
  int _best_index;  //indice di riga che indica il cromosoma con miglior distanza (minima) a una data generazione
  int _n_generation;  //numero di generazioni in una run
  vector<double> _remember_best_distance; //vettore contenente distanza migliore di ogni generazione
  vector <int> _breeders;  //vettore con gli indici che si riferiscono ai cromosomi (riga) fertili nella matrice _StringsOfCities
  vector<vector<double>> _best_distances50; //matrice della media e std delle migliori distanze ottenute da 50% popolazione

public:
  Salesman(int NumberOfCity, int NumberOfParents, int n_generation); //costruttore che crea generazione genitoriale
  ~Salesman();
  void Square(void); // riempie _Positions in posizioni (x,y) generate casualmente in un quadrato
  void Circle(void);// riempie _Positions in posizioni (x,y) generate casualmente su una circonferenza
  void Fitness(void); //calcolo  della distanza percorsa per ciascuna stringa (cromosoma) contenuto nella matrice _StringsOfCities
  void Crossover(int parent1, int parent2); //metodo di crossover(prede in input indici dei genitori che generano)
  void Mutation(double probability_mutation); //metodo di mutazione a una data probabilità fissata dall'esterno
  void Breeders_selection(void); //seleziona i cromosomi che si accoppiano (riempie il vettore _breeders)
  int Selection(void);  //pesca randomicamente tra i cromosomi fertili
  void Run(double probability_mutation);  //esecuzione algoritmo
  void Print(void); //stampa a schermo la matrice _StringsOfCities
  void Update_Distance50(int); //calcolo media e std delle migliori distanze del 50% popolazione
  void Best_path(const char * file_name); //scrive su file il percorso migliore tra tutte le generazioni
  void Best_string(void); //valuta il cromosoma con miglior distanza (minima) a una data generazione (riempie _best_index)
  void Best_distance(const char * file_name); //scrive su file le migliori distanze di ogni generazione
  void Print_Best_distance50(const char * file_name);//scrive su file le migliori distanze del 50% popolazione di ogni generazione con annesso errore
  void Restart(void);
};

#endif
