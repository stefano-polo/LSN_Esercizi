#ifndef _Metropolis_h_
#define _Metropolis_h_

#include "random.h"
#include "Funzioni.h"
#include <vector>

using namespace std;

class Metropolis{
private:
	Random _rnd;
	FunzioneBase * _f;   //distribuzione da campionare
	int _N;  //numero di passi dell'algoritmo
	double _passo; //passo di movimento dell'algoritmo
	vector<vector<double>> _posizioni; //Matrice _Nx_Ndim per salvare le posizioni (x,y,z) 
	int _Ndim;  //numero di dimensioni dello spazio dei parametri(in questo caso 3)
	int _accepted; //Numero di passi accettati durante il Metropolis
	const char * _metodo; //metodo con cui faccio il passo dell'algoritmo (uniforme o gaussiano)
	public:
		
	Metropolis(vector<double> posizioni0, int N, double passo,const char * metodo, int N_dimension,FunzioneBase * f);
	~Metropolis();
	void Run();  //algoritmo con passo uniforme o gaussiano a seconda della scelta (metodo)
	void Restart();
	void Equilibrate(int nstep);
	void SetPasso(double passo){_passo=passo;}
	void Print(const char * namefile); //stampa i punti campionati
	double Integrate(FunzioneBase *hamiltonian); //calcola l'integrale della funzione hamiltonian
	void PrintAccepted();
	
	
};











#endif //_Metropolis_h_
