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
	int _N;  //numero di iterazioni dell'algoritmo
	double _passo; //passo di movimento dell'algoritmo
	vector<vector<double>> _posizioni; //Matrice _Nx_Ndim per salvare le posizioni (x,y,z)
	int _Ndim;  //numero di dimensioni dello spazio dei parametri(in questo caso 3)
	int _accepted; //contatore per il calcolo del acceptance ratio

	public:

	Metropolis(vector<double> posizioni0, int N, double passo,int N_dimension,FunzioneBase * f);
	~Metropolis();
	void Run(const char * metodo);  //algoritmo con passo uniforme o gaussiano a seconda della scelta (metodo)
	void Restart(); //restart del campionamento alla posizione iniziale
	void Equilibrate(int nstep, const char * metodo); //metodo di equilibrazione dell'algoritmo
	void PrintAccepted(); //stampa a video il valore dell'acceptance ratio
	void SetPasso(double passo){_passo=passo;} //impostare l'ampiezza del passo dell'algoritmo
	void Print(const char * namefile); //stampare su file i punti campionati
	vector<double> Raggi(); //restituisce il vettore dei raggi calcolati con i punti x,y,z campionati


};











#endif //_Metropolis_h_
