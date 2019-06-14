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
	double _passo; //passo di tansizione dell'algoritmo
	vector<vector<double>> _posizioni; //Matrice _Nx_Ndim per salvare le posizioni (x)
	int _Ndim;  //numero di dimensioni dello spazio dei parametri(in questo caso 1)
	int _accepted; //Numero di passi accettati durante il Metropolis
	const char * _metodo; //metodo con cui faccio il passo dell'algoritmo (uniforme o gaussiano)

	//Variabili per l'istogramma della |psi|^2
	int _Nbins = 200; //numero di bin dell'istogramma della psi^2
	double _a, _b; //estremi dell'intervallo che divido in bin
	vector<double> _walker;
	vector<vector<double>> _hist; //vettore dell'istogramma

	public:
	Metropolis(vector<double> posizioni0, int N, double passo,const char * metodo, int N_dimension, FunzioneBase * f, double a, double b);
	~Metropolis();
	void Run();  //algoritmo con passo uniforme o gaussiano a seconda della scelta (metodo)
	void Restart();  //elimina tutti i punti campionati eccetto quello di partenza
	void Equilibrate(int nstep); //fase di equilibrazione per nstep
	void SetPasso(double passo){_passo=passo;} //settare dimensione del passo algoritmo
	void Print(const char * namefile); //stampa i punti campionati su file
	double Integrate(FunzioneBase *hamiltonian); //calcola l'integrale della funzione hamiltonian
	void PrintAccepted(); //stampa a video l'acceptance ratio
	void UpdateHistogram();//calcolo un nuovo istrogramma
	void BlockingHistogram(int Nblock, const char * filename); //data blocking su istogramma
	void PrintHistogram(); //stampo a video l'istogramma
	void Normalization(double somma);  //normalizzazione per avere istogramma normalizzato


};











#endif //_Metropolis_h_
