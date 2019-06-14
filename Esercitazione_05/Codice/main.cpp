#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "random.h"
#include "Funzioni.h"
#include "TLC.h"
#include "Metropolis.h"

using namespace std;

int main (int argc, char *argv[]){
	int N_dimensioni = 3; //campiono x,y,z
	double passo; //passo Metropolis
	int N_passi = 1e6; //passi di campionamento
	int N_equilibration = 1000; //passi di equilibrazione
	vector<double> r; //vettore che raccoglie i raggi
	cout<<endl<<endl<<"Codice lavora in unità del raggio di Bohr r = 0.0529 nm"<<endl<<endl;
	//orbitale 1s
	cout<<"ORBITALE 1S"<<endl;
	passo = 1.2;
  cout<<"Passo uniforme a = "<<passo<<endl;
	FunzioneBase * f1 = new S_uno(); //distribuzione di probabilità da campionare
	vector<double> posizioni0(3,1);  //posizione iniziale (1,1,1)
	Metropolis algoritmo(posizioni0, N_passi, passo, N_dimensioni, f1);
  algoritmo.Equilibrate(N_equilibration,"uniform"); //equilibrazione con passo uniforme
	algoritmo.Run("uniform"); //esecuzione algoritmo con passo uniforme
  algoritmo.PrintAccepted(); //stampo a video acceptance ratio
	algoritmo.Print("../Risultati/../Risultati/Risultati_S.out"); //stampo su file i punti campionati x,y,z
	r = algoritmo.Raggi(); //calcolo i raggi
	TLC(100, r, "../Risultati/Raggi_s.out"); //100 sono i numeri di blocchi per il data bloching
	passo = 0.75; //passo per l'algoritmo gaussiano in modo da fare 50% di passi
	cout<<"Passo gaussiano a = "<<passo<<endl;
	algoritmo.SetPasso(passo);
	algoritmo.Equilibrate(N_equilibration,"gauss");
	algoritmo.Run("gauss");
  algoritmo.PrintAccepted();
	r = algoritmo.Raggi();
	TLC(100, r, "../Risultati/Raggi_s_gauss.out");

  cout<<endl<<endl;
	//orbitale 2p
	cout<<"ORBITALE 2P"<<endl;
	passo = 2.8;
  cout<<"Passo uniforme a = "<<passo<<endl;
	FunzioneBase * f2 = new P_due();
	posizioni0={2.,3.,2.};
	Metropolis algoritmo2(posizioni0, N_passi, passo, N_dimensioni, f2);
  algoritmo2.Equilibrate(N_equilibration,"uniform");
	algoritmo2.Run("uniform");
  algoritmo2.PrintAccepted();
	algoritmo2.Print("../Risultati/Risultati_P.out");
	r = algoritmo2.Raggi();
	TLC(100, r, "../Risultati/Raggi_p.out");
	passo = 1.85; //passo per l'algoritmo gaussiano in modo da fare 50% di passi
	cout<<"Passo gaussiano a = "<<passo<<endl;
	algoritmo2.SetPasso(passo);
  algoritmo2.Equilibrate(N_equilibration,"gauss");
  algoritmo2.Run("gauss");
  algoritmo2.PrintAccepted();
	r = algoritmo2.Raggi();
	TLC(100, r, "../Risultati/Raggi_p_gauss.out");


 return 0;
}
