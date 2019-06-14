#ifndef _TLC_h_
#define _TLC_h_

#include "random.h"
#include <vector>
using namespace std;

//void Start(Random *rnd);
double Sum(vector<double> v, int start, int end);
double Mean(vector<double> v);
double Error(vector<double> v);
void TLC(int N, vector<double> v, const char*namefile);  //esegue data blocking su un vettore di elementi v
void AnalisiIstogramma(int N, int Nbins, double binsize, double a, vector<vector<double>> v, const char* namefile);
//esegue data blocking per istogramma (una matrice) di Nbins con dimensione binsize, con estremo sinistro a

  // N = numero di blocchi in cui divido le v.size() simulazioni
    // v = vettore in cui ci saranno i numeri casuali
    // namefile = nome del file in cui scriverà i risultati


#endif
