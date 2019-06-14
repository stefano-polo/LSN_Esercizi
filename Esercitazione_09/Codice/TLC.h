#ifndef _TLC_h_
#define _TLC_h_

#include "random.h"
#include <vector>
using namespace std;

void Start(Random *rnd);
double Sum(vector<double> v, int start, int end);
double Mean(vector<double> v);
double Error(vector<double> v);
void TLC(int N, vector<double> v, const char*namefile);
  
    // N = numero di blocchi in cui divido le M simulazioni
    // v = vettore in cui ci saranno i numeri casuali
    // namefile = nome del file in cui scriverà i risultati
void AnalysisMatrix(int N, vector<vector <double>> v, double scale_factor, double index_traslation, const char* namefile); //data blocking di una matrice stampando media e deviazione standard della media realizzata con tutti i blocchi



#endif
