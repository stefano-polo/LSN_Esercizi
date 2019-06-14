#ifndef _TLC_h_
#define _TLC_h_

#include <vector>
#include "random.h"
using namespace std;


double Error(vector<double> AV, vector<double> AV2, int n);
void TLC(int N, vector<double> s, const char * filename); //calcolo data blocking con N blocchi e scrittura su file
void TLC_modified(int N, vector<double> s, const char * filename); //per calcolare il secondo integrale e eseguire direttamente data blocking
void Chi2(int M, int N_tot, vector<double> s, const char * filename);

#endif //_TLC_h_
