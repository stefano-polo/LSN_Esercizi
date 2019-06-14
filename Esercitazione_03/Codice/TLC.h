#ifndef _TLC_h_
#define _TLC_h_

#include <vector>

using namespace std;




double Somma(vector<double> v, int i, int j);
double Media(vector<double> v);
double Varianza(vector<double> v);
double Error(vector<double> AV, vector<double> AV2, int n,int i);
void TLC(int N, int M, vector<double> s, const char * filename);

#endif //_TLC_h_
