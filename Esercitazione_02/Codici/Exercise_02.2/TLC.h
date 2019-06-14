#ifndef _TLC_h_
#define _TLC_h_

#include <vector>

using namespace std;

double Sum(vector<double> v, int i, int j); //somma gli elementi di v dall'i-esimo incluso, al j-esimo escluso
double Media(vector<double> v);
double Varianza(vector<double> v);
vector<double>  TLC(int N, vector <double> v); //fa datablocking restituendo valore medio e errore con 100 blocchi
void Print(vector<vector<double>> matrix, const char * filename); //stampa su file i risultati del random walk con una matrice contenente valore medio e errore

#endif
