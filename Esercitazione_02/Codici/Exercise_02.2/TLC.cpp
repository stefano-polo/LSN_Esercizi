#include <vector>
#include <cmath>
#include<iostream>
#include<fstream>
#include<string>

#include "TLC.h"

using namespace std;

double Sum(vector<double> v, int i, int j){

	double somma=0;
	for(int a=i; a<j; a++){
		somma += v.at(a);
	}
	return somma;
}

double Media(vector<double> v) {
	double somma=0;
	for(unsigned int i=0; i<v.size(); i++) {
		somma+=v.at(i);
	}

	return somma/v.size();

}

double Varianza(vector<double> v) {
	double somma=0;
	for(unsigned int i=0; i<v.size(); i++){
		somma+=pow(Media(v)-v.at(i), 2.);
	}

	return somma/(v.size()-1) ;
 }


vector<double> TLC(int N, vector<double> v){
	unsigned int M = v.size();
  int L = int(M/N);   //Numero di esperimenti per blocco
  vector<double> vmean(N,0);   //Vettore che mi contiene le medie dei singoli blocchi
  vector<double> vmean2(N,0);   //Vettore che contiene i quadrati delle medie dei singoli blocchi
  vector<double> mean_cum(N,0);  //Vettore che mi contiene la media cumulativa
  vector<double> mean2_cum(N,0);   //Vettore che contiene i quadrati delle medie cumulativi
  vector<double> err_cum(N,0);      //vettore che mi contiene  gli errori

  for(int i=0; i<N; i++){
 	 double temp = Sum(v,i*L,(i+1)*L)/L;
 	 vmean.at(i) = temp;
 	 vmean2.at(i) = pow(temp , 2);
 	 }

  for(int i=0; i<N; i++){
 	 mean_cum.at(i) = Sum(vmean, 0, i+1)/(i+1);
 	 mean2_cum.at(i) = Sum(vmean2, 0, i+1)/(i+1);
 	 err_cum.at(i) = sqrt( ( mean2_cum.at(i) - pow(mean_cum.at(i),2) ) / (i+1)  );
 	 }

	vector<double> output(2); //vettore contenente media e deviazione standard
	output[0] = sqrt(mean_cum.at(N-1));
	output[1] = 0.5*(1/sqrt(output[1]))*err_cum.at(N-1);  //propagazione degli errori
	return output;
}

void Print(vector<vector<double>> matrix, const char * filename) {
	ofstream Risultati;
	Risultati.open(filename);
	for(unsigned int i =0; i<matrix.size();i++){
		Risultati<<i<<" "<<matrix[i][0]<<" "<<matrix[i][1]<<endl;
	}

	Risultati.close();
}
