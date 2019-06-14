#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <fstream>
#include "TLC.h"

using namespace std;


double Somma(vector<double> v, int i, int j){

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


double Error(vector<double> AV, vector<double> AV2,int n, int i){
	if (n==0){
        	return 0;
	}
	else{
        	return sqrt((AV2[n] - pow(AV[n],2))/i);
	}
}
void TLC(int N, int M, vector<double> s, const char * namefile) {
// riempio il vettore mean con le medie ottenuti da blocchi di M/N numeri e std con le medie al quadrato

	int L = int(M/N);
	vector<double> ai(N);
	vector<double> ai_2(N);
 	double a; // variabile di appoggio

  	for(int i=0; i<N;i++){
		a= Somma(s,i*L,(i+1)*L)/L;
		ai.at(i) = a;
		ai_2.at(i) = pow(a,2);
	}

//calcolo la media progressiva del vettore mean e le deviazione standard della media progressiva
 	vector<double> mean_prog(N);
	vector<double> mean_2_prog(N);
	vector<double> mean_std(N);

	for(int i=0;i<N;i++){
		mean_prog.at(i) = Somma(ai,0,i+1)/(i+1);
		mean_2_prog.at(i) = Somma(ai_2,0,i+1)/(i+1);
		mean_std.at(i) = sqrt( ( mean_2_prog.at(i) - pow(mean_prog.at(i),2) ) / (i+1)  );
	}

// Stampo su file i risultati nel seguente formato: numero di dati usati per il calcolo di media e deviazione standard della media, media progressiva, deviazione standard della media progressiva
	ofstream Risultati;
	Risultati.open(namefile);

	for(int i=1;i<=N;i++){
		if (Risultati.is_open()){
   			Risultati <<i*L<< " "<< mean_prog.at(i-1) << " " <<mean_std.at(i-1)<< endl;
   		}
		else cerr << "PROBLEM: Unable to open the file.out" << endl;

	}

Risultati.close();
}
