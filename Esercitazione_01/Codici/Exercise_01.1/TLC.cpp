#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <fstream>
#include "TLC.h"
#include "random.h"
using namespace std;


double Error(vector<double> AV, vector<double> AV2, int n,int i){
	if (n==0){
        	return 0;
	}
	else{
        	return sqrt((AV2[n] - pow(AV[n],2))/i);
	}
}
void TLC(int N, vector<double> s, const char * filename) {
	ofstream Risultati;
	Risultati.open(filename);
	unsigned int M = s.size();
	int L = int(M/N);
	vector<double> ave(N);
	vector<double> av2(N);
	vector<double> sum_prog(N);
	vector<double> su2_prog(N);
	vector<double> err_prog(N);

	double sum;
	int k = 0;
	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			k = j+i*L;
			sum+= s[k];
		}
		ave[i] = sum/L;
		av2[i] = pow(ave[i],2);
	}


	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);
		su2_prog[i] /= (i+1);
		err_prog[i] = Error(sum_prog,su2_prog,i,i+1);
	}

	for(int i=0;i<N;i++){


	if (Risultati.is_open()){
   		Risultati <<i*L<< " "<< sum_prog[i] << " " <<err_prog[i]<< endl;
   }
	else cerr << "PROBLEM: Unable to open random.out" << endl;

   }
	Risultati.close();


}


void TLC_modified(int N, vector<double> s, const char * filename) {  //metodo per il calcolo del secondo integrale richiesto dall'esercizio
	ofstream Risultati;
	Risultati.open(filename);
	unsigned int M = s.size();
	int L = int(M/N);
	vector<double> ave(N);
	vector<double> av2(N);
	vector<double> sum_prog(N);
	vector<double> su2_prog(N);
	vector<double> err_prog(N);

	double sum;
	int k = 0;
	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			k = j+i*L;
			sum+= pow(-0.5+s[k],2);
		}
		ave[i] = sum/L;
		av2[i] = pow(ave[i],2);
	}


	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i] /= (i+1);
		su2_prog[i] /= (i+1);
		err_prog[i] = Error(sum_prog,su2_prog,i,i+1);
	}

	for(int i=0;i<N;i++){


	if (Risultati.is_open()){
   		Risultati <<i*L<< " "<< sum_prog[i] << " " <<err_prog[i]<< endl;
   }
	else cerr << "PROBLEM: Unable to open random.out" << endl;

   }
	Risultati.close();


}


void Chi2(int M, int N_tot, vector<double> s, const char * filename){
	//numero di dati generati casualmente sui cui calcolo i chi quadrati singoli
	int n = int(N_tot/M);

	//vettore contenente i chi quadrati singoli
	vector<double> chi(M);

	//vettore che contiene 0 e l'estremo di ogni intervallo in cui è diviso [0,1)
	vector<double> intervalli(M+1);
	for(int i=0; i<M+1; i++){
		intervalli.at(i) = (1.*i)/M;
	}


	//il vettore n_i serve a contenere il numero di punti che cadadono in ognuno degli M intervalli.
 	vector<double> n_i(M);

	//Segue il calcolo dei chi quadrati singoli
	for(int k=0; k<M; k++){

		//Ad ogni iterazione n_i viene riempito di zeri per iniziare il conteggio del numero di punti in ogni intervallo
		for(int b=0; b<M; b++){
			n_i.at(b) = 0.;
		}

		//calcolo del numero di punti che cadono in ognuno degli M intervalli
		  for(int i=k*n; i<(k+1)*n; i++){
			for(int j=0; j<M+1; j++){
		    		if( (s.at(i) >= intervalli.at(j) ) and ( s.at(i) < intervalli.at(j+1) ) ){
					n_i.at(j) = n_i.at(j)+1  ;
				}
			}
		  }



		//calcolo del chi quadrato di ognuno dei blocchi da n dati
		double somma= 0;
		double attesi = ((1.*n)/(1.*M));
		for(int a=0; a<M; a++){
			somma +=  pow( n_i.at(a)- attesi ,2) / attesi;
		}
		chi.at(k) = somma;

	}
//stampo su file i risultati
	ofstream Risultati;
	Risultati.open(filename);

	for(unsigned int i=0; i<chi.size(); i++){
		Risultati<< i << " " << chi.at(i) << endl;
	}



}
