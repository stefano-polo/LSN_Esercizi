#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <fstream>
#include "TLC.h"

using namespace std;



double Error(vector<double> AV, vector<double> AV2,int n, int i){
	if (n==0){
        	return 0;
	}
	else{
        	return sqrt((AV2[n] - pow(AV[n],2))/i);
	}
}
void TLC(int N,  vector<double> s, const char * filename) {
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
