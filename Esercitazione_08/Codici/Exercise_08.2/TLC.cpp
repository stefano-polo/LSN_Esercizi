#include "TLC.h"
#include "random.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double Sum(vector<double> v, int start, int end){

	double sum = 0;
	for(int i=start; i<end; i++) sum += v.at(i);
	return sum;
}

double Mean(vector<double> v) {
	double sum = 0;
	for(unsigned int i=0; i<v.size(); i++) sum += v.at(i);
	return sum/v.size();
}

double Error(vector<double> v) {
	double sum = 0;
	for(unsigned int i=0; i<v.size(); i++) sum += pow(Mean(v)-v.at(i), 2.);
	return sum/(v.size()-1);
}

void TLC(int N, vector<double> v, const char*namefile){         // Mi permette di trovare le medie e le deviazioni standard della media attrsaverso il metodo dei blocchi
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

  ofstream results;
  results.open(namefile);

  for(int i=0; i<N; i++){

   if(results.is_open()) results << i << " " << mean_cum.at(i) << " " << err_cum.at(i) << endl;
   else cerr << "PROBLEM: Unable to open result.txt" << endl;
  }

  results.close();

}

void AnalisiIstogramma(int N, int Nbins, double binsize, double a, vector<vector<double>> f, const char* namefile) {
	int M = f.size();
	int L = int(M/N);
	vector<double> appo(M);
	ofstream results;
	results.open(namefile);

	for(int r=0; r<Nbins; r++){
		for(int i=0; i<M; i++){
			appo[i] = f[i][r];
		}

		vector<double> vmean(N,0);   //Vettore che mi contiene le medie dei singoli blocchi
  		vector<double> vmean2(N,0);   //Vettore che contiene i quadrati delle medie dei singoli blocchi
  		vector<double> mean_cum(N,0);  //Vettore che mi contiene la media cumulativa
  		vector<double> mean2_cum(N,0);   //Vettore che contiene i quadrati delle medie cumulativi
  		vector<double> err_cum(N,0);      //vettore che mi contiene  gli errori

		for(int i=0; i<N; i++){
	  		double temp = Sum(appo,i*L,(i+1)*L)/L;
    			vmean.at(i) = temp;
    			vmean2.at(i) = pow(temp , 2);
	 	}

		for(int i=0; i<N; i++){
    			mean_cum.at(i) = Sum(vmean, 0, i+1)/(i+1);
    			mean2_cum.at(i) = Sum(vmean2, 0, i+1)/(i+1);
    			err_cum.at(i) = sqrt( ( mean2_cum.at(i) - pow(mean_cum.at(i),2) ) / (i+1)  );
    		}

		if(results.is_open()) results << a+(r)*binsize << " " << mean_cum.at(N-1) << " " << err_cum.at(N-1) << endl;
   		else cerr << "PROBLEM: Unable to open result.txt" << endl;


	}

	results.close();
}
