#include "Statistica.h"
#include <vector>
#include <cmath>
#include "random.h"
#include <iostream>
#include <fstream>

using namespace std;


double Media(vector<double> v) {	
	double somma=0;
	for(unsigned int i=0; i<v.size(); i++) {
		somma+=v.at(i);
	}

	return somma/v.size();

}

double Varianza(vector<double> v) {
	double sum=0;
	for(unsigned int i=0; i<v.size(); i++){
		sum+=pow(Media(v)-v.at(i), 2.);
	}

	return sum/(v.size()-1) ;
 }
