#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Salesman.h"
#include "TLC.h"
#include <algorithm>

using namespace std;

int main (){

    double mutation_probability=0.015; //probabilità di mutazione
    Salesman Walker(30,1000,1000);  //popolazione iniziale di 1000 e 1000 generazioni
    //Citta su circonferenza di raggio 1
    cout<<"Citta su circonferenza di raggio 1"<<endl;
    Walker.Circle();
    Walker.Run(mutation_probability);
    Walker.Best_path("../Risultati/path_circ.dat");
    Walker.Best_distance("../Risultati/distance_circ.dat");
    Walker.Print_Best_distance50("../Risultati/distance50_circ.dat");

   //Città dentro quadrato di lato 1
    Salesman Walker2(30,1000,1000);//popolazione iniziale di 1000 e 1000 generazioni
    mutation_probability=0.02; //
    cout<<"Città dentro quadrato di lato 1"<<endl;
    Walker2.Square();
    Walker2.Run(mutation_probability);
    Walker2.Best_path("../Risultati/path_sq.dat");
    Walker2.Best_distance("../Risultati/distance_sq.dat");
    Walker2.Print_Best_distance50("../Risultati/distance50_sq.dat");

  return 0;

}
