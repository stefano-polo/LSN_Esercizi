#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){

	Input(); //Inizialization
	for(int istep=1; istep <= nstep/1e2; ++istep){ //fase di equilibrazione
			Move(metro);

	}
	cout<<endl<<endl<<"Ising 1D for h = "<<h<<endl<<endl;
	do{
		beta = 1./temp;   //aggiornamento del parametro beta
		for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
			accepted = 0;
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep){
				Move(metro);
				Measure();
				Accumulate(); //Update block averages
			}
			if(iblk ==nblk){
			cout << "Block number " << iblk << endl;
			cout<<"Temperature " << temp <<endl;
		//	cout << "Acceptance rate " << accepted/attempted << endl << endl;  //per il Gibbs non ha senso impostare un acceptance rate dato che accetta sempre
			}
			Averages(iblk);   //Print results for current block--->calcola tutte le grandezze medie eccetto la magnetizzaione
			if(iblk ==nblk){
			cout << "----------------------------" << endl << endl;}
		}

		temp += 0.05;  //update della temperatura
	}
	while(temp <= 2.1);//2.1

	//Input(); //inizializzo nuovamente il sistema per simulare la magnetizzazione a campo h = 0.02
	temp = 0.5; //Inizializzo nuovamente la temperatura
	h = 0.02; //Accendo il campo
	beta = 1./temp;   //aggiornamento del parametro beta
	for(int istep=1; istep <= nstep/1e2; ++istep){  //fase di equilibrazione
				Move(metro);

			}
	cout<<endl<<endl<<"Ising 1D for h = "<<h<<endl<<endl;

	while(temp <= 2.1){
		beta = 1./temp;   //aggiornamento del parametro beta
		for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep){
				Move(metro);
				Measure();
				Accumulate(); //Update block averages
			}
			if(iblk ==nblk){
			cout << "Block number " << iblk << endl;
			cout<<"Temperature " << temp <<endl;
		//  cout << "Acceptance rate " << accepted/attempted << endl << endl;
			}
			Mag(iblk);   //Stampo risultati magnetizzazione media
			if(iblk ==nblk){
			cout << "----------------------------" << endl << endl;}
		}
		temp += 0.05;  //update della temperatura

	}

	ConfFinal(); //Write final configuration*/
	return 0;
}


void Input(void){

  ifstream ReadInput;

  cout << "Classic 1D Ising model         " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
	ReadInput >> mode;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  iu2 = 4;

  n_props = 5; //Number of observables

	if(mode == 1) {  //lettura configurazione spin da un file scelto
		cout << "Read old configuration from file old.0 " << endl << endl;
		ifstream ReadConf;
		ReadConf.open("old.0");

		for (int i=0; i<nspin; ++i){
			ReadConf >> s[i] ;

		}
		ReadConf.close();

	}
	else if(mode == 0){
//initial configuration                                   //Configurazione iniziale a T=infinito, cioè spin random
	cout<<"Initial configuration with T = infinity"<<endl;
	for (int i=0; i<nspin; ++i){

    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
}
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro){
  int o;

  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis

      Metropolis(o);
    }
    else{ //Gibbs sampling

      Gibbs(o);
    }
  }
}

double Boltzmann(int sm, int ip){
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure(){

  double u = 0.0, s_tot = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i){
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);   //Energia interna (H)
    s_tot += s[i];

  }

  walker[iu] = u;
	walker[iu2] = pow(u,2);
  walker[im] = s_tot;  //Magnetizzazione
  walker[ix] = pow(s_tot,2);    //Suscettività

}


void Reset(int iblk){ //Reset block averages

  if(iblk == 1){
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i){
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


void Accumulate(void){ //Update block averages

  for(int i=0; i<n_props; ++i){
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Energy(int iblk) {
	ofstream Ene, ResultsU;

	Ene.open("../../Risultati/Gibbs/output.ene.0",ios::app);
  	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  	glob_av[iu]  += stima_u;
  	glob_av2[iu] += stima_u*stima_u;
  	err_u=Error(glob_av[iu],glob_av2[iu],iblk);
	energy = glob_av[iu]/(double)iblk;
  	Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  	Ene.close();
	if(iblk == nblk){
   	 	ResultsU.open("../../Risultati/Gibbs/ene.result.0",ios::app);  //Stampo i valori piu precisi per il grafico in funzione della temperatura
    	ResultsU << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
			ResultsU.close();
	}

}

void Energy2(int iblk) {
	stima_u2 = blk_av[iu2]/blk_norm/(double)nspin; //Energy
  glob_av[iu2]  += stima_u2;
	energy2 = glob_av[iu2]/(double)iblk;

}

void Heath(int iblk) {

	ofstream Heat, ResultsC;

	Heat.open("../../Risultati/Gibbs/output.heat.0",ios::app);
  	//stima_u2 = blk_av[ic]/blk_norm/(double)nspin; // <H^2>
  	stima_c = pow(beta,2) * (energy2*(double)nspin-pow(energy,2)*(double)nspin*nspin)/(double)nspin; //Colore specifico
  	glob_av[ic]  += stima_c;
  	glob_av2[ic] += stima_c*stima_c;
  	err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  	Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  	Heat.close();
	if(iblk == nblk){
   	 	ResultsC.open("../../Risultati/Gibbs/heat.result.0",ios::app);  //Stampo i valori piu precisi per il grafico in funzione della temperatura
    		ResultsC << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
		ResultsC.close();
	}

}

void Mag(int iblk){
	ofstream Mag, ResultsM;

	Mag.open("../../Risultati/Gibbs/output.mag.0",ios::app);
	stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetizzazione
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);
	Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
	Mag.close();
	if(iblk == nblk){
   	 	 ResultsM.open("../../Risultati/Gibbs/mag.result.0",ios::app);  //Stampo i valori piu precisi per il grafico in funzione della temperatura
		 ResultsM << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
		 ResultsM.close();
	}


}

void Chi(int iblk){
	ofstream Chi, ResultsX;

	Chi.open("../../Risultati/Gibbs/output.chi.0",ios::app);
  	stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Suscettività
  	glob_av[ix]  += stima_x;
  	glob_av2[ix] += stima_x*stima_x;
  	err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  	Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  	Chi.close();
	if(iblk == nblk){
   	 	 ResultsX.open("../../Risultati/Gibbs/chi.result.0",ios::app); //Stampo i valori piu precisi per il grafico in funzione della temperatura
    		 ResultsX << setw(wd) << temp <<  setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
		 ResultsX.close();
	}
}

void Averages(int iblk){ //Print results for current block
	Energy(iblk);
	Energy2(iblk);
	Heath(iblk);
	Chi(iblk);
}


void ConfFinal(void){
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i){  //Algorithm for periodic boundary conditions

  if(i >= nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk){
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void Metropolis(int o){ //Algoritmo Metropolis
  double deltaE = 2.* J * s[o] * ( s[Pbc(o-1)] + s[Pbc(o+1)] ) +h*s[o]*2;

  if(deltaE < 0) s[o] = s[o]*(-1), accepted++;
  else{
    double a = rnd.Rannyu();
    double boltzweight = exp(-1/temp * deltaE);
	if( a <= boltzweight){
		s[o] = s[o]*(-1);
		accepted++;
         }
  }
	attempted++;
}


void Gibbs(int o){  //Algoritmo di Gibbs

 double deltaE = 2.* J * beta * ( s[Pbc(o-1)] + s[Pbc(o+1)] ) + 2 *beta*h;
 double p_up = 1./(1.+exp(-deltaE));  //probabilita spin up
 double a = rnd.Rannyu();
 if (a < p_up) s[o] = 1;
 else s[o] = -1;

}
