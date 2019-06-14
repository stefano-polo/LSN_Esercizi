#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <ostream>
#include <iomanip>
#include <vector>
#include "MolDyn_NVE.h"
#include "random.h"

using namespace std;

int main(){
  Input(0);   //Inizialization

  int n_termalizzazione = 4000;

  ofstream Term;
  Term.open("../../../Risultati/NVE_Verlet/Liquid/output_restart.dat",ios::app);

  cout<<endl<<endl<<"Thermalization of the system"<<endl;
   //termalizzazione del sistema con metodo di partire da una temperatura più alta (con il metodo Termalizzazione())
   for(int i = 0; i<n_termalizzazione;i++){
      Move();

      if(i%10 == 0){
         Termalizzazione();
         Term << stima_temp<< endl;

     }
  }
  Term.close();
  cout << "Temperature after thermalization: " << stima_temp << endl;

  cout << endl<<endl;

  cout<<"Starting simulation"<<endl;
  int nconf = 1;
  // evoluzione del sistema all'equilibrio
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
     cout <<endl<< "Block number " << iblk << endl;
     Reset(iblk);   //Reset block averages
     for(int istep=1; istep <= nstep; ++istep){
        Move();           //Move particles with Verlet algorithm
        if(istep%10 == 0){
          Measure();   //Properties measurement
          Accumulate(); //Update block averages
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
          nconf += 1;
        }
    }
    Averages(iblk);   //Print results for current block
   }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(int mode){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
 // double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

  //Lettura dei valori di input
  ReadInput.open("input.liquid"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblk;

  cout <<" Number of blocks for data blocking = "<<nblk<<endl;

  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();
  T_target = temp;
//Prepare array for measurements
 //Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
  it = 2;
  //commento gli indici per il calcolo dell'energia cinetica e energia totale poichè non sono richieste dall'esercizio
  //ik = 3;
  //itot = 4;
  n_props = 3; //Number of observables

//measurement of g(r)
  igofr = 3;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;


  //If you want the old configuration too, read it
  if(mode == 1) {
    cout << "Read old configuration from file old.0 " << endl << endl;
    ReadConf.open("old.0");

    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }

    ReadConf.close();


  }
  else{
  //Read initial configuration fcc
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  }

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu() - 0.5;
     vy[i] = rnd.Rannyu() - 0.5;
     vz[i] = rnd.Rannyu() - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = x[i] - vx[i] * delta;
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta;


   }


   return;
}

void Restart(void) {
	// r(t-dt)= xold2
        // r(t) =xold[i]
	//r(t+dt) = x[i]

   double sumv2, fs;
   Velocity();
   Square_Velocity();
   sumv2 = (t/0.5)/(double)npart;

   stima_temp = (2.0 / 3.0) * t/(double)npart;
     //double compare_T = abs(T_target-stima_temp);
     fs = sqrt(3 * T_target / sumv2);   // fs = velocity scale factor
     for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = x[i] - vx[i] * delta;   //ho messo io PBC e il 2
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta;
  }

     return;


}


void Termalizzazione(){     //processo che termalizza il sistema
  double sumv[3] = {0.0, 0.0, 0.0};

   for(int i=0; i<npart; ++i){ //Verlet integration scheme

     vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);
     vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
     vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }

   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   Square_Velocity();
   temp = (2.0 / 3.0) * t/(double)npart;
   stima_temp=temp;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;
   }



   return;

}

void Velocity(void){
    for(int i=0; i<npart; ++i){
    vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);    //r(t+dt)-r(t-dt)----> guarda la funzione Move();
    vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
    vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
}
  return;
}

void Square_Velocity(void) {    //funzione utile per il calcolo della temperatura e dell'energia cineticas
  t=0.;
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    //r(t-dt) lo salvo per il calcolo delle velocità nella funzione Velocity
    xold2[i] = xold[i];
    yold2[i] = yold[i];
    zold2[i] = zold[i];

    //r(t+dt)
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );


    //r(t)
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;   //r(t+dt)
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure()
{

  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     walker[igofr+int(dr/bin_size)]+=2; //ho dr/binsize mi identifica in quale bin mi trovo

     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;
  Velocity();
  Square_Velocity();
  Temperature();
//  Kinetic_Energy(); non mi interessa il calcolo dell'energia cinetica per questo esercizio
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();

  return;
}

void Temperature(void) {
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  walker[it] = stima_temp;
  return;
}

void Kinetic_Energy(void) {
  walker[ik] = t/(double)npart;

}



void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
  rnd.SaveSeed();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void Reset(int iblk) //Reset block averages
{

   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{

   double r,deltaV;
   ofstream Gofr, Gave, Epot, Pres,Temp,Ekin,Etot;
   const int wd=12;
    Epot.open("../../../Risultati/NVE_Verlet/Liquid/output.epot.0",ios::app);
    Pres.open("../../../Risultati/NVE_Verlet/Liquid/output.pres.0",ios::app);
    Gofr.open("../../../Risultati/NVE_Verlet/Liquid/output.gofr.0",ios::app);
    Gave.open("../../../Risultati/NVE_Verlet/Liquid/output.gave.0",ios::app);
    Temp.open("../../../Risultati/NVE_Verlet/Liquid/output.temp.0",ios::app);
    //commento il calcolo dell'energia cinetica e energia totale (medie) poichè non sono richieste dall'esercizio
  //  Ekin.open("../../../Risultati/NVE_Verlet/Solid/ave_ekin.out",ios::app);
  //  Etot.open("../../../Risultati/NVE_Verlet/Solid/ave_etot.out",ios::app);


    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_pres = rho * stima_temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);

    temperatura = blk_av[it]/blk_norm; //Temperatura
    glob_av[it] += temperatura;
    glob_av2[it] += temperatura*temperatura;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);


    /*stima_ekin = blk_av[ik]/blk_norm; //Energia Cinetica
    glob_av[ik] += stima_ekin;
    glob_av2[ik] += stima_ekin*stima_ekin;
    err_ekin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = stima_ekin+stima_pot; //Energia totale
    glob_av[itot] += stima_etot;
    glob_av2[itot] += stima_etot*stima_etot;
    err_etot=Error(glob_av[itot],glob_av2[itot],iblk);*/

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Pressure
    Pres << setw(wd) << iblk <<  setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;
//Average Temperature
   Temp << setw(wd) << iblk <<  setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Average kinetic energy
//   Ekin << setw(wd) << iblk <<  setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_ekin << endl;
//Average total energy
  // Etot << setw(wd) << iblk <<  setw(wd) << glob_av[itot]/(double)iblk << setw(wd) << err_etot << endl;


//g(r)   Radial function distribution

  for (int k=igofr; k<igofr+nbins; ++k) {
	   r = (k-1)*bin_size; //faccio variare r
	   deltaV = (4./3.)*M_PI*(pow(r+bin_size,3)-pow(r,3));  //calcolo il delta(r)
	   stima_gdir = (1./(rho*deltaV*(double)npart))*(blk_av[k]/blk_norm);
	   glob_av[k] += stima_gdir;
     glob_av2[k] += stima_gdir*stima_gdir;
     err_gdir=Error(glob_av[k],glob_av2[k],iblk);
	   Gofr   <<  setw(wd) << r << setw(wd) << stima_gdir<<endl;  //stampo i valori di g(r) di ciascun blocco
	   if(iblk == nblk){  //prendo il valore del blocco finale
		     Gave << setw(wd) << r <<  setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err_gdir<< endl; //stampo l'ultimo blocco di g(r) al variare del raggio con rispettivo errore
    }

  }

    cout <<endl<< "----------------------------" << endl << endl;

    Epot.close();
    Pres.close();
    Gofr.close();
    Gave.close();
    Temp.close();
  //  Ekin.close();
    //Etot.close();
}
