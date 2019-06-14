#ifndef __MolDyn_NVE__
#define __MolDyn_NVE__


#include "random.h"
//parameters, observables
const int m_props=3000;
int n_props, iv, iw, it,ik,itot,igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];
double T_target;

//generator
Random rnd;
int seed[4];

//numero di blocchi
int nblk;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];
double xold[m_part],yold[m_part],zold[m_part];
double xold2[m_part],yold2[m_part],zold2[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double t; //variabile per il calcolo della velocità quadratica media del sistema

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,err_pot,err_press,stima_gdir,err_gdir,stima_temp,temperatura, err_temp,stima_ekin,err_ekin,stima_etot,err_etot;
// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

//pigreco
const double pi=3.1415927;

// simulation
int nstep, iprint;
double delta;

//functions
void Input(int);
void Move(void);
void Restart(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Termalizzazione(void);
void Square_Velocity(void);
void Velocity(void);
void Kinetic_Energy(void);
void Total_Energy(void);
void Temperature(void);
void Print_Istant(void);
//blocking
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);



#endif // __MolDyn_NVE__
