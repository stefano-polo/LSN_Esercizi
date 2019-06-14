#ifndef __fluid_
#define __fluid_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

const int wd=12; //per lo spazio tra i dati nelle tabelle

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig,iu2;
double nbins;
double walker[m_props];
double Kb = 8.62e-5;


// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g,stima_u2;
double err_u,err_c,err_m,err_x,err_g;
double energy, energy2;
//configuration
const int m_spin=50;
double s[m_spin];
int mode; //se è 0 parto da spin random mentre se è 1 da una vecchia configurazione su file
// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk, metro;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);
void Metropolis(int o);
void Gibbs(int o);
void Energy(int);
void Energy2(int);
void Chi(int);
void Heath(int);
void Mag(int);

#endif
