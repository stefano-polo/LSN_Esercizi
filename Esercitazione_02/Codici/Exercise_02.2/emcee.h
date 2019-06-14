#ifndef _emcee_h_
#define _emcee_h_


#include "random.h"


class Walker:public Random {
private:
	Random _rnd;
	double _x;
	double _y;
	double _z;
	double _a;
	
	


	public:
	Walker(double a);
	~Walker();	
	void Passo();
	void Restart(){_x=0;_y=0;_z=0;}
	double RndWalk(int N);
	double Diffusione(int N);
	void Isos();
	
	
	
};


#endif
