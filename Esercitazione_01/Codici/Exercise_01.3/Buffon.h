#ifndef _Buffon_h_
#define _Buffon_h_


#include "random.h"


using namespace std;

class Buffon:public Random {
private:
	Random _rnd;
	double _L;  //lunghezza ago
	double _d;  //distanza tra le righe
	
	
	public:
		
	Buffon(double L, double d);
	~Buffon();
	double Experiment(int N_throws);
	double Sin_Angle();
};


#endif










