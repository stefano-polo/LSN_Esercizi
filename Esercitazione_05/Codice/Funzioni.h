#ifndef _Funzioni_h_
#define _Funzioni_h_

#include <cmath>
#include <vector>

using namespace std;

class FunzioneBase {
	public:
		virtual double Eval(vector<double> v) const=0;

};



class S_uno: public FunzioneBase {

	virtual double Eval (vector<double> v) const;


};

class P_due: public FunzioneBase{
	
	virtual double Eval (vector<double> v) const;
};




#endif //_Funzioni_h_

