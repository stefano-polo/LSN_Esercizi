#ifndef _Funzioni_h_
#define _Funzioni_h_

#include <cmath>



class FunzioneBase {
	public:
		virtual double Eval(double x) const=0;

};



class Coseno: public FunzioneBase {

	virtual double Eval (double x) const;


};

class G: public FunzioneBase{
	
	virtual double Eval (double x) const;
};




#endif //_Funzioni_h_

