#ifndef _Funzioni_h_
#define _Funzioni_h_

#include <cmath>
#include <vector>

using namespace std;

class FunzioneBase {
	public:
		virtual double Eval(vector<double> v) const=0;

};


class Psi: public FunzioneBase {   //è la funzione d'onda di trial dipendente da due parametri (non il suo modulo quadro)
	private:
		double _mu;
		double _sigma;
	public:
		Psi(double mu, double sigma);
		~Psi();
		virtual double Eval (vector<double> v) const;
		void SetMu(double parametro) {_mu = parametro;}//mi permette di settare i parametri mu e sigma della psi
		void SetSigma(double parametro) {_sigma = parametro;}
		double GetMu(){return _mu;}
		double GetSigma(){return _sigma;}

};

class Hamiltoniana: public FunzioneBase {   //è l'integranda ossia  (H*psi)/psi
	private:
		Psi * _psi;
	public:
		Hamiltoniana(Psi * psi){_psi = psi;}
		~Hamiltoniana();
		virtual double Eval (vector<double> v) const;

};


#endif //_Funzioni_h_
