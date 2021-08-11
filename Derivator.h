#ifndef __Derivator__
#define __Derivator__
#include"Func1D.h"
#include"TF1.h"

class Derivator: public Func1D {
public:
 Derivator(TF1*fp=NULL) :Func1D(fp) {;}
 ~Derivator(){;}
 void FDerivative(double,double,double&);
 void SDerivative(double, double, double&);
 void FourthDerivative(double, double, double&); 
 void IDerivative(double,double&,double*,double*,int);
 TF1* GetfSDerivative(double, double);
 TF1* GetfFDerivative(double,double);
 TF1* GetfFourthDerivative(double,double);
 private:
 double fSDerivative(double *fx,double *par);
 double fFDerivative(double *fx,double *par);
 double fFourthDerivative(double *fx,double *par);
 
					   
};
#endif
