#ifndef __Integrator__
#define __Integrator__
#include"Func1D.h"

class Integrator: public Func1D {
public:

 Integrator();
 Integrator(double fx0, double fx1, TF1*fp) :x0(fx0), x1(fx1), Func1D(fp) {;}
 ~Integrator(){;}
 void TrapezoidalRule1(int n, double& Integral, double& Error);
 void TrapezoidalRule2(int n, double& Integral, double& Error);
 void SimpsonRule(int n, double& Integral, double& Error);
 
 protected:
 double x0;
 double x1;
};
#endif
