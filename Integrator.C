#include"Func1D.h"
#include"Derivator.h"
#include"Integrator.h"
#include"TF1.h"
#include<iostream>
#include<cmath>
#include<stdio.h>

using namespace std;

//#define DEBUG


Integrator::Integrator(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
x0=0;
x1=0;

}

void Integrator::TrapezoidalRule1(int n,double &Integral,double &Error){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  double h=(x1-x0)/n;
  double*X=new double[n+1];
  X[0]=x0;
  X[n]=x1;
  for(int i=1;i<n;++i){X[i]=X[i-1]+h;}
  double*Y=new double[n+1];
  for(int i=0;i<n+1;++i){Y[i]=p->Eval(X[i]);}
 
  Integral=0;
  Error=0;
  Derivator D(p);
  TF1*Fd=D.GetfSDerivative(x0,x1);
  double d=Fd->GetMaximum(x0,x1);
  //cout<<d<<endl;
  Error=(x1-x0)*(x1-x0)*(x1-x0)*d/(12*n*n);
  //cout<<"ola"<<endl;
  for(int k=0;k<n;++k){
    Integral=Integral+0.5*h*(Y[k]+Y[k+1]);
    // cout<<Integral<<" "<<Error<<endl;
  }
}

void Integrator::TrapezoidalRule2(int n,double &Integral,double &Error){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  for(int k=0;k<n;++k){
    double temp=0;
    double pow2=pow(2,k-2);
    double pow1=pow(2,k-1);
    for(int i=1;i<pow2+1;++i){
      temp=temp+p->Eval(x0+(2*i-1)*(x1-x0)/pow1);
    }
    Error=Integral;
    Integral=Integral/2+(x1-x0)*temp/pow1;
    Error=Error-Integral;
  }
}

void Integrator::SimpsonRule(int n,double &Integral, double &Error){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  if(n%2!=0){n=n-1;}//n tem de ser numero par
    double h=(x1-x0)/n;
    double*X=new double[n+1];
    X[0]=x0;
    X[n]=x1;
    for(int i=1;i<n;++i){X[i]=X[i-1]+h;}
    double*Y=new double[n+1];
    for(int i=0;i<n+1;++i){Y[i]=p->Eval(X[i]);}
    Integral=Y[0]+Y[n];
    Derivator D(p);
    TF1*Fd=D.GetfFourthDerivative(x0,x1);
    double d=Fd->GetMaximum(x0,x1);
    //cout <<d<<endl;
    Error=(x1-x0)*(x1-x0)*(x1-x0)*(x1-x0)*(x1-x0)*d/(180*n*n*n*n);
    for(int i=1;i<n;++i){
      if(i%2==0){
	Integral=Integral+2*Y[i];
      }
      else{
	Integral=Integral+4*Y[i];
      }
    }  
    Integral=Integral*h/3;
}
