#include"Func1D.h"
#include"Derivator.h"
#include"Vec1.h"
#include"FCmatrix.h"
#include"FCmatrixFull.h"
#include"EqSolver.h"
#include"TF1.h"
#include<stdio.h>
#include<iostream>
using namespace std;

//#define DEBUG

void Derivator::FDerivative(double x,double h,double &Derivative){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  
  if(h!=0){
  	Derivative=(1/(12*h))*((p->Eval(x-2*h)+8*p->Eval(x+h))-(8*p->Eval(x-h)+p->Eval(x+2*h)));
  }
}

void Derivator::SDerivative(double x,double h,double &Derivative){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  if(h!=0){
    	Derivative=(1/(12*h*h))*((16*p->Eval(x-h)+16*p->Eval(x+h))-(p->Eval(x-2*h)+30*p->Eval(x)+p->Eval(x+2*h)));}
}

void Derivator::FourthDerivative(double x,double h,double &Derivative){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  if(h!=0){
    Derivative=((p->Eval(x+2*h)+6*p->Eval(x)+p->Eval(x-2*h))-(4*p->Eval(x+h)+4*p->Eval(x-h)))/(h*h*h*h);
  }
}

void Derivator::IDerivative(double fx, double &Derivative,double *x,double*y, int N){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  
  double **a=new double*[N-2];
  for(int i=0;i<N-2;++i){
    	a[i]=new double[N-2];
    	for(int j=0;j<N-2;++j){
      		if(i==j){
			a[i][j]=2*(x[i]-x[i+2]);
      		}
      		else if(i==j+1){
			a[i][j]=(x[i]-x[i+1]);
      		}
      		else if(j==i+1){
			a[i][j]=x[j]-x[j+1];
      		}
      		else{
			a[i][j]=0;
      		}
    	}
  }
  FCmatrixFull M(a,N-2,N-2);
  double *v=new double[N-2];
  for(int i=0;i<N-2;++i){
  	v[i]=6*((y[i]-y[i+1])/(x[i]-x[i+1])-(y[i+1]-y[i+2])/(x[i+1]-x[i+2]));
    	//cout<<v[i]<<endl;
  }
  Vec V(N-2,v);
  EqSolver E(M,V);
  Vec C=E.LUdecompositionSolver();
  double *K=new double[N];
	for(int i=1;i<N-1;++i){
    		K[i]=C[i-1];
    		//cout<<"K["<<i<<"]="<<K[i]<<endl;
    	}
  K[0]=0;
  K[N-1]=0;
  int i;
  for(i=1;i<N;++i){
  	if((fx-x[i])<0.){break;}     
  }
  Derivative=(K[i]/6)*(3*(fx-x[i+1])*(fx-x[i+1])/(x[i]-x[i+1])-(x[i]-x[i+1]))-(K[i+1]/6)*(3*(fx-x[i])*(fx-x[i])/(x[i]-x[i+1])-(x[i]-x[i+1]))+(y[i]-y[i+1])/(x[i]-x[i+1]);
  //segunda derivada=K[i]*(fx-x[i+1])/(x[i]-x[i+1])-K[i+1]*(fx-x[i])/(x[i]-x[i+1]);
}

double Derivator::fSDerivative(double *fx, double*par){
  double d;
  SDerivative(fx[0],0.00001,d);
  return d;
}

double Derivator::fFDerivative(double *fx, double*par){
  double d;
  FDerivative(fx[0],0.00001,d);
  return d;
}

double Derivator::fFourthDerivative(double *fx, double*par){
  double d;
  FourthDerivative(fx[0],0.01,d);
  return d;
}

TF1* Derivator::GetfSDerivative(double x0, double x1){
  TF1*Fd=new TF1("Fd",this, &Derivator::fSDerivative,x0,x1,0,"Derivator","SDerivative");
  return Fd;
}

TF1* Derivator::GetfFDerivative(double x0,double x1){
  TF1*Fd=new TF1("Fd",this, &Derivator::fFDerivative,x0,x1,0,"Derivator","SDerivative");
  return Fd;
}

TF1* Derivator::GetfFourthDerivative(double x0, double x1){
  TF1*Fd=new TF1("Fd",this, &Derivator::fFourthDerivative,x0,x1,0,"Derivator","SDerivative");
  return Fd;
}
