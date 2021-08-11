#include"Func1D.h"
#include"TF1.h"
#include"cFCgraphics.h"
#include<iostream>
#include <stdio.h>
using namespace std;

//#define DEBUG

Func1D::Func1D(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}


Func1D::Func1D(TF1* fp){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

  p=new TF1(*fp);

}

Func1D::~Func1D(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  delete p;
}

void Func1D::Draw(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  cFCgraphics G;
  TPad* pad1=G.CreatePad("pad1");
  G.AddObject(p,"pad1");
  G.AddObject(pad1);
  G.Draw();
  getchar();
}

double Func1D::Evaluate(double x){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  double f=p->Eval(x);
  return f;
}
  
