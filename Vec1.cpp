#include "Vec1.h"
#include <cmath>
#include <stdio.h>
#include <iostream>
using namespace std;

//#define DEBUG

//////////////////////////////////////////////////////////////////
///////////Constructors and Destructor///////////////////////////
/////////////////////////////////////////////////////////////////

Vec::Vec(){
	entries= new double[1];
	entries[0]=0;
	N=1;
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	}

Vec::Vec(int n): N(n){ 
	entries = new double[n];
        for (int i=0; i<n; i++){
	entries[i]=0;
	}
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif	
}

Vec::Vec(int n, double ini): N(n){ 
	entries = new double[n];
        for (int i=0; i<n; i++){
	entries[i]=ini;
	}
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}

Vec::Vec(int n, double* ini): N(n){ 
	entries = new double[n];
        for (int i=0; i<n; i++){
	entries[i]=ini[i];
	}
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}

Vec::Vec(const Vec& v){ 
        N=v.N;
	entries = new double[N];
        for (int i=0; i<N ;i++)
        entries[i] = v.entries[i]; 
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}


Vec::~Vec(){ 
	delete[] entries;
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}
////////////////////////////////////////////////////////////////
//////Operators
///////////////////////////////////////////////////////////////
const Vec& Vec::operator=(const Vec& v){
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif

  if(this!=&v)
    {
      N=v.N;
      delete[] entries;
      entries=new double[v.N];
      for(int i=0; i<N; i++)
	entries[i]=v.entries[i];
    }
return *this;
}

//-----------------------------------------------
Vec Vec:: operator+(const Vec& v)const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if(N!=v.N){
    cout<<"Erro-tamanhos diferentes"<<endl;
    Vec novo(N);
    return novo;}
  else{
    Vec novo(N);
    for(int i=0;i<N;++i){
      novo.entries[i]=entries[i]+v.entries[i];
    }
    return novo;
    }
}
//-------------------------------------------------
Vec Vec:: operator+=(const Vec& v){
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if(N!=v.N){
    cout<<"Erro-tamanhos diferentes"<<endl;
    return *this;}
  else{
    for(int i=0;i<N;++i){
      entries[i]=entries[i]+v.entries[i];
    }
    return *this;
    }
}
//------------------------------------------------
Vec Vec:: operator-(const Vec& v)const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if(N!=v.N){
    cout<<"Erro-tamanhos diferentes"<<endl;
    Vec novo(N);
    return novo;}
  else{
    Vec novo(N);
    for(int i=0;i<N;++i){
      novo.entries[i]=entries[i]-v.entries[i];
    }
    return novo;
    }
}
//-------------------------------------------------
Vec Vec:: operator-=(const Vec& v){
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if(N!=v.N){
    cout<<"Erro-tamanhos diferentes"<<endl;
    return *this;}
  else{
    for(int i=0;i<N;++i){
      entries[i]=entries[i]-v.entries[i];
    }
    return *this;
    }
}
//------------------------------------------------
double& Vec:: operator[](int i){
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
return entries[i];
}
//------------------------------------------------
Vec Vec:: operator-()const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
    Vec novo(N);
    for(int i=0;i<N;++i){
      novo.entries[i]=-entries[i];}
    return novo;
}
//------------------------------------------------
Vec Vec:: operator+()const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
    Vec novo(N);
    for(int i=0;i<N;++i){
      novo.entries[i]=abs(entries[i]);}
    return novo;
}

//-------------------------------------------------
Vec Vec:: operator*(const Vec& v)const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if(N!=v.N){
    cout<<"Erro-tamanhos diferentes"<<endl;
    Vec novo(N);
    return novo;}
  else{
    Vec novo(N);
    for(int i=0;i<N;++i){
      novo.entries[i]=entries[i]*v.entries[i];
    }
    return novo;
    }
}
//-------------------------------------------------
Vec Vec::operator*(double escalar)const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
Vec novo(N);
for(int i=0;i<N;++i){
    novo.entries[i]=escalar*entries[i];
    }
return novo;


}
 
////////////////////////////////////////////////////////////////
//////Metodos
///////////////////////////////////////////////////////////////
void Vec::Print(){
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
for(int i=0; i<N;i++)
cout<<entries[i]<<"  ";
cout<<endl;}


void Vec::SetEntries(int n, double* a){ 
        N=n;
	delete[] entries;
	entries = new double[N];
        for (int i=0; i<n ;i++)
        entries[i] = a[i]; 
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}

double Vec::At(int n)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	double valor=0;
	if (n<N){valor=entries[n];}
	else cout<<"erro-dimensao"<<endl;
	return valor;
}


int Vec::size()const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
return N;
}


double Vec::dot(const Vec& a)const{
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
if (N != a.N){
	printf("erro-produto interno-dimensao vetores\n");
	return 999;}
double soma=0;
for (int i=0; i<N; i++)
	soma+=entries[i]*a.entries[i];
return soma;
}


void Vec::swap(int a, int b){
#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
if(a>N-1 || b>N-1 ||a<0||b<0){cout<<"Erro-swap nao feito-valor inexistente"<<endl;}
  else{
  double temp;
  temp=entries[a];
  entries[a]=entries[b];
  entries[b]=temp;
  }
}


