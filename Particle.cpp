#include"Particle.h"

#include"TVector3.h"
#include<stdio.h>
#include<string>
#include<iostream>
using namespace std;

//Default Constructor////////////////////////////////////////////////////////////////////
Particle::Particle(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
  N="nada";
  vivo=1;
#endif
}

//Constructor///////////////////////////////////////////////////////////////////////////
Particle::Particle(string n,double energia, TVector3 mom, TVector3 inicio): N(n), En(energia),mo(mom), xi(inicio){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
  vivo=1;
#endif
}

//Set Energy///////////////////////////////////////////////////////////////////////////
void Particle::SetEnergy(double a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  En=a;
}

//Set Vivo///////////////////////////////////////////////////////////////////////////
void Particle::SetVivo(int a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  vivo=a;
}

//Set Linear Momentum///////////////////////////////////////////////////////////////////////////
void Particle::SetMomentum(TVector3 a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  mo=a;
}

//Set Initial Position///////////////////////////////////////////////////////////////////////////
void Particle::SetInitialP(TVector3 a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  xi=a;
}

//Set Final Position///////////////////////////////////////////////////////////////////////////
void Particle::SetFinalP(TVector3 a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  xf=a;
  vivo=0;
}

//Set Name///////////////////////////////////////////////////////////////////////////
void Particle::SetName(string a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  N=a;
}

//Get Vivo///////////////////////////////////////////////////////////////////////////
int Particle::GetVivo(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return vivo;
}

//Get Energy///////////////////////////////////////////////////////////////////////////
double Particle::GetEnergy(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return En;
}

//Get Linear Momentum///////////////////////////////////////////////////////////////////////////
TVector3 Particle::GetMomentum(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return mo;
}

//Get Initial Position///////////////////////////////////////////////////////////////////////////
TVector3 Particle::GetInitialP(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return xi;
}

//Get Final Position///////////////////////////////////////////////////////////////////////////
TVector3 Particle::GetFinalP(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return xf;
}

//Get Name///////////////////////////////////////////////////////////////////////////
string Particle::GetName(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return N;
}

//Energy///////////////////////////////////////////////////////////////////////////
double& Particle::E(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return En;
}

//Get Linear Momentum///////////////////////////////////////////////////////////////////////////
TVector3& Particle::p(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return mo;
}

//Initial Position///////////////////////////////////////////////////////////////////////////
TVector3& Particle::Xi(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return xi;
}

//Final Position///////////////////////////////////////////////////////////////////////////
TVector3& Particle::Xf(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return xf;
}

//N///////////////////////////////////////////////////////////////////////////
string& Particle::Name(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  return N;
}

//Print///////////////////////////////////////////////////////////////////////////////////
void Particle::Print(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  cout<<endl<<" Name="<<N<<endl;
  printf(" Energy=%lfMeV\n p=(%lf,%lf,%lf)\n xi=(%lf,%lf,%lf)\n xf=(%lf,%lf,%lf)\n\n",En,mo[0],mo[1],mo[2],xi[0],xi[1],xi[2],xf[0],xf[1],xf[2]);
}

//m//////////////////////////////////////////////////////////////////////////////////////
double Particle::m(){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif
  if(N=="Electron" || N=="Positron"){
    return 0.510998946131;
  }
  else{
    return 0;
  }
}


//copy constructor///////////////////////////////////////////////////////////////////////////
Particle::Particle(const Particle& a){
#if defined(DEBUG)
  printf("[%s]\n", __PRETTY_FUNCTION__);
#endif  
  En=a.En;
  mo=a.mo;
  xi=a.xi;
  xf=a.xf;
  N=a.N;
}

//operator =///////////////////////////////////////////////////////////////////////////////
const Particle Particle::operator=(const Particle& a){
  if(this!=&a)
    {
      En=a.En;
      mo=a.mo;
      xi=a.xi;
      xf=a.xf;
      N=a.N;
    }
  return *this;
}
    

    
  



