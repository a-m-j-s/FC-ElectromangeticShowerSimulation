#include"Particle.h"
#include"Material.h"
#include"Propagator.h"
#include"Funcoes.h"
#include"cFCgraphics.h"
#include"TH1F.h"
#include"TLine.h"
#include"TVector3.h"
#include<iostream>
#include<string>
#include<vector>
#include<fstream>

using namespace std;

vector<Particle> vPH;
vector<Particle> vPO;
vector<Particle> vEL;

void propagate(Particle P,material Al){
  if(P.Name()=="Photon"){
    Propagator prob(P,Al);
    vector<Particle> V1=prob.PropagatePH();
    Particle Paux=prob.GetParticle();
    if(V1[0].E()>50e-6){
      propagate(V1[0],Al);
    }
    if(V1[1].E()>50e-6){
      propagate(V1[1],Al);
    }
    vPH.push_back(Paux);
  }
  else if(P.Name()=="Electron"){
    Propagator prob(P,Al);
    vector<Particle> V2=prob.PropagateEL();
    Particle Paux=prob.GetParticle();
    for(int i=0;i<V2.size();++i){
      if(V2[i].E()>5){
	propagate(V2[i],Al);
      }
    }
    vEL.push_back(Paux);
  }
  else if(P.Name()=="Positron"){
    Propagator prob(P,Al);
    vector<Particle> V3=prob.PropagatePO();
    Particle Paux=prob.GetParticle();
    for(int i=0;i<V3.size();++i){
      if(V3[i].E()>5){
	propagate(V3[i],Al);
      }
    }
    vPO.push_back(Paux);				
  }
}

int main(){
  
  material Al={13,2.7,26.9815386,-4.24,0.1708,3.01,0.0802,3.63,"Aluminio"};
  
  double E=5000;
  
  TVector3 xi(0,0,0);//posicao inicila
  TVector3 p(0,0,1);//direcao 
  double me=0.5109989461; //Mev/(c*c)
  double mm1=sqrt((E+me)*(E+me)-me*me);
  p.SetMag(mm1);
  Particle Pi("Photon",E,p,xi);

  //obter particulas
  propagate(Pi,Al);

  Funcoes More(vPH,vPO,vEL);
  More.Ficheiros();
  More.Desenhar();
  More.Hist(1000);

  return 0;
}

  
