#ifndef __Propagator__
#define __Propagator__

#include"Particle.h"
#include"Material.h"
#include<vector>
using namespace std;

class Propagator{
 public:
  Propagator(){flag=1;};
  Propagator(Particle&,material);
  ~Propagator(){;};
  
  Particle& GetParticle();
  int Getflag();
  
  vector<Particle> PropagatePH();//fotao
  vector<Particle> PropagatePO();//fotao
  vector<Particle> PropagateEL();//eletrao
  
 private:
  material Al;
  Particle P;
  int flag;
  ///Criacao de pares
  double ParR ();
  void ParE(double Eg, double& Eele, double& Epos);
  void ParAng(double Epart, double& anga, double& angb);
  
  ///Brem
  double BremR(double E);
  void BremE(double E, double& Efotao, double& Enovo);
  void BremAng(double E, double& angf, double& azi);
  
  //Aniquilaçao de positrao
  double AnPR();
  double AnPu();
  void AnPE(double u, double& E1, double& E2);
  void AnPang(double u, double& teta1,double &teta2,double &fi1,double &fi2);
  
  //perda de energia
  void PE(double l);//muda a energia da particula do construtor se já estiver definida a xf
  
  //Auxi
  TVector3 MC(TVector3 S, TVector3 PS);
  
};
#endif
