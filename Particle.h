#ifndef __Particle__
#define __Particle__

#include"TVector3.h"
#include<string>
using namespace std;

class Particle{
 public:
  Particle();
  Particle(string,double,TVector3,TVector3);
  Particle(const Particle&);
  const Particle operator=(const Particle&);
  ~Particle(){;}
  
  void SetEnergy(double);
  void SetVivo(int);
  void SetMomentum(TVector3);
  void SetInitialP(TVector3);
  void SetFinalP(TVector3);
  void SetName(string);

  double GetEnergy();
  int GetVivo();
  TVector3 GetMomentum();
  TVector3 GetInitialP();
  TVector3 GetFinalP();
  string GetName();

  double& E();
  TVector3& p();
  TVector3& Xi();
  TVector3& Xf();
  string& Name();

  void Print();
  double m();

 protected:
  double En;//energia
  TVector3 mo;//momento linear
  TVector3 xi;//posicao inicial
  TVector3 xf;//posicao final
  string N;
  int vivo;
};
#endif
