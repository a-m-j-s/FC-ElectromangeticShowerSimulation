#ifndef __Funcoes__
#define __Funcoes__

#include"Particle.h"
#include"Material.h"
#include<vector>

class Funcoes{
public:
	Funcoes(){;}
	Funcoes(vector<Particle>,vector<Particle>,vector<Particle>);
	~Funcoes(){;}
  	void Desenhar();
 	void Hist(double);
 	void Ficheiros();

private:

 vector<Particle> vPH;
 vector<Particle> vPO;
 vector<Particle> vEL;
 
					   
};
#endif
