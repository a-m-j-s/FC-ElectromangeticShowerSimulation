#include"Propagator.h"
#include"Material.h"
#include"Particle.h"
#include"Integrator.h"

#include"cFCgraphics.h"
#include"TF1.h"
#include"TRandom.h"
#include"TVector3.h"

#include<cmath>
#include<vector>
#include<iostream>
using namespace std;

///constructor/////////////////////////////////////////////////////////////////////////////////////
Propagator::Propagator(Particle &p,material m){
  P=p;
  Al=m;
  flag=1;
}

///Get Particle////////////////////////////////////////////////////////////////////////////////////
Particle& Propagator::GetParticle(){
  return P;
}

///Get Flag///////////////////////////////////////////////////////////////////////////////////////
int Propagator::Getflag(){
  return flag;
}

///ParR///////////////////////////////////////////////////////////////////////////////////////////
double Propagator::ParR(){
  delete gRandom;
  gRandom = new TRandom(0);
  
  double alpha=1/137.035999;
  double Na=6.02214085774e23;
  double re=2.81794e-13;
  double Xo=Al.A/(4*alpha*Na*Al.Z*Al.Z*re*re*log(183*pow(Al.Z, 1/3)));
  double ppar=(Al.ro/Xo)*7/9;
  
  //cout<<ppar<<endl;
  
  TF1 *rpar = new TF1("rpar", "-log(x)/[0]", 0 , 1);
  
  rpar->SetParameter(0, ppar);
  
  double posi=0;
  
  double aleat=gRandom->Uniform();
  posi=rpar->Eval(aleat);
  
  return posi;
}

///ParE////////////////////////////////////////////////////////////////////////////////////////////
void Propagator:: ParE(double Eg, double& Eele, double& Epos){

  delete gRandom;
  gRandom = new TRandom(0);
  
  double me=0.5109989461; //Mev/(c*c)
  int c=1;
  
  double Emin=me*c*c/Eg;
  double Emax=1-Emin;
  double k=1/Emin;
  double ct=44.503;
  
  double alpha=1/137.035999;
  double a=alpha*Al.Z;
  double Fo=(-0.1774-12.10*a+11.18*a*a)*pow((2/k), 0.5) + (8.523+73.26*a-44.41*a*a)*(2/k) - (13.52+121.1*a-96.41*a*a)*pow((2/k),3/2) + (8.946+ 62.05*a -63.41*a*a)*(2/k)*(2/k);
  double Fc=a*a*(1/(1+a*a)+0.202059-0.03693*a*a+0.00835*a*a*a*a-0.00201*a*a*a*a*a*a+0.00049*a*a*a*a*a*a*a*a-0.00012*a*a*a*a*a*a*a*a*a*a+0.00003*a*a*a*a*a*a*a*a*a*a*a*a);
  double gk=4*log(ct)-4*Fc+Fo;
  
  TF1* ener = new TF1("ener","(2*(0.5-x)^2*(7./3.-2*log(1+([0]/(x*(1-x)))^2)-6*[0]/(x*(1-x))*atan(x*(1-x)/[0])-([0]/(x*(1-x)))^2*(4-4*[0]/(x*(1-x))*atan(x*(1-x)/[0])-3*log(1+([0]/(x*(1-x)))^(-2)))+[1])+11./6.-2*log(1+([0]/(x*(1-x)))^2)-3*[0]/(x*(1-x))*atan(x*(1-x)/[0])+0.5*([0]/(x*(1-x)))^2*(4-4*[0]/(x*(1-x))*atan(x*(1-x)/[0])-3*log(1+x*x*(1-x)*(1-x)/([0]*[0])))+[1])/[2]",Emin,Emax);
  ener->SetParameter(0,ct/(2*k));
  ener->SetParameter(1, gk);
  ener->SetParameter(2, 1);
  double max=ener->GetMaximum(Emin, Emax);
  //cout<<"max: "<<max<<endl;
  ener->SetParameter(2, max);
  
  int teste=1;
  double ep=0;
  double R2=0;
  
  while (teste){
    ep=(Emax-Emin)*gRandom->Uniform()-Emin;
    R2=gRandom->Uniform();
    if (R2<ener->Eval(ep)){teste=0;}
  }
  
  Eele = ep*Eg-me*c*c;
  Epos = Eg- Eele- 2*me*c*c;
  //cout<<"Eele: "<<Eele<<" Epos: "<<Epos<<endl;	
  
}

///Par Ang/////////////////////////////////////////////////////////////////////////////////////////
void Propagator::ParAng(double Epart, double& anga, double& angb){

  delete gRandom;
  gRandom = new TRandom(0);
 
 double me=0.5109989461; //Mev/(c*c)
 int c=1;
 
 double beta=sqrt(Epart*(Epart+2*me*c*c))/(Epart+me*c*c);
 
 TF1* ang = new TF1("ang","((1-[0]*cos(x))^(-2)*sin(x))/[1]",0, M_PI);
 ang->SetParameter(0, beta);
 ang->SetParameter(1, 1);
 double maxa=ang->GetMaximum(0, M_PI);
 ang->SetParameter(1, maxa);
 
 int teste=1;
 double R2;
 
 while (teste){
   anga=M_PI*gRandom->Uniform();
   R2=gRandom->Uniform();
   if (R2<ang->Eval(anga)){teste=0;}
 }
 
 angb=2*M_PI*gRandom->Uniform()-M_PI;
 //angb=2*M_PI*gRandom->Uniform();
 
}

///BremR///////////////////////////////////////////////////////////////////////////////////////////
double Propagator::BremR(double E){

  delete gRandom;
  gRandom = new TRandom(0);
  
  double re=2.81794e-13;
  double alpha=1/137.035999;
  double me=0.5109989461; //Mev/(c*c)
  int c=1;
  double Na=6.02214085774e23;
  
  TF1* rBrem = new TF1("rBrem","([0]/x*((x/[1])^2*(4*log([2])+2-log(1+([2]*[3]*(x/[1])/(1-x/[1]))^2)-4*([2]*[3]*(x/[1])/(1-x/[1]))*atan(([2]*[3]*(x/[1])/(1-x/[1]))^(-1)))+4/3*(1-x/[1])*(4*log([2])+7/2-2*log(1+([2]*[3]*(x/[1])/(1-x/[1]))^2)-6*([2]*[3]*(x/[1])/(1-x/[1]))*atan(([2]*[3]*(x/[1])/(1-x/[1]))^(-1))-([2]*[3]*(x/[1])/(1-x/[1]))^2*(4-4*([2]*[3]*(x/[1])/(1-x/[1]))*atan(([2]*[3]*(x/[1])/(1-x/[1]))^(-1))-3*log(1+([2]*[3]*(x/[1])/(1-x/[1]))^(-2))))))/[4]",50E-6,E);
  rBrem->SetParameter(0, re*re*alpha*Al.Z*Al.Z);
  rBrem->SetParameter(1, E+me*c*c);
  rBrem->SetParameter(2, 44.503);
  rBrem->SetParameter(3, 1/(1+(E/(me*c*c)))/2);
  rBrem->SetParameter(4, 1);
  
  double iinicial=50E-6;
  double ifinal=E; 
  double sigma=0, erro=0;
  
  Integrator fsigma(iinicial, ifinal, rBrem);
  fsigma.SimpsonRule(10000, sigma, erro);
  //cout<<"sigma: "<<sigma<<" e: "<<erro<<endl;
  
  double psigma=Na*Al.ro/Al.A*sigma;
  
  //cout<<psigma<<endl;
  
  TF1 *rsigma = new TF1("rsigma", "-log(x)/[0]", 0 , 1);
  
  rsigma->SetParameter(0, psigma);
  
  double posi=0;
 
  double aleat=gRandom->Uniform();
  posi=rsigma->Eval(aleat);
  
 return posi;
}

///BremE///////////////////////////////////////////////////////////////////////////////////////////
void Propagator::BremE(double E, double& Efotao, double& Enovo){

  delete gRandom;
  gRandom = new TRandom(0);
  
  double re=2.81794e-13;
  double alpha=1/137.035999;
  double me=0.5109989461; //Mev/(c*c)
  double Na=6.02214085774e23;
  int c=1;
  
  TF1* rBrem = new TF1("rBrem","([0]/x*((x/[1])^2*(4*log([2])+2-log(1+([2]*[3]*(x/[1])/(1-x/[1]))^2)-4*([2]*[3]*(x/[1])/(1-x/[1]))*atan(([2]*[3]*(x/[1])/(1-x/[1]))^(-1)))+4/3*(1-x/[1])*(4*log([2])+7/2-2*log(1+([2]*[3]*(x/[1])/(1-x/[1]))^2)-6*([2]*[3]*(x/[1])/(1-x/[1]))*atan(([2]*[3]*(x/[1])/(1-x/[1]))^(-1))-([2]*[3]*(x/[1])/(1-x/[1]))^2*(4-4*([2]*[3]*(x/[1])/(1-x/[1]))*atan(([2]*[3]*(x/[1])/(1-x/[1]))^(-1))-3*log(1+([2]*[3]*(x/[1])/(1-x/[1]))^(-2))))))/[4]",50E-6,E);
  rBrem->SetParameter(0, re*re*alpha*Al.Z*Al.Z);
  rBrem->SetParameter(1, E+me*c*c);
  rBrem->SetParameter(2, 44.503);
  rBrem->SetParameter(3, 1/(1+(E/(me*c*c)))/2);
  rBrem->SetParameter(4, 1);
  
  double iinicial=50E-6;
  double ifinal=E; 
  double max=rBrem->GetMaximum(iinicial, ifinal);
  //cout<<"max: "<<max<<endl;
  rBrem->SetParameter(4, max);

  double teste=1;
  double R2=0;

  while (teste){
    Efotao=(ifinal-iinicial)*gRandom->Uniform()+iinicial;//-+
    R2=gRandom->Uniform();	
    if (R2<rBrem->Eval(Efotao)){teste=0;}
  }

  Enovo=E-Efotao;
  //cout<<"Eele: "<<Efotao<<" Enovo: "<<Enovo<<endl;	
}

///BremAng////////////////////////////////////////////////////////////////////////////////////////
void Propagator::BremAng(double E, double& angf, double& azi){

  delete gRandom;
  gRandom = new TRandom(0);
  
  double me=0.5109989461; //Mev/(c*c)
  int c=1;
  
  double betapart=sqrt(E*(E+2*me*c*c))/(E+me*c*c);
  
  TF1* ang = new TF1("ang","((3/(16*[0])*(1+((cos(x)-[1])/(1-[1]*cos(x)))^2)*([2])^2*(1-[1]*cos(x))^(-2))*sin(x))/[3]",0, M_PI);
  ang->SetParameter(0, M_PI);
  ang->SetParameter(1, betapart);
  ang->SetParameter(2, me*c*c/(E+me*c*c));
  ang->SetParameter(3, 1);
  double maxa=ang->GetMaximum(0, M_PI);
  //ang->SetParameter(3, maxa);
  
  int teste=1;
  double R2;

  while (teste){
    angf=M_PI*gRandom->Uniform();
    R2=gRandom->Uniform();
    if (R2<ang->Eval(angf)){teste=0;}
  }
  //azi=2*M_PI*gRandom->Uniform();
  azi=M_PI*gRandom->Uniform()-M_PI;
}

//AnPR////////////////////////////////////////////////////////////////////////////////////////////
double Propagator::AnPR(){

  delete gRandom;
  gRandom=new TRandom(0);
  
  double re=2.8179403227e-13;
  double NA=6.022140858e23;
  double gamma=1+P.E()/P.m();
  double sigma=((gamma*gamma+4*gamma+1)*log(gamma+pow(gamma*gamma-1,0.5))-(3+gamma)*pow(gamma*gamma-1,0.5))*M_PI*re*re/((gamma+1)*(gamma*gamma-1));
  double prob=NA*Al.Z*sigma*Al.ro/(Al.A);
  
  double rint=-log(1-gRandom->Uniform())/prob;
  return rint;
}

//AnPu/////////////////////////////////////////////////////////////////////////////////////////////
double Propagator::AnPu(){
  delete gRandom;
  gRandom=new TRandom(0);
  
  double re=2.8179403227e-13;
  double gamma=1+P.E()/P.m();
  double umax=0.5;
  double umin=1/(1+gamma+pow(gamma*gamma-1,0.5));
  TF1 *f=new TF1("funcao","[0]*(2*[1]+[2]*(1/x)-(1/(x*x))+[2]*(1/(1-x))-(1/((1-x)*(1-x))))");
  double p0=M_PI*re*re/((gamma+1)*(gamma*gamma-1));
  double p1=-(gamma+1)*(gamma+1);
  double p2=gamma*gamma+4*gamma+1;
  f->SetParameters(p0,p1,p2);
  double fmax=f->GetMaximum(umin,umax);
  double u=0;
  int teste=0,a=0;
  while(teste!=1 && a<100000){
    double er=umin+(umax-umin)*gRandom->Uniform();
    double R=gRandom->Uniform();
    double val=(f->Eval(er))/fmax;
    if(R<val){
      u=er;
      teste=1;
    }
    ++a;
  }
  return u;
}

///AnPE////////////////////////////////////////////////////////////////////////////////////////////
void Propagator::AnPE(double u,double &E1, double &E2){
  E1=u*(P.E()+2*P.m());
  E2=P.E()+2*P.m()-E1;
}

//AnPang///////////////////////////////////////////////////////////////////////////////////////////
void Propagator::AnPang(double u, double &teta1,double &teta2,double &fi1,double &fi2){
  delete gRandom;
  gRandom=new TRandom(0);
  
  double gamma=1+P.E()/P.m();
  teta1=acos(pow(gamma*gamma-1,-0.5)*(gamma+1-1/u));
  teta2=acos(pow(gamma*gamma-1,-0.5)*(gamma+1-1/(1-u)));
  fi1=2*M_PI*gRandom->Uniform()-M_PI;
  //fi1=2*M_PI*gRandom->Uniform();
  fi2=fi1+M_PI;
}

//PE///////////////////////////////////////////////////////////////////////////////////////////////
void Propagator::PE(double l){
  double beta=sqrt(P.E()*(P.E()+2*P.m()))/(P.E()+P.m());
  double gamma=1+P.E()/P.m();
  double B=log10(beta*gamma);
  double f=0;
  if(P.Name()=="Electron"){
    f=1-beta*beta-(2*gamma-1)*log(2)/(gamma*gamma)+(gamma-1)*(gamma-1)/(8*gamma*gamma);
  }
  else if(P.Name()=="Positron"){
    f=2*log(2)-beta*beta*(23+14/(gamma+1)+10/((gamma+1)*(gamma+1))+4/((gamma+1)*(gamma+1)*(gamma+1)))/12;
  }
  double I=0;
  if(Al.Z<13){
    I=12*Al.Z+7;
  }
  else{
    I=9.76*Al.Z +58.8*pow(Al.Z,-1.19)*Al.Z;
  }
  double delta=0;
  if(B<Al.U0){
    delta=0;
  }
  else if(B>Al.U0 && B<Al.U1){
    delta=4.6052*B+Al.C0+Al.a*pow(Al.U1-B,Al.m);
  }
  else if(B>Al.U1){
    delta=4.6052*B+Al.C0;
  }

  double C=1e-6*I*I*((0.422377/(beta*beta*gamma*gamma)+0.0304043/(beta*beta*beta*beta*gamma*gamma*gamma*gamma)-0.00038106/(beta*beta*beta*beta*beta*beta*gamma*gamma*gamma*gamma*gamma*gamma))+(3.850190/(beta*beta*gamma*gamma)+0.1667989/(beta*beta*beta*beta*gamma*gamma*gamma*gamma)-0.00157955/(beta*beta*beta*beta*beta*beta*gamma*gamma*gamma*gamma*gamma*gamma))*I*1e-3);

  double dE=l*0.1532*Al.ro*Al.Z*(log(P.E()*1e12*P.E()*(gamma+1)/(I*I*2))+f-delta-2*C/Al.Z)/(Al.A*beta*beta);
  P.E()=P.E()-dE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//PropagatePH////////////////////////////////////////////////////////////////////////////////////
vector<Particle> Propagator::PropagatePH(){
  
  cout<<"---PH---"<<endl;
  
  Particle pos,ele;
  double me=0.5109989461; //Mev/(c*c)
  //R
  double R;
  R=ParR();
  //E
  double Eg=P.E();
  double Eele;
  double Epos;
  ParE(Eg, Eele, Epos);
  //Angpos
  double angpos;
  double azipos;
  ParAng(Epos, angpos, azipos);
  //Angele
  double angele;
  double aziele;
  ParAng(Eele, angele, aziele);
  //
  TVector3 Mi=P.p(); 
  Mi.SetMag(1);
  TVector3 Pi=P.Xi();
  TVector3 Pf=Pi+R*Mi;
  //pos
  pos.SetEnergy(Epos);
  pos.SetInitialP(Pf);
  double mmp=sqrt((Epos+me)*(Epos+me)-me*me);

  TVector3 V1, V2;
  V1.SetMagThetaPhi(1,angpos,azipos);
  V2=MC(Mi, V1);
  V2.SetMag(mmp);
  pos.SetMomentum(V2);
  pos.SetName("Positron");
  
  //ele
  ele.SetEnergy(Eele);
  ele.SetInitialP(Pf);
  
  double mme=sqrt((Eele+me)*(Eele+me)-me*me);
  
  TVector3 V3, V4;
  V3.SetMagThetaPhi(1,angele,aziele);
  V4=MC(Mi, V3);
  V4.SetMag(mme);
  ele.SetMomentum(V4);
  ele.SetName("Electron");


  flag=0;

  //alterar particula original
  P.SetFinalP(Pf);

  //vetor
  vector<Particle> V;
  V.push_back(ele);
  V.push_back(pos);
  
  /*
  cout<<"Propagate PH"<<endl;
  P.Print();
  pos.Print();
  ele.Print();
  */
  
  return V;
}

//PropagateEL//////////////////////////////////////////////////////////////////////////////////////
vector<Particle> Propagator::PropagateEL(){

  cout<<"---EL---"<<endl;

  vector<Particle> V;
  double r=0;

  while(P.E()>50e-6){//50e-6

    Particle ph;
    
    //mais de um brem
    double rintRB=BremR(P.E());
    r+=rintRB;
    PE(rintRB);
    
    if(P.E()<0){
      P.E()=0;
      //cout<<"OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"<<endl;
      break;
    }
  
    double En=0,Ep=0,teta=0,fi=0;
    // cout<<"Olá1"<<endl;
    BremE(P.E(),Ep,En);
    BremAng(P.E(),teta,fi);
    // cout<<"Olá2"<<endl;
    P.E()=En;
  
    double me=0.5109989461; //Mev/(c*c)
    double mme=sqrt((Ep+me)*(Ep+me)-me*me);
    TVector3 V3, V4;
    V3.SetMagThetaPhi(1,teta,fi);
    V4=MC(P.p(),V3);
    V4.SetMag(mme);
    
    if(Ep>5){
      double aux1=r/P.p().Mag();
      TVector3 xinicial=P.Xi()+P.p()*aux1;
      
      ph.SetMomentum(V4);
      ph.SetEnergy(Ep);
      ph.SetInitialP(xinicial);
      ph.SetName("Photon");
      
      V.push_back(ph);
      //ph.Print();
    }
    //cout<<"E="<<P.E()<<endl;
  }
  double aux=r/(P.p().Mag());
  TVector3 xfinal=P.Xi()+P.p()*aux;
  P.SetFinalP(xfinal);
  //cout<<"Ola2"<<endl;
  
  /* cout<<"Propagate El"<<endl;
  P.Print();
  for(int i=0;i<V.size();++i){
    V[i].Print();
    }*/
  
  flag=0;
  return V;
   
}

//PropagatePO/////////////////////////////////////////////////////////////////////////////////////
vector<Particle> Propagator::PropagatePO(){
  
  cout<<"---PO---"<<endl;
  
  double r=0;
  vector<Particle> V;

  while(P.E()>50e-6){//50e-6
  	double rintRB=BremR(P.E());
    	// cout<<rintRB<<endl;
    	double rintAP=AnPR();
         // cout<<rintAP<<endl;
    	if(rintRB<rintAP){
    		r+=rintRB;
      		Particle ph;
      		PE(rintRB);
      		if(P.E()<0){
			P.E()=0;
			//cout<<"Olá"<<endl;
			break;
      		}

      		// cout<<"Ola1"<<endl;
      		double En=0,Ep=0,teta=0,fi=0;
      		BremE(P.E(),Ep,En);
      		BremAng(P.E(),teta,fi);
      		// cout<<"ola2"<<endl;
      		P.E()=En;
  
	      	double me=0.5109989461; //Mev/(c*c)
	      	double mme=sqrt((Ep+me)*(Ep+me)-me*me);
	      	TVector3 V3, V4;
	      	V3.SetMagThetaPhi(1,teta,fi);
	      	V4=MC(P.p(),V3);
	      	V4.SetMag(mme);
      		
      		if(Ep>5){
			double aux1=r/P.p().Mag();
			TVector3 xinicial=P.Xi()+P.p()*aux1;
	
			ph.SetMomentum(V4);
			ph.SetEnergy(Ep);
			ph.SetInitialP(xinicial);
			ph.SetName("Photon");
	
			V.push_back(ph);
			//cout<<"po brem"<<endl;
       			//ph.Print();
      		}
      		//cout<<"E="<<P.E()<<endl;
    	}
    	else{
    		//cout<<"QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQq"<<endl;
      		Particle p1,p2;
      		double l=rintAP;//-r+rintRB;
      		r+=rintAP;
      		PE(l);
      		if(P.E()<0){
			P.E()=0;
			//cout<<"Olá"<<endl;
			break;
      		}
     	 	double u=AnPu();
      		double E1=0,E2=0,teta1=0,teta2=0,fi1=0,fi2=0;
      		AnPE(u,E1, E2);
      		AnPang(u,teta1,teta2,fi1,fi2);
      		TVector3 V1, V2;
      		TVector3 V3, V4;
      
      		double me=0.5109989461; //Mev/(c*c)
      		double mm1=sqrt((E1+me)*(E1+me)-me*me);
      		double mm2=sqrt((E2+me)*(E2+me)-me*me);

      		double aux1=r/P.p().Mag();
      		TVector3 xinicial=P.Xi()+P.p()*aux1;
      
      		V1.SetMagThetaPhi(1,teta1,fi1);
		V2=MC(P.p(), V1);
	     	V2.SetMag(mm1);
	      	p1.SetMomentum(V2);
	      
	      	p1.SetEnergy(E1);
	      	p1.SetInitialP(xinicial);
	      	p1.SetName("Photon");
	      
	      	V3.SetMagThetaPhi(1,teta2,fi2);
	      	V4=MC(P.p(), V3);
	 	V4.SetMag(mm2);
	      	p2.SetMomentum(V4);
	      
	      	p2.SetEnergy(E2);
	      	p2.SetInitialP(xinicial);
	      	p2.SetName("Photon");
	      	//cout<<"Propagate Po an"<<endl;
	      	//P.Print();
	      	//p1.Print();
	      	// p2.Print();
	      	V.push_back(p1);
	      	V.push_back(p2);
       		double aux=r/(P.p().Mag());
       		TVector3 xfinal=P.Xi()+P.p()*aux;
       		P.SetFinalP(xfinal);
       
       		flag=0;
       		return V;
  
    	}
  }
  double aux=r/(P.p().Mag());
  TVector3 xfinal=P.Xi()+P.p()*aux;
  P.SetFinalP(xfinal);

  //cout<<"Propagate po"<<endl;
  /* P.Print();
  for(int i=0;i<V.size();++i){
    V[i].Print();
    }*/

  flag=0;
  
  //cout<<"Ola2"<<endl;

  return V;
}
      


//////////////////////////////////////////////////////////////////////

TVector3 Propagator::MC(TVector3 S, TVector3 PS){
  
double teta=S.Theta();
double phi=S.Phi();
  
vector<TVector3> M;
TVector3 V1(cos(teta)*cos(phi),-sin(phi),sin(teta)*cos(phi));
TVector3 V2(cos(teta)*sin(phi),cos(phi),sin(teta)*sin(phi));
TVector3 V3(-sin(teta),0,cos(teta));
M.push_back(V1);
M.push_back(V2);
M.push_back(V3);

double *Paux=new double[3];
for(int i=0;i<3;++i){
	Paux[i]=M[i][0]*PS[0]+M[i][1]*PS[1]+M[i][2]*PS[2];
}
TVector3 PL(Paux[0],Paux[1],Paux[2]);
//cout<<"V2(Lab)=("<<Pi[0]<<","<<Pi[1]<<","<<Pi[2]<<");"<<endl;
delete []Paux;

return PL;
}


