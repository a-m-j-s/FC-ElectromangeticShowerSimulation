#include"Particle.h"
#include"Material.h"
#include"Propagator.h"
#include"Funcoes.h"
#include"cFCgraphics.h"
#include"TPolyLine3D.h"
#include"TH1F.h"
#include"TH2F.h"
#include"TLine.h"
#include"TVector3.h"
#include"TAxis3D.h"
#include<iostream>
#include<vector>
#include<fstream>
#include"Funcoes.h"

using namespace std;

Funcoes::Funcoes(vector<Particle> a,vector<Particle> b,vector<Particle> c){
  vPH=a;
  vPO=b;
  vEL=c;
}

void Funcoes::Ficheiros(){

  ofstream F;
  
  F.open("PH");
  for(int i=0; i<vPH.size();++i){
    vPH[i].Print();
    F<<" E="<<vPH[i].E()<<endl
     <<" xi=("<<vPH[i].Xi()[0]<<", "<<vPH[i].Xi()[1]<<", "<<vPH[i].Xi()[2]<<") "<<endl
     <<" xf=("<<vPH[i].Xf()[0]<<", "<<vPH[i].Xf()[1]<<", "<<vPH[i].Xf()[2]<<")"<<endl
     <<" p=("<<vPH[i].p()[0]<<","<<vPH[i].p()[1]<<","<<vPH[i].p()[2]<<")"<<endl<<endl;
      
  }
  F.close();
  
  F.open("PO");
  for(int i=0; i<vPO.size();++i){
    vPO[i].Print();
    F<<" E="<<vPO[i].E()<<endl
     <<" xi=("<<vPO[i].Xi()[0]<<", "<<vPO[i].Xi()[1]<<", "<<vPO[i].Xi()[2]<<") "<<endl
     <<" xf=("<<vPO[i].Xf()[0]<<", "<<vPO[i].Xf()[1]<<", "<<vPO[i].Xf()[2]<<")"<<endl
     <<" p=("<<vPO[i].p()[0]<<","<<vPO[i].p()[1]<<","<<vPO[i].p()[2]<<")"<<endl<<endl;
  }
  F.close();
  
  F.open("EL");
  for(int i=0; i<vEL.size();++i){
    vEL[i].Print();
    F<<" E="<<vEL[i].E()<<endl
     <<" xi=("<<vEL[i].Xi()[0]<<", "<<vEL[i].Xi()[1]<<", "<<vEL[i].Xi()[2]<<") "<<endl
     <<" xf=("<<vEL[i].Xf()[0]<<", "<<vEL[i].Xf()[1]<<", "<<vEL[i].Xf()[2]<<")"<<endl
     <<" p=("<<vEL[i].p()[0]<<","<<vEL[i].p()[1]<<","<<vEL[i].p()[2]<<")"<<endl<<endl;
  }



}


void Funcoes::Hist(double E){
  
  double zmax=0,zmin=0,xmax=0,xmin=0,ymax=0,ymin=0;

  for(int i=0;i<vPH.size();++i){
    if(vPH[i].Xf()[2]>zmax){
      zmax=vPH[i].Xf()[2];
    }
    if(vPH[i].Xi()[2]<zmin){
      zmin=vPH[i].Xi()[2];
    }
  }
  for(int i=0;i<vPO.size();++i){
    if(vPO[i].Xf()[2]>zmax){
      zmax=vPO[i].Xf()[2];
    }
    if(vPO[i].Xi()[2]<zmin){
      zmin=vPO[i].Xi()[2];
    }
  }
  for(int i=0;i<vEL.size();++i){
    if(vEL[i].Xf()[2]>zmax){
      zmax=vEL[i].Xf()[2];
    }
    if(vEL[i].Xi()[2]<zmin){
      zmin=vEL[i].Xi()[2];
    }
   
  }

  cout<<zmin<<" "<<zmax<<endl;

  cFCgraphics G2;
  TPad* pad11=G2.CreatePad("pad11");
  TPad* pad12=G2.CreatePad("pad12");
  TPad* pad13=G2.CreatePad("pad13");
  TPad* pad14=G2.CreatePad("pad14");
  
  TH1F *his11=new TH1F("his11","Histograma N Particulas",(int)(zmax-zmin),zmin,zmax);
  TH1F *his12=new TH1F("his12","Histograma N Fotoes",(int)(zmax-zmin),zmin,zmax);
  TH1F *his13=new TH1F("his13","Histograma N Eletroes",(int)(zmax-zmin),zmin,zmax);
  TH1F *his14=new TH1F("his14","Histograma N Positroes",(int)(zmax-zmin),zmin,zmax);

  for(int i=zmin;i<zmax;++i){
    for(int j=0;j<vPH.size();++j){
      double xi=vPH[j].Xi()[2];
      double xf=vPH[j].Xf()[2];
      if( (xi<i &&xf>i+1) || (xi<i+1 &&xi>i) || (xf<i+1 &&xf>i)){
	his11->AddBinContent(i);
	his12->AddBinContent(i);
      }
    }
    for(int j=0;j<vPO.size();++j){
      double xi=vPO[j].Xi()[2];
      double xf=vPO[j].Xf()[2];
      if( (xi<i &&xf>i+1) || (xi<i+1 &&xi>i) || (xf<i+1 &&xf>i)){
	his11->AddBinContent(i);
	his14->AddBinContent(i);
      }
    }
    for(int j=0;j<vEL.size();++j){
      double xi=vEL[j].Xi()[2];
      double xf=vEL[j].Xf()[2];
      if( (xi<i &&xf>i+1) || (xi<i+1 &&xi>i) || (xf<i+1 &&xf>i)){
	his11->AddBinContent(i);
	his13->AddBinContent(i);
      }
    }
  }

  int N=0;
  int Z=0;
  for(int i=zmin*1000;i<zmax*1000;++i){
    double Naux=0;
    for(int j=0;j<vPH.size();++j){
      double xi=vPH[j].Xi()[2];
      double xf=vPH[j].Xf()[2];
      if( (xi<i/1000 &&xf>i/1000+0.001) || (xi<i/1000+0.001 &&xi>i/1000) || (xf<i/1000+0.001 &&xf>i/1000)){
        Naux+=1;
      }
    }
    for(int j=0;j<vPO.size();++j){
      double xi=vPO[j].Xi()[2];
      double xf=vPO[j].Xf()[2];
      if( (xi<i/1000 &&xf>i/1000+0.001) || (xi<i/1000+0.001 &&xi>i/1000) || (xf<i/1000+0.001 &&xf>i/1000)){
        Naux+=1;
      }
    }
    for(int j=0;j<vEL.size();++j){
      double xi=vEL[j].Xi()[2];
      double xf=vEL[j].Xf()[2];
      if( (xi<i/1000 &&xf>i/1000+0.001) || (xi<i/1000+0.001 &&xi>i/1000) || (xf<i/1000+0.001 &&xf>i/1000)){
	Naux+=1;
      }
    }
    if(Naux>N){
      N=Naux;
      Z=i/1000+0.0005;
    }
  }

  double Ne=E/43;
  cout<<"---- N esperado="<<Ne<<"  N="<<N<<"  Z="<<Z<<"  -----"<<endl;


  G2.AddObject(his11,"pad11");
  G2.AddObject(his12,"pad12");
  G2.AddObject(his13,"pad13");
  G2.AddObject(his14,"pad14");
  G2.AddObject(pad11);
  G2.AddObject(pad12);
  G2.AddObject(pad13);
  G2.AddObject(pad14);
  G2.Draw();
  G2.Print("Histogramas.pdf");
}


void Funcoes::Desenhar(){
  
  double *aux=new double[6];

  for(int i=0;i<vPH.size();++i){
    for(int j=0;j<3;++j){
      if(vPH[i].Xi()[j]>aux[2*j+1]){
	aux[2*j+1]=vPH[i].Xi()[j];
      }
      if(vPH[i].Xf()[j]>aux[2*j+1]){
	aux[2*j+1]=vPH[i].Xf()[j];
      }
      if(vPH[i].Xi()[j]<aux[2*j]){
	aux[2*j]=vPH[i].Xi()[j];
      }
      if(vPH[i].Xf()[j]<aux[2*j]){
	aux[2*j]=vPH[i].Xf()[j];
      }		       
    }
  }
  for(int i=0;i<vPO.size();++i){
    for(int j=0;j<3;++j){
      if(vPO[i].Xi()[j]>aux[2*j+1]){
	aux[2*j+1]=vPO[i].Xi()[j];
      }
      if(vPO[i].Xf()[j]>aux[2*j+1]){
	aux[2*j+1]=vPO[i].Xf()[j];
      }
      if(vPO[i].Xi()[j]<aux[2*j]){
	aux[2*j]=vPO[i].Xi()[j];
      }
      if(vPO[i].Xf()[j]<aux[2*j]){
	aux[2*j]=vPO[i].Xf()[j];
      }		       
    }
  }
  for(int i=0;i<vEL.size();++i){
    for(int j=0;j<3;++j){
      if(vEL[i].Xi()[j]>aux[2*j+1]){
	aux[2*j+1]=vEL[i].Xi()[j];
      }
      if(vEL[i].Xf()[j]>aux[2*j+1]){
	aux[2*j+1]=vEL[i].Xf()[j];
      }
      if(vEL[i].Xi()[j]<aux[2*j]){
	aux[2*j]=vEL[i].Xi()[j];
      }
      if(vEL[i].Xf()[j]<aux[2*j]){
	aux[2*j]=vEL[i].Xf()[j];
      }		       
    }
  }

  cout<<aux[0]<<" "<<aux[1]<<" "<<aux[2]<<" "<<aux[3]<<" "<<aux[4]<<" "<<aux[5]<<endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  cFCgraphics G1;
  TPad *pad1=G1.CreatePad("pad1");
  pad1->Range(aux[4]-1,aux[2]-1,aux[5]+1,aux[3]+1);
  TH1F* frame1 = pad1->DrawFrame(aux[4]-1,aux[2]-1, aux[5]+1,aux[3]+1);
  frame1->GetXaxis()->SetTitle("z (cm)");
  frame1->GetYaxis()->SetTitle("y (cm)");
  frame1->SetLineColor(0);

  G1.AddObject(frame1,"pad1");
 for(int i=0;i<vPH.size();++i){
    TLine *l1=new TLine(vPH[i].Xi()[2],vPH[i].Xi()[1],vPH[i].Xf()[2],vPH[i].Xf()[1]);
    l1->SetLineColor(kRed+2);
    l1->SetLineStyle(9);
    l1->SetLineWidth(2);
    G1.AddObject(l1,"pad1");
    // G1.DrawPadFlush("pad1");
  }
  for(int i=0;i<vEL.size();++i){
    TLine *l1=new TLine(vEL[i].Xi()[2],vEL[i].Xi()[1],vEL[i].Xf()[2],vEL[i].Xf()[1]);
    l1->SetLineColor(kBlue+2);
    l1->SetLineStyle(1);
    l1->SetLineWidth(2);
    G1.AddObject(l1,"pad1");
    //G1.DrawPadFlush("pad1");
  }
  for(int i=0;i<vPO.size();++i){
    TLine *l1=new TLine(vPO[i].Xi()[2],vPO[i].Xi()[1],vPO[i].Xf()[2],vPO[i].Xf()[1]);
    l1->SetLineColor(kGreen+2);
    l1->SetLineStyle(1);
    l1->SetLineWidth(2);
    G1.AddObject(l1,"pad1");
    // G1.DrawPadFlush("pad1");
  }
  G1.AddObject(pad1);
  G1.ListObjects();
  G1.Draw();
  G1.Print("Desenhozy.pdf");
  
  G1.Clear();

  ////////////////////////////////////////////////////////////////////////

  cFCgraphics G2;
  TPad *pad2=G2.CreatePad("pad2");
  pad2->Range(aux[4]-1,aux[0]-1,aux[5]+1,aux[1]+1);
  TH1F* frame2 = pad2->DrawFrame(aux[4]-1, aux[0]-1, aux[5]+1, aux[1]+1);
  frame2->GetXaxis()->SetTitle("z (cm)");
  frame2->GetYaxis()->SetTitle("x (cm)");
  frame2->SetLineColor(0);

 
  G2.AddObject(frame2,"pad2");

  for(int i=0;i<vPH.size();++i){
     TLine *l0=new TLine(vPH[i].Xi()[2],vPH[i].Xi()[0],vPH[i].Xf()[2],vPH[i].Xf()[0]);
    l0->SetLineColor(kRed+2);
    l0->SetLineStyle(9);
    l0->SetLineWidth(2);
    G2.AddObject(l0,"pad2");
    //G1.DrawPadFlush("pad2");
  }
  for(int i=0;i<vEL.size();++i){
    TLine *l0=new TLine(vEL[i].Xi()[2],vEL[i].Xi()[0],vEL[i].Xf()[2],vEL[i].Xf()[0]);
    l0->SetLineColor(kBlue+2);
    l0->SetLineStyle(1);
    l0->SetLineWidth(2);
    G2.AddObject(l0,"pad2");
    //G1.DrawPadFlush("pad2");
  }
  for(int i=0;i<vPO.size();++i){
    TLine *l0=new TLine(vPO[i].Xi()[2],vPO[i].Xi()[0],vPO[i].Xf()[2],vPO[i].Xf()[0]);
    l0->SetLineColor(kGreen+2);
    l0->SetLineStyle(1);
    l0->SetLineWidth(2);
    G2.AddObject(l0,"pad2");
    // G1.DrawPadFlush("pad2");
  }
  G2.AddObject(pad2);
  G2.ListObjects();
  G2.Draw();
  G2.Print("Desenhoxz.pdf");
 
  G2.Clear();
  
  //////////////////////////////////////////////////////////////////////////
  
  cFCgraphics G3;
  TPad *pad3=G3.CreatePad("pad3");
  TAxis3D *Ax=new TAxis3D();
  Ax->SetAxisColor(1,"X");
  Ax->SetAxisColor(1,"Y");
  Ax->SetAxisColor(1,"Z");
  Ax->SetLabelColor(1,"X");
  Ax->SetLabelColor(1,"Y");
  Ax->SetLabelColor(1,"Z");
  Ax->SetXTitle("x (cm)");
  Ax->SetYTitle("y (cm)");
  Ax->SetZTitle("z (cm)");
  Ax->SetTitleOffset(2,"X");
  Ax->SetTitleOffset(2,"Y");
  Ax->SetTitleOffset(1.5,"Z");
  Ax->SetAxisRange(aux[0],aux[1],"X");
  Ax->SetAxisRange(aux[2],aux[3],"Y");
  Ax->SetAxisRange(aux[4],aux[5],"Z");
  G3.AddObject(Ax,"pad3");
 
  for(int i=0;i<vPH.size();++i){
  	double** valores= new double*[3];
  	for(int j=0; j<3; j++){
  		valores[j]= new double [2];
  		valores[j][0]=vPH[i].Xi()[j];
  		valores[j][1]=vPH[i].Xf()[j];
  	}
  
  	TPolyLine3D* t0 = new TPolyLine3D(2, valores[0], valores[1], valores[2]);
  
  	t0->SetLineColor(kRed+2);
  	t0->SetLineStyle(9);
  	t0->SetLineWidth(2);
  	G3.AddObject(t0,"pad3");
  }
  //..
  for(int i=0;i<vEL.size();++i){
  	double** valores= new double*[3];
  	for(int j=0; j<3; j++){
  		valores[j]= new double [2];
  		valores[j][0]=vEL[i].Xi()[j];
  		valores[j][1]=vEL[i].Xf()[j];
  	}
  	
  	TPolyLine3D* t0 = new TPolyLine3D(2, valores[0], valores[1], valores[2]);
  
  	t0->SetLineColor(kBlue+2);
    	t0->SetLineStyle(1);
    	t0->SetLineWidth(2);
  	G3.AddObject(t0,"pad3");
  }
  //..
  for(int i=0;i<vPO.size();++i){
  	double** valores= new double*[3];
  	for(int j=0; j<3; j++){
  		valores[j]= new double [2];
  		valores[j][0]=vPO[i].Xi()[j];
  		valores[j][1]=vPO[i].Xf()[j];
  	}
  	
  	TPolyLine3D* t0 = new TPolyLine3D(2, valores[0], valores[1], valores[2]);
  
   	t0->SetLineColor(kGreen+2);
    	t0->SetLineStyle(1);
    	t0->SetLineWidth(2);
  	G3.AddObject(t0,"pad3");
  }
  
  G3.AddObject(pad3);
  G3.ListObjects();
  G3.Draw();
  G3.Print("3D.pdf");
  
  G3.Clear();

  ///////////////////////////////////////////////////////////////////////////

  delete []aux;
}
