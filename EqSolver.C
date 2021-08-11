#include "EqSolver.h"
#include "Vec1.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include <iostream>
#include <stdio.h>
#include <vector>

#define DEBUG

EqSolver:: EqSolver(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}

EqSolver::EqSolver(const FCmatrix& Min, const Vec& vec){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
 b=vec;
 M= new FCmatrixFull(Min);
 //M->Print();

}

EqSolver:: ~EqSolver(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
delete M;
}


void EqSolver::Print(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
M->Print();
b.Print();
}

void EqSolver::SetConstants(const Vec& vec){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
b=vec;
}

void EqSolver::SetMatrix(const FCmatrix& mat){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
delete M;
 M= new FCmatrixFull(mat);
}
/////////////////////////////////////////////////////////////////////////////

void EqSolver::GaussElimination(FCmatrixFull& Mre, Vec& vre, int& t){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
//cout<<"inicio Gauss elimination"<<endl<<endl<<endl;
Mre=*M;
//Mre.Print();
vre=b;
 if((Mre.GetRow(0)).size()==(Mre.GetCol(0)).size() && Mre.Determinant()!=0){
	//if ((Mre.GetRow(0)).size()==1){
	//	vre[0]=b[0]/Mre[0][0];
	//}
	//else{
		t=1;
		//cout<<"entra condicao"<<endl<<endl<<endl;
		//tornar a matriz o mas dominante possivel
	    	for(int j=0;j<(Mre.GetRow(0)).size();++j){
	      		double pos=Mre.GetColMax(j);
	      		if(pos>j){//nao trocar com as linhas acima para que o primeiro pivot seja maior e !=0
				Mre.swapRows(j,pos);
				vre.swap(j,pos);
	      		}
	    	}
		//cout<<"fim pivot"<<endl<<endl<<endl;
	      	//Mre.Print();
	    	//gauss elimination  
	    	for(int k=0;k<(Mre.GetRow(0)).size()-1;++k){
	      		for(int i=k+1;i<(Mre.GetRow(0)).size();++i){
				double la=((Mre.GetRow(i)).At(k))/((Mre.GetRow(k)).At(k));
				Mre[i]=Mre[i]-(Mre[k]*la);
				vre[i]=vre[i]-(vre[k]*la);
	      		}	
	    	}
	    	Mre.Print();
	    	vre.Print();
	//}
  }
else{cout<<"gauss NAO efetuada"<<endl;}
}


Vec EqSolver::GaussEliminationSolver(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

Vec V1;
FCmatrixFull M1;
int t=0;
GaussElimination(M1, V1, t);

//cout<<"inicio solver"<<endl;
//M.Print();
//cout<<"vetor"<<endl;
//b.Print();
//cout<<"bora la"<<endl;
if (t){
for(int i=(M1.GetCol(0)).size()-1;i>=0;--i){
    double sum =0;
    for(int j=i+1;j<(M1.GetCol(0)).size();++j){
      sum=sum+V1[j]*M1[i][j];
      //cout<<sum<<endl;
    }
    V1[i]=(V1[i]-sum)/M1[i][i];
  }
  //V1.Print();
}
return V1;
}

///////////////////////////////////////////////////////


void EqSolver::LUdecomposition(FCmatrixFull& Mre, vector<int> &index, int& t){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
  Mre=*M;
  int n =(Mre.GetRow(0)).size();
  for(int i=0;i<n;++i){
    index.push_back(i);
  }
//cout<<"teste 1"<<endl;
  if((Mre.GetRow(0)).size()==(Mre.GetCol(0)).size() && Mre.Determinant()!=0){
	//if ((Mre.GetRow(0)).size()==1){
		//index[0]=b[0]/Mre[0][0];
	
	//}
	//else{
    		t=1;
		//tornar a matriz o mas dominante possivel
		for(int j=0;j<(Mre.GetRow(0)).size();++j){
			double pos=Mre.GetColMax(j);
			if(pos>j){//nao trocar com as linhas acima para que o primeiro pivot seja maior e !=0
				Mre.swapRows(j,pos);
				int temp=index[j];
				index[j]=index[pos];
				index[pos]=temp;
		      	}
    		}
		//Mre.Print();

		//gauss elimination  
		double *la=new double[(int)(n*n/2)];
		int a=0;
		for(int k=0;k<n-1;++k){
      			for(int i=k+1;i<n;++i){
				Mre.Print();
				la[a]=Mre[i][k]/Mre[k][k];
				cout<<Mre[i][k]<<" "<<Mre[k][k]<<endl;
				Mre[i]=Mre[i]-(Mre[k]*la[a]);
				//cout<<" la="<<la<<" "<<k<<endl;
				//printf("M[%d,%d]=%f  %f\n",i,k,Mre[i][k],la[a]);
				++a;
      			}
    		}
   		//Mre.Print();
   		a=0;
   		for(int k=0;k<n-1;++k){
      			for(int i=k+1;i<n;++i){
				Mre[i][k]=la[a];
				++a;
      			}
    		}
  		//Mre.Print();
  		delete []la;
  	//}
}
else{
	cout<<"LUdecomposition NAO efetuada"<<endl;}
}

Vec EqSolver::LUdecompositionSolver(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
Vec a;

int t=0;
  vector<int> v;
  FCmatrixFull M1;
  LUdecomposition(M1,v, t);

cout<<"inicio solver!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1"<<endl;
if(t){
  Vec V1((M1.GetRow(0)).size());
  Vec y((M1.GetRow(0)).size());

for(int i=0;i<(M1.GetRow(0)).size();++i){
    y[i]=b[v[i]];
  }

  for(int k=0;k<(M1.GetRow(0)).size();++k){
    double sum=0;
    for(int i=0;i<k;++i){
      sum=sum+y[i]*M1[k][i];
    }
    y[k]=y[k]-sum;
  }
  V1=y;
  for(int k=(M1.GetRow(0)).size()-1;k>=0;k--){
    double sum =0;
    for(int i=k+1;i<(M1.GetRow(0)).size();++i){
      sum=sum+V1[i]*M1[k][i];
    }
    V1[k]=(y[k]-sum)/M1[k][k];
    }
return V1;
}
   


return a;
}
