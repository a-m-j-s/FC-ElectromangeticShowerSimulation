#include "FCmatrixFull.h"
#include "FCmatrix.h"
#include "Vec1.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>

//#define DEBUG

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FCmatrixFull:: FCmatrixFull(){
classname="FCmatrixFull";
nc=0;
nr=0;
rowindices=new int[1];
rowindices[0]=0;
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}


FCmatrixFull::FCmatrixFull(double** fM, int fm, int fn){//matrix fm x fn
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
classname="FCmatrixFull";
nc=fn;
nr=fm;
rowindices=new int[nr];
rowindices[0]=0;
for(int i=0; i<fm; i++)
{	rowindices[i]=i;
 	Vec linha(fn, fM[i]);
 	//linha.Print();
 	M.push_back(linha);

}

}



FCmatrixFull::FCmatrixFull(double* fM, int fm, int fn){//matrix fm x fn
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
classname="FCmatrixFull";
nc=fn;
nr=fm;
rowindices=new int[nr];
for(int i=0; i<fm; i++)
{	rowindices[i]=i;
 	double temp[fn];
 	for (int k=0; k<fn;k++){
 		temp[k]=fM[i*fn+k];
	}
 	Vec linha(fn, temp);
 	//linha.Print();
 	M.push_back(linha);
}

}

FCmatrixFull::FCmatrixFull(vector<Vec> matrix){//matrix fm x fn
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
classname="FCmatrixFull";
M=matrix;
nr=matrix.size();
nc=M[0].size();
rowindices=new int[nr];
for (int i=0; i<nr; i++)
{ rowindices[i]=i;
}

//cout<<"1  "<<nr<<"    2  "<<nc<<endl;
}

FCmatrixFull:: ~FCmatrixFull(){
	delete[]rowindices;
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}

FCmatrixFull::FCmatrixFull(const FCmatrix& outra){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
classname="FCmatrixFull";

M.clear();

for(int i=0;i<(outra.GetCol(0)).size();++i){
       M.push_back(outra.GetRow(i));
      }

nr=M.size();
nc=M[0].size();

rowindices=new int[nr];
for (int i=0; i<nr; i++)
{ rowindices[i]=i;}

}

FCmatrixFull::FCmatrixFull(const FCmatrixFull& outra){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
classname="FCmatrixFull";
rowindices=new int[nr];
nr=outra.nr;
nc=outra.nc;
M=outra.M;
for (int i=0; i<outra.M.size(); i++)
{ rowindices[i]=outra.rowindices[i];
  //cout<<i<<endl;
}
}

/////////////////////////////////////////////////////////////////////////////////

Vec FCmatrixFull::GetRow(int i)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
if (i<nr&&i>-1)
{Vec linha = M[rowindices[i]];
return linha;}
else 
{printf("linha inexistente\n");
Vec linha;
return linha;}
} 

Vec FCmatrixFull::GetCol(int i)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
if (i<nc&&i>-1){
double *col1 = new double[nr];
for (int k=0; k<nr; k++){
col1[k]=M[rowindices[k]].At(i);
}
Vec col2(nr, col1);
return col2;}
else
{printf("coluna inexistente\n");
Vec col;
return col;}
} 

void FCmatrixFull::swapRows(int a,int b){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
if (a<nr&& b<nr && a>-1 &&b>-1){
int temp=rowindices[a];
rowindices[a]=rowindices[b];
rowindices[b]=temp;}
else cout<<"linhas erradas"<<endl;
}

void FCmatrixFull::Print(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
for (int i=0; i<nr; i++)
{
	M[rowindices[i]].Print();
}
cout<<endl;
}


int FCmatrixFull::GetRowMax(int i=0){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
int imax=0;
if (i<nr&&i>-1)
{
	double max= abs(M[rowindices[i]].At(0));
	for (int k=0; k<nc; k++){
		double valor = abs(M[rowindices[i]].At(k));
		if (valor>max)
 			{max=valor;
  			imax=k;}
	}
return imax;}
else 
{printf("linha inexistente\n");
return imax;}
}


int FCmatrixFull::GetColMax(int i=0){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
int imax=0;
if (i<nc&&i>-1)
{
	double max= abs(M[rowindices[0]].At(i));
	for (int k=0; k<nr; k++){
		double valor = abs(M[rowindices[k]].At(i));
		if (valor>max)
 			{max=valor;
  			imax=k;}
	}
return imax;}
else 
{printf("coluna inexistente\n");
return imax;}
}

double FCmatrixFull::Determinant(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
	int t=0;
	FCmatrixFull M1(M);
	if((M1.GetRow(0)).size()==(M1.GetCol(0)).size()){
	  for(int j=0;j<(M1.GetRow(0)).size();++j){
	    double pos=M1.GetColMax(j);
	    if(pos>j){
	      M1.swapRows(j,pos);
	      M1.Print();
	      ++t;
	    }
	  }
	}

	for(int k=0;k<(M1.GetRow(0)).size()-1;++k){
	  for(int i=k+1;i<(M1.GetRow(0)).size();++i){
	    double la=((M1.GetRow(i)).At(k))/((M1.GetRow(k)).At(k));
	    M1[i]=M1[i]-(M1[k]*la);       		
	  }	
	}

	double det=1;
	for(int a=0;a<(M1.GetRow(0)).size();++a){
	  det=det*M1[a][a];
	}
	cout<<det<<endl;
	if(t%2!=0){det=-1*det;}
	return det;
	  
    	//algoritmo alternativo---------------------	
	/*  	if(M.size()!=M[0].size()){
  	cout<<"matriz nao quadrada"<<endl;
  	return 0;
  }
  if(M.size()==1)
  	return M[0][0];

  if(M.size()==2 && M[0].size()==2){
    	double det=M[rowindices[0]].At(0)*M[rowindices[1]].At(1)-M[rowindices[0]].At(1)*M[rowindices[1]].At(0);
    	return det;
  }
  else{
    	double det=0;
    	for(int i=0;i<M[0].size();++i){
      		vector<Vec>Mat;
      		for(int it1=1; it1<M.size();++it1){
			double* d1=new double[M[0].size()-1];
			for(int it2=0;it2<M[0].size();++it2){
	  			if(it2<i){d1[it2]=M[rowindices[it1]].At(it2);}
	  			else if(it2>i){d1[it2-1]=M[rowindices[it1]].At(it2);}
			}
			Vec a(M[0].size()-1,d1);
			Mat.push_back(a);
			delete []d1;
     		}
      		FCmatrixFull Matriz(Mat);
      		//Matriz.Print();
      		if(i%2==0){det=det+M[rowindices[0]].At(i)*Matriz.Determinant();}
      		else{det=det-M[rowindices[0]].At(i)*Matriz.Determinant();}
    	}
   	 return det;
	 }*/
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Vec& FCmatrixFull::operator[] (int i){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
if (i<nr&&i>-1)
return M[rowindices[i]];
else 
{printf("linha inexistente\n");
return M[rowindices[0]];}
}


FCmatrixFull& FCmatrixFull::operator=(const FCmatrixFull& a){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif

	nr=a.nr;
	nc=a.nc;
    	M=a.M;
    	classname=a.classname;
	delete[]rowindices;
	rowindices=new int[nr];
	for (int i=0; i<nr; i++)
		{rowindices[i]=a.rowindices[i];}
    return *this;
}




FCmatrixFull FCmatrixFull::operator+(const FCmatrix& a)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
    if(nc != (a.GetRow(0)).size() || nr != (a.GetCol(0)).size()){
    cout<<"Erro-tamanhos diferentes"<<endl;
    FCmatrixFull novo;
    return novo;
    }
    else{
    vector<Vec> soma;
    for(int i=0;i<nr;++i){
    Vec novo1 = M[rowindices[i]]+a.GetRow(i);
    soma.push_back(novo1);
    }
    FCmatrixFull novo(soma);
    return novo;}
}


FCmatrixFull FCmatrixFull::operator-(const FCmatrix& a)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
    if(nc != (a.GetRow(0)).size() || nr != (a.GetCol(0)).size()){
    cout<<"Erro-tamanhos diferentes"<<endl;
    FCmatrixFull novo;
    return novo;
    }
    else{
    vector<Vec> soma;
    for(int i=0;i<nr;++i){
    Vec novo1 = M[rowindices[i]]-a.GetRow(i);
    soma.push_back(novo1);
    }
    FCmatrixFull novo(soma);
    return novo;}
}



FCmatrixFull FCmatrixFull::operator*(double c)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
    vector<Vec> soma;
    for(int i=0;i<nr;++i){
    Vec novo1 = M[rowindices[i]]*c;
    soma.push_back(novo1);
    }
    FCmatrixFull novo(soma);
    return novo;
}

FCmatrixFull FCmatrixFull::operator*(const FCmatrix& a)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
    if(nc != (a.GetCol(0)).size()){
    cout<<"Erro-dimensoes incompativeis"<<endl;
    FCmatrixFull novo;
    return novo;
    }
    else{
    double **temp=new double*[nr];
	for (int i=0; i<nr;i++){
		temp[i]=new double[(a.GetRow(0)).size()];
		for (int j=0; j<(a.GetRow(0)).size(); j++){
			temp[i][j]=M[rowindices[i]].dot(a.GetCol(j));
		}
	}
    FCmatrixFull novo(temp, nr,(a.GetRow(0)).size() );
    for (int i=0; i<nr; i++){delete temp[i];}
    delete[] temp;
    return novo;}
}

Vec FCmatrixFull::operator*(const Vec& a)const{
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
    if(nc != a.size()){
    cout<<"Erro-tamanhos diferentes"<<endl;
    Vec novo;
    return novo;
    }
    else{
    Vec mult(nr);
    double *temp=new double[nr];
    for(int i=0;i<nr;++i){
    	temp[i]= M[rowindices[i]].dot(a);
    }
    mult.SetEntries(nr, temp);
    delete[] temp;
    return mult;}
}


