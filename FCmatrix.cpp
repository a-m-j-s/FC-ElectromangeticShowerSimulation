#include "FCmatrix.h"
#include "Vec1.h"
#include <iostream>
#include <stdio.h>
//#define DEBUG
using namespace std;


FCmatrix:: FCmatrix(){
  classname="FCmatrix";
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}


FCmatrix::FCmatrix(double** fM, int fm, int fn){//matrix fm x fn
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
for(int i=0; i<fm; i++)
{
 Vec linha(fn, fM[i]);
 //linha.Print();
 M.push_back(linha);

}

}


FCmatrix::FCmatrix(double* fM, int fm, int fn){//matrix fm x fn
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
for(int i=0; i<fm; i++)
{
 //cout<<"teste "<<i<<endl;
 double temp[fn];
 	for (int k=0; k<fn;k++){
 		temp[k]=fM[i*fn+k];
	}
 	Vec linha(fn, temp);
 	//linha.Print();
 	M.push_back(linha);
}

}

FCmatrix::FCmatrix(vector<Vec> matrix){//matrix fm x fn
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
for (int i=0; i<matrix.size(); i++)
{
//matrix[i].Print();
M.push_back(matrix[i]);
}

}


FCmatrix:: ~FCmatrix(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
}
////////////////////////////////////////////////////////////////////////////////////////////
void FCmatrix::Print(){
	#if defined(DEBUG)
        printf("[%s]\n", __PRETTY_FUNCTION__);
	#endif
for (int i=0; i<M.size(); i++)
{
	M[i].Print();
}
cout<<endl;
}

