#ifndef __FCmatrixFull__
#define __FCmatrixFull__

#include "FCmatrix.h"
#include "Vec1.h"

class FCmatrixFull : public FCmatrix {
public:

// constructors
FCmatrixFull();
FCmatrixFull(double**, int, int); //matrix fm x fn
FCmatrixFull(double*, int, int);
FCmatrixFull(vector<Vec>);
~FCmatrixFull();

// copy constructor
FCmatrixFull(const FCmatrix&);
FCmatrixFull(const FCmatrixFull&);

// operators
FCmatrixFull& operator=(const FCmatrixFull&);
FCmatrixFull operator+(const FCmatrix&)const; // add 2 matrices of any kind
FCmatrixFull operator-(const FCmatrix&)const; // sub 2 matrices of any kind
FCmatrixFull operator*(const FCmatrix&)const; // mul 2 matrices of any kind
FCmatrixFull operator*(double) const; // mul matrix of any kind by scalar
Vec operator*(const Vec&)const; // mul matrix by Vec
Vec& operator[] (int);

// virtual inherited
Vec GetRow(int i) const; // retrieve row i
Vec GetCol(int i) const; // retrieve column i

// row max element index
int GetRowMax(int);

// row max element index
int GetColMax(int);

double Determinant();

void swapRows(int,int);
void Print();

private:
 int nr;//rows
 int nc; //col
 int *rowindices; // row indices (0,1,2,..
};

#endif
