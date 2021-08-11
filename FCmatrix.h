#ifndef __FCmatrix__
#define __FCmatrix__

#include "Vec1.h"
#include <vector>
#include <string>

using namespace std;

class FCmatrix {

public:
FCmatrix();
FCmatrix(double** fM, int fm, int fn); //matrix fm x fn
FCmatrix(double* fM, int fm, int fn);
FCmatrix(vector<Vec>);
virtual ~FCmatrix();

// operators
virtual Vec& operator[] (int) = 0;

// methods
virtual Vec GetRow(int i) const = 0; // retrieve row i
virtual Vec GetCol(int i) const = 0; // retrieve column i
virtual double Determinant() = 0;
virtual void Print();

// row max element index
virtual int GetRowMax(int) = 0;

// row max element index
virtual int GetColMax(int) = 0;

virtual void swapRows(int, int) = 0; // swap rows i,j

protected:
vector<Vec> M;
string classname;
};

#endif
