#include "Vec1.h"
#include "FCmatrix.h"
#include "FCmatrixFull.h"
#include <vector>


using namespace std;

class EqSolver{

public:
EqSolver();
EqSolver(const FCmatrix&, const Vec&); // matriz M e vector de constantes B
~EqSolver();

void Print();

// set
void SetConstants(const Vec&);
void SetMatrix(const FCmatrix&);

//eliminação de Gauss:
//resolução do sistema pelo método de eliminação de Gauss
Vec GaussEliminationSolver();
Vec LUdecompositionSolver();

private:
//decomposição LU com |L|=1
void LUdecomposition(FCmatrixFull&, vector<int>&, int&);
/* return triangular matrix and changed vector of constants */
void GaussElimination(FCmatrixFull&, Vec&, int&);

FCmatrix *M; //matriz de coeffs

Vec b; //vector de constantes
};
