#ifndef __material__
#define __material__
#include<string>
using namespace std;

struct material{
  double Z;
  double ro;
  double A;
  double C0;
  double U0;
  double U1;
  double a;
  double m;
  string name;
};
#endif
