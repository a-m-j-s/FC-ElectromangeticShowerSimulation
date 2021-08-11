#ifndef __Func1D__
#define __Func1D__

#include"TF1.h"
#include"cFCgraphics.h"

class Func1D{

 public:
  Func1D();
  Func1D(TF1 * fp);
  ~Func1D();
  void Draw();
  double Evaluate(double);
  
 protected:
  TF1*p;
};
#endif
