#include "math.h"
 
double leftRightGaussians(double *x, double *par){
  double xx =x[0];
  double f;
  if(xx<par[0])
    f=par[1]*exp(-0.5*((xx-par[0])/par[2])*((xx-par[0])/par[2]));
  else
    f=par[3]*exp(-0.5*((xx-par[0])/par[4])*((xx-par[0])/par[4]));
  return f;
}
