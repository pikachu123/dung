#include "population.h"

class Discretization {

  int      n;
  double *lo;
  double *hi;

  Discretization(Population *p, discretization Grain);

  virtual void discretize(Population *p)=0;
  virtual void undiscretize(Population *p)=0;

  ~Discretization();
};
