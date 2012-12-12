#include "population.h"

typedef struct DiscretizationParams {
  int discretizationType;
  int bitsPerVariable;
};

class Discretization {

 private:
  int    bitsPerVariable;
  long **index;
  long ***index2;
  int    binSize;
  long **clusterSize;
  long **fhhBound;
  long **fhhStart;
  long **fhhEnd;

 public:

  int numContinuous;
  int numDiscrete;
  int numBins;

  Population *continuousPopulation;
  Population *discretePopulation;

  double **lo, **hi;

  Discretization();
  ~Discretization();

  void discretize(Population *p, Population *q, DiscretizationParams *discretizationParams);
  void undiscretize(Population *p, Population *q, DiscretizationParams *params);

  void copyDiscreteVariables(char *dest, char *src, int n);


  void allocateBounds();
  void computeBoundsFHH(Population *p);
  void freeBounds();
  void discretizeFHH();
  void discretizeFWH();
  void undiscretizeFHH(Individual *from, Individual *to);
  void undiscretizeFWH(Individual *from, Individual *to);
  void discretizeContinuousVariablesFHH(char *dest, double *x, int n);

  void cluster(double *centers, long *centerSize, int k, int numClusters);
  void discretizeKMeans();
  void undiscretizeKMeans(Individual *from, Individual *to);

  void setCompare(Population *p, int whichVariable);

/*   Population getDiscretePopulation(); */
/*   Population getContinuousPopulation(); */
};

int discretizationCompare(const void *x, const void *y);
int doubleCompare(const void *a, const void *b);
