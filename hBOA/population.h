#ifndef _population_h_                                                               
#define _population_h_

#include "individual.h"

typedef struct {

  long N;
  int  n;
  int  numDiscrete;
  int  numContinuous;
  Individual *individual;

  char evaluated;

  long   best;
  long   worst;
  double avgFitness;

  long numOptimal;
  
  long *clusterSize;

  char *buffer;

} Population;

int allocatePopulation(Population *population, long N, int numDiscrete, int numContinuous);
int freePopulation(Population *population);
int generatePopulation(Population *population);
int evaluatePopulation(Population *population, Fitness *fitness);
int reevaluatePopulation(Population *population, Fitness *fitness);
int recomputeFitnessInfo(Population *population, Fitness *fitness);
int printPopulation(FILE *output, Population *population);
int swapIndividuals(Population *population, long i, long j);
int replaceIndividual(Population *population, long i, Individual *individual);
int shufflePopulation(Population *population);

int allocateUMF(double **p0, double **p1, int n);
int freeUMF(double **p0, double **p1, int n);
int calculateUMF(Population *population, double *p0, double *p1);
int calculateProportionateUMF(Population *population, double *p0, double *p1);
double calculateOrdering(double *p0, double *p1, int n);
int UMFCloserThanEpsilon(double *p0, double *p1, int n, double epsilon);

int allocateBMF(double ***p00, double ***p01, double ***p10, double ***p11, int n);
int freeBMF(double ***p00, double ***p01, double ***p10, double ***p11, int n);
int calculateBMF(Population *population, double **p00, double **p01, double **p10, double **p11);

int kMeansClusterPopulation(Population *population, int k, int numRestarts, char phenotypic);
int calculateClusterUMF(Population *p, int numClusters, double **p0, double **p1);
double fitnessVariance(Population *p);
void swapPopulations(Population *p, Population *q);

#endif
