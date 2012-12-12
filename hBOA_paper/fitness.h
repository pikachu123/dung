#ifndef _fitness_h_
#define _fitness_h_

#define HEAVY 1
#define LIGHT 0

void *getFitness(int n);
char *getFitnessDesc(int n);
int initRealValuedFitness(double *params);
double realValue(char *x, int n);
char isRealValued();

int priorConnected(int i, int j);

int resetFitnessCalls(void);
long fitnessCalled(void);
long fitnessCalled(long numCalls);
long getFitnessCalls(void);

typedef double FitnessFunction(char *x, int numDiscrete, double *c, int numContinuous);
typedef char   IsBest(char *x, int numDiscrete, double *c, int numContinuous, char type, double f);
typedef double GoodBBs(char *x, int numDiscrete, double *c, int numContinuous);
typedef int    InitFitness(double *params);
typedef int    LoadParameters(char *filename, double *params);
typedef int    DoneFitness(void);
//typedef int    InitPriorNetwork(double *params);
typedef int    InjectGoodGuys(void *p);
typedef long   GetNumberOfOptima(int numDiscrete, int numContinuous);
typedef long   WhichOptimum(char *x, int numDiscrete, double *c, int numContinuous);

int setFitnessNoiseVariance(double variance);
double getFitnessNoiseVariance();
double additionalFitnessNoise();

typedef struct {
  char              *description;

  FitnessFunction   *fitness;
  IsBest            *isBest;
  GoodBBs           *goodBBs;
  InitFitness       *init;
  DoneFitness       *doneFitness;
  LoadParameters    *load;
  InjectGoodGuys    *injectTheGood;
  GetNumberOfOptima *getNumberOfOptima;
  WhichOptimum      *whichOptimum;
} Fitness;

#endif
