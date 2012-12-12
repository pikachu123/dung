#ifndef _individual_h_
#define _individual_h_

#include <stdio.h>

#include "reordering.h"
#include "fitness.h"

// some macros

//#define getChromosome(guy,d) ((guy)->chromosome[d])
//#define getDigit(guy,i,j) ((guy)->chromosome[i][j])
//#define setDigit(guy,i,j,value) ((guy)->chromosome[i][j] = char(value))

// other definitions

class Individual {

 public:

  double f;
  char   fCalculated;
  char   *chromosome;
  double *continuous;

  double *mutation;

  int    cluster;
  double goodBBs;

  char operator< (Individual x)
    { 
      return (f<x.f);
    };

  char operator<= (Individual x)
    {
      return (f<=x.f);
    };
  
  char operator> (Individual x)
    {
      return (f>x.f);
    };
  
  char operator>= (Individual x)
    { 
      return (f>=x.f);
    };

  char operator== (Individual x)
    {
      return (f==x.f);
    };
};

int allocateIndividual(Individual *individual, int numDiscrete, int numContinuous);
int freeIndividual(Individual *individual);
int generateIndividual(Individual *individual, int numDiscrete, int numContinuous);
int copyIndividual(Individual *dest, Individual *source, int numDiscrete, int numContinuous);
int evaluateIndividual(Individual *individual, Fitness *fitness, int numDiscrete, int numContinuous);
int isBestIndividual(Individual *individual, Fitness *fitness, int numDiscrete, int numContinuous, int type);
//char getDigit(Individual *individual, int i, int j);
//char setDigit(Individual *individual, int i, int j, char value);
int printIndividual(FILE *output, Individual *individual, int numDiscrete, int numContinuous);
int readIndividual(FILE *input, Individual *individual, int numDiscrete, int numContinuous);
//char *getChromosome(Individual *individual);
char equalIndividuals(Individual *a, Individual *b, int numDiscrete, int numContinuous);

double hammingDistance(Individual *a, Individual *b, int n);

#endif
