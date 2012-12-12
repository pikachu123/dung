#include <string.h>

#include "random.h"
#include "memalloc.h"
//#include "fitness.h"
//#include "individual.h"
#include "population.h"

#define randomDigit() flipCoin()

//================================================================

int allocateIndividual(Individual *individual, int numDiscrete, int numContinuous)
{
  if (numDiscrete>0)
    individual->chromosome = new char[numDiscrete];
  else
    individual->chromosome = NULL;
  
  if (numContinuous>0)
    {
      individual->continuous = new double[numContinuous];
      individual->mutation   = new double[numContinuous];
    }
  else
    {
      individual->continuous = NULL;
      individual->mutation   = NULL;
    }
  
  individual->fCalculated=0;
  individual->cluster=-1;
  
  return 0;
}

//================================================================

int freeIndividual(Individual *individual)
{
  ///  printf("Freeing an individual\n");

  if (individual->chromosome)
    delete[] individual->chromosome;

  if (individual->continuous)
    {
      delete[] individual->continuous;
      delete[] individual->mutation;
    }

  ///  printf("Freed the individual\n");

  return 0;
}

//================================================================

int generateIndividual(Individual *individual, int numDiscrete, int numContinuous)
{
  for (int i=0; i<numDiscrete; i++)
    individual->chromosome[i]=randomDigit();

  for (int i=0; i<numContinuous; i++)
    {
      individual->continuous[i]=drand();
      individual->mutation[i]=0.05;
    }
  
  individual->fCalculated = 0;
  individual->cluster = -1;

  return 0;
}

//================================================================

int copyIndividual(Individual *dest, Individual *source, int numDiscrete, int numContinuous)
{
  if (numDiscrete>0)
    memcpy((void*)dest->chromosome,(void*)source->chromosome,numDiscrete);

  if (numContinuous>0)
    {
      memcpy((void*)dest->continuous,(void*)source->continuous,sizeof(double)*numContinuous);
      memcpy((void*)dest->mutation,(void*)source->mutation,sizeof(double)*numContinuous);
    };
  
  dest->f           = source->f;
  dest->goodBBs     = source->goodBBs;
  dest->fCalculated = source->fCalculated;
  dest->cluster     = source->cluster;
  
  return 0;
}

//================================================================

int evaluateIndividual(Individual *individual, Fitness *fitness, int numDiscrete, int numContinuous)
{
  if (individual->fCalculated)
    return 0;

  individual->goodBBs = 0;

  individual->f       = fitness->fitness(individual->chromosome,numDiscrete,individual->continuous,numContinuous) + additionalFitnessNoise();

  if (fitness->goodBBs)
    individual->goodBBs = fitness->goodBBs(individual->chromosome,numDiscrete,individual->continuous,numContinuous);
  
  fitnessCalled();
  
  individual->fCalculated=1;
    
  return 0;
}

//================================================================

int isBestIndividual(Individual *individual, Fitness *fitness, int numDiscrete, int numContinuous, int type)
{
  if (fitness->isBest==NULL)
    return 0;
  else
    return fitness->isBest(individual->chromosome,numDiscrete,individual->continuous,numContinuous, type,individual->f);
};

//================================================================

int printIndividual(FILE *output, Individual *individual, int numDiscrete, int numContinuous)
{
  for (int l=0; l<numDiscrete; l++)
    fprintf(output,"%u",int(individual->chromosome[l]));
  fprintf(output," ");
  for (int l=0; l<numContinuous; l++)
    fprintf(output,"%f ",individual->continuous[l]);
  fprintf(output,"\n");

  return 0;
}

//================================================================

int readIndividual(FILE *input, Individual *individual, int numDiscrete, int numContinuous)
{
  char c;

  for (int l=0; l<numDiscrete; l++)
    {
      do { c=fgetc(input); } while ((c!='0')&&(c!='1'));
      individual->chromosome[l]=c-'0';
    }
  
  printf("readIndividual is not fully implemented yet! exiting.\n");
  exit(-1);

  return 0;
}

//================================================================

//  char *getChromosome(Individual *individual)
//    {
//      return individual->chromosome[d];
//    }

char equalIndividuals(Individual *a, Individual *b, int numDiscrete, int numContinuous)
{
  int i;
  int result;

  result = 1;

  for (i=0; ((i<numDiscrete)&&(result)); i++)
    result = (a->chromosome[i]==b->chromosome[i]);
 
  for (i=0; ((i<numContinuous)&&(result)); i++)
    result = (a->continuous[i]==b->continuous[i]);
 
  return result;
};

double hammingDistance(Individual *a, Individual *b, int n)
{
  int j;
  double dist;

  dist=0;
  for (j=0; j<n; j++)
    if (a->chromosome[j]!=b->chromosome[j])
      dist++;

  return dist;
};
