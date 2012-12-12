#include <stdlib.h>

#include "population.h"
#include "recombination.h"
#include "boa.h"
#include "decisionGraphBoa.h"
#include "copy.h"

#define numRecombinations 3

Recombination recombinationInfo[numRecombinations]={
  {"BOA",&boaRecombination},
  {"Hierarchical BOA",&boaWithDecisionGraphsRecombination},
  {"Copy (reproduce all)",&copyRecombination}
};

//==========================================================

void *getRecombination(int n)
{
  if ((n<0)||(n>=numRecombinations))
    {
      fprintf(stderr,"ERROR: Unknown recombination method.");
      exit(-1);
    }
  
  return &(recombinationInfo[n]);
}

//==========================================================

char *getRecombinationDesc(int which)
{
  return recombinationInfo[which].description;
}


// =====================================================================
/*
int setPriorNetwork(Fitness *fitness, RecombinationParams *params)
{
  params->priorNetworkConnected = priorConnected;

  return 0;
};
*/
