#ifndef _priors_h_
#define _priors_h_

#include "population.h"
#include "fitness.h"
#include "graph.h"

typedef int PriorNetworkConstructor(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N); 

typedef struct {

  char                    *description;
  PriorNetworkConstructor *constructor;

} PriorNetworkSourceDescription;

int setPriorNetworkSource(int source);
char *getPriorNetworkDescription(int source);
char isPriorNetworkDefined();
int allocatePriorNetwork(OrientedGraph ***priorNetwork, int n, int d);
int freePriorNetwork(OrientedGraph ***priorNetwork, int n, int d);
int constructPriorNetwork(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N);

int emptyPriorNetwork(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N);
int LINCPriorNetwork(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N);

#endif
