#ifndef _recombination_h_
#define _recombination_h_

#include "graph.h"

#define BOA_RECOMBINATION      0
#define DECISION_TREE_BOA_RECOMBINATION 1
#define COPY_ALL 2

typedef struct {
  char   displayDependencies;
  double dependencyTreshold;

  int    maxIncoming;

  int    crossoverMethod;
  double crossoverRate;
  double mutationRate;

  int    numClusters;
  int    numRestarts;
  char   fitnessProportionalClusterReproduction;
  char   phenotypicClustering;

  char          priorNetwork;
  double        logKappa;
  double        logBeta;
  char          injectGoodGuys;

  int           allowAdditions;
  int           allowRemovals;
  int           allowReversals;
  int           allowJoints;
  int           useDefaultTables;

  int           dBOAMetricN;

  int           maxGroupSize;

  AcyclicOrientedGraph *G;
  double        recombinationBoltzmannBeta;
  char          useBoltzmannFrequencies;

  int           blockSize;
} RecombinationParams;

typedef int RecombinationFunction (Population *parents, Population *children, long M, RecombinationParams *recombinationParams);

typedef struct {
  char                  *description;
  RecombinationFunction *recombination;
} Recombination;

void *getRecombination(int n);
char *getRecombinationDesc(int which);

//int setPriorNetwork(Fitness *fitness, RecombinationParams *recombinationParams);

#endif
