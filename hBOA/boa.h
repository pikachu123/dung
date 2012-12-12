#ifndef _boa_h_
#define _boa_h_

#include "population.h"
#include "recombination.h"
#include "graph.h"

#define K2_METRIC  0
#define MDL_METRIC 1

typedef int BoaComputeAdditionGains( int i, 
				  double oldContribution,
				  double **gain, 
				  int *updateIdx, 
				  int numUpdated, 
				  int *parentList,
				  int numParents, 
				  Population *P,
				  RecombinationParams *params);

typedef int BoaComputeRemovalGains( int i, 
				 double oldContribution,
				 double **removalGain, 
				 int *parentList,
				 int numParents, 
				 Population *P,
				 RecombinationParams *params);

typedef int BoaComputeReversalGains();

typedef double BoaComputeIsolatedNodeContribution( int i,
						Population *P);

typedef struct {

  char                            *description;
  BoaComputeAdditionGains            *boaComputeAdditionGains;
  BoaComputeRemovalGains             *boaComputeRemovalGains;
  BoaComputeReversalGains            *boaComputeReversalGains;
  BoaComputeIsolatedNodeContribution *boaComputeIsolatedNodeContribution;

} MetricDescription;

int boaRecombination(Population *parents, 
		     Population *children, 
		     long M, 
		     RecombinationParams *params);

int boaRecomputeAdditionGains(int from,
			   int to,
			   double gain,
			   double *nodeContribution,
			   double **gainAddition, 
			   double **gainRemoval, 
			   double **gainReversal, 
			   AcyclicOrientedGraph *G, 
			   Population *population,
			   RecombinationParams *params);

int boaRecomputeRemovalGains(int from,
			  int to,
			  double gain,
			  double *nodeContribution,
			  double **gainAddition, 
			  double **gainRemoval, 
			  double **gainReversal, 
			  AcyclicOrientedGraph *G, 
			  Population *population,
			  RecombinationParams *params);

int boaRecomputeReversalGains(int from,
			   int to,
			   double gain,
			   double *nodeContribution,
			   double **gainAddition, 
			   double **gainRemoval, 
			   double **gainReversal, 
			   AcyclicOrientedGraph *G, 
			   Population *population,
			   RecombinationParams *params);

int boaUpdateAfterAddition(int newFrom, 
			int newTo, 
			char **additionPossible, 
			char **removalPossible, 
			char **reversalPossible, 
			AcyclicOrientedGraph *G, 
			char *full,
			int maxIncoming);

int boaUpdateAfterRemoval(int newFrom, 
		       int newTo, 
		       char **additionPossible, 
		       char **removalPossible, 
		       char **reversalPossible, 
		       AcyclicOrientedGraph *G, 
		       char *full,
		       int maxIncoming);

int boaUpdateAfterReversal(int newFrom, 
			int newTo, 
			char **additionPossible, 
			char **removalPossible, 
			char **reversalPossible, 
			AcyclicOrientedGraph *G, 
			char *full,
			int maxIncoming);

int boaRecomputeNodeContribution(int i, double *nodeContribution, AcyclicOrientedGraph *G, Population *parents);

double boaGetPenalty(int from, int to);

//double computeGain(int k, int l, Population *parents, AcyclicOrientedGraph *G, RecombinationParams *params);
int boaSetMetric(int n);
char *boaGetMetricDescription(int n);

#endif
