#include "population.h"
#include "recombination.h"

double logS(int i, Population *P, AcyclicOrientedGraph *G);
double logGain(int i, int *parents, int n, Population *P);

int K2ComputeLogGains( int i, 
		       double oldContribution,
		       double **gain, 
		       int *updateIdx, 
		       int numUpdated, 
		       int *parentList,
		       int numParents, 
		       Population *P,
		       RecombinationParams *params);

double K2ComputeIsolatedNodeContribution(int i, Population *P);
