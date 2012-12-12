#ifndef _bic_h_
#define _bic_h_

//#include "hBoa.h"
#include "population.h"
#include "recombination.h"

int bicComputeAdditionGains( int i, 
			     double oldContribution,
			     double **gain, 
			     int *updateIdx, 
			     int numUpdated, 
			     int *parentList,
			     int numParents, 
			     Population *P,
			     RecombinationParams *params);

int bicComputeRemovalGains( int i, 
			    double oldContribution,
			    double **removalGain, 
			    int *parentList,
			    int numParents, 
			    Population *P, 
			    RecombinationParams *params);

double bicComputeIsolatedNodeContribution(int i, Population *P);

#endif
