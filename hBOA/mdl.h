#ifndef _mdl_h_
#define _mdl_h_

#include "hBoa.h"
#include "population.h"
#include "recombination.h"

int mdlComputeAdditionGains( int i, 
			     double oldContribution,
			     double **gain, 
			     int *updateIdx, 
			     int numUpdated, 
			     int *parentList,
			     int numParents, 
			     Population *P,
			     RecombinationParams *params);

int mdlComputeRemovalGains( int i, 
			    double oldContribution,
			    double **removalGain, 
			    int *parentList,
			    int numParents, 
			    Population *P, 
			    RecombinationParams *params);

double mdlComputeIsolatedNodeContribution(int i, Population *P);

int mdlComputeGroupAdditionGains( int i,
				  int numVars,
				  GroupInformation *groupInformation,
				  OperatorGain     *operatorGain,
				  AcyclicOrientedGraph *G,
				  int *updateIdx, 
				  int numUpdated, 
				  int *parentList,
				  int numParents, 
				  Population *P,
				  RecombinationParams *params);

int mdlComputeGroupRemovalGains( int i,
				 int numVars,
				 GroupInformation *groupInformation,
				 OperatorGain     *operatorGain,
				 AcyclicOrientedGraph *G,
				 int *parentList,
				 int numParents, 
				 Population *P,
				 RecombinationParams *params);

double mdlComputeGroupContribution(int i, 
				   int numVars,
				   GroupInformation *groupInformation,
				   AcyclicOrientedGraph *G,
				   Population *P, 
				   RecombinationParams *params);

double buildDefaultTable(int n, long N, double **p, double *tp, long numParentConfigurations, int numParents);

#endif
