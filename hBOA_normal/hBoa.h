#ifndef _hBoa_h_
#define _hBoa_h_

#include "population.h"
#include "recombination.h"
#include "graph.h"
#include "frequencyTree.h"

typedef int HboaComputeAdditionGains();
typedef int HboaComputeRemovalGains();
typedef int HboaComputeReversalGains();
typedef int HboaComputeJointGains();
typedef int HboaComputeIsolatedNodeContribution();

typedef struct {
  int              numGroups;
  int              maxGroupSize;

  int             *groupSize;
  int            **groupIndex;

  int             *indexSize;
  int            **index;
  int             *parentIndexSize;

  FrequencyTree  **frequencyTree;
  int             *numInstances;
  char          ***instances;
  double         **frequencies;
  char            *full;

  double *contribution;
} GroupInformation;

typedef struct {
  char **addition;
  char **removal;
  char **reversal;
  char **joint;
} OperatorApplicability;

typedef struct {
  double **addition;
  double **removal;
  double **reversal;
  double **joint;
} OperatorGain;

int hboaRecombination(Population *parents, 
		      Population *children, 
		      long M, 
		      RecombinationParams *params);

int hboaRecomputeAdditionGains(int from,
			       int to,
			       int numVars,
			       GroupInformation *groupInformation,
			       OperatorGain *operatorGain,
			       OperatorApplicability *operatorApplicability,
			       AcyclicOrientedGraph *G, 
			       Population *population,
			       RecombinationParams *params);

int hboaRecomputeRemovalGains(int from,
			      int to,
			      int numVars,
			      GroupInformation *groupInformation,
			      OperatorGain *operatorGain,
			      OperatorApplicability *operatorApplicability,
			      AcyclicOrientedGraph *G, 
			      Population *population,
			      RecombinationParams *params);

int hboaRecomputeReversalGains(int from,
			       int to,
			       int numVars,
			       GroupInformation *groupInformation,
			       OperatorGain *operatorGain,
			       OperatorApplicability *operatorApplicability,
			       AcyclicOrientedGraph *G, 
			       Population *population,
			       RecombinationParams *params);

int hboaRecomputeJointGains(int from,
			    int to,
			    int numVars,
			    GroupInformation *groupInformation,
			    OperatorGain *operatorGain,
			    OperatorApplicability *operatorApplicability,
			    AcyclicOrientedGraph *G, 
			    Population *population,
			    RecombinationParams *params);

int hboaRecomputeAdditionGain(int from,
			      int to,
			      int numVars,
			      GroupInformation *groupInformation,
			      OperatorGain *operatorGain,
			      OperatorApplicability *operatorApplicability,
			      AcyclicOrientedGraph *G, 
			      Population *population,
			      RecombinationParams *params);

int hboaRecomputeRemovalGain(int from,
			     int to,
			     int numVars,
			     GroupInformation *groupInformation,
			     OperatorGain *operatorGain,
			     OperatorApplicability *operatorApplicability,
			     AcyclicOrientedGraph *G, 
			     Population *population,
			     RecombinationParams *params);

int hboaRecomputeReversalGain(int from,
			      int to,
			      int numVars,
			      GroupInformation *groupInformation,
			      OperatorGain *operatorGain,
			      OperatorApplicability *operatorApplicability,
			      AcyclicOrientedGraph *G, 
			      Population *population,
			      RecombinationParams *params);

// ---------------------------------------------------------------------------------------

int hboaUpdateAfterAddition(int newFrom, 
			    int newTo, 	
			    GroupInformation *groupInformation,
			    OperatorApplicability *operatorApplicability,
			    AcyclicOrientedGraph *G, 
			    RecombinationParams *params);

int hboaUpdateAfterRemoval(int newFrom, 
			   int newTo, 
			   GroupInformation *groupInformation,
			   OperatorApplicability *operatorApplicability,
			   AcyclicOrientedGraph *G, 
			   RecombinationParams *params);

int hboaUpdateAfterReversal(int newFrom, 
			    int newTo, 
			    GroupInformation *groupInformation,
			    OperatorApplicability *operatorApplicability,
			    AcyclicOrientedGraph *G, 
			    RecombinationParams *params);

int hboaUpdateAfterJoint(int newFrom, 
			 int newTo, 
			 GroupInformation *groupInformation,
			 OperatorApplicability *operatorApplicability,
			 AcyclicOrientedGraph *G, 
			 RecombinationParams *params);

int joinGroups( int a, 
		int b,   
		GroupInformation *groupInformation,
		OperatorApplicability *operatorApplicability,
		AcyclicOrientedGraph *G, 
		AcyclicOrientedGraph *newG,
		Population *P);

int hboaRecomputeNodeContribution(int i, GroupInformation *groupInformation, AcyclicOrientedGraph *G, Population *parents);

double hboaGetPenalty(int from, int to);

int joinGroupIndexes(GroupInformation *groupInformation, int *index, int n, int *destination);
//int joinIndexes(int *index1, int n1, int *index2, int n2, int *destination);
int indexedMatch(char *x, char *y, int *index, int n);
int copyIndexed(char *dest, char *src, int *index, int n);

#endif
