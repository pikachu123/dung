#ifndef _decisionGraphBoa_h_
#define _decisionGraphBoa_h_

#include "population.h"
#include "frequencyDecisionGraph.h"
#include "recombination.h"

struct Operator {
  float gain;
  int   type;
  int   where;
  int   label;
  FrequencyDecisionGraph *t;
  LabeledTreeNode *node;
  LabeledTreeNode *node2;
};

struct MergeOperator {
  float gain;
  LabeledTreeNode *a;
  LabeledTreeNode *b;
};

int boaWithDecisionGraphsRecombination(Population *parents, Population *children, long M, RecombinationParams *params);

int learnDecisionGraphBayesianNetwork(AcyclicOrientedGraph *G,
				     FrequencyDecisionGraph **T,
				     Population *P,
				     RecombinationParams *params);

int simulateDecisionGraphBayesianNetwork(AcyclicOrientedGraph *G,
					FrequencyDecisionGraph **T,
					Population *P);

int recomputeDecisionGraphSplitGains(FrequencyDecisionGraph *t,
				     LabeledTreeNode *x, 
				     AcyclicOrientedGraph *G,
				     int maxParents,
				     int node,
				     int n,
				     long N,
				     RecombinationParams *params);

int recomputeDecisionGraphMergeGains(FrequencyDecisionGraph *t,
				     AcyclicOrientedGraph *G,
				     MergeOperator **merge,
				     int *numMerges,
				     int node,
				     int n,
				     long N);

double nodeContribution(LabeledTreeNode *x, long N);

int resetOperator(Operator *x);
int updateBestOperator(Operator *best, Operator *x);
int updateBestNodeOperator(Operator *x, FrequencyDecisionGraph *t, MergeOperator *merge, int numMerges, int node, int n);
int applyOperator(Operator *x, AcyclicOrientedGraph *G);
int updateGainsAfterOperator(Operator *x, AcyclicOrientedGraph *G, int maxParents, long N, RecombinationParams *params);
int topologicalOrdering(AcyclicOrientedGraph *G, int *index);
int generateInstance(char *x,
		     AcyclicOrientedGraph *G,
		     int *index,
		     FrequencyDecisionGraph **T);

int operatorApplicable(Operator *x, AcyclicOrientedGraph *G);
int deleteOperator(FrequencyDecisionGraph **T, Operator *x);

int recomputeBoltzmannFrequencies(FrequencyDecisionGraph **T,
				  Population *P,
				  RecombinationParams *params);

char *dBOAGetMetricDescription(int n);

#endif
