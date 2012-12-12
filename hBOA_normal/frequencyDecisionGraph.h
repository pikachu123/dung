#ifndef _frequencyDecisionGraph_h_
#define _frequencyDecisionGraph_h_

#include "decisionGraph.h"
#include "population.h"

// -------------------------------------

class FrequencyDecisionGraph: public DecisionGraph
{
 private:

  Population *p;
  int         myPosition;

  double univariateFrequency(int k);
  double instanceFrequency(int *index, char *x, int n);

  LabeledTreeNode** myNode;
  
 public:
  FrequencyDecisionGraph(Population *p);
  FrequencyDecisionGraph(Population *p, int myPosition);
  ~FrequencyDecisionGraph();
  
  int merge(LabeledTreeNode *x, LabeledTreeNode *y);
  int split(LabeledTreeNode *x, int label);

  int computeFrequencies();
  int computeBoltzmannFrequencies(double beta, double coef);
  int computeSplitFrequencies(LabeledTreeNode *x, Value *left0, Value *left1, Value *right0, Value *right1);
};

#endif
