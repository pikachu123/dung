#ifndef _decisionGraph_h_
#define _decisionGraph_h_

#include <stdio.h>

#include "labeledTreeNode.h"

#define LEAF  0
#define SPLIT 1

// -------------------------------------

struct NodeListItem {
  LabeledTreeNode *x;

  NodeListItem *previous;
  NodeListItem *next;
};

// -------------------------------------

class DecisionGraph {

 private:
  LabeledTreeNode *root;

  int              numLeaves;
  NodeListItem    *leaves;

  LabeledTreeNode *iterator;
  NodeListItem    *leafIterator;

  int deleteSubtree(LabeledTreeNode *x);
  int deleteNodeList(NodeListItem *x);

  int updateLeavesAfterSplit(LabeledTreeNode *x);
  int updateLeavesBeforeMerge(LabeledTreeNode *x, LabeledTreeNode *y);

  int recursivePrint(FILE *out, LabeledTreeNode *x, int shift);

 protected:
  LabeledTreeNode *newNode();
  int deleteNode(LabeledTreeNode *x);

 public:
  DecisionGraph();
  ~DecisionGraph();

  // get the root

  LabeledTreeNode *getRoot();
  Value           setValue(LabeledTreeNode *x, Value value, char which);
  Value           getValue(LabeledTreeNode *x, char which);

  // node operations

  int split(LabeledTreeNode *x, int label);
  int merge(LabeledTreeNode *x, LabeledTreeNode *y);

  // functions on the leaves

  int           getNumLeaves();
  NodeListItem *getLeaves();

  // iterator functions

  LabeledTreeNode *getIterator();
  Value            getIteratorValue(char which);
  int              getIteratorLabel();
  int              iteratorGoToRoot();
  int              iteratorGoTo(LabeledTreeNode *x);
  int              iteratorGoLeft();
  int              iteratorGoRight();
  int              iteratorFollowInstance(char *x);
  int              iteratorFollowInstanceFromRoot(char *x);

  // leaf iterator functions

  NodeListItem    *getLeafIterator();
  double           getLeafIteratorValue(char which);
  LabeledTreeNode *getLeafIteratorNode();
  int              resetLeafIterator();
  int              leafIteratorNext();
  int              leafIteratorPrevious();

  // print out the tree, starting from the root

  int print(FILE *out);
};
  
#endif
