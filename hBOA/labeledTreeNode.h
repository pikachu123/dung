#ifndef _labeledTreeNode_h_
#define _labeledTreeNode_h_

#include <stdio.h>

typedef double Value;

// -------------------------------------
// well must polish this at some point

class LabeledTreeNode {

 public:
  
  int       which;
  
  int       depth;
  int       label;
  Value     value[2];
  
  void *leavesItem;
  void *leafParentsItem;
  
  char *parentLabelCoincidenceVector;
  
  Value *dArrayTmp;
  Value *rightValue0Array;
  Value *rightValue1Array;
  Value *leftValue0Array;
  Value *leftValue1Array;
  
  int             numParents;
  LabeledTreeNode **parent;
  
  LabeledTreeNode *left;
  LabeledTreeNode *right;
  
  LabeledTreeNode(int type);
  ~LabeledTreeNode();
  
  int freeTemporaryArrays();
  int allocateInstanceIndexEntries(int n);
  int reallocateInstanceIndexEntries(int n);
};

#endif
