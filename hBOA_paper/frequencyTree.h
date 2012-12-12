#ifndef _frequencyTree_h_
#define _frequencyTree_h_

#include "population.h"

// ----------------------------------

struct TreeNode {

  int       which;

  double    value;
  
  TreeNode *left;
  TreeNode *right;
};

// ----------------------------------

class FrequencyTree {

 private:

  TreeNode  *root;
  long       numLeaves;
  int        n;
  int        listSize;

  // some temporary variables not to waste memory

  long       globalK;
  double     globalDivisor;

 public:

  FrequencyTree();
  FrequencyTree(const FrequencyTree& frequencyTree);
  FrequencyTree(const FrequencyTree& frequencyTree, int *index, int indexLength);
  ~FrequencyTree();

  int  reset();
  int  computeFrequencies(Population *p);
  int  computeIndexedFrequencies(Population *p, int *index, int indexLength);
  int  computeIndexedFrequenciesWithSplit(Population *p, int *index, int indexLength, int **list, int *listLength, int listSize);
  
  long getNumInstances();
  long getNumInstancesWithSplit(int l);
  long getInstances(char **x);
  double getFrequency(char *x, int l);
  long getFrequencyList(double *x);
  long getInstancesAndFrequencies(char **x, double *f);
  long getInstancesAndFrequenciesWithSplit(char **x, double *f, int l);
  long getShiftedInstancesAndFrequencies(char **x, double *f, int shift);

 private:

  TreeNode *copySubtree(TreeNode *node);
  int       deleteSubtree(TreeNode *node);

  TreeNode *addInstance(char *x);
  TreeNode *addIndexedInstance(char *x, int *index);
  TreeNode *addIndexedInstanceWithSplit(char *x, int *index, int **list, int *listLength);

  int       divideNodes(double x);
  int       divideSubtreeNodes(TreeNode *node);

  long      getSubtreeNumInstances(TreeNode *node);
  long      getSubtreeNumInstancesWithSplit(TreeNode *node, int l);
  long      getSubtreeInstances(TreeNode *node, char **x, int depth, long i);
  long      getSubtreeFrequencyList(TreeNode *node, double *x, int depth, long i);
  long      getSubtreeInstancesAndFrequencies(TreeNode *node, char **x, double *f, int depth, long i);
  long      getSubtreeInstancesAndFrequenciesWithSplit(TreeNode *node, char **x, double *f, int l, int depth, long i);

  int       reduceByIndex(TreeNode *source, int *index, int indexLength);
  int       reduceSubtreeByIndex(TreeNode **currentNode, TreeNode *sourceTree, int *index, int indexLength, int depth);
  int       addTree(TreeNode **destination, TreeNode *source);

  long      recomputeNumLeaves(TreeNode *node);
};

// ----------------------------------

struct TreeSplit:TreeNode {
  int             which;

  int             numSubtrees;
  FrequencyTree **subtree;
};


#endif
