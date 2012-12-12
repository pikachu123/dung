#include <stdlib.h>

#include "labeledTreeNode.h"

LabeledTreeNode::LabeledTreeNode(int type)
{
  // assign the variables
  
  which=type; 
  
  left=NULL; 
  right=NULL; 
  
  label=-1; 
  value[0]=0; 
  value[1]=0;

  numParents=0;
  parent=NULL; 
  depth=0; 
  
  leavesItem=NULL;
  leafParentsItem=NULL;

  parentLabelCoincidenceVector=NULL;

  dArrayTmp = NULL;
  rightValue0Array = NULL;
  rightValue1Array = NULL;
  leftValue0Array  = NULL;
  leftValue1Array  = NULL;
};

LabeledTreeNode::~LabeledTreeNode()
{

  if (parentLabelCoincidenceVector)
    free(parentLabelCoincidenceVector);
  if (dArrayTmp)
    free(dArrayTmp);
  if (rightValue0Array)
    free(rightValue0Array);
  if (rightValue1Array)
    free(rightValue1Array);
  if (leftValue0Array)
    free(leftValue0Array);
  if (leftValue1Array)
    free(leftValue1Array);
};

int LabeledTreeNode::freeTemporaryArrays()
{
  if (dArrayTmp)
    free(dArrayTmp);
  if (rightValue0Array)
    free(rightValue0Array);
  if (rightValue1Array)
    free(rightValue1Array);
  if (leftValue0Array)
    free(leftValue0Array);
  if (leftValue1Array)
    free(leftValue1Array);

  dArrayTmp = NULL;
  rightValue0Array = NULL;
  rightValue1Array = NULL;
  leftValue0Array  = NULL;
  leftValue1Array  = NULL;

  return 0;
};

int LabeledTreeNode::allocateInstanceIndexEntries(int n)
{
  parentLabelCoincidenceVector = (char*) calloc(n,sizeof(char));

  // get back

  return 0;
};
