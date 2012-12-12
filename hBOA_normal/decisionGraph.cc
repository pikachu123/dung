#include <stdlib.h>

#include "decisionGraph.h"
#include "memalloc.h"
#include "random.h"

#define SHIFT_STEP 3

// -----------------------------------------
// constructor

DecisionGraph::DecisionGraph()
{
  // create the root node

  root = new LabeledTreeNode(LEAF);
  
  // one leaf at the moment (the root)

  numLeaves = 1;

  // add the root to the list of leaf nodes

  leaves           = new NodeListItem;
  leaves->x        = root;
  leaves->previous = leaves->next = NULL;

  // add pointer to the root's item in the leaf list to the root

  root->leavesItem = leaves;

  // root has no parents

  root->numParents = 0;
  root->parent     = NULL;

  // reset the iterators

  iteratorGoToRoot();
  resetLeafIterator();
};

// -----------------------------------------
// destructor

DecisionGraph::~DecisionGraph()
{
  deleteSubtree(root);
  deleteNodeList(leaves);
};

// -----------------------------------------
// deletes the subtree from memory

int DecisionGraph::deleteSubtree(LabeledTreeNode *x)
{
  int i;

  // is the tree NULL?

  if (x==NULL)
    return 0;

  // delete the subtrees (careful about identical children)

  if (x->left!=x->right)
    deleteSubtree(x->left);
  deleteSubtree(x->right);

  // cut the node from the parents (not to be deleted multiple times)

  for (i=0; i<x->numParents; i++)
    {
      if (x->parent[i]->left==x)
	x->parent[i]->left=NULL;
      if (x->parent[i]->right==x)
	x->parent[i]->right=NULL;
    }

  // delete the list of the parents of the leaf

  if (x->parent)
    Free(x->parent);

  // delete the node itself

  delete x;

  // get back

  return 0;
};

// -----------------------------------------
// deletes the leaf list from memory

int DecisionGraph::deleteNodeList(NodeListItem *x)
{
  if (x==NULL)
    return 0;

  deleteNodeList(x->next);
  delete x;

  return 0;
};

// -----------------------------------------

LabeledTreeNode *DecisionGraph::getRoot()
{
  // return the root

  return root;
};

// -----------------------------------------

Value DecisionGraph::setValue(LabeledTreeNode *x, Value value, char which)
{
  // sets the value of a node and returns it

  return x->value[which]=value;
};

// -----------------------------------------

Value DecisionGraph::getValue(LabeledTreeNode *x, char which)
{
  // return the value of a node

  return x->value[which];
};

// -----------------------------------------
// splits a node

int DecisionGraph::split(LabeledTreeNode *x, int label)
{
  // assign the label and type

  x->label = label;
  x->which = SPLIT;

  // split the node

  x->left  = new LabeledTreeNode(LEAF);
  x->right = new LabeledTreeNode(LEAF);

  // allocate a parent for each of the children

  x->left->numParents=1;
  x->left->parent = (LabeledTreeNode**) malloc(sizeof(LabeledTreeNode*));

  x->right->numParents=1;
  x->right->parent = (LabeledTreeNode**) malloc(sizeof(LabeledTreeNode*));

  // assign the parent
  
  x->left->parent[0] = x->right->parent[0] = x;
  
  // assign the depths

  x->left->depth = x->right->depth = x->depth+1;
  
  // replace the split node by its children
  
  updateLeavesAfterSplit(x);

  // get back

  return 0;
};

// -----------------------------------------
//XXXXXXXXXXXXXXXXX

int DecisionGraph::merge(LabeledTreeNode *x, LabeledTreeNode *y)
{
  int i,j,k;
  int newParents;
  int found;

  // update the leaf list

  updateLeavesBeforeMerge(x,y);

  // redirect the y's parents to x instead of y

  for (i=0; i<y->numParents; i++)
    if (y->parent[i]->left==y)
      y->parent[i]->left=x;
    else
      y->parent[i]->right=x;
  
  // assign all parents from y to x (only those that are new to x)

  newParents=0;
  for (i=0; i<y->numParents; i++)
    {
      found=0;
      for (j=0; j<x->numParents; j++)
	if (x->parent[j]==y->parent[i])
	  found=1;

      if (!found)
	newParents++;
    };

  x->parent=(LabeledTreeNode**) realloc(x->parent,(x->numParents+newParents)*sizeof(LabeledTreeNode*));

  k=x->numParents;
  for (i=0; i<y->numParents; i++)
    {
      found=0;
      for (j=0; j<x->numParents; j++)
	if (x->parent[j]==y->parent[i])
	  found=1;

      if (!found)
	x->parent[k++] = y->parent[i];
    };

  x->numParents+=newParents;

//   printf("Mergin nodes x and y, adding %u new parents.\n",newParents);

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::updateLeavesAfterSplit(LabeledTreeNode *x)
{
  NodeListItem *item;
  NodeListItem *newItem;

  // where is the leaf located in the list?

  item = (NodeListItem*) x->leavesItem;

  // replace the original item for x by an item for its left child
  
  item->x = x->left;
  x->left->leavesItem = item;

  // create a new list item and insert it between item and item->next

  newItem = new NodeListItem;
 
  if (item->next)
    item->next->previous = newItem;     // redirect the right neighbour

  newItem->next     = item->next;       // redirect the new item
  newItem->previous = item;             // 
  item->next        = newItem;          // redirect the left neighbour

  // new item points to the right child

  newItem->x = x->right;
  x->right->leavesItem = newItem;

  // increase the number of leaves by 1 (we replaced one by two)

  numLeaves++;

  // x's got no items in the list anymore

  x->leavesItem = NULL;

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::updateLeavesBeforeMerge(LabeledTreeNode *x, LabeledTreeNode *y)
{
  NodeListItem *xItem, *yItem;

  // assign x to the list item of the left child

  xItem = (NodeListItem*) x->leavesItem;
  yItem = (NodeListItem*) y->leavesItem;

  x->leavesItem = xItem;
  xItem->x = x;
  
  // delete the item for the right child

  if (yItem->previous)
    yItem->previous->next=yItem->next;

  if (yItem->next)
    yItem->next->previous=yItem->previous;

  if (yItem==leaves)
    leaves=yItem->next;

  Free(yItem);

  // decrement the number of leaves

  numLeaves--;

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::getNumLeaves()
{
  // return the number of leaves

  return numLeaves;
};

// -----------------------------------------

NodeListItem *DecisionGraph::getLeaves()
{
  // return the list of leaves

  return leaves;
};

// -----------------------------------------

LabeledTreeNode *DecisionGraph::getIterator()
{
  // returns the iterator

  return iterator;
};

// -----------------------------------------

Value DecisionGraph::getIteratorValue(char which)
{
  // returns the iterator

  return iterator->value[which];
};

// -----------------------------------------

int DecisionGraph::getIteratorLabel()
{
  // returns the iterator

  return iterator->label;
};

// -----------------------------------------

int DecisionGraph::iteratorGoToRoot()
{
  // move the iterator to the root

  iteratorGoTo(root);

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::iteratorGoTo(LabeledTreeNode *x)
{
  // move the iterator to the node x

  iterator = x;

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::iteratorGoLeft()
{
  // move the iterator to the left child, if possible

  iterator = iterator->left;

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::iteratorGoRight()
{
  // move the iterator to the right child, if possible

  iterator = iterator->right;

  // get back

  return 0;
};

// ------------------------------------------------------

int DecisionGraph::iteratorFollowInstance(char *x)
{
  while (iterator->which!=LEAF)
    {
      if (x[iterator->label]==0)
	iteratorGoLeft();
      else
	iteratorGoRight();
    };

  return 0;
};

// ------------------------------------------------------

int DecisionGraph::iteratorFollowInstanceFromRoot(char *x)
{
  iterator=root;

  while (iterator->which!=LEAF)
    {
      if (x[iterator->label]==0)
	iteratorGoLeft();
      else
	iteratorGoRight();
    };

  return 0;
};

// -----------------------------------------

int DecisionGraph::print(FILE *out)
{
  // call the recursive print function

  recursivePrint(out,root,0);

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::recursivePrint(FILE *out, LabeledTreeNode *x, int shift)
{
  int i;

  // if the tree is NULL, do nothing

  if (x==NULL)
    return 0;

  // print shift spaces

  for (i=0; i<shift; i++)
    fputc(' ',out);

  // for leaves, print the value, for splits, print the label

  if (x->which==LEAF)
    fprintf(out,"%1.2f  %1.2f (%u)\n",x->value[0],x->value[1],x->numParents);
  if (x->which==SPLIT)
    fprintf(out,"%i (%u)\n",x->label,x->numParents);

  // print the children sub-trees

  recursivePrint(out,x->left,shift+SHIFT_STEP);
  recursivePrint(out,x->right,shift+SHIFT_STEP);

  // get back

  return 0;
};

// -----------------------------------------

NodeListItem *DecisionGraph::getLeafIterator()
{
  // returns the leaf iterator

  return leafIterator;
};

// -----------------------------------------

int DecisionGraph::resetLeafIterator()
{
  // set the leaf iterator to the first leaf

  leafIterator = leaves;

  // get back

  return 0;
};

// -----------------------------------------

int DecisionGraph::leafIteratorNext()
{
  // can't do anything with NULL, error

  if (leafIterator==NULL)
    return -1;

  // move the leaf iterator to the next leaf

  leafIterator = leafIterator->next;
  
  // get back
  
  return 0;
};

// -----------------------------------------

int DecisionGraph::leafIteratorPrevious()
{
  // can't do anything with NULL, error

  if (leafIterator==NULL)
    return -1;

  // move the leaf iterator to the next leaf

  leafIterator = leafIterator->previous;
  
  // get back
  
  return 0;
};

// -----------------------------------------

double DecisionGraph::getLeafIteratorValue(char which)
{
  // return the value of the leaf corresponding to a current value of the leaf iterator

  return leafIterator->x->value[which];
};

// -----------------------------------------

LabeledTreeNode *DecisionGraph::getLeafIteratorNode()
{
  // return the leaf corresponding to a current value of the leaf iterator

  return leafIterator->x;
};
