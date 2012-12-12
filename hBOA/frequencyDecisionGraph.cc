#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "decisionGraph.h"
#include "frequencyDecisionGraph.h"

// --------------------------------------------------------------------

FrequencyDecisionGraph::FrequencyDecisionGraph(Population *p):DecisionGraph()
{
  // set the active population and dimension

  this->p          = p;
  this->myPosition = -1;

  // allocate the instance index stuff (we don't do this in traditional decision trees)

  getRoot()->allocateInstanceIndexEntries(p->n);

  // the probability in the root is 1 (no position associated with the tree)

  setValue(getRoot(),1,0);

  // myNode is for nothing until we have some position

  myNode = NULL;
};

// --------------------------------------------------------------------

FrequencyDecisionGraph::FrequencyDecisionGraph(Population *p, int myPosition):DecisionGraph()
{
  double f;

  // set the active population and dimension

  this->p          = p;
  this->myPosition = myPosition;

  // allocate the instance index stuff (we don't do this in traditional decision trees)

  getRoot()->allocateInstanceIndexEntries(p->n);

  // the probability in the root is the univariate frequency in the position

  f = univariateFrequency(myPosition);
  setValue(getRoot(),1-f,0);
  setValue(getRoot(),f,1);

  // everyone is satisfied by the empty tree (which has only root node)

  myNode = (LabeledTreeNode**) calloc(p->N,sizeof(LabeledTreeNode*));
  for (int i=0; i<p->N; i++)
    myNode[i] = getRoot();
};

// --------------------------------------------------------------------

FrequencyDecisionGraph::~FrequencyDecisionGraph()
{
  if (myNode!=NULL)
    free(myNode);
  myNode=NULL;
};

// --------------------------------------------------------------------
// merge the children of the node (with update of the frequency values)

int FrequencyDecisionGraph::merge(LabeledTreeNode *x, LabeledTreeNode *y)
{
  // increase the value in x by the value in y

  x->value[0] += y->value[0];
  if (myPosition>=0)
    x->value[1] += y->value[1];

  // must check the myNode entries and merge them

  for (int i=0; i<p->N; i++)
    if (myNode[i]==y)
      myNode[i]=x;

  // merge the children of x now
 
  DecisionGraph::merge(x,y);
  
  // get back
  
  return 0;
};

// --------------------------------------------------------------------
// merge the children of the node (with update of the frequency values)

int FrequencyDecisionGraph::split(LabeledTreeNode *x, int label)
{
  // split the node x now

  DecisionGraph::split(x,label);

  // allocate instance/frequency related stuff for both children

  x->left->allocateInstanceIndexEntries(p->n);
  x->right->allocateInstanceIndexEntries(p->n);
 
  // it is sufficient to compute frequencies for one of the two new instances,
  // because the second one can be computed by using the value in x (which is
  // the sum of the two), so that's what we're gonna do

  // copy the instance information from the parent to the new children (if any)
  
  if (x->depth>0)
    {
      memcpy(x->left->parentLabelCoincidenceVector,x->parentLabelCoincidenceVector,p->n);
      memcpy(x->right->parentLabelCoincidenceVector,x->parentLabelCoincidenceVector,p->n);
    };
  
  // add the new information to each child
  
  x->left->parentLabelCoincidenceVector[x->label]  = 1;
  x->right->parentLabelCoincidenceVector[x->label] = 1;
  
  // should recompute frequencies?!!!
  
  // get back
  
  return 0;
};

// ------------------------------------------------------

double FrequencyDecisionGraph::univariateFrequency(int k)
{
  long i;
  long N;
  long count;

  // assign helper variables

  N = p->N;

  // compute the frequency of x_k=1

  count=0;
  for (i=0; i<N; i++)
    count += p->individual[i].chromosome[k];

  // return the resulting frequency

  return double(count)/N;
};

// ------------------------------------------------------

int FrequencyDecisionGraph::computeFrequencies()
{
  long i;
  int  j;
  int numLeaves;
  long N;
  LabeledTreeNode *node;
  char *s;
  NodeListItem *leafItem;

  // set helper variables

  N         = p->N;
  numLeaves = getNumLeaves();

  // set the value in all leaves to 0

  leafItem = getLeaves();
  for (j=0; j<numLeaves; j++, leafItem=leafItem->next)
    leafItem->x->value[0]=leafItem->x->value[1]=0;
  
  // compute the new values

  for(i=0; i<N; i++)
    {
      s=p->individual[i].chromosome;
      //      node = getRoot();
      node = myNode[i];

      if (node->which!=LEAF)
	{
	  while (node->which!=LEAF)
	    {
	      if (s[node->label]==0)
		node = node->left;
	      else
		node = node->right;
	    };
	  myNode[i] = node;
	};
      
      if (myPosition>=0)
	node->value[s[myPosition]]++;
      else
	node->value[0]++;
    };
  
  // divide the values in all leaves by N

  leafItem = getLeaves();
  if (myPosition>=0)
    for (j=0; j<numLeaves; j++, leafItem=leafItem->next)
      {
	leafItem->x->value[0]/=N;
	leafItem->x->value[1]/=N;
      }
  else
    for (j=0; j<numLeaves; j++, leafItem=leafItem->next)
      leafItem->x->value[0]/=N;
    
  // get back
 
  return 0;
};

// ------------------------------------------------------

int FrequencyDecisionGraph::computeBoltzmannFrequencies(double beta, double coef)
{
  long i;
  int  j;
  int numLeaves;
  long N;
  LabeledTreeNode *node;
  char *s;
  NodeListItem *leafItem;
  double total;

  printf("Computing boltzmann frequencies\n");

  // set helper variables

  N         = p->N;
  numLeaves = getNumLeaves();

  // set the value in all leaves to 0

  leafItem = getLeaves();
  for (j=0; j<numLeaves; j++, leafItem=leafItem->next)
    leafItem->x->value[0]=leafItem->x->value[1]=0;
  
  // compute the new values

  total=0;
  for(i=0; i<N; i++)
    {
      s=p->individual[i].chromosome;
      node = getRoot();

      while (node->which!=LEAF)
	{
	  if (s[node->label]==0)
	    node = node->left;
	  else
	    node = node->right;
	};
      
      if (myPosition>=0)
	{
	  node->value[s[myPosition]]+=exp(beta*p->individual[i].f/coef);
	  total += node->value[s[myPosition]];
	}
      else
	{
	  node->value[0]+=exp(beta*p->individual[i].f/coef);
	  total += node->value[0];
	};
    };


  // divide the values in all leaves by N

  leafItem = getLeaves();
  if (myPosition>=0)
    for (j=0; j<numLeaves; j++, leafItem=leafItem->next)
      {
	leafItem->x->value[0]/=total;
	leafItem->x->value[1]/=total;
      }
  else
    for (j=0; j<numLeaves; j++, leafItem=leafItem->next)
      {
	leafItem->x->value[0]/=total;
      };

  // get back
 
  return 0;
};

// ------------------------------------------------------

int FrequencyDecisionGraph::computeSplitFrequencies(LabeledTreeNode *x, Value *left0, Value *left1, Value *right0, Value *right1)
{
  long i;
  int  n;
  long N;
  LabeledTreeNode *node;
  char *s;
  int label;

  // set helper variables

  N         = p->N;
  n         = p->n;

  // reset the counts

  for (label=0; label<n; label++)
    left0[label] = left1[label] = right0[label] = right1[label] = 0;
  
  // compute the frequencies
  
  for(i=0; i<N; i++)
    {
      s=p->individual[i].chromosome;
      //      node = rootNode;

      //      while (node->which!=LEAF)
      //	{
      //	  if (s[node->label]==0)
      //	    node = node->left;
      //	  else
      //	    node = node->right;
      //	};

      node = myNode[i];
      if (node->which!=LEAF)
	{
	  if (s[node->label])
	    node = node->left;
	  else
	    node = node->right;
	  myNode[i]=node;
	};
      
      if (node==x)
	if (myPosition>=0)
	  {
	    for (label=0; label<n; label++)
	      if (s[label]==0)
		if (s[myPosition])
		  left1[label]++;
		else
		  left0[label]++;
	      else
		if (s[myPosition])
		  right1[label]++;
		else
		  right0[label]++;
	  }
	else
	  {
	    for (label=0; label<n; label++)
	      if (s[label]==0)
		left0[label]++;
	      else
		right0[label]++;
	  };
    };

  // divide the values in all leaves by N

  if (myPosition>=0)
    for (label=0; label<n; label++)
      {
	left0[label]  /= N;
	left1[label]  /= N;
	right0[label] /= N;
	right1[label] /= N;
      }
  else
    for (label=0; label<n; label++)
      {
	left0[label]  /= N;
	right0[label] /= N;
      };

  // get back
 
  return 0;
};
