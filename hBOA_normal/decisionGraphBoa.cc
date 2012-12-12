#include <stdio.h>
#include <string.h>

#include "decisionGraphBoa.h"

#include "mymath.h"
#include "memalloc.h"
#include "individual.h"
#include "population.h"
#include "recombination.h"
#include "graph.h"
#include "binary.h"
#include "utils.h"
#include "BDe.h"
//#include "mda.h"
#include "mdl.h"
#include "priors.h"
#include "random.h"
#include "counts.h"
#include "frequencyDecisionGraph.h"

#define YES 1
#define NO  0

#define OPERATOR_NONE      -1
#define OPERATOR_SPLIT_NODE 0
#define OPERATOR_MERGE_NODE 1

#define FIXED_THRESHOLD 0

#define EPSILON 20

//#define PRINT_SPLITS
//#define PRINT_GROUPS
//#define PRINT_OPERATORS
//#define PRINT_MAX_DEPTH

//#define STORE_SPLIT_FREQUENCIES

long numAppliedSplits,numAppliedMerges;
int  maxDepth;
int metricToUse;
int penalty;

#define BAYESIAN     0
#define MDL          1

char *dBOAMetricDesc[3] = { "K2 metric with penalty",
			    "K2 metric without penalty",
			    "MDL metric" };

int boaWithDecisionGraphsRecombination(Population *parents, Population *children, long M, RecombinationParams *params)
{
  int k;
  int  n;
  AcyclicOrientedGraph   *G;
  FrequencyDecisionGraph **T;

  switch (params->dBOAMetricN) {
  case 0:
    metricToUse = BAYESIAN;
    penalty     = 1; 
    break;
    
  case 1:
   metricToUse = BAYESIAN;
   penalty     = 0; 
   break;

  case 2:
    metricToUse = MDL;
    penalty     = 1;
    break;

  default:
    fprintf(stderr,"ERROR: No such metric in dBOA...exiting.\n");
    exit(-1);
  };
 
  // set helper variables

  n = parents->n;
 
  // allocate the children population (just do it before we do the rest)

  allocatePopulation(children,M,parents->numDiscrete,parents->numContinuous);

  // allocate the memory for frequency trees
  
  T = (FrequencyDecisionGraph**) Calloc(n,sizeof(FrequencyDecisionGraph*));

  // create a new graph
  
  G = new AcyclicOrientedGraph(n);
  
      // reset the graph

      G->removeAllEdges();

      // create a frequency decision tree for each node in the graph

      for (k=0; k<n; k++)
	T[k] = new FrequencyDecisionGraph(parents,k);

      // learn the Bayesian network for the population of parents

      learnDecisionGraphBayesianNetwork(G,T,parents,params);

      // temporary!!!!!!!!!!!

      if (params->useBoltzmannFrequencies)
	recomputeBoltzmannFrequencies(T,parents,params);

      // and use it to generate new points

      simulateDecisionGraphBayesianNetwork(G,T,children);
      
      // delete the decision trees

      for (k=0; k<n; k++)
	delete T[k];

  // put the resulting graph in the parameters (might be needed later)

  if (params->G)
    delete params->G;
  
  params->G = G;
  
  // free the array for the trees
  
  Free(T);

  // the created offspring has not been evaluated yet

  children->evaluated = 0;
  
  // get back

  return 0;
};

// ====================================================================

int recomputeBoltzmannFrequencies(FrequencyDecisionGraph **T,
				  Population *P,
				  RecombinationParams *params)
{
  double variance;
  double max;
 
  long i,N;

  N = P->N;

  variance = fitnessVariance(P);

  max=0;
  for (i=0; i<N; i++)
    if (P->individual[i].f>max)
      max=P->individual[i].f;
  
  if (params->recombinationBoltzmannBeta*max/variance>23)
    variance=(max*params->recombinationBoltzmannBeta)/23;

  for (i=0; i<P->n; i++)
    T[i]->computeBoltzmannFrequencies(params->recombinationBoltzmannBeta,variance);

  return 0;
};

// ====================================================================

int learnDecisionGraphBayesianNetwork(AcyclicOrientedGraph *G,
				     FrequencyDecisionGraph **T,
				     Population *P,
				     RecombinationParams *params)
{
  int    i,j;
  int    n;
  long   N;
  int    maxIncoming;
  MergeOperator **merge;
  int           *numMerges;

  char **canAdd;

  long     operatorsLeft;
  Operator bestOperator;
  Operator *bestNodeOperator;

  // no operators applied

  numAppliedSplits=numAppliedMerges=0;
  maxDepth=0;

  // set the helper variables

  N           = P->N;
  n           = P->n;
  maxIncoming = params->maxIncoming;

  // allocate all needed parameters

  canAdd = (char**) calloc(n,sizeof(char*));
  for (i=0; i<n; i++)
    canAdd[i] = (char*) malloc(n);

  // allocate the operators for each node

  bestNodeOperator = (Operator*) calloc(n,sizeof(Operator));
  for (i=0; i<n; i++)
    resetOperator(&(bestNodeOperator[i]));

  // well allocate some memory for the info on merge operators

  numMerges = (int*) calloc(n,sizeof(int));
  merge = (MergeOperator**) calloc(n,sizeof(MergeOperator*));
  for (i=0; i<n; i++)
    {
      merge[i]=NULL;
      numMerges[i]=0;
    };

  // remove all edges (just in case...)

  G->removeAllEdges();

  // ---- reset all needed parameters

  // all edges can be added (unless we are UMDA which we are not)
  
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      canAdd[i][j]=YES;


  // no best operator at the moment

  resetOperator(&bestOperator);

  // allocate and recompute the gains (not to worry about freeing them, it is done automatically)

  for (i=0; i<n; i++)
    {
      T[i]->getLeaves()->x->dArrayTmp=(Value*) calloc(n,sizeof(Value));
#ifdef STORE_SPLIT_FREQUENCIES
      T[i]->getLeaves()->x->rightValue0Array = (Value*) calloc(G->size(),sizeof(Value));
      T[i]->getLeaves()->x->rightValue1Array = (Value*) calloc(G->size(),sizeof(Value));
      T[i]->getLeaves()->x->leftValue0Array = (Value*) calloc(G->size(),sizeof(Value));
      T[i]->getLeaves()->x->leftValue1Array = (Value*) calloc(G->size(),sizeof(Value));
#endif
      recomputeDecisionGraphSplitGains(T[i],T[i]->getLeaves()->x,G,maxIncoming,i,n,N,params);
      updateBestNodeOperator(&(bestNodeOperator[i]),T[i],merge[i],numMerges[i],i,n);
      updateBestOperator(&bestOperator,&(bestNodeOperator[i]));
    };

  // all operators are left now, we can split in each node once on any variable

  operatorsLeft = n*n;

  // while we gained and can do some more do splits

  while ((bestOperator.gain>0)&&(operatorsLeft>0))
    {
      if (operatorApplicable(&bestOperator,G))
	{
	  // apply the operator

	  applyOperator(&bestOperator,G);
	 
	  // update the gains
 
	  updateGainsAfterOperator(&bestOperator,G,maxIncoming,N,params);
	  recomputeDecisionGraphMergeGains(T[bestOperator.where],
					   G,
					   merge,
					   &(numMerges[bestOperator.where]),
					   bestOperator.where,
					   n,
					   N);
	}
      else
	// delete the operator (we can't apply it for some reason)

	deleteOperator(T,&bestOperator);

      // update the best operator in the node we've touched

      updateBestNodeOperator(&(bestNodeOperator[bestOperator.where]),
			     T[bestOperator.where],
			     merge[bestOperator.where],
			     numMerges[bestOperator.where],
			     bestOperator.where,
			     n);

      // everything alright, now we must get the best gain again
      
      resetOperator(&bestOperator);
      
      for (i=0; i<n; i++)
	updateBestOperator(&bestOperator,&(bestNodeOperator[i]));
    };

  // free it all

  for (i=0; i<n; i++)
    free(canAdd[i]);
  free(canAdd);

  free(bestNodeOperator);

  free(numMerges);
  for (i=0; i<n; i++)
    if (merge[i]!=NULL)
      free(merge[i]);
  free(merge);

#ifdef PRINT_OPERATORS
  printf("Number of splits performed: %lu\n",numAppliedSplits);
  printf("Number of merges performed: %lu\n",numAppliedMerges);
#endif

#ifdef PRINT_MAX_DEPTH
  printf("Maximum depth: %u\n",maxDepth);
#endif

  // get back

  return 0;
};

// ====================================================================

int simulateDecisionGraphBayesianNetwork(AcyclicOrientedGraph *G,
					FrequencyDecisionGraph **T,
					Population *P)
{
  long N;
  int  n;
  int *index;

  // assign the helper variables
  
  N = P->N;
  n = P->n;

  // allocate index for topological ordering

  index = (int*) calloc(n,sizeof(int));

  // order the vertices topologically

  topologicalOrdering(G,index);

  // print the groups (just for debugging)

#ifdef PRINT_GROUPS

  printf("\n\n=======================================================\nGroups:\n");

  printf("Topological order: ");
  for (int i=0; i<n; i++)
    printf("%2u ",index[i]);
  printf("\n\n");
  for (int i=0; i<n; i++)
    {
      printf("%u -> ",index[i]);
      for (int j=0; j<G->getNumIn(index[i]); j++)
	printf("%2u ",G->getParentList(index[i])[j]);
      printf("\n");
    };
  getchar();

#endif

  // generate the new stuff

  for (long i=0; i<N; i++)
    generateInstance(P->individual[i].chromosome,G,index,T);

  // free memory

  free(index);

  // get back

  return 0;
};

// ====================================================================

int generateInstance(char *x,
		     AcyclicOrientedGraph *G,
		     int *index,
		     FrequencyDecisionGraph **T)
{
  int i;
  int n;
  int position;
  int value;
  double p0,p1;

  // assign the helper variables

  n = G->size();

  // generate the stuff

  for (i=0; i<n; i++)
    {
      position = index[i];
      T[position]->iteratorFollowInstanceFromRoot(x);
      
      p0 = T[position]->getIterator()->value[0];
      p1 = T[position]->getIterator()->value[1];
      p0 = p0/(p0+p1);
      p1 = 1-p0;

      value=(drand()<p1)? 1:0;

      x[position]=value;
    };

  // get back

  return 0;
};

// ====================================================================

int topologicalOrdering(AcyclicOrientedGraph *G,
			int *index)
{
  int i,j;
  int n;
  int numAdded;
  char *added;
  int numParents;
  int  *parentList;
  char canAdd;

  // assign the helper variables

  n = G->size();

  // allocate memory for the added array

  added = (char*) malloc(n);
  for (i=0; i<n; i++)
    added[i]=0;

  // order

  numAdded=0;
  while (numAdded<n)
    for (i=0; i<n; i++)
      if (!added[i])
	{
	  canAdd=1;
	  parentList=G->getParentList(i);
	  numParents=G->getNumIn(i);
	  for (j=0; (j<numParents)&&(canAdd); j++)
	    {
	      if (!added[parentList[j]])
		canAdd=0;
	    };
	  
	  if (canAdd)
	    {
	      index[numAdded]=i;
	      numAdded++;
	      added[i]=1;
	    };
	};
	
  // free memory

  free(added);

  // get back

  return 0;
};

// ====================================================================

int recomputeDecisionGraphSplitGains(FrequencyDecisionGraph *t,
				     LabeledTreeNode *x, 
				     AcyclicOrientedGraph *G,
				     int maxParents,
				     int node,
				     int n,
				     long N,
				     RecombinationParams *params)
{
  int    label;
  double scoreBefore;
  double scoreAfter;
  double gain;
  LabeledTreeNode *dummy;
  Value *leftValue0,*leftValue1;
  Value *rightValue0,*rightValue1
;
  double leftContribution;
  double rightContribution;

  // if this bit is almost fixed now, ignore it

  if ((x->value[0]<=FIXED_THRESHOLD)||(x->value[1]<=FIXED_THRESHOLD))
    {
      for (label=0; label<n; label++)
	x->dArrayTmp[label]=-1;
    };

  // we'll need this dummy node

  dummy = new LabeledTreeNode(LEAF);

  // depth is too big - we can't split

  if (x->depth>=maxParents)
    {
      for (label=0; label<n; label++)
	x->dArrayTmp[label]=-1;
      return 0;
    };

  // allocate memory for the frequencies of the splits

  leftValue0  = (Value*) calloc(n,sizeof(Value));
  leftValue1  = (Value*) calloc(n,sizeof(Value));
  rightValue0 = (Value*) calloc(n,sizeof(Value));
  rightValue1 = (Value*) calloc(n,sizeof(Value));

  t->computeSplitFrequencies(x,leftValue0,leftValue1,rightValue0,rightValue1);

  // compute basic contribution of this node (before split)

  scoreBefore = nodeContribution(x,N);

  // try all splits we can do (must not be our parent, and must not
  // create a cycle)

  for (label=0; label<n; label++)
    {
      if ((node!=label)&&
      	  (x->parentLabelCoincidenceVector[label]==0)&&
      	  ((G->connected(label,node))||
	   ((G->getNumIn(node)<maxParents)&&(G->canAddEdge(label,node)))))
	{
	  // compute frequencies if we did split
	  
	  dummy->value[0] = leftValue0[label];
	  dummy->value[1] = leftValue1[label];
	  leftContribution = nodeContribution(dummy,N);
	  
          ///printf("dummy: %f %f %f\n",dummy->value[0],dummy->value[1],leftContribution);

	  dummy->value[0] = rightValue0[label];
	  dummy->value[1] = rightValue1[label];
	  rightContribution = nodeContribution(dummy,N);
	  
          ///printf("dummy: %f %f %f\n",dummy->value[0],dummy->value[1],rightContribution);

	  // new contribution
	  
	  scoreAfter = leftContribution+rightContribution;
	  
	  // compute the gain
	  
	  if (params->priorNetwork)//!!!&&(!G->connected(label,node)))
	    if (priorConnected(label,node))
	      {
		gain = scoreAfter-scoreBefore;
		///	      printf("Using prior network edge %u - %u\n",label,node);
	      }
	    else
	      gain = scoreAfter-scoreBefore - log2(N)/2-params->logKappa;
	  else
	    if (metricToUse==BAYESIAN)
	      if (penalty)
		gain = scoreAfter-scoreBefore - log2(N)/2;
	      else
		gain = scoreAfter-scoreBefore;
	    else
	      {
		gain = scoreAfter-scoreBefore-log2(N)/2;
		//		printf("Node %u: Split on %u has a gain %f (%f,%f,%f)\n",node,label,gain,scoreBefore,scoreAfter,log2(N)/2); 
		//getchar(); 
	      };
	  
	  // update the contribution
	  
	  x->dArrayTmp[label]=gain;
	  
	  // update the values we computed (for future reference)
	  
#ifdef STORE_SPLIT_FREQUENCIES
	  x->rightValue0Array[label]=rightValue0[label];
	  x->rightValue1Array[label]=rightValue1[label];
	  x->leftValue0Array[label]=leftValue0[label];
	  x->leftValue1Array[label]=leftValue1[label];
#endif
	}
      else
	x->dArrayTmp[label]=-1;
    }
  
  // delete the dummy node
  
  delete dummy;

  // free the memory used by the frequencies of the splits

  free(leftValue0);
  free(leftValue1);
  free(rightValue0);
  free(rightValue1);
  
  // get back
  
  return 0;
};


// ====================================================================

int recomputeDecisionGraphMergeGains(FrequencyDecisionGraph *t,
				     AcyclicOrientedGraph *G,
				     MergeOperator **merge,
				     int *numMerges,
				     int node,
				     int n,
				     long N)
{
  double scoreBefore;
  double scoreAfter;
  double gain;
  LabeledTreeNode *dummy;
  NodeListItem *a, *b;
  int numLeaves;

  return 0; // this prohibits merge operation!!!!!!!! check and do the corrections!!!!!!

  // allocate a dummy node 

  dummy = new LabeledTreeNode(LEAF);

  // assign helper variables

  numLeaves = t->getNumLeaves();

  // allocate the merge operator array

  if (merge[node]!=NULL)
    free(merge[node]);
  
  merge[node]   = (MergeOperator*) calloc(numLeaves*(numLeaves+1)/2,sizeof(MergeOperator));
  *numMerges    = 0;

  // compute merge gains for all pair of nodes
 
  for (a=t->getLeaves(); a->next!=NULL; a=a->next)
    for (b=a->next; b!=NULL; b=b->next)
      {
	// compute the score before the merge

	scoreBefore = nodeContribution(a->x,N)+nodeContribution(b->x,N);

	// update the dummy node (its frequencies are equal the sum of the frequencies 
	//                        of the two merged nodes)

	dummy->value[0] = a->x->value[0]+b->x->value[0];
	dummy->value[1] = a->x->value[1]+b->x->value[1];

	// compute the score after the merge (without really having to do the merge)

	scoreAfter = nodeContribution(dummy,N);

	// compute the gain

	if (metricToUse==BAYESIAN)
	  if (penalty)
	    gain = scoreAfter-scoreBefore+log2(N)/2;
	  else
	    gain = scoreAfter-scoreBefore;
	else
	  gain = scoreAfter-scoreBefore+log2(N)/2;
	
	if (gain>0)
	  {
	    // compute the score after the merge
	    
	    merge[node][*numMerges].gain = gain;
	    merge[node][*numMerges].a    = a->x;
	    merge[node][*numMerges].b    = b->x;
	    
	    (*numMerges)++;
	  };
      };

  // delete the dummy node

  delete dummy;

  // get back

  return 0;
};

// ====================================================================

double logGamma(long x)
{
  if (x==1)
    return 0;
  else
    return getPrecomputedCummulativeLog(1,x-1);
};

// ====================================================================

double nodeContribution(LabeledTreeNode *x, long N)
{
  double p_x_0;
  double p_x_1;
  long   m_x_0;
  long   m_x_1;
  double p_p;
  long   m_p;
  double score;

  // compute the probability and frequency of the node x 

  p_x_0 = x->value[0];
  m_x_0 = (long) ((double)p_x_0*N+0.5);
  p_x_1 = x->value[1];
  m_x_1 = (long) ((double)p_x_1*N+0.5);

  ///printf("p0,p1 = %f,%f\n",p_x_0,p_x_1);
  ///printf("m0,m1 = %lu %lu\n",m_x_0,m_x_1);

  // is this node reachable? if not get the hell out of here

  if ((p_x_0==0)&&(p_x_1==0))
    return -10E10;

  if (metricToUse==BAYESIAN)
    {
      // alright, let's go for it, compute the probability and frequency of the parents
      
      p_p = p_x_0 + p_x_1;
      m_p = (long) ((double)p_p*N +0.5);
      
      // compute basic contribution of this node 
      
      score = -logGamma(2+m_p) + logGamma(1+m_x_0) + logGamma(1+m_x_1);
    }
  else
    {
      score = m_x_0*log2(p_x_0/(p_x_0+p_x_1))+m_x_1*log2(p_x_1/(p_x_0+p_x_1));
      
/*       printf("x1 = %f\n",log2(p_x_0/(p_x_0+p_x_1))); */
/*       printf("x2 = %f\n",log2(p_x_1/(p_x_0+p_x_1))); */
    }

/*   printf("Final score: %f\n",score); */

  // return the score
  
  return score;
};

// ====================================================================

int resetOperator(Operator *x)
{
  // reset the operator (empty and stuff)

  x->gain  = -1;
  x->type  = OPERATOR_NONE;
  x->where = -1;
  x->label = -1;
  x->node  = NULL;
  x->node2 = NULL;

  // get back

  return 0;
};

// ====================================================================

int updateBestOperator(Operator *best, Operator *x)
{
  // take the one with a better gain

  if (x->gain>best->gain)
    memcpy(best,x,sizeof(Operator));

  // get back
  
  return 0;
};

// ====================================================================

int updateBestNodeOperator(Operator *x, FrequencyDecisionGraph *t, MergeOperator *merge, int numMerges, int node, int n)
{
  int i,j;
  int numLeaves;
  NodeListItem *leafItem;

  //  printf("Updating best gain for node %u (%u + %u)\n",node,t->getNumLeaves(),numMerges);

  // assign some helper variables

  numLeaves = t->getNumLeaves();

  // reset the operator

  resetOperator(x);

  // try all splits, and update the best oeprator for this node

  for (i=0, leafItem=t->getLeaves(); (i<numLeaves); i++, leafItem=leafItem->next)
    for (j=0; j<n; j++)
      if (x->gain<leafItem->x->dArrayTmp[j])
	{
	  x->gain  = leafItem->x->dArrayTmp[j];
	  x->type  = OPERATOR_SPLIT_NODE;
	  x->where = node;
	  x->label = j;
	  x->node  = leafItem->x;
	  x->node2 = NULL;
	  x->t     = t;
	};

  // try all merges

  for (i=0; i<numMerges; i++)
    if (merge[i].gain>x->gain)
      {
	x->gain  = merge[i].gain;
	x->type  = OPERATOR_MERGE_NODE;
	x->where = node;
	x->label = 0;
	x->node  = merge[i].a;
	x->node2 = merge[i].b;
	x->t     = t;
      };

  // get back
  
  return 0;
};

// ====================================================================

int applyOperator(Operator *x, AcyclicOrientedGraph *G)
{
  // apply the operator

  switch (x->type) {
  case OPERATOR_SPLIT_NODE: 
    x->t->split(x->node,x->label);
    x->node->left->dArrayTmp = (Value*) calloc(G->size(),sizeof(Value));
    x->node->right->dArrayTmp = (Value*) calloc(G->size(),sizeof(Value));
#ifdef STORE_SPLIT_FREQUENCIES
    x->node->left->rightValue0Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->left->rightValue1Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->left->leftValue0Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->left->leftValue1Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->right->rightValue0Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->right->rightValue1Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->right->leftValue0Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->right->leftValue1Array = (Value*) calloc(G->size(),sizeof(Value));
    x->node->left->value[0]=x->node->leftValue0Array[x->label];
    x->node->left->value[1]=x->node->leftValue1Array[x->label];
    x->node->right->value[0]=x->node->rightValue0Array[x->label];
    x->node->right->value[1]=x->node->rightValue1Array[x->label];
    x->node->freeTemporaryArrays();
#else
  x->t->computeFrequencies();
#endif

    if (x->node->depth+1>maxDepth)
      maxDepth=x->node->depth+1;
    G->addEdge(x->label,x->where);
    numAppliedSplits++;
#ifdef PRINT_SPLITS
    printf("Split %u on %u applied, with a gain %1.2f\n",x->label,x->where,x->gain);
    //    getchar();
#endif
    break;

  case OPERATOR_MERGE_NODE:
    x->t->merge(x->node,x->node2);
    numAppliedMerges++;
    break;

  default:
    printf("ERROR: Operator not defined (applyOperator)! Exiting!\n");
  };

  // get back

  return 0;
};

// ====================================================================

int updateGainsAfterOperator(Operator *x, AcyclicOrientedGraph *G, int maxParents, long N, RecombinationParams *params)
{
   switch (x->type) {
   case OPERATOR_SPLIT_NODE: 
     recomputeDecisionGraphSplitGains(x->t,x->node->left,G,maxParents,x->where,G->size(),N,params);
     recomputeDecisionGraphSplitGains(x->t,x->node->right,G,maxParents,x->where,G->size(),N,params);
     break;

   case OPERATOR_MERGE_NODE:
     recomputeDecisionGraphSplitGains(x->t,x->node,G,maxParents,x->where,G->size(),N,params);
     break;


  default:
    printf("ERROR: Operator not defined (updateGainsAfterOperator)! Exiting!\n");
  };

  // get back

  return 0;
};

// ====================================================================

int operatorApplicable(Operator *x, AcyclicOrientedGraph *G)
{
  // we can split only when we can add a parent or it's already there

  if (x->type==OPERATOR_SPLIT_NODE)
    return ((G->connected(x->label,x->where))||(G->canAddEdge(x->label,x->where)));
  else
    if (x->type==OPERATOR_MERGE_NODE)
      return 1;

  // get back

  return 0;
};

// ====================================================================

int deleteOperator(FrequencyDecisionGraph **T, Operator *x)
{
  // well, if it's the split, set the corresponding gain to -1

  if (x->type==OPERATOR_SPLIT_NODE)
    x->node->dArrayTmp[x->label]=-1;

  // get back
  
  return 0;
};

char *dBOAGetMetricDescription(int n)
{
  return dBOAMetricDesc[n];
};
