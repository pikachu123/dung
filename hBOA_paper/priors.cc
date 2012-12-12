#include <math.h>
#include "priors.h"
#include "memalloc.h"

#define NOPRIORNETWORK    -1
#define EMPTYPRIORNETWORK  0
#define LINCPRIORNETWORK   1

#define EPSILON 1E-03

#define numPriorNetworkSourceDescriptions 3

//#define DEBUG

char                    *priorNetworkDescriptor=NULL;
PriorNetworkConstructor *priorNetworkConstructor=NULL;
char                     priorNetworkDefined=0;

PriorNetworkSourceDescription priorNetworkSourceDescription[numPriorNetworkSourceDescriptions] = {
  {"None",NULL},
  {"Empty",&emptyPriorNetwork},
  {"LINC",&LINCPriorNetwork}
};

//==========================================================

int setPriorNetworkSource(int source)
{
  // set the source

  priorNetworkConstructor = priorNetworkSourceDescription[source].constructor;
  priorNetworkDescriptor  = priorNetworkSourceDescription[source].description;
  priorNetworkDefined     = (priorNetworkConstructor!=NULL);

  // get back

  return 0;
};

//==========================================================

char *getPriorNetworkDescription(int source)
{
  //  return priorNetworkDescriptor;
  return priorNetworkSourceDescription[source].description;
};

//==========================================================

char isPriorNetworkDefined()
{
  return priorNetworkDefined;
};

//==========================================================

int allocatePriorNetwork(OrientedGraph ***priorNetwork, int n, int d)
{
  int i;

  // are we using a prior network?
 
  if (priorNetworkDefined)
    {
      // allocate the memory for the array of prior networks

      (*priorNetwork) = (OrientedGraph**) Calloc(d,sizeof(OrientedGraph**));
      
      // initialize each of them
      
      for (i=0; i<d; i++)
	(*priorNetwork)[i] = new OrientedGraph(n);
    };

  // get back

  return 0;
};

//==========================================================

int freePriorNetwork(OrientedGraph ***priorNetwork, int n, int d)
{
  int i;

  // are we using a prior network?

  if (priorNetworkConstructor)
  {

    // call destructor of each prior network

    for (i=0; i<d; i++)
      delete ((*priorNetwork)[i]);

    // free the array of networks

    Free(*priorNetwork);
  }

  // get back

  return 0;
};

//==========================================================

int constructPriorNetwork(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N)
{
  // if there is assigned a non-null constructor, use it and construct the prior network

  if (priorNetworkConstructor)
    (*priorNetworkConstructor)(P,fitness,priorNetwork,N);
  
  // get back

  return 0;
};

//==========================================================

int emptyPriorNetwork(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N)
{
  // empty the prior network

  priorNetwork->removeAllEdges();
  
  // get back
  
  return 0;
};

//==========================================================

int LINCPriorNetwork(Population *P, Fitness *fitness, OrientedGraph *priorNetwork, long N)
{
  long i;
  int k,l,numDiscrete,numContinuous;
  Individual *individual;
  char *x;
  double f;
  double f1,f2,f12;
  double nonlinearity;

  // initialize some variables
  
  numDiscrete  = P->numDiscrete;
  numContinuous  = P->numContinuous;
 
  // create the network

  priorNetwork->removeAllEdges();

  for (i=0; i<N; i++)
    {
      individual = &(P->individual[i]);
      x = individual->chromosome;
      f = individual->f;

      for (k=0; k<numDiscrete; k++)
	for (l=0; l<k; l++)
	  {
	    x[k] = 1-x[k];
	    individual->fCalculated = 0;
	    evaluateIndividual(individual,fitness,numDiscrete,numContinuous);
	    f1 = individual->f-f;

	    x[l] = 1-x[l];
	    individual->fCalculated = 0;
	    evaluateIndividual(individual,fitness,numDiscrete,numContinuous);
	    f12 = individual->f-f;
		
	    x[k] = 1-x[k];
	    individual->fCalculated = 0;
	    evaluateIndividual(individual,fitness,numDiscrete,numContinuous);
	    f2 = individual->f-f;

	    x[l] = 1-x[l];

	    nonlinearity = f12-f1-f2;
#ifdef DEBUG
	    printf("Nonlinearity = %f\n",nonlinearity);
#endif

	    if ((fabs(nonlinearity)>=EPSILON) && (!priorNetwork->connected(k,l)))
	      {
#ifdef DEBUG
		printf("Adding edge (%u,%u) to prior network (%f + %f != %f)\n",k,l,f1,f2,f12);
		getchar();
#endif
		priorNetwork->addEdge(k,l);
		priorNetwork->addEdge(l,k);
	      }
	    else
	      {
#ifdef DEBUG
		printf("Not adding (%u %u)\n",k,l);
		getchar();
#endif
	      };
	  }

      individual->f = f;
      individual->fCalculated = 1;
    }

  // get back

  return 0;
};
