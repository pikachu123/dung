#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "memalloc.h"
#include "population.h"
#include "random.h"
#include "mymath.h"
#include "userBOA.h"

#ifndef M_PI
	#define M_PI (3.14159265358979323846)
#endif

//#define RANDOM_CLAUSE_INJECTION
//#define ORDER_1_INJECTION
#define HILLCLIMBING_INJECTION

//#define DEBUG

#define numFitness 11

#define RATIO     0.2
#define HIGH_PEAK 1
#define LOW_PEAK  0.9

#define BISECTION_PENALTY 0.1

#define newedge(x,y) { node1[numEdges]=(x); node2[numEdges]=(y); numEdges++; coincidence[x][y]=coincidence[y][x]=1; }

static double onemax(char *x, int numDiscrete, double *c, int numContinuous);
static double trap4(char *x, int numDiscrete, double *c, int numContinuous);
static double trap5(char *x, int numDiscrete, double *c, int numContinuous);
static double realvalued(char *x, int numDiscrete, double *c, int numContinuous);
static double F1(char *x, int numDiscrete, double *c, int numContinuous);
static double HIFF(char *x, int numDiscrete, double *c, int numContinuous);
static double hierarchicalTrap3(char *x, int numDiscrete, double *c, int numContinuous);
static double hierarchicalDeceptiveTrap3(char *x, int numDiscrete, double *c, int numContinuous);
static double continuousDeceptiveTwoPeaks(char *x, int numDiscrete, double *continuous, int numContinuous);

static double httpEvaluator(char *x, int numDiscrete, double *continuous, int numContinous);

static char areAllGenesOne(char *x, int numDiscrete, double *c, int numContinuous, char type, double f);
static char areAllGenesZeroOrOne(char *x, int numDiscrete, double *c, int numContinuous, char type, double f);
static char realValuedIsBest(char *x, int numDiscrete, double *c, int numContinuous, char type, double f);
static char continuousDeceptiveTwoPeaksIsBest(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double f);

static double onemaxGoodBBs(char *x, int numDiscrete, double *c, int numContinuous);
static double trap4GoodBBs(char *x, int numDiscrete, double *c, int numContinuous);
static double trap5GoodBBs(char *x, int numDiscrete, double *c, int numContinuous);
static double proportionOfOnes(char *x, int numDiscrete, double *c, int numContinuous);

static long getNumOptimaHIFF(int numDiscrete, int numContinuous);
static long whichOptimumHIFF(char *x, int numDiscrete, double *c, int numContinuous);

static int trap(char *x, int n, int order);

static double peak(double x, double center, double width, double height);
static double generalTrap(char *x, int k, double leftPeak, double rightPeak);

double  additionalFitnessNoiseVariance;
double  additionalFitnessNoiseDeviation;
char realValued;
int numOptima;

int numNodes;
int numEdges;
int *node1;
int *node2;
char **coincidence=NULL;

double twomaxFactor;

// --------------------------------
// params desc:
//
// 0...length of chromosome
// 1...left bound
// 2...right bound
// 3...two-max bias to the right
// 4...number of optima for real valued (with variable number of optima)
// --------------------------------

double userEvaluator(char *x, int numDiscrete, double *continuous, int numContinous);
char userIsBest(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double f);
double userGoodBBs(char *x, int numDiscrete, double *continuous, int numContinuous);

static Fitness fitness[numFitness] = {
  {"ONEMAX",&onemax,&areAllGenesOne,onemaxGoodBBs,NULL,NULL,NULL,NULL,NULL,NULL},
  {"4-ORDER TRAP",&trap4,&areAllGenesOne,trap4GoodBBs,NULL,NULL,NULL,NULL,NULL,NULL},
  {"5-ORDER TRAP",&trap5,&areAllGenesOne,trap5GoodBBs,NULL,NULL,NULL,NULL,NULL,NULL},
  {"REAL-VALUED (with variable number of optima)",&realvalued,&realValuedIsBest,NULL,&initRealValuedFitness,NULL,NULL,NULL,NULL,NULL},
  {"F1 (Goldberg&Richardson)",&F1,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL},
  {"HIFF (Watson, Pollack)",&HIFF,&areAllGenesZeroOrOne,NULL,NULL,NULL,NULL,NULL,&getNumOptimaHIFF,&whichOptimumHIFF},
  {"Hierarchical trap 3 (d=0.8,0.9)",&hierarchicalTrap3,&areAllGenesOne,NULL,NULL,NULL,NULL,NULL,NULL,NULL},
  {"Hierarchical deceptive trap 3 (d=0.8,0.9)",&hierarchicalDeceptiveTrap3,&areAllGenesOne,NULL,NULL,NULL,NULL,NULL,NULL,NULL},
  {"Continuous deceptive made from twoPeaks",&continuousDeceptiveTwoPeaks,&continuousDeceptiveTwoPeaksIsBest,NULL,NULL,NULL,NULL,NULL,NULL,NULL},
  {"Remote Evaluator", &httpEvaluator, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
  {"User Evaluator", &userEvaluator, &userIsBest, &userGoodBBs, NULL,NULL,NULL,NULL,NULL,NULL}
};

double a,b;
double coef;
long   fitnessCalls_;

double *NK_coeff;
int **NK_neighbors;
int NK_k;
int NK_k1;
int NK_n=-1;

char **priorCoincidence;

// ============================================
// User Evaluator
double userEvaluator(char *x, int numDiscrete, double *continuous, int numContinous) 
{
	if (fitnessDefinition.fitness)
		return fitnessDefinition.fitness(numDiscrete, x);
	else
		return 0.0;
}

// ==============================================
// User IsBest

char userIsBest(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double f)
{
	if (fitnessDefinition.isBest)
		return fitnessDefinition.isBest(numDiscrete, x);
	else
		return 0;
}

// ==============================================

double userGoodBBs(char *x, int numDiscrete, double *continuous, int numContinuous)
{
	if (fitnessDefinition.goodBBs)
		return fitnessDefinition.goodBBs(numDiscrete, x);
	else
		return 0.0;
}

// ============================================
// Remote Evaluator
double httpEvaluator(char *x, int numDiscrete, double *continuous, int numContinous) 
{
  return 0.0;
}

// ============================================
// ONEMAX fitness function

double onemax(char *x, int numDiscrete, double *continuous, int numContinuous)
{
   int s;

   s=0;

   for (register int i=0; i<numDiscrete; i++)
       s += x[i];
   return (double) s;
}

// ============================================
// order-4 trap function

double trap4(char *x, int numDiscrete, double *continuous, int numContinuous)
{
   return (double) trap(x,numDiscrete,4)/*+n/5*/;
}

double trap5(char *x, int numDiscrete, double *continuous, int numContinuous)
{
   return (double) trap(x,numDiscrete,5)/*+n/5*/;
}

// ============================================
// real valued with variable number of optima

double realvalued(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  double f,v;

  v = realValue(x,numDiscrete);

  f = cos(M_PI*v*numOptima)*(1-v*v*0.05);

  return f;
}

// ==========================================================
// F1 fitness function (Mahfoud) 5 equal peaks in (0,1)
// maxima = 0.1,0.3,0.5,0.7,0.9

double F1(char *x, int numDiscrete, double *continuous, int numContinuous)
{
   double d;
   int k;

   d=0;
   for (k=0; k<numDiscrete; k++)
     {
       d *= 2;
       d += x[k];
     };
   
   d /= (double) ((1<<(numDiscrete))-1);

   return (double) pow(sin(5*M_PI*d),6);
}

// ==============================================

double HIFF(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int    i;
  char   *transform;
  char   nullValue = 3;
  int    level;
  int    numLevels;
  int    last;
  double result;
  int    bonus;

  // set the variables

  numLevels = (int) ((double)log2(numDiscrete)+0.5); 
  last      = numDiscrete;
  bonus     = 2;
  result    = numDiscrete;   // to make it consistent with the definition by Watson...

  // allocate some memory

  transform = (char*) Malloc(numDiscrete);

  // copy the string to the transform array

  for (i=0; i<numDiscrete; i++)
    transform[i]=x[i];
  
  // process all levels (add reinforcements and interpret the string to a higher level)

  for (level=0; level<numLevels; level++)
    {
      for (i=0; i<last; i+=2)
	{
	  if ((transform[i]==1)&&(transform[i+1]==1))
	    {
	      result += bonus;
	      transform[i>>1]=1;
	    }
	  else
	    if ((transform[i]==0)&&(transform[i+1]==0))
	      {
		result += bonus;
		transform[i>>1]=0;
	      }
	    else
	      transform[i>>1]=nullValue;
	};
      
      // increase the bonus on a higher level

      bonus<<=1;

      // we've got only half interpretations at this point

      last>>=1;
    };

  // free the memory

  Free(transform);

  // get back

  return result;
};

long getNumOptimaHIFF(int numDiscrete, int numContinuous)
{
  return 2;
}

long whichOptimumHIFF(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int s;
  
  s=0;
  
  for (int i=0; i<numDiscrete; i++)
    s += x[i];
  
  if (s==0)
    return 0;
  else
    if (s==numDiscrete)
      return 1;
    else
      return -1;
}

// ==============================================
// are all genes one?

char areAllGenesOne(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double f)
{
  register int i;

  for (i=0; (i<numDiscrete)&&(x[i]==1); i++);
    
  return (i==numDiscrete);
}

// ==============================================

double hierarchicalTrap3(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int    i;
  char   *transform;
  char   nullValue = 3;
  int    level;
  int    numLevels;
  int    last;
  double result;
  int    bonus;
  int    tmp;

  // set the variables

  numLevels = (int) ((double)log3(numDiscrete)+0.5); 
  last      = numDiscrete;
  bonus     = 1;
  result    = 0;

  // allocate some memory

  transform = (char*) Malloc(numDiscrete);

  // copy the string to the transform array

  for (i=0; i<numDiscrete; i++)
    transform[i]=x[i];
  
  // process all levels (add reinforcements and interpret the string to a higher level)

  for (level=0; level<numLevels; level++)
    {
//       printf("Level %u\n",level);
//       for (i=0; i<last; i++)
// 	{
// 	  if ((i!=0)&&(i%3==0))
// 	    printf(" ");
// 	  printf("%u",transform[i]);
// 	};
//       printf("\n\n");

      for (i=0; i<last; i+=3)
	{
	  if (level<numLevels-1)
	    result += bonus*generalTrap(transform+i,3,1,1);
	  else
	    result += bonus*generalTrap(transform+i,3,0.9,1);

	  tmp = transform[i]+transform[i+1]+transform[i+2];

	  if (tmp==0)
	    transform[i/3]=0;
	  else
	    if (tmp==3)
	      transform[i/3]=1;
	    else
	      transform[i/3]=nullValue;
	};
       
      // increase the bonus on a higher level

      bonus*=3;

      // we've got only half interpretations at this point

      last/=3;
    };

//   getchar();

  // free the memory

  Free(transform);

  // get back

  return result;
};

// ===================================================================

double hierarchicalDeceptiveTrap3(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int    i;
  char   *transform;
  char   nullValue = 3;
  int    level;
  int    numLevels;
  int    last;
  double result;
  int    bonus;
  int    tmp;

  // set the variables

  numLevels = (int) ((double)log3(numDiscrete)+0.5); 
  last      = numDiscrete;
  bonus     = 1;
  result    = 0;

  // allocate some memory

  transform = (char*) Malloc(numDiscrete);

  // copy the string to the transform array

  for (i=0; i<numDiscrete; i++)
    transform[i]=x[i];
  
  // process all levels (add reinforcements and interpret the string to a higher level)

  for (level=0; level<numLevels; level++)
    {
//       printf("Level %u\n",level);
//       for (i=0; i<last; i++)
// 	{
// 	  if ((i!=0)&&(i%3==0))
// 	    printf(" ");
// 	  printf("%u",transform[i]);
// 	};
//       printf("\n\n");

      for (i=0; i<last; i+=3)
	{
	  if (level<numLevels-1)
	    {
// 	      printf("Giving bonus of %f to the string %i%i%i\n",generalTrap(transform+i,3,1,1-0.2/(numLevels)),transform[i],transform[i+1],transform[i+2]);
// 	      getchar();
	      result += bonus*generalTrap(transform+i,3,1+0.1/numLevels,1);
	    }
	  else
	    result += bonus*generalTrap(transform+i,3,0.9,1);

	  tmp = transform[i]+transform[i+1]+transform[i+2];

	  if (tmp==0)
	    transform[i/3]=0;
	  else
	    if (tmp==3)
	      transform[i/3]=1;
	    else
	      transform[i/3]=nullValue;
	};
       
      // increase the bonus on a higher level

      bonus*=3;

      // we've got only one third of the interpretations at this point

      last/=3;
    };

//   getchar();

  // free the memory

  Free(transform);

  // get back

  return result;
};

double continuousDeceptiveTwoPeaks(char *x, int numDiscrete, double *continuous, int n)
{
  double f=0;
  for (int i=0; i<n; i+=2)
    {
      double a= continuous[i];
      double b= continuous[i+1];
      double v=sqrt((dsqr(a)+dsqr(b))/2.0);
      
      if (v<0)
	f+=0;
      else
	if (v<RATIO)
	  f+=peak(v,RATIO/2,RATIO,HIGH_PEAK);
	else
	  if (v<1)
	    f+=peak(v,RATIO+(1-RATIO)/2,1-RATIO,LOW_PEAK);
	  else
	    f+=0;
    };

  return f;
};

// ==============================================
// are all genes one?

char areAllGenesZeroOrOne(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double f)
{
  register int i;

  for (i=0; (i<numDiscrete)&&(x[i]==1); i++);

  if (i<numDiscrete)
    for (i=0; (i<numDiscrete)&&(x[i]==0); i++);    
    
  return (i==numDiscrete);
}

// ==============================================
// returns the number of ones

double proportionOfOnes(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  register int i;
  int s;

  s=x[0];

  for (i=1; i<numDiscrete; i++)
    s+=x[i];
    
  return (double) s/numDiscrete;
}

// =============================================

char realValuedIsBest(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double fit)
{
  double f;

  f = realValue(x,numDiscrete);

  return (fabs(f)<=1E-4);
};

// ==============================================

double onemaxGoodBBs(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int i,k;

  k=0;

  for (i=0; i<numDiscrete; i++)
    {
      k += x[i];
    }

  return double(k)/numDiscrete;
};

// ==============================================

double trap4GoodBBs(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int i,k,s;

  k=0;

  for (i=0; i<numDiscrete;)
    {
      s  = x[i++];
      s += x[i++];
      s += x[i++];
      s += x[i++];

      if (s==4)
	k++;
    }

  return double(k)*4/numDiscrete;
};


double trap5GoodBBs(char *x, int numDiscrete, double *continuous, int numContinuous)
{
  int i,k,s;

  k=0;

  for (i=0; i<numDiscrete;)
    {
      s  = x[i++];
      s += x[i++];
      s += x[i++];
      s += x[i++];
      s += x[i++];
      if (s==5)
	k++;
    }

  return double(k)*5/numDiscrete;
};

char continuousDeceptiveTwoPeaksIsBest(char *x, int numDiscrete, double *continuous, int numContinuous, char type, double fit)
{
  double f;
  
  f = 0;
  for (int i=0; i<numContinuous; i+=2)
    {
      double a= continuous[i];
      double b= continuous[i+1];
      double v=sqrt((dsqr(a)+dsqr(b))/2.0);
      f+=dsqr(v-RATIO/2);
    };
  
  return (f<1E-2);
};

// ==============================================

void *getFitness(int n)
{
  realValued = 0;

  if ((n>=0) && (n<numFitness))
    return (void*) &(fitness[n]);
  else
    {
      fprintf(stderr,"ERROR: Specified fitness function doesn't exist (%u)!",n);
      exit(-1);
    }

  return NULL;
}

char *getFitnessDesc(int n)
{
  return fitness[n].description;
};

char isRealValued()
{
  return realValued;
};

int initRealValuedFitness(double *params)
{
   int n = (int) params[0];
   a = params[1];
   b = params[2];
   coef = (double) (b-a)/(1L<<n);
   
   realValued = 1;
   numOptima  = (int) params[4];

   return 0;
}


double realValue(char *x, int n)

{
   register int i;
   double k,p;

   k=0;
   p=1;

   for (i=n-1; i>=0; i--)
   {
      if (x[i])
	 k += p;
      p*=2;
   }


   return (double) a+coef*k;
}

int resetFitnessCalls(void)
{
  return fitnessCalls_=0;
}

long fitnessCalled(void)
{
  return fitnessCalls_++;
}

long fitnessCalled(long numCalls)
{
  return fitnessCalls_+=numCalls;
}

long getFitnessCalls(void)
{
  return fitnessCalls_;
}

int trap(char *x, int n, int order)
{
   int order1,s;
   register int i,j,k;

   order1 = order-1;

   s = 0;
   i = 0;
   
   do
   {
     k=0;
     for (j=0; j<order; j++)
         k+=x[i++];

     if (k==order)
        s+=order;
     else
        s+=order1-k;
   } while (i<n);

   return s;
}

// returns true if the two nodes are connected (for prior network stuff)

int priorConnected(int i, int j)
{
  if (coincidence==NULL)
    {
      fprintf(stderr,"ERROR: Prior network is not defined and can't be used!\n");
      exit(-1);
    };

  return (coincidence[i][j]);
};

// injects good solutions to the initial population (for incorporating prior info into the initial population)
//
// to be used after random generation - because it actually generates the guys at random (at least here)
// also must have the fitness evaluated!

int bisectionInjectClusters(void *P)
{
  long i;
  int ii,j,l,m;
  int n;
  long N;
  int *selected;
  char *best;
  int numBest;
  int *edgesToSelected;
  int *otherEdges;
  int picked;
  int clusterSize;

  Population *p;

  // this is the population (we had to use a void pointer)

  p = (Population*) P;

  // assign the helper variables

  N = p->N;
  n = p->n;

  // allocate some memory

  selected = (int*) Calloc(n,sizeof(int));
  best     = (char*) Malloc(n);
  edgesToSelected = (int*) Calloc(n,sizeof(int));
  otherEdges = (int*) Calloc(n,sizeof(int));

  // how big clusters to insert

  clusterSize=(int)((double)0.05*n);

  // go for it

  for (i=0; i<N; i++)
    {
      // reset the individual

      for (l=0; l<n; l++)
	p->individual[i].chromosome[l]=0;

	  // select the first random position in the string at random

      selected[0]=intRand(n);

      // select other positions according with most edges into the selected guys, when a tie encountered
      // we select the one with the least other edges. if that's a tie too, we pick a random one from the
      // best candidates.

      for (l=1; l<clusterSize; l++)
	{
	  int maxEdgesToSelected;
	  int minOtherEdges;

	  for (m=0; m<n; m++)
	    {
	      int alreadySelected=0;

	      for (j=0; j<l; j++)
		if (selected[j]==m)
		  alreadySelected=1;

	      if (!alreadySelected)
		{
		  edgesToSelected[m] = 0;
		  otherEdges[m]=0;

		  for (j=0; j<numEdges; j++)
		    if ((node1[j]==m)||(node2[j]==m))
		      {
			int toSelected=0;

			for (ii=0; ii<l; ii++)
			  if ((node1[j]==selected[ii])||
			      (node2[j]==selected[ii]))
			    toSelected=1;

			if (toSelected)
			  edgesToSelected[m]++;
			else
			  otherEdges[m]++;	
		  
			// find the best candidate 
		      };
		}
	      else
		edgesToSelected[m]=-1;
	    };

	  maxEdgesToSelected=-1;
	  minOtherEdges=-1;

	  for (m=0; m<n; m++)
	    if ((edgesToSelected[m]>maxEdgesToSelected)||
		((edgesToSelected[m]==maxEdgesToSelected)&&
		 (otherEdges[m]<minOtherEdges)))
	      {
		maxEdgesToSelected = edgesToSelected[m];
		minOtherEdges      = otherEdges[m];
	      };
	      
	  // count all guys that are just as good
	      
	  numBest=0;
	  for (m=0; m<n; m++)
	    if ((edgesToSelected[m]==maxEdgesToSelected)&&
		(otherEdges[m]==minOtherEdges))
	      {
		best[m]=1;
		numBest++;
	      }
	    else
	      best[m]=0;
	      
	  // pick one of the best
	      
	  picked = intRand(numBest);
	      
	  // get him
	      
	  for (m=0; m<n; m++)
	    if (best[m]==1)
	      {
		if (picked==0)
		  {
		    selected[l]=m;
		    picked=-1;
		  }
		else
		  picked--;
	      };
	};
	  

      for (m=0; m<n/2; m++)
	p->individual[i].chromosome[selected[m]]=1;
    };

  // free the memory

  Free(selected);
  Free(best);
  Free(edgesToSelected);
  Free(otherEdges);

  // get back

  return 0;
};

// ---------------------------------------------------------
// set the additional noise to the fitness

int setFitnessNoiseVariance(double variance)
{
  additionalFitnessNoiseVariance  = variance;
  additionalFitnessNoiseDeviation = sqrt(variance);

  return 0;
};

// ---------------------------------------------------------
// get the additional noise to the fitness

double getFitnessNoiseVariance()
{
  return additionalFitnessNoiseVariance;
};

// ---------------------------------------------------------
// return the noisy variable (if none, return 0)

double additionalFitnessNoise()
{
  if (additionalFitnessNoiseDeviation>0)
    return gaussianRandom(0,additionalFitnessNoiseDeviation);
  else
    return 0;
};

// makes a trap with specified parameters
// here d is the ratio of deceptive and isolated attractors,
// k is the length of the block
// and x is input
// the bigger peak is always 1
// -> allows nullValue (for hierarchical stuff)

double generalTrap(char *x, int k, double leftPeak, double rightPeak)
{
  int s;
  
  s=0;
  for (int i=0; i<k; i++)
    if (x[i]==1)
      s++;
    else
      if(x[i]!=0)
	return 0;
  
  if (s==k)
    return rightPeak;
  else
    return leftPeak*(k-1-s)/(k-1);
};

// returns a peak function with a center and width (based on cos)

double peak(double x, double center, double width, double height)
{
  return height*(cos(2/width*M_PI*(x-center))+1)/2;
};
