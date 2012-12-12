#include <stdio.h>
#include <malloc.h>
#include <math.h>

#include "population.h"
#include "select.h"
#include "memalloc.h"
#include "random.h"
#include "utils.h"

static char *selectionDesc[8]={
			    "Proportional",
			    "Tournament with replacement",
                            "Truncation",
			    "Tournament without replacement",
			    "No selection",
			    "Tournament without replacement and continuous sharing",
			    "Normalized Boltzmann selection",
                            "Tournament (w/o repl.) for exp. scaled dec3",
			  };

int tournSize=2;
int treshold=2;
int nn;
int minHamm;
double beta;

inline int moreThan(Individual a, Individual b, int numDiscrete, int numContinuous)
{
  return (a.f>b.f);
};

inline int lessThan(Individual a, Individual b, int numDiscrete, int numContinuous)
{
  return (a.f<b.f);
};

/*

inline int moreThanExpDec3(Individual a, Individual b, int numDiscrete, int numContinuous)
{
  char *x = a.chromosome;
  char *y = b.chromosome;
 
  for (int i=0; i<numDiscrete; i+=3)
    {
      if ((x[i]!=y[i])||
	  (x[i+1]!=y[i+1])||
	  (x[i+2]!=y[i+2]))
	{
	  double fx = f3deceptive(x+i,3,NULL,0);
	  double fy = f3deceptive(y+i,3,NULL,0);
	  if (fx>fy)
	    return 1;
	  else
	    return 0;
	};
    };
  
  return 0;
};

inline int lessThanExpDec3(Individual a, Individual b, int numDiscrete, int numContinuous)
{
  char *x = a.chromosome;
  char *y = b.chromosome;
 
  for (int i=0; i<numDiscrete; i+=3)
    {
      if ((x[i]!=y[i])||
	  (x[i+1]!=y[i+1])||
	  (x[i+2]!=y[i+2]))
	{
	  double fx = f3deceptive(x+i,3,NULL,0);
	  double fy = f3deceptive(y+i,3,NULL,0);
	  if (fy>fx)
	    return 1;
	  else
	    return 0;
	};
    };
  
  return 0;
};

*/

//==============================================================

int selectionProportional(Population *population, Population *parents, long M)
{
  register long i,j;
  long N1,N;
  int numDiscrete,numContinuous;
  double *fS;
  double r;

  N  = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;
  N1 = N-1;

  allocatePopulation(parents,M,numDiscrete,numContinuous);

  fS = (double*) Calloc(N,sizeof(*fS));

  fS[0]=population->individual[0].f;
  for (i=1; i<N; i++)
    fS[i] = fS[i-1]+population->individual[i].f;
  for (i=0; i<N; i++)
    fS[i] /= fS[N1];

  for (j=0; j<M; j++)
  {
    r = drand();

    for (i=0; fS[i]<r; i++);

    copyIndividual(&(parents->individual[j]),&(population->individual[i]),numDiscrete,numContinuous);
  }

  Free(fS);

  return 0;
}

//==============================================================
int selectionTournamentWithReplacement(Population *population, Population *parents, long M)
{
  register long i,j;
  long N,max;
  int numDiscrete,numContinuous;
  long  *chosen;
  char generated;

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  shufflePopulation(population);

  allocatePopulation(parents,M,numDiscrete,numContinuous);

  chosen = (long*) Calloc(tournSize,sizeof(long));

  for (j=0; j<M; j++)
  {
    chosen[0]=longRand(N);
    max = chosen[0];

    for (i=1; i<tournSize; i++)
    {
       do {
	 generated=1;
	 chosen[i]=longRand(N);
	 for (int k=0; k<i; k++)
	     if (chosen[i]==chosen[k])
		generated=0;
       } while (!generated);

       if (moreThan(population->individual[chosen[i]],population->individual[max],numDiscrete,numContinuous))
	 max  = chosen[i];
    }

    copyIndividual(&(parents->individual[j]),&(population->individual[max]),numDiscrete,numContinuous);
  }
  
  Free(chosen);

  return 0;
}

//==============================================================

int selectionTournamentWithoutReplacement(Population *population, Population *parents, long M)
{
  register long i,j;
  long N,max;
  int numDiscrete,numContinuous;
  int k;
  Individual x;

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  allocatePopulation(parents,M,numDiscrete,numContinuous);
  allocateIndividual(&x,numDiscrete,numContinuous);

  i = 0;

  for (j=0; j<M; j++)
  {
    if (i+tournSize>N)
      {
	shufflePopulation(population);
	i=0;
      };
    
    max=i;
    copyIndividual(&x,&(population->individual[max]),numDiscrete,numContinuous);
    
    for (k=1; k<tournSize; k++)
      {
	i++;
	
	if (moreThan(population->individual[i],population->individual[max],numDiscrete,numContinuous))
	  {
	    max  = i;
	    copyIndividual(&x,&(population->individual[max]),numDiscrete,numContinuous);
	  };
      };

    copyIndividual(&(parents->individual[j]),&x,numDiscrete,numContinuous);
    i++;
  }

  freeIndividual(&x);

  return 0;
}

//==============================================================
/*
int selectionTournamentWithoutReplacementExpDec3(Population *population, Population *parents, long M)
{
  register long i,j;
  long N,max;
  int numDiscrete,numContinuous;
  int k;
  Individual x;

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  allocatePopulation(parents,M,numDiscrete,numContinuous);
  allocateIndividual(&x,numDiscrete,numContinuous);

  i = 0;

  for (j=0; j<M; j++)
  {
    if (i+tournSize>N)
      {
	shufflePopulation(population);
	i=0;
      };
    
    max=i;
    copyIndividual(&x,&(population->individual[max]),numDiscrete,numContinuous);
    
    for (k=1; k<tournSize; k++)
      {
	i++;
	
	if (moreThanExpDec3(population->individual[i],population->individual[max],numDiscrete,numContinuous))
	  {
	    max  = i;
	    copyIndividual(&x,&(population->individual[max]),numDiscrete,numContinuous);
	  };
      };

    copyIndividual(&(parents->individual[j]),&x,numDiscrete,numContinuous);
    i++;
  }

  freeIndividual(&x);

  return 0;
}
*/

int initSelectionTournament(int tournamentSize_)
{
  tournSize = tournamentSize_;
  return 0;
}

//==============================================================

int selectionTruncation2(Population *population, Population *parents, long M)
{
  register long i,j;
  long N;
  int numDiscrete,numContinuous;

  long max;

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  allocatePopulation(parents,M,numDiscrete,numContinuous);

  for (j=0; j<M; j++)
  {
    max = j;

    for (i=j+1; i<N; i++)
    {
      if (moreThan(population->individual[i],population->individual[max],numDiscrete,numContinuous))
	 max  = i;
    }

    if (max!=j)
      swapIndividuals(population,j,max);

    copyIndividual(&(parents->individual[j]),&(population->individual[j]),numDiscrete,numContinuous);
  }

  ///  printf("Selection performed..\n");

  return 0;
}

//==============================================================

int selectionTruncation(Population *population, Population *parents, long M)
{
  register long i;
  long N;
  int numDiscrete,numContinuous;

  // initialize some variables

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  // allocate memory for selected individuals

  allocatePopulation(parents,M,numDiscrete,numContinuous);

  // shuffle the individuals in the population

  shufflePopulation(population);

  // perform quick-sort but only so that first M are better than the rest - no need for more

  divideBest(population,0,N-1,numDiscrete, numContinuous,M);

  /*
  for (i=0; i<N; i++)
    {
      if (i==M)
	printf("------------------\n");

      printf("fitness[%3lu]= %1.2f\n",i,population->individual[i].f);
    };
  getchar();
  */

  // copy M best to the population of selected strings

  for (i=0; i<M; i++)
    copyIndividual(&(parents->individual[i]),&(population->individual[i]),numDiscrete, numContinuous);

  // get back

  return 0;
}

//==============================================================

int selectionNone(Population *population, Population *parents, long M)
{
  register long i;
  long N;
  int numDiscrete,numContinuous;

  // initialize some variables

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  // allocate memory for selected individuals

  allocatePopulation(parents,M,numDiscrete, numContinuous);

  // shuffle the individuals in the population

  shufflePopulation(population);

  // copy the guys directly

  for (i=0; i<M; i++)
    copyIndividual(&(parents->individual[i]),&(population->individual[i%(population->N)]),numDiscrete, numContinuous);

  // get back

  return 0;
};

//==================================================
//==================================================

char *getSelectionDesc(int n)
{
  return selectionDesc[n];
}

//==================================================

int commonOnes(char *a, char *b, int n)
{
  register int i;
  int S=0;

  for (i=0; i<n; i++)
      S+=(a[i])&&(b[i]);

  return S;
}

//==================================================
// do divide and conquer until first M individuals
// get better or equal than the rest

int divideBest(Population *population, long left, long right, int numDiscrete, int numContinuous, long M)
{
  long l,r,lr2;
  Individual pivot;

  l = left;
  r = right;
  lr2 = (l+r)/2;

  allocateIndividual(&pivot,numDiscrete, numContinuous);

  copyIndividual(&pivot,&((population->individual)[lr2]), numDiscrete, numContinuous);

  while (l<=r)
    {
      while ((l<right)&&(moreThan(population->individual[l],pivot,numDiscrete,numContinuous))) l++;
      while ((r>left)&&(lessThan(population->individual[r],pivot,numDiscrete,numContinuous))) r--;

      if (l<=r)
	{
	  if (l!=r)
	    swapIndividuals(population,l,r);

	  l++;
	  r--;
	}
    };
  
  freeIndividual(&pivot);

  if ((l==M)||(r==(M-1)))
    return 0;

  if ((r>=M)&&(left<r))
    divideBest(population,left,r,numDiscrete,numContinuous,M);
  
  if ((l<M)&&(l<right))
    divideBest(population,l,right,numDiscrete,numContinuous,M);

  return 0;
};

//==============================================================

int selectionTournamentWithReplacementAndContinuousSharing(Population *population, Population *parents, long M)
{
  register long i,j,l;
  long N,max;
  int numDiscrete, numContinuous;
  long  *chosen;
  char generated;
  long  *id;
  long  *numerosity;
  long  maxID;

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  allocatePopulation(parents,M,numDiscrete, numContinuous);

  chosen     = (long*)   Calloc(tournSize,sizeof(long));
  id         = (long*)   Calloc(N,sizeof(long));
  numerosity = (long*)   Calloc(N,sizeof(long));
  
  maxID  = -1;

  for (i=0; i<N; i++)
    {
      id[i]=-1;
      numerosity[i]=1;
    };

  for (i=0; i<N; i++)
    if (id[i]<=0)
      {
	maxID++;
	
	id[i]=maxID;
	
	for (j=i+1; j<N; j++)
	  if (equalIndividuals(&(parents->individual[i]),&(parents->individual[j]),numDiscrete,numContinuous))
	    id[j]=maxID;
      };

  for (j=0; j<M; j++)
  {
    chosen[0]=longRand(N);
    max = chosen[0];

    for (i=1; i<tournSize; i++)
    {
       do {
	 generated=1;
	 chosen[i]=longRand(N);
	 for (int k=0; k<i; k++)
	     if (chosen[i]==chosen[k])
		generated=0;
       } while (!generated);

       if (population->individual[chosen[i]].f/numerosity[chosen[i]]>population->individual[max].f/numerosity[max])
	 max=chosen[i];	 
       
       //!!! moreThan(population->individual[chosen[i]],population->individual[max],numDiscrete,numContinuous))
       // 	 max  = chosen[i];
    }

    copyIndividual(&(parents->individual[j]),&(population->individual[max]),numDiscrete,numContinuous);
    
    // update the numerosities
    
    for (l=0; l<N; l++)
      if (id[l]==id[max])
	numerosity[l]++;
  }
  
  Free(chosen);
  Free(numerosity);
  Free(id);

  return 0;
}

//==============================================================

// int initSelectionTournamentWithReplacementAndCoevolutionarySharing(int tS, int nBM, int minH, int n, int d)
// {
//   tournSize   = tS;
//   minHamm  = minH;

//   initBusinessMenPopulation(nBM,n,d);

//   return 0;
// }

// int initBusinessMenPopulation(int nBM, int n, int d)
// {
//   int i,j;
//   char *x;
//   int good;

//   printf("Initializing BM pop. (nbm=%u, n=%u)\n",nBM,n);

//   allocatePopulation(&businessMen,nBM,n,d);

//   x = (char*) Malloc(n);

//   for (i=0; i<businessMen.N; i++)
//     {      
//       do {
// 	// randomly initialize x
	
// 	for (j=0; j<n; j++)
// 	  x[j]=flipCoin();
	
// 	for (j=0; j<n; j++)
// 	  printf("%u",x[j]);
// 	printf("\n");
// 	getchar();

// 	good = 1;
// 	for (j=0; j<i; j++)
// 	  if (hamming(x,businessMen.individual[j].chromosome[0],n)<minHamm)
// 	    good = 0;

//       } while (good==0);
      
//       for (j=0; j<n; j++)
// 	businessMen.individual[i].chromosome[0][j]=x[j];

//       printf("Got the %uth businessman.\n",i+1);
//     };

//   Free(x);

//   printf("Done initializing BM pop.\n");

//   // get back

//   return 0;
// };

// int doneBusinessMenPopulation()
// {
//   freePopulation(&businessMen);

//   return 0;
// };

// int selectionTournamentWithReplacementAndCoevolutionarySharing(Population *population, Population *parents, long M)
// {
//   register long i,j;
//   long N,max;
//   int n,d;
//   long  *chosen;
//   char generated;
//   long  *id;
//   long  *numerosity;
//   long  maxID;
//   long  min;
//   int   minDistance;

//   N = population->N;
//   n = population->n;
//   d = population->d;

//   allocatePopulation(parents,M,n,d);

//   chosen     = (long*)   Calloc(tournSize,sizeof(long));
//   id         = (long*)   Calloc(N,sizeof(long));
//   numerosity = (long*)   Calloc(businessMen.N,sizeof(long));
  
//   maxID  = -1;

//   for (j=0; j<businessMen.N; j++)
//     numerosity[j]=1;
  
//   for (i=0; i<N; i++)
//     {
//       min = 0;
//       minDistance = hamming(parents->individual[i].chromosome[0],businessMen.individual[0].chromosome[0],n);
      
//       for (j=1; j<businessMen.N; j++)
// 	if (hamming(parents->individual[i].chromosome[0],businessMen.individual[j].chromosome[0],n)<minDistance)
// 	  {
// 	    min = j;
// 	    minDistance = hamming(parents->individual[i].chromosome[0],businessMen.individual[j].chromosome[0],n);
// 	  };

//       id[i]=min;
//     };
  
//   for (j=0; j<M; j++)
//     {
//       chosen[0]=longRand(N);
//       max = chosen[0];
      
//       for (i=1; i<tournSize; i++)
// 	{
// 	  do {
// 	    generated=1;
// 	    chosen[i]=longRand(N);
// 	    for (int k=0; k<i; k++)
// 	      if (chosen[i]==chosen[k])
// 		generated=0;
// 	  } while (!generated);
	  
// 	  if (population->individual[chosen[i]].f/pow(numerosity[id[chosen[i]]],2)>population->individual[max].f/pow(numerosity[id[max]],2))
// 	    max=chosen[i];	 
// 	}
      
//       copyIndividual(&(parents->individual[j]),&(population->individual[max]),n,d);
      
//       // update the numerosities
      
//       numerosity[id[max]]++;
//     }
  
//   Free(chosen);
//   Free(numerosity);
//   Free(id);
  
//   return 0;
// }

// //==============================================================


int selectionBoltzmann(Population *population, Population *parents, long M)
{
  register long i,j;
  long N1,N;
  int numDiscrete, numContinuous;
  double *cP;
  double r;
  double variance;
  double max;

  // set the helper variables

  N  = population->N;
  N1 = N-1;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  // allocate the new population

  allocatePopulation(parents,M,numDiscrete,numContinuous);

  // compute the fitness variance

  variance = fitnessVariance(population);

  // allocate the cummulative probabilities

  cP = (double*) Calloc(N,sizeof(double));

  // compute the maximal fitness and adjust the variance so that we don't get huge numbers

  max=0;
  for (i=0; i<N; i++)
    if (population->individual[i].f>max)
      max=population->individual[i].f;
  
  if (beta*max/variance>23)
    variance=(max*beta)/23;
  
  // compute the cummulative probabilities
  
  max=0;
  for (i=0; i<N; i++)
    {
      cP[i] = exp(beta*population->individual[i].f/variance);
      if (cP[i]>max)
	max=cP[i];
    };
  printf("Max cp:   %f\n",max);
  printf("Variance: %f\n",variance);

  for (i=1; i<N; i++)
    cP[i] += cP[i-1];
  for (i=0; i<N; i++)
    cP[i] /= cP[N1];

  for (j=0; j<M; j++)
  {
    r = drand();

    for (i=0; cP[i]<r; i++);

    copyIndividual(&(parents->individual[j]),&(population->individual[i]),numDiscrete,numContinuous);
  }

  // free the cummulative probabilities

  Free(cP);

  // get back

  return 0;
};

// =======================================================================

int initSelectionBoltzmann(double new_beta)
{
  // set the beta

  beta = new_beta;

  // get back

  return 0;
};
