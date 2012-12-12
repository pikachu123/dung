#include "mymath.h"
#include "memalloc.h"
#include "random.h"
#include "population.h"
#include "distance.h"
#include "binary.h"

#define NUM_SAMPLES 20

//#define PRINT_CLUSTERS

#define CENTER_TOURNAMENT_SIZE 1

//======================================================================

int allocatePopulation(Population *population, long N, int numDiscrete, int numContinuous)
{
  register long i;

  population->N = N;
  population->n = numDiscrete+numContinuous;
  population->numDiscrete = numDiscrete;
  population->numContinuous = numContinuous;
 
  // allocate memory for array of individuals of length N and dxn chromosome for each of them

  population->individual = (Individual*) Calloc(N,sizeof(Individual));

  for (i=0; i<N; i++)
    allocateIndividual(&(population->individual[i]),numDiscrete,numContinuous);

  population->evaluated = 0;
  population->numOptimal = 0;
  population->clusterSize = NULL;

  // get back
  
  return 0;
}

//======================================================================

int freePopulation(Population *population)
{
  long N;
  int  n;

  N = population->N;
  n = population->n;

  // free from the population all the memory allocated by allocatePopulation

  for (int i=0; i<N; i++)
    freeIndividual(&(population->individual[i]));

  Free(population->individual);
  population->N=0;

  if (population->clusterSize)
    Free(population->clusterSize);

  // get back

  return 0;
}

//======================================================================

int generatePopulation(Population *population)
{
  register long i;

  long N;
  int  numDiscrete,numContinuous;

  N = population->N;
  numDiscrete = population->numDiscrete;
  numContinuous = population->numContinuous;

  // randomly generate all individuals

  for (i=0; i<N; i++)
    generateIndividual(&(population->individual[i]),numDiscrete,numContinuous);

  population->evaluated = 0;
  population->numOptimal = 0;

  // get back

  return 0;
}

//======================================================================

int evaluatePopulation(Population *population, Fitness *fitness)
{
  long i;
  long N,best,worst;
  int numDiscrete,numContinuous;
  double bestF,worstF;
  double sum;
  
  // initialize some variables

  N = population->N;
  numDiscrete = population->numDiscrete;
  numContinuous = population->numContinuous;

  population->numOptimal = 0;

  // evaluate just the first individual

  evaluateIndividual(&(population->individual[0]),fitness,numDiscrete,numContinuous);
  
  // if best increase counter

  if (isBestIndividual(&(population->individual[0]),fitness,numDiscrete,numContinuous,LIGHT))
    population->numOptimal++;
  
  // prepare best and worst individuals identifiers

  best=0;
  worst=0;
  sum = bestF = worstF = population->individual[0].f;
  
  // for the rest of individuals do...
  
  for (i=1; i<N; i++)
    {
      // evaluate the individual

      evaluateIndividual(&(population->individual[i]),fitness,numDiscrete,numContinuous);

      // if optimal increase counter

      if (isBestIndividual(&(population->individual[i]),fitness,numDiscrete,numContinuous,LIGHT))
	population->numOptimal++;
      
      // update sum of fitnesses

      sum += population->individual[i].f;

      // update best and worst individuals

      if (population->individual[i].f>bestF)
	{
	  best  = i;
	  bestF = population->individual[i].f;
	}
      else
	if (population->individual[i].f<worstF)
	  {
	    worst  = i;
	    worstF = population->individual[i].f;
	  }
     }

  // set the best and worst chromosomes, and average fitness
  
  population->best       = best;
  population->worst      = worst;
  population->avgFitness = sum/N;

  population->evaluated  = 1;

  // get back

  return 0;
}

//======================================================================

int reevaluatePopulation(Population *population, Fitness *fitness)
{
  long i;

  printf("\nReevaluating parents\n");

  for (i=0; i<population->N; i++)
    population->individual[i].fCalculated=0;

  evaluatePopulation(population,fitness);

  return 0;
};
//======================================================================

int recomputeFitnessInfo(Population *population, Fitness *fitness)
{
  long i;
  long N,best,worst;
  int numDiscrete,numContinuous;
  double bestF,worstF;
  double sum;
  
  // initialize some variables

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  population->numOptimal = 0;
  
  // if best increase counter

  if (isBestIndividual(&(population->individual[0]),fitness,numDiscrete,numContinuous,LIGHT))
    population->numOptimal++;

  // prepare best and worst individuals identifiers

  best=0;
  worst=0;
  sum = bestF = worstF = population->individual[0].f;
  
  // for the rest of individuals do...

  for (i=1; i<N; i++)
    {
      // if optimal increase counter

      if (isBestIndividual(&(population->individual[i]),fitness,numDiscrete,numContinuous,LIGHT))
	population->numOptimal++;

      // update sum of fitnesses

      sum += population->individual[i].f;
      
      // update best and worst individuals

      if (population->individual[i].f>bestF)
	{
	  best  = i;
	  bestF = population->individual[i].f;
	}
      else
	if (population->individual[i].f<worstF)
	  {
	    worst  = i;
	    worstF = population->individual[i].f;
	  }
     }

  // set the best and worst chromosomes, and average fitness

  population->best       = best;
  population->worst      = worst;
  population->avgFitness = sum/N;

  // get back

  return 0;
}

//===========================================================

int swapIndividuals(Population *population, long i, long j)
{
  Individual tmp;

  allocateIndividual(&tmp,population->numDiscrete,population->numContinuous);
  copyIndividual(&tmp,&(population->individual[i]),population->numDiscrete,population->numContinuous);
  copyIndividual(&(population->individual[i]),&(population->individual[j]),population->numDiscrete,population->numContinuous);
  copyIndividual(&(population->individual[j]),&tmp,population->numDiscrete,population->numContinuous);
  freeIndividual(&tmp);
  
  if (population->evaluated)
    {
      if (population->best==i)
	population->best=j;
      if (population->worst==i)
	population->worst=j;
    }

  return 0;
}

//===========================================================

int replaceIndividual(Population *population, long i, Individual *individual)
{
  population->evaluated = 0;
  
  copyIndividual(&(population->individual[i]),individual,population->numDiscrete,population->numContinuous);

  return 0;
};

//===========================================================

int printPopulation(FILE *output, Population *population)
{
  register long i;
  long N;
  int numDiscrete,numContinuous;

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  for (i=0; i<N; i++)
    printIndividual(output, &(population->individual[i]), numDiscrete, numContinuous);
  
  return 0;
}

//===========================================================

int shufflePopulation(Population *population)
{
  long N,i,j;

  // set the variables

  N = population->N;

  // shuffle the individuals a little
  
  for (i=0; i<N; i++)
    {
      j = (long) ((double) drand()*N);

      if (i!=j)
	swapIndividuals(population,i,j);
    };

  // get back

  return 0;
}

//================================================================

int allocateUMF(double **p0, double **p1, int n)
{
  if (n>0)
    {
      *p0 = (double*) Calloc(n,sizeof(double));
      *p1 = (double*) Calloc(n,sizeof(double));
    };

  return 0;
}

//=================================================================

int freeUMF(double **p0, double **p1, int n)
{
  if (n>0)
    {
      Free(*p0);
      Free(*p1);
    };

  return 0;
}

//=================================================================

int calculateUMF(Population *population, double *p0, double *p1)
{
  long i;
  int l;
  long N;
  int  n;

  N = population->N;
  n = population->numDiscrete;

  for(l=0; l<n; l++)
    p1[l]=0;

  for (i=0; i<N; i++)
    for(l=0; l<n; l++)
      p1[l]+=population->individual[i].chromosome[l];

  for(l=0; l<n; l++)
    {
      p1[l] /= N;
      p0[l] = 1-p1[l];
    }

  return 0;
}

//=================================================================

int calculateProportionateUMF(Population *population, double *p0, double *p1)
{
  register long i;
  register int l;
  long N;
  int  n;
  long n0,n1;
  double f0,f1;

  N = population->N;
  n = population->n;

  for(l=0; l<n; l++)
    p1[l]=0;
  
  for(l=0; l<n; l++)
    {
      n0=n1=0;
      f0=f1=0;

      for (i=0; i<N; i++)
	if (population->individual[i].chromosome[l]==1)
	  {
	    n1++;
	    f1+=population->individual[i].f;
	  }
	else
	  {
	    n0++;
	    f0+=population->individual[i].f;
	  };
      
      if (n0==0)
	{
	  p0[l]=1;
	  p1[l]=0;
	}
	else
	  if (n1==0)
	    {
	      p0[l]=0;
	      p1[l]=1;
	    }
	  else
	    {
	      f0 /= n0;
	      f1 /= n1;
	      
	      f0 = exp(1*f0);
	      f1 = exp(1*f1);
	
	      p0[l] = f0/(f0+f1);
	      p1[l] = f1/(f0+f1);
	    };
    };
 
  return 0;
}

//=================================================================

double calculateOrdering(double *p0, double *p1, int n)
{
  double o;

  o=0;

  for (int l=0; l<n; l++)
    o+=fsqr(p1[l]-0.5);

  o*=4;
  o/=n;

  return o;
}

//=================================================================

int UMFCloserThanEpsilon(double *p0, double *p1, int n, double epsilon)
{
  int isCloser;

  isCloser=1;

  for (int l=0; l<n; l++)
    if ((p0[l]>=epsilon)&&(p1[l]>=epsilon))
      isCloser=0;

  return isCloser;
}

//==================================================================

int allocateBMF(double ***p00, double ***p01, double ***p10, double ***p11, int n)
{
  (*p00) = (double**) Calloc(n,sizeof(double*));
  (*p01) = (double**) Calloc(n,sizeof(double*));
  (*p10) = (double**) Calloc(n,sizeof(double*));
  (*p11) = (double**) Calloc(n,sizeof(double*));
  
  for (int j=0; j<n; j++)
    {
      (*p00)[j] = (double*) Calloc(n,sizeof(double));
      (*p01)[j] = (double*) Calloc(n,sizeof(double));
      (*p10)[j] = (double*) Calloc(n,sizeof(double));
      (*p11)[j] = (double*) Calloc(n,sizeof(double));
    }

  return 0;
}

//==================================================================

int freeBMF(double ***p00, double ***p01, double ***p10, double ***p11, int n)
{
  for (int j=1; j<n; j++)
    {
      Free((*p00)[j]);
      Free((*p01)[j]);
      Free((*p10)[j]);
      Free((*p11)[j]);
    }

  Free(*p00);
  Free(*p01);
  Free(*p10);
  Free(*p11);

  return 0;
}

//==================================================================

int calculateBMF(Population *population, double **p00, double **p01, double **p10, double **p11)
{
  register long i;
  register int k,l;
  long N;
  int n;
  char valueK,valueL;
  long n00,n01,n10,n11;

  N = population->N;
  n = population->n;
  
  for (k=1; k<n; k++)
    for (l=0; l<k; l++)
      {
	n00=n01=n10=n11=0;
	
	for (i=0; i<N; i++)
	  {
	    //	      printf("(%lu %u %u %u)\n",i,j,k,l);
	    valueK = population->individual[i].chromosome[k];
	    valueL = population->individual[i].chromosome[l];
	    // 	      valueK = getDigit(&(population->individual[i]),j,k);
	    // 	      valueL = getDigit(&(population->individual[i]),j,l);
	    if (valueK)
	      if (valueL)
		n11++;
	      else
		n10++;
	    else
	      if (valueL)
		n01++;
	  }
	
	n00 = N-n01-n10-n11;
	
	p00[k][l] = p00[l][k] = (double) n00/N;
	p01[k][l] = p10[l][k] = (double) n01/N;
	p10[k][l] = p01[l][k] = (double) n10/N;
	p11[k][l] = p11[l][k] = (double) n11/N;
      }
  
  return 0;
}

//==================================================================

int calculateClusterUMF(Population *p, int numClusters, double **p0, double **p1)
{
  int i,j,n;
  long l;
  long **n0,**n1;
  char *x=NULL;

  // assign some variables

  n=p->n;

  // allocate the memory and initialize

  n0 = (long**) Calloc(numClusters,sizeof(long*));
  n1 = (long**) Calloc(numClusters,sizeof(long*));

  for (i=0; i<numClusters; i++)
    {
      n0[i]=(long*) Calloc(n,sizeof(long));
      n1[i]=(long*) Calloc(n,sizeof(long));

      for (j=0; j<n; j++)
	n0[i][j]=n1[i][j]=0;
    };

  // compute the counts

  for (l=0; l<p->N; l++)
    {
      x=p->individual[l].chromosome;
      
      for (i=0; i<p->n; i++)
	if (x[i]==0)
	  n0[p->individual[l].cluster][i]++;
	else
	  n1[p->individual[l].cluster][i]++;
    };
 
  // compute the probabilities

  for (i=0; i<numClusters; i++)
    for (j=0; j<n; j++)
      {
	p0[i][j]=double(n0[i][j])/(n0[i][j]+n1[i][j]);
	p1[i][j]=1-p0[i][j];
      };

  // free the memory

  for (i=0; i<numClusters; i++)
    {
      Free(n0[i]);
      Free(n1[i]);
    };

  Free(n0);
  Free(n1);

  // get back

  return 0;
};

//==================================================================

int kMeansClusterPopulation(Population *population, int k, int numRestarts, char phenotypic)
{
  int i,j,pass;
  int cluster=0;
  int oldcluster=0;
  long l;
  int changed;

  double **center=NULL;
  double **bestCenter=NULL;
  double *phenotypeCenter=NULL;
  double *bestPhenotypeCenter=NULL;
  double  maxPhenotypeCenter;
  double *phenotype=NULL;

  double distance, dummy, totalDistance, bestDistance=0;

  char *x=NULL;

  //   printf("Clustering...\n");

  // allocate memory for cluster sizes

  if (population->clusterSize)
    Free(population->clusterSize);

  population->clusterSize   = (long*) Calloc(k,sizeof(long));

  // allocate memory for centers

  if (!phenotypic)
    {
      center = (double**) Calloc(k,sizeof(double*));
      bestCenter = (double**) Calloc(k,sizeof(double*));
      for (i=0; i<k; i++)
	{
	  center[i] = (double*) Calloc(population->n,sizeof(double));
	  bestCenter[i] = (double*) Calloc(population->n,sizeof(double));
	};
    }
  else
    {
      phenotypeCenter = (double*) Calloc(k,sizeof(double));
      bestPhenotypeCenter = (double*) Calloc(k,sizeof(double));

      maxPhenotypeCenter=1;
      for (i=1; i<population->n; i++)
	{
	  maxPhenotypeCenter*=2;
	  maxPhenotypeCenter++;
	};

      phenotype = (double*) Calloc(population->N,sizeof(double));
      for (l=0; l<population->N; l++)
	phenotype[l] = binaryToDouble(population->individual[l].chromosome,population->n)/maxPhenotypeCenter;
    };

  // all but the first cluster sizes are 0

  for (i=0; i<k; i++)
    population->clusterSize[i]=0;
  population->clusterSize[0]=population->N;

  // all guys into the first cluster

  for (l=0; l<population->N; l++)
    population->individual[l].cluster = 0;

  // perform numRestart passes

  for (pass=0; pass<numRestarts; pass++)
    {
      
      //      printf("Starting pass %2u\n",pass);

      // generate centers
      
      if (!phenotypic)
	{
	  for (i=0; i<k; i++)
	    {
	      l = longRand(population->N);
	      for (j=0; j<population->n; j++)
		center[i][j]=population->individual[l].chromosome[j];
	    };
	}
      else
	{	  
 	  for (i=0; i<k; i++)
 	    {
// 	      phenotypeCenter[i]=0;

// 	      for (j=0; j<1; j++)
// 		{
// 		  l=longRand(population->N);
// 		  phenotypeCenter[i] += phenotype[l];
// 		};

// 	      phenotypeCenter[i]/=1;
	      
 	      l = longRand(population->N);
 	      phenotypeCenter[i] = drand();//phenotype[l];
	      
// 	      if (i==0)
// 		{
// 		  guy = longRand(population->N);
// 		  phenotypeCenter[i] = phenotype[guy];
// 		}
// 	      else
// 		{
// 		  max=-1;
// 		  for (m=0; m<CENTER_TOURNAMENT_SIZE; m++)
// 		    {
// 		      l=longRand(population->N);
// 		      bestDistance = fabs(phenotype[l]-phenotypeCenter[0]);
// 		      for (j=1; j<i; j++)
// 			{
// 			  dummy = fabs(phenotype[l]-phenotypeCenter[j]);
// 			  if (dummy<bestDistance)
// 			    bestDistance=dummy;
// 			};
		      
// 		      if (bestDistance>max)
// 			{
// 			  max=bestDistance;
// 			  guy=l;
// 			};
// 		    };

// 		  phenotypeCenter[i]=phenotype[guy];
// 		};
	    };
	};
      
      do {
	
	// choose a cluster for each guy and count the changes, updating the centers on the fly
	
	changed       = 0;
	totalDistance = 0;

	for (l=0; l<population->N; l++)
	  {
	    //	    printf("   assigning individual %lu with phenotype %f\n",l,phenotype[l]);

	    if (!phenotypic)
	      {
		x          = population->individual[l].chromosome;
		oldcluster = population->individual[l].cluster;
		
		// find out the closest cluster center
		
		cluster=0;
		distance=euclidBinaryDouble(x,center[0],population->n);
		
		for (i=1; i<k; i++)
		  {
		    dummy = euclidBinaryDouble(x,center[i],population->n);
		    if (dummy<distance)
		      {
			cluster=i;
			distance=dummy;
		      };
		  };
	      }
	    else
	      {
		oldcluster = population->individual[l].cluster;

		cluster=0;
		distance=fabs(phenotypeCenter[0]-phenotype[l]);

		for (i=1; i<k; i++)
		  {
		    dummy = fabs(phenotypeCenter[i]-phenotype[l]);
		    if (dummy<distance)
		      {
			cluster = i;
			distance = dummy;
		      };
		  };
	      };

	    // cluster changed?

	    if (cluster!=oldcluster)
	      {
		// update the sizes

		population->clusterSize[oldcluster]--;
		population->clusterSize[cluster]++;

		// change the cluster information in the individual
	    
		changed++;
		population->individual[l].cluster = cluster;
	      };

	    totalDistance += distance;
	  };

	// recompute the centers

	if (!phenotypic)
	  {
	    for (i=0; i<k; i++)
	      for (j=0; j<population->n; j++)
		center[i][j]=0;
	    
	    for (l=0; l<population->N; l++)
	      {
		x          = population->individual[l].chromosome;
		oldcluster = population->individual[l].cluster;
		
		for (j=0; j<population->n; j++)
		  center[oldcluster][j] += x[j];
	      };

	    for (i=0; i<k; i++)
	      {
		for (j=0; j<population->n; j++)
		  {
		    if (population->clusterSize[i]>0)
		      center[i][j]/=population->clusterSize[i];
		    else
		      center[i][j]=0;
		  };
	      };
	  }
	else
	  {
	    for (i=0; i<k; i++)
	      phenotypeCenter[i]=0;

	    for (l=0; l<population->N; l++)
	      phenotypeCenter[population->individual[l].cluster] += phenotype[l];

	    for (i=0; i<k; i++)
	      {
		phenotypeCenter[i]/=population->clusterSize[i];
		//		printf("   center[%u] = %f   (%4lu guys)\n",i,phenotypeCenter[i],population->clusterSize[i]);
	      }
	  };

	//	printf("Clustering...changes = %u\n",changed);
	//	getchar();

      } while (changed>0);
  
      //      printf("New distance: %f\n",totalDistance);

      if ((pass==0)||(totalDistance<bestDistance))
	{
	  if (!phenotypic)
	    {
	      for (i=0; i<k; i++)
		for (j=0; j<population->n; j++)
		  bestCenter[i][j]=center[i][j];
	    }
	  else
	    for (i=0; i<k; i++)
	      bestPhenotypeCenter[i] = phenotypeCenter[i];

	  bestDistance=totalDistance;
	};
      
      //      printf("Pass %2u: distance = %f\n",pass,totalDistance);
    };

  // cluster the guys according to the best seen so far

  for (i=0; i<k; i++)
    population->clusterSize[i]=0;

  for (l=0; l<population->N; l++)
    {
      if (!phenotypic)
	{
	  x          = population->individual[l].chromosome;
	  
	  // find out the closest cluster center
	  
	  cluster=0;
	  distance=euclidBinaryDouble(x,bestCenter[0],population->n);
	  
	  for (i=1; i<k; i++)
	    {
	      dummy = euclidBinaryDouble(x,bestCenter[i],population->n);
	      if (dummy<distance)
		{
		  cluster=i;
		  distance=dummy;
		};
	    };
	}
      else
	{
	  cluster = 0;
	  distance=fabs(bestPhenotypeCenter[0]-phenotype[l]);
	  
	  for (i=1; i<k; i++)
	    {
	      dummy = fabs(bestPhenotypeCenter[i]-phenotype[l]);
	      if (dummy<distance)
		{
		  cluster = i;
		  distance = dummy;
		};
	    };
	};

      population->individual[l].cluster = cluster;
      population->clusterSize[cluster]++;
    };
  
#ifdef PRINT_CLUSTERS

  if (!phenotypic)
    for (i=0; i<k; i++)
      {
 	printf("Cluster %2u\n",i);
 	printf("   Size   : %lu\n",population->clusterSize[i]);
 	printf("   Center : ");
	for (j=0; j<population->n; j++)
	  printf(" %f",bestCenter[i][j]);
	printf("\n\n");
      }
  else
    for (i=0; i<k; i++)
      printf("   center[%u] = %f   (%4lu guys)\n",i,bestPhenotypeCenter[i],population->clusterSize[i]);
  //getchar();
#endif
    
    // free memory used by the centers
  
    if (!phenotypic)
      {
	for (i=0;i<k; i++)
	  {
	    Free(center[i]);
	    Free(bestCenter[i]);
	  };
	
	Free(center);
	Free(bestCenter);
      }
    else
      {
	Free(phenotypeCenter);
	Free(bestPhenotypeCenter);
	Free(phenotype);
      };

    //   printf("Leaving clustering.\n");
  
    // get back
  
    return 0;
};

// ==============================================================

double fitnessVariance(Population *p)
{
  long   i;
  long   N;
  double sum;
  double mean;
  double var;

  // helper variables

  N = p->N;

  // compute the mean first

  sum=0;
  for (i=0; i<N; i++)
    sum+=p->individual[i].f;
  mean = sum/N;

  // now compute the variance

  sum=0;
  for (i=0; i<N; i++)
    sum+=fsqr(p->individual[i].f-mean);
  var=sum/(N-1);

  // and return it

  return var;
};


// ==============================================================

void swapPopulations(Population *p, Population *q)
{
  Population aux;

  aux = *p;
  *p  = *q;
  *q  = aux;
};
