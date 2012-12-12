#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "binary.h"
#include "discretization.h"
#include "random.h"

//#define DEBUG
//#define RANDOM_UNDISCRETIZE
#define STORE_HISTOGRAM

#ifdef DEBUG
#define debug(statement) {statement;fflush(stdout);}
#else
#define debug(statement)
#endif

Discretization::Discretization()
{
  continuousPopulation = discretePopulation = NULL;
  bitsPerVariable=0;
  numBins=0;
  numContinuous=numDiscrete=0;
  index=NULL;
  lo=hi=NULL;
  clusterSize=NULL;
  index2=NULL;
};

Discretization::~Discretization()
{
  if (lo)
    freeBounds();
  
  if (clusterSize)
    {
      for (int i=0; i<numContinuous; i++)
	free(clusterSize[i]);
      free(clusterSize);
    };
  
  if (index)
    {
      for (int i=0; i<numContinuous; i++)
	free(index[i]);
      free(index);
    };
};

void Discretization::allocateBounds()
{
  if (lo)
    freeBounds();

  lo = (double**) calloc(numContinuous,sizeof(double*));
  hi = (double**) calloc(numContinuous,sizeof(double*));

  for (int i=0; i<numContinuous; i++)
    {
      lo[i] = (double*) calloc(numBins,sizeof(double));
      hi[i] = (double*) calloc(numBins,sizeof(double));
    };
};

void Discretization::freeBounds()
{
  for (int i=0; i<numContinuous; i++)
    {
      free(lo[i]);
      free(hi[i]);
    };

  free(lo);
  free(hi);
};

void Discretization::discretizeFHH()
{
  // allocate memory for the bounds

  debug(printf("    - going to allocate the bounds\n"));
  allocateBounds();
 
  // now let's compute the bounds
  debug(printf("    - going to compute the bounds\n"));
  computeBoundsFHH(continuousPopulation);

  // now we must allocate the target population and create it, really

  debug(printf("    - going to allocate the target population\n"));
  allocatePopulation(discretePopulation,continuousPopulation->N,continuousPopulation->numContinuous*bitsPerVariable,0);

  // now we must copy the discrete variables and fit the rest depending on the values
  // of continuous variables and the created histogram

  // the rest of the discrete variables is created by discretization itself

  for (long i=0; i<continuousPopulation->N; i++)
    {
      copyDiscreteVariables(discretePopulation->individual[i].chromosome,
			    continuousPopulation->individual[i].chromosome,
			    continuousPopulation->numDiscrete);
      discretizeContinuousVariablesFHH(discretePopulation->individual[i].chromosome+continuousPopulation->numDiscrete,
				       continuousPopulation->individual[i].continuous,
				       continuousPopulation->numContinuous);
    };
};


void Discretization::discretizeFWH()
{
  if (index2)
    {
      for (int i=0; i<numContinuous; i++)
	{
	  for (int j=0; j<numBins; j++)
	    if (index2[i][j])
	      free(index2[i][j]);
	  free(index2[i]);
	};
      free(index2);
    };
  
  index2 = (long***) calloc(numContinuous,sizeof(long**));
  for (int i=0; i<numContinuous; i++)
    {
      index2[i] = (long**) calloc(numBins,sizeof(long*)); // will finish allocating as we know the sizes of clusters
    };
 
  if (clusterSize)
    {
      for (int i=0; i<numContinuous; i++)
	free(clusterSize[i]);
      free(clusterSize);
    };

  clusterSize = (long**) calloc(numContinuous,sizeof(long*));
  for (int i=0; i<numContinuous; i++)
    clusterSize[i] = (long*) calloc(numBins,sizeof(long));

  allocatePopulation(discretePopulation,continuousPopulation->N,continuousPopulation->numContinuous*bitsPerVariable,0);

  for (long i=0; i<continuousPopulation->N; i++)
    copyDiscreteVariables(discretePopulation->individual[i].chromosome,
			  continuousPopulation->individual[i].chromosome,
			  continuousPopulation->numDiscrete);

  for (int k=0; k<numContinuous; k++)
    {
      for (int bin=0; bin<numBins; bin++)
	clusterSize[k][bin]=0;

      for (long i=0; i<continuousPopulation->N; i++)
	{
	  double d=0;
	  double dStep=double(1)/numBins;
	  int bin=0;
	  
	  while (d+dStep<continuousPopulation->individual[i].continuous[k])
	    {
	      d+=dStep;
	      bin++;
	    };

	  if (bin>=numBins)
	    bin=numBins-1;
	  
	  clusterSize[k][bin]++;
	  continuousPopulation->individual[i].cluster=bin;
	  
	  intToBinary(discretePopulation->individual[i].chromosome+continuousPopulation->numDiscrete+k*bitsPerVariable,bin,bitsPerVariable);
	};

      for (int i=0; i<numBins; i++)
	if (clusterSize[k][i]>0)
	  index2[k][i]=(long*)calloc(clusterSize[k][i],sizeof(long));
	else
	  index2[k][i]=NULL;

      ///      printf("Cluster size of %u is %u\n",9,clusterSize[k][9]);
      
      for (int i=0; i<numBins; i++)
	clusterSize[k][i]=0;
      
      for (long i=0; i<continuousPopulation->N; i++)
	{
	  int cluster=continuousPopulation->individual[i].cluster;
	  ///	  printf("%u %u %u\n",k,cluster,clusterSize[k][cluster]);
	  index2[k][cluster][clusterSize[k][cluster]]=i;
	  clusterSize[k][cluster]++;
	};
    };
};


int doubleCompare(const void *a, const void *b)
{
  double x = *((double*)a);
  double y = *((double*)b);

  if (x>y)
    return 1;
  else
    if (x<y)
      return -1;
    else
      return 0;
}

void Discretization::discretizeKMeans()
{
  debug(printf("Inside discretizeKMeans\n"););
  
  if (index2)
    {
      for (int i=0; i<numContinuous; i++)
	{
	  for (int j=0; j<numBins; j++)
	    if (index2[i][j])
	      free(index2[i][j]);
	  free(index2[i]);
	};
      free(index2);
    };
  
  index2 = (long***) calloc(numContinuous,sizeof(long**));
  for (int i=0; i<numContinuous; i++)
    {
      index2[i] = (long**) calloc(numBins,sizeof(long*)); // will finish allocating as we know the sizes of clusters
    };
 
  double *centers = (double*) calloc(numBins,sizeof(double));

  if (clusterSize)
    {
      for (int i=0; i<numContinuous; i++)
	free(clusterSize[i]);
      free(clusterSize);
    };

  clusterSize = (long**) calloc(numContinuous,sizeof(long*));
  for (int i=0; i<numContinuous; i++)
    clusterSize[i] = (long*) calloc(numBins,sizeof(long));

  allocatePopulation(discretePopulation,continuousPopulation->N,continuousPopulation->numContinuous*bitsPerVariable,0);

  for (int k=0; k<numContinuous; k++)
    {
      for (int i=0; i<numBins; i++)
	centers[i]=drand();

      qsort(centers,numBins,sizeof(double),&doubleCompare);
    
      cluster(centers,clusterSize[k],k,numBins);

      for (int i=0; i<numBins; i++)
	if (clusterSize[k][i]>0)
	  {
	    index2[k][i] = (long*) calloc(clusterSize[k][i],sizeof(long));
	    clusterSize[k][i]=0;
	  }
	else
	  index2[k][i]=NULL;

      for (long i=0; i<continuousPopulation->N; i++)
	{
	  int cluster = continuousPopulation->individual[i].cluster;
	  index2[k][cluster][clusterSize[k][cluster]++]=i;
	  intToBinary(discretePopulation->individual[i].chromosome+continuousPopulation->numDiscrete+k*bitsPerVariable,cluster,bitsPerVariable);
	};
    };
  
  for (long i=0; i<continuousPopulation->N; i++)
    {
      copyDiscreteVariables(discretePopulation->individual[i].chromosome,
			    continuousPopulation->individual[i].chromosome,
			    continuousPopulation->numDiscrete);
    };

  free(centers);
  
  debug(printf("Outside discretizeKMeans\n"););
};

void Discretization::cluster(double *centers, long *centerSize, int k, int numClusters)
{
  double *newCenters = (double*) calloc(numClusters,sizeof(double));
  long    change=continuousPopulation->N;
  long    previousChanges=0;

  for (long i=0; i<continuousPopulation->N; i++)
    continuousPopulation->individual[i].cluster=-1;

  while (change>0) 
    {
      change=0;

      for (int j=0; j<numClusters; j++)
	{
	  centerSize[j]=0;
	  newCenters[j]=0;
	};
    
      for (long i=0; i<continuousPopulation->N; i++)
	{
	  int    minIndex = continuousPopulation->individual[i].cluster;
	  if (minIndex<0) minIndex=0;
	  double min = fabs(centers[minIndex]-continuousPopulation->individual[i].continuous[k]);
	
	  for (int j=0; j<numBins; j++)
	    {
	      double d = fabs(centers[j]-continuousPopulation->individual[i].continuous[k]);
	  
	      if (d<min)
		{
		  min      = d;
		  minIndex = j;
		};
	    };

	  if (continuousPopulation->individual[i].cluster!=minIndex)
	    {
	      /// had to add this other condition on the differences, points were repeatedly migrating sometimes (rarely but still)
	      double difference = fabs(fabs(continuousPopulation->individual[i].continuous[k]-centers[minIndex])-fabs(continuousPopulation->individual[i].continuous[k]-centers[continuousPopulation->individual[i].cluster]));
	      if (centers[continuousPopulation->individual[i].cluster]!=centers[minIndex])
		difference /= fabs(centers[continuousPopulation->individual[i].cluster]-centers[minIndex]);
	      
	      if (difference>1E-5)
		{
		  change++;
#ifdef DEBUG
		  if (previousChanges==1)
		    {
		      printf("Changing from %u to %u\n",continuousPopulation->individual[i].cluster,minIndex);
		      printf("Centers: %6.3f %6.3f\n",centers[continuousPopulation->individual[i].cluster],centers[minIndex]);
		      printf("Point: %6.3f\n",continuousPopulation->individual[i].continuous[k]);
		      printf("Difference: %f\n",fabs(continuousPopulation->individual[i].continuous[k]-centers[minIndex])-fabs(continuousPopulation->individual[i].continuous[k]-centers[continuousPopulation->individual[i].cluster]));
		    };
#endif
		};
	    };
	  
	  centerSize[minIndex]++;
	  continuousPopulation->individual[i].cluster = minIndex;
      
	  newCenters[minIndex] += continuousPopulation->individual[i].continuous[k];
	};
    
      for (int j=0; j<numClusters; j++)
	newCenters[j] /= centerSize[j];

#ifdef DEBUG
      printf("Old center: ");
      for (int j=0; j<numClusters; j++)
	{
	  printf("%f (%lu), ",newCenters[j],centerSize[j]);
	};
      printf("\n");
      printf("New center: ");
      for (int j=0; j<numClusters; j++)
	{
	  printf("%f (%lu), ",newCenters[j],centerSize[j]);
	};
      printf("\n");
#endif
    
      memcpy(centers,newCenters,numClusters);

      debug(printf("Number of changes = %3lu\n",change););
      previousChanges=change;
#ifdef DEBUG
      printf("New centers: ");
      for (int j=0; j<numClusters; j++)
	printf("%5.2f ",newCenters[j]);
      printf("\n");
#endif
    };
  
  free(newCenters);
};

void Discretization::discretize(Population *p, Population *q, DiscretizationParams *discretizationParams)
{
  continuousPopulation=p;
  discretePopulation=q;

  debug(printf("-> Discretize entered.\n")); 

  // first assign the important variables 

  numDiscrete = p->numDiscrete;
  numContinuous = p->numContinuous;
  bitsPerVariable = discretizationParams->bitsPerVariable;
  numBins = 1L<<bitsPerVariable;
  
  switch (discretizationParams->discretizationType) 
    {
    case 1: 
      discretizeFHH();
      break;

    case 2:
      discretizeKMeans();
      break;

    case 3:
      discretizeFWH();
      break;

    default: 
      fprintf(stderr,"ERROR: Unknown discretization method!\n");
      exit(-1);
    };
  

/*   printf("\n\nDiscretized population:\n"); */
/*   printPopulation(stdout,discretePopulation); */
/*   getchar(); */

  // all done

  return;
};


void Discretization::copyDiscreteVariables(char *dest, char *src, int n)
{
  memcpy(dest,src,n);
};

void Discretization::discretizeContinuousVariablesFHH(char *dest, double *x, int n)
{
  for (int k=0; k<n; k++)
    {
      int which;

      for (which=0; (hi[k][which]<x[k]); which++);

      intToBinary(dest+k*bitsPerVariable,which,bitsPerVariable);
      
/*       printf("Translated %u = ",which); */
/*       for (int i=0; i<bitsPerVariable; i++) */
/* 	printf("%u",(dest+k*bitsPerVariable)[i]? 1:0); */
/*       printf("\n");getchar(); */
    };
};

void Discretization::computeBoundsFHH(Population *p)
{
  if (index)
    {
      for (int i=0; i<numContinuous; i++)
	free(index[i]);
      free(index);
    };

  index = (long**) calloc(numContinuous,sizeof(long*));
  for (int i=0; i<numContinuous; i++)
    index[i] = (long*) calloc(p->N,sizeof(long));

  binSize = p->N/numBins;

  debug(printf("       : computing the bounds, ready to iterate\n"));

  for (int k=0; k<numContinuous; k++)
    {
      debug(printf("       : iteration %u\n",k));
      debug(printf("       : setting up index\n"));

      for (long i=0; i<p->N; i++)
	index[k][i]=i;

      debug(printf("       : setting up compare\n"));

      setCompare(p,k);

      debug(printf("       : sorting\n"));

      qsort(index[k],p->N,sizeof(long),&discretizationCompare);

/*       for (long i=0; i<p->N; i++) */
/* 	{ */
/* 	  printf("%f\n",p->individual[index[k][i]].continuous[k]); */
/* 	}; */
/*       getchar(); */

      debug(printf("       : computing the bounds finally\n"));

      // compute the bounds now

      for (int bin=0; bin<numBins; bin++)
	{
	  if (bin>0)
	    lo[k][bin]=hi[k][bin-1];
	  else
	    lo[k][bin]=0; //!! might do this differently

	  if (bin<numBins-1)
	    {
	      long i1 = index[k][binSize*(bin+1)-1];
	      long i2 = index[k][binSize*(bin+1)];

	      double x = continuousPopulation->individual[i1].continuous[k];
	      double y = continuousPopulation->individual[i2].continuous[k];

	      hi[k][bin] = 0.5*(x+y);
	    }
	  else
	    hi[k][bin]=1; //!! might do this differently

	  debug(printf("         : new interval (%f,%f)\n",lo[k][bin],hi[k][bin]));
	};
    };
  
#ifdef STORE_HISTOGRAM
  FILE *out = fopen("hist","w");
  for (int bin=0; bin<numBins; bin++)
    {
	{
	  fprintf(out,"%f %f\n",hi[0][bin],0.0);
	  fprintf(out,"%f %f\n",hi[0][bin],1.0);
	  fprintf(out,"%f %f\n",hi[0][bin],0.0);
	};
    };
  fclose(out);

  out = fopen("hist2","w");
  for (int bin=0; bin<numBins; bin++)
    {
	{
	  fprintf(out,"%f %f\n",0.0,hi[1][bin]);
	  fprintf(out,"%f %f\n",1.0,hi[1][bin]);
	  fprintf(out,"%f %f\n",0.0,hi[1][bin]);
	};
    };
  fclose(out);
#endif
  
};

Population *comparePopulation;
long       *compareIndex;
int         compareWhichVariable;

void Discretization::setCompare(Population *p, int whichVariable)
{
  comparePopulation = p;
  compareWhichVariable = whichVariable;
};

int discretizationCompare(const void *x, const void *y)
{
  long i= *((long*)x);
  long j= *((long*)y);

  double a=comparePopulation->individual[i].continuous[compareWhichVariable];
  double b=comparePopulation->individual[j].continuous[compareWhichVariable];

  if (a>b)
    return +1;
  else
    if (a<b)
      return -1;
    else
      return 0;
};


void Discretization::undiscretizeFHH(Individual *from, Individual *to)
{
  for (int k=0; k<numContinuous; k++)
    {
      char *x;
      
      x = from->chromosome+numDiscrete+k*bitsPerVariable;
	  
      int which = binaryToInt(x,bitsPerVariable);
	  
#ifdef RANDOM_UNDISCRETIZE
      to->continuous[k] = lo[k][which]+(hi[k][which]-lo[k][which])*drand();
#else
      long l = binSize*which+longRand(binSize);
      
      /* 	  printf("Binsize: %u\n",binSize); */
      /* 	  printf("popsize: %u\n",continuousPopulation->N); */
      /* 	  printf("which:   %u\n",which); */
      /* 	  printf("picked:  %u\n",l); */
      
      /* 	  printf("[%u] => %u\n",k,which); */
      
      if (l>=continuousPopulation->N)
	l=continuousPopulation->N-1;
      
      l = index[k][l];
      
      debug(printf("Original param: %f (%f %f)\n",continuousPopulation->individual[l].continuous[k],lo[k][which],hi[k][which]));
      
      to->continuous[k] = continuousPopulation->individual[l].continuous[k];
      to->mutation[k]   = continuousPopulation->individual[l].mutation[k];
#endif
    };
};

void Discretization::undiscretizeKMeans(Individual *from, Individual *to)
{
  for (int k=0; k<numContinuous; k++)
    {
      char *x;
      
      x = from->chromosome+numDiscrete+k*bitsPerVariable;
   
      int which = binaryToInt(x,bitsPerVariable);
      
      long l = longRand(clusterSize[k][which]);

      if (clusterSize[k][which]==0)
	l=longRand(continuousPopulation->N);
      else
	l = index2[k][which][l];

      to->continuous[k] = continuousPopulation->individual[l].continuous[k];
      to->mutation[k]   = continuousPopulation->individual[l].mutation[k];
    };
};

void Discretization::undiscretize(Population *p,Population *q, DiscretizationParams *params)
{
  debug(printf("   -> Undiscretization called.\n"));
  debug(printf("        : going to allocate memory for the new population\n"));
  
  allocatePopulation(q,p->N,numDiscrete,numContinuous);

  debug(printf("        : ready to iterate\n"));

/*   printf("\n\nDiscrete population:\n"); */
/*   printPopulation(stdout,p); */
/*   getchar(); */

  for (long i=0; i<p->N; i++)
    {
      debug(printf("        : copying discrete variables\n"));

      copyDiscreteVariables(q->individual[i].chromosome,p->individual[i].chromosome,numDiscrete);

      debug(printf("        : going to generate continuous variables\n"));

      // !!! must put a switch here
      
      switch (params->discretizationType)
	{
	case 1:
	  undiscretizeFHH(&(p->individual[i]),&(q->individual[i]));
	  break;

	case 2:
	  undiscretizeKMeans(&(p->individual[i]),&(q->individual[i]));
	  break;

	case 3:
	  undiscretizeKMeans(&(p->individual[i]),&(q->individual[i]));//fwh is the same
	  break;
	  
	default:
	  fprintf(stderr,"ERROR: Unknown discretization type.\n");
	  exit(-1);
	};

      double tau=4;
      double coef = exp(tau*gaussianRandom(0,1)/sqrt(numContinuous));
      
      for (int k=0; k<numContinuous; k++)
	{
	  q->individual[i].mutation[k]*=coef;
	  q->individual[i].continuous[k] += gaussianRandom(0,q->individual[i].mutation[k]);

	  if (q->individual[i].continuous[k]>1)
	    q->individual[i].continuous[k]=1;
	  if (q->individual[i].continuous[k]<0)
	    q->individual[i].continuous[k]=0;
	};
    };
};
