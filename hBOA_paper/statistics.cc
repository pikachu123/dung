#include <stdio.h>
#include <math.h>

#include "mymath.h"
#include "memalloc.h"
#include "population.h"
#include "statistics.h"
#include "distance.h"

//#define SHOW_SELECTION_INTENSITY
//#define SHOW_SELECTED_UMF
//#define LOG_DECODED_POPULATION
//#define LOG_HISTOGRAM_OF_DECODED_POPULATION
//#define LOG_BIT_SUM
//#define PRINT_DISTANCE_FROM_ORIGIN
//#define DISPLAY_ONEMAX_SIGNAL
//#define PRINT_EXTREMES
//#define STORE_EXTREMES

#define STEPS 50

//==================================================================================

int prepareFiles(StatisticsParams *statisticsParams, char *fileRoot, char *fileExtension)
{
  char s[200];

  // output of fitness info required?

  if ((statisticsParams->OutputFitness) && (fileRoot))
    {
      sprintf(s,"%s.fitness.%s",fileRoot,fileExtension);
      statisticsParams->outputFitness = fopen(s,"w");
    }
  else
    statisticsParams->outputFitness = NULL;

  // output of UMFs required?

  if ((statisticsParams->OutputUMF) && (fileRoot))
    {
      sprintf(s,"%s.UMF.%s",fileRoot,fileExtension);
      statisticsParams->outputUMF = fopen(s,"w");
    }
  else
    statisticsParams->outputUMF = NULL;

  // get back

  return 0;
};

//==================================================================================

int closeFiles(StatisticsParams *statisticsParams)
{
  // the file for fitness info open?

  if (statisticsParams->outputFitness)
    fclose(statisticsParams->outputFitness);

  // the file for UMFs open?

  if (statisticsParams->outputUMF)
    fclose(statisticsParams->outputUMF);

  // get back

  return 0;
};

//==================================================================================

int statistics(Population *population, long t, Fitness *fitness, StatisticsParams *statisticsParams)
{
  long N;
  register int l;
  int  n;
  double *p0,*p1;

  // initialize variables

  N = population->N;
  n = population->n;

  // if needed, allocate memory for univariate marginal frequencies and calculate them

  if ((statisticsParams->displayOrdering) ||
      (statisticsParams->displayGuidance) || 
      (statisticsParams->displayUMF) || 
      (statisticsParams->OutputUMF))
    {
      allocateUMF(&p0,&p1,n);
      calculateUMF(population,p0,p1);
    }

  // display epoch information and the number of fitness calls

  //comment by pikachu123 printf("Epoch (fCalc) %4li (%6lu)",t,getFitnessCalls());

  // if required, display fitness info
  
  if (statisticsParams->displayFitness)
    printf(" # fitness (worst/best/avg) = (%1.8f,%1.8f,%1.8f)",
	   population->individual[population->worst].f,
	   population->individual[population->best].f,
	   population->avgFitness);
    
  // if the function is real valued....

  if ((statisticsParams->displayFitness)&&(isRealValued()))
    {
      double f;

      printf(" # Best in (");
      
      f=realValue(population->individual[population->best].chromosome,n);
      printf(" %f",f);
      printf(")");
    }

#ifdef PRINT_DISTANCE_FROM_ORIGIN
  if (isRealValued())
    {
      double *origin=(double*) Calloc(n,sizeof(double));
      double *point=(double*) Calloc(n,sizeof(double));
      double current;
      double best = -1;

      printf(" # Best dist. from the origin ");

      for (long i=0; i<N; i++)
	{
	  for (int j=0; j<d; j++)
	    point[j]=realValue(population->individual[i].chromosome[j],n);

	  current = euclid(origin,point,n);

	  if ((best==-1)||(current<best))
	    best=current;
	};
      
      printf("%1.9f",best);

      Free(origin);
      Free(point);
    };
#endif

  // if required display number of optimal solutions found so far

  if ((fitness->isBest!=NULL)&&(statisticsParams->displayNumOptimal))
    {
      //printf(" # NumOpt = %4lu ",population->numOptimal);

#ifdef PRINT_EXTREMES
      
      if (fitness->getNumberOfOptima)
	{
	  long *numExtremes;
	  
	  numExtremes = new long[fitness->getNumberOfOptima(population->n)];

	  for (int i=0; i<fitness->getNumberOfOptima(population->n); i++)
	    numExtremes[i]=0;
	  
	  for (long i=0; i<population->N; i++)
	    {
	      long which;

	      which = fitness->whichOptimum(population->individual[i].chromosome[0],population->n);

	      if (which>=0)
		numExtremes[which]++;
	    };

	  printf("\n");
	  printf("Extremes:\n");
	  for (long i=0; i<fitness->getNumberOfOptima(population->n); i++)
 	    printf("%lu %lu\n",i,numExtremes[i]);
	  printf("\n");

	  delete numExtremes;
	};

      {
	long numOnes=0;
	long numZeros=0;

	for (i=0; i<population->N; i++)
	  {
	    int ones=0;
	    
	    for (int j=0; j<population->n; j++)
	      ones+=population->individual[i].chromosome[0][j];
	    
	    if (ones==population->n)
	      numOnes++;
	    if (ones==0)
	      numZeros++;
	  };

	printf("\n");
	printf("# 111..1 = %lu\n",numOnes);
	printf("# 000..0 = %lu\n",numZeros);
	printf("\n");
      };

#endif

#ifdef STORE_EXTREMES
      
      FILE *extremes = fopen("extremes","a");

      if (fitness->getNumberOfOptima)
	{
	  long *numExtremes;
	  
	  numExtremes = new long[fitness->getNumberOfOptima(population->numDiscrete,population->numContinuous)];

	  for (int i=0; i<fitness->getNumberOfOptima(population->numDiscrete,population->numContinuous); i++)
	    numExtremes[i]=0;
	  
	  for (long i=0; i<population->N; i++)
	    {
	      long which;

	      which = fitness->whichOptimum(population->individual[i].chromosome,population->numDiscrete,population->individual[i].continuous,population->numContinuous);

	      if (which>=0)
		numExtremes[which]++;
	    };

	  fprintf(extremes,"%li ",t);
	  for (long i=0; i<fitness->getNumberOfOptima(population->numDiscrete,population->numContinuous); i++)
 	    fprintf(extremes,"%lu ",numExtremes[i]);
	  fprintf(extremes,"\n");

	  delete numExtremes;
	};


      fclose(extremes);

#endif

    }

  // if required, display guidance of the search
  
  if (statisticsParams->displayGuidance)
    {
      double X;

      X = statisticsParams->guidanceTreshold+0.5;

      printf(" # Guided to: ");
      printf(" ");
      for (l=0; l<n; l++)
	if (p0[l]>X)
	  printf("0");
	else
	  if (p1[l]>X)
	    printf("1");
	  else
	    printf(".");
    }

  // if required, display ordering parameter

  if (statisticsParams->displayOrdering)
    {
      double ordering;

      ordering = calculateOrdering(p0,p1,n);

      printf(" # Ordering = %1.4f",ordering);
    }

  // if required, display univariate marginal frequencies

  if (statisticsParams->displayUMF)
    {
      printf("\n\nUMF: ");
      for (l=0; l<n; l++)
	printf("%1.2f ",p1[l]);
      printf("\n");
    }

  if (statisticsParams->pause)
    getchar();
  else//comment by pikchu123
   ; //printf("\n");

  // if needed, free allocated memory

  ///  printf("Checkpoint 0\n");

  if ((statisticsParams->displayOrdering)||(statisticsParams->displayGuidance)||(statisticsParams->displayUMF))
    freeUMF(&p0,&p1,n);

  // --------------------------------------------------------

  // output fitness info (to file)?

  ///  printf("Checkpoint 1\n");

  if (statisticsParams->outputFitness)
    fprintf(statisticsParams->outputFitness,
	    "%li %f %f %f\n",
	    getFitnessCalls(),
	    population->individual[population->worst].f,
	    population->individual[population->best].f,
	    population->avgFitness);
  
  // output UMFs (to file)?

  ///  printf("Checkpoint 2\n");

  if (statisticsParams->outputUMF)
    {
      ///      printf("Checkpoint 2.5\n");

      fprintf(statisticsParams->outputUMF,"%li  ",t);

      for (l=0; l<n; l++)
	fprintf(statisticsParams->outputUMF,"%1.2f ",p1[l]);
      fprintf(statisticsParams->outputUMF,"\n");
    }

#ifdef LOG_DECODED_POPULATION

  {
    FILE *out;
    char s[100];

    sprintf(s,"decodedPop.%lu",t);
    out = fopen(s,"w");
    
    for (i=0; i<population->N; i++)
      {
	char *x;
	double d;
	
	x=population->individual[i].chromosome[0];
	
	d=0;
	for (k=0; k<n; k++)
	  {
	    d *= 2;
	    d += x[k];
	  };
	
	d /= (double) ((1<<(n))-1);

	fprintf(out,"%f %f\n",d,pow(sin(5*M_PI*d),6));
      };

    fclose(out);
  };

#endif

#ifdef LOG_HISTOGRAM_OF_DECODED_POPULATION
  {
    FILE *out;
    char s[100];
    double *fq;
    int bin;

    fq=(double*)Calloc(STEPS+1,sizeof(double));
    
    sprintf(s,"histogramOfDecodedPop.%li",t);
    out = fopen(s,"w");
    
    for (i=0; i<population->N; i++)
      {
	char *x;
	double d;
	
	x=population->individual[i].chromosome[0];
	
	d=0;
	for (k=0; k<n; k++)
	  {
	    d *= 2;
	    d += x[k];
	  };
	
	d /= (double) ((1<<(n))-1);
	d += 0.5/STEPS;

	bin = (int) double(d*STEPS);
	
	fq[bin]++;
      };

    for (k=0; k<STEPS; k++)
      {
	//	fprintf(out,"%f %f\n",k*double(1)/STEPS,fq[k]/population->N);
	//	fprintf(out,"%f %f\n",(k+1)*double(1)/STEPS,fq[k]/population->N);
	fprintf(out,"%f %f\n",k*double(1)/STEPS+0.5/STEPS,fq[k]/population->N);
      };

    fclose(out);
    Free(fq);

  };

#endif

#ifdef DISPLAY_ONEMAX_SIGNAL

  {
    double *f1,*f0;
    long n0,n1;
    double f;
    int l;
    double max,min;

    f1=(double*) Calloc(n,sizeof(double));
    f0=(double*) Calloc(n,sizeof(double));
    
    for (int k=0; k<n; k++)
      {
	f1[k]=0;
	f0[k]=0;
	n0=0;
	n1=0;
	
	for (i=0; i<population->N; i++)
	  {
	    if (population->individual[i].chromosome[0][k]==1)
	      {
		f1[k]+=population->individual[i].f;
		n1++;
	      }
	    else
	      {
		f0[k]+=population->individual[i].f;
		n0++;
	      };
	  };
	
	if (n0*n1>0)
	  {
	    f0[k]/=n0;
	    f1[k]/=n1;
	  }
	else
	  {
	    f0[k]=-1;
	    f1[k]=-1;
	  };
      };

    printf("\n\n");
    printf("Signal: ");
    for (int k=0; k<n; k++)
      if (f0[k]!=-1)
	printf("%1.2f ",f1[k]-f0[k]);
      else
	printf("<conv> ");
    f=0;l=0;
    max=-1E+08;
    min=1E+08;
    for (int k=0; k<n; k++)
      if (f0[k]!=-1)
	{
	  f+=f1[k]-f0[k];
	  l++;

	  if (f1[k]-f0[k]<min)
	    min=f1[k]-f0[k];

	  if (f1[k]-f0[k]>max)
	    max=f1[k]-f0[k];
	};
    f/=l;
    printf("\nAverage signal: %1.2f\n",f);
    printf("Min: %1.2f\n",min);
    printf("Max: %1.2f\n",max);
    getchar();
    
    Free(f0);
    Free(f1);

  }

#endif

#ifdef LOG_BIT_SUM

  {
    FILE *out;
    char s[100];

    sprintf(s,"bitSum.%li",t);
    out = fopen(s,"w");
    
    for (i=0; i<population->N; i++)
      {
	char *x;
	int d;
	
	x=population->individual[i].chromosome[0];
	
	d=0;
	for (k=0; k<n; k++)
	    d += x[k];

	fprintf(out,"%u\n",d);
      };

    fclose(out);
  };

#endif

//   printf("  best guy: ");
//   printIndividual(stdout,&(population->individual[population->best]),population->numDiscrete,population->numContinuous);
 
  ///  printf("Checkpoint 3\n");

  // get back

  return 0;
}

//==================================================================================

int runStatistics(Population *population, 
		  Fitness *fitness,
		  long generations,
		  double epsilon, 
		  double maxOptimal, 
		  int maxFailures,
		  RunLog *runLog, int run)
{
  long  N;
  int numDiscrete,numContinuous;
  double *p0, *p1;
  register int j;
  char converged;

  // initialize variables

  N = population->N;
  numDiscrete  = population->numDiscrete;
  numContinuous = population->numContinuous;

  (*runLog)[run].numDiscrete = numDiscrete;
  (*runLog)[run].numContinuous = numContinuous;

  // allocate memory for and calculate univariate marginal frequencies

  allocateUMF(&p0,&p1,numDiscrete);
  calculateUMF(population,p0,p1);

  // fill in things in run log array

  allocateIndividual(&((*runLog)[run].best),numDiscrete,numContinuous);
  allocateIndividual(&((*runLog)[run].worst),numDiscrete,numContinuous);

  copyIndividual(&((*runLog)[run].best),&(population->individual[population->best]),numDiscrete,numContinuous);
  copyIndividual(&((*runLog)[run].worst),&(population->individual[population->worst]),numDiscrete,numContinuous);

  // compute the number of good BBs in the best guy (wrt to number of good BBs)

  if (fitness->goodBBs)
    {
      double best;
      double active;
      
      best = population->individual[0].goodBBs;
      for (long i=1; i<population->N; i++)
	{
	  active = population->individual[i].goodBBs;
	  if (active>best)
	    best=active;
	};
      
      (*runLog)[run].goodBBs = best;
      /*comment by pikachu123
      printf("Best is            %lu\n",population->best);
      printf("BBs of the best:   %f\n",fitness->goodBBs(population->individual[population->best].chromosome,numDiscrete,population->individual[population->best].continuous,numContinuous));
      printf("which is supp. :   %f\n",population->individual[population->best].goodBBs);
      printf("BBs of the winner: %f\n",best);
      */
    }
  else
    (*runLog)[run].goodBBs = -1;

  (*runLog)[run].avgFitness   = population->avgFitness;
  (*runLog)[run].ordering     = calculateOrdering(p0,p1,numDiscrete);

  if (fitness->isBest==NULL)
    {
      (*runLog)[run].bestDefined = 0;
      (*runLog)[run].bestFound   = 0;
    }
  else
    {
      (*runLog)[run].bestDefined = 1;
      (*runLog)[run].bestFound   = isBestIndividual(&((*runLog)[run].best),fitness,numDiscrete,numContinuous,HEAVY);
      if ((*runLog)[run].bestFound)
	(*runLog)[run].bestFoundIn = getFitnessCalls();
      else
	(*runLog)[run].bestFoundIn = 0;
    };
  
  converged = 1;
  for (j=0; j<numDiscrete; j++)
    if ((p1[j]>=epsilon)&&(p0[j]>=epsilon))
      converged=0;

  if (converged)
    {
      Individual x;

      allocateIndividual(&x,numDiscrete,numContinuous);

      for (j=0; j<numDiscrete; j++)
	if (p1[j]>0.5)
	  x.chromosome[j]=1;
	  else
	    x.chromosome[j]=0;

      (*runLog)[run].converged = 1;
      (*runLog)[run].timeToConvergence = getFitnessCalls();
      (*runLog)[run].generationsToConvergence = generations;

      if (fitness->isBest)
	{
	  (*runLog)[run].convergedToTheBest = isBestIndividual(&x,fitness,numDiscrete,numContinuous,HEAVY);
	  (*runLog)[run].failed = !(*runLog)[run].convergedToTheBest;
	}  
      else
	(*runLog)[run].convergedToTheBest = 0;
      
      freeIndividual(&x);
    }
  else
    {
      (*runLog)[run].converged = 0;
      (*runLog)[run].convergedToTheBest = 0;
      (*runLog)[run].failed = 1;
      
      if (fitness->isBest)
	if (isBestIndividual(&(population->individual[population->best]),fitness,numDiscrete,numContinuous,HEAVY))
	  (*runLog)[run].failed=0;
    }
  
  if (((*runLog)[run].bestFound)&&(population->numOptimal==0))
    population->numOptimal=1;
    

  ///  printf("Failed check.... -> %u\n",(*runLog)[run].failed);
  // reached max optimal?
  
  if (maxOptimal>0)
    if ((double)population->numOptimal/N>=maxOptimal)
      {
	(*runLog)[run].maxOptimalReached = 1;
	(*runLog)[run].maxOptimalReachedIn = getFitnessCalls();
	(*runLog)[run].failed = 0;
      }
    else
      {
        printf("Did not reach number of optima (%lu/%lu<%f)!!!\n",population->numOptimal,N,maxOptimal);
	(*runLog)[run].maxOptimalReached = 0;
	(*runLog)[run].maxOptimalReachedIn = 0;
	(*runLog)[run].failed = 1;
      }

   // print gotten results out to standard output
//printf("Results of the gen # %3u\n",generations);
//comment by pikachu123
/*
  printf("---i
   fitnessDefinition.fitness=myfit2;
   fitnessDefinition.goodBBs=NULL;
   fitnessDefinition.isBest=mybest;----------------------------------------------------------\n");
  printf("Results of the run # %3u\n",run);
  printf("\n");

  printf("Best fitness         : %f\n",(*runLog)[run].best.f);
  printf("Worst fitness        : %f\n",(*runLog)[run].worst.f);
  printf("Avg. fitness         : %f\n",(*runLog)[run].avgFitness);

  if ((*runLog)[run].goodBBs!=-1)
    printf("Good BBs             : %f\n",(*runLog)[run].goodBBs);

  if ((*runLog)[run].bestDefined)
    printf("Found optimum?       : %s\n",((*runLog)[run].bestFound)? "Yes":"No");

  printf("Converged?           : %s\n",((*runLog)[run].converged)? "Yes":"No");

  if (((*runLog)[run].bestDefined) && ((*runLog)[run].converged))
    printf("Converged to optimum : %s\n",((*runLog)[run].convergedToTheBest)? "Yes":"No");
  if (maxOptimal>0)
    printf("Max optimal reached  : %s\n",((*runLog)[run].maxOptimalReached)? "Yes":"No");

  printf("--------------------------------------------------------------\n");
*/
  // free memory used by univariate marginal frequencies

  freeUMF(&p0,&p1,numDiscrete);

  //  get back

  return 0;
}

//==================================================================================

int finalStatistics(RunLog *runLog, int numRuns, FILE *outputFile)
{
  register int i;
  
  double bestBest,worstBest,avgBest,devBest;
  double bestWorst,worstWorst,avgWorst,devWorst;
  double bestAvg,worstAvg,avgAvg,devAvg;

  double avgGoodBBs,devGoodBBs;

  double success,failure,convergence,successfulConvergence;

  double avgTimeToConvergence,avgTimeToSuccessfulConvergence;
  double devTimeToConvergence,devTimeToSuccessfulConvergence;
  double avgMaxOptimalReachedIn,devMaxOptimalReachedIn;

  double avgGenerationsToConvergence,devGenerationsToConvergence;
  double avgGenerationsToSuccessfulConvergence,devGenerationsToSuccessfulConvergence;

  double maxOptimalReached;

  FILE *output;

  // calculate all data

  bestBest   = (*runLog)[0].best.f;
  worstBest  = (*runLog)[0].best.f;
  avgBest    = (*runLog)[0].best.f;
  bestWorst  = (*runLog)[0].worst.f;
  worstWorst = (*runLog)[0].worst.f;
  avgWorst   = (*runLog)[0].worst.f;
  bestAvg    = (*runLog)[0].avgFitness;
  worstAvg   = (*runLog)[0].avgFitness;
  avgAvg     = (*runLog)[0].avgFitness;
  avgGoodBBs = (*runLog)[0].goodBBs;

  success    = 0;
  failure    = 0;
  successfulConvergence = convergence = 0;
  maxOptimalReached = avgMaxOptimalReachedIn = 0;
  double avgBestFoundIn = 0;

  if ((*runLog)[0].bestDefined)
    if ((*runLog)[0].bestFound)
      {
	avgBestFoundIn += (*runLog)[0].bestFoundIn;
	success++;
      }
    else
      failure++;

  if ((*runLog)[0].converged)
    {
      convergence++;
      avgTimeToConvergence = (*runLog)[0].timeToConvergence;
      avgGenerationsToConvergence = (*runLog)[0].generationsToConvergence;
    }
  else
    {
      avgTimeToConvergence = 0;
      avgGenerationsToConvergence = 0;
    }

  if ((*runLog)[0].convergedToTheBest)
    {
      successfulConvergence++;
      avgTimeToSuccessfulConvergence = (*runLog)[0].timeToConvergence;
      avgGenerationsToSuccessfulConvergence = (*runLog)[0].generationsToConvergence;
    }
  else
    {
      avgTimeToSuccessfulConvergence = 0;
      avgGenerationsToSuccessfulConvergence = 0;
    }

  if ((*runLog)[0].maxOptimalReached)
    {
      maxOptimalReached++;
      avgMaxOptimalReachedIn = (*runLog)[0].maxOptimalReachedIn;
    }

  for (i=1; i<numRuns; i++)
    {
      if ((*runLog)[i].best.f>bestBest)
	bestBest = (*runLog)[i].best.f;
      else
	if ((*runLog)[i].best.f<worstBest)
	  worstBest = (*runLog)[i].best.f;

      avgBest += (*runLog)[i].best.f;

      if ((*runLog)[i].worst.f>bestWorst)
	bestWorst = (*runLog)[i].worst.f;
      if ((*runLog)[i].worst.f<worstWorst)
	worstWorst = (*runLog)[i].worst.f;

      avgWorst += (*runLog)[i].worst.f;

      if ((*runLog)[i].avgFitness>bestAvg)
	bestAvg = (*runLog)[i].avgFitness;
      else
	if ((*runLog)[i].avgFitness<worstAvg)
	  worstAvg = (*runLog)[i].avgFitness;

      avgAvg += (*runLog)[i].avgFitness;

      avgGoodBBs += (*runLog)[i].goodBBs;

      if ((*runLog)[i].bestDefined)
	if ((*runLog)[i].bestFound)
	  {
	    avgBestFoundIn += (*runLog)[i].bestFoundIn;
	    success++;
	  }
	else
	  failure++;

      if ((*runLog)[i].converged)
	{
	  convergence++;
	  avgTimeToConvergence += (*runLog)[i].timeToConvergence;
	  avgGenerationsToConvergence += (*runLog)[i].generationsToConvergence;
	}

      if ((*runLog)[i].convergedToTheBest)
	{
	  successfulConvergence++;
	  avgTimeToSuccessfulConvergence += (*runLog)[i].timeToConvergence;
	  avgGenerationsToSuccessfulConvergence += (*runLog)[i].generationsToConvergence;
	}

      if ((*runLog)[i].maxOptimalReached)
	{
	  maxOptimalReached++;
	  avgMaxOptimalReachedIn += (*runLog)[i].maxOptimalReachedIn;
	}
    }

  avgBest  /= numRuns;
  avgWorst /= numRuns;
  avgAvg   /= numRuns;

  avgGoodBBs /= numRuns;

  if (success)
    avgBestFoundIn /= success;

  if (convergence>0)
    {
      avgTimeToConvergence        /= convergence;
      avgGenerationsToConvergence /= convergence;
    };

  if (successfulConvergence)
    {
      avgTimeToSuccessfulConvergence        /= successfulConvergence;
      avgGenerationsToSuccessfulConvergence /= successfulConvergence;
    };

  if (maxOptimalReached>0)
    avgMaxOptimalReachedIn /= maxOptimalReached;
  else
    avgMaxOptimalReachedIn = 0;  

  devBest = devWorst = devAvg = 0;

  devGoodBBs = 0;

  devTimeToConvergence = devTimeToSuccessfulConvergence = 0;
  devGenerationsToConvergence = devGenerationsToSuccessfulConvergence = 0;

  devMaxOptimalReachedIn = 0;

  double devBestFoundIn = 0;

  for (i=0; i<numRuns; i++)
    {
      devBest  += fsqr((*runLog)[i].best.f-avgBest);
      devWorst += fsqr((*runLog)[i].worst.f-avgWorst);
      devAvg   += fsqr((*runLog)[i].avgFitness-avgAvg);

      devGoodBBs += fsqr((*runLog)[i].goodBBs-avgGoodBBs);

      if ((*runLog)[i].bestFound)
	devBestFoundIn += fsqr((*runLog)[i].bestFoundIn-avgBestFoundIn);

      if ((*runLog)[i].converged)
	{
	  devTimeToConvergence += dsqr((*runLog)[i].timeToConvergence-avgTimeToConvergence);
	  devGenerationsToConvergence += dsqr((*runLog)[i].generationsToConvergence-avgGenerationsToConvergence);
	};

      if ((*runLog)[i].convergedToTheBest)
	{
	  devTimeToSuccessfulConvergence += dsqr((*runLog)[i].timeToConvergence-avgTimeToSuccessfulConvergence);
	  devGenerationsToSuccessfulConvergence += dsqr((*runLog)[i].generationsToConvergence-avgGenerationsToSuccessfulConvergence);
	};

      if ((*runLog)[i].maxOptimalReached)
	devMaxOptimalReachedIn += dsqr((*runLog)[i].maxOptimalReachedIn-avgMaxOptimalReachedIn);
    }

  if (numRuns>1)
    {
      devBest  = sqrt(devBest/(numRuns-1));
      devWorst = sqrt(devWorst/(numRuns-1));
      devAvg   = sqrt(devAvg/(numRuns-1));
      devGoodBBs = sqrt(devGoodBBs/(numRuns-1));
      if (success*numRuns>1)
	devBestFoundIn = sqrt(devBestFoundIn/(success*numRuns-1));
      else
	devBestFoundIn = 0;
    }
  else
    devBest = devWorst = devAvg = devGoodBBs = 0;

  if (convergence>1)
    {
      devTimeToConvergence = sqrt(devTimeToConvergence/(convergence-1));
      devGenerationsToConvergence = sqrt(devGenerationsToConvergence/(convergence-1));
    }
  else
    {
      devTimeToConvergence = 0;
      devGenerationsToConvergence = 0;
    };
      

  if (successfulConvergence>1)
    {
      devTimeToSuccessfulConvergence = sqrt(devTimeToSuccessfulConvergence/(successfulConvergence-1));
      devGenerationsToSuccessfulConvergence = sqrt(devGenerationsToSuccessfulConvergence/(successfulConvergence-1));
    }
  else
    {
      devTimeToSuccessfulConvergence = 0;
      devGenerationsToSuccessfulConvergence = 0;
    };

  if (maxOptimalReached>1)
    devMaxOptimalReachedIn = sqrt(devMaxOptimalReachedIn/(maxOptimalReached-1));
  else
    devMaxOptimalReachedIn = 0;

  success  /= numRuns;
  success  *= 100;

  failure  /= numRuns;
  failure  *= 100;

  convergence /= numRuns;
  convergence *= 100;

  successfulConvergence /= numRuns;
  successfulConvergence *= 100;

  maxOptimalReached /= numRuns;
  maxOptimalReached *= 100;

  // output them all to standard output and output file (if any)

  output=stdout;
  do
    {
      fprintf(output,"==================================================\n");
      fprintf(output,"FINAL STATISTICS\n");
      fprintf(output,"\n");
      fprintf(output,"Number of runs          : %u\n",numRuns);
      fprintf(output,"\n");
      fprintf(output,"Best fitness ever       : %f\n",bestBest);
      fprintf(output,"Worst of best found     : %f\n",worstBest);
      fprintf(output,"Avg. best fitness       : %f\n",avgBest);
      if (numRuns>1)
	fprintf(output,"Dev. of best fitness    : %f\n",devBest);
      fprintf(output,"\n");
      fprintf(output,"Best of worse found     : %f\n",bestWorst);
      fprintf(output,"Worst fitness ever      : %f\n",worstWorst);
      fprintf(output,"Avg. worst fitness      : %f\n",avgWorst);
      if (numRuns>1)
	fprintf(output,"Dev. of worst fintess   : %f\n",devWorst);
      fprintf(output,"\n");
      fprintf(output,"Best average fitness    : %f\n",bestAvg);
      fprintf(output,"Worst average fitness   : %f\n",worstAvg);
      fprintf(output,"Avg. average fitness    : %f\n",avgAvg);
      if (numRuns>1)
	fprintf(output,"Dev. of average fitness : %f\n",devAvg);
      fprintf(output,"\n");
      if ((*runLog)[0].bestDefined)
	{
	  fprintf(output,"Found optimum in         : %6.2f",success);fprintf(output," %s\n","%");
	  fprintf(output,"Failed to find it in     : %6.2f",failure);fprintf(output," %s\n","%");
	  fprintf(output,"\n");

	  fprintf(output,"Reached num. of optimal solutions in : %1.2f",maxOptimalReached);fprintf(output,"%s\n","%");
	  fprintf(output,"Average time to reach this was       : %f\n",avgMaxOptimalReachedIn);
	  fprintf(output,"Std. Dev. of time to reach it        : %f\n",devMaxOptimalReachedIn);
	  fprintf(output,"\n");
	}

      fprintf(output,"Found optimum in               : %1.2f",success);fprintf(output,"%s\n","%");
      fprintf(output,"Average time to find it was    : %f\n",avgBestFoundIn);
      fprintf(output,"Std. Dev. of time to find it   : %f\n",devBestFoundIn);
      fprintf(output,"\n");
      
      if (avgGoodBBs>0)
	{
	  fprintf(output,"Avg proportion of good BB: %f\n",avgGoodBBs);
	  if (numRuns>1)
	    fprintf(output,"Dev. of prop. of good BBs: %f\n",devGoodBBs);
	  fprintf(output,"\n");
	};

      fprintf(output,"----------------------------------\n");
      fprintf(output,"Converged in             : %1.2f",convergence);fprintf(output," %s\n","%");

      if (convergence)
	{
	  fprintf(output,"\nFitness evaluations until convergence\n");
	  fprintf(output,"Average  : %f\n",avgTimeToConvergence);
	};

      if (convergence>1)
	fprintf(output,"Std dev. : %f\n",devTimeToConvergence);

      if (convergence)
	{
	  fprintf(output,"\nGenerations until convergence\n");
	  fprintf(output,"Average  : %f\n",avgGenerationsToConvergence);
	};

      if (convergence>1)
	fprintf(output,"Std dev. : %f\n",devGenerationsToConvergence);

      if ((*runLog)[0].bestDefined)
	{
	  fprintf(output,"\n");
	  fprintf(output,"----------------------------------\n");
	  fprintf(output,"Converged to optimum in  : %1.2f",successfulConvergence);fprintf(output," %s\n","%");

	  if (successfulConvergence)
	    {
	      fprintf(output,"\nFitness evaluations until successful convergence\n");
	      fprintf(output,"Average  : %f\n",avgTimeToSuccessfulConvergence);
	    };

	  if (successfulConvergence>1)
	    fprintf(output,"Std dev  : %f",devTimeToSuccessfulConvergence);

	  fprintf(output,"\n");

	  if (successfulConvergence)
	    {
	      fprintf(output,"\nGenerations until successful convergence\n");
	      fprintf(output,"Average  : %f\n",avgGenerationsToSuccessfulConvergence);
	    };

	  if (successfulConvergence>1)
	    fprintf(output,"Std dev  : %f",devGenerationsToSuccessfulConvergence);

	  fprintf(output,"\n");
  
	}

      fprintf(output,"==================================================\n\n");

      if (output==stdout)
	{
	  output = outputFile;
	}
      else
	{
	  output = NULL;
	};
    } while (output!=NULL);
     
  // get back

  return 0;
}

//==================================================================================

int initRunLog(RunLog *runLog, int numRuns)
{
  // allocate the memory needed for the run log

  (*runLog) = (RunLog) Calloc(numRuns, sizeof(RunLogItem));

  // get back

  return 0;
}

//==================================================================================

int doneRunLog(RunLog *runLog, int numRuns)
{
  int i;

  // free the memory used by best and worst guys in the run logs

  for (i=0; i<numRuns; i++)
    {
      freeIndividual(&((*runLog)[i].best));
      freeIndividual(&((*runLog)[i].worst));
    }

  // free the memory used by the run log
  
  Free(*runLog);
  
  // get back
  
  return 0;
}

//==================================================================================

int selectionStatistics(Population *old, Population *selected)
{

#ifdef SHOW_SELECTION_INTENSITY
  {
    long i;
    double oldFAvg;
    double selectionIntensity;
    double oldFStd;
    double selectedFAvg;
    
    oldFAvg = 0;
    for (i=0; i<old->N; i++)
      oldFAvg += old->individual[i].f;
    oldFAvg /= old->N;
    
    oldFStd = 0;
    for (i=0; i<old->N; i++)
      oldFStd += fsqr(oldFAvg-old->individual[i].f);
    oldFStd = sqrt(oldFStd/(old->N-1));
    
    selectedFAvg = 0;
    for (i=0; i<selected->N; i++)
      selectedFAvg += selected->individual[i].f;
    selectedFAvg /= selected->N;
    
    printf("Old population avg. fitness = %1.3f\n",oldFAvg);
    printf("Sel population avg. fitness = %1.3f\n",selectedFAvg);
    
    selectionIntensity = (selectedFAvg-oldFAvg)/oldFStd;
    
    printf("Selection Intensity = %1.3f\n",selectionIntensity);
  }
#endif 

#ifdef SHOW_SELECTED_UMF

  {
    int k,l;
    double **p0,**p1;

    allocateUMF(&p0,&p1,selected->n,selected->d);
    calculateUMF(selected,p0,p1);
    
    printf("\n\nUMF of the selected:\n");
    for (k=0; k<selected->d; k++)
      {
	printf(" ");
	for (l=0; l<selected->n; l++)
	  printf("%1.2f ",p1[k][l]);
      }
    printf("\n");
  }

#endif

  return 0;
} 
