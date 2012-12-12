#include <stdio.h>
#include <string.h>

#include "graph.h"
#include "random.h"
#include "population.h"
#include "memalloc.h"
#include "select.h"
#include "recombination.h"
#include "replace.h"
#include "statistics.h"
//#include "mda.h"
#include "priors.h"
#include "discretization.h"

#define randomAllele randomDigit

/* #define STORE_PARENTS */
/* #define STORE_CHILDREN */
/* #define STORE_CUMULATIVE_VAR1 */
//#define DEBUG

#ifdef DEBUG
#define debug(statement) statement;
#else
#define debug(statement)
#endif

//=====================================================================

int done(Population *population, 
	 Fitness *fitness, 
	 long t, 
	 long tMax, 
	 long maxFC, 
	 double epsilon, 
	 char stopWhenFoundOptimum, 
	 double maxAverageFitness,
	 double maxOptimal)
{
  int finito;

  // ----------------------------

  finito = (t>=tMax);

  // ----------------------------

  if (!finito)
      finito = ((maxFC!=-1) && (getFitnessCalls()>=maxFC));

  // ----------------------------

  if (!finito)
    finito = (stopWhenFoundOptimum) && (isBestIndividual(&(population->individual[population->best]),fitness,population->numDiscrete,population->numContinuous,HEAVY));
  
  // ----------------------------

  if ((!finito)&&(epsilon>=0))
    {
      double *p0,*p1;

      allocateUMF(&p0,&p1,population->numDiscrete);

      calculateUMF(population,p0,p1);

      finito=UMFCloserThanEpsilon(p0,p1,population->numDiscrete,epsilon);

      freeUMF(&p0,&p1,population->numDiscrete);
    }

  if ((maxAverageFitness>0)&&(maxAverageFitness<=population->avgFitness))
    finito = 1;

  if (maxOptimal>0)
    {
      ///      printf("NumOptimal = %u\n",population->numOptimal);
      if (((double)population->numOptimal)/population->N>=maxOptimal)
	finito = 1;
    }

  return finito;
}

//======================================================================

int ea( long N,
	long offspringSize,
	long M,
	int numDiscrete, 
	int numContinuous,

	DiscretizationParams *discretizationParams,

  	long *tPerformed,
        long tMax, 
	  
        Fitness *fitness, 

	long   maxFC,
	double epsilon,
	char   stopWhenFoundOptimum,
	double maxAverageFitness,
	double maxOptimal,

        Selection *selection,

	Recombination *recombination,
	RecombinationParams *recombinationParams,
	
	Replacement *replace,

	StatisticsParams *statisticsParams,

	char *fileRoot,
	char *fileExtension,

	Population *theLastPopulation)
{
  long t;
  Population population;
  Population parents;
  Population children;
  Discretization *discretization = new Discretization();
  Population tmpPopulation;

  // prepare output files

  prepareFiles(statisticsParams,fileRoot,fileExtension);

  // reset the number of fitness calls

  resetFitnessCalls();

  // create initial population (allocate, generate, and evaluate)

  debug(printf("Going to allocate the population...\n"));

  allocatePopulation(&population, N, numDiscrete, numContinuous);
  generatePopulation(&population);
  //printf("size:%d\n",N);

  if ((recombinationParams->injectGoodGuys)&&(fitness->injectTheGood))
    fitness->injectTheGood((void*)&population);

  debug(printf("Going to evaluate the population...\n"));

  evaluatePopulation(&population, fitness); 

  debug(printf("Going to iterate the EA...\n"));

  // do it!

  t=-1;

  while (!done(&population,fitness,t,tMax,maxFC,epsilon,stopWhenFoundOptimum,maxAverageFitness,maxOptimal))
    {

      ///      printf("Population:\n");
      ///      printPopulation(stdout,&population);

      // output information about population

      debug(printf("Going to compute the statistic stuff...\n"););

      statistics(&population,t,fitness,statisticsParams);

      // select the set of selected parents

      debug(printf("Going to select the parents.\n");)

      (*selection)(&population,&parents,M);

#ifdef STORE_PARENTS
      FILE *tmpFile = fopen("parents","w");
      printPopulation(tmpFile,&parents);
      fclose(tmpFile);
#endif

#ifdef STORE_CUMULATIVE_VAR1
      FILE *tmpFile3 = fopen("cumm1","w");
      double *x      = (double*) calloc(parents.N,sizeof(double));

      for (int i=0; i<parents.N; i++)
	x[i] = parents.individual[i].continuous[1];

      qsort(x,parents.N,sizeof(double),&doubleCompare);

      for (int i=0; i<parents.N; i++)
	fprintf(tmpFile3,"%f %f\n",x[i],double(i)/parents.N);

      fclose(tmpFile3);
      free(x);
#endif

      if (discretizationParams->discretizationType>0)
	{
	  debug(printf(" -> I am about to discretize\n"));

	  discretization->discretize(&parents,&tmpPopulation,discretizationParams);
	  
	  swapPopulations(&parents,&tmpPopulation);
	};

      // selection statistics

      selectionStatistics(&population,&parents);

      // recombination the population

#ifdef DEBUG
      printf("Before recombination:\n");  
/*        printPopulation(stdout,&parents);  */
/*        getchar();  */
#endif

       (*recombination->recombination)(&parents,&children,offspringSize,recombinationParams);

#ifdef DEBUG
      printf("After recombination:\n");  
/*         printPopulation(stdout,&children);  */
/*         getchar();  */
#endif

      if (discretizationParams->discretizationType>0)
	{
	  Population tmpPopulation2;
	  swapPopulations(&parents,&tmpPopulation);

	  debug(printf(" -> I am about to undiscretize\n"));

	  // this undiscretizes children, fixes up everything about discretization

	  discretization->undiscretize(&children,&tmpPopulation2,discretizationParams);

	  debug(printf(" -> undiscretized.\n"));

	  swapPopulations(&children,&tmpPopulation2);
	  swapPopulations(&parents,&tmpPopulation);
	  
	  freePopulation(&tmpPopulation);
	  freePopulation(&tmpPopulation2);
	};

      debug(printf("Recombined the selected.\n"));

#ifdef STORE_CHILDREN
      FILE *tmpFile2 = fopen("children","w");
      printPopulation(tmpFile2,&children);
      fclose(tmpFile2);
#endif

      // evaluate the children

      evaluatePopulation(&children,fitness);
      
      // if noisy, we must evaluate the parents too
      
      if ((getFitnessNoiseVariance()>0)&&
	  (parents.N!=children.N))
	reevaluatePopulation(&parents,fitness);
      
      debug(printf("Evaluated the children.\n"));

      // add offspring to the population

      (*replace)(&population,&children);

      debug(printf("Replaced the old.\n"));

      // recompute the best, worst, and average fitnesses

      recomputeFitnessInfo(&population, fitness);

      debug(printf("Recomputed the fitness information.\n"));

      // free the children population 	  

      freePopulation(&children);

      debug(printf("Freed the children population.\n"));

      // free the parent population

      freePopulation(&parents);

      // increase the epoch

      t++;

      debug(printf("Epoch done\n"));
    }

  // output final statistics over the final population

  statistics(&population,t,fitness,statisticsParams);

  // put the last population into theLastPopulation variable serving for run statistics later

  memcpy(theLastPopulation,&population,sizeof(Population));

  // close output files (if any open)

  closeFiles(statisticsParams);

  // the number of generations performed is equal to the current value of t

  *tPerformed = t;

  // get rid of the discretization

  delete discretization;

  // get back
  int fdopt=0;
  if(done(&population,fitness,t,tMax,maxFC,epsilon,stopWhenFoundOptimum,maxAverageFitness,maxOptimal))
     fdopt=1;
  //return 0;
    return fdopt;
}
