#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include <unistd.h>
#include <time.h>

//#include <sys/resource.h>
#include "utils.h"
#include "getFileArgs.h"
#include "args.h"
#include "population.h"
#include "statistics.h"
#include "select.h"
#include "recombination.h"
#include "replace.h"
#include "ea.h"
#include "random.h"
#include "boa.h"
#include "memalloc.h"
#include "mymath.h"
#include "priors.h"
#include "decisionGraphBoa.h"
#include "hBOAmain.h"

//---------------------------

int numRuns;                // the number of runs of algorithm to perform

int numDiscrete;           // number of discrete variables (not including the continuous variables)
int numContinuous;         // number of continuous variables 
int d;                     // dimension of the problem
long N;                    // size of the population
long offspringSize;        // size of offspring
long M;                    // size of the parent population
double percentM;            // how much percent of N should M be?
double offspringPercentage; // size of offspring in % of N
long tMax;                 // number of generations (if all other term. criteria not satisfied)

DiscretizationParams discretizationParams; // parameters for the discretization

char stopWhenFoundOptimum; // stop if optimum has been found (if possible with fitness..)?
double maxAverageFitness;  // stop when fitness exceeds maxAverageFitness (-1...ignore)
double maxOptimal;         // stop when the rate of optimal and nonoptimal individuals reaches this value
int maxFailures;           // stop when the number of failed runs reaches this value

int fitnessN;	           // the number of used fitness function
double fitnessNoiseVariance;// variance of the additional noise for fitness function
Fitness *fitness;          // fitness function description
double fitnessParams[10];  // fitness parameters
long maxFC;                // maximal number of fitness calls
char *fitnessFilename;     // filename of the file with fitness parameters

int reorderN;              // the number of reordering operator
double reorderingParams[1]; // params for reordering operator

double epsilon;             // when univariate frequencies are closer than epsilon to 0 or 1 terminate

int selectionN;            // number of the used selection method
Selection *selection;      // selection function
int tournamentSize;        // tournament size (for tournament selection)
double boltzmannBeta;      // beta parameter for the Boltzmann selection

int recombinationN;                      // number of the used recombination method
Recombination *recombination;            // recombination method
RecombinationParams recombinationParams; // parameters for recombination
double dependencyTreshold;               // dep. treshold for BMDA
char displayDeps;                        // display dependencies?

int replacementN;        // number of the used replacement method
Replacement *replace;    // replacement method

char *outputFileName;    // name of output file (if any)
char *paramFile;         // name of the file with parameters

FILE *outputGeneral;     // handle for the file with general info (parameters)

long randSeed;		 // random seed

//int  numBusinessMen;
//int  minHamming;

int metricN;             // a metric to use

StatisticsParams statisticsParams; // parameters for statistics

//char useBisection;

//---------------------------

ParamStruct params[] = {
	{PARAM_INT,
		"numRuns",
		&numRuns,
		"1",
		"The number of runs to perform",
		NULL},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_INT,
		"recombination",
		&recombinationN,
		"1",
		"The number of recombination method",
		&getRecombinationDesc},

	{PARAM_INT,
		"blockSize",
		&recombinationParams.blockSize,
		"1",
		"Size of the block for UMDA to combine",
		NULL},

	{PARAM_DOUBLE,
		"dependencyTreshold",
		&recombinationParams.dependencyTreshold,
		"3.84",
		"The treshold for dependencies in BMDA(take only greater)",
		NULL},
	{PARAM_INT,
		"maxIncoming",
		&recombinationParams.maxIncoming,
		"50",
		"Maximal number of incoming edges in dep. graph for XBMDA",
		NULL},

	{PARAM_INT,
		"metric",
		&metricN,
		"2",
		"Metric to use in the BOA",
		&boaGetMetricDescription},

	{PARAM_INT,
		"dBOAMetric",
		&recombinationParams.dBOAMetricN,
		"0",
		"Metric to use in the dBOA",
		&dBOAGetMetricDescription},

	{PARAM_INT,
		"allowAdditions",
		&recombinationParams.allowAdditions,
		"1",
		"Allow edge additions when constructing a network?",
		&yesNoDescriptor},

	{PARAM_INT,
		"allowRemovals",
		&recombinationParams.allowRemovals,
		"1",
		"Allow edge removals when constructing a network?",
		&yesNoDescriptor},

	{PARAM_INT,
		"allowReversals",
		&recombinationParams.allowReversals,
		"0",
		"Allow edge reversals when constructing a network?",
		&yesNoDescriptor},

	{PARAM_INT,
		"allowJoints",
		&recombinationParams.allowJoints,
		"0",
		"Allow variable grouping when constructing a network (in hBOA only)?",
		&yesNoDescriptor},

	{PARAM_INT,
		"useDefaultTables",
		&recombinationParams.useDefaultTables,
		"0",
		"Use default tables in MDL metric?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"priorNetwork",
		&recombinationParams.priorNetwork,
		"0",
		"Use prior network?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"injectGoodGuys",
		&recombinationParams.injectGoodGuys,
		"0",
		"Inject good guys in the initial population (if possible)?",
		&yesNoDescriptor},

	{PARAM_DOUBLE,
		"logKappa",
		&recombinationParams.logKappa,
		"0",
		"Logarithm of the kappa penalty for edges that do not match the prior network.",
		NULL},

	{PARAM_DOUBLE,
		"logBeta",
		&recombinationParams.logKappa,
		"5",
		"Logarithm of the beta reward for edges that match the prior network.",
		NULL},

	{PARAM_CHAR,
		"useBoltzmannFrequencies",
		&recombinationParams.useBoltzmannFrequencies,
		"0",
		"Use boltzmann frequencies?",
		NULL},

	{PARAM_DOUBLE,
		"recombinationBoltzmannBeta",
		&recombinationParams.recombinationBoltzmannBeta,
		"1.5",
		"Beta for the boltzmann frequencies",
		NULL},

	{PARAM_LONG,
		"populationSize",
		&N,
		"100",
		"Size of the population",
		NULL},

	{PARAM_DOUBLE,
		"offspringPercentage",
		&offspringPercentage,
		"100",
		"Size of offspring to create in % from population",
		NULL},

	{PARAM_DOUBLE,
		"parentsPercentage",
		&percentM,
		"100",
		"The number of parents to select (% from population)",
		NULL},

	{PARAM_INT,
		"problemSize",
		&numDiscrete,
		"50",
		"number of discrete variables (not continuous!)",
		NULL},

	{PARAM_INT,
		"numContinuous",
		&numContinuous,
		"0",
		"Number of continuous variables in a problem",
		NULL},

	{PARAM_INT,
		"discretizationType",
		&discretizationParams.discretizationType,
		"0",
		"Type of discretization (none=0, FHHt=1)",
		NULL},

	{PARAM_INT,
		"bitsPerVariable",
		&discretizationParams.bitsPerVariable,
		"0",
		"Number of bits per variable",
		NULL},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},  

	{PARAM_LONG,
		"maxNumberOfGenerations",
		&tMax,
		"50",
		"Maximal Number of Generations to Perform",
		NULL},

	{PARAM_LONG,
		"maxFitnessCalls",
		&maxFC,
		"-1",
		"Maximal Number of Fitness Calls (-1 when unbounded)",
                NULL},

	{PARAM_DOUBLE,
		"epsilon",
		&epsilon,
		"-1",
		"Termination treshold for distance of univ. freq. from (0,1)",
		NULL},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_INT,
		"selection",
		&selectionN,
		"1",
		"Selection method",
		&getSelectionDesc},

	{PARAM_INT,
		"tournamentSize",
		&tournamentSize,
		"2",
		"Tournament size (when tournament selection)",
		NULL},

	{PARAM_DOUBLE,
		"boltzmannBeta",
		&boltzmannBeta,
		"1.5",
		"Beta parameter for the Boltzmann selection (makes sense from 0.5-2.5)",
		NULL},

	{PARAM_INT,
		"numClusters",
		&(recombinationParams.numClusters),
		"1",
		"Number of clusters (if applicable)",
		NULL},

	{PARAM_CHAR,
		"phenotypicClustering",
		&(recombinationParams.phenotypicClustering),
		"0",
		"Clusters points according to the phenotype?",
		NULL},

	{PARAM_INT,
		"numRestarts",
		&(recombinationParams.numRestarts),
		"1",
		"Number of restarts for k-means (if applicable)",
		NULL},

	{PARAM_CHAR,
		"clusterFitnessProportional",
		&(recombinationParams.fitnessProportionalClusterReproduction),
		"1",
		"Size of cluster kids proportional to its avg fitness?",
		&yesNoDescriptor},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_INT,
		"replacement",
		&replacementN,
		"4",
		"Replacement method",
		&getReplacementDesc},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_INT,
		"fitnessFunction",
		&fitnessN,
		"9",
		"Number of fitness function to use",
		&getFitnessDesc},

	{PARAM_DOUBLE,
		"noiseVariance",
		&fitnessNoiseVariance,
		"0",
		"Variance of the additional noise to fitness",
		NULL},

	{PARAM_DOUBLE,
		"lowerBound",
		&(fitnessParams[1]),
		"0",
		"The lower bound for real valued fitness function",
		NULL},

	{PARAM_DOUBLE,
		"upperBound",
		&(fitnessParams[2]),
		"1",
		"The upper bound for real valued fitness function",
		NULL},

	{PARAM_DOUBLE,
		"numOptima",
		&(fitnessParams[4]),
		"21",
		"Number of optima for the real valued function with a variable number of optima.",
		NULL},

	{PARAM_DOUBLE,
		"twomaxFactor",
		&(fitnessParams[3]),
		"1.0",
		"The factor to bias two-max to the right",
		NULL},

	{PARAM_DOUBLE,
		"k",
		&fitnessParams[3],
		"1",
		"k for NK fitness function",
		NULL},

	{PARAM_STRING,
		"fitnessFile",
		&fitnessFilename,
		NULL,
		"name of the file with fitness parameters (if allowed)",
		NULL},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_INT,
		"reordering",
		&reorderN,
		"0",
		"Number of reordering to use",
		&getReorderingDescription}, 

	{PARAM_DOUBLE,
		"disorderK",
		&(reorderingParams[0]),
		"2",
		"k for disorder-k reordering operator",
		NULL},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_STRING,
		"outputFile",
		&outputFileName,
		NULL,
		"Output file name",
		NULL},

	{PARAM_CHAR,
		"outputFitness",
		&statisticsParams.OutputFitness,
		"1",
		"Output fitness info to file?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"outputUMF",
		&statisticsParams.OutputUMF,
		"0",
		"Output univariate marginal frequencies to file?",
		&yesNoDescriptor},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_CHAR,
		"displayFitness",
		&statisticsParams.displayFitness,
		"1",
		"Display fitness information?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"displayOrdering",
		&statisticsParams.displayOrdering,
		"0",
		"Display ordering parameter?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"displayGuidance",
		&statisticsParams.displayGuidance,
		"0",
		"Display guidance of the search?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"displayUMF",
		&statisticsParams.displayUMF,
		"0",
		"Display univariate marginal frequencies",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"displayNumOptimal",
		&statisticsParams.displayNumOptimal,
		"1",
		"Display the number of opt. solutions found so far",
		&yesNoDescriptor},

	{PARAM_DOUBLE,
		"guidanceTreshold",
		&statisticsParams.guidanceTreshold,
		"0.1",
		"How far can UMF from 0.5 be to keep the bit undecided?",
		NULL},

	{PARAM_CHAR,
		"displayDependencies",
		&recombinationParams.displayDependencies,
		"0",
		"Display dependencies in BMDA?",
		&yesNoDescriptor},

	{PARAM_CHAR,
		"pause",
		&statisticsParams.pause,
		"0",
		"Wait for enter after displaying statistics?",
		&yesNoDescriptor},

	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	{PARAM_CHAR,
		"stopWhenFoundOptimum",
		&stopWhenFoundOptimum,
		"1",
		"Stop if the optimum was found?", 
		&yesNoDescriptor},

	{PARAM_DOUBLE,
		"maxAverageFitness",
		&maxAverageFitness,
		"-1",
   "Stop when avg. fitness exceeds this number (-1 for ignoring)",
		NULL},

	{PARAM_DOUBLE,
		"maxOptimal",
		&maxOptimal,
		"-1",
		"Stops when the rate of opt. and nonopt. ind. reaches this value (-1 is ignore)",
		NULL},

	{PARAM_INT,
		"maxFailures",
		&maxFailures,
		"-1",
		"Stops when the number of failures is higher than this value",
		NULL},

	{PARAM_LONG,
		"randSeed",
		&randSeed,
		"time",
		"Random Seed",
		NULL},


	{PARAM_DIVIDER,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL},

	//  {PARAM_CHAR,
	//   "useBisection",
	//   &useBisection,
	//   "0",
	//   "Use bisection to size populations until success", 
	//   &yesNoDescriptor},

	{PARAM_END,
		NULL,
		NULL,
		NULL,
		NULL}
};

//=======================================================================

int help(char what) {
	//  printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");

	printf("------------------------------------------\n");
	printf("Bivariate Marginal Frequencies Algorithm\n");
	printf("Copyright (c) 1998 Martin Pelikan\n");
	printf("------------------------------------------\n\n");

	if (what==0) {
		printf("Command line parameters:\n");
		printf("-h                   display this help screen\n");
		printf("-config=<filename>   configuration file name\n");
		printf("-configDesc          print out the description of configuration file\n");
	}
	else {
		printParamsDescription(stdout,params);
	}

	return 0;
}

//=======================================================================

int generalInfo(FILE *output) {

	//------------------------------------------------------
	// display intro and the values for all the parameters
//comment by pikachu123
	fprintf(output,"----------------------------------------------------\n");
	fprintf(output,"    Hierarchical Bayesian Optimization Algorithm\n");
	fprintf(output,"    (c) 1999 Martin Pelikan\n");
	fprintf(output,"----------------------------------------------------\n");

	if (paramFile)
		printf("Using parameter file %s...\n\n",paramFile);
	else
		printf("Using defaults for all parameters (no parameter file specified)...\n\n");

	//  printParamValues(output,params);
	//  fprintf(output,"\n");

	return 0;
}

//==============================================================

void changedPopulationSize() {
	freePrecomputedCummulativeLogarithms();
	if (recombinationN==BOA_RECOMBINATION) {
		if (metricN==K2_METRIC)
			precomputeCummulativeLogarithms(N+3);

		// for mdl we need to precompute combination numbers

		if (metricN==MDL_METRIC)
			precomputeCummulativeLogarithms(numDiscrete+1);
	}
	if (recombinationN==DECISION_TREE_BOA_RECOMBINATION) {
		precomputeCummulativeLogarithms(N+3);
	}

	// set the size of parent population (having given % of population)

	M = (long) ((double) percentM/100*N);

	// set the size of offspring

	offspringSize = (long) ((double) offspringPercentage/100*N);

	// set random seed

	setSeed(randSeed);
}

//==============================================================

void writeFinalInfo(RunLog *runLog, int numRuns) {
	// write down the general info to the output file (if any)

	if (outputFileName) {
		char s[100];
		printf("\n\n\n\n\n\n\n\n\n\n\n\nRunning with output to %s\n",outputFileName);
		sprintf(s,"%s.generalInfo",outputFileName);
		outputGeneral=fopen(s,"w");
		generalInfo(outputGeneral);
	}
	else
		outputGeneral=NULL;

	finalStatistics(runLog,numRuns,outputGeneral);
	fclose(outputGeneral);
	outputGeneral=NULL;
}

//==============================================================

int startUp(int argc, char **argv,int pop,int ell) {
	// set the priority to the lowest value

	//setpriority(PRIO_PROCESS, getpid(), 19);

	// help requested?

	if (isArg("-h",argc,argv)) {
		help(1);
		exit(0);
	}

	// description of parameters requested?

	if (isArg("-configDesc",argc,argv)) {
		help(1);
		exit(0);
	}

	// read the paramters file name and read it all

	paramFile = getArgString("-config=",argc,argv);
	getParamsFromFile(paramFile,params);
        N=pop;
        numDiscrete=ell;
	// set used recombination method
	// This is hBOA, please do not change it, everything will collapse!!!
	recombinationN = 1; 
	recombination = (Recombination*) getRecombination(recombinationN);

	boaSetMetric(metricN); // well this line is not nice but necessary !!!!!!

	// cut here the stuff for initstuff (re:metrics)

	// prepare, assign, and initialize fitness function

	// check for size limitation 
	switch(fitnessN) {
	case 0:
		// OneMax: 250 bits
		if (numDiscrete > 250) {
			printf("\n");
			printf("The problem size limitation for OneMax is 250 bits.\n");
			printf("\n");
			exit(-1);
		}
		break;
	case 1:
		// deceptive: 80 bits (and % 4 = 0)
		if ((numDiscrete > 80) || ((numDiscrete % 4) != 0)) {
			printf("\n");
			printf("The problem size for deceptive can only be multiples of 4 and less than or equal to 80 bits.\n");
			printf("\n");
			exit(-1);
		}
		break;
	case 2:
		// realValued (sphere): 120 bits (and % 12 = 0)
		if ((numDiscrete > 120) || ((numDiscrete % 12) != 0)) {
			printf("\n");
			printf("The problem size for sphere can only be multiples of 12 and less than or equal to 120 bits.\n");
			printf("\n");
			exit(-1);
		}
		break;
	case 3:
		// hTrap: 9, 27, 81, or 243 bits
		if ((numDiscrete != 9) &&
		  (numDiscrete != 27) &&
		  (numDiscrete != 81) &&
		  (numDiscrete != 243)) {
			printf("\n");
			printf("The problem size for hTrap can only be 9, 27, 81, or 243 bits.\n");
			printf("\n");
			exit(-1);
		}
		break;
	case 4:
		// hDeceptive: 9, 27, 81, or 243 bits
		if ((numDiscrete != 9) &&
		  (numDiscrete != 27) &&
		  (numDiscrete != 81) &&
		  (numDiscrete != 243)) {
			printf("\n");
			printf("The problem size for hDeceptive can only be 9, 27, 81, or 243 bits.\n");
			printf("\n");
			exit(-1);
		}
		break;
	default:
		// User defined function: 100 bits
		if (numDiscrete > 100) {
			printf("\n");
			printf("The problem size limitation for user defined functions is 100 bits.\n");
			printf("\n");
			exit(-1);
		}
		break;
	}
	// Always user's fitness function
	fitnessN = 9;//user_define_fit 
	fitnessParams[0] = numDiscrete;
	fitness = (Fitness*) getFitness(fitnessN);

	if (fitnessFilename){
		if (fitness->load) {
			fitness->load(fitnessFilename,fitnessParams);
			numDiscrete = (int) fitnessParams[0];
			numContinuous = 0; // !!! update for continuous variables !!!
		}
		else {
			fprintf(stderr,"ERROR: Fitness filename set though the fitness has no load function!");
			exit(-1);
		}
	}
	else{
		if (fitness->init)
			fitness->init(fitnessParams);
	}

	// set the additional noise to the fitness

	setFitnessNoiseVariance(fitnessNoiseVariance);

	// init reordering method !!! works only for discrete variables !!! (i don't care now)

	((Reordering*) getReorderingOperator(reorderN))->init(numDiscrete,reorderingParams); 

	// set used selection method
	// Always let users use tournament selection with replacement
	selectionN = 1;
	switch (selectionN) {
	case 0: 
		selection = &selectionProportional; 
		break;

	case 1: 
		selection = &selectionTournamentWithReplacement;
		initSelectionTournament(tournamentSize);
		break;

	case 2: 
		selection = &selectionTruncation;
		break;

	case 3: 
		selection = &selectionTournamentWithoutReplacement;
		initSelectionTournament(tournamentSize);
		break;

	case 4: selection = &selectionNone;
			break;

	case 5: 
			selection = &selectionTournamentWithReplacementAndContinuousSharing;
			initSelectionTournament(tournamentSize);
			break;

	case 6:
			selection = &selectionBoltzmann;
			initSelectionBoltzmann(boltzmannBeta);
			break;

			/*
			   case 7:
			   selection=&selectionTournamentWithoutReplacementExpDec3;
			   initSelectionTournament(tournamentSize);
			   break;
			 */
	default: fprintf(stderr,"ERROR: Couldn't find specified selection method (%u)!",selectionN);
			 exit(-1);
	}

	// set used replacement method
	replacementN = 4;
	switch (replacementN)
	{
	case 0: 
		replace = &replaceWorst; 
		break;

	case 1:
		replace = &replaceAny;
		break;

	case 2:
		replace = &lambdaCommaMju;
		break;

	case 3:
		replace = &lambdaPlusMju;
		break;

	case 4:
		replace = &restrictedTournament;
		break;

	case 5:
		replace = &crowding;
		break;

	default: fprintf(stderr,"ERROR: Couldn't find specified replacement method (%u)!",replacementN);
			 exit(-1);
	}

	// write the general info to stdout

	//fprintf(stdout,"\n");
	generalInfo(stdout);
	//comment by pikachu123 fprintf(stdout,"\n------------------------------------------------------------\n\n");

	// set the prior network (if requested)

	//  if (recombinationParams.priorNetwork)
	//    setPriorNetwork(fitness,&recombinationParams);

	changedPopulationSize();

	// get back

	return 0;
}

int done() {

	// close general output file

	//if (outputGeneral!=NULL)
	//  fclose(outputGeneral);

	// fitness won't be used anymore

	if (fitness->doneFitness)
		fitness->doneFitness();

	// free the memory used by precomputed cummulative logarithms
	if (recombinationN==BOA_RECOMBINATION) { 
		if (metricN==K2_METRIC)
			freePrecomputedCummulativeLogarithms();

		if (metricN==MDL_METRIC)
			freePrecomputedCummulativeLogarithms();
	}

	// free the memory used by the reordering operator (always even when "normal order" is used)

	doneReordering();

	doneParams(params);

	return 0;
}

//=======================================================================

void progressReport(char *s) {
	char filename[200];

	sprintf(filename,"%s.progress",outputFileName);
	FILE *progress = fopen(filename,"a");

	fprintf(progress,s);

	fclose(progress);
}

//=======================================================================
/*
   int bisection() {
   char s[300];
   double accuracy = 0.1;

   int failed;

   long Nmax;
   long Nmin;

   printf("BISECTION: Started.\n");
   progressReport("BISECTION: Started\n");

   printf("BISECTION: Going to determine initial bounds.\n");
   progressReport("BISECTION: Going to determine initial bounds.\n");

   printf("BISECTION:   Simulating N=%lu\n",N);
   sprintf(s,"BISECTION:   Simulating N=%lu\n",N);
   progressReport(s);

   failed = oneTime();

   if (failed) {
   printf("BISECTION:     Result = failed\n");
   progressReport("BISECTION:     Result = failed\n"); do {
   N*=2;
   printf("BISECTION:   Simulating N=%lu\n",N);
   sprintf(s,"BISECTION:   Simulating N=%lu\n",N);
   progressReport(s);
   changedPopulationSize();
   failed=oneTime();
if (failed) {
	printf("BISECTION:     Result = failed\n");
	progressReport("BISECTION:     Result = failed\n");
}
else {
	printf("BISECTION:     Result = succeeded\n");
	progressReport("BISECTION:     Result = succeeded\n");
}
} while (failed);
Nmax=N;
Nmin=Nmax/2;
}
else {
	printf("BISECTION:     Result = succeeded\n");
	progressReport("BISECTION:     Result = succeeded\n"); do {
		N/=2;
		printf("BISECTION:   Simulating N=%lu\n",N);
		sprintf(s,"BISECTION:   Simulating N=%lu\n",N);
		progressReport(s);
		changedPopulationSize();
		failed=oneTime();
		if (failed) {
			printf("BISECTION:     Result = failed\n");
			progressReport("BISECTION: Result = failed\n");
		}
		else {
			printf("BISECTION:     Result = succeeded\n");
			progressReport("BISECTION: Result = succeeded\n");
		}
	} while (!failed);
	Nmax=N*2;
	Nmin=Nmax/2;
}

printf("BISECTION: Initial bounds %lu %lu\n",Nmin,Nmax);
sprintf(s,"BISECTION: Initial bounds %lu %lu\n",Nmin,Nmax);
progressReport(s);
do {
	N = (Nmin+Nmax)/2;
	printf("BISECTION:   Simulating N=%lu\n",N);
	sprintf(s,"BISECTION:   Simulating N=%lu\n",N);
	progressReport(s);
	changedPopulationSize();
	failed=oneTime();
	if (failed) {
		printf("BISECTION:     Result = failed\n");
		progressReport("BISECTION:     Result = failed\n");
	}
	else {
		printf("BISECTION:     Result = succeeded\n");
		progressReport("BISECTION:     Result = succeeded\n");
	}
	if (failed)
		Nmin = N;
	else
		Nmax = N;
} while (((((double)(Nmax-Nmin))/Nmax)>accuracy)&&(Nmax>10));

printf("BISECTION: Finished with N=%lu\n",Nmax);
sprintf(s,"BISECTION: Finished with N=%lu\n",Nmax);
progressReport(s);

return 1;
}
*/

//=======================================================================

int runHeader(int run) {
	printf("\n");
	printf("##################\n");
	printf("###  RUN %3u   ###\n",run);
	printf("##################\n");
	printf("\n");

	return 0;
}

//========================================================================

int oneTime() {
        
	double startTime,endTime;
	Population theLastPopulation;
	int run;
	RunLog runLog;
	int totalFailures;
	long generations;

	// init the timing

	startTime = getTime();

	// init the log array for all runs

	initRunLog(&runLog,numRuns);

	// no failures yet

	totalFailures = 0;

	// perform all the runs

	//  for (run=0; (run<numRuns)&&((maxFailures==-1)||(totalFailures<=maxFailures)); run++)
	run = 0; 
	
	{
		// output the header for a run

		//    runHeader(run);

		// reset fitness calls

		resetFitnessCalls();

		// run evolutionary algorithm
                
		ea(N,offspringSize,M,numDiscrete,numContinuous,&discretizationParams,&generations,tMax,fitness,maxFC,epsilon,stopWhenFoundOptimum,maxAverageFitness,maxOptimal,selection,recombination,&recombinationParams,replace,&statisticsParams,outputFileName,"0",&theLastPopulation);

		// output the statistics for the last run

		runStatistics(&theLastPopulation,fitness,generations,epsilon,maxOptimal,maxFailures,&runLog,run);

		// failed?

		/*
		   if ((fitness->isBest!=NULL)&&(runLog[run].failed)) {
		   totalFailures++;
		   }
		 */
		// free the population used for statistics (allocated by evolutionary algorithm)

		freePopulation(&theLastPopulation);
	}

	// output final statistics

	//  int failed=0;
	/*
	   if (run==numRuns)
	   failed=0;
	   else
	   failed=1;
	 */
	//    if (!useBisection) 
	{
		if (run==numRuns)
			finalStatistics(&runLog,numRuns,outputGeneral);
		else {
			if (outputGeneral) {
				/* WRITE SOMETHING ELSE HERE */
				/*
				   fprintf(outputGeneral,"==================================================\n");
				   fprintf(outputGeneral,"\nThe computation was interrupted since the maximal number of failures was reached!\n");
				   fprintf(outputGeneral,"Runs performed: %u\n",run);
				   fprintf(outputGeneral,"Failed in     : %u runs\n",totalFailures);
				 */
			}
			/*
			   printf("==================================================\n");
			   printf("\nThe computation was interrupted since the maximal number of failures was reached!\n");
			   printf("Runs performed: %u\n",run);
			   printf("Failed in     : %u runs\n",totalFailures);      
			   failed=1;
			 */
		}
	}
	/*
	   else
	   {
	   if (!failed)
	   writeFinalInfo(&runLog,numRuns);

	   printf("BISECTION: Writing final info.\n");
	   }
	 */
	endTime = getTime();
	//printf("   Time spent %f seconds\n",endTime-startTime);

	// run log is not needed anymore
        int fdopt=0;
        if(runLog[run].bestFound)
           fdopt=1;

	doneRunLog(&runLog,numRuns);

	//return 0;
	return fdopt;
	// return failed;
}

//=======================================================================

int hBOAmain(int argc, char **argv,int ell,int pop) {
	// set the priority to the lowest value

	//setpriority(PRIO_PROCESS, getpid(), 19);

	// initialize parameters
        //N=pop;
        //numDiscrete=ell;
	startUp(argc, argv,pop,ell);

	// now go ahead with the bisection or one 

	//  if (useBisection)
	//  bisection();
	//else
	int fdopt=oneTime();

	// things done

	done();

	if (outputGeneral)
		fclose(outputGeneral);

	// get back

	//return 0;
        return fdopt;
}
