#include "discretization.h"

int ea( long N,
	long offspringSize,
	long M,
	int n, 
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

	Population *theLastPopulation);
