#ifndef _userBOA_h_
#define _userBOA_h_

typedef double typeFitness(int iLength, char *pChromosome);
typedef char   typeIsBest(int iLength, char *pChromosome);
typedef double typeGoodBBs(int iLength, char *pChromosome);

typedef struct	tagFitnessDefinition {
	typeFitness	*fitness;
	typeIsBest	*isBest;
	typeGoodBBs	*goodBBs;
} FitnessDefinition;

extern FitnessDefinition fitnessDefinition;

// hBOA entry point
int hBOAmain(int argc, char **argv,int,int);

#endif
