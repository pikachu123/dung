#ifndef _replace_h_
#define _replace_h_

typedef int Replacement(Population *population, Population *offspring);

int replaceWorst(Population *population, Population *offspring);
int replaceAny(Population *population, Population *offspring);
int lambdaCommaMju(Population *population, Population *offspring);
int lambdaPlusMju(Population *population, Population *offspring);
int restrictedTournament(Population *population, Population *offspring);
int crowding(Population *population, Population *offspring);

char *getReplacementDesc(int n);
int divideWorst(Population *population, long left, long right, int numDiscrete, int numContinuous, long M);

#endif
