#ifndef _select_h_
#define _select_h_

int selectionProportional(Population *population, Population *parents, long M);
int initSelectionTournament(int tournamentSize_);
int selectionTournamentWithReplacement(Population *population, Population *parents, long M);
int selectionTournamentWithoutReplacement(Population *population, Population *parents, long M);
//int selectionTournamentWithoutReplacementExpDec3(Population *population, Population *parents, long M);
int selectionTresholding(Population *population, Population *parents, long M);
int initSelectionTresholding(int n_,int t);
int selectionTruncation(Population *population, Population *parents, long M);
int selectionTruncation2(Population *population, Population *parents, long M);
int selectionNone(Population *population, Population *parents, long M);
int selectionTournamentWithReplacementAndContinuousSharing(Population *population, Population *parents, long M);
int selectionTournamentWithReplacementAndCoevolutionarySharing(Population *population, Population *parents, long M);
int selectionBoltzmann(Population *population, Population *parents, long M);
int initSelectionBoltzmann(double new_beta);

//int initBusinessMenPopulation(int nBM, int n, int d);

char *getSelectionDesc(int n);
int commonOnes(char *a, char *b, int n);
int divideBest(Population *population, long left, long right, int numDiscrete, int numContinuous, long M);

typedef int Selection(Population *population, Population *parents, long M);

#endif
