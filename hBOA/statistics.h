#ifndef _statistics_h_
#define _statistics_h_

typedef struct {

    char displayFitness;
    char displayOrdering;
    char displayGuidance;
    char displayUMF;
    char displayNumOptimal;

    char OutputFitness;
    char OutputUMF;
    FILE *outputFitness;
    FILE *outputUMF;

    char pause;

    double guidanceTreshold;

} StatisticsParams;

typedef struct {

    int        numDiscrete;
    int        numContinuous;

    Individual best;
    Individual worst;
    double     avgFitness;
    double     ordering;
    char       bestDefined;
    char       bestFound;
    long       bestFoundIn; 
    char       converged;
    char       convergedToTheBest;
    long       timeToConvergence;
    long       generationsToConvergence;
    char       failed;

    float      goodBBs;

    char       maxOptimalReached;
    long       maxOptimalReachedIn;

} RunLogItem;

typedef RunLogItem *RunLog;

int prepareFiles(StatisticsParams *statisticsParams, char *fileRoot, char *fileExtension );
int closeFiles(StatisticsParams *statisticsParams);

int statistics(Population *population, long t, Fitness *fitness, StatisticsParams *statisticsParams);

int runStatistics(Population *population, 
		  Fitness *fitness, 
		  long generations,
		  double epsilon, 
		  double maxOptimal, 
		  int maxFailures,
		  RunLog *runLog, 
		  int run);
int finalStatistics(RunLog *runLog, int numRuns, FILE *outputFile);

int initRunLog(RunLog *runLog, int numRuns);
int doneRunLog(RunLog *runLog, int numRuns);

int selectionStatistics(Population *old, Population *selected);

#endif
