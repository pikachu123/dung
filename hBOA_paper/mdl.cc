#include "mdl.h"

#include "counts.h"
#include "population.h"
#include "memalloc.h"
#include "mymath.h"
#include "utils.h"
#include "recombination.h"
#include "index.h"
#include "frequencyTree.h"

//#define DISPLAY_THE_INSTANCES
//#define DEBUG

int mdlComputeAdditionGains( int i, 
			     double oldContribution,
			     double **additionGain, 
			     int *updateIdx, 
			     int numUpdated, 
			     int *parentList,
			     int numParents, 
			     Population *P,
			     RecombinationParams *params)
{
  long  j;
  int   k,l;
  long  **count;
  double **pwi,*tpwi;
  long  numParentConfigurations;
  long  numParentConfigurations2;
  long  numParentConfigurations4;
  double result;
  double structureDescriptionLength;
  double tableAndDataDescriptionLength = 0;
  double entropyN;

#ifdef DEBUG
  printf("in the mdl (additions)!\n");
#endif

  // initialize some variables
  
  numParentConfigurations  = (long) 1<<numParents;
  numParentConfigurations2 = (long) numParentConfigurations<<1;
  numParentConfigurations4 = (long) numParentConfigurations2<<1;

  // allocate memory for counts of all n-ths

  if (numUpdated>0)
    count = (long**) Calloc(numUpdated,sizeof(long*));
  else
    count = NULL;
  
  for (l=0; l<numUpdated; l++)
    count[l] = (long*) Calloc(numParentConfigurations4,sizeof(long));

  // allocate memory for p's and so (need only one of these)

  tpwi  = (double*)  Calloc(numParentConfigurations2,sizeof(double));
  pwi   = (double**) Calloc(numParentConfigurations2,sizeof(double*));

  for (j=0; j<numParentConfigurations2; j++)
    pwi[j] = (double*) Calloc(2,sizeof(double));

  // compute the counts -----------------------------------------------

  computeCountsForList(i,updateIdx,numUpdated,parentList,numParents,P,count);

  // compute the frequencies for the edge and its parents before the addition

  structureDescriptionLength  = log2(P->n)+getPrecomputedCummulativeLog(P->n-numParents,P->n)-getPrecomputedCummulativeLog(1,numParents+1);
  
  // for each element of the nodes to be updated update the gain
  
  for (l=0; l<numUpdated; l++)
    {      
      ///printf("Computing gain for (%u,%u)\n",updateIdx[l],i);
      
      // compute the N's
      
      for (j=0; j<numParentConfigurations2; j++)
	{
	  tpwi[j] = 0;
	  
	  for (k=0; k<2; k++)
	    {
	      pwi[j][k] = (double) (count[l][(j<<1)+k])/P->N;
	      tpwi[j]  += pwi[j][k];
	    };
	  
	  ///	  if (i==1)
	  ///	    printf("parent = %u, tpwi[%u] = %f\n",updateIdx[l],j,tpwi[j]);
	};
      

      // compute the new entropy
      
      if (!params->useDefaultTables)
	{
	  entropyN = 0;
      
	  for (j=0; j<numParentConfigurations2; j++)
	    for (k=0; k<2; k++) 
	      if (pwi[j][k]>0)
		entropyN -= pwi[j][k]*log2(pwi[j][k]/tpwi[j]);
      
	  entropyN *= P->N;

	  tableAndDataDescriptionLength = entropyN+double(numParentConfigurations)*log2(P->N)/2;
	  ///	  printf("table complexity: %lf\n",tableAndDataDescriptionLength);
	  ///	  getchar();
	  ///	  printf("Entropy = %lf\n tableDescription = %lf\n",entropyN,double(numParentConfigurations)*log2(P->N));
	};

      // build a default table ?

      if (params->useDefaultTables)
	{
	  if (numParents>=20)
	    printf("Default table for %u->%u\n",updateIdx[l],i);
	  tableAndDataDescriptionLength = buildDefaultTable(P->n,P->N,pwi,tpwi,numParentConfigurations2,numParents+1);
	};
      
      // compute the resulting gain for the addition of an edge from updateIdx[l] to i
      
      result = -oldContribution-structureDescriptionLength-tableAndDataDescriptionLength;
      
      // update the gain
       
      ///printf("Addition gain for (%u,%u)...%f (%f, %f, %f) with %u parents\n",updateIdx[l],i,result,oldContribution,structureDescriptionLength,tableAndDataDescriptionLength,numParents);
      ///getchar();
      
      additionGain[updateIdx[l]][i]=result;
    };
  
  // free the memory
  
  for (l=0; l<numUpdated; l++)
    Free(count[l]);
  
  if (numUpdated>0)
    Free(count);
  
  for (j=0; j<numParentConfigurations2; j++)
    Free(pwi[j]);
  
  Free(pwi);
  Free(tpwi);
  
  // get back
  
#ifdef DEBUG
  printf("Exiting mdl...\n");
#endif

  return 0;
};

// ====================================================================

int mdlComputeRemovalGains( int i, 
			    double oldContribution,
			    double **removalGain, 
			    int *parentList,
			    int numParents, 
			    Population *P,
			    RecombinationParams *params)
{
  long  j;
  int   k,l;
  long  *count;
  double **poi,*tpoi;
  long  numParentConfigurations;
  long  numParentConfigurations2;
  long  numParentConfigurationsD2;
  double result;
  double structureDescriptionLength;
  double tableAndDataDescriptionLength = 0;
  double entropyN;
  long lExp;
  long lTmp;

  ///  printf("in the mdl (removals)!\n");

  // any parents to remove?
  
  if (numParents==0)
    return 0;

  // initialize some variables
  
  numParentConfigurations   = (long) 1<<numParents;
  numParentConfigurations2  = (long) numParentConfigurations<<1;
  numParentConfigurationsD2 = (long) numParentConfigurations>>1;

  // allocate memory for counts of all n-ths

  count = (long*) Calloc(numParentConfigurations2,sizeof(long));

  // allocate memory for p's and so (need only one of these)

  tpoi  = (double*)  Calloc(numParentConfigurationsD2,sizeof(double));
  poi   = (double**) Calloc(numParentConfigurationsD2,sizeof(double*));

  for (j=0; j<numParentConfigurationsD2; j++)
    poi[j] = (double*) Calloc(2,sizeof(double));

  // add the current node into the list of parents to compute the counts (at the end)

  parentList[numParents]=i;

  // compute the counts -----------------------------------------------

  computeCounts(parentList,numParents+1,P,count);

  // compute the model description length decrease (the same for all) as 
  // a sum of the increase in graph size in bits + size of tables in bits

  structureDescriptionLength  = log2(P->n)+getPrecomputedCummulativeLog(P->n-numParents+2,P->n)-getPrecomputedCummulativeLog(1,numParents-1);

  // for each element of the nodes to be removed update the gain

  for (l=0; l<numParents; l++)
    {      
      ///      printf("Compute removal gain for (%u, %u)\n",parentList[l],i);
      
      lExp = (long) 1<<(numParents-l);

      // compute the po's
      
      for (j=0; j<numParentConfigurationsD2; j++)
	{
	  tpoi[j] = 0;
	  
	  lTmp = (((j>>(numParents-l-1))<<(numParents-l))+(j&(((long)1<<(numParents-l-1))-1)))<<1;

	  ///	  printf("lExp = %lu, lTmp = %lu\n",lExp,lTmp);

	  for (k=0; k<2; k++)
	    {
	      poi[j][k] = (double) (count[lTmp+k]+count[lTmp+lExp+k])/P->N;
	      tpoi[j]  += poi[j][k];
	      ///	      if (i==1)
	      ///		printf("parent = %u, poi[%u][%u] = %f\n",parentList[l],j,k,poi[j][k]);
	    };
	};

      // compute the new entropy

      if (!params->useDefaultTables)
	{
	  entropyN = 0;

	  for (j=0; j<numParentConfigurationsD2; j++)
	    for (k=0; k<2; k++) 
	      if (poi[j][k]>0)
		entropyN -= poi[j][k]*log2(poi[j][k]/tpoi[j]);
	  
	  entropyN *= P->N;

	  tableAndDataDescriptionLength = entropyN+double(0.25)*numParentConfigurations*log2(P->N);
	};

      // build a default table ?

      if (params->useDefaultTables)
	{
	  ///	  printf("Default table for %u->%u\n",updateIdx[l],i);
	  tableAndDataDescriptionLength = buildDefaultTable(P->n,P->N,poi,tpoi,numParentConfigurationsD2,numParents-1);
	};

      // compute the resulting gain for the removal of an edge from updateIdx[l] to i

      result = -oldContribution-structureDescriptionLength-tableAndDataDescriptionLength;
      
      // update the gain

      ///      printf("Removal gain for (%u,%u)...%f (%f + %f + %f)\n",parentList[l],i,result,oldContribution,structureDescriptionLength,tableAndDataDescriptionLength);
      ///      getchar();
     
      removalGain[parentList[l]][i]=result;
    };

  // free the memory

  Free(count);

  for (j=0; j<numParentConfigurationsD2; j++)
    Free(poi[j]);
  
  Free(poi);
  Free(tpoi);

  // get back

  ///  printf("Exiting mdl...\n");

  return 0;
};

// =======================================================================================

double mdlComputeIsolatedNodeContribution(int i, Population *P)
{
  double result;
  int    idx[1];
  long   counts[2];

  // compute the frequency

  idx[0]=i;
  computeCounts(idx,1,P,counts);

  // network complexity is simply log n

  result = log2(P->n);

  // table complexity is log N / 2

  result += log2(P->N)/2;

  // entropy is like this

  result -= counts[0]*log2(double(counts[0])/P->N);
  result -= counts[1]*log2(double(counts[1])/P->N);

  ///  printf("Contribution[%u] = %lf (%lf,%lf)\n",i,result,log2(P->n)+log2(P->N)/2,counts[0]*log2(double(counts[0])/P->N)+counts[1]*log2(double(counts[1])/P->N));

  // get back with the result

  return -result;
  
};

// =====================================================================================

double buildDefaultTable(int n, long N, double **p, double *tp, long numParentConfigurations, int numParents)
{
  int k;
  long j;
  long   *index;
  long   numDefaultLines;
  double defaultProb[2], defaultTProb;
  double newDefaultProb[2], newDefaultTProb;
  double maxEntropyDecrease, change, entropyDecrease;
  long   maxEntropyDecreaseIdx;
  double entropy,tableComplexity;
  double oldTableComplexity;
  long   numRows;

  index = (long*) Calloc(numParentConfigurations,sizeof(long));

  defaultTProb   = 0;

  for (j=0; j<numParentConfigurations; j++)
    {
      defaultTProb += tp[j];
      ///      printf("adding %f\n",tp[j]);
    }
  
  for (k=0; k<2; k++)
    {
      defaultProb[k] = 0;
      
      for (j=0; j<numParentConfigurations; j++)
	defaultProb[k] += p[j][k];
    };

  ///  printf("defaultTProb = %lf\n",defaultTProb);

  entropy = 0;

  for (j=0; j<numParentConfigurations; j++)
    for (k=0; k<2; k++)
      if (p[j][k]>0)
	{
	  if (log2(numParentConfigurations)/log2(2)>=20)
	    printf("Prob[%3lu][%3u]...%f, %f, %f\n",j,k,p[j][k],tp[j],p[j][k]/tp[j]);
	  entropy -= p[j][k]*log2(p[j][k]/tp[j]);
	};

  oldTableComplexity = entropy*N+double(0.5)*numParentConfigurations*log2(N);

  // print out old table complexity (without defaults)

  if (log2(numParentConfigurations)/log2(2)>=20)
    printf("Old table complexity: %f (%f,%f)\n",oldTableComplexity, entropy*N, double(0.5)*numParentConfigurations*log2(N));

  entropy = 0;

  for (k=0; k<2; k++)
    if (defaultProb[k]>0)
      {
	if (log2(numParentConfigurations)/log2(2)>=20)
	  printf("Prob...%f, %f, %f\n",defaultProb[k],defaultTProb,defaultProb[k]/defaultTProb);
	entropy -= defaultProb[k]*log2(defaultProb[k]/defaultTProb);
      }

  tableComplexity = N*entropy + log2(N)/2 + log2(numParentConfigurations);

  ///  printf("Only default lines...complexity: %f (%f, %f, %f)\n",tableComplexity,N*entropy,log2(N)/2,log2(numParentConfigurations));

  numDefaultLines = numParentConfigurations;
  numRows         = 0;

  for (j=0; j<numParentConfigurations; j++)
    index[j]=j;

  if (numParentConfigurations>1)
    do {

      maxEntropyDecrease    = -1;
      maxEntropyDecreaseIdx = -1;

      for (j=0; j<numDefaultLines; j++)
	{
	  entropyDecrease = 0;
	
	  newDefaultTProb = defaultTProb - tp[index[j]];

	  for (k=0; k<2; k++)
	    {
	      newDefaultProb[k] = defaultProb[k]-p[index[j]][k];

	      ///	    printf("Prob...%f, %f\n",p[index[j]][k],tp[index[j]]);
	      if (p[index[j]][k]>0)
		entropyDecrease -= p[index[j]][k]*log2(p[index[j]][k]/tp[index[j]]);
	    
	      ///	    printf("Prob...%f, %f\n",newDefaultProb[k],newDefaultTProb);
	      if (newDefaultProb[k]>0)
		entropyDecrease -= newDefaultProb[k]*log2(newDefaultProb[k]/newDefaultTProb);
	    };

	  ///	printf("entropy decrease = %f\n",entropyDecrease);

	  if ((entropyDecrease<maxEntropyDecrease)||(maxEntropyDecrease<0))
	    {
	      maxEntropyDecrease    = entropyDecrease;
	      maxEntropyDecreaseIdx = j;
	    };
	};

      for (k=0; k<2; k++)
	if (defaultProb[k]>0)
	  maxEntropyDecrease += defaultProb[k]*log2(defaultProb[k]/defaultTProb);

      ///    printf("max entropy decrease = %f\n",maxEntropyDecrease);

      //    printf("Num rows = %lu, Parent configs = %lu\n",numRows,numParentConfigurations);

      change = N*maxEntropyDecrease+log2(N)/2+log2(double(numParentConfigurations-numRows)/(numRows+1));

      ///    printf("Change: %f (%f, %f)\n",change,N*maxEntropyDecrease,log2(N)/2);
      ///    getchar();
    
      if (change<0)
	{
	  defaultTProb -= tp[index[maxEntropyDecreaseIdx]];
	
	  for (k=0; k<2; k++)
	    defaultProb[k] -= p[index[maxEntropyDecreaseIdx]][k];
	
	  if ((numDefaultLines-1)!=maxEntropyDecreaseIdx)
	    swapLong(&(index[numDefaultLines-1]),&(index[maxEntropyDecreaseIdx]));
	
	  tableComplexity += change;

	  //	printf("Extracting %u from a default line (%f, %f, %f)\n",maxEntropyDecreaseIdx,N*maxEntropyDecrease,log2(N)/2,log2(double(numParentConfigurations-numRows)/(numRows+1)));

	  numDefaultLines--;
	  numRows++;
	  ///	getchar();
	};

    } while ((change<0)&&(numDefaultLines>1));

  // compute the final table complexity

  entropy = 0;

  for (j=numDefaultLines; j<numParentConfigurations; j++)
    for (k=0; k<2; k++)
      if (p[index[j]][k]>0)
	{
	  if (log2(numParentConfigurations)/log2(2)>=20)
	    printf("Prob...%f, %f, %f\n",p[index[j]][k],tp[index[j]],p[index[j]][k]/tp[index[j]]);
	  entropy -= p[index[j]][k]*log2(p[index[j]][k]/tp[index[j]]);
	}
  if (log2(numParentConfigurations)/log2(2)>=20)
    printf("Default line\n");

  for (k=0; k<2; k++)
    if (defaultProb[k]>0)
      {
	if (log2(numParentConfigurations)/log2(2)>=20)
	  printf("Prob...%f, %f, %f\n",defaultProb[k],defaultTProb,defaultProb[k]/defaultTProb);
	entropy -= defaultProb[k]*log2(defaultProb[k]/defaultTProb);
      }

  if (log2(numParentConfigurations)/log2(2)>=20)
    {
      printf("New table complexity: %f (%f,%f)\n",tableComplexity, entropy*N, double(0.5)*(numParentConfigurations-numDefaultLines+1)*log2(N));
      ///      getchar();
    }

  if (log2(numParentConfigurations)/log2(2)>=20)
    if (numDefaultLines>0)
      { 
	printf("Created %lu default lines (saved proportion %f)...\n",numDefaultLines,double(numDefaultLines-1)/numParentConfigurations);
	if (numParentConfigurations>2)
	  getchar();
      };

  ///  printf("Exiting build-default-table.\n");

  Free(index);

  return tableComplexity;
};

// =====================================================================================

int mdlComputeGroupAdditionGainsBackup( int i,
					int numVars,
					GroupInformation *groupInformation,
					OperatorGain     *operatorGain,
					AcyclicOrientedGraph *G,
					int *updateIdx, 
					int numUpdated, 
					int *parentList,
					int numParents, 
					Population *P,
					RecombinationParams *params)
{
  long           j;
  int            k,l;

  int            parentIndexSize;
  int           *parentIndex;

  int            allIndexSize;
  int           *allIndex;
  int            allIndexTmpSize;
  int           *allIndexTmp;

  double         result;
  double         structureDescriptionLengthWithout;
  double         structureDescriptionLengthWith;
  double         descriptionLengthWithout;
  double         descriptionLengthWith;
 
  int             numInstances;
  char          **instance;
  double         *frequency;

  FrequencyTree  *treeWithout;
  FrequencyTree  *treeWith;

  double         entropyWithout;
  double         entropyWith;

  int            numGroups;

  numGroups = groupInformation->numGroups;

  // compute the structure description length without the addition (only the fields that do change by the addition)

  structureDescriptionLengthWithout  = 
    log2(numGroups-1)+                                                   // storing the number of parents
    getPrecomputedCummulativeLog(numGroups-numParents+1,numGroups)-      // storing the set of parents
    getPrecomputedCummulativeLog(1,numParents);                

  // initialize the frequency tree

  treeWithout = new FrequencyTree();

  // join the indexes of the parents

  parentIndexSize = 0;

  for (k=0; k<numParents; k++)
    parentIndexSize += groupInformation->groupSize[parentList[k]];

  parentIndex = (int*) Calloc(parentIndexSize+1,sizeof(int));

  joinGroupIndexes(groupInformation,parentList,numParents,parentIndex);

  // compute the instances and the frequencies WITHOUT the new guys

  allIndexTmpSize = groupInformation->groupSize[i]+parentIndexSize;
  allIndexTmp     = (int*) Calloc(allIndexTmpSize, sizeof(int));
  
  joinIndexes(parentIndex,
	      parentIndexSize,
	      groupInformation->groupIndex[i],
	      groupInformation->groupSize[i],
	      allIndexTmp);

  treeWithout->computeIndexedFrequencies(P,allIndexTmp,allIndexTmpSize);

  numInstances = treeWithout->getNumInstances();
    
  instance = (char**) Calloc(numInstances,sizeof(char*));
  for (j=0; j<numInstances; j++)
    instance[j] = (char*) Malloc(allIndexTmpSize);
 
  frequency = (double*) Calloc(numInstances,sizeof(double));
    
  treeWithout->getInstancesAndFrequencies(instance,frequency);

  // compute the entropy without the new guys

  entropyWithout = 0;
  for (j=0; j<numInstances; j++)
    entropyWithout -= frequency[j]*log2(frequency[j]/treeWithout->getFrequency(instance[j],parentIndexSize)); 

  // compute the description length without the addition

  descriptionLengthWithout = 
    structureDescriptionLengthWithout+              // storing the number of parents and the set of the parents
    P->N*entropyWithout+                            // storing the compressed date (entropy*N)
    log2(numInstances)+                             // storing the number of instances we have
    0.5*log2(P->N)*numInstances+                    // storing the table length for the frequencies
    numInstances*allIndexTmpSize;                   // storing the instances this group and its parents have in the population

   printf("\nThe description length without the addition:\n");
   printf("Structure description length: %f\n",structureDescriptionLengthWithout);
   printf("Entropy                     : %f\n",entropyWithout*P->N);
   printf("Frequencies                 : %f\n",0.5*log2(P->N)*numInstances);
   printf("Number of instances         : %u\n",numInstances);
   printf("Instances                   : %f\n",(float)numInstances*allIndexTmpSize);
   printf("Total                       : %f\n",descriptionLengthWithout);

//   printf("\nGroup Contribution: %f\n",mdlComputeGroupContribution(i,d,numVars,groupInformation,G,P,params)-log2(numVars));
   getchar();

#ifdef DISPLAY_THE_INSTANCES

  for (j=0; j<numInstances; j++)
    {
      printf("Instance %2u: ");
      for (l=0; l<allIndexTmpSize; l++)
	printf("%u",instance[j][l]);
      printf("  %5.2f\n",frequency[j]);
    };

#endif

  // free the memory

  for (j=0; j<numInstances; j++)
    Free(instance[j]);
  Free(instance);
  Free(frequency);
  Free(parentIndex);
  Free(allIndexTmp);

  delete treeWithout;

  // for each element of the nodes to be updated update the gain

  for (l=0; l<numUpdated; l++)
    {   
      // allocate the memory for the tree

      treeWith = new FrequencyTree();
   
      // compute the structure description length with the addition
      
      structureDescriptionLengthWith  = 
	log2(numGroups-1)+                                              // storing the number of parents
	getPrecomputedCummulativeLog(numGroups-numParents,numGroups)-        // storing the set of parents
	getPrecomputedCummulativeLog(1,numParents+1);

      // prepare the index for the frequencies, etc.

      allIndexSize = allIndexTmpSize+groupInformation->groupSize[updateIdx[l]];
      allIndex     = (int*) Calloc(allIndexSize,sizeof(int));
      
      joinIndexes(groupInformation->groupIndex[updateIdx[l]],
		  groupInformation->groupSize[updateIdx[l]],
		  allIndexTmp,
		  allIndexTmpSize,
		  allIndex);
	    
      // compute the frequencies of the tree WITH the new guy updateIdx[l]

      treeWith->computeIndexedFrequencies(P,allIndex,allIndexSize);
	    
      // total number of instances?

      numInstances = treeWith->getNumInstances();
      
      // allocate and get the instances

      instance = (char**) Calloc(numInstances,sizeof(char*));
      for (j=0; j<numInstances; j++)
	instance[j] = (char*) Malloc(allIndexSize);

      frequency = (double*) Calloc(numInstances,sizeof(double));

      treeWith->getInstancesAndFrequencies(instance,frequency);

      // compute the new entropy
	    
      printf("Parent index size: %u\n",parentIndexSize);
      printf("Added group size : %u\n",groupInformation->groupSize[updateIdx[l]]);
       	    printf("AllIndexTmpSize  : %u\n",allIndexTmpSize);
       	    printf("index: ");
       	    for (j=0; j<allIndexSize; j++)
       	      printf("%u ",allIndex[j]);
       	    printf("\n");

      entropyWith = 0;
      for (j=0; j<numInstances; j++)
	entropyWith -= frequency[j]*log2(frequency[j]/treeWith->getFrequency(instance[j],
									     parentIndexSize+
									     groupInformation->groupSize[updateIdx[l]]));
      
      // compute the description length with an added edge

      descriptionLengthWith = 
	structureDescriptionLengthWith+               // storing the number and the set extended parents
	P->N*entropyWith+                             // storing the compressed data
	log2(numInstances)+                           // storing the number of instances
	0.5*log2(P->N)*numInstances+                  // storing the frequencies
	numInstances*allIndexSize;                    // storing the instances themselves

//       printf("\nThe description length with the addition (%u,%u):\n",updateIdx[l],i);
//       printf("Structure description length: %f\n",structureDescriptionLengthWith);
//       printf("Entropy                     : %f\n",entropyWith*P->N);
//       printf("Frequencies                 : %f\n",0.5*log2(P->N)*numInstances);
//       printf("Number of instances         : %u\n",numInstances);
//       printf("Instances                   : %u\n",numInstances*allIndexSize);
//       printf("Total                       : %f\n",(float)descriptionLengthWith);
//       G->addEdge(updateIdx[l],i);
//       printf("\nGroup Contribution: %f\n",mdlComputeGroupContribution(i,d,numVars,groupInformation,G,P,params)-log2(numVars));
//       G->removeEdge(updateIdx[l],i);
//       getchar();
     	
#ifdef DISPLAY_THE_INSTANCES

      for (j=0; j<numInstances; j++)
	{
	  printf("Instance %2u: ");
	  for (k=0; k<allIndexSize; k++)
	    printf("%u",instance[j][k]);
	  printf("  %5.2f\n",frequency[j]);
	};

      getchar();
      printf("Key pressed.\n");

#endif

      // compute the resulting change in score and assign it where it belongs

      result = descriptionLengthWithout-descriptionLengthWith;

      // assign the gain to the right spot in the array of addition gains

      operatorGain->addition[updateIdx[l]][i]=result;  
      printf("Gain for addition (%u,%u) = %f\n",updateIdx[l],i,result);
      getchar();
      
      // free the memory

      for (j=0; j<numInstances; j++)
	Free(instance[j]);
      Free (instance);
      Free (frequency);
      Free (allIndex);
	    
      delete treeWith;
    };

  // free the memory

  free(allIndexTmp);

  // get back
  
  return 0;
};

// =====================================================================================

int mdlComputeGroupAdditionGains( int i,
				  int numVars,
				  GroupInformation *groupInformation,
				  OperatorGain     *operatorGain,
				  AcyclicOrientedGraph *G,
				  int *updateIdx, 
				  int numUpdated, 
				  int *parentList,
				  int numParents, 
				  Population *P,
				  RecombinationParams *params)
{
  int            l;

  double         result;
  double         contributionBefore;
  double         contributionAfter;
 
  // compute the structure description length without the addition (only the fields that do change by the addition)

  contributionBefore  = groupInformation->contribution[i];              

  // for each element of the nodes to be updated update the gain

  for (l=0; l<numUpdated; l++)
    {   

      // try to add edge and compute the new contribution of the group
   
      G->addEdge(updateIdx[l],i);

      contributionAfter  = mdlComputeGroupContribution(i,numVars,groupInformation,G,P,params);

      // remove the temporary edge

      G->removeEdge(updateIdx[l],i);

      // compute the resulting change in score and assign it where it belongs

      result = contributionAfter-contributionBefore;

      // assign the gain to the right spot in the array of addition gains

      operatorGain->addition[updateIdx[l]][i]=result;
//       printf("Before: %f\n",contributionBefore);
//       printf("After:  %f\n",contributionAfter);
//       printf("Gain for addition (%u,%u) = %f\n",updateIdx[l],i,result);
//       getchar();
     };

  // get back
  
  return 0;
};

// =====================================================================================

int mdlComputeGroupRemovalGains( int i,
				 int numVars,
				 GroupInformation *groupInformation,
				 OperatorGain     *operatorGain,
				 AcyclicOrientedGraph *G,
				 int *parentList,
				 int numParents, 
				 Population *P,
				 RecombinationParams *params)
{
  int            l;

  double         result;
  double         descriptionLengthBefore;
  double         descriptionLengthAfter;
 
  //  printf("Going to compute the removals into the node %u\n",i);

  // compute the structure description length without the addition (only the fields that do change by the addition)

  descriptionLengthBefore = groupInformation->contribution[i];              

  // for each element of the nodes to be updated update the gain

  for (l=0; l<numParents; l++)
    {   

      // try to remove an edge and compute the new contribution of the group
   
      G->removeEdge(parentList[l],i);
      
      descriptionLengthAfter = mdlComputeGroupContribution(i,numVars,groupInformation,G,P,params);
      
      // remove the temporary edge
      
      G->addEdge(parentList[l],i);

      // compute the resulting change in score and assign it where it belongs (the lower description, the greater result)
      
      result = descriptionLengthAfter-descriptionLengthBefore;
      
      // assign the gain to the right spot in the array of addition gains
      
//       operatorGain->removal[parentList[l]][i]=result;  
//       printf("Gain for removal (%u,%u) = %f\n",parentList[l],i,result);
//       getchar();
    };

  // get back
  
  return 0;
 
};

// =====================================================================================

int mdlComputeGroupJointGains( int i,
			       int numVars,
			       GroupInformation *groupInformation,
			       OperatorGain     *operatorGain,
			       int *updateIdx, 
			       int numUpdated, 
			       int *parentList,
			       int numParents, 
			       Population *P,
			       RecombinationParams *params)
{
  //  AcyclicOrientedGraph *newG;

  printf("ERROR: mdlComputeGroupJointGains not implemented yet!\n");
  exit(-1);

  return 0;
};

// =====================================================================================

double mdlComputeGroupContribution(int i, 
				   int numVars,
				   GroupInformation *groupInformation,
				   AcyclicOrientedGraph *G,
				   Population *P, 
				   RecombinationParams *params)
{
  int            j,k;

  double         descriptionLength;

  int           *parentList;
  int            numParents;
  int            numGroups;

  FrequencyTree *tree;
  double         entropy;

  int           *parentIndex;
  int            parentIndexSize;

  int           *allIndexTmp;
  int            allIndexTmpSize;

  char         **instance;
  double        *frequency;
  int            numInstances;

  double         tableDescription, tableDescription1, tableDescription2;

  numGroups  = groupInformation->numGroups;
  
  // initialize the parent list

  parentList = G->getParentList(i);
  numParents = G->getNumIn(i);

  // size of the memory we need to store the group (the size and which variables are in this group)
  
  descriptionLength = log2(groupInformation->groupSize[i])+
    getPrecomputedCummulativeLog(numVars-groupInformation->groupSize[i]+1,numVars)-
    getPrecomputedCummulativeLog(1,groupInformation->groupSize[i]);
  
  // size of the memory we need to store the set of parents of the group

  descriptionLength += 
    log2(numGroups-1)+                                                   // storing the number of parents
    getPrecomputedCummulativeLog(numGroups-numParents+1,numGroups)-      // storing the set of parents
    getPrecomputedCummulativeLog(1,numParents);                

  // initialize the frequency tree

  tree = new FrequencyTree();
  
  // join the indexes of the parents
  
  parentIndexSize = 0;

  for (k=0; k<numParents; k++)
    parentIndexSize += groupInformation->groupSize[parentList[k]];

  parentIndex = (int*) Calloc(parentIndexSize+1,sizeof(int));

  joinGroupIndexes(groupInformation,parentList,numParents,parentIndex);

  // compute the instances and the frequencies WITHOUT the new guys

  allIndexTmpSize = groupInformation->groupSize[i]+parentIndexSize;
  allIndexTmp     = (int*) Calloc(allIndexTmpSize, sizeof(int));
  
  joinIndexes(parentIndex,
	      parentIndexSize,
	      groupInformation->groupIndex[i],
	      groupInformation->groupSize[i],
	      allIndexTmp);

  tree->computeIndexedFrequencies(P,allIndexTmp,allIndexTmpSize);

  numInstances = tree->getNumInstances();
    
  instance = (char**) Calloc(numInstances,sizeof(char*));
  for (j=0; j<numInstances; j++)
    instance[j] = (char*) Malloc(allIndexTmpSize);
 
  frequency = (double*) Calloc(numInstances,sizeof(double));
    
  tree->getInstancesAndFrequencies(instance,frequency);

  // compute the entropy without the new guys

  entropy = 0;
  for (j=0; j<numInstances; j++)
    entropy -= frequency[j]*log2(frequency[j]/tree->getFrequency(instance[j],parentIndexSize)); 

  // table description

  tableDescription1 = 1+log2(numInstances)+0.5*log2(P->N)*numInstances+numInstances*allIndexTmpSize;
  tableDescription2 = 1+0.5*log2(P->N)*(long(1)<<allIndexTmpSize);
  
  tableDescription=(tableDescription1>tableDescription2)? tableDescription2:tableDescription1;

//   if (tableDescription2>tableDescription1)
//     {
//       printf("Tables: %f %f\n",tableDescription1,tableDescription2);
//       getchar();
//     };

  // compute the description length without the addition

  descriptionLength += 
    P->N*entropy+                                   // storing the compressed date (entropy*N)
    tableDescription;                               // table of frequencies and the instances

  // free the memory

  for (j=0; j<numInstances; j++)
    Free(instance[j]);
  Free(instance);
  Free(frequency);
  Free(allIndexTmp);
  Free(parentIndex);

  // destruct the tree

  delete tree;

  // return the result

  return -descriptionLength;
};
