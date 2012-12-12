#include "index.h"

#include "utils.h"

// =======================================================================================

int joinGroupIndexes(int **groupIndex, int *groupSize, int *index, int n, int *destination)
{
  int i,j,k;

  for (i=0,k=0; i<n; i++)
    for (j=0; j<groupSize[index[i]]; j++)
	destination[k++] = groupIndex[index[i]][k];
   
  return 0;
};


// =======================================================================================

int joinIndexes(int *index1, int n1, int *index2, int n2, int *destination)
{
  int i,j;

  for (i=0; i<n1; i++)
    destination[i] = index1[i];

  for (i=0,j=n1; i<n2; i++,j++)
    destination[j] = index2[i];

  return 0;
};

// =======================================================================================

int sortIndex(int *index, int n)
{
  register int k,l;
  int      minValue,minIndex;

  for (k=0; k<n; k++)
    {
      minIndex = k;
      minValue = index[minIndex];

      for (l=k+1; l<n; l++)
	if (index[l]<minValue)
	  {
	    minIndex = l;
	    minValue = index[minIndex];
	  };

      if (minIndex!=k)
	swapInt(&(index[minIndex]),&(index[k]));
    };

  return 0;
};
