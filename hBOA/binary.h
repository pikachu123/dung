#ifndef _binary_h_
#define _binary_h_

#include <stdio.h>

int binaryToInt(char *x, char n);
//int indexedBinaryToInt(char *x, int *index, int n);
int intToBinary(char *x, int number, int numBits);

double binaryToDouble(char *x, int n);

int longToBinary(long a, char *x, char n);

int printBinary(FILE *output, char *x, char n);

int unitation(char *x, char n);

inline int indexedBinaryToInt(char *x, int *index, int n)
{
    int val;
    int *idx;
    int *indexN;
    
    val=0;
    indexN = index+n;
    
    for (idx=index; idx<indexN; idx++)
	{
	    val <<= 1;
	    val +=  x[*idx];
	}
    
    return val;
};


#endif
