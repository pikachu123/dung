#include <stdio.h>

#include "binary.h"

int binaryToInt(char *x, char n)
{
  register int i,val;
  
  val=0;
  
  for (i=0; i<n; i++)
    {
      val <<= 1;
      val += x[i];
/*       printf("%u",x[i]); */
    }
/*   printf(" = %u\n",val); */
/*   getchar(); */

  return val;
}

double binaryToDouble(char *x, int n)
{
  register int i;
  double   val;
  
  val=0;

  for (i=0; i<n; i++)
    {
      val *= 2;
      val += x[i];
    }

  return val;
}

//  int indexedBinaryToInt(char *x, int *index, int n)
//  {
//    int val;
//    int *idx;
//    int *indexN;

//    val=0;
//    indexN = index+n;

//    for (idx=index; idx<indexN; idx++)
//      {
//        val <<= 1;
//        val +=  x[*idx];
//      }
 
//    return val;
//  };

int longToBinary(long a, char *x, char n)
{
  for (int i=n-1; i>=0; i--)
    {
      x[i] = a%2;
      a>>=1;
    };

  // get back

  return 0;
}

int printBinary(FILE *output, char *x, char n)
{
  // print out the binary string

  for (int i=0; i<n; i++)
    fprintf(output,"%u",x[i]);
  
  // get back

  return 0;
}

int unitation(char *x, char n)
{
  int u;

  // sum the bits in the string

  u=0;
  for (int i=0; i<n; i++)
    u+=x[i];

  // return the value 

  return u;
}

int intToBinary(char *x, int number, int numBits)
{
  int tmp = number;
  int i=numBits-1;

  while (tmp>0)
    {
      x[i]=tmp%2;
      tmp >>= 1;
      i--;
    };

  for (; i>=0; i--)
    x[i]=0;

  return 0;
};
