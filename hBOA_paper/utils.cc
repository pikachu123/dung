#include <sys/time.h>
#include <stdlib.h>

//==============================================================

int swapInt(int *a, int *b)
{
  int aux;
  
  aux = *a;
  *a  = *b;
  *b  = aux;

  return 0;
};

//==============================================================

double swapDoubles(double *a, double *b)
{
  double aux;
  
  aux = *a;
  *a  = *b;
  *b  = aux;

  return 0;
};

//==============================================================

int swapLong(long *a, long *b)
{
  long aux;
  
  aux = *a;
  *a  = *b;
  *b  = aux;

  return 0;
};

//==============================================================

void swapPointers(void **a, void **b)
{
  void *aux;

  aux = *a;
  *a  = *b;
  *b  = aux;
}


//==============================================================

// double string2double(char *string)
// {
//   unsigned int i;
//   double get;
//   double ret,j;
  
//   j=0.;
//   ret=0.;

//   for(i=0; i<strlen(string); i++)
//       {
// 	  if(string[i]=='.')
// 	      j=10.;
// 	  else
// 	      {
// 		  get=string[i]-'0';      
// 		  if(j==0.)
// 		      {
// 			  ret=ret*10.;
// 			  ret=ret+get;
// 		      }
// 		  else
// 		      {
// 			  ret=ret+(get/j);
// 			  j=j*10.;
// 		      }
// 	      }
//       }

//   return ret;
// }

//============================================================

double getTime()
{
  timeval tv;

  gettimeofday(&tv,NULL);

  return (tv.tv_sec*1000000+tv.tv_usec)/1000000.0;
}
