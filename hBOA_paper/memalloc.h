#ifndef _memalloc_h_
#define _memalloc_h_

#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>

// inline functions

inline void *Malloc(int x)
{
   void *p;

   if (x==0)
     {
       printf("ERROR: Cannot allocate empty block (malloc) !\n");
       exit(-1);
     };

   p=malloc(x);

   if (p==NULL)
   {
      printf("ERROR: Not enough memory (for a block of size %u)\n",x);
      exit(-1);
   }

   return p;
}

inline void *Calloc(int x, int s)
{
   void *p;

   if (x==0)
     {
       printf("ERROR: Cannot allocate empty block (calloc) !\n");
       exit(-1);
     };

   p=calloc(x,s);

   if (p==NULL)
   {
      printf("ERROR: Not enough memory. (for a block of size %u)\n",x*s);
      exit(-1);
   }

///   cout << "Allocated..." << flush

   return p;
}

inline void Free(void *x)
{
  free(x);
}

#endif
