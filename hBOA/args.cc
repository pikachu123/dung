#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "args.h"
#include "arrays.h"

int isInt(char *s)
{
   char is = 1;
   int  i=0;

   while ((s[i]!=0) && (is))
   {
      if ((s[i]<'0')||(s[i]>'9')) is=0;
      i++;
   }

   return is;
}

int isFloat(char *s)
{
   char is = 1;

   s = s;

   return is;
}

int isArg(char *s, int argc, char **argv)
{

  for (int i=0; i<argc; i++)  
      if (!strcmp(s,argv[i]))
	 return 1;

  return 0;
}

int getArgInt(char *s, int def, int argc, char **argv)
{
  int n  = def;
  int sL = strlen(s);

  for (int i=0; i<argc; i++)
      if ((!strncmp(s,argv[i],sL)) && (isInt(argv[i]+sL)))
         sscanf(argv[i]+sL,"%u",&n);

  return n;
}

long getArgLong(char *s, int def, int argc, char **argv)
{
  long n  = def;
  int sL = strlen(s);

  for (int i=0; i<argc; i++)
      if ((!strncmp(s,argv[i],sL)) && (isInt(argv[i]+sL)))
         sscanf(argv[i]+sL,"%lu",&n);

  return n;
}


float getArgFloat(char *s, float def, int argc, char **argv)
{
  float n  = def;
  int   sL = strlen(s);

  for (int i=0; i<argc; i++)
      if ((!strncmp(s,argv[i],sL)) && (isFloat(argv[i]+sL)))
         sscanf(argv[i]+sL,"%f",&n);

  return n;
}

char *getArgString(char *s, int argc, char **argv)
{
  int sL = strlen(s);

  for (int i=0; i<argc; i++)
      if (!strncmp(s,argv[i],sL))
         return argv[i]+sL;

  return NULL;
}


int getArgPar_(char *s, int *n, float *value, int argc, char **argv)
{
    int sL = strlen(s);
  char *c,*e;

  c=NULL;

  for (int i=0; i<argc; i++)
      if (!strncmp(s,argv[i],sL))
	  c=argv[i]+sL;   

  e = strchr(c,'=');

  if (e==NULL)
    {
       printf("ERROR: '=' missing in parameter declaration %s!\n",s);
       exit(9);
    }

  *e = 0;

  sscanf(c,"%u",n);
  sscanf(e+1,"%f",value);

  *(c-1)=0;

  return 0;
}

int getArgCoupleArray(char *s, void *a_, int argc, char **argv)
{
  char *c=NULL,*p,*q,*r;
  int sL;
  int x,y;
  CoupleArray *a;
  int i;

  a=(CoupleArray*) a_;

  a->clear();

  sL = strlen(s);

  for (i=0; i<argc; i++)
      if (!strncmp(s,argv[i],sL))
         c=argv[i]+sL;   

  if (c!=NULL)
  {
    c++;

    while (strchr(c,'('))
    {
      p = strchr(c,'(');
      q = strchr(c,',');
      r = strchr(c,')');
  
      *q = 0;
      *r = 0;
      sscanf(p+1,"%u",&x);
      sscanf(q+1,"%u",&y);
      a->add(x,y);
      c=r+1;
    }
  }

  printf("\nCouples: %s",s);
  for (i=0; i<a->getSize(); i++)
  {
    a->get(&x,&y,i);
    printf("(%u,%u) ",x,y);
  }
  printf("\n");
  getchar();

  return 0;
}
