#include <stdio.h>
#include <stdlib.h>

#include "stack.h"
#include "memalloc.h"

//==========================================================

IntStack::IntStack(int max)
{
  maxSize = max;
  marker  = 0;
  size    = 0;

  s = (int*) Calloc(max,sizeof(int));
}

//==========================================================

IntStack::~IntStack()
{
  Free(s);
}

//==========================================================

int IntStack::push(int x)
{
  if (size>=maxSize)
    {
      fprintf(stderr,"ERROR: push method called for a full stack!\n");
      exit(-1);
    }

  return s[size++]=x;
}

//==========================================================

int IntStack::pop()
{
  if (size>0)
    return s[--size];
  else
    {
      fprintf(stderr,"ERROR: pop method called for an empty stack!\n");
      exit(-1);
    }

  return 0;
}

//==========================================================

int IntStack::setMarker()
{
  return marker = size;
}

//==========================================================

int IntStack::goToMarker()
{
  return size = marker;
}

//==========================================================

int IntStack::empty()
{
  return (size==0);
}

//==========================================================

int IntStack::notEmpty()
{
  return size;
}

//==========================================================

int IntStack::full()
{
  return (size==maxSize);
}

//==========================================================

int IntStack::getSize()
{
  return size;
}

//==========================================================
