#ifndef _reordering_h_
#define _reordering_h_

typedef int InitReorder(int n, double *params);

typedef struct {
  char            *description;
  InitReorder     *init;
} Reordering;

void *getReorderingOperator(int n);
char *getReorderingDescription(int n);

int normalOrdering(int n, double *params);
int disorderK(int n, double *params);

int reorder(char *y, char *x, int n);

int doneReordering();

#endif
