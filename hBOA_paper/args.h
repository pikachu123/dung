#ifndef _args_h_
#define _args_h_

int isArg(char *s, int argc, char **argv);
int getArgInt(char *s, int def, int argc, char **argv);
long getArgLong(char *s, int def, int argc, char **argv);
float getArgFloat(char *s, float def, int argc, char **argv);
char *getArgString(char *s, int argc, char **argv);
int getArgPar_(char *s, int *n, float *value, int argc, char **argv);
int getArgCoupleArray(char *s, void *a_, int argc, char **argv);

#endif
