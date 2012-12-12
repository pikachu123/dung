#ifndef _getFileArgs_h_
#define _getFileArgs_h_

#define PARAM_CHAR          0
#define PARAM_INT           1
#define PARAM_LONG          2
#define PARAM_DOUBLE        3
#define PARAM_STRING        4
#define PARAM_COUPLEARRAY   5
#define PARAM_DIVIDER     100
#define PARAM_END         101

typedef char *GetValueDescription(int n);

typedef struct {
  char type;
  char *identifier;
  void *where;
  char *defValue;
  char *description;
  GetValueDescription *getValueDescription;
} ParamStruct;

int getFirstString(FILE *f, char *s);
int getParamsFromFile(char *filename, ParamStruct params[]);
int printParamsDescription(FILE *out, ParamStruct params[]);
int printParamValues(FILE *out, ParamStruct params[]);
int doneParams(ParamStruct params[]);

char *yesNoDescriptor(int i);

#endif
