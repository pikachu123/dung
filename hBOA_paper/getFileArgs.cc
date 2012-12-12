#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "getFileArgs.h"
#include "memalloc.h"
#include "arrays.h"
#include "utils.h"

char *typeDesc[30] = 
{ "char", 
  "int", 
  "long", 
  "float", 
  "char*",
  "(int,int)*"
};

char *yesString = "Yes";
char *noString = "No";

//================================================================

char *yesNoDescriptor(int i)
{
    if (i)
	return yesString;
    else
	return noString;
}

//================================================================

int getFirstString(FILE *f, char *s, char *restC)
{
    register int i;
    char c;

    do 
	{
	    c = fgetc(f);
	} while ((!feof(f)) && 
		 (((c<'a') || (c>'z')) &&
		  ((c<'A') || (c>'Z')) &&
		  ((c<'0') || (c>'9')) &&
		  (c!='_') &&
		  (c!='-')));
    if (!feof(f))
	{
	    i=0;
	    while ((c!=10) && (!feof(f)) && (c!=' ') && (c!='='))
		{
		    s[i++] = c;
		    c = fgetc(f);
		}
	    s[i]=0;
	}

    if (restC)
	*restC=c;

    return (!feof(f));
}

//================================================================

char firstChar(FILE *f)
{
    char c;
    do
	{
	    c=fgetc(f);
	}
    while ((!feof(f)) && ((c==10) || (c==' ')));

    if ((c!=10) && (c!=' '))
	return c;
    else
	return 0;
}

//================================================================

int setParamValue(FILE *f, char *filename, ParamStruct *param)
{
    char s[100];
    int iTmp;
    float fTmp;
    time_t t;

    // set the value of paremeter from the file, if the file is NULL then the value is chosen
    // to be the default value given as the param's attribute

    switch (param->type) 
	{
	case PARAM_CHAR:
	    if (f)
		getFirstString(f,s,NULL);
	    else
		strcpy(s,param->defValue);

	    sscanf(s,"%i",&iTmp);

	    *((char*)param->where)=(char)iTmp;

	    break;
      
	case PARAM_INT: 
	    if (f)
		getFirstString(f,s,NULL);
	    else
		strcpy(s,param->defValue);

	    sscanf(s,"%i",(int*) param->where);

	    break;
      
	case PARAM_LONG: 
	    if (f)
		getFirstString(f,s,NULL);
	    else
		strcpy(s,param->defValue);

	    if (strcmp(s,"time"))
		sscanf(s,"%li",(long*)param->where);
	    else
		*((long*)param->where) = time(&t);

	    break;
      
	case PARAM_DOUBLE:
	    if (f)
		getFirstString(f,s,NULL);
	    else
		strcpy(s,param->defValue);

	    sscanf(s,"%f",&fTmp);
	    *((double*)param->where) = (double) fTmp;
//	    printf("Got the double. (%f)\n",fTmp);
	    //      *((double*)param->where) = string2double(s);
	    break;

	case PARAM_STRING:
	    if (f)
		{
		    getFirstString(f,(char*) s,NULL);
		    *((char**)param->where) = strdup(s);
		}
	    else
		{
		    if (param->defValue)
			*((char**)param->where) = strdup((char*)param->defValue);
		    else
			*((char**)param->where) = (char*) NULL;
		}
	    break;
      

	case PARAM_COUPLEARRAY:
	    if ((f)&&(!feof(f)))
		{
		    char finito;
		    int a,b;
		    char c;
		    char i;

		    *((CoupleArray**)param->where) = (new CoupleArray());

		    // try to read the beginning of the block with pairs of integers (identified by '{')

		    c=firstChar(f);

		    // found '{'?

		    if (c!='{')
			{
			    fprintf(stderr,"In config file %s in definition of %s is \'{\' missing!\n",filename, param->identifier);
			    printf("FOUND => %c\n",c);
			    exit(-1);
			}	
	
		    // until the end of the block with (int,int)'s found (id. by '}'), read all the pairs

		    finito = 0;

		    while (!finito)
			{
			    c=firstChar(f);
	    
			    // end of the block found?

			    if (c=='}')
				finito=1;
			    else

				// beginning of pair found?

				if (c=='(')
				    {

					// read the pair

					i=0;
					c=fgetc(f);
					while ((c!=',') && (!feof(f)))
					    {
						s[i++]=c;
						c=fgetc(f);
					    }
					s[i]=0;
		  
					// missing ',' after the first integer?

					if (c!=',')
					    {
						fprintf(stderr,"In definition %s in file %s, \',\' expected!\n",s,filename);
						exit(-1);
					    }
					sscanf(s,"%i",&a);
					c=fgetc(f);
					i=0;
					while ((c!=')') && (!feof(f)))
					    {
						s[i++]=c;
						c=fgetc(f);
					    }
					s[i]=0;

					// missing the ')' after the second integer?

					if (c!=')')
					    {
						fprintf(stderr,"In definition %s in file %s, \')\' expected!\n",s,filename);
						exit(-1);
					    }
					sscanf(s,"%i",&b);

					// add the last readed pair to the array of pairs

					(*((CoupleArray**) param->where)) -> add(a,b);
				    }
				else
				    {

					// nothing from what we wanted to find found....

					fprintf(stderr,"Unexpected character %c in the definition of %s in file %s\n",c,param->identifier,filename);
					exit(-1);
				    }
	      
			}
		}
	    else
		param->where = NULL;

	    break;
	}

    // get back 
    return 0;
}

//================================================================

int getParamsFromFile(char *filename, ParamStruct params[])
{
    FILE *f=NULL;
    char s[200];
    register int i;
    int numParams;
    char *defined;
    int which;
    char c;

    numParams=0;
    while (params[numParams].type!=PARAM_END) numParams++;

    // allocate memory for the array saying what was defined and what hasn't been defined yet

    defined = (char*) Malloc(numParams);

    // set this array values to 0 (nothing has been defined yet)

    for (i=0; i<numParams; i++)
	defined[i]=0;

    // if there's configuration file name given, this file has got to exist

    if (filename)
	{
	    f = fopen(filename,"r");

	    if (f==NULL)
		{
		    fprintf(stderr,"ERROR: File %s doesn't exist!\n",filename);
		    exit(-1);
		}
	}

    // read configuration file...

    if (f)
      {
	while (!feof(f))
	    {
		// if it is possible, read the first identifier

		if (getFirstString(f,s,&c))
		    {
			which=-1;
			for (i=0; (i<numParams) && (which==-1); i++)
			    if ((params[i].type!=PARAM_DIVIDER) && (!strcmp(params[i].identifier,s)))
				which=i;

			// does identifier not exist?

			if (which==-1)
			    {
				fprintf(stderr,"ERROR: Parameter %s in file %s is not understood!\n",s,filename);
				exit(-1);
			    }

			// defined twice?

			if (defined[which])
			    {
				fprintf(stderr,"ERROR: Parameter %s in file %s was redefined!\n",s,filename);
				exit(-1);
			    }

			// missing '=' after identifier?

			while ((!feof(f)) && (c!='=')) c = fgetc(f);
			if ((feof(f))||(c!='='))
			    {
				fprintf(stderr,"ERROR: Parameter identifier %s in file %s is not followed by '='!\n",params[which].identifier,filename);
				exit(-1);
			    }

			// read the parameter value 

			setParamValue(f,filename,&(params[which]));

			defined[which]=1;
		    }
	    }

	fclose(f);
      };

    // set default values for the rest of parameters (that has not been define in configuration
    // file)

    for (i=0; i<numParams; i++)
	if (!defined[i])
	    setParamValue(NULL,filename,&(params[i]));

    // free allocated memory

    free(defined);

    // get back

    return 0;
}

//================================================================

int printParamsDescription(FILE *out, ParamStruct params[])
{
    int i=0;

    // print the header (description of information pieces that follow)

    printf("Configuration file description:\n\n");
    printf("%-58s %-26s %-10s %-8s\n","Description","Identifier","Type","Default");
    printf("----------------------------------------------------------------------------------------------------------------\n");

    // print the description of all parameters

    while (params[i].type!=PARAM_END)
	{
	    if (params[i].type==PARAM_DIVIDER)
		printf("\n");
	    else
		fprintf(out,"%-58s %-26s %-10s %-8s\n",params[i].description,params[i].identifier, typeDesc[params[i].type], params[i].defValue);
	    i++;
	}

    // get back

    return 0;
}

//================================================================

int printParamValues(FILE *out, ParamStruct params[])
{

    int i=0;

    // print out the header (the description of information that follows)

    fprintf(out,"Parameter Values:\n\n");
    fprintf(out,"%-58s %-26s %-10s %-8s\n","Description","Identifier","Type","Value");
    fprintf(out,"----------------------------------------------------------------------------------------------------------------\n");

    // print out the descriptions and the values of all parameters

    while (params[i].type!=PARAM_END)
	{
	    if (params[i].type!=PARAM_DIVIDER)
		fprintf(out,"%-58s %-26s %-10s %-8s",params[i].description,params[i].identifier, typeDesc[params[i].type], params[i].defValue);

	    switch (params[i].type)
		{
		case PARAM_CHAR:
		    if (params[i].getValueDescription)
			fprintf(out,"%i (%s)",*((char*)params[i].where),(*((GetValueDescription*)params[i].getValueDescription))(*((char*)params[i].where)));
		    else
			fprintf(out,"%i",*((char*)params[i].where));

		    break;

		case PARAM_INT:
		    if (params[i].getValueDescription)
			fprintf(out,"%i (%s)",*((char*)params[i].where),(*((GetValueDescription*)params[i].getValueDescription))(*((int*)params[i].where)));
		    else
			fprintf(out,"%i",*((int*)params[i].where));
		    break;

		case PARAM_LONG:
		    fprintf(out,"%li",*((long*)params[i].where));
		    break;

		case PARAM_DOUBLE:
		    fprintf(out,"%f",*((double*)params[i].where));
		    break;

		case PARAM_STRING:
		    fprintf(out,"%s",*((char**)params[i].where));
		    break;

		case PARAM_COUPLEARRAY:
		    if (params[i].where)
			{
			    int a,b;

			    fprintf(out,"{\n");
			    for (int j=0; j<(*((CoupleArray**)params[i].where))->getSize(); j++)
				{
				    (*((CoupleArray**)params[i].where))->get(&a,&b,j);
				    fprintf(out,"%-58s %-26s %-10s   (%3i,%3i)\n","","","",a,b);
				}	 
			    fprintf(out,"%-58s %-26s %-10s }","","","");
			}
		    break;

		case PARAM_DIVIDER:
		    fprintf(out,"\n");
		    break;
		}

	    fprintf(out,"\n");

	    i++;
	}

    // get back
  
    return 0;
}


//================================================================

int doneParams(ParamStruct params[])
{
  int i;

  i=0;

  while (params[i].type!=PARAM_END)
    {
      if ((params[i].type==PARAM_STRING)&&((*((char**)params[i].where))!=NULL))
	  Free((void*) *((char**)params[i].where));

      i++;
    };

    // get back
  
    return 0;
}
