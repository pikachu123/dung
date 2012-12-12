#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define CGITMAX 500
#define CGEPS 1.0E-10
#define MNBGOLD 1.618034
#define MNBGLIMIT 100.0
#define MNBTINY 1.0E-20
#define MNBSHIFT(a,b,c,d) (a)=(b); (b)=(c); (c)=(d);
#define BRITMAX 200
#define BRCGOLD 0.3819660
#define BRZEPS 1.0E-10
#define BRTOL 2.0E-4

#define SI_Rcutoff 1.8
#define SI_p 4
#define SI_q 0
#define SW_SI_lambda 21.0
#define G_SI_lambda 25.0
#define SI_gamma 1.20
#define SI_A 7.049556227
#define SI_B 0.6022245584
#define SI_c0 -0.5
#define SI_c1 0.45

void conjugateGradient(double*, double*, int, int, double, int*);
int lineMinimization(double*, double*, int, int, double*);
double f1dim(double, int, int, double*, double*, double*, int*); 
void mnbrak(double*, double*, double*, double*, double*, double*, 
	    int, int, double*, double*, double*, int*);
void brent(double, double, double, double*, double*, int, int, 
	   double*, double*, double*, int*);
  		     
double StillWeberPotnl(int, double*, double*, double*, double*);
double GongPotnl(int, double*);
void DGongPotnl(int, double*, double*);
