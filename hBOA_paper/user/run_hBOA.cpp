#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "userBOA.h"

using namespace std;
// User-defined fitness functions
// See the type definition in userBOA.h

// Sample function 1: OneMax
double myOneMax(int iLength, char *pChromosome);
char myOneMaxIsBest(int iLength, char *pChromosome);
double myOneMaxGoodBBs(int iLength, char *pChromosome);

// Sample function 2: order-4 trap
double myTrap4(int iLength, char *pChromosome);

// Sample function 3: realValued (sphere)
double mySphere(int iLength, char *pChromosome);

// Sample function 4: hTrap
double myhTrap(int iLength, char *pChromosome);

// Sample function 5: hDeceptive
double myhDeceptive(int iLength, char *pChromosome);

double dungFunc(int iLength, char *pChromosome);

int main(int argc, char **argv) {
	//
	// Construct the fitness definition for BOA
	//

	// Using sample function 1
	fitnessDefinition.fitness	= &myOneMax;
	fitnessDefinition.isBest	= NULL;
	fitnessDefinition.goodBBs	= NULL;
	// 
	// Only the fitness function is necessary.
	// These two functions are for research purpose.
	//
	//fitnessDefinition.isBest	= &myOneMaxIsBest;
	//fitnessDefinition.goodBBs	= &myOneMaxGoodBBs;

	// Using sample function 2
	fitnessDefinition.fitness	= &myTrap4;
	fitnessDefinition.isBest	= NULL;
	fitnessDefinition.goodBBs	= NULL;

	// Using sample function 3
	fitnessDefinition.fitness	= &mySphere;
	fitnessDefinition.isBest	= NULL;
	fitnessDefinition.goodBBs	= NULL;

	// Using sample function 4
	fitnessDefinition.fitness	= &myhTrap;
	fitnessDefinition.isBest	= NULL;
	fitnessDefinition.goodBBs	= NULL;

	// Using sample function 5
	fitnessDefinition.fitness	= &myhDeceptive;
	fitnessDefinition.isBest	= NULL;
	fitnessDefinition.goodBBs	= NULL;

	//Using dung function
	fitnessDefinition.fitness	= &dungFunc;
	fitnessDefinition.isBest	= NULL;
	fitnessDefinition.goodBBs	= NULL;

	printf("\nThis is a sample user code for using the hBOA library.\n");

	// Run hBOA
	return hBOAmain(argc, argv);
}

double dungFunc(int iLength, char *pChromosome){
	for(int i=0;i<iLength;++i){
		cout<<(int)pChromosome[i]<<" ";
	}
	cout<<endl;
	return 0;

}

// Sample function 1: OneMax
// Only the fitness function is necessary.
// The other two functions are for research purpose.
double myOneMax(int iLength, char *pChromosome) {
	int s = 0;

	for (int i = 0 ; i < iLength ; i ++)
		s += pChromosome[i];

	return ((double) s);
}

char myOneMaxIsBest(int iLength, char *pChromosome) {
	int i;

	for (i = 0 ; (i < iLength) && (pChromosome[i] == 1) ; i ++);

	return (i == iLength);
}

double myOneMaxGoodBBs(int iLength, char *pChromosome) {
	int k = 0;

	for (int i = 0 ; i < iLength ; i ++)
		k += pChromosome[i];

	return (((double )k)/((double )iLength));
}

// Sample function 2: order-4 trap
int trap(char *x, int n, int order);

double myTrap4(int iLength, char *pChromosome) {
	return (double) trap(pChromosome, iLength, 4);
}

int trap(char *x, int n, int order) {
	int order1,s;
	register int i,j,k;

	order1 = order-1;

	s = 0;
	i = 0;

	do {
		k=0;
		for (j=0; j<order; j++)
			k+=x[i++];

		if (k==order)
			s+=order;
		else
			s+=order1-k;
	} while (i<n);

	return s;
}

// Sample function 3: realValued (sphere)
// iLength bits -> [dLower..dUpper]
double decodeRealValue(int iLength, char *pChromosome, double dLower, double dUpper);

double mySphere(int iLength, char *pChromosome) { int iCount;
	int iNumVar;
	int iBitsPerVar;
	double dValue, dTemp;

	iBitsPerVar = 12;
	iNumVar = iLength/iBitsPerVar;

	dValue = 0.0;

	for (iCount = 0 ; iCount < iNumVar ; iCount ++) {
		dTemp = decodeRealValue(iBitsPerVar, &pChromosome[iBitsPerVar*iCount], -1.0, 1.0);
		dValue += (dTemp*dTemp);
	}

	// minimization
	return (-1.0)*dValue;
}

double decodeRealValue(int iLength, char *pChromosome, double dLower, double dUpper) {
	int i;
	double k,p;

	k=0;
	p=1;

	for (i = iLength-1; i>=0; i--) {
		if (pChromosome[i])
			k += p;
		p *= 2.0;
	}

	return (double) dLower + ((dUpper-dLower)/p)*k;
}

// Sample function 4: hTrap
double myGeneralTrap(char *x, int k, double leftPeak, double rightPeak);

double myhTrap(int iLength, char *pChromosome) {
	int    i;
	char   *transform;
	char   nullValue = 3;
	int    level;
	int    numLevels;
	int    last;
	double result;
	int    bonus;
	int    tmp;

	// set the variables
	numLevels = (int) ((double)(log(iLength)/log(3))+0.5); 
	last      = iLength;
	bonus     = 1;
	result    = 0;

	// allocate some memory
	transform = (char *) malloc(iLength);

	// copy the string to the transform array
	for (i=0; i<iLength; i++)
		transform[i]=pChromosome[i];

	// process all levels (add reinforcements and interpret the string to a higher level)
	for (level=0; level<numLevels; level++) {
		for (i=0; i<last; i+=3) {
			if (level<numLevels-1)
				result += bonus*myGeneralTrap(transform+i,3,1,1);
			else
				result += bonus*myGeneralTrap(transform+i,3,0.9,1);

			tmp = transform[i]+transform[i+1]+transform[i+2];

			if (tmp==0)
				transform[i/3]=0;
			else
				if (tmp==3)
					transform[i/3]=1;
				else
					transform[i/3]=nullValue;
		}

		// increase the bonus on a higher level
		bonus*=3;

		// we've got only half interpretations at this point
		last/=3;
	}

	// free the memory
	free(transform);

	// get back
	return result;
}

double myGeneralTrap(char *x, int k, double leftPeak, double rightPeak) {
	int s;

	s=0;
	for (int i=0; i<k; i++)
		if (x[i]==1)
			s++;
		else
			if(x[i]!=0)
				return 0;

	if (s==k)
		return rightPeak;
	else
		return leftPeak*(k-1-s)/(k-1);
}

// Sample function 5: hDeceptive
double myhDeceptive(int iLength, char *pChromosome) {
	int    i;
	char   *transform;
	char   nullValue = 3;
	int    level;
	int    numLevels;
	int    last;
	double result;
	int    bonus;
	int    tmp;

	// set the variables

	numLevels = (int) ((double)(log(iLength)/log(3))+0.5); 
	last      = iLength;
	bonus     = 1;
	result    = 0;

	// allocate some memory

	transform = (char*) malloc(iLength);

	// copy the string to the transform array

	for (i=0; i<iLength; i++)
		transform[i]=pChromosome[i];

	// process all levels (add reinforcements and interpret the string to a higher level)
	for (level=0; level<numLevels; level++) {

		for (i=0; i<last; i+=3) {
			if (level<numLevels-1)
				result += bonus*myGeneralTrap(transform+i,3,1+0.1/numLevels,1);
			else
				result += bonus*myGeneralTrap(transform+i,3,0.9,1);

			tmp = transform[i]+transform[i+1]+transform[i+2];

			if (tmp==0)
				transform[i/3]=0;
			else
				if (tmp==3)
					transform[i/3]=1;
				else
					transform[i/3]=nullValue;
		}

		// increase the bonus on a higher level

		bonus*=3;

		// we've got only one third of the interpretations at this point

		last/=3;
	}

	// free the memory
	free(transform);

	// get back
	return result;
}
