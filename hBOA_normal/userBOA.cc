#include "userBOA.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <ctime>
#define SHOW_BISECTION 1
using namespace std;
FitnessDefinition fitnessDefinition;
extern double evaluate_2(vector<char>, vector<char>, int );
extern int subevaluate(vector<char>, vector<char>, int, int, int );
extern bool checkflag(vector<bool> );
extern double evaluate_3(vector<char> s1, vector<char> s2, vector<char> s3,vector<char> s4, int numchannel);
//extern double evaluate_3(vector<char> s1, vector<char> s2, vector<char> s3, int numchannel);
#define seqNum 4
#define maxSize 11
#define splen 20
//int seqSize[seqNum]={11,13,14};
int seqSize[seqNum]={7,8,9,11};
//int channelopt[seqNum]={2,4,4};
//char pattern[seqNum]={0,1,2};
char pattern[seqNum]={0,1,2,3};
/*char seqPattern[seqNum][maxSize]={
   {0,0,0,1,2,-1,-1,-1},
   {1,1,0,1,2,2,0,-1},
   {2,2,2,1,2,0,0,1}
};*/
char seqPattern[seqNum][maxSize]={
        {0,0,2,0,3,1,1,-1,-1,-1,-1},
        {1,1,1,0,1,2,0,3,-1,-1,-1},
        {2,2,2,3,2,0,3,1,0,-1,-1},
        {3,3,3,2,0,3,1,1,0,2,2}
};
/*char seqPattern3[seqNum][maxSize]={
        {0,0,0,1,2,0,1,1,2,1,2,-1,-1,-1},
        {1,1,0,1,2,0,2,2,0,1,2,2,0,-1},
        {2,2,2,2,0,1,1,2,1,0,1,0,1,0}
};*/
int max(int *list,int count,int ch)
{
   int max_num=-1;
   int max_index=-1;
   for(int i=0;i<count;i++)
   {
      if(i==ch)
         continue;
      else if(max_num < list[i])
      {
         max_num=list[i];
         max_index=i;
      }
   }
   return max_index;
}
int strtoint(char *x,int base,int len)
{
   int res=0;
   int p=(int )pow(base,len-1);
   while(*x!=2)
   {
      res+=(*x)*p;
      p/=base;
      x++;
   }
   return res;
}
double myfit2(int len,char *x)
{
   int digits=0,num=seqNum-1;
   
   while(num>0)
   {
      digits++;
      num/=2;
   }
   //printf("digits:%d\n",digits);
   //getchar();
   int *seq=(int *)malloc(sizeof(int)*splen);
   char *temp=(char*)malloc(sizeof(char)*(digits+1));
   for(int i=0;i<len/2;i++)
   {
      memcpy(temp,x,digits);
      x+=digits;
      temp[digits]=2;
      seq[i]=strtoint(temp,2,digits);
      //printf("%d",seq[i]);
   }
   vector<char> sequence[seqNum];
   vector<int> seqVaryIndex[seqNum];
   int index=0;     
   for(int i=0;i<seqNum;i++)
   {
      for(int j=0;j<seqSize[i];j++)
      {
         
         if(seqPattern[i][j]!=i)
         {
            seqVaryIndex[i].push_back(j);
            sequence[i].push_back(pattern[seq[index]]);
            index++;
         }
         else
            sequence[i].push_back(seqPattern[i][j]);
      }
   }
   int invalid=0;
   int atleastone[seqNum][seqNum];
   for(int i=0;i<seqNum;i++)
   {
      for (int j = 0; j < seqNum; j++) {
         atleastone[i][j]=0;
      }
   }
   for (int i = 0; i < seqNum; i++) {
      for (int j = 0; j < seqSize[i]; j++) {
         atleastone[i][sequence[i][j]]++;
         if(seqPattern[i][j]!=i)
         {
            if(sequence[i][j]==i||sequence[i][j]>=seqNum)
               invalid++;
         }
      }
   }
   for(int i=0;i<seqNum;i++){
      for (int j = 0; j < seqNum; j++) {
         if(atleastone[i][j]==0)
            invalid++;
      }
   }
   free(seq);
   free(temp);
   seq=NULL;
   temp=NULL;
   if(invalid==0)
   {
      //return -1*evaluate_3(sequence[0],sequence[1],sequence[2],sequence[3],4);
      //for(int i=0;i<seqNum;i++)
      //   printf("%d ",atleastone[i]);
      /*for (int i = 0; i < seqNum; i++) {
         for (int j = 0; j < (int)sequence[i].size(); j++) {
            printf("%d ",(int)sequence[i][j]);
         }
         printf("\n");
      }
      getchar();*/
       //return  -1*evaluate_3(sequence[0],sequence[1],sequence[2],3);
      return  -1*evaluate_3(sequence[0],sequence[1],sequence[2],sequence[3],4);
   }
   else
   {
      /*
      for (int i = 0; i < seqNum; i++) {
         for (int j = 0; j < (int)sequence[i].size(); j++) {
            printf("%d ",(int)sequence[i][j]);
         }
         printf("\n");
      }
      getchar();*/
      return -1*pow(2,7+invalid);
   }


}
/*
double myfit(int len,char* x)
{

   int digits=0,num=seqNum;
   
   while(num>0)
   {
      digits++;
      num/=2;
   }
   int *seq=(int *)malloc(sizeof(int)*splen);
   char *temp=(char*)malloc(sizeof(char)*(digits+1));
   for(int i=0;i<len/2;i++)
   {
      memcpy(temp,x,2);
      x+=2;
      temp[digits]='\0';
      seq[i]=strtoint(temp,2);
   }
   vector<char> sequence[seqNum];
   vector<int> seqVaryIndex[seqNum];
   int index=0;
   for(int i=0;i<seqNum;i++)
   {
      for(int j=0;j<seqSize[i];j++)
      {
         
         if(seqPattern[i][j]!=i)
         {
            seqVaryIndex[i].push_back(j);
            if(seq[index]==i)
               sequence[i].push_back((pattern[seq[index]]+1)%seqNum);
            else if(seq[index]>=seqNum)
            {
               if(seq[index]%seqNum==i)
                 sequence[i].push_back(pattern[((seq[index]%seqNum)+1)%seqNum]);
               else
                 sequence[i].push_back(pattern[seq[index]%seqNum]);
            }
            else
               sequence[i].push_back(pattern[seq[index]]);
            index++;
         }
         else
            sequence[i].push_back(seqPattern[i][j]);
      }
   }
   srand((unsigned int)time(NULL));
   
   for(int i=0;i<seqNum;i++)
   {
      vector<int> invalid_channel;
      bool invalid=false;
      int atleast_one[seqNum]={0};
      for(int j=0;j<seqSize[i];j++)
      {
         atleast_one[sequence[i][j]]++;
      }
      for(int j=0;j<seqNum;j++)
      {
         if(atleast_one[j]==0)
         {
            invalid=true;
            invalid_channel.push_back(j);
         }
      }
      int maxch=0;
      if(invalid)
      {
         while(invalid_channel.size()>0)
         {
            vector<int> insert_index;
            for(int j=0;j<seqNum;j++)
               atleast_one[j]=0;
            for(int j=0;j<seqSize[i];j++)
            {
               atleast_one[sequence[i][j]]++;
            }
            maxch=max(atleast_one,seqNum,i);
            for(int j=0;j<seqSize[i];j++)
               if(sequence[i][j]==maxch)
                  insert_index.push_back(j);
            int inserted=rand()%(atleast_one[maxch]);
            sequence[i][insert_index[inserted]]=invalid_channel.back();
            invalid_channel.pop_back();

         }
      }

   }   
   
   //for(int i=0;i<3;i++)
   //{
   //   for(int j=0;j<channelopt[i];j++)
   //   {
   //      if(seq[index]==i || seq[index]>=seqNum)
   //         sequence[i].push_back(seqNum-1);
   //      else
   //         sequence[i].push_back(pattern[seq[index]]);
   //      
   //   }
   //   index+=channelopt[i];
   //   
   //}
   //
   //for (int i = 0; i < seqNum; i++) {
   //   for (int j = 0; j < (int)sequence[i].size(); j++) {
   //      printf("%d ",(int)sequence[i][j]);
   //   }
   //   printf("\n");
   //}
   //getchar();
   
   free(temp);
   free(seq);
   temp=NULL;
   seq=NULL;
   //delete [] sequence;
   //return 1;
   //printf("%f\n",evaluate_3(sequence[0],sequence[1],sequence[2],sequence[3],4));
   
   
   //for (int i = 0; i < seqNum; i++) {
   //   for (int j = 0; j < (int)sequence[i].size(); j++) {
   //      printf("%d ",(int)sequence[i][j]);
   //   }
   //   printf("\n");
   //}
   //printf("%f\n",evaluate_3(sequence[0],sequence[1],sequence[2],sequence[3],4));

   //getchar();
    
   return -1*evaluate_3(sequence[0],sequence[1],sequence[2],sequence[3],4);

}
*/
char mybest(int len,char *x)
{
   if(myfit2(len,x)>=-23.95800265)
      return 1;
   else
      return 0;
}
/*int main(int argc,char **argv)
{
   fitnessDefinition.fitness=myfit2;
   fitnessDefinition.goodBBs=NULL;
   fitnessDefinition.isBest=mybest;
 
   if (argc != 6)
     {
       printf ("GA ell numConvergence lower upper loop\n");
       return -1;
     }

   int ell = atoi (argv[1]);// problem size
   int numConvergence = atoi (argv[2]);	// consecutive succsess
   int lower = atoi(argv[3]);//max pop
   int upper = atoi(argv[4]);//min pop
   int loop = atoi(argv[5]);//trial times
   int middle_average=0;
   for(int it=0; it<loop;it++)
   {
      int nInitial = (int) (0.5 * ell * log((double)ell) / log(2.71828));  // initial population size

      //int nInitial = 10;


      int j;


      int left, right, middle;


      int populationSize = nInitial/2;
      bool foundOptima;

      if (lower < 0 || upper < 0) {

         if (SHOW_BISECTION) printf("Bisection phase 1\n");

         do {

            populationSize *= 2;

            if (SHOW_BISECTION) printf("[%d]: ", populationSize);

            foundOptima = true;


            for (j=0; j<numConvergence; j++) {

               
      


               if (hBOAmain(argc,argv,ell,populationSize)!=1) {

                  foundOptima = false;

                  if (SHOW_BISECTION) {
                     printf("-");
                     fflush(NULL);
                  }
                  break;
               }

               if (SHOW_BISECTION) {
                  printf("+");
                  fflush(NULL);
               }
            }

      if (SHOW_BISECTION) printf("\n");

      } while (!foundOptima);

      left = populationSize/2;
      right = populationSize;
   }
      else {
         left = lower;
         right = upper;
      }


   middle = (left + right)/2;

   if (SHOW_BISECTION) printf("Bisection phase 2\n");

   while ((right > 1.05 * left) && right > left + 2) {

      middle = (left + right) / 2;

      if (SHOW_BISECTION) printf("[%d] left[%d] right[%d]: ", middle,left,right);

      foundOptima = true;


               for (j=0; j<numConvergence; j++) {

                  if (hBOAmain(argc,argv,ell,middle)!=1) {
                     foundOptima = false;
                     if (SHOW_BISECTION) {
                        printf("-");
                        fflush(NULL);
                     }
                     break;
                  }

                  if (SHOW_BISECTION) {
                     printf("+");
                     fflush(NULL);
                  }
               }

               if (foundOptima)
                  right = middle;
               else
                  left = middle;

               if (SHOW_BISECTION) printf("\n");
            

            

            middle = (left + right) / 2;

            printf("===============\n%d\n", middle);
            middle_average += middle;
   }  
}
           printf("===============\n middle average %d",middle_average/loop);

   return 0;   
}*/

int main(int argc,char **argv)
{
 
   fitnessDefinition.fitness=myfit2;
   fitnessDefinition.goodBBs=NULL;
   fitnessDefinition.isBest=mybest;
   for(int i=0;i<10;i++)
   {
   if(hBOAmain(argc,argv,40,5000)==1)
      printf("found opt\n");
   else 
      printf("not found opt\n");
   } 
}
