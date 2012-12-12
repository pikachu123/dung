#include <iostream>
#include <vector>
using namespace std;

#define Debug

extern double evaluate_2(vector<char>, vector<char>, int );
extern int subevaluate(vector<char>, vector<char>, int, int, int );
extern bool checkflag(vector<bool> );
extern double evaluate_3(vector<char> s1, vector<char> s2, vector<char> s3, int numchannel);


#define seqNum 3
#define maxSize 8

int seqSize[seqNum]={5,7,8};
char seqPattern[seqNum][maxSize]={
	{0,0,0,1,2,-1,-1.-1},
	{1,1,0,1,0,2,0,-1},
	{2,2,2,1,2,0,0,0}
};

vector<char>* sequence;
vector<int>* seqVaryIndex;

void initSeq(){
	sequence=new vector<char>[seqNum];
	seqVaryIndex=new vector<int>[seqNum];
	for(int i=0;i<seqNum;++i){
		for(int j=0;j<seqSize[i];++j){
			sequence[i].push_back(seqPattern[i][j]);
			if(seqPattern[i][j]!=i){
				seqVaryIndex[i].push_back(j);
			}
		}
	}
}

int main(){
	initSeq();
//	for(int i=0;i<seqNum;++i){
//		for(int j=0;j<seqSize[i];++j){
//			cout<<sequence[i][j];
//		}
//		cout<<endl;
//	}
        for(int i=0;i<seqNum;++i)
        {
          for(int j=0;j<seqSize[i];++j)
             cout<<(int)sequence[i][j];
          cout<<endl;
        }
	
	cout<<evaluate_3(sequence[0],sequence[1],sequence[2],3)<<endl;
	return 0;
}
