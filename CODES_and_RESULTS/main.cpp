#include <iostream>
#include <vector>
#include "evaluate.cpp"
using namespace std;

#define Debug

int main(){
//	int sequence0[] = {0,0,0,1,2};
    int sequence0[] = {0,0,2,0,3,1,1};
	vector<int> s0(sequence0, sequence0 + sizeof(sequence0)/sizeof(int));
	//int sequence1[] = {1,1,0,1,2,2,0};
	int sequence1[] = {1,1,1,0,1,2,0,3};
	vector<int> s1(sequence1, sequence1 + sizeof(sequence1)/sizeof(int));
	//int sequence2[] = {2,2,2,1,2,0,0,1};
	int sequence2[] = {2,2,2,3,2,0,3,1,0};
	vector<int> s2(sequence2, sequence2 + sizeof(sequence2)/sizeof(int));
	int sequence3[] = {3,3,3,2,0,3,1,1,0,2,2};
	vector<int> s3(sequence3, sequence3 + sizeof(sequence3)/sizeof(int));
#if 0
	for(int i=0; i<s0.size(); i++)
		cout<<s0[i]<<" ";
	cout<<endl;
	for(int i=0; i<s1.size(); i++)
		cout<<s1[i]<<" ";
	cout<<endl;
	for(int i=0; i<s2.size(); i++)
		cout<<s2[i]<<" ";
	cout<<endl;
#endif
	cout<<evaluate_3(s0,s1,s2,s3,4)<<endl;
	system("pause");
	return 0;
}
