#include <iostream>
#include <vector>

using namespace std;
double evaluate_2(vector<char>, vector<char>, int );
int subevaluate(vector<char>, vector<char>, int, int, int );
bool checkflag(vector<bool> );


double evaluate_3(vector<char> s1, vector<char> s2, vector<char> s3,vector<char> s4, int numchannel){
	double a12 = evaluate_2(s1,s2,numchannel);//cout<<"a12 = "<<a12<<endl;
	double a13 = evaluate_2(s1,s3,numchannel);//cout<<"a13 = "<<a13<<endl;
	double a14 = evaluate_2(s1,s4,numchannel);//cout<<"a23 = "<<a23<<endl;
	double a23 = evaluate_2(s2,s3,numchannel);//cout<<"a12 = "<<a12<<endl;
	double a24 = evaluate_2(s2,s4,numchannel);//cout<<"a13 = "<<a13<<endl;
	double a34 = evaluate_2(s3,s4,numchannel);//cout<<"a23 = "<<a23<<endl;
	return (a12+a13+a14+a23+a24+a34)/6;
}

double evaluate_2(vector<char> s1, vector<char> s2, int numchannel){ //numchannel is # of channel; value = 2 means that channel = 0, 1;
	int size_1 = s1.size();
	int size_2 = s2.size();
	int sum = 0;
	double result=0;
	for(int i=0; i<size_1; i++){
		for(int j=0; j<size_2; j++){
			sum = sum + subevaluate(s1,s2,numchannel,i,j);
			}
		}
	//cout<<"sum = "<<sum<<endl;
	result=double(sum)/(size_1*size_2);
	return result;
}
//Given initial position of two sequence, than this ft returns the time that all channel has already meets at.
//start1 => initial position of sequence 1;
int subevaluate(vector<char> s1, vector<char> s2, int numchannel, int start1, int start2){
	vector<bool> flag(numchannel);
	int step = 0;	
	while(checkflag(flag) == 0){
		step++;
		int d1 = s1[start1];
		int d2 = s2[start2];
		start1++;
		start2++;
		if(start1 >= s1.size() )start1 = start1 - s1.size();
		if(start2 >= s2.size() )start2 = start2 - s2.size();
		if(d1 == d2){
			flag[d1] = 1; 
			}
		}
	//cout<<"step = "<<step<<endl;
	return step;
}

bool checkflag(vector<bool> flag){//if elements of flag don't have any 0, output = 1;
	int size = flag.size();
	for(int i=0; i<size; i++){
		if(flag[i] == 0)return 0;
		}
	return 1;
}
