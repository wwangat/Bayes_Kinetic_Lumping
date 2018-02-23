#include<iostream>
#include <math.h>
using namespace std;
//part 1: max of an array
double max(double *array, int size){
	int i;
	double max=array[0];
	for (i=1;i<size;i++){
		if (array[i]>max)
			{max = array[i];}
	}
	return max;
}

double stirling(double n){
    //calculate log(n!) approximation
    double result;
    result = 0.5*log(2*3.1415926535)+0.5*log(n)+n*log(n)-n+1.0/(12*n);
//    result = n*log(n)-n;
    return result;
}

//part 2: generate random number
int randi(int min, int max){
	return rand() % (max-min+1)+min;
} //generate random integer between min and max
double randd(){
	double randnum;
	randnum = rand()/(RAND_MAX+1.0);
	return randnum;
} //generate pseudo random number between 0 and 1


double posterior_v2(double **tCount, int **nCount, double alpha_nd,double beta_nd, double alpha_dd, double beta_dd, int sizeMacro){
//calculate gamma conjugate prior posterior probability
	double sum = 0.0, temp;
	int i, j;
	for (i = 0; i<sizeMacro; i++){
        sum += lgamma(tCount[i][i]+alpha_dd)-(tCount[i][i]+alpha_dd)*log(nCount[i][i]+beta_dd);
		for (j = i+1; j<sizeMacro; j++){
			//in this situation, since tCount starts from 0, we should add them by one
			sum += lgamma(tCount[i][j]+alpha_nd)-(tCount[i][j]+alpha_nd)*log(nCount[i][j]+beta_nd);
			sum += lgamma(tCount[j][i]+alpha_nd)-(tCount[j][i]+alpha_nd)*log(nCount[j][i]+beta_nd);
		}
	}
	return sum;
}


/*metastabilityz*/
double metastability(double **matrix, int sizeMacro)
{
    int i,j;
    double temp, temp_meta=0.0;
    for (i=0;i<sizeMacro;i++){
        temp=0.0;
        for (j=0;j<sizeMacro;j++){
            temp=temp+matrix[j][i];   //colsum, bug fix on Mar 15, 2017, should be column summation
        }
        temp_meta = temp_meta+matrix[i][i]/temp;
    }
    return temp_meta;
}

char* num2str(int num)
{
	static char str[10];
	int rem;
	int n, len=0, i;
	n = num;
	while (n!=0)
	{
		len++;
		n /= 10;
	}
	for (i=0;i<len;i++)
	{
		rem = num %10;
		num = num/10;
		str[len -(i+1)] = rem + '0';
	}
	str[len] = '\0';
	return str;
}
/*monte carlo*/
int MonteCarlo(double *array, int size){
	int index = 0, i;
	double randd();
	double randnum, sum = array[0];
	for (i = 1;i<size;i++){
		sum += array[i];
	}
	randnum = sum*randd();
	if (randnum<=array[0])
		index = 0;
	else
	{
		sum = 0.0;
		for (i = 0; i<size-1; i++){
			if(sum+array[i]<randnum && sum+array[i]+array[i+1]>=randnum){//bug fix on July 3, 2016. add +array[i]
				index = i+1;
				break;
			}
			else
			{
				sum = sum+array[i];
			}

		}
	}
	return index;
}
