/*
 * ============================================================================
 *       Filename:  main_v18.cpp
 *    Description:  Bayes lumping algorithm
 *			 Function:  Input MicroAssignment, Output MacroAssignment.  Doing lumping based on Bayes posterior probability,
										select parameters based on metastability criteria.  Can deal with different population level.
										Finally select the assignment with biggest metastability
 *        Created:  2016-2-19 00:58
 *       Modified:  2016-12-29 16:19 based on main_v9.cpp, new function: add jumping window
 *       Modified:  2017-1-1 4:33am, remove empty states as well as select most populated states based on the transition count matrix
 *		Modified: 2017-2-23 12:58pm, add the diagonal elements of the block, do not symmetrize the transition count matrix
 *    Modified: 2017-03-20: use the average of estimated transition in this version to decide whether to throw away the lumping
 *          Author:  Wei WANG        (wwangat@gmail.com)
 * ============================================================================
*/

///////////////////////////////////////////////////////////////////////////////
//compile: g++ main_v20.cpp -o main.o 
//usage: ./main.o inputfile outputdir nMicro nMacro pop
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "allocate.cpp"
#include "read.cpp"
#include "operator.cpp"
#include "msm_clean_v2.cpp"
#include "basic.cpp"
#include <iostream>
#include "gibbs_process_v6.cpp"
using namespace std;
const int MCstep = 200;
const int numBeta = 10;   //this number should be changed, as well as the max_beta and
const int RunTime = 100;
const int RunTime_test = 20;
const int traj_num = 100;//number of trajectories
const double ratio_mean_dd_nd = 1e10; //1e10 //beta_dd/beta_nd, for general example
const double ratio_var_dd_nd = 1e4;//1e4, ratio_mean_dd_nd and ratio_var_dd_nd are related to the free energy landscape

int main(int argc, char *argv[]){
	//input parameters.........	................////////////////////////////
	int lagtime=80; //set as microstate Markovian lagtime
    int jump_step=1;//for jumping window, these two parameters should revised
    double max_beta = 1e10, min_beta= 1e6;   //nolinear inteval

	//inputs: MicroAssignment filename, outputdir
	char *inputfile, *outputdir;
	inputfile = argv[1];
	outputdir = argv[2];
	struct stat st={0};
	if (stat(outputdir, &st) == -1){
		mkdir(outputdir, 0700);
	}
	int nMicro = atoi(argv[3]);
	int	nMacro = atoi(argv[4]);
	double pop = atof(argv[5]);   //we can only use the most populated states for analysis

	clock_t start, end;
	double cpu_time_used;
	//variable section
	int nline = countline(inputfile);
	int *micro = alloarray_int(nline), *macro=alloarray_int(nline);
	readarray(inputfile, micro); //read microassignment to array

	int temp1, temp2;
	temp1 = max(micro, nline);
	if(temp1 != nMicro-1){
			cout<<"max of microstate trajectory index is "<<temp1<<", regulating"<<endl;
			temp2 = min(micro, nline);
			for(int jj = 0; jj<nline; jj++){
					micro[jj] = (nMicro-1)*(micro[jj]-temp2)/(temp1-temp2);
			}
	}//now the microstate trajectory index is starting from 0 to nMicro-1

	////////////////////////////////calculating microstate transition count matrix in order to remove empty states as well as select most populated states///
	double **old_CountMatrix;
	old_CountMatrix = (double **)malloc(nMicro*sizeof(double *));
	for (int i=0;i<nMicro;i++){
		old_CountMatrix[i] = (double *)malloc(nMicro*sizeof(double));
	}
	//calculate CountMatrix
	int *traj_len = alloarray_int(traj_num);
	if (traj_num != 1){
		printf("Please indicate the trajectory length file\n");
		char file1[4096];
		scanf("%s", file1);
		FILE *fp=fopen(file1, "r");
		for (int j=0; j<traj_num;j++){
			fscanf(fp, "%d", &traj_len[j]);
		}
		fclose(fp);
	}
	else{
		traj_len[0] = 1; //actually no use in this case
	}

	transCount(old_CountMatrix, nMicro, lagtime, nline, micro, traj_len, traj_num, jump_step);   //get the updated CountMatrix
	//this is used to calculate the stationary population

	///////////////////////check whether contains empty states by checking the microstate CountMatrix////////////////////
	int *select_index;
	//define temp variables
	int *temp_micro = alloarray_int(nline), temp_int;
	for(int j=0;j<nline;j++){
		temp_micro[j] = -1;
	}
	///////////////////////////////////////////////added on Jan 1, 2017 by Wang Wei
	double *temp_sum = alloarray_double(nMicro);
	for(int j=0;j<nMicro;j++){
		temp_sum[j] = 0.0;
	}
	for(int j=0;j<nMicro;j++){
		for (int k =0; k<nMicro; k++){
			temp_sum[j] = temp_sum[j]+old_CountMatrix[j][k];
		}
	}

	double cutpop=temp_sum[0];
	for (int j=1;j<nMicro;j++){
		cutpop = cutpop+temp_sum[j];
	}
	cutpop = cutpop*pop;   ///overall conformations
	int *sorted_index;
	sorted_index = alloarray_int(nMicro);
	for(int j=0;j<nMicro;j++){
		sorted_index[j] = j;
	}
	quicksort(temp_sum, nMicro, sorted_index); //since we want to find most populated states, we need to sort in 'descend order', temp_sum is the population of each microstate.
	//cummulate summation
	cumsum(temp_sum, nMicro);
	//cut according to cutpop
	int temp = 0;
	for(int j=1;j<nMicro;j++){
		if(temp_sum[j]>=cutpop){
			temp = j; //number of most populated states
			break;
		}
	}
	temp = temp+1;
	for(int j=0;j<nMicro;j++){
		free(old_CountMatrix[j]);
	}
	free(old_CountMatrix);

	cout<<"contain "<<temp<<" states due to selection of most populated states or remove empty states"<<endl;
	nMicro = temp;
	select_index = alloarray_int(nMicro);
	for(int j=0;j<nMicro;j++){
		select_index[j] = sorted_index[j];
	}
	quicksort_int(select_index, nMicro);
	reverse(select_index, nMicro); //original select_index is from populated to less populated; we put them in 'ascend' order of index
	for(int j=0;j<nMicro;j++){
		temp_int = select_index[j];
		for(int k=0;k<nline;k++){
			if(micro[k] == temp_int){
				temp_micro[k] = j;
			}
		}
	}
	for(int j=0;j<nline;j++){
		micro[j] = temp_micro[j];
	}

	free(temp_sum);
	free(select_index);
	//Now begin Bayes Lumping
	double evaluation;  //evaluation methods can be metastability, posterior function or modularity
	double totaledge, alpha_nd=10, optibeta, **tCount;
	//add
	double alpha_dd = ratio_mean_dd_nd*ratio_mean_dd_nd/ratio_var_dd_nd*alpha_nd;
	double beta_dd;
	//add end
	int i,j, k, m, n, row, col, success_iter, microIndex, flag, ITER,*mapping, **nCount, *best_mapping;//success_iter:successful trial gibbs sampling
	double optimodu=-1e20, tempdouble;
	int *find_resultSize, **find_resultArray;
	double average;
	double beta_nd[numBeta];
	double ratio=pow(max_beta/min_beta, 1.0/numBeta);
	for (int j = 0;j<numBeta;j++){
			 beta_nd[j] = floor(min_beta*pow(ratio, j));
	}

	srand(time(NULL));
	start=clock();
	double **CountMatrix;//the new microstate transition count matrix
	CountMatrix = (double **)malloc(nMicro*sizeof(double *));
	for (int j=0;j<nMicro;j++){
		CountMatrix[j] = (double *)malloc(nMicro*sizeof(double));
	}
	transCount(CountMatrix, nMicro, lagtime, nline, micro, traj_len, traj_num, jump_step, 0);
	//0 means not symmetrize, 1 means to do symmetrize
	cout<<"first and last elements in new CountMatrix"<<endl;
	cout<<CountMatrix[0][0]<<'\t'<<CountMatrix[nMicro-1][nMicro-1]<<endl;
	mapping = (int *)malloc(nMicro*sizeof(int));
	best_mapping = (int *)malloc(nMicro*sizeof(int));
	tCount = (double **)malloc(nMacro*sizeof(double *));
	nCount = (int **)malloc(nMacro*sizeof(int *));
	for (k=0;k<nMacro;k++){
		tCount[k] = (double *)malloc(nMacro*sizeof(double));
		nCount[k] = (int *)malloc(nMacro*sizeof(int));
	}
	find_resultSize = (int *)malloc(nMacro*sizeof(int));
	find_resultArray = (int **)malloc(nMacro*sizeof(int *));
	for (i = 0; i<nMacro; i++){
		find_resultArray[i] = (int *)malloc(nMicro*sizeof(int));
	}
	//scan paramter sets and use modularity to choose beta
	optibeta=beta_nd[0];
	cout<<"Now begin to scan the parameters"<<endl;
	double optimodu_temp;
	for (i = 0; i<numBeta; i++){
		printf("begin a new parameter %lf\n", beta_nd[i]);
		beta_dd =ratio_mean_dd_nd/ratio_var_dd_nd*beta_nd[i];
		cout<<"non diag alpha part:"<<alpha_nd<<", diag alpha is:"<<alpha_dd<<endl;
		cout<<"non diag beta part:"<<beta_nd[i]<<", diag beta is:"<<beta_dd<<endl;
		j=1;success_iter=0;
		average = 0.0;
		optimodu_temp = 0;
		double *meta_list = alloarray_double(RunTime_test);
		while (j<=RunTime_test && success_iter<1000){
			for (k = 0; k<nMacro; k++){
				find_resultSize[k] = 0;
			}
			flag=0;
			initialization(mapping, nMicro, tCount, nCount, nMacro, CountMatrix, find_resultArray, find_resultSize);
			double oldsum = 100;
			int old_index;
			for (ITER = 1;ITER<=MCstep;ITER++){
				for (microIndex = 0; microIndex<nMicro; microIndex++){
					update(mapping, nMicro, tCount, nCount, nMacro, CountMatrix, microIndex, alpha_nd, beta_nd[i], alpha_dd, beta_dd, find_resultArray, find_resultSize);
				}

				for (k = 0; k<nMacro; k++){
					if (find_resultSize[k] == 0){
						flag = 1; break;
					}
				}
				if (flag == 1) break;
				double newsum = 0.0;
				for (k = 0;k<nMacro;k++){
					newsum += tCount[k][k];
				}
				if (fabs(newsum-oldsum)<1.0){
					break;
				}
				else	oldsum = newsum;
			}
			//finished iteration
			j++;
			//judge whether the result is okay according to the average count of transition
			for (int k = 0;k<nMacro; k++){
				double temp_max=tCount[0][k];
				for(int m = 1;m<nMacro; m++){
					//baceuse column normalized
					if(tCount[m][k]>temp_max){
						temp_max = tCount[m][k];
					}
				}
				if(tCount[k][k]<temp_max)   //then not metastable
					{flag = 1; break;}
			}
			if (flag==0){
				success_iter++;
				evaluation = metastability(tCount, nMacro);
				meta_list[success_iter-1]=evaluation;
			}
		}
		double *temp_list = alloarray_double(success_iter);
		for(int j = 0;j<success_iter;j++){
			temp_list[j] = meta_list[j];
		}
		free(meta_list);
		quicksort(temp_list, success_iter); //since we want to find most populated states, we need to sort in 'descend order', temp_sum is the population of each microstate.
		optimodu_temp = temp_list[int(success_iter/2)];
		cout<<"success for "<<success_iter<<" iterations, "<< "median metastability is: "<<optimodu_temp<<endl;
		if (optimodu_temp>optimodu & success_iter>0){
			optimodu = optimodu_temp;
			optibeta = beta_nd[i];
		}
		free(temp_list);
	}
	cout <<"END SCANNING PARAMETERS, the optimal parameter is beta_nd=" << optibeta << ", with average metastability to be" <<optimodu << endl;
	cout <<"Now begin to run at the optimal parameters"<<endl;

	optimodu = -1e20;
	for (i = 0; i<1; i++){
		printf("begin a new parameter %lf\n", optibeta);
		beta_dd =ratio_mean_dd_nd/ratio_var_dd_nd*optibeta;
		cout<<"non diag alpha part:"<<alpha_nd<<", diag alpha is:"<<alpha_dd<<endl;
		cout<<"non diag beta part:"<<optibeta<<", diag beta is:"<<beta_dd<<endl;
		j=1;success_iter=0;
		average = 0.0;
		while (j<=RunTime && success_iter<100){
			for (k = 0; k<nMacro; k++){
				find_resultSize[k] = 0;
			}
			initialization(mapping, nMicro, tCount, nCount, nMacro, CountMatrix, find_resultArray, find_resultSize);
			cout<<"initial posterior:"<<fixed<<posterior_v2(tCount, nCount, alpha_nd, optibeta, alpha_dd, beta_dd, nMacro)<<endl;
			//get initial mapping, tCount, nCount
			//updating by System-scan gibbs sampling
			flag = 0;
			double oldsum = 100;
			int old_index;
			for (ITER = 1;ITER<=MCstep;ITER++){
				for (microIndex = 0; microIndex<nMicro; microIndex++){
					update(mapping, nMicro, tCount, nCount, nMacro, CountMatrix, microIndex, alpha_nd, optibeta, alpha_dd, beta_dd, find_resultArray, find_resultSize);
				}

				for (k = 0; k<nMacro; k++){
					if (find_resultSize[k] == 0){
						flag = 1; break;
					}
				}
				if (flag == 1) break;
				double newsum = 0.0;
				for (k = 0;k<nMacro;k++){
					newsum += tCount[k][k];
				}
				if (fabs(newsum-oldsum)<1.0){
					break;
				}
				else	oldsum = newsum;
			}
			j++;
			//judge whether the result is okay according to the average count of transition
			for (int k = 0;k<nMacro; k++){
				double temp_max=tCount[0][k];
				for(int m = 1;m<nMacro; m++){
					//because column normalized
					if(tCount[m][k]>temp_max){
						temp_max = tCount[m][k];
					}
				}
				if(tCount[k][k]<temp_max)   //then not metastable
					{flag = 1; break;}
			}
			if (flag==0){
				success_iter++;
				evaluation = posterior_v2(tCount, nCount, alpha_nd, optibeta, alpha_dd, beta_dd, nMacro);
				if (evaluation>optimodu){
						optimodu = evaluation;
						for (int k=0;k<nMicro;k++){
								best_mapping[k] = mapping[k];
						}
				}
				printf("successful run with posterior to be %lf :\n", evaluation);
				printf("corresponding metastability to be %lf and its final block count matrix is as follows:\n", metastability(tCount, nMacro));

				for (int kk = 0;kk<nMacro;kk++){
					for (int jj=0;jj<nMacro;jj++){
						cout<<tCount[kk][jj]<<'\t';
					}
					cout<<'\n';
				}

				for (int kk = 0;kk<nMacro;kk++){
					for (int jj=0;jj<nMacro;jj++){
						cout<<nCount[kk][jj]<<'\t';
					}
					cout<<'\n';
				}
			}
		}
	}

	for(int j=0;j<nline;j++){
		 macro[j] = -1;
	}
	for(int j=0;j<nMicro;j++){
		for(int k=0;k<nline;k++){
			if(micro[k] == j){
				macro[k] = best_mapping[j];
			}
		}
	}

	char fn[512];
	strcpy(fn, outputdir);
	strcat(fn, "MacroAssignment.txt");
		FILE *outputfile = fopen(fn, "w");
	if (!outputfile) exit(0);
	for (int k = 0; k<nline; k++){
		fprintf(outputfile, "%d\n", macro[k]);
	}
	fclose(outputfile);
	char fn1[512];
	strcpy(fn1, outputdir);
	strcat(fn1, "mapping.txt");
		FILE *outputfile1 = fopen(fn1, "w");
	if (!outputfile1) exit(0);
	for (k = 0; k<nMicro; k++){
		fprintf(outputfile1, "%d\n", best_mapping[k]);
	}
	fclose(outputfile1);
	free(find_resultSize);
	for (k = 0; k<nMacro; k++){
		free(tCount[k]);
		free(nCount[k]);
		free(find_resultArray[k]);
	}
	free(find_resultArray);
	free(tCount);
	free(nCount);
	free(mapping);
		free(best_mapping);
	free(micro);
	free(macro);
	free(temp_micro);
	printf("program finished, best result with posterior to be %lf, now post-processing by yourselves\n", optimodu);
	for(i=0;i<nMicro;i++){
		free(CountMatrix[i]);
	}
	free(CountMatrix);
	end=clock();
	cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
	printf("use a total of %lf seconds\n", cpu_time_used);
	return 0;
}
