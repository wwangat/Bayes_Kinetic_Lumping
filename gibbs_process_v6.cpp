/*
name:gibbs_process_v5.cpp
revised on Feb23 2017
fast versionß
 */

//in this version, all the elements in the microstate TCM are used, and we don't care whether TCM is symmetrized or not.
int initialization(int *mapping, int sizeMicro, double **macrotCount, int **macronCount, int sizeMacro, double **data, int **find_array, int *find_size){
	int i, j, row, col;
	for (i=0;i<sizeMicro;i++){
		mapping[i] = randi(1, sizeMacro)-1;
	}
	//calculate macro macromacrotCount and macromacronCount
	for (i=0;i<sizeMacro; i++){
		for (j=0;j<sizeMacro;j++){
			macrotCount[i][j] = 0.0;
		}
	}
	for (i = 0; i<sizeMicro; i++){
		j = mapping[i];
		find_array[j][find_size[j]] = i;
		find_size[j]++;
	}
	for (i = 0; i<sizeMacro; i++){
		for(j = 0;j<sizeMacro;j++){
			macronCount[i][j] = find_size[i]*find_size[j];
			for(row=0; row<find_size[i]; row++){
				for(col=0; col<find_size[j];col++){
					macrotCount[i][j] += data[find_array[i][row]][find_array[j][col]];
				}
			}
		}
	}
	return 0;
}

int update(int *mapping, int sizeMicro, double **macrotCount, int **macronCount, int sizeMacro, double **data, int microIndex, double alpha_nd, double beta_nd, double alpha_dd, double beta_dd, int **find_resultArray, int *find_resultSize){
	//diagonal: alpha_dd, beta_dd; non-diagonal:alpha_nd, beta_nd
	int k, **old_macronCount, old_index, new_index, macroIndex, i, j, m, n, l, row, col, tempint, *B1, *B2;
	double sum, *A1, *A2;
	int MC_index;
	double temp_new, temp_old, deduct_cond;
	double **old_macrotCount, *cond_prob, temp;
	old_macrotCount = (double **)malloc(sizeMacro*sizeof(double *));
	old_macronCount = (int **)malloc(sizeMacro*sizeof(int *));
	cond_prob = (double *)malloc(sizeMacro*sizeof(double));
	A1 = (double *)malloc(sizeMacro*sizeof(double));
	B1 = (int *)malloc(sizeMacro*sizeof(int));
	A2 = (double *)malloc(sizeMacro*sizeof(double));
	B2 = (int *)malloc(sizeMacro*sizeof(int));
	for (i=0;i<sizeMacro;i++){
		old_macrotCount[i] = (double *)malloc(sizeMacro*sizeof(double));
		old_macronCount[i] = (int *)malloc(sizeMacro*sizeof(int));
	}
	old_index = mapping[microIndex];
	//
	deduct_cond = 0.0;
	//posterior, old-index and new-index part
	for(l=0;l<sizeMacro;l++){
		if (l != old_index){
			deduct_cond -= lgamma(macrotCount[old_index][l]+alpha_nd)-(macrotCount[old_index][l]+alpha_nd)*log(macronCount[old_index][l]+beta_nd);
			deduct_cond -= lgamma(macrotCount[l][old_index]+alpha_nd)-(macrotCount[l][old_index]+alpha_nd)*log(macronCount[l][old_index]+beta_nd);
			//thus we have considered all the elements involved
		}
		else{
			deduct_cond -= lgamma(macrotCount[old_index][old_index]+alpha_dd)-(macrotCount[old_index][old_index]+alpha_dd)*log(macronCount[old_index][old_index]+beta_dd);
			//we also removed the diagonal parts by using its corresponding alpha and beta
		}
	}
	macrotCount[old_index][old_index] -= data[microIndex][microIndex];
	macronCount[old_index][old_index]--;
	for (l = 0; l<sizeMacro; l++)
	{
		if (l == old_index){
			B1[l] = find_resultSize[l]-1;
			sum = 0.0;
			for (m = 0; m<find_resultSize[l];m++){
				if (find_resultArray[l][m] != microIndex)	sum += data[microIndex][find_resultArray[l][m]];
			}
			A1[l] = sum;
		}
		else{
			B1[l] = find_resultSize[l];//all the microstates index that belongs to this macrostate
			if (B1[l] != 0){
				sum = 0.0;
				for (m = 0; m<find_resultSize[l];m++){
					sum += data[microIndex][find_resultArray[l][m]];
				}
				A1[l] = sum;
			}
			else{
				A1[l] = 0.0;
			}
		}	//link table, here A is the summation of the elements from microstate 1 to other microstates that belongs to macrostate 2
		//add on Feb23, 2017
		macrotCount[old_index][l] -= A1[l];
		macronCount[old_index][l] -= B1[l];
	}
	/////////////////add on Feb23, 2017
	for (l = 0; l<sizeMacro; l++)
	{
		if (l == old_index){
			B2[l] = find_resultSize[l]-1;
			sum = 0.0;
			for (m = 0; m<find_resultSize[l];m++){
				if (find_resultArray[l][m] != microIndex)	sum += data[find_resultArray[l][m]][microIndex];
			}
			A2[l] = sum;
		}
		else{
			B2[l] = find_resultSize[l];//all the microstates index that belongs to this macrostate
			if (B2[l] != 0){
				sum = 0.0;
				for (m = 0; m<find_resultSize[l];m++){
					sum += data[find_resultArray[l][m]][microIndex];
				}
				A2[l] = sum;
			}
			else{
				A2[l] = 0.0;
			}
		}	//link table, here A is the summation of the elements from microstate 1 to other microstates that belongs to macrostate 2
		//add on Feb23, 2017
		macrotCount[l][old_index] -= A2[l];
		macronCount[l][old_index] -= B2[l];
	}

	for (m = 0; m<sizeMacro; m++){
		for (n = 0; n<sizeMacro; n++){
			old_macrotCount[m][n] = macrotCount[m][n];
			old_macronCount[m][n] = macronCount[m][n];
		}
	}
	 //having changed part of old_index, now not symmetric temporatory

	for(l=0;l<sizeMacro;l++){
		if (l != old_index){
			deduct_cond += lgamma(macrotCount[old_index][l]+alpha_nd)-(macrotCount[old_index][l]+alpha_nd)*log(macronCount[old_index][l]+beta_nd);
			deduct_cond += lgamma(macrotCount[l][old_index]+alpha_nd)-(macrotCount[l][old_index]+alpha_nd)*log(macronCount[l][old_index]+beta_nd);
			//thus we have considered all the elements involved
		}
		else{
			deduct_cond += lgamma(macrotCount[old_index][old_index]+alpha_dd)-(macrotCount[old_index][old_index]+alpha_dd)*log(macronCount[old_index][old_index]+beta_dd);
			//we also removed the diagonal parts by using its corresponding alpha and beta
		}
	}
////////////////having dealt with the old index parts
	cond_prob[old_index] = 0.0;
	for (macroIndex = 0;macroIndex<sizeMacro;macroIndex++){
		new_index = macroIndex;
		//updating
		if (old_index != new_index){
			cond_prob[new_index] = deduct_cond;
			for (l = 0; l<sizeMacro;l++){
				if (new_index != l){
					cond_prob[new_index] -= lgamma(macrotCount[new_index][l]+alpha_nd)-(macrotCount[new_index][l]+alpha_nd)*log(macronCount[new_index][l]+beta_nd);
					cond_prob[new_index] -= lgamma(macrotCount[l][new_index]+alpha_nd)-(macrotCount[l][new_index]+alpha_nd)*log(macronCount[l][new_index]+beta_nd);
				}
				else{
					cond_prob[new_index] -= lgamma(macrotCount[new_index][new_index]+alpha_dd)-(macrotCount[new_index][new_index]+alpha_dd)*log(macronCount[new_index][new_index]+beta_dd);
				}
			}
			macrotCount[new_index][new_index] += data[microIndex][microIndex];
			macronCount[new_index][new_index]++;
			for (l = 0; l<sizeMacro; l++)
			{
				if (new_index!=l){
					macrotCount[new_index][l] += A1[l];
					macronCount[new_index][l] += B1[l];
					macrotCount[l][new_index] += A2[l];
					macronCount[l][new_index] += B2[l];
				}
				else{
					macrotCount[new_index][l] += A1[l]+A2[l];
					macronCount[new_index][l] += B1[l]+B2[l];   //need to check on this
				}
			}

			for(l=0; l<sizeMacro;l++){
				if (new_index != l){
					cond_prob[new_index] += lgamma(macrotCount[new_index][l]+alpha_nd)-(macrotCount[new_index][l]+alpha_nd)*log(macronCount[new_index][l]+beta_nd);
					cond_prob[new_index] += lgamma(macrotCount[l][new_index]+alpha_nd)-(macrotCount[l][new_index]+alpha_nd)*log(macronCount[l][new_index]+beta_nd);
				}
				else{
					cond_prob[new_index] += lgamma(macrotCount[new_index][new_index]+alpha_dd)-(macrotCount[new_index][new_index]+alpha_dd)*log(macronCount[new_index][new_index]+beta_dd);
					//we also removed the diagonal parts by using its corresponding alpha and beta
				}
			}
			for (l = 0; l<sizeMacro; l++){
				macrotCount[new_index][l] = old_macrotCount[new_index][l];
				macronCount[new_index][l] = old_macronCount[new_index][l];
				macrotCount[l][new_index] = old_macrotCount[l][new_index];
				macronCount[l][new_index] = old_macronCount[l][new_index];
			}  //this is for changing the microstate index to other macrostates
		}

		//updating finished
		//calculate posterior probability
	}
	//updating through Monte Carlo
	temp = max(cond_prob, sizeMacro);
	for (k = 0;k<sizeMacro; k++){
		cond_prob[k] = cond_prob[k]-temp;
		cond_prob[k] = exp(cond_prob[k]);
	}
	MC_index = MonteCarlo(cond_prob, sizeMacro);
	free(cond_prob);
	mapping[microIndex] = MC_index;
	for(int m = 0; m<sizeMacro; m++){
		for (int n = 0; n<sizeMacro; n++){
			macrotCount[m][n] = old_macrotCount[m][n];
			macronCount[m][n] = old_macronCount[m][n];
		}
	}
	macrotCount[MC_index][MC_index] += data[microIndex][microIndex];
	macronCount[MC_index][MC_index]++;

	for (l = 0; l<sizeMacro; l++)
	{

		if (MC_index != l)
		{
			macrotCount[MC_index][l] += A1[l];
			macronCount[MC_index][l] += B1[l];
			macrotCount[l][MC_index] += A2[l];
			macronCount[l][MC_index] += B2[l];
		}
		else
		{
			macrotCount[l][MC_index] += A1[l]+A2[l];
			macronCount[l][MC_index] += B1[l]+B2[l];
		}
	}

	for (i = 0; i<sizeMacro; i++){
		free(old_macrotCount[i]);
		free(old_macronCount[i]);
	}
	free(old_macrotCount);
	free(old_macronCount);
	free(A1);
	free(B1);
	free(A2);
	free(B2);
	//updating find_resultArray and find_resultSize
	if (MC_index != old_index){
		int *temp_array;
		temp_array = (int *)malloc(sizeMicro*sizeof(int *));
		j = 0;
		for(int i = 0; i<find_resultSize[old_index]; i++){
			if (find_resultArray[old_index][i] != microIndex){
				temp_array[j] = find_resultArray[old_index][i];
				j++;
			}
		}
		find_resultSize[old_index]--;
		for(int i = 0; i<find_resultSize[old_index];i++){
			find_resultArray[old_index][i] = temp_array[i];
		}
		find_resultSize[MC_index]++;
		find_resultArray[MC_index][find_resultSize[MC_index]-1] = microIndex;
		free(temp_array);
	}
	return 0;
}
