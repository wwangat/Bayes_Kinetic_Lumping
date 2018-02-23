#for the 3-basin, 300 microstate example, run the algorithm by $bash run_3state_numerical.sh in the terminal.
g++ main_v20.cpp 
resultdir=test_on_rugged_potential/;echo traj_len.txt | ./a.out MicroAssignment.txt $resultdir 300 3 1.0 | tee $resultdir/logfile
