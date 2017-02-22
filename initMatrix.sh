#!/bin/bash
# This script will store the data from Accumulo table to TFS.
# the input matrix using $1 machines; $2 processors in the cluster, $3 is the number of iterations;
# usage: ./ParallelReading.sh 16 32 10 >> parallelReading_2_18_105edges.log
cd /home/yhuang9/my_workspace/AlluxioEigenSolver_v2;
matlab -nodesktop -nosplash -r "initMatrix($1,$2,$3,$4); exit;" &>> myinitmatrix.out &
