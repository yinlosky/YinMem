#!/bin/bash
cd /home/yhuang9/my_workspace/AlluxioEigenSolver_v2;
matlab -nodesktop -nosplash -r "initRandVectorB($1,$2,$3); exit;" &>> initBlog.out &

