#!/bin/bash
# vector_calc.sh data_in_file data_out_file

export VECTOR_DATA_IN=$1
export VECTOR_DATA_OUT=$2

# /testbed/tools/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -nojvm -batch "run('/home/ec2-user/code_kim/VECTOR/api/matlab_io.m');exit;"
/testbed/tools/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -nojvm -batch "run('/opt/python/current/app/matlab_io.m');exit;"

