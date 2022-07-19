#/bin/bash

# set the OPENMP cores for parallelization
export OMP_NUM_THREADS=32
# read input XYZ file
quip_inp_coords=$1
# give the path for quip command.
/home/marescaf/software/QUIP_OPENMP/build/linux_x86_64_gfortran_openmp/quip \
output_file=a.out \
verbosity=ANALYSIS \
atoms_filename=${quip_inp_coords} \
calc_args="local_gap_variance" \
param_filename=../../GAP_model/gap_Fe_test.xml \
E=T \ 
F=T \
L=T \
local_gap_variance=T 

