#!/usr/bin/bash
# lei.zhang@rug.nl, 14 March 2022,
# 

# Controlling script of active learning process.
# Please read README.md before using the script.

# set the error criterion for stopping the simulation
criterion="0.01"
# loop the whole active learning process
# Will be stopped once the desired convergence is achieved.
iter=0
eaddress="lei.zhang@rug.nl"
while [ 1 -gt 0 ]
do
# increase iter index
let "iter++"
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
##########################################################
# Section: 0
# Folder: 0.K_test            
# Run crack simulation with current GAP potential
# Code: LAMMPS
##########################################################
#---------------------------------------------------------
# Create 'potential.in' for current GAP
# The potential xml_label need to be extracted from GAP potential file
#---------------------------------------------------------
# Get the label of current GAP model
gaplabstr=`awk 'NR==1 {print $1}' ./GAP_model/gap_Fe_test.xml`
gaplabel="${gaplabstr:1:-1}"

cat > ./0.K_test/potential.in <<EOF
pair_style      quip
pair_coeff      * * ../../GAP_model/gap_Fe_test.xml   "Potential xml_label=$gaplabel" 26 26
EOF
# Set lammps executable path
# You can run the job on local machine or 
# Customise your own SLURM script to submit a job: ./0.K_test/submit
#--------------------------------------------
# Submit LAMMPS JOBs
#--------------------------------------------
cd ./0.K_test/100_011
sbatch submit > crack_jobID1
# get the job IDs -> to set the job dependencies
c_jobID1=`awk 'NR==1 {print $4}' crack_jobID1`
cd ../110_001
sbatch submit > crack_jobID2
# get the job IDs -> to set the job dependencies
c_jobID2=`awk 'NR==1 {print $4}' crack_jobID2`
# back to main folder
cd ../../
#=============================================

##########################################################
# Section: 1
# Folder: 1.GAP_evaluation
# GAP predicted error evaluation.
# code: QUIP + PYTHON
# The crack tip is extracted once the maximum error exceed a preset threshold (Here we set it to 10 meV).
# A random small displacements with a Gaussian distribution is applied to all atoms along three directions. 
##########################################################
# Evaluate every dump file produced in Section-0.
# Evaluation needs to be done for two crack systems.
# Start the job after crack simulation is finished.
cd ./1.GAP_evaluation/100_011

sbatch --dependency=afterany:${c_jobID1}:${c_jobID2} submit > quip_jobID_1
q_jobID1=`awk 'NR==1 {print $4}' quip_jobID_1`
cd ../110_001/
sbatch --dependency=afterany:${c_jobID1}:${c_jobID2} submit > quip_jobID_2
q_jobID2=`awk 'NR==1 {print $4}' quip_jobID_2`
cd ../../
#=============================================


# set a flag to detect whether the quip evaluation is finished.
# If all jobs are finished, jump out of the while loop.
stop_flag1=0
echo "Waiting for crack simulation and GAP evaluation."
while [ ${stop_flag1} == 0 ]
do
    squeue -j ${q_jobID1} > quip_jobstate1
    squeue -j ${q_jobID2} > quip_jobstate2
    numlineq1=`wc -l < quip_jobstate1`
    numlineq2=`wc -l < quip_jobstate1`
    
    if [ ${numlineq1} == 2 ] && [ ${numlineq2} == 2 ]; then
        sleep 60
    else
        stop_flag1=1
    fi
done

##########################################################
# Section: 2
# Folder: 2.DFT_calculation
# Run DFT calculation of small crack tip
# code: Quantum Espresso + BASH
##########################################################
##########################################################
# Flag used to break the active learning. -> stop the script
# once the desired accuracy is achieved.
##########################################################
#=============================================
if [ -f ./1.GAP_evaluation/100_011/gap_error/extrapolated_frame ] && [ -f ./1.GAP_evaluation/110_001/gap_error/extrapolated_frame ]; then
    cd ./2.DFT_calculation/100_011
    sbatch submit > dft_jobID_1
    d_jobID1=`awk 'NR==1 {print $4}' dft_jobID_1`	
    cd ../110_001/
    sbatch submit > dft_jobID_2
    d_jobID2=`awk 'NR==1 {print $4}' dft_jobID_2`
    cd ../../
    dftjobnumber=2
elif [ -f ./1.GAP_evaluation/100_011/gap_error/extrapolated_frame ]; then
    cd ./2.DFT_calculation/100_011
    sbatch submit > dft_jobID
    d_jobID=`awk 'NR==1 {print $4}' dft_jobID`
    dftjobnumber=1
    cd ../../
elif [ -f ./1.GAP_evaluation/110_001/gap_error/extrapolated_frame ]; then
    cd ./2.DFT_calculation/110_001
    sbatch submit > dft_jobID
    d_jobID=`awk 'NR==1 {print $4}' dft_jobID`
    dftjobnumber=1
    cd ../../
else
    # Jump out of the controlling loop
    # Stop the training iterations.
    dftjobnumber=0
mail -s "Active learning is Stopped!" "${eaddress}"<<EOF
GAP is converged to the desired accuracy: GAP uncertainty less than ${criterion} eV.
EOF
    break
fi
#=============================================
##########################################################
# Section: 3
# Folder: ./3.NEWGAP_training
# Train new GAP with newly added crack tip DFT data.
# code: QUIP
##########################################################
if [ ${dftjobnumber} == 2 ]; then
    cd ./3.NEWGAP_training
    sbatch --dependency=afterany:${d_jobID1}:${d_jobID2} submit > training_jobID
else
    cd ./3.NEWGAP_training
    sbatch --dependency=afterany:${d_jobID} submit > training_jobID
fi

t_jobID=`awk 'NR==1 {print $4}' training_jobID`
cd ..
# set a flag to detect whether the quip evaluation is finished.
# If all jobs are finished, jump out of the while loop.
# Check the status of both simulation
echo "Waiting for DFT calculation."
stop_flag2=0
while [ ${stop_flag2} == 0 ]
do
    squeue -j ${t_jobID} > training_jobstate
    numlineq=`wc -l < training_jobstate`
    if [ ${numlineq} == 2 ]; then
        sleep 60
    else
        stop_flag2=1
    fi
done

##########################################################
# Section: 4
# Folder: ./DATA
# code: BASH
# This section is used to backup all the data, clear the folders
# and initialize to restart the process. 
##########################################################

savedata="./DATA/${iter}"
mkdir ${savedata}
# copy LAMMPS dump file
cp -r ./0.K_test/100_011/dump/ ${savedata}/100_011_lmp_dump
cp -r ./0.K_test/110_001/dump/ ${savedata}/110_001_lmp_dump
mv ./0.K_test/100_011/slurm* ${savedata}/100_011_lmp_dump
mv ./0.K_test/110_001/slurm* ${savedata}/110_001_lmp_dump
# clean folder
rm ./0.K_test/100_011/dump/*
rm ./0.K_test/110_001/dump/*
rm ./0.K_test/potential.in

# copy GAP evaluation
cp -r ./1.GAP_evaluation/100_011/gap_error ${savedata}/100_011_gap_error
cp -r ./1.GAP_evaluation/110_001/gap_error ${savedata}/110_001_gap_error
cp -r ./1.GAP_evaluation/100_011/qe_inp ${savedata}/100_011_qe_inp
cp -r ./1.GAP_evaluation/110_001/qe_inp ${savedata}/110_001_qe_inp
mv ./1.GAP_evaluation/100_011/slurm* ${savedata}/100_011_gap_error
mv ./1.GAP_evaluation/110_001/slurm* ${savedata}/110_001_gap_error
# Those folders are already cleaned!

# copy DFT data
cp -r ./2.DFT_calculation/100_011 ${savedata}/100_011_dft
cp -r ./2.DFT_calculation/110_001 ${savedata}/110_001_dft
rm ./2.DFT_calculation/100_011/slurm*
rm ./2.DFT_calculation/100_011/crack*
rm ./2.DFT_calculation/100_011/DFT_crack_100_011
rm ./2.DFT_calculation/110_001/slurm*
rm ./2.DFT_calculation/110_001/crack*
rm ./2.DFT_calculation/110_001/DFT_crack_110_001

# copy GAP model: copy the latest version to ./GAP_model
# move current one as back up
mkdir ${savedata}/GAP_model_old
mv ./GAP_model/gap_Fe_test* ${savedata}/GAP_model_old
cp ./3.NEWGAP_training/gap_Fe_test* ./GAP_model
cp -r ./3.NEWGAP_training ${savedata}/GAP_model_new
# clean folder
rm ./3.NEWGAP_training/slurm*
rm ./3.NEWGAP_training/gap_Fe_test*
rm ./3.NEWGAP_training/train.xyz.idx

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
done
