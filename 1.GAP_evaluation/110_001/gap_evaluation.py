# Evaluate LAMMPS dump file

from ovito.io import import_file, export_file
import ovito.data
import sys
import os
import numpy as np
import math

# Read path and name for the dump file.
try:
    dumpfile = sys.argv[1]
    sif = sys.argv[2]
    criterion = float(sys.argv[3])

except IndexError:
    print("Please enther the dump file and the corresponding K value!")
    print("****************************")
    print("Example usage: python file.dump sif")
    print("****************************")
    exit()

###################################################################
# Define functions that are used to evaluate GAP uncertainty.
###################################################################
#------------------------------------------------------------------
# Convert Dump file to XYZ format, multiple frame is allowed.
#------------------------------------------------------------------
def dump2xyz(dumpfile):
    pipeline = import_file(dumpfile)
    export_file(pipeline, './coords_temp/coords.xyz.*', 'xyz', columns=['Particle Type', 'Particle Identifier', 'Position.X', 'Position.Y', 'Position.Z'], multiple_frames=True)
    
    return pipeline.source.num_frames

#------------------------------------------------------------------
# Convert XYZ format to Extended XYZ that is readable for QUIP 
# --> to evaluate the GAP local variance.
#------------------------------------------------------------------
def xyz2gapExtxyz(xyzfile, Extxyzfile):
    with open(xyzfile, 'r') as f:
        line1 = f.readline()
        line2 = f.readline()
        lines = f.readlines()
        f.close()
    lines = [ 'Fe'+line[1:] for line in lines]
    with open(Extxyzfile, 'w') as f:
        f.writelines(line1)
        f.writelines(line2)
        f.writelines(lines)
        f.close()

#------------------------------------------------------------------
# Evaluation GAP predicted uncertainty.
#------------------------------------------------------------------
def gapevaluation(inpxyzfile, quipoutfile):
    # assess the error index by QUIP
    cmd='./local_variances.sh ' + inpxyzfile +' > quip_out'
    os.system(cmd)
    # output .xyz file with gap error
    lineout = []
    with open('quip_out', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == 'A':
                line = line[3:]
                lineout.append(line)
        f.close()

    with open(quipoutfile, 'w') as f:
        f.writelines(lineout)
        f.close()
#------------------------------------------------------------------     
# Function returns N largest elements
#------------------------------------------------------------------
def Nmaxelements(list1, N):
    final_list = []
    for i in range(0, N): 
        max1 = 0      
        for j in range(len(list1)):     
            if list1[j] > max1:
                max1 = list1[j];  
        list1.remove(max1);
        final_list.append(max1)
    
    return final_list

#------------------------------------------------------------------
# Find ten max error index for every frame of crack simulation
# For the first frame of each K value.
#------------------------------------------------------------------
def findmaxerror(quipoutfile, outfile, output_error_all, output_error_k, sif, framenum):
    cutoff = 40 # Only the atoms inside the cutoff will be considered for the error evaluation.
    gap_error = []
    gap_errorWithID = []
    output = open(outfile, 'a')
    f = open(quipoutfile, 'r')
    num_lines1 = f.readline()
    num_lines = int(num_lines1)
    sec_lin = f.readline().split()
    sec_lin[-1] = 'Properties=species:S:1:id:I:1:pos:R:3:local_gap_variance:R:1:sqrt_gap_variance:R:1'
    sec_inf = ' '.join([str(item) for item in sec_lin]) + '\n'
    output.write(num_lines1)
    output.write(sec_inf)
    for i in range(0,num_lines):
        line = f.readline().split()
        pos_x = float(line[1])
        pos_y = float(line[2])
        pos_z = float(line[3])
        atomid = int(line[4])
        gap_var = float(line[-1])
        gap_sqrtVar = np.sqrt(gap_var)
        temp = 'Fe' + ' ' + line[4] + ' ' + line[1] + ' ' + line[2] + ' ' + line[3]  + ' ' + line[-1]  + ' ' + str(gap_sqrtVar) + '\n'
        output.write(temp)
        if pos_x**2 + pos_y**2 < cutoff**2: 
            gap_error.append(gap_sqrtVar)
    f.close()
    
    # Here we can find the atomID with the largest error and construct DFT cell.
    
    
    # Save max ten to output error file.
    Maxten = Nmaxelements(gap_error,10)
    # 
    if Maxten[0] > criterion:
        f = open('./gap_error/extrapolated_frame', 'a')
        f.write(sif +' ' + str(framenum) + ' ' + str(Maxten[0]) + ' ')
        f.write('\n')
        f.close() 
    # save max ten to output error file.    
    f = open(output_error_all, 'a') 
    for i in range(0,10):
        f.write(str(Maxten[i]) + ' ')
    f.write('\n')
    f.close()
    if framenum == 0:
        f = open(output_error_k, 'a') 
        for i in range(0,10):
            f.write(str(Maxten[i]) + ' ')
        f.write('\n')
        f.close()        

###################################################################
# Execute the code 
###################################################################
# Convert LAMMPS dump to xyz file, store them in ./coords_temp temporarily
# Get number of frames
nframes = dump2xyz(dumpfile)

# Define the name for final output files
# Output xyz file with GAP error
outfile = './gap_error/gap_erro_xyz.' + sif
# Output fisrt ten maximum errors around crack tip for each frame
output_error_all = './gap_error/output_error_all'
# Output fisrt ten maximum errors around crack tip for each K
output_error_k = './gap_error/output_error_k'
for i in np.arange(0,nframes):
    print('Evaluating for K=', int(sif)*0.01 , 'frame', i)
    xyzfile = './coords_temp/coords.xyz.' + str(i)
    # convert to extended XYZ for QUIP
    xyz2gapExtxyz(xyzfile, 'EXT_xyz')
    # GAP evaluation
    gapevaluation('EXT_xyz', 'quipout_xyz')
    # Prepare the final data
    findmaxerror('quipout_xyz', outfile, output_error_all, output_error_k, sif, i)
