from ovito.io import import_file, export_file
import ovito.data
import sys
import os
import numpy as np
import math

# set the error criterion for stopping the simulation
criterion = 0.015

# Read path and name for the dump file.
try:
    dumpfile = sys.argv[1]
    outfile = sys.argv[2]

except IndexError:
    print("Please enther the dump file and the output file!")
    print("****************************")
    print("Example usage: python file.dump gap_variance.xyz")
    print("****************************")
    exit()

###################################################################
# Define functions that are used to evaluate GAP uncertainty.
# 
###################################################################
# Convert Dump file to XYZ format 
def dump2xyz(dumpfile):
    pipeline = import_file(dumpfile)
    export_file(pipeline, './coords_temp/coords.xyz.*', 'xyz', columns=['Particle Type','Position.X', 'Position.Y', 'Position.Z'], multiple_frames=True)

# Convert XYZ format to Extended XYZ that is readable for QUIP --> to evaluate the GAP local variance.
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
        
# Function returns N largest elements
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

def gapevaluation(inpxyzfile):
    # assess the error index by QUIP
    cmd='./local_variances.sh ' + inpxyz +' > quip_out'
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

    with open('gap_error.xyz', 'w') as f:
        f.writelines(lineout)
        f.close()

# Find ten max error index.
def findmaxerror(outfile):
    cutoff = 40 # Only the atoms inside the cutoff will be considered for the error evaluation.
    gap_error = []
    output = open(outfile, 'w')
    with open('gap_error.xyz', 'r') as f:
        num_lines1 = f.readline()
        num_lines = int(num_lines1)
        sec_lin = f.readline().split()
        sec_lin[-1] = 'Properties=species:S:1:pos:R:3:local_gap_variance:R:1:sqrt_gap_variance:R:1'
        sec_inf = ' '.join([str(item) for item in sec_lin]) + '\n'
        output.write(num_lines1)
        output.write(sec_inf)
        for i in range(0,num_lines):
            line = f.readline().split()
            pos_x = float(line[1])
            pos_y = float(line[2])
            pos_z = float(line[3])
            gap_var = float(line[-1])
            gap_sqrtVar = np.sqrt(gap_var)
            temp = 'Fe' + ' ' + line[1] + ' ' + line[2] + ' ' + line[3] + ' ' + line[-1]  + ' ' + str(gap_sqrtVar) + '\n'
            output.write(temp)
            if pos_x**2 + pos_y**2 < cutoff**2: 
                gap_error.append(gap_sqrtVar)
        f.close()
        
#####################################################
# Define functions that are used to extract crack tips
# and constructing Quantum Espresso input file
#####################################################
# function used to read single frame dump file. 
# only atom ID and positions are extracted.
# Total atoms, id, positions are returned.
def readdump(dumpfile):
    data = open(dumpfile,'r')
    lines = data.readlines()
    # check the number of frames contained in one file.
    # commentted, need to upgrade the code of multiple frame in dump file is used.
    '''
    frame = 0
    for line in lines:
        line = line.split()
        if line[0] == 'ITEM:':
            frame = frame + 1
    frame = int(frame/4)
    '''
    numlines = list(range(0, len(lines)))
    id_all = []
    position_all = []
    for i, line in zip(numlines, lines):
        try:
            line = line.split()
        except AttributeError:
            return
        # third line: number of atoms
        if i == 3:
            natoms = int(line[0])
        # box size
        if i == 5:
            xmin = float(line[0])
            xmax = float(line[1])
        if i == 6:
            ymin = float(line[0])
            ymax = float(line[1])
        if i == 7:
            zmin = float(line[0])
            zmax = float(line[1])
        if i == 8:
            properties = len(line) - 2
        # read information for all atoms
        if i > 8:
            id_i = int(line[0])
            position = []
            position.append(float(line[2]))
            position.append(float(line[3]))
            position.append(float(line[4]))
            id_all.append(id_i)
            position_all.append(position)
            boundaries = ([xmin,xmax,ymin,ymax,zmin,zmax])
            
    return natoms, id_all, position_all, boundaries

# Store the id in the extracted crack tip region
def find_id(atominfo, region):    
    natoms = atominfo[0]
    natoms_lst = list(range(0, natoms))
    pbc_id = []
    for i in natoms_lst:
        id_i = atominfo[1][i]
        coords = atominfo[2][i]
        
        # select atoms in side a certain region.
        # store the atomic ID of all atoms.
        if coords[0] > region[0] and coords[0] < region[1] and \
        coords[1] > region[2] and coords[1] <region[3]:
            pbc_id.append(id_i)
    return pbc_id

def convert_pbc(atominfo, id_tip, sif):
    natoms = atominfo[0]
    natoms_lst = list(range(0, natoms))
    pbc_coords = []
    # select all atoms
    for i in natoms_lst:
        id_i = atominfo[1][i]
        coords = atominfo[2][i]
        for id_j in id_tip:
            if id_i == id_j:
                pbc_coords.append(coords)
    # move the box ----> (minx,miny) = (0,0)
    # measure the distance along X, Y and Z
    pbc_coords = np.array(pbc_coords)
    minx = min(pbc_coords[:, 0])
    maxx = max(pbc_coords[:, 0])
    miny = min(pbc_coords[:, 1])
    maxy = max(pbc_coords[:, 1])
    minz = min(pbc_coords[:, 2])
    maxz = max(pbc_coords[:, 2])
    xspacing = maxx - minx
    yspacing = maxy - miny
    zspacing = maxz - minz
    rbt_coords = []
    for everyele in pbc_coords:
        rbt_x = everyele[0] - minx
        rbt_y = everyele[1] - miny
        rbt_z = everyele[2] - minz
        rbt = ([rbt_x, rbt_y, rbt_z])
        rbt_coords.append(rbt)
    # copy image and rotate
    # create an symmetry image, rotation symmetry, counterclockwise by pi.
    # can be represented by R=[-1,0;0,-1]
    # [x,y] --> R.*[x,y]=[-x,-y]
    # displace the image rigidly: move up by the distance of (Ymax-Ymin-one_lattice)
    pbc_coords_img = []
    for coords in rbt_coords:
        rot_x = -coords[0] + xspacing + 0.2
        rot_y = -coords[1] + 2*yspacing + (140-sif) * 0.011
        rot_z = coords[2]
        pbc_coords_img.append([rot_x, rot_y, rot_z])
    # Tune the size of boundaries
    xlo = 0
    xhi = xspacing
    ylo = 0
    yhi = 2*yspacing 
    zlo = atominfo[3][4]
    zhi = atominfo[3][5]
    boundaries = ([xlo,xhi,ylo,yhi,zlo,zhi])
    
    # combine both images
    pbc_coords = rbt_coords + pbc_coords_img
    
    return boundaries, pbc_coords

# Write coordinates to Quantum Espresso input
# Can be further customized.
def write_qe(final_data, out_filename, sif):
    # calculate K points according to k_spacing = 0.03
    xlength = final_data[0][1] - final_data[0][0] + 1.6
    ylength = final_data[0][3] - final_data[0][2] + (140-sif) * 0.018
    zlength = final_data[0][5] - final_data[0][4]
    K1 = math.ceil((1 / xlength) / 0.03)
    K2 = math.ceil((1 / ylength) / 0.03)
    K3 = math.ceil((1 / zlength) / 0.03)
    # write QE input file
    qe_inp = open(out_filename,'w')
    qe_inp.write('&CONTROL\n')
    qe_inp.write("    calculation= 'scf'\n")
    qe_inp.write("    verbosity= 'high'\n")
    qe_inp.write("    restart_mode= 'from_scratch'\n")
    qe_inp.write("    prefix= 'Fe_crack_100_010_%i'\n" % sif)
    qe_inp.write("    outdir= './out'\n")
    qe_inp.write("    pseudo_dir= './pseudo'\n")
    qe_inp.write("    tprnfor=.true.\n")
    qe_inp.write("    tstress=.true.\n")
    qe_inp.write("    disk_io='none'\n")
    qe_inp.write("    wf_collect= .false.\n")
    qe_inp.write('/\n\n')

    qe_inp.write('&SYSTEM\n')
    qe_inp.write('    ibrav= 0\n')
    qe_inp.write('    nat= %i\n' % len(final_data[1]))
    qe_inp.write('    ntyp= 1\n')
    qe_inp.write('    ecutwfc= 90\n')
    qe_inp.write('    ecutrho= 1080\n')
    qe_inp.write("    occupations= 'smearing'\n")
    qe_inp.write("    smearing= 'marzari-vanderbilt'\n")
    qe_inp.write('    degauss= 0.01\n')
    qe_inp.write('    nspin= 2\n')
    qe_inp.write('    nosym= .true.\n')
    qe_inp.write('    starting_magnetization(1)= 0.32\n')
    qe_inp.write('/\n\n');

    qe_inp.write('&ELECTRONS\n')
    qe_inp.write('electron_maxstep= 800\n')
    qe_inp.write('conv_thr= 1e-7\n')
    qe_inp.write('mixing_beta= 0.05\n')
    qe_inp.write('/\n\n')

    qe_inp.write('ATOMIC_SPECIES\n')
    qe_inp.write('Fe 55.845 Fe.pbe-spn-rrkjus_psl.0.2.1.UPF\n\n')

    qe_inp.write('K_POINTS automatic\n')
    qe_inp.write('%i %i %i %i %i %i\n\n' % (K1,K2,K3,1,1,1))
    
    qe_inp.write('CELL_PARAMETERS angstrom\n')
    qe_inp.write('%12.8f %12.8f %12.8f \n' % (xlength, 0, 0))
    qe_inp.write('%12.8f %12.8f %12.8f \n' % (0, ylength, 0))
    qe_inp.write('%12.8f %12.8f %12.8f \n\n' % (0, 0, zlength))
    
    qe_inp.write('ATOMIC_POSITIONS angstrom\n')
    for coords in final_data[1]:
        qe_inp.write('%s %12.8f %12.8f %12.8f\n' % ('Fe', coords[0], coords[1], coords[2]))
    print("****************************")
    print("Write Quantum Espresso input file successful!")
    print("QE input file name: %s" % outfile)
    print("****************************")
    
    return

###################################################################
# Execute the code 
###################################################################

#------------------------------------------------------------------
# First part 
# GAP error evaluation
#------------------------------------------------------------------
Maxerror = Nmaxelements(gap_error,1)

# Check if the maximum error is larger than the preset criterion.
if Maxerror > criterion:
    os.system('touch %s' % ('../flag_stop_1'))
    
    
#------------------------------------------------------------------
# Second part 
# Crack tip convertion
#------------------------------------------------------------------
# define the region
region = list([-5, 11, -8, 10])

refdata = readdump(dumpfile_ref)
curdata = readdump(dumpfile_cur)
pbc_id = find_id(refdata, region)
pbc_final = convert_pbc(curdata, pbc_id, sif_k)
write_qe(pbc_final, outfile, sif_k)
