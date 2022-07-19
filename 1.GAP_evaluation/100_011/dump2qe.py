"""
-------------------------------
This script is used to extract the extrapolated atomic environment around crack tip.
It will create a periodic cell that contains the extrapolated atomic environment.
This small periodic cell is used for DFT calculation.
-------------------------------
The reasoning behind the code:
1. Read the crack configuration with GAP predicted per atom error -> find the ID of the 
atom with largest error around crack tip (within a cutoff radius of 40 Angstrom). 
2. Read a reference configuration (undeformed) that has the same atom ID assignment. (This file is dumped
in the same simulation before K loading.)
3. Find the atom ID of a rectangle region that is centered at the atom with the largest error.
3. Using previous obtained IDs to indentify the crack tip in the crack configuration.
4. Create PBC cell for DFT: replicate + rotate + rigid displacement.
5. Write Quantum Espresso input file.
-------------------------------
Some variables related to the crystallography needs to be specified for different crack systems.
-------------------------------
@Author: lei.zhang@rug.nl
@Date: Mon 24 Jan 2022
-------------------------------
"""
import numpy as np
import math
import sys

# Dump/xyz file name: at least two files are required.
# 1. dumpfile_ref -> Reference configuration: undeformed (with the same assignment of IDs).
# 2. xyzfile_cur -> Crack configuration: the one needs to be extracted.
# 3. frame_num -> number of frame in xyzfile (xyz file is a multi-frame data.)
# 4. outfile -> The name for QE input file (with path).
# This is to make sure the small window is correct.
try:
    dumpfile_ref = sys.argv[1]
    xyzfile_cur = sys.argv[2] 
    frame_num = sys.argv[3]
    outfile = sys.argv[4]
except IndexError:
    print("Please enther the dump file of reference, xyz file of current configuration, frame number and output file name!")
    print("****************************")
    print("Example usage: python dump2xyz_pbc.py ref_dump(full path) current_dump(full path) frame number output_QE_name")
    print("****************************")
    exit()

#####################################################
# Define functions that are used to extract crack tips
# and constructing Quantum Espresso input file
#####################################################

# Read data from multi frame XYZ file for a single frame
def readxyz(xyzfile, frameNum):
    data = open(xyzfile,'r')
    N = int(frameNum)
    lines = data.readlines()
    natoms = int(lines[0])
    # the data for Frame number = N is between
    # line number (natoms+2)*N:(natoms+2)*(N+1)
    # excluding the first two lines 
    linesN = lines[(natoms+2)*N+2:(natoms+2)*(N+1)] 
    id_all = []
    position_all = []
    gap_sqrt_variance = []
    for line in linesN:
        line = line.split()
        id_i = int(line[1])
        position = []
        position.append(float(line[2]))
        position.append(float(line[3]))
        position.append(float(line[4]))
        variance = float(line[6])
        id_all.append(id_i)
        position_all.append(position)
        gap_sqrt_variance.append(variance)
            
    return natoms, id_all, position_all, gap_sqrt_variance

# Read LAMMPS dump file
# Total atoms, id, positions, boundaries are returned.
def readdump(dumpfile):
    data = open(dumpfile,'r')
    lines = data.readlines()
    numlines = np.arange(len(lines))
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
    
# find the id of the atom with largest gap error
# around crack tip (within a cutoff of 40 Angstrom 
def find_largest_error_id(natoms, id_all, position_all, gap_sqrt_variance):
    cutoff = 20
    largest_variance = 0.0000000001
    for atom in np.arange(natoms):
        atomid = id_all[atom]
        pos_x = position_all[atom][0]
        pos_y = position_all[atom][1]
        sqrt_variance = gap_sqrt_variance[atom] 
        if pos_x**2 + pos_y**2 < cutoff**2: 
            if sqrt_variance > largest_variance:
                largest_variance = sqrt_variance
                id_largest = atomid
                posx_largest = pos_x  
                posy_largest = pos_y
    
    return id_largest

# find the tip region to be extracted
# in the reference frame
def find_tip_id(largesterrorid, natoms, id_all, position_all, region_ylo, region_yhi, xlatspacing, ylatspacing):
    ## Find the region first
    # find the coordinates of the largesterrorid
    for i in np.arange(natoms):
        id_i = id_all[i]
        if id_i == largesterrorid:
            posx_largest = position_all[i][0]
            posy_largest = position_all[i][1]
    # Along X axis, the position is determined by the atom with 
    # largest error
    # left part: + 3 layers; right part: + 4 layers
    # atomic distance along x <110> direction is
    xlo = posx_largest - 3 * xlatspacing - 0.5 # Minus/Plus an error tolenrence of 1 due to distortion
    xhi = posx_largest + 4 * xlatspacing + 0.5
    # Along Y axis, it is defined by a preview of the structure.
    ylo = region_ylo 
    yhi = region_yhi
    region = [xlo, xhi, ylo, yhi]    
# Store the id in the extracted crack tip region
    pbc_id = []
    for i in np.arange(natoms):
        id_i = id_all[i]
        coords = position_all[i]        
        # select atoms in side a certain region.
        # store the atomic ID of all atoms.
        if coords[0] > region[0] and coords[0] < region[1] and \
        coords[1] > region[2] and coords[1] <region[3]:
            pbc_id.append(id_i)
            
    return pbc_id

# Convert the small region to a PBC cell.
def convert_pbc(atominfo, boundaries, id_tip, xlatspacing, ylatspacing):
    newid = 0
    origin_coords = []
    natoms = atominfo[0]
    id_all = atominfo[1]
    position_all = atominfo[2]
    for i in np.arange(natoms):
        id_i = id_all[i]
        coords = position_all[i]
        for id_j in id_tip:
            if id_i == id_j:
                newid = newid + 1
                newcoords = (newid,coords[0],coords[1],coords[2])
                origin_coords.append(newcoords)
    # measure the distance along X, Y and Z
    origin_coords = np.array(origin_coords)
    rbt_coords = []
    # move the lower left atom to origin (0,0)
    for everyele in origin_coords:
        rbt_x = everyele[1] - min(origin_coords[:, 1])
        rbt_y = everyele[2] - min(origin_coords[:, 2])
        rbt_z = everyele[3]
        rbt = ([everyele[0], rbt_x, rbt_y, rbt_z])
        rbt_coords.append(rbt)
    rbt_coords = np.array(rbt_coords)
    # find coordinates of left upper and bottom atom 
    upl = rbt_coords[np.argmax(rbt_coords[:,2])]
    lol = rbt_coords[np.argmin(rbt_coords[:,2])]
    upl_index = np.where(rbt_coords[:,2] == upl[2])
    lol_index = np.where(rbt_coords[:,2] == lol[2])
    upl_coords = rbt_coords[upl_index][0]
    lol_coords = rbt_coords[lol_index][0]
    # find coordinates of right upper and bottom atom 
    # find the eight atoms in the right of the crack tip
    sorted_x = np.argsort(rbt_coords[:,1])
    max8x = rbt_coords[sorted_x[-8:]]  
    # right upper and bottom atom 
    upr_coords = np.array([0,0,-10000,0])
    lor_coords = np.array([0,0,10000,0])
    for ele in max8x:
        if ele[2] > upr_coords[2]:
            upr_coords = ele
        if ele[2] < lor_coords[2]:
            lor_coords = ele
    # copy image and rotate
    # create an symmetry image, rotation symmetry, counterclockwise by pi.
    # can be represented by R=[-1,0;0,-1]
    # [x,y] --> R.*[x,y]=[-x,-y]
    # displace the image rigidly: move up by the distance of (Ymax-Ymin-one_lattice)
    
    pbc_coords_img = []
    for coords in rbt_coords:
        atomid = int(coords[0])
        rot_x = -coords[1] 
        rot_y = -coords[2] 
        rot_z = coords[3]
        pbc_coords_img.append([len(rbt_coords) + atomid,rot_x, rot_y, rot_z])
    # the left bottom corner atom of rotated image is
    # the right upper of origin before rotation
    # find the coordinate of that atom
        if atomid == int(upr_coords[0]):
            upr_coords_image = [len(rbt_coords) + atomid,rot_x, rot_y, rot_z]
        if atomid == int(lor_coords[0]):
            lor_coords_image = [len(rbt_coords) + atomid,rot_x, rot_y, rot_z]
    # find the rigid displacement of the rotated image
    # in order to align with the original one
    # displacement along X
    disp_x = lol_coords[1] - lor_coords_image[1]
    disp_y = upl_coords[2] - upr_coords_image[2] + ylatspacing
    # align image -> move the lower left atom to origin (0,0)
    all_coords = rbt_coords.tolist()
    for everyele in pbc_coords_img:
        rbt_x = everyele[1] + disp_x
        rbt_y = everyele[2] + disp_y
        rbt_z = everyele[3]
        rbt = ([everyele[0], rbt_x, rbt_y, rbt_z])
        all_coords.append(rbt)
    all_coords = np.array(all_coords)
    pbc_coords = all_coords[:,1:4]
    # add Gaussin noisy (0,0.05)
    noise = np.random.normal(0,0.05,[len(pbc_coords),3])
    pbc_coords = (pbc_coords + noise).tolist()
    # Tune the size of boundaries
    xlo = 0
    xhi = np.min(max8x[:,1]) + xlatspacing
    ylo = 0
    # find the coordinates of upper left image
    up_lim = all_coords[np.where(all_coords[:,0] == lor_coords_image[0])]
    yhi = up_lim[0][2] + ylatspacing
    zlo = boundaries[4]
    zhi = boundaries[5]
    boundaries = ([xlo,xhi,ylo,yhi,zlo,zhi])
    
    return boundaries, pbc_coords

# Write structure information to Quantum Espresso input file
def write_qe(final_data, out_filename):
    
    # calculate K points according to k_spacing=0.03
    xlength = final_data[0][1] - final_data[0][0]
    ylength = final_data[0][3] - final_data[0][2]
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
    qe_inp.write("    prefix= 'Fe_crack_110_010'\n")
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
# define some variables related to the crystallography
latconst = 2.834
xlatspacing = latconst * 0.5 * np.sqrt(2)
ylatspacing = latconst * 0.5 * np.sqrt(1)
region_ylo = -5
region_yhi = 6
#----------------------
# Read current configuration
cur_config = readxyz(xyzfile_cur, frame_num)
# find_largest_error_id
largesterrorid = find_largest_error_id(cur_config[0], cur_config[1], cur_config[2], cur_config[3])
# Read reference configuration
ref_config = readdump(dumpfile_ref)
# Find id of atoms around crack tip in a square shape
tip_ids = find_tip_id(largesterrorid, ref_config[0], ref_config[1], ref_config[2], region_ylo, region_yhi, xlatspacing, ylatspacing)
boundaries = ref_config[3]
# Convert crack tip to a PBC cell
pbc_final = convert_pbc(cur_config, boundaries, tip_ids, xlatspacing, ylatspacing)
# Write input file for QE 
write_qe(pbc_final, outfile)
