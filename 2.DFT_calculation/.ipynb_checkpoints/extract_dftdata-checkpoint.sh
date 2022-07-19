#/bin/bash
# Extract the DFT data (total energy, cell information, atomic positions and forces, k-points, PBC)
# Output it in extendXYZ format

# Define QE output file.
outfile=$1
iter=$2

# Define the configuration name and output file
config_name='small_crack'
output='DFT_crack.'${iter}

# Define constants
Ry_to_eV='13.60569301'
Bohr_to_Ang='0.529177208'
RyBohr_to_eVAngs='25.71104309541616'
RyBohr3_to_eVAngs3=`echo "scale=10;  $Ry_to_eV /( $Bohr_to_Ang*$Bohr_to_Ang*$Bohr_to_Ang)" | bc`

# Number of atoms in the configuration
alat=`echo "scale=10;  $Bohr_to_Ang * $(grep 'lattice parameter (alat)' ${outfile} | awk '{print $5}')" | bc`
Natom=$(grep 'number of atoms/cell' ${outfile} | tr -dc '0-9')
# Cutoff energy
KineticCut=$(grep 'kinetic-energy cutoff' ${outfile} | awk '{print $4}')
RhoCut=$(grep 'charge density cutoff' ${outfile} | awk '{print $5}')
# Box size
a1=`echo "scale=10;  $alat * $(grep 'a(1) =' ${outfile} | awk '{print $4}')" | bc`
a2=`echo "scale=10;  $alat * $(grep 'a(1) =' ${outfile} | awk '{print $5}')" | bc`
a3=`echo "scale=10;  $alat * $(grep 'a(1) =' ${outfile} | awk '{print $6}')" | bc`
b1=`echo "scale=10;  $alat * $(grep 'a(2) =' ${outfile} | awk '{print $4}')" | bc`
b2=`echo "scale=10;  $alat * $(grep 'a(2) =' ${outfile} | awk '{print $5}')" | bc`
b3=`echo "scale=10;  $alat * $(grep 'a(2) =' ${outfile} | awk '{print $6}')" | bc`
c1=`echo "scale=10;  $alat * $(grep 'a(3) =' ${outfile} | awk '{print $4}')" | bc`
c2=`echo "scale=10;  $alat * $(grep 'a(3) =' ${outfile} | awk '{print $5}')" | bc`
c3=`echo "scale=10;  $alat * $(grep 'a(3) =' ${outfile} | awk '{print $6}')" | bc`
echo "======="
echo $b1
# atomic positions
# get the number of the line before the line containing atomic positions
line_pos=$(grep -n 'site n.     atom                  positions (alat units)' ${outfile} | awk -F  ":" '{print $1}')
line_for=$(grep -n 'Forces acting on atoms' ${outfile} | awk -F  ":" '{print $1}')
line_for=$(echo "$(($line_for + 1))")

# Loop to get atomic positions and forces
for k in $(seq 1 1 $Natom)
do
line_1=$(echo "$(($line_pos + $k))")
# Read Atomic coordinates
x1[$k]=`echo "scale=10;  $alat * $(awk 'NR=='$line_1' {print $7}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
x2[$k]=`echo "scale=10;  $alat * $(awk 'NR=='$line_1' {print $8}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
x3[$k]=`echo "scale=10;  $alat * $(awk 'NR=='$line_1' {print $9}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`

line_2=$(echo "$(($line_for + $k))")
# Read the forces on atoms
f1[$k]=`echo "scale=10;  $RyBohr_to_eVAngs * $(awk 'NR=='$line_2' {print $7}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
f2[$k]=`echo "scale=10;  $RyBohr_to_eVAngs * $(awk 'NR=='$line_2' {print $8}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
f3[$k]=`echo "scale=10;  $RyBohr_to_eVAngs * $(awk 'NR=='$line_2' {print $9}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
#f3[$k]=$(awk 'NR=='$line_2' {print $9}' ${outfile})
done
# get the total energy
total_E=`echo "scale=8;  $Ry_to_eV * $(grep '!' ${outfile} | awk '{print $5}')" | bc`
echo $total_E
# stress on box
line_str=$(grep -n 'total   stress  (Ry/bohr\*\*3)' ${outfile} | awk -F  ":" '{print $1}')
line_str1=$(echo "$(($line_str + 1))")
line_str2=$(echo "$(($line_str + 2))")
line_str3=$(echo "$(($line_str + 3))")
sigma_xx=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str1' {print $1}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_xy=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str1' {print $2}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_xz=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str1' {print $3}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_yx=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str2' {print $1}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_yy=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str2' {print $2}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_yz=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str2' {print $3}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_zx=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str3' {print $1}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_zy=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str3' {print $2}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`
sigma_zz=`echo "scale=10; $RyBohr3_to_eVAngs3 * $(awk 'NR=='$line_str3' {print $3}' ${outfile})" | bc | awk '{printf "%.10f\n", $0}'`

echo $Natom |tee -a ${output}
echo 'Lattice=''"'$a1 $a2 $a3 $b1 $b2 $b3 $c1 $c2 $c3'"' 'Properties=species:S:1:pos:R:3:force:R:3:Z:I:1' 'config_type='$config_name 'energy='$total_E 'pbc="T T T"' |tee -a ${output}
for l in $(seq 1 1 $Natom)
do
echo 'Fe' ${x1[$l]} ${x2[$l]} ${x3[$l]} ${f1[$l]} ${f2[$l]} ${f3[$l]} '26'|tee -a ${output}
done
