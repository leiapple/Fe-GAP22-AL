# The script includes all necessary inputs/codes for GAP active learning.

The active learning steps are listed below:

0. Crack simulation: K test (LAMMPS)
1. GAP evaluation (QUIP) + Crack tip extraction (PYTHON)
2. DFT calculation (Quantum Espresso)
3. GAP training (QUIP) -> 0.K test


The other folders are:

- GAP_model: current GAP that is used for current atomistic simulation
- DATA: data produced in each iteration

Prerequisites are needed:

- Python module ovito/matplotlib are needed.
- The path for LAMMPS/QE executable needs to be customised by the user in slurm submit file.
- The path for QUIP need to be defined in *1.GAP_evaluation/100_011/local_variances.sh* ,  *1.GAP_evaluation/110_001/local_variances.sh* and *3.NEWGAP_training/submit*.
- 