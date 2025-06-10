# Molecular Dynamics Simulation Analysis Script
- Author: iawnix

# Usage
1. bin/aaname:
    - bash: `aaname N`
    - -> ASN
    - bash: `aaname ASN`
    - -> N
2. bin/AlterBfactorPDB
    - bash: `AlterBfactorPDB pro.pdb rmsf.txt pro_out`
    - -> This script can add the rmsf values about residues in rmsf.txt as bfactors to pro.pdb and generate a new file pro_out.pdb, Then you con run `spectrum b, blue_white_red, minimum=min_value, maximum=max_value` in Pymol
3. bin/MDplot
    - This script is used to draw a simple analysis chart
4. bin/MolVolume.py
    - bash: `MolVolume.py -cpu xx -top xxx.prmtop -pdb xxx.pdb -box "x,y,z,x1,y1,z1" -space xx > xx.log`
    - This script is used to estimate the volume of molecules using Monte Carlo method based on the atomic vdw radius in the top file
5. bin/RESP.py
    - This script is used to assist in the construction of the Amber position for small molecules, and the charge is calculated using Gaussian's RESP charge
