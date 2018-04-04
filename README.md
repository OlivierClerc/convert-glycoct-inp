# convert-glycoct-inp

This script in perl converts Glycosaminoglycans (GAGs) in GlycoCT format to INP format.
The output is then used by POLYS tool to create 3D structure in PDB format.

Takes as args a csv file and a glycoct file.
Print on stdout the corresponding inp format for POLYS tool.

Example of use:
perl convertor1.0.pl phi_psi_all_disacc.csv HP_dp02_0001.glycoct > HP_dp02_0001.inp
