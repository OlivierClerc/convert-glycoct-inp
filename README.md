# convert-glycoct-inp

This work is part of [MatrixDB](http://matrixdb.univ-lyon1.fr) project.

This script in perl converts Glycosaminoglycans (GAGs) in GlycoCT format to INP format.
The output is then used by [POLYS](http://glycan-builder.cermav.cnrs.fr) tool to create 3D structure in PDB format.

This script takes as args a CSV file and a GlycoCT file.
It prints on stdout the corresponding INP format for POLYS tool.

## Important to know
* This script doesn't support ALT, UND and ISO sections of GlycoCT
* This script support REP sections but not nested repetitions
* This script only converts GAG's monosaccharides
* In order to build a 3D structure, POLYS need to have all information about monosaccharides: alpha/beta, L/D, etc...
* Iduronic acid can have many conformations (1C4, 4C1, 2S0). We have not this information in GlycoCT format but we need it for POLYS. Thus we decided to choose 1C4 conformation by default.

Example of use:
```bash
    perl convertor2.0.pl phi_psi_all_disacc.csv HP_dp02_0001.glycoct > HP_dp02_0001.inp
```

* convertor2.0.pl is the script in perl
* phi_psi_all_disacc.csv is the CSV file in which we have phi and psi values of GAG's disaccharides
* HP_dp02_0001.glycoct is a GlycoCT file of disaccharide of heparin
