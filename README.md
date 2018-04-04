# convert-glycoct-inp

This work is part of [MatrixDB](http://matrixdb.univ-lyon1.fr) project.

This script in perl converts Glycosaminoglycans (GAGs) in GlycoCT format to INP format.
The output is then used by [POLYS](http://glycan-builder.cermav.cnrs.fr) tool to create 3D structure in PDB format.

This script takes as args a CSV file and a GlycoCT file.
It prints on stdout the corresponding INP format for POLYS tool.

Example of use:
```bash
    perl convertor1.0.pl phi_psi_all_disacc.csv HP_dp02_0001.glycoct > HP_dp02_0001.inp
```

* convertor1.0.pl is the script in perl
* phi_psi_all_disacc.csv is the CSV file in which we have phi and psi values of GAGs disaccharides
* HP_dp02_0001.glycoct is the GlycoCT file
