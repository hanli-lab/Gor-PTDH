This is the Rosetta docking protocol used in the manuscript from Professor Han Li lab at UC Irvine and Professor Justin Siegel lab at UC Davis. The Rosetta files listed here are for docking the ligand NMN and PS-PTDH LY-7 variant. Rosetta version: 2018.24.post.dev+17.master.450949e 450949e481542459ae6534e867a61fe9709846be

In order to successfully run the simulation, all the files are required:
1.  PDB file with the target ligand NMN bound in the active site. (PtDH_WT.pdb)
2.  NMN parameter file that is readable by Rosetta and the corresponding conformers' library pdb. (X00.params, X00_conformers.pdb)
3.  Constraint file that restricts the movement of active site catalytic residues and relative position between ligand and active site residues (ptdh_nmn.cst)
4.  Rosetta simulation option flag file (dock_flags)
5.  Rosetta ligand docking protocol xml file (dock.xml)
6.  Slurm submission file, only executed for cluster jobs running (submit_cluster_LY7.bash)

Run the following command to perform Rosetta docking locally in MacOS system:
Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -database /Rosetta/main/database @dock_flags -out:pdb_gz -out:path:all ./

Run the following command to perform Rosetta docking in cluster:
Sbatch submit_cluster_LY7.bash
