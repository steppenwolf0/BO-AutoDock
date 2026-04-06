# BO-AutoDock
Evolutionary-Driven Bayesian Optimization for Automated Molecular Docking with AutoDock Vina

This repository implements an automated molecular docking pipeline using Bayesian Optimization and AutoDock Vina.
--------------------------------------------------------------------
Project Structure
In the base folder:

- Receptor:
  The receptor file corresponds to PDB ID 1N26 and is already prepared for AutoDock Vina.
  https://www.rcsb.org/structure/1N26

- Ligands:
  res1 to res20 are the ligands to be tested.
  Ligands were obtained from the ZINC database:
  https://zinc.docking.org/

- AutoDock Vina executables:
  vina.exe
  vina_split.exe

- Perl script:
  runVina.pl

AutoDock Vina executables can be downloaded from:
https://vina.scripps.edu/downloads/
--------------------------------------------------------------------
Requirements

- Windows 10 or Windows 11
- Python
- AutoDock Vina executables (vina.exe and vina_split.exe) located in the base folder

Note:
Installing AutoDock Vina via pip does NOT include vina.exe or vina_split.exe.
You must download the standalone binaries.
--------------------------------------------------------------------
Instructions

1. First, run the following script to verify that everything is working correctly:

   python pipeline0.py

2. Once the pipeline check is successful, run the Bayesian Optimization docking:

   python BO-AutoDock.py
--------------------------------------------------------------------
Notes

- All docking calculations are performed using AutoDock Vina.
- The code has been tested on Windows 10 and Windows 11.
