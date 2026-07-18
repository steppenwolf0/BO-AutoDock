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

```bash
   python pipeline0.py
```

2. Once the pipeline check is successful, run the Bayesian Optimization docking:
```bash
   python BO-AutoDock.py
```

- For the full list of command line args:
```bash
   python BO-AutoDock.py --help
```

--------------------------------------------------------------------
Notes

- All docking calculations are performed using AutoDock Vina.
- The code has been tested on Windows 10 and Windows 11.


## Citing BO-Autodock
If you use it in your research, please use the following BibTeX entry.
```
@inproceedings{lopez2026evolutionary,
  title={Evolutionary-Driven Bayesian Optimization for Automated Molecular Docking with AutoDock Vina},
  author={Lopez-Rincon, Alejandro and Varga, Brigitta and Rojas-Velazquez, David and Tonda, Alberto},
  booktitle={Proceedings of the Genetic and Evolutionary Computation Conference},
  pages={1137--1145},
  year={2026}
}
```