import subprocess
import sys
from pathlib import Path

LIGAND_FILE = "ligands.txt"
CONFIG_FILE = "conf_vs.txt"
VINA_CMD = "vina.exe"

if not Path(LIGAND_FILE).is_file():
    print(f"Cannot open file: {LIGAND_FILE}")
    sys.exit(1)

with open(LIGAND_FILE, "r") as f:
    ligands = [line.strip() for line in f.readlines() if line.strip()]
    for ligand in ligands:
        print(ligand, flush=True)

    for ligand in ligands:
        if ligand:
            print(ligand, flush=True)
            subprocess.run([
                VINA_CMD,
                "--config", CONFIG_FILE,
                "--ligand", ligand,
                "--log", f"{ligand}_log.log"
            ])
