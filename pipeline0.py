
import os
import sys
import shutil
import subprocess
import re
from pathlib import Path


def run_command(cmd, cwd=None):
    """
    Runs a command, captures stdout/stderr, raises with context if it fails.
    Forces UTF-8 for Python subprocess output (safe on Windows).
    """
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUTF8"] = "1"

    print("\n>>> Running: " + " ".join(cmd))
    if cwd:
        print("    (cwd: " + str(cwd) + ")")

    cp = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        check=False
    )

    if cp.stdout:
        print("---- stdout ----")
        print(cp.stdout)
    if cp.stderr:
        print("---- stderr ----")
        print(cp.stderr)

    if cp.returncode != 0:
        raise RuntimeError(
            "Command failed with exit code {code}: {cmd}\n"
            "cwd={cwd}\n"
            "stderr:\n{stderr}".format(
                code=cp.returncode,
                cmd=" ".join(cmd),
                cwd=cwd,
                stderr=cp.stderr
            )
        )
    return cp


def copy_base_to_run(base_dir: str, run_dir: str) -> Path:
    """Copies base_dir into run_dir. Deletes run_dir first if it exists."""
    base = Path(base_dir)
    run = Path(run_dir)

    if not base.exists() or not base.is_dir():
        raise FileNotFoundError("Base folder not found: " + str(base.resolve()))

    if run.exists():
        print("Removing existing folder: " + str(run.resolve()))
        shutil.rmtree(run)

    print("Copying {src} -> {dst}".format(src=base.resolve(), dst=run.resolve()))
    shutil.copytree(base, run)
    return run


def write_conf_vs(run_path: Path,
                  receptor="receptor.pdbqt",
                  center_x=22.208, center_y=48.77, center_z=76.016,
                  size_x=40, size_y=40, size_z=40,
                  num_modes=10, energy_range=4, exhaustiveness=16,
                  filename="conf_vs.txt") -> Path:
    """Writes an AutoDock Vina-style config file into run_path/filename."""
    conf_path = run_path / filename
    content = (
        f"receptor = {receptor}\n\n"
        f"center_x = {center_x}\n"
        f"center_y = {center_y}\n"
        f"center_z = {center_z}\n\n"
        f"size_x = {size_x}\n"
        f"size_y = {size_y}\n"
        f"size_z = {size_z}\n\n"
        f"num_modes = {num_modes}\n"
        f"energy_range = {energy_range}\n"
        f"exhaustiveness = {exhaustiveness}\n"
    )
    conf_path.write_text(content, encoding="utf-8")
    print("Generated config: " + str(conf_path.resolve()))
    return conf_path


def extract_mean_from_affinity_summary(summary_path: Path) -> float:
    """
    Reads affinity_summary.txt and extracts the numeric value from a line like:
        Mean: -6.340
    """
    if not summary_path.exists():
        raise FileNotFoundError("Summary file not found: " + str(summary_path.resolve()))

    lines = summary_path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        m = re.search(r"^\s*Mean:\s*([-+]?\d+(?:\.\d+)?)\s*$", line)
        if m:
            return float(m.group(1))

    tail = "\n".join(lines[-25:])
    raise ValueError(
        "Could not find 'Mean:' line in {f}.\nLast lines:\n{tail}".format(
            f=summary_path.name,
            tail=tail
        )
    )


def mainPipeline(base_dir="base",
         run_id=0,
         run_prefix="run",
         vina_script="runVina.py",
         read_results_script="readResults.py",
         conf_name="conf_vs.txt",
         summary_name="affinity_summary.txt",
         mean_outfile="mean.txt",
         receptor="receptor.pdbqt",
         center_x=22.208, center_y=48.77, center_z=76.016,
         size_x=40, size_y=40, size_z=40,
         num_modes=10,
         energy_range=4,
         exhaustiveness=16):

    run_dir = f"{run_prefix}{run_id}"
    run_path = copy_base_to_run(base_dir, run_dir)

    write_conf_vs(
        run_path,
        receptor=receptor,
        center_x=center_x, center_y=center_y, center_z=center_z,
        size_x=size_x, size_y=size_y, size_z=size_z,
        num_modes=num_modes,
        energy_range=energy_range,
        exhaustiveness=exhaustiveness,
        filename=conf_name
    )

    # Optional sanity check
    receptor_path = run_path / receptor
    if not receptor_path.exists():
        print("Warning: receptor file not found: " + str(receptor_path.resolve()))

    # Run vina script
    vina_path = run_path / vina_script
    if not vina_path.exists():
        raise FileNotFoundError("Vina script not found: " + str(vina_path.resolve()))
    run_command([sys.executable, vina_script], cwd=str(run_path))

    # Run readResults.py (creates affinity_summary.txt)
    rr_path = run_path / read_results_script
    if not rr_path.exists():
        raise FileNotFoundError("readResults.py not found: " + str(rr_path.resolve()))
    run_command([sys.executable, read_results_script], cwd=str(run_path))

    # Parse mean from affinity_summary.txt
    summary_path = run_path / summary_name
    mean_val = extract_mean_from_affinity_summary(summary_path)

    print("\nMean best affinity: " + str(mean_val))

    # Write mean to file
    out_path = run_path / mean_outfile
    out_path.write_text(str(mean_val) + "\n", encoding="utf-8")
    print("Wrote mean to: " + str(out_path.resolve()))

    # IMPORTANT: return mean so askPoint.py can use it
    return mean_val


if __name__ == "__main__":
    run_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    base_dir = sys.argv[2] if len(sys.argv) > 2 else "base"
    mainPipeline(base_dir=base_dir, run_id=run_id)
