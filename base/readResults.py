
import os
from statistics import mean


def extract_affinities_from_log(filepath):
    """
    Extracts affinity values from a Vina-like results table.

    Expected table row format (whitespace-separated):
        mode  affinity  rmsd_lb  rmsd_ub

    Example:
        1   -7.6   0.000   0.000

    Returns:
        List[float] of affinities found in the file.
    """
    affinities = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue

            # parts[0] should be integer mode; parts[1] affinity float
            try:
                int(parts[0])           # mode
                affinity = float(parts[1])
                # Optional: validate rmsd columns are floats
                float(parts[2])
                float(parts[3])
            except ValueError:
                continue

            affinities.append(affinity)

    return affinities


def summarize_logs(n=20, folder=".", out_txt="affinity_summary.txt"):
    results = []  # (index, filename, best_affinity or None, status)

    for i in range(1, n + 1):
        filename = f"res{i}.pdbqt_log.log"
        path = os.path.join(folder, filename)

        if not os.path.exists(path):
            results.append((i, filename, None, "MISSING"))
            continue

        affinities = extract_affinities_from_log(path)
        if not affinities:
            results.append((i, filename, None, "NO_AFFINITY_FOUND"))
            continue

        best_aff = min(affinities)  # most negative = best
        results.append((i, filename, best_aff, "OK"))

    valid = [r[2] for r in results if r[2] is not None]

    # Write summary file
    with open(out_txt, "w", encoding="utf-8") as out:
        out.write("Docking affinity summary (best = most negative)\n")
        out.write("=" * 60 + "\n")
        out.write("\n")
        out.write(f"{'Idx':>3} {'File':<25} {'BestAffinity':>12} Status\n")
        out.write("-" * 60 + "\n")

        for idx, fname, best, status in results:
            best_str = f"{best: .3f}" if best is not None else "   N/A"
            out.write(f"{idx:>3} {fname:<25} {best_str:>12} {status}\n")

        out.write("=" * 60 + "\n")
        out.write("Overall statistics (valid files only)\n")
        out.write("-" * 60 + "\n")

        if valid:
            out.write(f"Count valid: {len(valid)} / {n}\n")
            out.write(f"Best (min): {min(valid):.3f}\n")
            out.write(f"Worst (max): {max(valid):.3f}\n")
            out.write(f"Mean: {mean(valid):.3f}\n")
        else:
            out.write("No valid affinities were found.\n")

    # ASCII-only print (safe on Windows cp1252)
    print(f"Wrote summary to: {out_txt}")


if __name__ == "__main__":
    summarize_logs(n=20, folder=".", out_txt="affinity_summary.txt")
