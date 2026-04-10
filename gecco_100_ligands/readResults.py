
import os
import re
from statistics import mean

# Regex to match a Vina-like result row:
# e.g. "   1   -7.6   0.000   0.000"
ROW_RE = re.compile(
    r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$"
)

def extract_affinities_from_log(filepath):
    """
    Returns a list of affinity floats found in Vina-like tables.
    """
    affinities = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            m = ROW_RE.match(line)
            if m:
                # mode = int(m.group(1))  # not used, but available
                affinity = float(m.group(2))
                affinities.append(affinity)
    return affinities

def summarize_logs(n=100, folder=".", out_txt="affinity_summary.txt"):
    results = []  # list of tuples: (index, filename, best_affinity or None)

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

        best_aff = min(affinities)  # most negative is "best" in Vina convention
        results.append((i, filename, best_aff, "OK"))

    # Collect valid affinities for global stats
    valid = [r[2] for r in results if r[2] is not None]

    # Write summary
    with open(out_txt, "w", encoding="utf-8") as out:
        out.write("Docking affinity summary (best = most negative)\n")
        out.write("=" * 60 + "\n\n")
        out.write(f"{'Idx':>3}  {'File':<25}  {'BestAffinity':>12}  Status\n")
        out.write("-" * 60 + "\n")

        for i, fname, best, status in results:
            best_str = f"{best: .3f}" if best is not None else "   N/A"
            out.write(f"{i:>3}  {fname:<25}  {best_str:>12}  {status}\n")

        out.write("\n" + "=" * 60 + "\n")
        out.write("Overall statistics (valid files only)\n")
        out.write("-" * 60 + "\n")

        if valid:
            out.write(f"Count valid: {len(valid)} / {n}\n")
            out.write(f"Best (min):  {min(valid):.3f}\n")
            out.write(f"Worst (max): {max(valid):.3f}\n")
            out.write(f"Mean:        {mean(valid):.3f}\n")
        else:
            out.write("No valid affinities were found.\n")

    print(f"✅ Wrote summary to: {out_txt}")

if __name__ == "__main__":
    summarize_logs(n=100, folder=".", out_txt="affinity_summary.txt")
