#!/usr/bin/env python3
"""One command that reruns everything, from raw downloads to figures.

  python run_all.py            full study (hours; safe to stop and restart)
  python run_all.py --pilot    small real-data run to check your setup (about 10 minutes)

Stages: setup -> candidates -> scan -> calibration on the original's articles -> validation sheet -> analyze.
The manual validation labels (step 3) are the one thing you do by hand, see README.md.
"""
import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

# small PMC id slices used by --pilot (one per year), so no full-year search is needed
PILOT_SLICES = {2019: (6400000, 6430000), 2020: (7000000, 7030000), 2021: (8200000, 8230000),
                2022: (9000000, 9030000), 2023: (10100000, 10130000),
                2024: (10900000, 10930000), 2025: (11800000, 11830000)}


def run(*args):
    print("\n$", " ".join(str(a) for a in args), flush=True)
    r = subprocess.run([PY, *map(str, args)], cwd=HERE)
    if r.returncode:
        sys.exit("stopped: %s failed (exit %d)" % (args[0], r.returncode))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", action="store_true", help="small run, writes to data/pilot.sqlite and output_pilot/")
    ap.add_argument("--years", type=int, nargs="+")
    ap.add_argument("--target", type=int, help="eligible articles per year (default from config.py)")
    ap.add_argument("--procs", type=int, help="worker processes")
    ap.add_argument("--with-original-lists", action="store_true",
                    help="also rescan with the 2015 gene lists only (sensitivity run, doubles the time)")
    ap.add_argument("--skip", nargs="*", default=[], choices=["setup", "candidates", "scan", "calibrate", "analyze", "validation"])
    args = ap.parse_args()

    import config
    years = args.years or config.YEARS
    extra = []
    if args.procs:
        extra += ["--procs", args.procs]

    if "setup" not in args.skip:
        run("00_setup_reference.py")

    if args.pilot:
        target = args.target or 60
        for y in years:
            lo, hi = PILOT_SLICES[y]
            run("01_build_candidates.py", "--years", y, "--uid-range", lo, hi, "--tag", "_pilot")
        if "scan" not in args.skip:
            run("02_scan_articles.py", "--years", *years, "--tag", "_pilot", "--target", target,
                "--max-candidates", 500, "--db", "data/pilot.sqlite", *extra)
        if "analyze" not in args.skip:
            run("04_analyze.py", "--db", "data/pilot.sqlite", "--target", target, "--out", "output_pilot")
        print("\nPilot done. Numbers from the pilot are NOT results (tiny, non-random slices).")
        return

    target = args.target or config.TARGET_ELIGIBLE_PER_YEAR
    if "candidates" not in args.skip:
        run("01_build_candidates.py", "--years", *years)
    if "scan" not in args.skip:
        run("02_scan_articles.py", "--years", *years, "--target", target, *extra)
        if args.with_original_lists:
            run("02_scan_articles.py", "--years", *years, "--target", target, *extra,
                "--original-lists-only", "--db", "data/scan_original_lists.sqlite")
    if "calibrate" not in args.skip:
        run("05_calibrate_original.py")
    if "validation" not in args.skip:
        run("03_validation.py", "make", "--target", target)
    if "analyze" not in args.skip:
        run("04_analyze.py", "--target", target)
    print("\nDone. Next: label the validation sheet (README step 3), then rerun: "
          "python run_all.py --skip setup candidates scan calibrate validation")


if __name__ == "__main__":
    main()
