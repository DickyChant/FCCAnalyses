#!/usr/bin/env python3
"""Compare per-event thrust CSVs from `thrust_legacy.py` and `thrust_edm4hep.py`.

Exits 0 if every common event agrees within `--tol` on the gen-level
thrust (default 1e-6). Reco-level is reported but not asserted, because
its expected value depends on how the EDM4hep converter was invoked:
- `delphi_to_edm4hep --require-lvlock-zero` matches legacy's SKELANA
  track-quality cut → reco-thrust is bit-identical (max Δ = 0).
- without that flag the EDM4hep keeps every TracRaw_*, so the two
  writers feed different track multisets to the kernel — the reco delta
  reflects that input-selection difference, not an algorithm bug.
"""
import sys
import csv
import math
import argparse
import statistics


def load(fn):
    with open(fn) as f:
        return {(int(r["runNo"]), int(r["evtNo"])): r for r in csv.DictReader(f)}


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("legacy_csv")
    ap.add_argument("edm4hep_csv")
    ap.add_argument("--tol", type=float, default=1e-6)
    args = ap.parse_args()

    L = load(args.legacy_csv)
    E = load(args.edm4hep_csv)
    common = sorted(set(L) & set(E))
    print(f"common events: {len(common)}  (legacy {len(L)}, edm4hep {len(E)})")

    dr, dg = [], []
    n_eq = n_lt_tol = 0
    n_reco_eq = 0
    for k in common:
        rl, re = L[k], E[k]
        tr_l = float(rl["thrust_reco_charged"])
        tr_e = float(re["thrust_reco_charged"])
        tg_l = float(rl["thrust_gen_all"])
        tg_e = float(re["thrust_gen_all"])
        if not (math.isnan(tr_l) or math.isnan(tr_e)):
            d = abs(tr_l - tr_e)
            dr.append(d)
            if d == 0:
                n_reco_eq += 1
        if not (math.isnan(tg_l) or math.isnan(tg_e)):
            d = abs(tg_l - tg_e)
            dg.append(d)
            if d == 0:
                n_eq += 1
            if d < args.tol:
                n_lt_tol += 1

    print()
    print(f"=== Gen thrust (algorithm equivalence test) ===")
    if dg:
        print(f"  bit-identical (d=0):  {n_eq}/{len(dg)}")
        print(f"  |Δ| < {args.tol:g}:        {n_lt_tol}/{len(dg)}")
        print(f"  mean / median / max:  {statistics.mean(dg):.3e}  {statistics.median(dg):.3e}  {max(dg):.3e}")
    print()
    print(f"=== Reco thrust ===")
    if dr:
        print(f"  bit-identical (d=0):  {n_reco_eq}/{len(dr)}")
        print(f"  mean / median / max:  {statistics.mean(dr):.3e}  {statistics.median(dr):.3e}  {max(dr):.3e}")
        if n_reco_eq != len(dr):
            print(f"  (non-zero Δ → run delphi_to_edm4hep with --require-lvlock-zero")
            print(f"   to apply SKELANA's track-quality cut and match legacy's input set)")

    # Pass/fail on gen
    return 0 if (dg and n_lt_tol == len(dg)) else 1


if __name__ == "__main__":
    sys.exit(main())
