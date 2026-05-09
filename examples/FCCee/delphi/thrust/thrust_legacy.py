#!/usr/bin/env python3
"""Reference: per-event thrust on the legacy delphi-nanoaod TTrees.

Reads `t / tgenBefore` produced by jingyu's `delphi-nanoaod` (skelana-based
binary, takes the same DELPHI SDST as `delphi-raw-nanoaod`) and writes the
same CSV format as `thrust_edm4hep.py`. Exists only as the reference for
the algorithm-equivalence cross-check; for production use the FCCAnalyses
port `thrust_edm4hep.py` directly.
"""
import sys
import argparse
import numpy as np
import ROOT

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from common_functions import thrust_axis_fast_optimized

LONG_LIVED_NEUTRALS_ABS = {130, 3122, 310}


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("input")
    ap.add_argument("output_csv")
    args = ap.parse_args()

    f = ROOT.TFile.Open(args.input)
    t_reco, t_gen = f.Get("t"), f.Get("tgenBefore")
    n = t_reco.GetEntries()
    print(f"legacy nanoaod: {n} events", flush=True)

    with open(args.output_csv, "w") as out:
        out.write("runNo,evtNo,nReco_charged,nGen,thrust_reco_charged,thrust_gen_all\n")
        for i in range(n):
            t_reco.GetEntry(i)
            t_gen.GetEntry(i)
            run, evn = int(t_reco.RunNo), int(t_reco.EventNo)

            # Reco charged
            n_p = int(t_reco.nParticle)
            n_reco = 0
            thrust_reco = float("nan")
            if n_p > 0:
                q = np.fromiter((t_reco.charge[k] for k in range(n_p)), dtype=np.float32)
                sel = np.abs(q) > 0.1
                if sel.sum() >= 2:
                    px = np.fromiter((t_reco.px[k] for k in range(n_p)), dtype=np.float32)
                    py = np.fromiter((t_reco.py[k] for k in range(n_p)), dtype=np.float32)
                    pz = np.fromiter((t_reco.pz[k] for k in range(n_p)), dtype=np.float32)
                    p3 = np.stack([px[sel], py[sel], pz[sel]], axis=1)
                    _, T = thrust_axis_fast_optimized(p3, include_met=False)
                    thrust_reco = float(1.0 - T)
                n_reco = int(sel.sum())

            # Gen — drop pid==0 (string remnants) + long-lived neutrals
            n_g = int(t_gen.nParticle)
            n_gen = 0
            thrust_gen = float("nan")
            if n_g > 0:
                pid = np.fromiter((t_gen.pid[k] for k in range(n_g)), dtype=np.int32)
                keep = (pid != 0) & ~np.isin(np.abs(pid), list(LONG_LIVED_NEUTRALS_ABS))
                if keep.sum() >= 2:
                    px = np.fromiter((t_gen.px[k] for k in range(n_g)), dtype=np.float32)
                    py = np.fromiter((t_gen.py[k] for k in range(n_g)), dtype=np.float32)
                    pz = np.fromiter((t_gen.pz[k] for k in range(n_g)), dtype=np.float32)
                    p3 = np.stack([px[keep], py[keep], pz[keep]], axis=1)
                    _, Tg = thrust_axis_fast_optimized(p3, include_met=False)
                    thrust_gen = float(1.0 - Tg)
                n_gen = int(keep.sum())

            out.write(f"{run},{evn},{n_reco},{n_gen},{thrust_reco:.10g},{thrust_gen:.10g}\n")

    print(f"wrote {args.output_csv}", flush=True)


if __name__ == "__main__":
    main()
