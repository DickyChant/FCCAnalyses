#!/usr/bin/env python3
"""Per-event thrust on a DELPHI EDM4hep file (FCCAnalyses port of
jingyu's delphi-nanoaod analysis_thrust.py).

Same numpy thrust kernel as the original (see common_functions.py copied
from upstream); only the input layer is replaced — instead of reading
the legacy delphi-nanoaod TTrees `t / tgenBefore`, we read EDM4hep
`ReconstructedParticles` and `Particle` collections produced by
`delphi-improved-reco/ild/delphi_to_edm4hep`.

Validated against the legacy analyzer on a 100-event Z->bb Pythia run:
gen-level thrust matches at 1e-7 (float precision) on 100/100 events
after harmonizing the long-lived-neutral convention (PSHLUJ filters
K_L / Lambda / K_S as "decays in flight"; legacy tgenBefore keeps
them — both analyzers drop the same set here).
"""
import sys
import argparse
import numpy as np
import uproot

# common_functions lives next to this script
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent))
from common_functions import thrust_axis_fast_optimized

# Long-lived neutrals that PSHLUJ filters from the LUJETS record but
# the legacy SKELANA-based writer keeps. Drop in both analyzers for
# apples-to-apples comparison.
LONG_LIVED_NEUTRALS_ABS = {130, 3122, 310}


def get_run_evt(gp_keys_event, gp_vals_event):
    keys = list(gp_keys_event)
    vals = list(gp_vals_event)
    out = {}
    for k, v in zip(keys, vals):
        out[k] = int(v[0])
    return out.get("runNumber", -1), out.get("eventNumber", -1)


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("input")
    ap.add_argument("output_csv")
    args = ap.parse_args()

    ev = uproot.open(args.input)["events"]
    n = ev.num_entries
    print(f"EDM4hep: {n} events", flush=True)

    # Reco
    rp_px = ev["ReconstructedParticles/ReconstructedParticles.momentum.x"].array()
    rp_py = ev["ReconstructedParticles/ReconstructedParticles.momentum.y"].array()
    rp_pz = ev["ReconstructedParticles/ReconstructedParticles.momentum.z"].array()
    rp_q = ev["ReconstructedParticles/ReconstructedParticles.charge"].array()
    # Gen
    mc_st = ev["Particle/Particle.generatorStatus"].array()
    mc_pdg = ev["Particle/Particle.PDG"].array()
    mc_px = ev["Particle/Particle.momentum.x"].array()
    mc_py = ev["Particle/Particle.momentum.y"].array()
    mc_pz = ev["Particle/Particle.momentum.z"].array()
    # Run / event keys live in podio's GenericParameter int map
    gp_keys = ev["GPIntKeys"].array()
    gp_vals = ev["GPIntValues"].array()

    with open(args.output_csv, "w") as out:
        out.write(
            "runNo,evtNo,nReco_charged,nGen,thrust_reco_charged,thrust_gen_all\n"
        )
        for i in range(n):
            run, evn = get_run_evt(gp_keys[i], gp_vals[i])
            if evn == 0:
                # DELSIM run-start / housekeeping records
                continue

            # --- Reco (charged-only) -----------------------------------
            q_i = np.asarray(rp_q[i])
            sel = np.abs(q_i) > 0.1
            n_reco = int(sel.sum())
            if n_reco >= 2:
                p3 = np.stack(
                    [
                        np.asarray(rp_px[i])[sel],
                        np.asarray(rp_py[i])[sel],
                        np.asarray(rp_pz[i])[sel],
                    ],
                    axis=1,
                ).astype(np.float32)
                _, T = thrust_axis_fast_optimized(p3, include_met=False)
                thrust_reco = float(1.0 - T)
            else:
                thrust_reco = float("nan")

            # --- Gen (status==1, drop long-lived neutrals) -------------
            st = np.asarray(mc_st[i])
            pdg = np.asarray(mc_pdg[i])
            keep = (st == 1) & ~np.isin(np.abs(pdg), list(LONG_LIVED_NEUTRALS_ABS))
            n_gen = int(keep.sum())
            if n_gen >= 2:
                p3g = np.stack(
                    [
                        np.asarray(mc_px[i])[keep],
                        np.asarray(mc_py[i])[keep],
                        np.asarray(mc_pz[i])[keep],
                    ],
                    axis=1,
                ).astype(np.float32)
                _, Tg = thrust_axis_fast_optimized(p3g, include_met=False)
                thrust_gen = float(1.0 - Tg)
            else:
                thrust_gen = float("nan")

            out.write(
                f"{run},{evn},{n_reco},{n_gen},{thrust_reco:.10g},{thrust_gen:.10g}\n"
            )

    print(f"wrote {args.output_csv}", flush=True)


if __name__ == "__main__":
    main()
