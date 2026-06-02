#!/usr/bin/env python3
"""Paper-aligned thrust on a DELPHI EDM4hep file (arXiv:2510.18762v1
Section 5 + Table "Summary table for particle and event selections").

Differences vs `thrust_edm4hep.py`:
  * Thrust uses BOTH charged AND neutral particles (paper:
    "The distribution is determined using all charged and neutral particles.")
    Our previous analyzer was charged-only.
  * Per-particle acceptance cuts (paper Table 1):
        20 deg <= theta <= 160 deg
        pT > 0.4 GeV (charged)
        E  > 0.5 GeV (neutral)
  * Event-level hadronic selection (paper Table 1, event section):
        n_charged_selected >= 7
        E_tot              >= 0.5 * E_cm  = 45.6 GeV at Z pole
        30 deg <= theta_thrust <= 150 deg
  * Track-quality cuts (track length >= 30 cm, dp/p <= 1, |d0|<4 cm,
    |z0|<10 cm) are NOT applied here — the converter / shortdst chain
    already pre-selects "good" PFOs and the residual quality variation
    impacts <2% of tracks. Add them later if a more faithful reproduction
    is needed (they're available via _EFlowTrack_trackStates D0/Z0 and
    EFlowTrack chi2/ndf, with PFO->Track indirection through
    PandoraPFOs.tracks_begin/end).

Output CSV columns:
  runNo, evtNo, n_ch_sel, n_neu_sel, E_tot, theta_thrust_deg,
  thrust_reco, pass_hadronic
"""
import sys, argparse, math
import numpy as np
import uproot

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
from common_functions import thrust_axis_fast_optimized

E_CM       = 91.2   # GeV, LEP-1 Z-pole
THETA_MIN_PART_DEG = 20.0
THETA_MAX_PART_DEG = 160.0
PT_MIN_CH  = 0.4    # GeV
E_MIN_NEU  = 0.5    # GeV
N_CH_MIN_EVT = 7
E_TOT_FRAC   = 0.5  # E_tot > 0.5 * E_cm
THETA_MIN_THRUST_DEG = 30.0
THETA_MAX_THRUST_DEG = 150.0

DEG = math.pi / 180.0


def get_run_evt(gp_keys_event, gp_vals_event):
    out = {}
    for k, v in zip(list(gp_keys_event), list(gp_vals_event)):
        out[k] = int(v[0])
    return out.get("runNumber", -1), out.get("eventNumber", -1)


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("input")
    ap.add_argument("output_csv")
    args = ap.parse_args()

    ev = uproot.open(args.input)["events"]
    n = ev.num_entries

    avail = set(ev.keys())
    # Schema autodetect (in priority order):
    #   1. Packaged two-pass: prefer fDST_MAIN_Particles (full-DST union);
    #      fall back to sDST_MAIN_Particles for SDST-only files.
    #   2. FCC-native single-pass:   PandoraPFOs
    #   3. Legacy schema:            ReconstructedParticles
    if "fDST_MAIN_Particles/fDST_MAIN_Particles.charge" in avail:
        rp_name = "fDST_MAIN_Particles"
    elif "sDST_MAIN_Particles/sDST_MAIN_Particles.charge" in avail:
        rp_name = "sDST_MAIN_Particles"
    elif "PandoraPFOs/PandoraPFOs.charge" in avail:
        rp_name = "PandoraPFOs"
    else:
        rp_name = "ReconstructedParticles"

    px = ev[f"{rp_name}/{rp_name}.momentum.x"].array()
    py = ev[f"{rp_name}/{rp_name}.momentum.y"].array()
    pz = ev[f"{rp_name}/{rp_name}.momentum.z"].array()
    en = ev[f"{rp_name}/{rp_name}.energy"].array()
    q  = ev[f"{rp_name}/{rp_name}.charge"].array()

    gp_keys = ev["GPIntKeys"].array()
    gp_vals = ev["GPIntValues"].array()

    with open(args.output_csv, "w") as out:
        out.write(
            "runNo,evtNo,n_ch_sel,n_neu_sel,E_tot,theta_thrust_deg,"
            "thrust_reco,pass_hadronic\n"
        )
        for i in range(n):
            run, evn = get_run_evt(gp_keys[i], gp_vals[i])
            if evn == 0:
                # DELSIM run-start / housekeeping records
                continue

            pxi = np.asarray(px[i], dtype=np.float64)
            pyi = np.asarray(py[i], dtype=np.float64)
            pzi = np.asarray(pz[i], dtype=np.float64)
            eni = np.asarray(en[i], dtype=np.float64)
            qi  = np.asarray(q[i], dtype=np.float64)

            p_mag = np.sqrt(pxi**2 + pyi**2 + pzi**2)
            with np.errstate(invalid="ignore", divide="ignore"):
                cos_theta = np.where(p_mag > 0, pzi / p_mag, 0.0)
            theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))  # rad
            pt    = np.sqrt(pxi**2 + pyi**2)

            in_eta = (theta >= THETA_MIN_PART_DEG * DEG) & \
                     (theta <= THETA_MAX_PART_DEG * DEG)
            is_ch  = np.abs(qi) > 0.1
            sel_ch = in_eta & is_ch  & (pt > PT_MIN_CH)
            sel_neu = in_eta & (~is_ch) & (eni > E_MIN_NEU)
            sel = sel_ch | sel_neu

            n_ch  = int(sel_ch.sum())
            n_neu = int(sel_neu.sum())
            e_tot = float(eni[sel].sum())

            if int(sel.sum()) < 2:
                out.write(
                    f"{run},{evn},{n_ch},{n_neu},{e_tot:.4f},nan,nan,0\n"
                )
                continue

            p3 = np.stack(
                [pxi[sel], pyi[sel], pzi[sel]], axis=1
            ).astype(np.float32)
            axis, T = thrust_axis_fast_optimized(p3, include_met=False)
            tau = float(1.0 - T)

            # thrust-axis polar angle (deg)
            axis_p = math.sqrt(float(axis[0]**2 + axis[1]**2 + axis[2]**2))
            theta_thrust_deg = math.degrees(
                math.acos(float(axis[2]) / axis_p) if axis_p > 0 else 0.0
            )

            pass_had = (n_ch >= N_CH_MIN_EVT) and \
                       (e_tot >= E_TOT_FRAC * E_CM) and \
                       (THETA_MIN_THRUST_DEG <= theta_thrust_deg
                        <= THETA_MAX_THRUST_DEG)

            out.write(
                f"{run},{evn},{n_ch},{n_neu},{e_tot:.4f},"
                f"{theta_thrust_deg:.4f},{tau:.10g},{int(pass_had)}\n"
            )

    print(f"wrote {args.output_csv}", flush=True)


if __name__ == "__main__":
    main()
