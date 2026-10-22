#!/usr/bin/env python3
"""Checks the DELPHI momentum-scale correction on a delphi-edm4hep file, and
shows how to use it.

    python3 check_momentum_scale.py FILE.edm4hep.root CORRECTION.json.gz [-n N] [--variation V]

Applies track_momentum_scale twice: in RDataFrame with delphi_corrections.h
(apply_momentum_scale) and in numpy with correctionlib's Python package
(momentum_scale_factors). Requires both to give the same factor for every
particle, and the corrected particles to keep their direction, their order
and their mass. On simulation (negative run numbers) every particle must come
out unchanged. Exits 1 on any failure.
"""
import argparse
import sys

import numpy as np
import uproot

import delphi_corrections as dc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("input")
    ap.add_argument("correction")
    ap.add_argument("-n", type=int, default=2000, help="events to check (default 2000)")
    ap.add_argument("--variation", default="nominal")
    ap.add_argument("--particles", default="sDST_MAIN_Particles")
    ap.add_argument("--all-tracks", action="store_true",
                    help="correct every charged track, not only those in the calibration's selection")
    args = ap.parse_args()
    window = None if args.all_tracks else dc.CALIBRATED_WINDOW

    import ROOT
    df = ROOT.RDataFrame("events", args.input).Range(args.n)
    df = dc.apply_momentum_scale(df, args.correction, variation=args.variation,
                                 particles=args.particles, output="Corrected", window=window)
    P = args.particles
    # One column per field of the corrected particles. (ROOT 6.38 types a
    # VecOps::Map over a lambda as void here, hence the explicit loops.)
    field = lambda f: f"ROOT::VecOps::RVec<float> v; for (const auto& p : Corrected) v.push_back(p.{f}); return v;"
    cols = {
        "run": "Corrected_runNumber",
        "processing": "Corrected_processing",
        "r0": "Corrected_r0",
        "s": "Corrected_momentumScale",
        "px_c": field("momentum.x"),
        "py_c": field("momentum.y"),
        "pz_c": field("momentum.z"),
        "e_c": field("energy"),
        "m_c": field("mass"),
        "q_c": field("charge"),
    }
    for k, expr in cols.items():
        df = df.Define(f"chk_{k}", expr)
    rdf = df.AsNumpy([f"chk_{k}" for k in cols])
    rdf = {k: rdf[f"chk_{k}"] for k in cols}

    t = uproot.open(args.input)["events"]
    prefix = P.split("_", 1)[0]
    V = f"{prefix}_V0_V0Candidates"
    br = lambda name: t[name].array(entry_stop=args.n, library="np")
    px, py, pz = (br(f"{P}/{P}.momentum.{c}") for c in "xyz")
    q, e, m = br(f"{P}/{P}.charge"), br(f"{P}/{P}.energy"), br(f"{P}/{P}.mass")
    vx, vy = br(f"{V}/{V}.position.x"), br(f"{V}/{V}.position.y")
    vpb, vpe = br(f"{V}/{V}.particles_begin"), br(f"{V}/{V}.particles_end")
    vab, vae = br(f"{V}/{V}.parameters_begin"), br(f"{V}/{V}.parameters_end")
    vidx = br(f"_{V}_particles/_{V}_particles.index")
    vpar = br(f"_{V}_parameters")

    import correctionlib
    corr = correctionlib.CorrectionSet.from_file(args.correction)["track_momentum_scale"]

    failures = []
    n_data = n_mc = n_charged = n_v0 = n_corrected = 0
    smin, smax = np.inf, -np.inf
    for i in range(len(px)):
        run, processing = int(rdf["run"][i]), str(rdf["processing"][i])
        r0 = dc.v0_daughter_radius(len(px[i]), vx[i], vy[i], vpb[i], vpe[i], vidx[i], vab[i], vae[i], vpar[i])
        if not np.array_equal(r0, np.asarray(rdf["r0"][i])):
            failures.append(f"entry {i}: r0 differs")
        p0 = np.stack([px[i], py[i], pz[i]], 1).astype(np.float64)
        pc = np.stack([rdf["px_c"][i], rdf["py_c"][i], rdf["pz_c"][i]], 1).astype(np.float64)
        s_rdf = np.asarray(rdf["s"][i])
        if run < 0:
            n_mc += 1
            if not (np.array_equal(pc, p0) and np.all(s_rdf == 1)):
                failures.append(f"entry {i}: simulation (run {run}) was changed")
            continue
        n_data += 1
        s_py = dc.momentum_scale_factors(corr, processing, args.variation, px[i], py[i], pz[i], q[i], r0,
                                         window=window)
        if not np.allclose(s_rdf, s_py, rtol=1e-12, atol=0):
            failures.append(f"entry {i}: factors differ, max |diff| {np.max(np.abs(s_rdf - s_py)):.3g}")
        charged = q[i] != 0
        n_charged += int(charged.sum())
        n_v0 += int((r0 > 0).sum())
        n_corrected += int((s_py != 1).sum())
        if charged.any():
            smin, smax = min(smin, s_py[charged].min()), max(smax, s_py[charged].max())
        if not np.allclose(pc, p0 * s_py[:, None], rtol=2e-7, atol=0):
            failures.append(f"entry {i}: corrected momentum is not s * p")
        if not np.array_equal(pc[~charged], p0[~charged]) or not np.array_equal(np.asarray(rdf["e_c"][i])[~charged], e[i][~charged]):
            failures.append(f"entry {i}: a neutral particle changed")
        e_expect = np.sqrt((pc**2).sum(1) + m[i].astype(np.float64)**2)
        if not np.allclose(np.asarray(rdf["e_c"][i])[charged], e_expect[charged], rtol=2e-7):
            failures.append(f"entry {i}: energy is not sqrt(p^2 + m^2)")
        if not (np.array_equal(np.asarray(rdf["m_c"][i]), m[i]) and np.array_equal(np.asarray(rdf["q_c"][i]), q[i])):
            failures.append(f"entry {i}: mass or charge changed")

    print(f"{args.input}: {n_data} data events, {n_mc} simulation events; "
          f"{n_charged} charged particles in data, {n_v0} of them V0 daughters (r0 > 0)")
    if n_charged:
        print(f"{n_corrected} tracks corrected ({'every charged track' if window is None else 'window p, theta = ' + str(window)}); "
              f"factor s ({args.variation}): {smin:.6f} to {smax:.6f}")
    print(f"{len(failures)} failure(s)")
    for f in failures[:10]:
        print("  " + f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
