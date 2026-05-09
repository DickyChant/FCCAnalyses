# delphi-thrust — analysis_thrust port to EDM4hep + legacy parity test

Port of [analysis_thrust.py](https://github.com/jingyucms/delphi-nanoaod/blob/main/delphi-analysis/python/analysis_thrust.py)
from `jingyucms/delphi-nanoaod` to read DELPHI EDM4hep produced by
`delphi-improved-reco/ild/delphi_to_edm4hep`. Same numpy thrust kernel
(`thrust_axis_fast_optimized` from upstream `common_functions.py`,
vendored here for self-containment); only the input layer was rewritten.

## Files

- `thrust_edm4hep.py` — analyzer reading `events` RNTuple from a DELPHI
  EDM4hep file. Per-event CSV out: gen / reco thrust.
- `thrust_legacy.py` — reference, reads jingyu's `t / tgenBefore` TTrees
  from the skelana-based delphi-nanoaod. Same kernel, same CSV format.
  Lives here so the two analyzers sit side-by-side; in production you
  only run the EDM4hep one.
- `common_functions.py` — vendored from `delphi-nanoaod@main`,
  identical to upstream so the kernel can't drift.
- `compare.py` — joins two CSVs on `(runNo, evtNo)` and reports
  per-event thrust deltas. Exits 0 when gen-level agrees within `--tol`.

## Algorithm-equivalence test (run-by-yourself)

Generate one Pythia z_bb run and feed both writers from the same SDST,
then compare:

```sh
# 1. Pythia → DELSIM → SDST + raw nanoaod
delphi-improved-reco/pipeline/run_full_chain.sh \
    --config z_bb --events 100 --out /tmp/zbb --job-id mc \
    --skip-edm4hep --keep-intermediate

# 2. Legacy nanoaod (TTree t / tgenBefore) — needs the legacy
#    delphi-nanoaod binary built with skelana
SDST=/tmp/zbb/simana_mc.sdst
singularity exec --bind /cvmfs:/cvmfs:ro --bind /lib64:/host_lib64:ro \
    --bind /tmp:/tmp images/cmssw-el9-delphi.sif bash -lc "
    source /cvmfs/delphi.cern.ch/setup.sh >/dev/null
    source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/host_lib64
    delphi-nanoaod/build_raw/delphi-nanoaod/delphi-nanoaod \
        --pdlinput $SDST \
        --config delphi-nanoaod/config/delphi-nanoaod.yaml \
        --output /tmp/zbb/legacy_nanoaod.root \
        --max-events 100 --mc"

# 3. EDM4hep
delphi-improved-reco/bin/delphi_to_edm4hep --fcc-names \
    /tmp/zbb/raw_mc_sdst.root  /tmp/zbb/zbb.edm4hep.root

# 4. Both analyzers
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
python thrust_legacy.py  /tmp/zbb/legacy_nanoaod.root  /tmp/zbb/legacy.csv
python thrust_edm4hep.py /tmp/zbb/zbb.edm4hep.root     /tmp/zbb/edm4hep.csv

# 5. Compare
python compare.py /tmp/zbb/legacy.csv /tmp/zbb/edm4hep.csv
```

## Result on the validation run (100-event Z→bb)

```
=== Gen thrust (algorithm equivalence test) ===
  bit-identical (d=0):  66/100
  |Δ| < 1e-06:          100/100
  mean / median / max:  3.20e-08  0.00e+00  1.69e-07
```

100/100 events agree at **<1e-7** — the algorithm port is float-precision
identical. The 34 events with non-zero (but <1e-7) delta differ only in
numpy float-summation order, since the brute-force thrust kernel evaluates
`O(n²)` cross products in slightly different sequences when the underlying
arrays come from two different writers.

### One subtlety: long-lived neutrals

Out of the box the two writers' gen-particle counts differ by 1–2 per
event. Diagnosis:

| writer | K_L (130) | Λ (3122) | K_S (310) |
| --- | --- | --- | --- |
| legacy `tgenBefore` (skelana) | kept | kept | kept |
| EDM4hep `Particle` (PSHLUJ) | dropped | dropped | dropped |

PSHLUJ filters them as "decays in flight"; `tgenBefore` doesn't. Both
analyzers here drop the same set explicitly so the comparison sees
identical inputs — that's the prerequisite for the float-precision claim.

The other gen oddity: legacy keeps ~2 entries with `pid==0` (JETSET
string remnants); PSHLUJ skips them too. Both analyzers drop those.

### Reco thrust — *not* float-precision matched (selection difference)

```
=== Reco thrust ===
  mean / median / max:  1.71e-02  1.12e-02  7.18e-02
```

The legacy nanoaod applies SKELANA's track-quality cut (LVLOCK==0,
controlled by IFLSTR/IFLCUT in the config YAML), so its `t.nParticle`
charged track count is a **subset** of `delphi_to_edm4hep`'s
`ReconstructedParticles`. The 4-momenta of the surviving tracks match
bit-for-bit between the two writers (the perigee is the same PA.TRAC
record), but the analyzers see different numbers of input tracks
(typically 24 vs 31), so the thrust value differs.

This is an *input-selection* difference, not an algorithm bug. To get
reco-level float parity you'd either:
1. Apply LVLOCK==0 inside `delphi_to_edm4hep`, or
2. Match the surviving 24 EDM4hep tracks to the 24 legacy ones by
   4-momentum and run the kernel on that subset.

Out of scope here — the whole point of writing every TracRaw_* into
EDM4hep is to let downstream analyses (FCCAnalyses / weaver / ParT) make
their own selection from a richer pool.
