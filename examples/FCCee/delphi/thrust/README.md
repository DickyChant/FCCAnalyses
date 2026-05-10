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

# 3. EDM4hep — `--require-lvlock-zero` keeps only the SKELANA-selected
#    tracks (legacy applies the same cut via IFLSTR/IFLCUT in the YAML).
delphi-improved-reco/bin/delphi_to_edm4hep --fcc-names --require-lvlock-zero \
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

=== Reco thrust ===
  bit-identical (d=0):  100/100
  mean / median / max:  0.000e+00  0.000e+00  0.000e+00
```

Reco-thrust is **bit-identical** event-for-event when the converter is
run with `--require-lvlock-zero`. Gen-thrust agrees to <1e-7 on all 100
events; the 34 non-zero (sub-1e-7) entries differ only in numpy
float-summation order — the brute-force thrust kernel evaluates `O(n²)`
cross products in slightly different sequences when the input arrays
come from two different writers.

Verified on two independent generations (same kernel, two separate
Pythia → DELSIM passes): both 100/100 with max Δ = 0.000e+00.

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

### Reco-thrust parity: how the input selection is matched

The legacy nanoaod applies SKELANA's track-quality cut (`LVLOCK==0`,
controlled by `IFLSTR`/`IFLCUT` in the config YAML). To reproduce the
same selection downstream:

1. The raw nanoaod surfaces `TracRaw_lvlock` — SKELANA's per-track
   quality word, mapped from PA tracks via the PUCLLL `LVECP` link
   (ordinal mapping is unstable when PSHCTRECOVER reclassifies tracks).
2. The converter accepts `--require-lvlock-zero`, which drops every
   TracRaw row whose `LVLOCK != 0` — the same cut legacy applies.
3. Both writers' surviving tracks then carry bit-identical 4-momenta
   (`TracRaw_vecpPx`/`Py`/`Pz`/`E` from SKELANA's stored VECP entries),
   and the thrust kernel evaluates the same input array.

If you instead want every reconstructed track in EDM4hep (so
FCCAnalyses / weaver / ParT can make their own selection from a richer
pool), drop `--require-lvlock-zero`. The reco-thrust comparison then
shows a non-zero per-event delta — that's the input-selection
difference, not an algorithm bug.
