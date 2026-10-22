# DELPHI corrections (delphi-correctionlib)

Applies the corrections of
[delphi-correctionlib](https://github.com/delphi-fullDST-edm4hep/delphi-correctionlib)
to DELPHI EDM4hep files written by
[delphi-edm4hep](https://github.com/delphi-fullDST-edm4hep/delphi-edm4hep), in
RDataFrame analyses and in uproot/numpy scripts. So far there is one
correction: the charged-track momentum scale (`track_momentum_scale` in
`tracking/momentumScale_9495.json.gz`), for the 94C2 and 95D2 processings.

| File | What |
|---|---|
| `delphi_corrections.h` | C++ for RDataFrame: `TrackMomentumScale`, `v0DaughterRadius`, frame-parameter getters |
| `delphi_corrections.py` | `apply_momentum_scale(df, …)` for RDataFrame; `momentum_scale_factors` and `v0_daughter_radius` for numpy |
| `check_momentum_scale.py` | checks the two against each other on a file; also a usage example |

## Setup

correctionlib is not in the Key4hep release (2026-04-08), so install it next to
it:

```sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
python3 -m pip install --target ~/.local/correctionlib correctionlib==2.9.0
export PYTHONPATH=~/.local/correctionlib:$PYTHONPATH
```

The pip package carries the C++ library and headers; `delphi_corrections.py`
loads them into ROOT with `correctionlib.register_pyroot_binding()`. Get the
correction file from a clone of delphi-correctionlib, or download
`tracking/momentumScale_9495.json.gz` from it.

## Use

RDataFrame:

```python
import delphi_corrections
df = delphi_corrections.apply_momentum_scale(
    df, "momentumScale_9495.json.gz", variation="nominal",
    particles="sDST_MAIN_Particles", output="ReconstructedParticles")
```

This defines `ReconstructedParticles` as the corrected particles, in the same
order and with the same relations as `sDST_MAIN_Particles`, so index-based
helpers keep working. It also defines `ReconstructedParticles_momentumScale`
(the factor per particle) and `ReconstructedParticles_r0`. Run it once per
systematic variation: `nominal`, `scale_up`, `scale_down`, `alignment_up`,
`alignment_down`.

Wired in already: `delphi/thrust/thrust_edm4hep_arxiv.py`, with
`--momentum-scale FILE` (`--momentum-scale-variation`,
`--momentum-scale-all-tracks`).

## What it does

For each charged particle, the momentum is multiplied by
`s(dataset, variation, r0, theta, p, q)` and the energy recomputed from the
particle's mass. The direction is unchanged, and neutral particles are untouched.

- **Data only.** Events with a negative run number (`<pass>_EVT_runNumber`,
  i.e. simulation) pass through unchanged.
- **dataset** is the file's processing tag (`<pass>_EVT_dstProcessingTag`, e.g.
  `94C2`), unless given. A processing the correction does not cover stops the
  job with an error.
- **r0** is the xy radius of the V0 vertex for daughters of a DELPHI V0
  (`<pass>_V0_V0Candidates`) with tag ≥ 1, as in the calibration; 0 otherwise.
- **Only tracks in the calibration's selection are corrected:**
  0.4 < p < 50 GeV and 20° < θ < 160°, the cuts of `momentumScale_skim.py`.
  Outside it, the model extrapolates. Its energy-loss term grows as 1/p_T: on
  94C2 data it gives s < 0.5 for tracks of a few tens of MeV, and s < 0 for
  some, while inside the selection s stays within 0.96–1.03. Tracks outside
  keep their momentum. `window=None` (`--momentum-scale-all-tracks`) corrects
  every charged track.

## Validation

`check_momentum_scale.py FILE momentumScale_9495.json.gz` applies the
correction through RDataFrame and through numpy. It requires the same factor
to 1e-12, s·p for the corrected momenta, E = √(p² + m²), and unchanged
neutrals, masses and charges. On simulation, nothing may change. On
delphi-edm4hep output (3000 events each):

| Sample | Corrected tracks | s (nominal) | Result |
|---|---|---|---|
| 94C2 data, `short94_c2` | 15523 of 22897 charged | 0.968 – 1.027 | pass |
| 95D2 data, `long95p3_d2` | 2087 of 4172 | 0.951 – 1.083 | pass |
| 94C2 simulation, `sh_qqps_r94_2l_c2` | none (simulation) | — | pass |
| 99D1 data, `xsdst99_d1` | — | — | stops: processing not covered |

In `thrust_edm4hep_arxiv.py` on the 94C2 file, the correction moves 1 − T by
2 × 10⁻⁵ (median) in hadronic events, and by at most 10⁻³ when the selected
tracks stay the same.
