"""DELPHI corrections from delphi-correctionlib, for FCCAnalyses scripts.

RDataFrame analyses:

    import delphi_corrections
    df = delphi_corrections.apply_momentum_scale(df, '/path/to/momentumScale_9495.json.gz')

uproot / numpy analyses: momentum_scale_factors() and v0_daughter_radius()
evaluate the same correction with correctionlib's Python package.

Both need correctionlib (`pip install correctionlib`); see README.md.
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# The calibration's track selection (momentumScale_skim.py): p in GeV, theta in
# degrees. Tracks outside it keep their momentum; see TrackMomentumScale in
# delphi_corrections.h. None corrects every charged track.
CALIBRATED_WINDOW = (0.4, 50.0, 20.0, 160.0)

_registered = False
_instances = 0


def register():
    """Loads correctionlib into ROOT's interpreter and declares
    delphi_corrections.h. Safe to call more than once."""
    global _registered
    if _registered:
        return
    import correctionlib
    import ROOT
    correctionlib.register_pyroot_binding()
    if not ROOT.gInterpreter.Declare(f'#include "{HERE}/delphi_corrections.h"'):
        raise RuntimeError(f"could not compile {HERE}/delphi_corrections.h")
    _registered = True


def _prefix(collection):
    """'sDST' or 'fDST': the pass a collection like sDST_MAIN_Particles belongs to."""
    return collection.split("_", 1)[0]


def apply_momentum_scale(df, file, variation="nominal", dataset="",
                         particles="sDST_MAIN_Particles", output="ReconstructedParticles",
                         window=CALIBRATED_WINDOW):
    """Defines `output` as `particles` with the charged-track momentum scale
    applied (TrackMomentumScale in delphi_corrections.h), and the columns
    <output>_momentumScale (the factor per particle) and <output>_r0.

    Data only: simulation (negative run number) passes through unchanged. The
    dataset key is the file's processing tag (<pass>_EVT_dstProcessingTag)
    unless `dataset` is given. The V0s of the same pass give r0. Only tracks
    inside `window` (p_min, p_max in GeV, theta_min, theta_max in degrees) are
    corrected.
    """
    global _instances
    register()
    import ROOT
    p = _prefix(particles)
    name = f"delphiMomentumScale{_instances}"
    _instances += 1
    pmin, pmax, tmin, tmax = window or (0.0, 0.0, 0.0, 180.0)
    ROOT.gInterpreter.Declare(
        f'const FCCAnalyses::DelphiCorrections::TrackMomentumScale {name}'
        f'{{"{os.path.abspath(file)}", "{variation}", "{dataset}", {pmin!r}, {pmax!r}, {tmin!r}, {tmax!r}}};')
    ns = "FCCAnalyses::DelphiCorrections"
    run = f"{output}_runNumber"
    processing = f"{output}_processing"
    return (df
            .Define(run, f'{ns}::intParameter(GPIntKeys, GPIntValues, "{p}_EVT_runNumber")')
            .Define(processing, f'{ns}::stringParameter(GPStringKeys, GPStringValues, '
                                f'"{p}_EVT_dstProcessingTag")')
            .Define(f"{output}_r0", f"{ns}::v0DaughterRadius({particles}, {p}_V0_V0Candidates, "
                                    f"_{p}_V0_V0Candidates_particles, _{p}_V0_V0Candidates_parameters)")
            .Define(f"{output}_momentumScale", f"{name}.factors({particles}, {output}_r0, {run}, {processing})")
            .Define(output, f"{name}({particles}, {output}_r0, {run}, {processing})"))


def v0_daughter_radius(n_particles, v0_x, v0_y, v0_particles_begin, v0_particles_end,
                       v0_particle_index, v0_parameters_begin, v0_parameters_end, v0_parameters):
    """r0 (mm) per particle of one event, as v0DaughterRadius() in
    delphi_corrections.h: the V0 radius for daughters of a V0 with tag >= 1,
    else 0. The arguments are the event's V0 branches (position.x/y,
    particles_begin/end, parameters_begin/end) and its
    _<pass>_V0_V0Candidates_particles.index and _parameters arrays."""
    r0 = np.zeros(n_particles, dtype=np.float32)
    for x, y, pb, pe, ab, ae in zip(v0_x, v0_y, v0_particles_begin, v0_particles_end,
                                    v0_parameters_begin, v0_parameters_end):
        if ab >= ae or v0_parameters[ab] < 1:
            continue
        radius = np.float32(np.hypot(x, y))
        for k in range(pb, pe):
            i = v0_particle_index[k]
            if 0 <= i < n_particles and r0[i] == 0:
                r0[i] = radius
    return r0


def momentum_scale_factors(correction, dataset, variation, px, py, pz, q, r0,
                           window=CALIBRATED_WINDOW):
    """The factor s per particle (numpy arrays of one event or flattened), with
    correctionlib's Python package: 1 for neutral particles and for tracks
    outside `window`, as in delphi_corrections.h. `correction` is
    CorrectionSet.from_file(...)["track_momentum_scale"]. The caller decides
    whether the event is data; for simulation, don't apply it."""
    px, py, pz, q = (np.asarray(a, dtype=np.float64) for a in (px, py, pz, q))
    p = np.sqrt(px**2 + py**2 + pz**2)
    s = np.ones_like(p)
    theta = np.arctan2(np.hypot(px, py), pz)
    charged = (q != 0) & (p > 0)
    if window is not None:
        pmin, pmax, tmin, tmax = window
        charged &= (p > pmin) & (p < (pmax if pmax > 0 else np.inf))
        charged &= (theta > np.radians(tmin)) & (theta < np.radians(tmax))
    if charged.any():
        s[charged] = correction.evaluate(dataset, variation, np.asarray(r0, dtype=np.float64)[charged],
                                         theta[charged], p[charged], q[charged])
    return s
