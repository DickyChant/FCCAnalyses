# hadron_tagging — modern key4hep stack port (full feature set)

Modern-stack port of `hadron_tagging` from
[zuoxunwu/FCCAnalyses@hadronTagger_dev](https://github.com/zuoxunwu/FCCAnalyses/tree/hadronTagger_dev/examples/FCCee/flavour/hadron_tagging).
The original was on the centos7 / EDM4hep 0.7 stack; this branch runs on
the modern key4hep stack (FCCAnalyses 0.12.x, EDM4hep 1.0) and uses DELPHI
EDM4hep produced by `delphi-improved-reco`'s `delphi_to_edm4hep` converter
(branch `feature/edm4hep-pipeline`).

## Running

```sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
fccanalysis run examples/FCCee/flavour/hadron_tagging/analysis_stage1.py
```

The default `inputDir` points at our local DELPHI Z→bb output. `processList`
sets the file basename without `.edm4hep.root` extension; a symlink
`events_phA_bb_fcc.root → events_phA_bb_fcc.edm4hep.root` is required since
`fccanalysis` appends `.root`.

## What this does

Full-feature stage1 ntuple, restored from the legacy version. Event
selection is `EVT_hasPV==1` (~96% pass on Z→bb).

- **Gen-level B-hadron / Λ_b counts + thrust hemisphere assignment**
  (`n_genBs/Bu/Bd/Bc/Lb`, `genB?_thrustangle`, `label_B?_Emin/Emax`).
- **Reco P kinematics + helix params at IP** (`RP_e`, `RP_p{x,y,z}`,
  `RP_charge`, `RP_trk_d0/z0/phi/omega/tanLambda`).
- **TPC dE/dx per RP** (`RP_dndx`) — pulled from the modern
  `EFlowTrack_dNdx` (RecDqdx) collection via the explicit
  `_EFlowTrack_dNdx_track` index relation.
- **Reco→MC truth match** (`RP_MCidx`, `RP_nMC`, `RP_fromBs/Bu/Bd/Bc/Lb`)
  via the `MCRecoAssociations` collection produced by `delphi_to_edm4hep`'s
  helix matcher (greedy nearest-neighbour in (theta, phi, p), ~85% MC-side
  match rate on Z→bb).
- **Reco vertex object** (`Vertex_*`) built by `myUtils::get_VertexObject`
  using the MC vertex tree + RP-MC assoc + tracks.
  - Modern EDM4hep 1.0 stores per-track AtIP states in
    `_EFlowTrack_trackStates` (one per track in our converter), passed
    directly as the TrackState input — replaces the legacy `EFlowTrack_1`
    flat-table convention.
- **Thrust + hemisphere energy / multiplicity splits**
  (`EVT_thrust*`, `EVT_ThrustEmin/Emax_*`).
- **PV — reco vs. MC offset** (`PV_x/y/z`, `PV_x/y/z_offset`).

## Cross-check (Z→bb vs. Z→light, 100 events each)

| | Z→bb | Z→light |
| --- | --- | --- |
| events passing PV filter | 100/104 | 100/104 |
| total reco PFOs | 2987 | 2696 |
| RPs tagged from-Bs / Bd / Bu / Λb | 81 / 442 / 324 / 33 | 0 / 5 / 3 / 0 |
| dE/dx population (>0) | 78% | 79% |
| dE/dx mean | 1.85 | 1.90 |

The flavor labels cleanly separate Z→bb (heavy-flavor) from Z→light
(prompt-only). RP-level d0 distribution centered at 0 with hadronic-decay
tail.

## LVLOCK selection (DELPHI-native variant only)

`analysis_stage1_delphi_native.py` exposes SKELANA's per-PFO `Track_lvlock`
quality word as two parallel columns on the stage1 ntuple:

- `RP_lvlock` — the raw signed-int32 bitmask, in `RecoPartPIDAtVertex`
  order. `0` = pass, `1` = LVSELE fail, `INT32_MIN` = REMCLU
  calo-cluster overlap (bit 32), `INT32_MIN+1` = both, `-1` = nanoaod
  predates the field.
- `RP_passLvlock` — `0/1` mask, `1 ⇔ lvlock == 0 || lvlock == -1`.
  Apply this as a per-RP cut at the analysis output level (e.g., when
  filling per-RP histograms or selecting BDT inputs).

Do **not** filter `ReconstructedParticles` with `filterRPbyLvlock` before
calling `myUtils::get_VertexObject` / `myUtils::PID` / `getRP2MC_index`
— those helpers consume `MCRecoAssociations*.index` written against
the original `PandoraPFOs` ordering and segfault when the collection
is silently re-indexed under them (observed on
`/eos/experiment/eealliance/Users/zhangj/edm4hep/edm4hep/events_qq_*.edm4hep.root`,
which is the first sample to actually carry non-zero LVLOCK bits).
The `lvlockPassMask` helper added in `functions.h` produces the
parallel-to-PandoraPFOs flag without rewriting the collection.

Always use `RP_passLvlock == 1` rather than `RP_lvlock <= 0`: the
REMCLU class (bit 32) decodes as INT32_MIN which is numerically
negative, so `<=0` keeps it by mistake.

## Modern-stack quirks

Several FCCAnalyses helpers depend on EDM4hep 0.x field names that were
renamed / replaced in 1.0. The local `functions.h` works around the gap:

- `vertex.primary` → bit-encoded in `vertex.type`. The modern
  `VertexFitter` writes `type = Primary` literally instead of using the
  bit, so `myUtils::get_Vertex_isPV` and `myUtils::hasPV` always return
  false. Local replacements `get_Vertex_isPV` / `hasPV` use
  `FCCAnalysesVertex.mc_ind == 0` instead.
- `TrackData.dxQuantities_begin` (legacy dE/dx index) → modern
  `RecDqdxCollection` is its own collection with explicit
  `_EFlowTrack_dNdx_track.index` relation. Local `get_RP_dndx` walks
  this relation.
- `TrackerHitData` → `TrackerHit3DData` (renamed).
- `RP::type` → gone.

## Reduced features (vs. original) — still dropped

- **TMVA RBDT inference** (xzuo's afs path). Re-add when the model
  is available from a stable location.
- **`RP_mtof`** — DELPHI shortDST has no calorimeter timing.

## Reference

- Modern path output: `/eos/user/s/sqian/www/delphi_edm4hep/phaseA-fcc/`
  - web: https://sqian.web.cern.ch/sqian/delphi_edm4hep/phaseA-fcc/
- Legacy path (centos7 stack, no MCParticles): `/eos/user/s/sqian/www/delphi_edm4hep/phaseA-legacy/`
