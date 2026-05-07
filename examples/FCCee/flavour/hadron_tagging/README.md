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
