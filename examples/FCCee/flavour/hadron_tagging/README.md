# hadron_tagging — modern key4hep stack port

This is the modern-stack port of `hadron_tagging` from
[zuoxunwu/FCCAnalyses@hadronTagger_dev](https://github.com/zuoxunwu/FCCAnalyses/tree/hadronTagger_dev/examples/FCCee/flavour/hadron_tagging).
The original was on the centos7 / EDM4hep 0.7 stack and targeted the centrally
produced FCC IDEA winter2023 samples; this version runs on the modern key4hep
stack (FCCAnalyses 0.12.x, EDM4hep 1.0) and uses DELPHI EDM4hep produced by
`delphi-improved-reco`'s `delphi_to_edm4hep` converter (branch
`feature/edm4hep-pipeline`).

## Running

```sh
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
fccanalysis run examples/FCCee/flavour/hadron_tagging/analysis_stage1.py
```

The default `inputDir` points at our local DELPHI Z→bb output:

```
/eos/user/s/sqian/www/delphi_edm4hep/phaseA-fcc/zbb/events_phA_bb_fcc.root
```

## What this does

A reduced-feature smoke test that verifies the modern FCCAnalyses toolchain
runs end-to-end on DELPHI EDM4hep. It:

- counts gen B-hadron and Λb species (Bs, Bu, Bd, Bc, Λb)
- assigns per-thrust-hemisphere truth labels (label_Bs_Emin, label_Bu_Emin, …)
- computes basic reco kinematics (RP_e, RP_px, RP_py, RP_pz, RP_charge, …)
- computes thrust + hemisphere energy/multiplicity splits
- computes the legacy reco PV — MC PV offset (PV_x_offset, PV_y_offset, PV_z_offset)

Verified Z→bb vs Z→light separation at gen level (0 vs 2 B's/evt).

## Differences from the original (legacy) version

| feature | original | modern smoke test |
| --- | --- | --- |
| stack | spackages6 / centos7 / EDM4hep 0.7 | key4hep / el9 / EDM4hep 1.0 |
| input sample | FCC IDEA winter2023 SimDelphes | DELPHI shortDST → EDM4hep |
| BDT inference | TMVA RBDT (xzuo's afs path) | dropped |
| MCRecoAssociations features | yes | dropped (collection not in our converter output) |
| EFlowTrack_1/2 (track state, dNdx) | yes | dropped (modern schema; would need re-aliasing) |
| Vertex_* (recovered SVs, mass, d2PV, …) | yes | dropped (depends on get_VertexObject) |
| RP_trk_d0/z0/phi/omega/tanLambda | yes | dropped |
| RP_dndx / RP_mtof | yes | dropped (no calorimeter timing in DELPHI) |
| gen B counts + thrust labels | yes | yes |
| reco P4 + thrust hemispheres | yes | yes |

## To restore full feature set

The dropped features are recoverable but each depends on adding the missing
collection / association to the converter:

1. **MCRecoAssociations**: add to `delphi_to_edm4hep` a per-RP truth match
   built from the `_ReconstructedParticles_tracks` ↔ `_Particle` map. Many
   DELPHI samples link via the dominant TPC track of each RP back to a single
   MCParticle; PV and SV consumers need this.
2. **EFlowTrack_1/2 aliases**: in modern EDM4hep 1.0 the track state lives in
   `_EFlowTrack_trackStates` and dE/dx in `EFlowTrack_dNdx` (RecDqdx).
   The legacy FCC analyses' `EFlowTrack_1` / `EFlowTrack_2` aliases would
   need to be replaced with calls into modern `ReconstructedTrack::*` helpers.
3. **functions.h**: trimmed for this smoke test. The original had helpers
   keyed off EDM4hep 0.x field names (`vertex.primary`,
   `TrackData.dxQuantities_begin`, `TrackerHitData`, `RP::type`); restore
   piecewise as needed.
4. **BDT inference**: load the trained model from a robust location (not afs)
   and re-add the `computeModel1` line. This is decoupled from the schema
   migration above.
