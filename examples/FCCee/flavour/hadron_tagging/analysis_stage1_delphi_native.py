'''
hadron_tagging stage1 — DELPHI-native collection-name variant.

Same body as `analysis_stage1.py`, but reads EDM4hep files that came out of
`delphi-improved-reco/ild/delphi_to_edm4hep` *without* the `--fcc-names`
flag — i.e. with the DELPHI-native collection names

    MCParticles  PandoraPFOs  Tracks       _Tracks_trackStates
    _MCParticles_parents  _MCParticles_daughters

instead of the FCC names

    Particle     ReconstructedParticles  EFlowTrack  _EFlowTrack_trackStates
    _Particle_parents     _Particle_daughters

Bridge: a small block of RDataFrame `.Alias()` calls maps the file-side
branch names onto the FCC names that all the FCCAnalyses helpers expect.

Same coverage as the FCC-renamed variant once the converter emits
`Tracks_dNdx` (parallel RecDqdx form, alongside the rich
`ParticleID_dEdx`). The stage1 picks it up via the
`EFlowTrack_dNdx`/`_EFlowTrack_dNdx_track` aliases below and the
`get_RP_dndx` helper.

Per-PFO PID flags (`RP_isMu`, `RP_isEl`, `RP_hasRich`) are pulled from
the `_ParticleID_Muon_particle` / `_ParticleID_Electron_particle` /
`_ParticleID_HadronRich_particle` index relations via the local
`hasPIDLink` helper in `functions.h`.

Verified on smoke sample
  /eos/experiment/eealliance/Users/zhangj/edm4hep/edm4hep/events_bb_*.edm4hep.root
(~95% pass the EVT_hasPV==1 filter, ≈19/20 on a 20-event smoke).
'''

# Mandatory: list of processes
processList = {
    'events_bb_1778449473_00000.edm4hep': {'output': 'smoke_zbb_delphi_native'},
}

# Mandatory: production tag (None when running on local files)
prodTag = None

# Mandatory: input directory — change to your DELPHI-native EDM4hep dir
inputDir = '/eos/experiment/eealliance/Users/zhangj/edm4hep/edm4hep'

# Mandatory: output directory
outputDir = 'outputs/FCCee/flavour/hadron_tagging/stage1_delphi_native'

# Optional: include local helpers
includePaths = ['functions.h']

# Optional: number of threads
nCPUS = 4


class RDFanalysis():

    def analysers(df):
        # =====================================================================
        # Schema detection: the post-beff550 converter writes FCC-native
        # collection names (EFlowTrack / EFlowPhoton / EFlowNeutralHadron /
        # ParticleID_dEdx) and dropped the per-track lvlock + MC-Reco
        # association collections that the older zhangj samples carried.
        # We branch the alias/stub block on schema; the rest of the
        # analyzer (Defines / Filters / output) is shared.
        #
        # On new-schema files the MC-Reco link is gone, so per-RP truth
        # features (RP_MCidx, RP_fromB*, Vertex_fromB*) come out as zero
        # arrays. That's acceptable for tagger inference and for the
        # parton-flavour retrain (channel-based labels, no per-RP truth).
        # =====================================================================
        cols = set(df.GetColumnNames())
        new_schema = ('_EFlowTrack_trackStates' in cols
                      and '_Tracks_trackStates' not in cols)

        if new_schema:
            df = (df
                  # MCParticles / PandoraPFOs are named the same in both
                  # schemas; only their downstream FCC-name aliases change.
                  .Alias('Particle',               'MCParticles')
                  .Alias('ReconstructedParticles', 'PandoraPFOs')

                  # Track_lvlock dropped — stub to zeros so every PFO passes.
                  .Define('Track_lvlock',
                          'ROOT::VecOps::RVec<int>(EFlowTrack.size(), 0)')
                  .Define('RP_passLvlock_input',
                          'FCCAnalyses::ZHfunctions::lvlockPassMask('
                          'PandoraPFOs, Track_lvlock)')

                  # dN/dx moved into ParticleID_dEdx in the new schema.
                  # Stub the *legacy* (Tracks_dNdx + _Tracks_dNdx_track)
                  # names with empty vectors so the existing get_RP_dndx
                  # call returns all-zero dN/dx vectors.
                  .Define('EFlowTrack_dNdx',
                          'ROOT::VecOps::RVec<edm4hep::Quantity>{}')
                  .Define('_EFlowTrack_dNdx_track',
                          'ROOT::VecOps::RVec<podio::ObjectID>{}')

                  # MC-particle index relations — same field names in both.
                  .Alias('Particle0', '_MCParticles_parents.index')
                  .Alias('Particle1', '_MCParticles_daughters.index')

                  # _MCRecoAssociations_{from,to} dropped — stub empty so
                  # get_VertexObject / PID / getRP2MC_index get zero-length
                  # association vectors. Downstream features that depend on
                  # the link (RP_MCidx, RP_fromB*, Vertex_fromB*) come out
                  # as zero arrays of the appropriate length.
                  .Define('MCRecoAssociations0',
                          'ROOT::VecOps::RVec<int>{}')
                  .Define('MCRecoAssociations1',
                          'ROOT::VecOps::RVec<int>{}')
                  )
        else:
            df = (df
                  #############################################
                  ## Legacy DELPHI-native (zhangj) schema aliases.
                  #############################################
                  .Alias('Particle',                 'MCParticles')
                  .Alias('_EFlowTrack_trackStates',  '_Tracks_trackStates')

                  ## DELPHI-legacy track-quality cut: SKELANA `LVLOCK == 0`.
                  ## See repo CLAUDE.md for the convention; the writer emits
                  ## every PFO so the cut lives at the FCCAnalyser entry
                  ## point and is exposed as RP_passLvlock_input (don't
                  ## shrink the RP collection before vertexing — relations
                  ## point into the original PandoraPFOs ordering).
                  .Alias('ReconstructedParticles',   'PandoraPFOs')
                  .Define('RP_passLvlock_input',
                          'FCCAnalyses::ZHfunctions::lvlockPassMask('
                          'PandoraPFOs, Track_lvlock)')
                  .Alias('EFlowTrack_dNdx',          'Tracks_dNdx')
                  .Alias('_EFlowTrack_dNdx_track',   '_Tracks_dNdx_track')

                  #############################################
                  ## Index-relation aliases (file-side branch names)
                  #############################################
                  .Alias('Particle0',           '_MCParticles_parents.index')
                  .Alias('Particle1',           '_MCParticles_daughters.index')
                  .Alias('MCRecoAssociations0', '_MCRecoAssociations_from.index')
                  .Alias('MCRecoAssociations1', '_MCRecoAssociations_to.index')
                  )

        df2 = (
            df

            #############################################
            ## MC bookkeeping
            #############################################
            .Define('MC_PDG', 'FCCAnalyses::MCParticle::get_pdg(Particle)')
            .Define('MC_n',   'int(MC_PDG.size())')
            .Define('MC_M1',  'myUtils::get_MCMother1(Particle, Particle0)')
            .Define('MC_M2',  'myUtils::get_MCMother2(Particle, Particle0)')
            .Define('MC_D1',  'myUtils::get_MCDaughter1(Particle, Particle1)')
            .Define('MC_D2',  'myUtils::get_MCDaughter2(Particle, Particle1)')

            #############################################
            ## Gen B-hadron / b-quark counters
            #############################################
            .Define('genBottom',     'FCCAnalyses::MCParticle::sel_pdgID(5, true)(Particle)')
            .Define('n_genBottoms',  'FCCAnalyses::MCParticle::get_n(genBottom)')
            .Define('genBottom_px',  'FCCAnalyses::MCParticle::get_px(genBottom)')
            .Define('genBottom_py',  'FCCAnalyses::MCParticle::get_py(genBottom)')
            .Define('genBottom_pz',  'FCCAnalyses::MCParticle::get_pz(genBottom)')
            .Define('genBottom_pdg', 'FCCAnalyses::MCParticle::get_pdg(genBottom)')

            .Define('genBs', 'FCCAnalyses::MCParticle::sel_pdgID(531, true)(Particle)')
            .Define('n_genBs', 'FCCAnalyses::MCParticle::get_n(genBs)')
            .Define('genBs_px', 'FCCAnalyses::MCParticle::get_px(genBs)')
            .Define('genBs_py', 'FCCAnalyses::MCParticle::get_py(genBs)')
            .Define('genBs_pz', 'FCCAnalyses::MCParticle::get_pz(genBs)')

            .Define('genBu', 'FCCAnalyses::MCParticle::sel_pdgID(521, true)(Particle)')
            .Define('n_genBu', 'FCCAnalyses::MCParticle::get_n(genBu)')
            .Define('genBu_px', 'FCCAnalyses::MCParticle::get_px(genBu)')
            .Define('genBu_py', 'FCCAnalyses::MCParticle::get_py(genBu)')
            .Define('genBu_pz', 'FCCAnalyses::MCParticle::get_pz(genBu)')

            .Define('genBd', 'FCCAnalyses::MCParticle::sel_pdgID(511, true)(Particle)')
            .Define('n_genBd', 'FCCAnalyses::MCParticle::get_n(genBd)')
            .Define('genBd_px', 'FCCAnalyses::MCParticle::get_px(genBd)')
            .Define('genBd_py', 'FCCAnalyses::MCParticle::get_py(genBd)')
            .Define('genBd_pz', 'FCCAnalyses::MCParticle::get_pz(genBd)')

            .Define('genBc', 'FCCAnalyses::MCParticle::sel_pdgID(541, true)(Particle)')
            .Define('n_genBc', 'FCCAnalyses::MCParticle::get_n(genBc)')
            .Define('genBc_px', 'FCCAnalyses::MCParticle::get_px(genBc)')
            .Define('genBc_py', 'FCCAnalyses::MCParticle::get_py(genBc)')
            .Define('genBc_pz', 'FCCAnalyses::MCParticle::get_pz(genBc)')

            .Define('genLb', 'FCCAnalyses::MCParticle::sel_pdgID(5122, true)(Particle)')
            .Define('n_genLb', 'FCCAnalyses::MCParticle::get_n(genLb)')
            .Define('genLb_px', 'FCCAnalyses::MCParticle::get_px(genLb)')
            .Define('genLb_py', 'FCCAnalyses::MCParticle::get_py(genLb)')
            .Define('genLb_pz', 'FCCAnalyses::MCParticle::get_pz(genLb)')

            #############################################
            ## Gen C-hadron / c-quark counters
            ## (inclusive c-hemisphere truth tag)
            #############################################
            .Define('genCharm',     'FCCAnalyses::MCParticle::sel_pdgID(4, true)(Particle)')
            .Define('n_genCharms',  'FCCAnalyses::MCParticle::get_n(genCharm)')
            .Define('genCharm_px',  'FCCAnalyses::MCParticle::get_px(genCharm)')
            .Define('genCharm_py',  'FCCAnalyses::MCParticle::get_py(genCharm)')
            .Define('genCharm_pz',  'FCCAnalyses::MCParticle::get_pz(genCharm)')
            .Define('genCharm_pdg', 'FCCAnalyses::MCParticle::get_pdg(genCharm)')

            .Define('genD0', 'FCCAnalyses::MCParticle::sel_pdgID(421, true)(Particle)')
            .Define('n_genD0', 'FCCAnalyses::MCParticle::get_n(genD0)')
            .Define('genD0_px', 'FCCAnalyses::MCParticle::get_px(genD0)')
            .Define('genD0_py', 'FCCAnalyses::MCParticle::get_py(genD0)')
            .Define('genD0_pz', 'FCCAnalyses::MCParticle::get_pz(genD0)')

            .Define('genDp', 'FCCAnalyses::MCParticle::sel_pdgID(411, true)(Particle)')
            .Define('n_genDp', 'FCCAnalyses::MCParticle::get_n(genDp)')
            .Define('genDp_px', 'FCCAnalyses::MCParticle::get_px(genDp)')
            .Define('genDp_py', 'FCCAnalyses::MCParticle::get_py(genDp)')
            .Define('genDp_pz', 'FCCAnalyses::MCParticle::get_pz(genDp)')

            .Define('genDs', 'FCCAnalyses::MCParticle::sel_pdgID(431, true)(Particle)')
            .Define('n_genDs', 'FCCAnalyses::MCParticle::get_n(genDs)')
            .Define('genDs_px', 'FCCAnalyses::MCParticle::get_px(genDs)')
            .Define('genDs_py', 'FCCAnalyses::MCParticle::get_py(genDs)')
            .Define('genDs_pz', 'FCCAnalyses::MCParticle::get_pz(genDs)')

            .Define('genLc', 'FCCAnalyses::MCParticle::sel_pdgID(4122, true)(Particle)')
            .Define('n_genLc', 'FCCAnalyses::MCParticle::get_n(genLc)')
            .Define('genLc_px', 'FCCAnalyses::MCParticle::get_px(genLc)')
            .Define('genLc_py', 'FCCAnalyses::MCParticle::get_py(genLc)')
            .Define('genLc_pz', 'FCCAnalyses::MCParticle::get_pz(genLc)')

            #############################################
            ## MC vertex object + MC PV
            #############################################
            .Define('MCVertexObject', 'myUtils::get_MCVertexObject(Particle, Particle0)')
            .Define('MC_Vertex_x',    'myUtils::get_MCVertex_x(MCVertexObject)')
            .Define('MC_Vertex_y',    'myUtils::get_MCVertex_y(MCVertexObject)')
            .Define('MC_Vertex_z',    'myUtils::get_MCVertex_z(MCVertexObject)')
            .Define('MC_Vertex_PDG',         'myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)')
            .Define('MC_Vertex_PDGmother',   'myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)')
            .Define('MC_Vertex_PDGgmother',  'myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)')
            ## MC primary vertex. FCCAnalyses' stock
            ## `MCParticle::get_EventPrimaryVertexP4()` reads
            ## `MCParticles[0].vertex`, which on DELPHI-native EDM4hep is the
            ## incoming-beam particle (generatorStatus==21) sitting at (0,0,0).
            ## Use a status==1 production vertex instead — every status==1
            ## particle in an event shares the same production point.
            .Define('MC_PV_xyzt',     'FCCAnalyses::ZHfunctions::get_MC_PV_status1(Particle)')

            #############################################
            ## RECO vertex object — schema-conditional.
            ##  Legacy:  refit via myUtils::get_VertexObject (uses MC-Reco
            ##           link to label vertices).
            ##  New:     build RVec<FCCAnalysesVertex> directly from
            ##           file-side PrimaryVertex + SecondaryVertices via
            ##           ZHfunctions::buildFCCAnalysesVertexFromFileSide.
            ##           No refit, no MC-Reco needed; the converter already
            ##           wrote fitted positions, chi2, and the particles
            ##           relation (_*_particles -> PandoraPFOs).
            ##           Also fills the fitter-output RVec fields (track
            ##           momentum/parameters/phase/chi2) sized to reco_ind
            ##           so myUtils::get_RP_atVertex's `.at(i)` loops
            ##           don't go out of range.
            #############################################
            .Define('VertexObject',
                    ('FCCAnalyses::ZHfunctions::buildFCCAnalysesVertexFromFileSide('
                     'PrimaryVertex, SecondaryVertices, '
                     '_PrimaryVertex_particles, _SecondaryVertices_particles, '
                     'ReconstructedParticles)'
                     if new_schema else
                     'myUtils::get_VertexObject(MCVertexObject, '
                     'ReconstructedParticles, _EFlowTrack_trackStates, '
                     'MCRecoAssociations0, MCRecoAssociations1)'))

            #############################################
            ## PV bookkeeping
            #############################################
            ## hasPV(VertexObject) is the legacy fitter-based check; on
            ## new-schema files the VertexObject ends up empty because the
            ## MC-Reco associations were stubbed (the converter dropped
            ## them), so we OR in the file-side PrimaryVertex test that
            ## both schemas carry.
            .Define('EVT_hasPV',
                    'int(FCCAnalyses::ZHfunctions::hasPV(VertexObject) != 0 '
                    '    || PrimaryVertex.position.x.size() > 0)')
            .Define('EVT_NtracksPV',  'float(myUtils::get_PV_ntracks(VertexObject))')
            .Define('EVT_NVertex',    'float(VertexObject.size())')
            .Filter('EVT_hasPV==1')

            #############################################
            ## Missing energy
            #############################################
            .Define('missingEnergy',
                    'FCCAnalyses::ZHfunctions::missingEnergy(91.188, ReconstructedParticles)')
            .Define('recoEmiss_px', 'missingEnergy[0].momentum.x')
            .Define('recoEmiss_py', 'missingEnergy[0].momentum.y')
            .Define('recoEmiss_pz', 'missingEnergy[0].momentum.z')
            .Define('recoEmiss_e',  'missingEnergy[0].energy')

            #############################################
            ## Reco P with PID, evaluated at the assigned vertex
            #############################################
            .Define('RecoPartPID',
                    'myUtils::PID(ReconstructedParticles, MCRecoAssociations0, MCRecoAssociations1, Particle)')
            .Define('RecoPartPIDAtVertex',
                    'myUtils::get_RP_atVertex(RecoPartPID, VertexObject)')

            #############################################
            ## Vertex-level features
            #############################################
            .Define('Vertex_x',       'myUtils::get_Vertex_x(VertexObject)')
            .Define('Vertex_y',       'myUtils::get_Vertex_y(VertexObject)')
            .Define('Vertex_z',       'myUtils::get_Vertex_z(VertexObject)')
            .Define('Vertex_xErr',    'myUtils::get_Vertex_xErr(VertexObject)')
            .Define('Vertex_yErr',    'myUtils::get_Vertex_yErr(VertexObject)')
            .Define('Vertex_zErr',    'myUtils::get_Vertex_zErr(VertexObject)')
            .Define('Vertex_chi2',    'myUtils::get_Vertex_chi2(VertexObject)')
            .Define('Vertex_isPV',    'FCCAnalyses::ZHfunctions::get_Vertex_isPV(VertexObject)')
            # Per-vertex quality flags (expose-don't-cut, same pattern as
            # RP_passLvlock). isInDet=1 means the SV sits inside the
            # b/c-decay envelope (|d2PV|<50mm, |z|<150mm) — that's
            # where genuine heavy-flavor SVs live. isV0=1 means the
            # 2-track SV has K0_S or Λ-consistent mass — strangeness /
            # baryon signature. Both feed the model as soft features.
            .Define('Vertex_isInDet',
                    'FCCAnalyses::ZHfunctions::get_Vertex_isInDet(VertexObject, 50.0, 150.0)')
            .Define('Vertex_ntrk',    'myUtils::get_Vertex_ntracks(VertexObject)')
            .Define('Vertex_n',       'int(Vertex_x.size())')
            .Define('Vertex_mass',    'myUtils::get_Vertex_mass(VertexObject, RecoPartPIDAtVertex)')
            .Define('Vertex_isV0',
                    'FCCAnalyses::ZHfunctions::get_Vertex_isV0(VertexObject, Vertex_mass)')

            .Define('Vertex_d2PV',    'myUtils::get_Vertex_d2PV(VertexObject, -1)')
            .Define('Vertex_d2PVx',   'myUtils::get_Vertex_d2PV(VertexObject, 0)')
            .Define('Vertex_d2PVy',   'myUtils::get_Vertex_d2PV(VertexObject, 1)')
            .Define('Vertex_d2PVz',   'myUtils::get_Vertex_d2PV(VertexObject, 2)')
            .Define('Vertex_d2PVErr', 'myUtils::get_Vertex_d2PVError(VertexObject, -1)')
            .Define('Vertex_d2PVSig', 'Vertex_d2PV / Vertex_d2PVErr')

            #############################################
            ## Vertex 4-momentum from constituent RPs
            ## (`VertexObject.reco_ind` is filled by `myUtils::get_VertexObject`).
            ## Exposed as per-Vertex flat columns so downstream taggers can
            ## consume Vertex_px/py/pz/e the same way as RP_px/py/pz/e.
            #############################################
            .Define('Vertex_p4',
                    'FCCAnalyses::ZHfunctions::get_Vertex_p4(VertexObject, RecoPartPIDAtVertex)')
            .Define('Vertex_px',    'FCCAnalyses::ZHfunctions::get_p4_px(Vertex_p4)')
            .Define('Vertex_py',    'FCCAnalyses::ZHfunctions::get_p4_py(Vertex_p4)')
            .Define('Vertex_pz',    'FCCAnalyses::ZHfunctions::get_p4_pz(Vertex_p4)')
            .Define('Vertex_e',     'FCCAnalyses::ZHfunctions::get_p4_e(Vertex_p4)')
            .Define('Vertex_phi',   'FCCAnalyses::ZHfunctions::get_p4_phi(Vertex_p4)')
            .Define('Vertex_theta', 'FCCAnalyses::ZHfunctions::get_p4_theta(Vertex_p4)')

            #############################################
            ## Reco P kinematics + per-track helix params
            ##  (DROPPED: RP_dndx — DELPHI-native EDM4hep has no
            ##   EFlowTrack_dNdx collection)
            #############################################
            .Define('RP_e',       'ReconstructedParticle::get_e(RecoPartPIDAtVertex)')
            .Define('RP_px',      'ReconstructedParticle::get_px(RecoPartPIDAtVertex)')
            .Define('RP_py',      'ReconstructedParticle::get_py(RecoPartPIDAtVertex)')
            .Define('RP_pz',      'ReconstructedParticle::get_pz(RecoPartPIDAtVertex)')
            .Define('RP_eta',     'ReconstructedParticle::get_eta(RecoPartPIDAtVertex)')
            .Define('RP_phi',     'ReconstructedParticle::get_phi(RecoPartPIDAtVertex)')
            .Define('RP_theta',   'ReconstructedParticle::get_theta(RecoPartPIDAtVertex)')
            .Define('RP_charge',  'ReconstructedParticle::get_charge(RecoPartPIDAtVertex)')
            .Define('RP_m_true',  'ReconstructedParticle::get_mass(RecoPartPIDAtVertex)')
            .Define('RP_m_reco',  'ReconstructedParticle::get_mass(ReconstructedParticles)')
            .Define('RP_n',       'int(RP_e.size())')
            .Define('RP_fromPV',  'FCCAnalyses::ZHfunctions::get_RP_isfromPV(VertexObject, RecoPartPIDAtVertex)')
            .Define('RP_vert_ind','FCCAnalyses::ZHfunctions::get_RP_Vert_Ind(VertexObject, RecoPartPIDAtVertex)')

            .Define('RP_trk_d0',
                    'ReconstructedParticle2Track::getRP2TRK_D0(RecoPartPIDAtVertex, _EFlowTrack_trackStates)')
            .Define('RP_trk_z0',
                    'ReconstructedParticle2Track::getRP2TRK_Z0(RecoPartPIDAtVertex, _EFlowTrack_trackStates)')
            .Define('RP_trk_phi',
                    'ReconstructedParticle2Track::getRP2TRK_phi(RecoPartPIDAtVertex, _EFlowTrack_trackStates)')
            .Define('RP_trk_omega',
                    'ReconstructedParticle2Track::getRP2TRK_omega(RecoPartPIDAtVertex, _EFlowTrack_trackStates)')
            .Define('RP_trk_tanLambda',
                    'ReconstructedParticle2Track::getRP2TRK_tanLambda(RecoPartPIDAtVertex, _EFlowTrack_trackStates)')

            ## dE/dx per RP via the RecDqdx collection (Tracks_dNdx in
            ## DELPHI-native mode, aliased to EFlowTrack_dNdx above).
            .Define('RP_dndx',
                    'FCCAnalyses::ZHfunctions::get_RP_dndx(RecoPartPIDAtVertex, EFlowTrack_dNdx, _EFlowTrack_dNdx_track)')

            ## Per-PFO PID flags from the *_particle index relations.
            ## Resolved against the ORIGINAL ReconstructedParticles ordering
            ## (PandoraPFOs); RecoPartPIDAtVertex is a permutation that
            ## doesn't drop entries, so hasPIDLink_onRP can re-emit the flag
            ## in RecoPartPIDAtVertex order using a stable PFO index.
            .Define('RP_isMu',     'FCCAnalyses::ZHfunctions::hasPIDLink_onRP(ReconstructedParticles, RecoPartPIDAtVertex, _ParticleID_Muon_particle)')
            .Define('RP_isEl',     'FCCAnalyses::ZHfunctions::hasPIDLink_onRP(ReconstructedParticles, RecoPartPIDAtVertex, _ParticleID_Electron_particle)')
            .Define('RP_hasRich',  'FCCAnalyses::ZHfunctions::hasPIDLink_onRP(ReconstructedParticles, RecoPartPIDAtVertex, _ParticleID_HadronRich_particle)')

            ## SKELANA LVLOCK quality word per PFO. Same ordering as
            ## ReconstructedParticles (= PandoraPFOs in the input). The
            ## helper permutes it into RecoPartPIDAtVertex order so this
            ## column is parallel to RP_e/RP_px/etc.
            ##   -1 = input lacks the field (e.g. pre-LVLOCK nanoaod)
            ##    0 = passes legacy IFLSTR=11/IFLCUT=3 selection
            ##   >0 = locked. Apply `RP_lvlock == 0` for the standard cut
            ##         OR use the pre-computed `RP_passLvlock` flag (0/1).
            ## NOTE on bit-32 (REMCLU) values: SKELANA writes lvlock as a
            ## signed int32 bitmask, so the bit-32-only state shows up as
            ## INT32_MIN (-2147483648). Use `RP_passLvlock == 1` rather
            ## than `RP_lvlock <= 0` to capture the REMCLU class.
            .Define('RP_lvlock',
                    'FCCAnalyses::ZHfunctions::permuteIntOnRP(ReconstructedParticles, RecoPartPIDAtVertex, Track_lvlock)')
            .Define('RP_passLvlock',
                    'FCCAnalyses::ZHfunctions::permuteIntOnRP(ReconstructedParticles, RecoPartPIDAtVertex, RP_passLvlock_input)')

            .Define('RP_nMC',     'FCCAnalyses::ZHfunctions::getRP2MC_nMC(MCRecoAssociations0, MCRecoAssociations1, RecoPartPIDAtVertex)')
            # RP_MCidx: index into MCParticles per RP. The legacy
            # ReconstructedParticle2MC::getRP2MC_index relies on
            # MCRecoAssociations, which is dropped in new-schema files.
            # The angle-based matcher in functions.h reconstructs the link
            # by matching every RP to its closest stable MC particle in
            # (theta, phi) with a loose |p| consistency cut. On legacy
            # files MCRecoAssociations is populated, so the matcher's
            # output supersedes it harmlessly.
            .Define('RP_MCidx',
                    'FCCAnalyses::ZHfunctions::matchRPtoMCByAngle('
                    'ReconstructedParticles, Particle, 0.02, 0.20)')
            # Per-RP truth PDG (0 if no match). Useful for K/pi diagnostics
            # and for sanity-checking the per-hemisphere truth label
            # assignment in prep_bhadron.
            .Define('RP_truthPDG',
                    'FCCAnalyses::ZHfunctions::getRPMatchedPDG(RP_MCidx, Particle)')

            # ---- New-schema PID parameters per RP ------------------------
            # ParticleID_dEdx: 2 params per PID [value, sigma]
            #   - RP_dndx_clean = TPC dE/dx mean (replaces the broken -9 sentinel
            #                     that legacy get_RP_dndx returns on new schema)
            #   - RP_dndx_sigma = uncertainty
            # ParticleID_HadronRich: 18 params per PID, indices for the
            # RICH gas / liquid Cherenkov tags by particle hypothesis.
            # We expose the most useful ones; the full vector is too noisy
            # for direct training but the gas-pi/K/p triplet is the actual
            # PID discriminator the analyzers used historically.
            .Alias('PIDdEdx_params',          '_ParticleID_dEdx_parameters')
            .Alias('PIDdEdx_particle_idx',    '_ParticleID_dEdx_particle.index')
            .Alias('PIDRich_params',          '_ParticleID_HadronRich_parameters')
            .Alias('PIDRich_particle_idx',    '_ParticleID_HadronRich_particle.index')

            # PA.MTPC parameter layout (empirically verified):
            #   param 0 = number of TPC samples in the track (~7 typical)
            #   param 1 = mean dE/dx value -- THE actual ionization measurement
            # The converter README's "[value, sigma]" doc is incorrect.
            .Define('RP_dndx_nSamp',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_dEdx,'
                    ' PIDdEdx_params, PIDdEdx_particle_idx, 0)')
            .Define('RP_dndx_clean',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_dEdx,'
                    ' PIDdEdx_params, PIDdEdx_particle_idx, 1)')

            # ParticleID_HadronRich parameter layout (per ild/...
            # delphi_sdst_to_edm4hep.cpp):
            #   idx 0..2 : KHAID(4..6)  — gas RICH per-hypothesis tag codes
            #   idx 3..4 : KHAID(2,3)   — combined hadronic ID tags
            #   idx 5..6 : QHAID(7,8)   — liquid RICH per-hypothesis tags
            #   idx 7    : KHAID(9)     — auxiliary HAID tag
            #   idx 8..11: THEG, SIGG, NPHG, NEPG  — gas Cherenkov angle + n_photons
            #   idx 12+  : liquid Cherenkov angle + n_photons (mirror of 8..11)
            # The tag values are integer codes (not likelihoods) but the model
            # learns them fine as numerical inputs. We expose the discriminating
            # subset; the angle + n-photon measurements are the most
            # information-dense (continuous, momentum-aware).
            .Define('RP_HAID_gas0',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 0)')
            .Define('RP_HAID_gas1',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 1)')
            .Define('RP_HAID_gas2',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 2)')
            .Define('RP_HAID_liq0',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 5)')
            .Define('RP_HAID_liq1',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 6)')
            .Define('RP_RICHgas_theta',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 8)')
            .Define('RP_RICHgas_nphot',
                    'FCCAnalyses::ZHfunctions::getRPPIDParam('
                    'ReconstructedParticles, ParticleID_HadronRich,'
                    ' PIDRich_params, PIDRich_particle_idx, 10)')
            .Define('RP_fromBs',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(531,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromBu',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(521,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromBd',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(511,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromBc',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(541,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromLb',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(5122, true)(RP_MCidx, Particle, Particle1)')

            ## Inclusive C-hadron descendant tag — any RP that is a stable
            ## descendant of D⁰(421) / D⁺(411) / D_s(431) / Λc(4122).
            .Define('RP_fromD',
                    'FCCAnalyses::ZHfunctions::get_RP_isDescendantAny(std::vector<int>{421,411,431,4122}, true)(RP_MCidx, Particle, Particle1)')
            ## Convenience aggregate: any heavy-flavor (b OR c OR Λb).
            .Define('RP_fromHF',
                    'ROOT::VecOps::RVec<int>(((RP_fromBs+RP_fromBu+RP_fromBd+RP_fromBc+RP_fromLb+RP_fromD) > 0))')

            .Define('Vertex_fromBs', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBs)')
            .Define('Vertex_fromBu', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBu)')
            .Define('Vertex_fromBd', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBd)')
            .Define('Vertex_fromBc', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBc)')
            .Define('Vertex_fromLb', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromLb)')
            .Define('Vertex_fromD',  'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromD)')

            #############################################
            ## Thrust + hemisphere split
            #############################################
            .Define('EVT_thrustNP',
                    'Algorithms::minimize_thrust("Minuit2","Migrad")(RP_px, RP_py, RP_pz)')
            .Define('RP_thrustangleNP',
                    'Algorithms::getAxisCosTheta(EVT_thrustNP, RP_px, RP_py, RP_pz)')
            .Define('EVT_thrust',
                    'Algorithms::getThrustPointing(1.)(RP_thrustangleNP, RP_e, EVT_thrustNP)')
            .Define('RP_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, RP_px, RP_py, RP_pz)')
            .Define('EVT_thrust_phi',   'FCCAnalyses::ZHfunctions::getAxisPhi(EVT_thrust)')
            .Define('EVT_thrust_theta', 'FCCAnalyses::ZHfunctions::getAxisTheta(EVT_thrust)')

            #############################################
            ## Tagger-ready angular deltas + per-RP linked-vertex attrs.
            ## Depends on EVT_thrust_{phi,theta} + Vertex_p4 + RP_vert_ind.
            ## Δφ wrapped to (-π, π]; Δθ unwrapped (theta ∈ [0,π]).
            ## Vertex_thrustangle = cos(angle between Vertex momentum and
            ## thrust axis), parallel to RP_thrustangle.
            ## RP_vert_{e,mass} = per-RP lookup of the assigned vertex's
            ## energy/mass; -1 for RPs not assigned to any vertex.
            #############################################
            .Define('RP_Dphi',       'FCCAnalyses::ZHfunctions::delta_phi(RP_phi, EVT_thrust_phi)')
            .Define('RP_Dtheta',     'FCCAnalyses::ZHfunctions::delta_theta(RP_theta, EVT_thrust_theta)')
            .Define('Vertex_Dphi',   'FCCAnalyses::ZHfunctions::delta_phi(Vertex_phi, EVT_thrust_phi)')
            .Define('Vertex_Dtheta', 'FCCAnalyses::ZHfunctions::delta_theta(Vertex_theta, EVT_thrust_theta)')
            .Define('Vertex_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, Vertex_px, Vertex_py, Vertex_pz)')
            .Define('RP_vert_e',
                    'FCCAnalyses::ZHfunctions::get_RP_vert_attr(RP_vert_ind, Vertex_e)')
            .Define('RP_vert_mass',
                    'FCCAnalyses::ZHfunctions::get_RP_vert_attr(RP_vert_ind, Vertex_mass)')

            # Per-RP "is-V0-daughter" proxy: if the RP belongs to a vertex
            # whose mass is consistent with K0_S or Λ, flag it. This is the
            # stage1-level surrogate for SKELANA's LVLOCK V0-tag bit
            # (which the new-schema converter doesn't yet emit). Approximate
            # but the model can learn it as a soft strangeness/V0 feature.
            .Define('RP_vert_isV0_int',
                    'ROOT::VecOps::RVec<float>(Vertex_isV0.begin(), Vertex_isV0.end())')
            .Define('RP_isV0daughter_f',
                    'FCCAnalyses::ZHfunctions::get_RP_vert_attr(RP_vert_ind, RP_vert_isV0_int)')
            .Define('RP_isV0daughter',
                    'ROOT::VecOps::RVec<int>(RP_isV0daughter_f.begin(), RP_isV0daughter_f.end())')

            .Define('EVT_thrusthemis0_n', 'Algorithms::getAxisN(0)(RP_thrustangle, RP_charge)')
            .Define('EVT_thrusthemis1_n', 'Algorithms::getAxisN(1)(RP_thrustangle, RP_charge)')
            .Define('EVT_thrusthemis0_e', 'Algorithms::getAxisEnergy(0)(RP_thrustangle, RP_charge, RP_e)')
            .Define('EVT_thrusthemis1_e', 'Algorithms::getAxisEnergy(1)(RP_thrustangle, RP_charge, RP_e)')

            .Define('EVT_ThrustEmax_E',         'EVT_thrusthemis0_e.at(0)')
            .Define('EVT_ThrustEmax_Echarged',  'EVT_thrusthemis0_e.at(1)')
            .Define('EVT_ThrustEmax_Eneutral',  'EVT_thrusthemis0_e.at(2)')
            .Define('EVT_ThrustEmax_N',         'float(EVT_thrusthemis0_n.at(0))')
            .Define('EVT_ThrustEmax_Ncharged',  'float(EVT_thrusthemis0_n.at(1))')
            .Define('EVT_ThrustEmax_Nneutral',  'float(EVT_thrusthemis0_n.at(2))')

            .Define('EVT_ThrustEmin_E',         'EVT_thrusthemis1_e.at(0)')
            .Define('EVT_ThrustEmin_Echarged',  'EVT_thrusthemis1_e.at(1)')
            .Define('EVT_ThrustEmin_Eneutral',  'EVT_thrusthemis1_e.at(2)')
            .Define('EVT_ThrustEmin_N',         'float(EVT_thrusthemis1_n.at(0))')
            .Define('EVT_ThrustEmin_Ncharged',  'float(EVT_thrusthemis1_n.at(1))')
            .Define('EVT_ThrustEmin_Nneutral',  'float(EVT_thrusthemis1_n.at(2))')

            .Define('EVT_Thrust_Mag', 'EVT_thrust.at(0)')

            ###############################################################
            ## Durham k_T exclusive clustering -> EXACTLY 2 jets (dijet-like
            ## event topology). The paper arXiv:2603.06524 applies this same
            ## clustering and then cuts on:
            ##    y_3 < 0.13  (well-defined 2-jet event, no hard gluon)
            ##    jet |p| > 10 GeV
            ##    |cos theta_jet| < 0.65
            ## We expose all the needed quantities; the downstream prep
            ## step (prep_partonflavour.py --mode jet) applies the cuts and
            ## emits one row per jet.
            ##
            ## clustering_ee_kt(arg_exclusive=2, arg_cut=2, arg_sorted=1,
            ##                  arg_recombination=0):
            ##   - exclusive=2: cluster to EXACTLY arg_cut jets
            ##   - cut=2: N=2 jets
            ##   - sorted=1: order jets by energy (descending)
            ##   - recombination=0: E-scheme (standard for Durham)
            ##
            ## get_exclusive_dmerge(jet, 2) returns d_{2,3} - the merge
            ## distance at the 3->2 transition. Divide by E_cm^2 (91.2 GeV)
            ## to get the paper's y_3 normalization.
            ###############################################################
            .Define('PseudoJets_RP',
                    'JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)')
            .Define('FCCAnalysesJets_durham',
                    'JetClustering::clustering_ee_kt(2, 2, 1, 0)(PseudoJets_RP)')
            .Define('Jets_durham',
                    'JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_durham)')
            .Define('JetConstituents_durham',
                    'JetClusteringUtils::get_constituents(FCCAnalysesJets_durham)')

            ## y_3 = d_{2,3} / E_cm^2 (normalized as in the ALEPH paper).
            ## 91.2^2 = 8317.44.
            .Define('EVT_y3_durham',
                    'JetClusteringUtils::get_exclusive_dmerge(FCCAnalysesJets_durham, 2) / 8317.44f')

            ## Per-jet kinematics (jets are E-ordered: jet[0] = leading,
            ## jet[1] = subleading).
            .Define('Jets_e',     'JetClusteringUtils::get_e(Jets_durham)')
            .Define('Jets_p',     'JetClusteringUtils::get_p(Jets_durham)')
            .Define('Jets_pt',    'JetClusteringUtils::get_pt(Jets_durham)')
            .Define('Jets_theta', 'JetClusteringUtils::get_theta(Jets_durham)')
            .Define('Jets_phi',   'JetClusteringUtils::get_phi(Jets_durham)')

            ## Per-RP jet index: 0 = leading jet, 1 = subleading; -1 if not
            ## assigned (should not happen with exclusive clustering). Done
            ## inline via an RVec invert of the constituents list.
            .Define('RP_jet_idx',
                    'ROOT::VecOps::RVec<int> idx(static_cast<size_t>(RP_n), -1);'
                    ' for (size_t j = 0; j < JetConstituents_durham.size(); ++j) {'
                    '   for (auto rp : JetConstituents_durham[j]) {'
                    '     if (rp >= 0 && static_cast<int>(rp) < RP_n) idx[rp] = static_cast<int>(j);'
                    '   }'
                    ' }'
                    ' return idx;')

            .Define('genBs_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBs_px, genBs_py, genBs_pz)')
            .Define('genBu_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBu_px, genBu_py, genBu_pz)')
            .Define('genBd_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBd_px, genBd_py, genBd_pz)')
            .Define('genBc_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBc_px, genBc_py, genBc_pz)')
            .Define('genLb_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genLb_px, genLb_py, genLb_pz)')

            ## c-hadron thrust angles (used by hemisphere label_C_*)
            .Define('genCharm_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genCharm_px, genCharm_py, genCharm_pz)')
            .Define('genD0_thrustangle',    'Algorithms::getAxisCosTheta(EVT_thrust, genD0_px, genD0_py, genD0_pz)')
            .Define('genDp_thrustangle',    'Algorithms::getAxisCosTheta(EVT_thrust, genDp_px, genDp_py, genDp_pz)')
            .Define('genDs_thrustangle',    'Algorithms::getAxisCosTheta(EVT_thrust, genDs_px, genDs_py, genDs_pz)')
            .Define('genLc_thrustangle',    'Algorithms::getAxisCosTheta(EVT_thrust, genLc_px, genLc_py, genLc_pz)')

            #############################################
            ## Per-hemisphere truth labels
            #############################################
            .Define('n_Bs_Emin', 'int(genBs_thrustangle[genBs_thrustangle>0].size())')
            .Define('n_Bu_Emin', 'int(genBu_thrustangle[genBu_thrustangle>0].size())')
            .Define('n_Bd_Emin', 'int(genBd_thrustangle[genBd_thrustangle>0].size())')
            .Define('n_Bc_Emin', 'int(genBc_thrustangle[genBc_thrustangle>0].size())')
            .Define('n_Lb_Emin', 'int(genLb_thrustangle[genLb_thrustangle>0].size())')
            .Define('label_Bs_Emin', 'int(n_Bc_Emin==0 && n_Bs_Emin==1 && n_Bu_Emin==0 && n_Bd_Emin==0 && n_Lb_Emin==0)')
            .Define('label_Bu_Emin', 'int(n_Bc_Emin==0 && n_Bs_Emin==0 && n_Bu_Emin==1 && n_Bd_Emin==0 && n_Lb_Emin==0)')
            .Define('label_Bd_Emin', 'int(n_Bc_Emin==0 && n_Bs_Emin==0 && n_Bu_Emin==0 && n_Bd_Emin==1 && n_Lb_Emin==0)')
            .Define('label_Bc_Emin', 'int(n_Bc_Emin==1 && n_Bs_Emin==0 && n_Bu_Emin==0 && n_Bd_Emin==0 && n_Lb_Emin==0)')
            .Define('label_Lb_Emin', 'int(n_Bc_Emin==0 && n_Bs_Emin==0 && n_Bu_Emin==0 && n_Bd_Emin==0 && n_Lb_Emin==1)')

            ## Heavy-flavor count = any of the five B species (incl. Λb).
            .Define('n_B_Emin',   'int(n_Bs_Emin + n_Bu_Emin + n_Bd_Emin + n_Bc_Emin + n_Lb_Emin)')
            ## Inclusive C count: D⁰ + D⁺ + D_s + Λc with positive thrust cosine.
            .Define('n_D0_Emin', 'int(genD0_thrustangle[genD0_thrustangle>0].size())')
            .Define('n_Dp_Emin', 'int(genDp_thrustangle[genDp_thrustangle>0].size())')
            .Define('n_Ds_Emin', 'int(genDs_thrustangle[genDs_thrustangle>0].size())')
            .Define('n_Lc_Emin', 'int(genLc_thrustangle[genLc_thrustangle>0].size())')
            .Define('n_C_Emin',  'int(n_D0_Emin + n_Dp_Emin + n_Ds_Emin + n_Lc_Emin)')
            ## label_C: ≥1 C-hadron AND no B-hadron in this hemisphere.
            .Define('label_C_Emin',     'int(n_C_Emin >= 1 && n_B_Emin == 0)')
            ## label_udsg: no heavy-flavor hadron of any kind (light hemisphere).
            .Define('label_udsg_Emin',  'int(n_B_Emin == 0 && n_C_Emin == 0)')

            .Define('n_Bs_Emax', 'int(genBs_thrustangle[genBs_thrustangle<0].size())')
            .Define('n_Bu_Emax', 'int(genBu_thrustangle[genBu_thrustangle<0].size())')
            .Define('n_Bd_Emax', 'int(genBd_thrustangle[genBd_thrustangle<0].size())')
            .Define('n_Bc_Emax', 'int(genBc_thrustangle[genBc_thrustangle<0].size())')
            .Define('n_Lb_Emax', 'int(genLb_thrustangle[genLb_thrustangle<0].size())')
            .Define('label_Bs_Emax', 'int(n_Bc_Emax==0 && n_Bs_Emax==1 && n_Bu_Emax==0 && n_Bd_Emax==0 && n_Lb_Emax==0)')
            .Define('label_Bu_Emax', 'int(n_Bc_Emax==0 && n_Bs_Emax==0 && n_Bu_Emax==1 && n_Bd_Emax==0 && n_Lb_Emax==0)')
            .Define('label_Bd_Emax', 'int(n_Bc_Emax==0 && n_Bs_Emax==0 && n_Bu_Emax==0 && n_Bd_Emax==1 && n_Lb_Emax==0)')
            .Define('label_Bc_Emax', 'int(n_Bc_Emax==1 && n_Bs_Emax==0 && n_Bu_Emax==0 && n_Bd_Emax==0 && n_Lb_Emax==0)')
            .Define('label_Lb_Emax', 'int(n_Bc_Emax==0 && n_Bs_Emax==0 && n_Bu_Emax==0 && n_Bd_Emax==0 && n_Lb_Emax==1)')

            .Define('n_B_Emax',  'int(n_Bs_Emax + n_Bu_Emax + n_Bd_Emax + n_Bc_Emax + n_Lb_Emax)')
            .Define('n_D0_Emax', 'int(genD0_thrustangle[genD0_thrustangle<0].size())')
            .Define('n_Dp_Emax', 'int(genDp_thrustangle[genDp_thrustangle<0].size())')
            .Define('n_Ds_Emax', 'int(genDs_thrustangle[genDs_thrustangle<0].size())')
            .Define('n_Lc_Emax', 'int(genLc_thrustangle[genLc_thrustangle<0].size())')
            .Define('n_C_Emax',  'int(n_D0_Emax + n_Dp_Emax + n_Ds_Emax + n_Lc_Emax)')
            .Define('label_C_Emax',     'int(n_C_Emax >= 1 && n_B_Emax == 0)')
            .Define('label_udsg_Emax',  'int(n_B_Emax == 0 && n_C_Emax == 0)')

            #############################################
            ## PV — both reco PV and offset vs. MC truth
            #############################################
            ## FCC re-fit PV (from `myUtils::get_VertexObject`). Useful for
            ## reco-vs-reco closure between the FCC vertexing and the
            ## converter's PrimaryVertex, but it returns coordinates that are
            ## NOT in the same frame as the file PV (BS-prior subtraction +
            ## a different fitter), so don't use it for data/MC PV closure.
            .Define('PVfit_x',       'Vertex_x[Vertex_isPV==1]')
            .Define('PVfit_y',       'Vertex_y[Vertex_isPV==1]')
            .Define('PVfit_z',       'Vertex_z[Vertex_isPV==1]')
            .Define('PV_ntrk',       'Vertex_ntrk[Vertex_isPV==1]')

            ## File-side reco PV — the converter-written PrimaryVertex
            ## collection (one entry per event in our converter). This is
            ## the right PV for data/MC closure on PV offsets.
            .Define('PV_x',          'float(PrimaryVertex.position.x[0])')
            .Define('PV_y',          'float(PrimaryVertex.position.y[0])')
            .Define('PV_z',          'float(PrimaryVertex.position.z[0])')
            .Define('PV_x_offset',   'PV_x - float(MC_PV_xyzt.X())')
            .Define('PV_y_offset',   'PV_y - float(MC_PV_xyzt.Y())')
            .Define('PV_z_offset',   'PV_z - float(MC_PV_xyzt.Z())')

            #############################################
            ## Secondary vertices — exposed from the converter's
            ## SecondaryVertices collection (pair-seed-and-add fits from
            ## delphi-improved-reco/reco/vertex/sv_reco.cpp). Position,
            ## chi²/ndf, and displacement from PV are available. The
            ## constituent track relation `_SecondaryVertices_particles`
            ## is empty on this generation of files, so SV_ntrk == 0
            ## everywhere and SV_mass cannot be computed at the analysis
            ## level — that needs a converter-side fix. `SV_inDet` flags
            ## SVs whose position is plausibly inside the tracker
            ## (|d3D|<5cm, |z|<15cm) — pair-seed sometimes produces
            ## extrapolation outliers at O(1 m).
            #############################################
            .Define('SV_n',     'FCCAnalyses::ZHfunctions::get_SV_n(SecondaryVertices)')
            .Define('SV_x',     'FCCAnalyses::ZHfunctions::get_SV_x(SecondaryVertices)')
            .Define('SV_y',     'FCCAnalyses::ZHfunctions::get_SV_y(SecondaryVertices)')
            .Define('SV_z',     'FCCAnalyses::ZHfunctions::get_SV_z(SecondaryVertices)')
            .Define('SV_chi2',  'FCCAnalyses::ZHfunctions::get_SV_chi2(SecondaryVertices)')
            .Define('SV_ndf',   'FCCAnalyses::ZHfunctions::get_SV_ndf(SecondaryVertices)')
            .Define('SV_ntrk',  'FCCAnalyses::ZHfunctions::get_SV_ntrk(SecondaryVertices)')
            .Define('SV_d2PV',  'FCCAnalyses::ZHfunctions::get_SV_d2PV(SecondaryVertices, PV_x, PV_y, PV_z, -1)')
            .Define('SV_d2PVx', 'FCCAnalyses::ZHfunctions::get_SV_d2PV(SecondaryVertices, PV_x, PV_y, PV_z, 0)')
            .Define('SV_d2PVy', 'FCCAnalyses::ZHfunctions::get_SV_d2PV(SecondaryVertices, PV_x, PV_y, PV_z, 1)')
            .Define('SV_d2PVz', 'FCCAnalyses::ZHfunctions::get_SV_d2PV(SecondaryVertices, PV_x, PV_y, PV_z, 2)')
            .Define('SV_inDet', 'FCCAnalyses::ZHfunctions::get_SV_inDet(SecondaryVertices, PV_x, PV_y, PV_z)')
            .Define('SV_n_inDet', 'int(ROOT::VecOps::Sum(SV_inDet))')
        )
        return df2

    def output():
        return [
            'EVT_NVertex', 'EVT_NtracksPV', 'EVT_hasPV',
            'MC_n', 'n_genBottoms',
            'n_genBs', 'n_genBu', 'n_genBd', 'n_genBc', 'n_genLb',
            'genBs_px', 'genBs_py', 'genBs_pz',
            'genBu_px', 'genBu_py', 'genBu_pz',
            'genBd_px', 'genBd_py', 'genBd_pz',
            'genBc_px', 'genBc_py', 'genBc_pz',
            'genLb_px', 'genLb_py', 'genLb_pz',

            'n_genCharms', 'n_genD0', 'n_genDp', 'n_genDs', 'n_genLc',
            'genCharm_px', 'genCharm_py', 'genCharm_pz', 'genCharm_pdg',
            'genD0_px', 'genD0_py', 'genD0_pz',
            'genDp_px', 'genDp_py', 'genDp_pz',
            'genDs_px', 'genDs_py', 'genDs_pz',
            'genLc_px', 'genLc_py', 'genLc_pz',

            'recoEmiss_px', 'recoEmiss_py', 'recoEmiss_pz', 'recoEmiss_e',

            'RP_n', 'RP_e', 'RP_m_true', 'RP_m_reco',
            'RP_px', 'RP_py', 'RP_pz', 'RP_phi', 'RP_theta', 'RP_charge',
            'RP_thrustangle', 'RP_Dphi', 'RP_Dtheta',
            'RP_fromPV', 'RP_vert_ind', 'RP_vert_e', 'RP_vert_mass',
            'RP_isV0daughter',
            'RP_trk_d0', 'RP_trk_z0', 'RP_trk_phi', 'RP_trk_omega', 'RP_trk_tanLambda',
            'RP_dndx', 'RP_isMu', 'RP_isEl', 'RP_hasRich', 'RP_lvlock', 'RP_passLvlock',
            'RP_nMC', 'RP_MCidx', 'RP_truthPDG',
            'RP_dndx_clean', 'RP_dndx_nSamp',
            'RP_HAID_gas0', 'RP_HAID_gas1', 'RP_HAID_gas2',
            'RP_HAID_liq0', 'RP_HAID_liq1',
            'RP_RICHgas_theta', 'RP_RICHgas_nphot',
            'RP_fromBs', 'RP_fromBu', 'RP_fromBd', 'RP_fromBc', 'RP_fromLb',
            'RP_fromD', 'RP_fromHF',

            'Vertex_n', 'Vertex_x', 'Vertex_y', 'Vertex_z',
            'Vertex_xErr', 'Vertex_yErr', 'Vertex_zErr', 'Vertex_chi2',
            'Vertex_isPV', 'Vertex_isInDet', 'Vertex_isV0',
            'Vertex_ntrk', 'Vertex_mass',
            'Vertex_px', 'Vertex_py', 'Vertex_pz', 'Vertex_e',
            'Vertex_phi', 'Vertex_theta',
            'Vertex_Dphi', 'Vertex_Dtheta', 'Vertex_thrustangle',
            'Vertex_d2PV', 'Vertex_d2PVx', 'Vertex_d2PVy', 'Vertex_d2PVz',
            'Vertex_d2PVErr', 'Vertex_d2PVSig',
            'Vertex_fromBs', 'Vertex_fromBu', 'Vertex_fromBd', 'Vertex_fromBc', 'Vertex_fromLb',
            'Vertex_fromD',

            'EVT_Thrust_Mag', 'EVT_thrust_phi', 'EVT_thrust_theta',
            'EVT_ThrustEmin_E', 'EVT_ThrustEmin_Echarged', 'EVT_ThrustEmin_Eneutral',
            'EVT_ThrustEmin_N', 'EVT_ThrustEmin_Ncharged', 'EVT_ThrustEmin_Nneutral',
            'EVT_ThrustEmax_E', 'EVT_ThrustEmax_Echarged', 'EVT_ThrustEmax_Eneutral',
            'EVT_ThrustEmax_N', 'EVT_ThrustEmax_Ncharged', 'EVT_ThrustEmax_Nneutral',

            ## Durham k_T exclusive clustering (2 jets), for jet-level training.
            'EVT_y3_durham',
            'Jets_e', 'Jets_p', 'Jets_pt', 'Jets_theta', 'Jets_phi',
            'RP_jet_idx',

            'genBs_thrustangle', 'genBu_thrustangle', 'genBd_thrustangle',
            'genBc_thrustangle', 'genLb_thrustangle',
            'genCharm_thrustangle', 'genD0_thrustangle', 'genDp_thrustangle',
            'genDs_thrustangle', 'genLc_thrustangle',
            'n_Bs_Emin', 'n_Bu_Emin', 'n_Bd_Emin', 'n_Bc_Emin', 'n_Lb_Emin',
            'n_B_Emin', 'n_D0_Emin', 'n_Dp_Emin', 'n_Ds_Emin', 'n_Lc_Emin', 'n_C_Emin',
            'label_Bs_Emin', 'label_Bu_Emin', 'label_Bd_Emin', 'label_Bc_Emin', 'label_Lb_Emin',
            'label_C_Emin', 'label_udsg_Emin',
            'n_Bs_Emax', 'n_Bu_Emax', 'n_Bd_Emax', 'n_Bc_Emax', 'n_Lb_Emax',
            'n_B_Emax', 'n_D0_Emax', 'n_Dp_Emax', 'n_Ds_Emax', 'n_Lc_Emax', 'n_C_Emax',
            'label_Bs_Emax', 'label_Bu_Emax', 'label_Bd_Emax', 'label_Bc_Emax', 'label_Lb_Emax',
            'label_C_Emax', 'label_udsg_Emax',

            'PV_x', 'PV_y', 'PV_z', 'PV_ntrk',
            'PV_x_offset', 'PV_y_offset', 'PV_z_offset',
            'PVfit_x', 'PVfit_y', 'PVfit_z',

            'SV_n', 'SV_n_inDet',
            'SV_x', 'SV_y', 'SV_z',
            'SV_chi2', 'SV_ndf', 'SV_ntrk',
            'SV_d2PV', 'SV_d2PVx', 'SV_d2PVy', 'SV_d2PVz',
            'SV_inDet',
        ]
