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
        df2 = (
            df
            #############################################
            ## Collection-name aliases: DELPHI-native -> FCC names used
            ## throughout the FCCAnalyses helpers / Defines below.
            #############################################
            .Alias('Particle',                 'MCParticles')
            .Alias('_EFlowTrack_trackStates',  '_Tracks_trackStates')

            ## DELPHI-legacy track-quality cut: SKELANA `LVLOCK == 0`
            ## (with -1 allowed as the legacy-no-info value). The writer
            ## emits every PFO so the cut lives at the FCCAnalyser entry
            ## point. IMPORTANT: do NOT shrink the RP collection before
            ## vertexing — `MCRecoAssociations_to.index` from the writer
            ## points into the ORIGINAL PandoraPFOs ordering, and the
            ## `get_VertexObject` / `myUtils::PID` / `getRP2MC_index`
            ## helpers will segfault if we re-index out from under those
            ## relations. Instead, alias the unfiltered collection, expose
            ## the lvlock-pass flag as a parallel mask, and use it as a
            ## per-RP cut downstream (or as the BDT input cut later).
            .Alias('ReconstructedParticles',   'PandoraPFOs')
            .Define('RP_passLvlock_input',
                    'FCCAnalyses::ZHfunctions::lvlockPassMask(PandoraPFOs, Track_lvlock)')
            .Alias('EFlowTrack_dNdx',          'Tracks_dNdx')
            .Alias('_EFlowTrack_dNdx_track',   '_Tracks_dNdx_track')

            #############################################
            ## Index-relation aliases (file-side branch names)
            #############################################
            .Alias('Particle0',           '_MCParticles_parents.index')
            .Alias('Particle1',           '_MCParticles_daughters.index')
            .Alias('MCRecoAssociations0', '_MCRecoAssociations_from.index')
            .Alias('MCRecoAssociations1', '_MCRecoAssociations_to.index')

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
            ## MC vertex object + MC PV
            #############################################
            .Define('MCVertexObject', 'myUtils::get_MCVertexObject(Particle, Particle0)')
            .Define('MC_Vertex_x',    'myUtils::get_MCVertex_x(MCVertexObject)')
            .Define('MC_Vertex_y',    'myUtils::get_MCVertex_y(MCVertexObject)')
            .Define('MC_Vertex_z',    'myUtils::get_MCVertex_z(MCVertexObject)')
            .Define('MC_Vertex_PDG',         'myUtils::get_MCpdgMCVertex(MCVertexObject, Particle)')
            .Define('MC_Vertex_PDGmother',   'myUtils::get_MCpdgMotherMCVertex(MCVertexObject, Particle)')
            .Define('MC_Vertex_PDGgmother',  'myUtils::get_MCpdgGMotherMCVertex(MCVertexObject, Particle)')
            .Define('MC_PV_xyzt', 'FCCAnalyses::MCParticle::get_EventPrimaryVertexP4()(Particle)')

            #############################################
            ## RECO vertex object
            #############################################
            .Define('VertexObject',
                    'myUtils::get_VertexObject(MCVertexObject, ReconstructedParticles, _EFlowTrack_trackStates, MCRecoAssociations0, MCRecoAssociations1)')

            #############################################
            ## PV bookkeeping
            #############################################
            .Define('EVT_hasPV',      'FCCAnalyses::ZHfunctions::hasPV(VertexObject)')
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
            .Define('Vertex_ntrk',    'myUtils::get_Vertex_ntracks(VertexObject)')
            .Define('Vertex_n',       'int(Vertex_x.size())')
            .Define('Vertex_mass',    'myUtils::get_Vertex_mass(VertexObject, RecoPartPIDAtVertex)')

            .Define('Vertex_d2PV',    'myUtils::get_Vertex_d2PV(VertexObject, -1)')
            .Define('Vertex_d2PVx',   'myUtils::get_Vertex_d2PV(VertexObject, 0)')
            .Define('Vertex_d2PVy',   'myUtils::get_Vertex_d2PV(VertexObject, 1)')
            .Define('Vertex_d2PVz',   'myUtils::get_Vertex_d2PV(VertexObject, 2)')
            .Define('Vertex_d2PVErr', 'myUtils::get_Vertex_d2PVError(VertexObject, -1)')
            .Define('Vertex_d2PVSig', 'Vertex_d2PV / Vertex_d2PVErr')

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
            .Define('RP_MCidx',   'ReconstructedParticle2MC::getRP2MC_index(MCRecoAssociations0, MCRecoAssociations1, RecoPartPIDAtVertex)')
            .Define('RP_fromBs',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(531,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromBu',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(521,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromBd',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(511,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromBc',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(541,  true)(RP_MCidx, Particle, Particle1)')
            .Define('RP_fromLb',  'FCCAnalyses::ZHfunctions::get_RP_isDescendant(5122, true)(RP_MCidx, Particle, Particle1)')

            .Define('Vertex_fromBs', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBs)')
            .Define('Vertex_fromBu', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBu)')
            .Define('Vertex_fromBd', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBd)')
            .Define('Vertex_fromBc', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromBc)')
            .Define('Vertex_fromLb', 'FCCAnalyses::ZHfunctions::get_Vertex_containDescendant(VertexObject, RP_fromLb)')

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

            .Define('genBs_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBs_px, genBs_py, genBs_pz)')
            .Define('genBu_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBu_px, genBu_py, genBu_pz)')
            .Define('genBd_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBd_px, genBd_py, genBd_pz)')
            .Define('genBc_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genBc_px, genBc_py, genBc_pz)')
            .Define('genLb_thrustangle', 'Algorithms::getAxisCosTheta(EVT_thrust, genLb_px, genLb_py, genLb_pz)')

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

            #############################################
            ## PV — both reco PV and offset vs. MC truth
            #############################################
            .Define('PV_x',          'Vertex_x[Vertex_isPV==1]')
            .Define('PV_y',          'Vertex_y[Vertex_isPV==1]')
            .Define('PV_z',          'Vertex_z[Vertex_isPV==1]')
            .Define('PV_ntrk',       'Vertex_ntrk[Vertex_isPV==1]')
            .Define('PV_x_offset',   'PV_x - MC_PV_xyzt.X()')
            .Define('PV_y_offset',   'PV_y - MC_PV_xyzt.Y()')
            .Define('PV_z_offset',   'PV_z - MC_PV_xyzt.Z()')
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

            'recoEmiss_px', 'recoEmiss_py', 'recoEmiss_pz', 'recoEmiss_e',

            'RP_n', 'RP_e', 'RP_m_reco',
            'RP_px', 'RP_py', 'RP_pz', 'RP_phi', 'RP_theta', 'RP_charge',
            'RP_thrustangle', 'RP_fromPV', 'RP_vert_ind',
            'RP_trk_d0', 'RP_trk_z0', 'RP_trk_phi', 'RP_trk_omega', 'RP_trk_tanLambda',
            'RP_dndx', 'RP_isMu', 'RP_isEl', 'RP_hasRich', 'RP_lvlock', 'RP_passLvlock',
            'RP_nMC', 'RP_MCidx',
            'RP_fromBs', 'RP_fromBu', 'RP_fromBd', 'RP_fromBc', 'RP_fromLb',

            'Vertex_n', 'Vertex_x', 'Vertex_y', 'Vertex_z',
            'Vertex_xErr', 'Vertex_yErr', 'Vertex_zErr', 'Vertex_chi2',
            'Vertex_isPV', 'Vertex_ntrk', 'Vertex_mass',
            'Vertex_d2PV', 'Vertex_d2PVx', 'Vertex_d2PVy', 'Vertex_d2PVz',
            'Vertex_d2PVErr', 'Vertex_d2PVSig',
            'Vertex_fromBs', 'Vertex_fromBu', 'Vertex_fromBd', 'Vertex_fromBc', 'Vertex_fromLb',

            'EVT_Thrust_Mag', 'EVT_thrust_phi', 'EVT_thrust_theta',
            'EVT_ThrustEmin_E', 'EVT_ThrustEmin_Echarged', 'EVT_ThrustEmin_Eneutral',
            'EVT_ThrustEmin_N', 'EVT_ThrustEmin_Ncharged', 'EVT_ThrustEmin_Nneutral',
            'EVT_ThrustEmax_E', 'EVT_ThrustEmax_Echarged', 'EVT_ThrustEmax_Eneutral',
            'EVT_ThrustEmax_N', 'EVT_ThrustEmax_Ncharged', 'EVT_ThrustEmax_Nneutral',

            'genBs_thrustangle', 'genBu_thrustangle', 'genBd_thrustangle',
            'genBc_thrustangle', 'genLb_thrustangle',
            'n_Bs_Emin', 'n_Bu_Emin', 'n_Bd_Emin', 'n_Bc_Emin', 'n_Lb_Emin',
            'label_Bs_Emin', 'label_Bu_Emin', 'label_Bd_Emin', 'label_Bc_Emin', 'label_Lb_Emin',
            'n_Bs_Emax', 'n_Bu_Emax', 'n_Bd_Emax', 'n_Bc_Emax', 'n_Lb_Emax',
            'label_Bs_Emax', 'label_Bu_Emax', 'label_Bd_Emax', 'label_Bc_Emax', 'label_Lb_Emax',

            'PV_x', 'PV_y', 'PV_z', 'PV_ntrk',
            'PV_x_offset', 'PV_y_offset', 'PV_z_offset',
        ]
