'''
hadron_tagging stage1 — modern port (DELPHI EDM4hep 1.0 smoke test).

Adapted from the centos7-stack version on branch `hadronTagger_dev`. Major
differences vs. the legacy version:

- target stack: key4hep modern (FCCAnalyses 0.12.x, EDM4hep 1.0).
- target sample: DELPHI Z->bb EDM4hep produced by `delphi-improved-reco/ild/delphi_to_edm4hep`
  (under `feature/edm4hep-pipeline`). Schema differences from FCC IDEA:
    * no `MCRecoAssociations` collection (no truth-match map)
    * no `EFlowTrack_1`/`EFlowTrack_2` (modern schema; track state is in
      `_EFlowTrack_trackStates`, RecDqdx is its own collection `EFlowTrack_dNdx`)
    * `Particle` and `ReconstructedParticles` are present
    * `PrimaryVertex` is a single-entry Vertex collection
- features kept: gen-level B-hadron tagging labels, reco P4, thrust hemispheres
- features dropped (depend on missing collections): `RP_fromBs`/`RP_MCidx`,
  `Vertex_*`, `RP_trk_*`, `RP_dndx`, `RP_mtof`, BDT inference

This is a smoke test that the toolchain runs end-to-end on our DELPHI EDM4hep
files. Restoring the full feature set requires adding `MCRecoAssociations` and
the other missing branches to `delphi_to_edm4hep`.
'''

# Mandatory: List of processes
processList = {
    'events_phA_bb_fcc': {'output': 'fcc_smoke_zbb'},
}

# Mandatory: Production tag — set to None and use inputDir below.
prodTag = None

# Mandatory: Input directory (local DELPHI EDM4hep output).
inputDir = '/eos/user/s/sqian/www/delphi_edm4hep/phaseA-fcc/zbb'

# Mandatory: output directory
outputDir = 'outputs/FCCee/flavour/hadron_tagging/stage1'

# Optional: include local helpers
includePaths = ['functions.h']

# Optional: number of threads
nCPUS = 4

# Mandatory: RDFanalysis class
class RDFanalysis():

    def analysers(df):
        df2 = (
            df
            #############################################
            ## relations alias (Particle parents/daughters)
            #############################################
            .Alias('Particle0', '_Particle_parents.index')
            .Alias('Particle1', '_Particle_daughters.index')

            #############################################
            ## MC particles bookkeeping
            #############################################
            .Define('MC_PDG', 'FCCAnalyses::MCParticle::get_pdg(Particle)')
            .Define('MC_n', 'int(MC_PDG.size())')
            .Define('MC_M1', 'myUtils::get_MCMother1(Particle, Particle0)')
            .Define('MC_M2', 'myUtils::get_MCMother2(Particle, Particle0)')
            .Define('MC_D1', 'myUtils::get_MCDaughter1(Particle, Particle1)')
            .Define('MC_D2', 'myUtils::get_MCDaughter2(Particle, Particle1)')

            #############################################
            ## Gen B-hadron / b-quark counters
            #############################################
            .Define('genBottom', 'FCCAnalyses::MCParticle::sel_pdgID(5, true)(Particle)')
            .Define('n_genBottoms', 'FCCAnalyses::MCParticle::get_n(genBottom)')
            .Define('genBottom_px', 'FCCAnalyses::MCParticle::get_px(genBottom)')
            .Define('genBottom_py', 'FCCAnalyses::MCParticle::get_py(genBottom)')
            .Define('genBottom_pz', 'FCCAnalyses::MCParticle::get_pz(genBottom)')
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
            ## True PV (xyzt) from MCParticles
            #############################################
            .Define('MC_PV_xyzt', 'FCCAnalyses::MCParticle::get_EventPrimaryVertexP4()(Particle)')

            #############################################
            ## Reco kinematics — directly off ReconstructedParticles
            #############################################
            .Define('RP_e', 'ReconstructedParticle::get_e(ReconstructedParticles)')
            .Define('RP_m_reco', 'ReconstructedParticle::get_mass(ReconstructedParticles)')
            .Define('RP_px', 'ReconstructedParticle::get_px(ReconstructedParticles)')
            .Define('RP_py', 'ReconstructedParticle::get_py(ReconstructedParticles)')
            .Define('RP_pz', 'ReconstructedParticle::get_pz(ReconstructedParticles)')
            .Define('RP_eta', 'ReconstructedParticle::get_eta(ReconstructedParticles)')
            .Define('RP_phi', 'ReconstructedParticle::get_phi(ReconstructedParticles)')
            .Define('RP_theta', 'ReconstructedParticle::get_theta(ReconstructedParticles)')
            .Define('RP_charge', 'ReconstructedParticle::get_charge(ReconstructedParticles)')
            .Define('RP_n', 'int(RP_e.size())')

            #############################################
            ## Missing energy (vs. CM-frame)
            #############################################
            .Define('missingEnergy',
                    'FCCAnalyses::ZHfunctions::missingEnergy(91.188, ReconstructedParticles)')
            .Define('recoEmiss_px', 'missingEnergy[0].momentum.x')
            .Define('recoEmiss_py', 'missingEnergy[0].momentum.y')
            .Define('recoEmiss_pz', 'missingEnergy[0].momentum.z')
            .Define('recoEmiss_e', 'missingEnergy[0].energy')

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
            .Define('EVT_thrust_phi',
                    'FCCAnalyses::ZHfunctions::getAxisPhi(EVT_thrust)')
            .Define('EVT_thrust_theta',
                    'FCCAnalyses::ZHfunctions::getAxisTheta(EVT_thrust)')

            .Define('EVT_thrusthemis0_n',
                    'Algorithms::getAxisN(0)(RP_thrustangle, RP_charge)')
            .Define('EVT_thrusthemis1_n',
                    'Algorithms::getAxisN(1)(RP_thrustangle, RP_charge)')
            .Define('EVT_thrusthemis0_e',
                    'Algorithms::getAxisEnergy(0)(RP_thrustangle, RP_charge, RP_e)')
            .Define('EVT_thrusthemis1_e',
                    'Algorithms::getAxisEnergy(1)(RP_thrustangle, RP_charge, RP_e)')

            .Define('EVT_ThrustEmax_E', 'EVT_thrusthemis0_e.at(0)')
            .Define('EVT_ThrustEmax_Echarged', 'EVT_thrusthemis0_e.at(1)')
            .Define('EVT_ThrustEmax_Eneutral', 'EVT_thrusthemis0_e.at(2)')
            .Define('EVT_ThrustEmax_N', 'float(EVT_thrusthemis0_n.at(0))')
            .Define('EVT_ThrustEmax_Ncharged', 'float(EVT_thrusthemis0_n.at(1))')
            .Define('EVT_ThrustEmax_Nneutral', 'float(EVT_thrusthemis0_n.at(2))')

            .Define('EVT_ThrustEmin_E', 'EVT_thrusthemis1_e.at(0)')
            .Define('EVT_ThrustEmin_Echarged', 'EVT_thrusthemis1_e.at(1)')
            .Define('EVT_ThrustEmin_Eneutral', 'EVT_thrusthemis1_e.at(2)')
            .Define('EVT_ThrustEmin_N', 'float(EVT_thrusthemis1_n.at(0))')
            .Define('EVT_ThrustEmin_Ncharged', 'float(EVT_thrusthemis1_n.at(1))')
            .Define('EVT_ThrustEmin_Nneutral', 'float(EVT_thrusthemis1_n.at(2))')

            .Define('EVT_Thrust_Mag', 'EVT_thrust.at(0)')
            .Define('EVT_Thrust_X', 'EVT_thrust.at(1)')
            .Define('EVT_Thrust_Y', 'EVT_thrust.at(3)')
            .Define('EVT_Thrust_Z', 'EVT_thrust.at(5)')

            #############################################
            ## Gen B thrust angles + per-hemisphere counts
            #############################################
            .Define('genBs_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, genBs_px, genBs_py, genBs_pz)')
            .Define('genBu_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, genBu_px, genBu_py, genBu_pz)')
            .Define('genBd_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, genBd_px, genBd_py, genBd_pz)')
            .Define('genBc_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, genBc_px, genBc_py, genBc_pz)')
            .Define('genLb_thrustangle',
                    'Algorithms::getAxisCosTheta(EVT_thrust, genLb_px, genLb_py, genLb_pz)')

            .Define('n_Bs_Emin', 'int(genBs_thrustangle[genBs_thrustangle>0].size())')
            .Define('n_Bu_Emin', 'int(genBu_thrustangle[genBu_thrustangle>0].size())')
            .Define('n_Bd_Emin', 'int(genBd_thrustangle[genBd_thrustangle>0].size())')
            .Define('n_Bc_Emin', 'int(genBc_thrustangle[genBc_thrustangle>0].size())')
            .Define('n_Lb_Emin', 'int(genLb_thrustangle[genLb_thrustangle>0].size())')

            .Define('label_Bs_Emin',
                    'int(n_Bc_Emin==0 && n_Bs_Emin==1 && n_Bu_Emin==0 && n_Bd_Emin==0 && n_Lb_Emin==0)')
            .Define('label_Bu_Emin',
                    'int(n_Bc_Emin==0 && n_Bs_Emin==0 && n_Bu_Emin==1 && n_Bd_Emin==0 && n_Lb_Emin==0)')
            .Define('label_Bd_Emin',
                    'int(n_Bc_Emin==0 && n_Bs_Emin==0 && n_Bu_Emin==0 && n_Bd_Emin==1 && n_Lb_Emin==0)')
            .Define('label_Bc_Emin',
                    'int(n_Bc_Emin==1 && n_Bs_Emin==0 && n_Bu_Emin==0 && n_Bd_Emin==0 && n_Lb_Emin==0)')
            .Define('label_Lb_Emin',
                    'int(n_Bc_Emin==0 && n_Bs_Emin==0 && n_Bu_Emin==0 && n_Bd_Emin==0 && n_Lb_Emin==1)')

            .Define('n_Bs_Emax', 'int(genBs_thrustangle[genBs_thrustangle<0].size())')
            .Define('n_Bu_Emax', 'int(genBu_thrustangle[genBu_thrustangle<0].size())')
            .Define('n_Bd_Emax', 'int(genBd_thrustangle[genBd_thrustangle<0].size())')
            .Define('n_Bc_Emax', 'int(genBc_thrustangle[genBc_thrustangle<0].size())')
            .Define('n_Lb_Emax', 'int(genLb_thrustangle[genLb_thrustangle<0].size())')

            .Define('label_Bs_Emax',
                    'int(n_Bc_Emax==0 && n_Bs_Emax==1 && n_Bu_Emax==0 && n_Bd_Emax==0 && n_Lb_Emax==0)')
            .Define('label_Bu_Emax',
                    'int(n_Bc_Emax==0 && n_Bs_Emax==0 && n_Bu_Emax==1 && n_Bd_Emax==0 && n_Lb_Emax==0)')
            .Define('label_Bd_Emax',
                    'int(n_Bc_Emax==0 && n_Bs_Emax==0 && n_Bu_Emax==0 && n_Bd_Emax==1 && n_Lb_Emax==0)')
            .Define('label_Bc_Emax',
                    'int(n_Bc_Emax==1 && n_Bs_Emax==0 && n_Bu_Emax==0 && n_Bd_Emax==0 && n_Lb_Emax==0)')
            .Define('label_Lb_Emax',
                    'int(n_Bc_Emax==0 && n_Bs_Emax==0 && n_Bu_Emax==0 && n_Bd_Emax==0 && n_Lb_Emax==1)')

            #############################################
            ## Reco PV (single-vertex collection)
            #############################################
            .Define('PV_x', 'PrimaryVertex.position.x')
            .Define('PV_y', 'PrimaryVertex.position.y')
            .Define('PV_z', 'PrimaryVertex.position.z')
            .Define('PV_chi2', 'PrimaryVertex.chi2')
            .Define('PV_x_offset', 'PV_x - MC_PV_xyzt.X()')
            .Define('PV_y_offset', 'PV_y - MC_PV_xyzt.Y()')
            .Define('PV_z_offset', 'PV_z - MC_PV_xyzt.Z()')
        )
        return df2

    def output():
        return [
            'MC_n', 'n_genBottoms',
            'n_genBs', 'n_genBu', 'n_genBd', 'n_genBc', 'n_genLb',
            'genBs_px', 'genBs_py', 'genBs_pz',
            'genBu_px', 'genBu_py', 'genBu_pz',
            'genBd_px', 'genBd_py', 'genBd_pz',
            'genBc_px', 'genBc_py', 'genBc_pz',
            'genLb_px', 'genLb_py', 'genLb_pz',
            'RP_n',
            'RP_e', 'RP_px', 'RP_py', 'RP_pz', 'RP_charge',
            'RP_thrustangle',
            'recoEmiss_px', 'recoEmiss_py', 'recoEmiss_pz', 'recoEmiss_e',
            'EVT_Thrust_Mag', 'EVT_Thrust_X', 'EVT_Thrust_Y', 'EVT_Thrust_Z',
            'EVT_thrust_phi', 'EVT_thrust_theta',
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
            'PV_x', 'PV_y', 'PV_z', 'PV_chi2',
            'PV_x_offset', 'PV_y_offset', 'PV_z_offset',
        ]
