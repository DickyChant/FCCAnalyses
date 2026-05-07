#ifndef ZHfunctions_H
#define ZHfunctions_H

// Restored helpers + DELPHI EDM4hep 1.0 schema adapters for hadron_tagging.
//
// Modern key4hep / EDM4hep 1.0 update of the centos7-stack helpers from
// hadronTagger_dev/functions.h. Removed parts that depended on EDM4hep 0.x
// field names that no longer exist (vertex.primary, TrackData.dxQuantities_*,
// TrackerHitData, RP::type) — those features should be rebuilt from the
// modern schema (Vertex.algorithmType, RecDqdxCollection, TrackerHit3DData,
// RP::PDG) once they're actually needed by the analysis. For the smoke
// test we only need the helpers below.

#include <cmath>
#include <cstdint>
#include <vector>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/RecDqdxData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "FCCAnalyses/MCParticle.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "FCCAnalyses/myUtils.h"

namespace FCCAnalyses { namespace ZHfunctions {

using rp = edm4hep::ReconstructedParticleData;
using Vec_rp = ROOT::VecOps::RVec<rp>;

// missing energy (vs. an assumed CM-frame total energy).
inline Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
    float px = 0, py = 0, pz = 0, e = 0;
    for (auto &p : in) {
        if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y) < p_cutoff) continue;
        px += -p.momentum.x;
        py += -p.momentum.y;
        pz += -p.momentum.z;
        e  += p.energy;
    }
    Vec_rp ret;
    rp res;
    res.momentum.x = px;
    res.momentum.y = py;
    res.momentum.z = pz;
    res.energy = ecm - e;
    ret.emplace_back(res);
    return ret;
}

inline float getAxisPhi(const ROOT::VecOps::RVec<float> & axis) {
    TLorentzVector tlv;
    tlv.SetXYZM(axis[1], axis[3], axis[5], 0);
    return tlv.Phi();
}

inline float getAxisTheta(const ROOT::VecOps::RVec<float> & axis) {
    TLorentzVector tlv;
    tlv.SetXYZM(axis[1], axis[3], axis[5], 0);
    return tlv.Theta();
}

// FCCAnalysesVertex helpers — restored verbatim from the legacy version.

inline ROOT::VecOps::RVec<TLorentzVector>
build_p4(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py,
         ROOT::VecOps::RVec<float> pz, ROOT::VecOps::RVec<float> mass) {
    ROOT::VecOps::RVec<TLorentzVector> p4;
    for (size_t i = 0; i < px.size(); ++i) {
        TLorentzVector tlv;
        tlv.SetXYZM(px[i], py[i], pz[i], mass[i]);
        p4.push_back(tlv);
    }
    return p4;
}

inline ROOT::VecOps::RVec<TLorentzVector>
get_Vertex_p4(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
              Vec_rp reco) {
    ROOT::VecOps::RVec<TLorentzVector> result;
    for (auto &p : vertex) {
        ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
        TLorentzVector tlv = myUtils::build_tlv(reco, reco_ind);
        result.push_back(tlv);
    }
    return result;
}

// Modern equivalent of the legacy `vertex.primary == 1` check. In EDM4hep 1.0
// `Vertex::primary` is gone; the type is encoded in a bitfield (BITPrimary
// = bit 1). FCCAnalyses' `VertexFitter` writes `vertex.type = Primary`
// literally, which is wrong for the bit convention, so we use mc_ind which is
// robustly set by `get_VertexObject` (mc_ind==0 is the MC primary vertex).
inline ROOT::VecOps::RVec<int>
get_Vertex_isPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
    ROOT::VecOps::RVec<int> result;
    result.reserve(vertex.size());
    for (auto &p : vertex) result.push_back(p.mc_ind == 0 ? 1 : 0);
    return result;
}

inline int hasPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
    for (auto &p : vertex) if (p.mc_ind == 0) return 1;
    return 0;
}

inline ROOT::VecOps::RVec<int>
get_RP_isfromPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                Vec_rp reco) {
    ROOT::VecOps::RVec<int> result;
    result.resize(reco.size(), -1);
    for (auto &p : vertex) {
        ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
        const int isPV = (p.mc_ind == 0) ? 1 : 2;
        for (size_t j = 0; j < reco_ind.size(); ++j)
            result[reco_ind.at(j)] = isPV;
    }
    return result;
}

inline ROOT::VecOps::RVec<int>
get_RP_Vert_Ind(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                Vec_rp reco) {
    ROOT::VecOps::RVec<int> result;
    result.resize(reco.size(), -1);
    for (size_t iv = 0; iv < vertex.size(); ++iv) {
        auto & p = vertex[iv];
        ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
        for (size_t ip = 0; ip < reco_ind.size(); ++ip) {
            result[reco_ind.at(ip)] = iv;
        }
    }
    return result;
}

inline ROOT::VecOps::RVec<int>
getRP2MC_nMC(ROOT::VecOps::RVec<int> recind,
             ROOT::VecOps::RVec<int> /*mcind*/,
             Vec_rp reco) {
    ROOT::VecOps::RVec<int> count;
    count.resize(reco.size(), 0);
    for (size_t i = 0; i < recind.size(); ++i) {
        count[recind.at(i)] += 1;
    }
    return count;
}

// Mark RPs whose truth-matched MCParticle is a stable descendant of any
// MC particle with the given PDG.
struct get_RP_isDescendant {
    int  m_pdg               = 13;
    bool m_chargeconjugate   = true;

    get_RP_isDescendant(int arg_pdg, bool arg_chargeconjugate)
        : m_pdg(arg_pdg), m_chargeconjugate(arg_chargeconjugate) {}

    ROOT::VecOps::RVec<int>
    operator()(ROOT::VecOps::RVec<int> reco_mcidx,
               ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
               ROOT::VecOps::RVec<int> ind) {
        // collect indices of all stable descendants of any matching mother
        ROOT::VecOps::RVec<int> descd;
        for (size_t i = 0; i < in.size(); ++i) {
            const auto &p = in[i];
            const bool match = m_chargeconjugate ? (std::abs(p.PDG) == std::abs(m_pdg))
                                                 : (p.PDG == m_pdg);
            if (match) {
                std::vector<int> rr =
                    MCParticle::get_list_of_stable_particles_from_decay(i, in, ind);
                descd.insert(descd.end(), rr.begin(), rr.end());
            }
        }
        ROOT::VecOps::RVec<int> result;
        result.resize(reco_mcidx.size(), 0);
        for (size_t i = 0; i < reco_mcidx.size(); ++i) {
            if (std::find(descd.begin(), descd.end(), reco_mcidx[i]) != descd.end())
                result[i] = 1;
        }
        return result;
    }
};

inline ROOT::VecOps::RVec<int>
get_Vertex_containDescendant(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
                             ROOT::VecOps::RVec<int> rp_isDescendant) {
    ROOT::VecOps::RVec<int> result;
    for (auto &p : vertex) {
        ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
        int contain = 0;
        for (size_t i = 0; i < reco_ind.size(); ++i) {
            contain += rp_isDescendant[reco_ind.at(i)];
        }
        result.push_back(contain);
    }
    return result;
}

// Per-RP dE/dx pulled from the modern EDM4hep RecDqdxCollection schema.
// dNdx is parallel to `_EFlowTrack_dNdx_track.index`, which gives the index
// into the EFlowTrack collection that this RecDqdx entry refers to. We index
// the dNdx by RP via the RP's tracks_begin (the first track of the RP is
// the canonical one in our DELPHI converter). Returned in DELPHI's native
// dE/dx 80%-truncated-mean units (~1.0 for MIPs, no MeV/cm conversion).
inline ROOT::VecOps::RVec<float>
get_RP_dndx(Vec_rp in,
            ROOT::VecOps::RVec<edm4hep::RecDqdxData> dNdx,
            ROOT::VecOps::RVec<podio::ObjectID> dNdx_track) {
    std::vector<int> trk_to_dndx(8192, -1);
    for (size_t k = 0; k < dNdx_track.size() && k < dNdx.size(); ++k) {
        const int trk_idx = dNdx_track[k].index;
        if (trk_idx < 0) continue;
        if (static_cast<size_t>(trk_idx) >= trk_to_dndx.size())
            trk_to_dndx.resize(trk_idx + 1, -1);
        trk_to_dndx[trk_idx] = static_cast<int>(k);
    }
    ROOT::VecOps::RVec<float> result;
    result.reserve(in.size());
    for (auto &p : in) {
        if (p.charge == 0 || p.tracks_begin >= static_cast<unsigned>(trk_to_dndx.size())
            || trk_to_dndx[p.tracks_begin] < 0) {
            result.push_back(-9.f);
            continue;
        }
        const int k = trk_to_dndx[p.tracks_begin];
        result.push_back(dNdx[k].dQdx.value);
    }
    return result;
}

}}  // namespace

#endif
