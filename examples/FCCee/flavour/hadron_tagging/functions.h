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

#include "FCCAnalyses/MCParticle.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "FCCAnalyses/myUtils.h"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/RecDqdxData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"

namespace FCCAnalyses {
namespace ZHfunctions {

using rp = edm4hep::ReconstructedParticleData;
using Vec_rp = ROOT::VecOps::RVec<rp>;

// missing energy (vs. an assumed CM-frame total energy).
inline Vec_rp missingEnergy(float ecm, Vec_rp in, float p_cutoff = 0.0) {
  float px = 0, py = 0, pz = 0, e = 0;
  for (auto &p : in) {
    if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y * p.momentum.y) <
        p_cutoff)
      continue;
    px += -p.momentum.x;
    py += -p.momentum.y;
    pz += -p.momentum.z;
    e += p.energy;
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

inline float getAxisPhi(const ROOT::VecOps::RVec<float> &axis) {
  TLorentzVector tlv;
  tlv.SetXYZM(axis[1], axis[3], axis[5], 0);
  return tlv.Phi();
}

inline float getAxisTheta(const ROOT::VecOps::RVec<float> &axis) {
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
  for (auto &p : vertex)
    result.push_back(p.mc_ind == 0 ? 1 : 0);
  return result;
}

inline int hasPV(ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex) {
  for (auto &p : vertex)
    if (p.mc_ind == 0)
      return 1;
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
    auto &p = vertex[iv];
    ROOT::VecOps::RVec<int> reco_ind = p.reco_ind;
    for (size_t ip = 0; ip < reco_ind.size(); ++ip) {
      result[reco_ind.at(ip)] = iv;
    }
  }
  return result;
}

inline ROOT::VecOps::RVec<int> getRP2MC_nMC(ROOT::VecOps::RVec<int> recind,
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
  int m_pdg = 13;
  bool m_chargeconjugate = true;

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
      const bool match = m_chargeconjugate
                             ? (std::abs(p.PDG) == std::abs(m_pdg))
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

inline ROOT::VecOps::RVec<int> get_Vertex_containDescendant(
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vertex,
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

// Per-PFO PID flag from a `_ParticleID_*_particle` index relation.
// Resolves the index against the ORIGINAL ReconstructedParticles ordering
// (the column the relation was written against), then re-emits the flag
// in the order of `rpAtVertex` (a permutation of the same particles).
//
// We rely on stable per-RP identity: every RP in rpAtVertex appears
// exactly once in rp, and we match by the byte-identical
// ReconstructedParticleData payload (cheap: track_begin is unique per RP
// in our converter).
inline ROOT::VecOps::RVec<int>
hasPIDLink_onRP(Vec_rp rp, Vec_rp rpAtVertex,
                ROOT::VecOps::RVec<podio::ObjectID> pid_to_part) {
  std::vector<bool> hasPID(rp.size(), false);
  for (auto &oid : pid_to_part) {
    if (oid.index >= 0 && static_cast<size_t>(oid.index) < rp.size())
      hasPID[oid.index] = true;
  }
  ROOT::VecOps::RVec<int> result;
  result.reserve(rpAtVertex.size());
  for (auto &p : rpAtVertex) {
    bool flag = false;
    for (size_t j = 0; j < rp.size(); ++j) {
      if (rp[j].tracks_begin == p.tracks_begin &&
          rp[j].tracks_end == p.tracks_end && rp[j].charge == p.charge &&
          rp[j].energy == p.energy) {
        flag = hasPID[j];
        break;
      }
    }
    result.push_back(flag ? 1 : 0);
  }
  return result;
}

// Filter a ReconstructedParticles collection by SKELANA's per-track
// LVLOCK bitmask (`Track_lvlock` UserData column, parallel to the
// PandoraPFOs/Tracks ordering on the file side). The legacy SKELANA
// convention is `pass <=> LVLOCK == 0` (all 32 bits zero).
//
//   lvlock == 0          all bits clear — keep (legacy selection)
//   lvlock == 1          bit 1: LVSELE quality fail — drop
//   lvlock == INT32_MIN  bit 32: REMCLU calo-cluster overlap — drop
//                                (energy already merged into a charged
//                                 track PFO; keeping these double-counts
//                                 calorimeter energy in thrust / EEC /
//                                 missing-pt sums)
//   lvlock == INT32_MIN+1 bit 32 + bit 1                       — drop
//   lvlock == -1         input nanoaod predates the field — keep
//                                (degrade gracefully on legacy files)
//
// Implementation: keep when `lvlock == 0 || lvlock == -1`. Any other
// value is a locked Part and MUST be dropped — even when "<= 0" looks
// inclusive (the INT32_MIN bit-32 class is < 0 numerically but is the
// most important reject class for calorimeter energy bookkeeping).
//
// CAVEAT: this function changes the size and indexing of the returned
// collection. `MCRecoAssociations_to.index` from the writer points into
// the *original* PandoraPFOs ordering, so passing the filtered output
// to any helper that consumes `MCRecoAssociations*` (notably
// `myUtils::get_VertexObject`, `myUtils::PID`, `getRP2MC_index`) reads
// past the end and segfaults. For stage1 use `lvlockPassMask` to get a
// per-RP `int` flag parallel to PandoraPFOs and apply the selection at
// the analysis output level, not before vertexing.
inline Vec_rp filterRPbyLvlock(Vec_rp rp, ROOT::VecOps::RVec<int> lvlock) {
  Vec_rp out;
  out.reserve(rp.size());
  for (size_t i = 0; i < rp.size(); ++i) {
    int lv = (i < lvlock.size()) ? lvlock[i] : -1;
    if (lv == 0 || lv == -1)
      out.push_back(rp[i]); // strict zero (or unknown) keeps
  }
  return out;
}

// Per-RP `passes-LVLOCK` flag, parallel to `rp` (= PandoraPFOs on file).
// 1 = keep (lvlock == 0 or -1), 0 = drop. Same semantics as
// `filterRPbyLvlock` but as a mask — preserves the original RP indexing
// so MCRecoAssociations / Vertex helpers stay valid.
inline ROOT::VecOps::RVec<int>
lvlockPassMask(Vec_rp rp, ROOT::VecOps::RVec<int> lvlock) {
  ROOT::VecOps::RVec<int> out;
  out.reserve(rp.size());
  for (size_t i = 0; i < rp.size(); ++i) {
    int lv = (i < lvlock.size()) ? lvlock[i] : -1;
    out.push_back((lv == 0 || lv == -1) ? 1 : 0);
  }
  return out;
}

// Full DELPHI-legacy track selection at the FCCAnalyser level, matching
// Jingyu's `apply_track_selection_delphi` in delphi-analysis/python.
// Applied to the *charged* RP subset (our converter is charged-only PFO):
//
//   sel_c = (lvlock <= 0)                                  # quality
//         & (charge != 0)                                  # charged
//         & (p > pT_min)                                   # momentum
//         & (|D0| < d0_max_mm)                             # transverse IP
//         & (|Z0 * sin θ| < z0_sin_max_mm)                 # longitudinal IP
//         & (theta_min_rad < θ < theta_max_rad)            # polar accept
//
// D0 / Z0 are read off the parallel `_Tracks_trackStates` column (mm).
// The PandoraPFOs.tracks_begin field indexes into the Tracks collection;
// since our converter emits one TrackState per Track at AtIP, the
// trackState index equals the track index, which equals the PFO index
// (1:1 charged-RP : Track : trackState mapping in delphi_to_edm4hep).
// Fall back to charge/p/θ only if the trackState lookup fails (no IP cut).
//
// Units: D0/Z0 in mm (EDM4hep convention), p in GeV. Defaults reproduce
// Jingyu's nominal cuts (pT>0.4, |d0|<4cm=40mm, |z0sinθ|<4cm=40mm,
// 20°<θ<160°).
inline Vec_rp
filterRP_delphi_legacy(Vec_rp rp, ROOT::VecOps::RVec<int> lvlock,
                       ROOT::VecOps::RVec<edm4hep::TrackState> trackStates,
                       float pT_min = 0.4f, float d0_max_mm = 40.0f,
                       float z0_sin_max_mm = 40.0f,
                       float theta_min_rad = 0.349065850f, // 20° = π/9
                       float theta_max_rad = 2.792526803f) // 160°
{
  Vec_rp out;
  out.reserve(rp.size());
  for (size_t i = 0; i < rp.size(); ++i) {
    int lv = (i < lvlock.size()) ? lvlock[i] : -1;
    // SKELANA LVLOCK convention: keep iff all 32 bits zero (or no info).
    // `lv > 0` alone misses the bit-32 REMCLU class (INT32_MIN < 0).
    if (!(lv == 0 || lv == -1))
      continue;
    const auto &p = rp[i];
    if (std::abs(p.charge) < 0.1f)
      continue; // require charged
    const float px = p.momentum.x;
    const float py = p.momentum.y;
    const float pz = p.momentum.z;
    const float pT2 = px * px + py * py;
    const float pMag = std::sqrt(pT2 + pz * pz);
    if (pMag < pT_min)
      continue;
    const float theta = std::atan2(std::sqrt(pT2), pz);
    if (theta < theta_min_rad || theta > theta_max_rad)
      continue;
    // IP cuts via trackState. In delphi_to_edm4hep the i-th PFO has
    // tracks_begin = i (1:1 PFO↔Track) and one TrackState at AtIP at
    // trackState index = i, so direct lookup is fine.
    if (i < trackStates.size()) {
      const auto &ts = trackStates[i];
      const float d0 = ts.D0;
      const float z0sin = ts.Z0 * std::sin(theta);
      if (std::abs(d0) > d0_max_mm)
        continue;
      if (std::abs(z0sin) > z0_sin_max_mm)
        continue;
    } else {
      continue; // no trackState → can't apply IP cuts, drop conservatively
    }
    out.push_back(p);
  }
  return out;
}

// Permute a per-PFO int32 UserData column (e.g., Track_lvlock) from the
// original ReconstructedParticles ordering into the RecoPartPIDAtVertex
// ordering used elsewhere in the analysis. Returns one entry per
// rpAtVertex element. Falls back to -1 if no matching original RP found.
inline ROOT::VecOps::RVec<int> permuteIntOnRP(Vec_rp rp, Vec_rp rpAtVertex,
                                              ROOT::VecOps::RVec<int> flat) {
  ROOT::VecOps::RVec<int> result;
  result.reserve(rpAtVertex.size());
  for (auto &p : rpAtVertex) {
    int value = -1;
    for (size_t j = 0; j < rp.size() && j < flat.size(); ++j) {
      if (rp[j].tracks_begin == p.tracks_begin &&
          rp[j].tracks_end == p.tracks_end && rp[j].charge == p.charge &&
          rp[j].energy == p.energy) {
        value = flat[j];
        break;
      }
    }
    result.push_back(value);
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
get_RP_dndx(Vec_rp in, ROOT::VecOps::RVec<edm4hep::RecDqdxData> dNdx,
            ROOT::VecOps::RVec<podio::ObjectID> dNdx_track) {
  std::vector<int> trk_to_dndx(8192, -1);
  for (size_t k = 0; k < dNdx_track.size() && k < dNdx.size(); ++k) {
    const int trk_idx = dNdx_track[k].index;
    if (trk_idx < 0)
      continue;
    if (static_cast<size_t>(trk_idx) >= trk_to_dndx.size())
      trk_to_dndx.resize(trk_idx + 1, -1);
    trk_to_dndx[trk_idx] = static_cast<int>(k);
  }
  ROOT::VecOps::RVec<float> result;
  result.reserve(in.size());
  for (auto &p : in) {
    if (p.charge == 0 ||
        p.tracks_begin >= static_cast<unsigned>(trk_to_dndx.size()) ||
        trk_to_dndx[p.tracks_begin] < 0) {
      result.push_back(-9.f);
      continue;
    }
    const int k = trk_to_dndx[p.tracks_begin];
    result.push_back(dNdx[k].dQdx.value);
  }
  return result;
}

} // namespace ZHfunctions
} // namespace FCCAnalyses

#endif
