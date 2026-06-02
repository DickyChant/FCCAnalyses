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
#include "FCCAnalyses/VertexingUtils.h"
#include "FCCAnalyses/myUtils.h"
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/RecDqdxData.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/VertexData.h"
#include "podio/ObjectID.h"

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

// True MC primary vertex from the MCParticles tree. The FCCAnalyses
// stock helper `MCParticle::get_EventPrimaryVertexP4()` reads
// `MCParticles[0].vertex`, but on DELPHI-native EDM4hep that slot is
// the incoming beam particle (generatorStatus == 21) at the
// conventional origin (0,0,0). The real MC PV lives in the production
// vertex of any final-state (status == 1) particle — all status==1
// particles in an event share the same production point — so we walk
// to the first one and read its `vertex` field.
inline TLorentzVector
get_MC_PV_status1(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  TLorentzVector tlv;
  for (auto &p : in) {
    if (p.generatorStatus == 1) {
      tlv.SetXYZT(p.vertex.x, p.vertex.y, p.vertex.z, p.time);
      return tlv;
    }
  }
  // Fallback: MCParticles[0].vertex (the FCCAnalyses default behaviour).
  if (!in.empty())
    tlv.SetXYZT(in[0].vertex.x, in[0].vertex.y, in[0].vertex.z, in[0].time);
  return tlv;
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

// Same as get_RP_isDescendant but accepts a list of PDGs and marks RPs
// whose truth-matched MCParticle descends from ANY of them. Used for
// inclusive c-hadron tagging (D⁰/D⁺/D_s/Λc) and for the "no heavy
// flavor" light-flavor label.
struct get_RP_isDescendantAny {
  std::vector<int> m_pdgs;
  bool m_chargeconjugate = true;

  get_RP_isDescendantAny(std::vector<int> arg_pdgs, bool arg_cc)
      : m_pdgs(std::move(arg_pdgs)), m_chargeconjugate(arg_cc) {}

  ROOT::VecOps::RVec<int>
  operator()(ROOT::VecOps::RVec<int> reco_mcidx,
             ROOT::VecOps::RVec<edm4hep::MCParticleData> in,
             ROOT::VecOps::RVec<int> ind) {
    ROOT::VecOps::RVec<int> descd;
    for (size_t i = 0; i < in.size(); ++i) {
      const auto &p = in[i];
      bool match = false;
      for (int q : m_pdgs) {
        const bool m =
            m_chargeconjugate ? (std::abs(p.PDG) == std::abs(q)) : (p.PDG == q);
        if (m) {
          match = true;
          break;
        }
      }
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
inline ROOT::VecOps::RVec<int> lvlockPassMask(Vec_rp rp,
                                              ROOT::VecOps::RVec<int> lvlock) {
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

// =================== Secondary Vertex accessors ====================
// Modern EDM4hep 1.0 stores secondary vertices in their own collection
// (`SecondaryVertices`) as `edm4hep::Vertex`. `delphi_to_edm4hep` writes
// position + chi² + ndf for every SV found by the pair-seed-and-add SV
// finder, but the constituent track relation
// (`_SecondaryVertices_particles`) is empty on this generation of files
// (`particles_end == particles_begin`). So SV_mass and per-track SV
// linkage are not derivable on the analysis side — they need a
// converter-side fix. Everything else (position, chi², displacement
// from PV, quality flag) is here.

inline int get_SV_n(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  return static_cast<int>(sv.size());
}

inline ROOT::VecOps::RVec<float>
get_SV_x(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(sv.size());
  for (auto &v : sv)
    out.push_back(v.position.x);
  return out;
}
inline ROOT::VecOps::RVec<float>
get_SV_y(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(sv.size());
  for (auto &v : sv)
    out.push_back(v.position.y);
  return out;
}
inline ROOT::VecOps::RVec<float>
get_SV_z(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(sv.size());
  for (auto &v : sv)
    out.push_back(v.position.z);
  return out;
}
inline ROOT::VecOps::RVec<float>
get_SV_chi2(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(sv.size());
  for (auto &v : sv)
    out.push_back(v.chi2);
  return out;
}
inline ROOT::VecOps::RVec<int>
get_SV_ndf(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  ROOT::VecOps::RVec<int> out;
  out.reserve(sv.size());
  for (auto &v : sv)
    out.push_back(v.ndf);
  return out;
}
// Number of constituent tracks: (particles_end - particles_begin).
// Will be 0 on inputs where the converter didn't fill the relation.
inline ROOT::VecOps::RVec<int>
get_SV_ntrk(ROOT::VecOps::RVec<edm4hep::VertexData> sv) {
  ROOT::VecOps::RVec<int> out;
  out.reserve(sv.size());
  for (auto &v : sv)
    out.push_back(static_cast<int>(v.particles_end) -
                  static_cast<int>(v.particles_begin));
  return out;
}

// Distance from the (file-side) PV. axis: -1 = 3D, 0 = x, 1 = y, 2 = z.
inline ROOT::VecOps::RVec<float>
get_SV_d2PV(ROOT::VecOps::RVec<edm4hep::VertexData> sv, float PVx, float PVy,
            float PVz, int axis = -1) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(sv.size());
  for (auto &v : sv) {
    const float dx = v.position.x - PVx;
    const float dy = v.position.y - PVy;
    const float dz = v.position.z - PVz;
    if (axis == 0)
      out.push_back(dx);
    else if (axis == 1)
      out.push_back(dy);
    else if (axis == 2)
      out.push_back(dz);
    else
      out.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
  }
  return out;
}

// In-detector quality flag: |SV-PV|_3D < max_mm AND |SV.z|<150 mm.
// DELPHI's sv_reco occasionally produces extrapolation outliers at
// O(1 m) — reject them with this flag before downstream cuts.
inline ROOT::VecOps::RVec<int>
get_SV_inDet(ROOT::VecOps::RVec<edm4hep::VertexData> sv, float PVx, float PVy,
             float PVz, float max_d2PV_mm = 50.0f, float max_absz_mm = 150.0f) {
  ROOT::VecOps::RVec<int> out;
  out.reserve(sv.size());
  for (auto &v : sv) {
    const float dx = v.position.x - PVx;
    const float dy = v.position.y - PVy;
    const float dz = v.position.z - PVz;
    const float d = std::sqrt(dx * dx + dy * dy + dz * dz);
    const bool ok = (d < max_d2PV_mm) && (std::abs(v.position.z) < max_absz_mm);
    out.push_back(ok ? 1 : 0);
  }
  return out;
}

// =================== Tagger-ready helpers ====================
// Extract per-component float vectors from a Vec<TLorentzVector>.
// `get_Vertex_p4` already gives us Vertex 4-momenta as TLorentzVectors built
// by summing constituent RP 4-momenta — these are the cheapest way to expose
// px/py/pz/e/phi/theta as flat per-Vertex columns the tagger expects.
inline ROOT::VecOps::RVec<float>
get_p4_px(ROOT::VecOps::RVec<TLorentzVector> p4) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(p4.size());
  for (auto &v : p4)
    out.push_back(static_cast<float>(v.Px()));
  return out;
}
inline ROOT::VecOps::RVec<float>
get_p4_py(ROOT::VecOps::RVec<TLorentzVector> p4) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(p4.size());
  for (auto &v : p4)
    out.push_back(static_cast<float>(v.Py()));
  return out;
}
inline ROOT::VecOps::RVec<float>
get_p4_pz(ROOT::VecOps::RVec<TLorentzVector> p4) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(p4.size());
  for (auto &v : p4)
    out.push_back(static_cast<float>(v.Pz()));
  return out;
}
inline ROOT::VecOps::RVec<float>
get_p4_e(ROOT::VecOps::RVec<TLorentzVector> p4) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(p4.size());
  for (auto &v : p4)
    out.push_back(static_cast<float>(v.E()));
  return out;
}
inline ROOT::VecOps::RVec<float>
get_p4_phi(ROOT::VecOps::RVec<TLorentzVector> p4) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(p4.size());
  for (auto &v : p4)
    out.push_back(static_cast<float>(v.Phi()));
  return out;
}
inline ROOT::VecOps::RVec<float>
get_p4_theta(ROOT::VecOps::RVec<TLorentzVector> p4) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(p4.size());
  for (auto &v : p4)
    out.push_back(static_cast<float>(v.Theta()));
  return out;
}

// Per-RP lookup into a per-Vertex attribute via `RP_vert_ind`. Used to expose
// `RP_vert_e` / `RP_vert_mass` (the energy / invariant mass of the vertex an
// RP is assigned to). For RPs not assigned to any vertex (`RP_vert_ind < 0`)
// the sentinel value is emitted instead.
inline ROOT::VecOps::RVec<float>
get_RP_vert_attr(ROOT::VecOps::RVec<int> rp_vert_ind,
                 ROOT::VecOps::RVec<float> vertex_attr,
                 float sentinel = -1.0f) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(rp_vert_ind.size());
  for (int idx : rp_vert_ind) {
    if (idx >= 0 && static_cast<size_t>(idx) < vertex_attr.size())
      out.push_back(vertex_attr[idx]);
    else
      out.push_back(sentinel);
  }
  return out;
}

// Wrap (-π, π]. Used by `delta_phi`.
inline float wrap_pi(double d) {
  const double TWOPI = 2.0 * M_PI;
  while (d > M_PI)
    d -= TWOPI;
  while (d <= -M_PI)
    d += TWOPI;
  return static_cast<float>(d);
}

// φ delta to a scalar reference angle (e.g. EVT_thrust_phi), wrapped to (-π,
// π].
inline ROOT::VecOps::RVec<float> delta_phi(ROOT::VecOps::RVec<float> phi,
                                           float ref) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(phi.size());
  for (float p : phi)
    out.push_back(wrap_pi(static_cast<double>(p) - ref));
  return out;
}

// θ delta to a scalar reference angle (e.g. EVT_thrust_theta).
inline ROOT::VecOps::RVec<float> delta_theta(ROOT::VecOps::RVec<float> theta,
                                             float ref) {
  ROOT::VecOps::RVec<float> out;
  out.reserve(theta.size());
  for (float t : theta)
    out.push_back(t - ref);
  return out;
}

// Build a VertexObject (RVec<FCCAnalysesVertex>) directly from file-side
// PrimaryVertex + SecondaryVertices on post-beff550 EDM4hep — bypasses
// myUtils::get_VertexObject which requires MC-Reco associations that the
// new converter dropped. The fitter inside get_VertexObject re-fits
// vertices from reco tracks AND uses MC-Reco to *label* them; we don't
// need either step because the converter already wrote the fitted PV and
// SVs. Mapping:
//
//   FCCAnalysesVertex.vertex     <- edm4hep::VertexData (copy)
//   FCCAnalysesVertex.ntracks    <- (particles_end - particles_begin)
//   FCCAnalysesVertex.reco_ind   <- indices into PandoraPFOs from
//                                   _PrimaryVertex_particles /
//                                   _SecondaryVertices_particles
//   FCCAnalysesVertex.mc_ind     <- -1 (no MC truth on new schema)
//
//   The fitter-output RVecs (updated_track_momentum_at_vertex,
//   updated_track_parameters, final_track_phases, reco_chi2) MUST be
//   sized equal to reco_ind — myUtils::get_RP_atVertex (and friends)
//   loop `for i in [0, reco_ind.size())` and do `.at(i)` on these
//   fields, so an empty RVec there triggers a bounds-checked crash.
//   We fill them with sensible defaults:
//     - updated_track_momentum_at_vertex[i] = TVector3(PFO_i px,py,pz)
//       (no refit; file-side PFO momentum is the post-fit value)
//     - updated_track_parameters[i] = TVectorD(5, zeros)   (sentinel; only
//                                     line 1412/1423 of myUtils.cc read
//                                     [0]/[3] for d0/z0 — those reads
//                                     belong to RP_d0_atVertex code that
//                                     we don't need for the tagger)
//     - final_track_phases[i] = 0
//     - reco_chi2[i] = 0
//
// PV goes first (so VertexObject[0] == PV, matches the convention the
// legacy fitter uses and what get_Vertex_isPV expects).
inline ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex>
buildFCCAnalysesVertexFromFileSide(
    ROOT::VecOps::RVec<edm4hep::VertexData> primaryVertex,
    ROOT::VecOps::RVec<edm4hep::VertexData> secondaryVertices,
    ROOT::VecOps::RVec<podio::ObjectID> primaryVertexParticles,
    ROOT::VecOps::RVec<podio::ObjectID> secondaryVerticesParticles,
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco) {

  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;

  auto fill = [&](const edm4hep::VertexData &vd,
                  const ROOT::VecOps::RVec<podio::ObjectID> &rel,
                  bool is_pv) {
    VertexingUtils::FCCAnalysesVertex v;
    v.vertex = vd;
    // ZHfunctions::get_Vertex_isPV flags entries with mc_ind == 0 as the
    // primary vertex. Tag PVs accordingly so Vertex_isPV downstream is
    // [1, 0, 0, ...] (PV first, then SVs).
    v.mc_ind = is_pv ? 0 : -1;

    const int beg = vd.particles_begin;
    const int end = vd.particles_end;
    for (int k = beg; k < end && k < static_cast<int>(rel.size()); ++k) {
      const int pfo_idx = rel[k].index;
      v.reco_ind.push_back(pfo_idx);

      // updated_track_momentum_at_vertex[i] = PFO_i momentum (no refit)
      if (pfo_idx >= 0 && pfo_idx < static_cast<int>(reco.size())) {
        const auto &m = reco[pfo_idx].momentum;
        v.updated_track_momentum_at_vertex.emplace_back(m.x, m.y, m.z);
      } else {
        v.updated_track_momentum_at_vertex.emplace_back(0.f, 0.f, 0.f);
      }

      // Sentinel-sized fitter outputs so .at(i) loops in myUtils.cc
      // (d0/z0 / chi2 / phase reads) don't go out of range. Values are
      // unused by the parton-flavour / B-hadron tagger pipeline.
      v.updated_track_parameters.emplace_back(5);  // TVectorD(5) zeros
      v.final_track_phases.push_back(0.f);
      v.reco_chi2.push_back(0.f);
    }
    v.ntracks = static_cast<int>(v.reco_ind.size());
    result.push_back(v);
  };

  for (const auto &vd : primaryVertex)
    fill(vd, primaryVertexParticles, /*is_pv=*/true);

  // Outer-envelope drop only: |d2PV| < 1000 mm AND |z| < 500 mm. This
  // removes pair-seed extrapolation outliers / beam-pipe + cryostat
  // interactions but KEEPS V0 candidates (K0_S, Λ) that the model needs
  // for strangeness / baryon tagging. The per-vertex isInDet / isV0
  // flags below let the model distinguish "B/D-like SV" from
  // "V0-like SV" — same expose-don't-cut principle as RP_passLvlock.
  const float MAX_D2PV_MM = 1000.0f;
  const float MAX_ABSZ_MM = 500.0f;
  float pvx = 0.f, pvy = 0.f, pvz = 0.f;
  if (!primaryVertex.empty()) {
    pvx = primaryVertex[0].position.x;
    pvy = primaryVertex[0].position.y;
    pvz = primaryVertex[0].position.z;
  }
  for (const auto &vd : secondaryVertices) {
    const float dx = vd.position.x - pvx;
    const float dy = vd.position.y - pvy;
    const float dz = vd.position.z - pvz;
    const float d = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (d >= MAX_D2PV_MM) continue;            // genuine beam-pipe / pair-seed junk
    if (std::abs(vd.position.z) >= MAX_ABSZ_MM) continue;
    fill(vd, secondaryVerticesParticles, /*is_pv=*/false);
  }

  return result;
}


// Per-vertex inDet flag (|d2PV|<50mm AND |z|<150mm), analogous to
// RP_passLvlock. Exposed as a feature, NOT used to filter the
// VertexObject collection.
inline ROOT::VecOps::RVec<int>
get_Vertex_isInDet(
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vobj,
    float max_d2PV_mm = 50.0f, float max_absz_mm = 150.0f) {
  ROOT::VecOps::RVec<int> out;
  out.reserve(vobj.size());
  if (vobj.empty()) return out;
  const auto &pv = vobj[0].vertex.position;
  for (auto &v : vobj) {
    const float dx = v.vertex.position.x - pv.x;
    const float dy = v.vertex.position.y - pv.y;
    const float dz = v.vertex.position.z - pv.z;
    const float d = std::sqrt(dx * dx + dy * dy + dz * dz);
    const bool inDet = (d < max_d2PV_mm) &&
                       (std::abs(v.vertex.position.z) < max_absz_mm);
    out.push_back(inDet ? 1 : 0);
  }
  return out;
}


// Per-vertex V0 candidate flag. A V0 (K0_S → π+π− or Λ → pπ−) is
// characterised by:
//   - displaced position (typical K0_S γcτ ~ 13 cm, Λ γcτ ~ 30 cm at LEP1)
//   - exactly 2 charged daughters
//   - invariant mass near m_K0S (497 MeV) or m_Lambda (1115 MeV)
// We flag any 2-track SV whose mass falls in the K0_S or Λ window. Both
// flags can be true simultaneously near the overlap region; downstream
// the model can use them as soft tags.
inline ROOT::VecOps::RVec<int>
get_Vertex_isV0(
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vobj,
    ROOT::VecOps::RVec<float> v_mass) {
  ROOT::VecOps::RVec<int> out;
  out.reserve(vobj.size());
  for (size_t i = 0; i < vobj.size(); ++i) {
    // PV first, V0s never have isPV=1
    if (vobj[i].mc_ind == 0) { out.push_back(0); continue; }
    if (vobj[i].ntracks != 2) { out.push_back(0); continue; }
    if (i >= v_mass.size())   { out.push_back(0); continue; }
    const float m = v_mass[i];
    // Loose K0_S mass window (442–552 MeV, ±2σ for DELPHI resolution)
    const bool isKs = (m > 0.442f && m < 0.552f);
    // Loose Λ mass window (1.103–1.127 MeV)
    const bool isLm = (m > 1.103f && m < 1.127f);
    out.push_back((isKs || isLm) ? 1 : 0);
  }
  return out;
}


// ============================================================================
// Angle-based RP -> MCParticle matcher.
//
// The new-schema EDM4hep production (delphi_sdst_to_edm4hep) doesn't write
// MCRecoAssociations, so the canonical helper getRP2MC_index returns all -1
// and RP_fromB* / Vertex_fromB* downstream come out as zero arrays. This
// helper reconstructs the link by matching every reco PFO to its closest
// stable MC particle in (eta, phi) with a loose |p| consistency cut. It is
// roughly equivalent to what an MC-truth-table emitting converter would
// produce, accurate enough for per-RP truth-PDG queries and per-hemisphere
// truth label assignment.
//
// Convention: return value is a per-RP index into MCParticles (i.e. into
// the `mc` argument); -1 means "no good match". This matches the signature
// the rest of the pipeline expects from RP_MCidx and lets get_RP_isDescendant
// / get_RP_isDescendantAny work unchanged.
inline ROOT::VecOps::RVec<int>
matchRPtoMCByAngle(
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
    double dR_max = 0.02,
    double p_rel_tol = 0.20) {
  ROOT::VecOps::RVec<int> out(reco.size(), -1);
  std::vector<char> taken(mc.size(), 0);

  // Pre-compute MC angles & momenta for stable (status==1) particles only.
  std::vector<int>    cand_idx;
  std::vector<double> cand_theta, cand_phi, cand_p;
  cand_idx.reserve(mc.size());
  cand_theta.reserve(mc.size());
  cand_phi.reserve(mc.size());
  cand_p.reserve(mc.size());
  for (size_t j = 0; j < mc.size(); ++j) {
    const auto &m = mc[j];
    if (m.generatorStatus != 1) continue;
    double pp = std::sqrt(m.momentum.x * m.momentum.x +
                          m.momentum.y * m.momentum.y +
                          m.momentum.z * m.momentum.z);
    if (pp <= 0.) continue;
    cand_idx.push_back(static_cast<int>(j));
    cand_theta.push_back(std::acos(m.momentum.z / pp));
    cand_phi.push_back(std::atan2(m.momentum.y, m.momentum.x));
    cand_p.push_back(pp);
  }

  for (size_t i = 0; i < reco.size(); ++i) {
    const auto &r = reco[i];
    double rp = std::sqrt(r.momentum.x * r.momentum.x +
                          r.momentum.y * r.momentum.y +
                          r.momentum.z * r.momentum.z);
    if (rp <= 0.) continue;
    double r_theta = std::acos(r.momentum.z / rp);
    double r_phi   = std::atan2(r.momentum.y, r.momentum.x);

    int best_k = -1;
    double best_score = dR_max * dR_max;   // squared dR threshold
    for (size_t k = 0; k < cand_idx.size(); ++k) {
      if (taken[cand_idx[k]]) continue;
      // Loose momentum consistency to avoid pathological 30 GeV photon →
      // 0.4 GeV pion matches that would otherwise pass dR.
      if (std::abs(rp - cand_p[k]) / cand_p[k] > p_rel_tol) continue;
      double dphi = cand_phi[k] - r_phi;
      while (dphi >  M_PI) dphi -= 2 * M_PI;
      while (dphi < -M_PI) dphi += 2 * M_PI;
      double dtheta = cand_theta[k] - r_theta;
      double dr2 = dtheta * dtheta + dphi * dphi;
      if (dr2 < best_score) {
        best_score = dr2;
        best_k = static_cast<int>(k);
      }
    }
    if (best_k >= 0) {
      out[i] = cand_idx[best_k];
      taken[cand_idx[best_k]] = 1;
    }
  }
  return out;
}


// Per-RP truth PDG, given the angle-matched index from matchRPtoMCByAngle.
// Returns 0 for unmatched RPs (sentinel different from any real PDG).
inline ROOT::VecOps::RVec<int>
getRPMatchedPDG(
    ROOT::VecOps::RVec<int> matchIdx,
    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {
  ROOT::VecOps::RVec<int> out(matchIdx.size(), 0);
  for (size_t i = 0; i < matchIdx.size(); ++i) {
    int j = matchIdx[i];
    if (j >= 0 && j < static_cast<int>(mc.size())) out[i] = mc[j].PDG;
  }
  return out;
}


// Per-RP value pulled from any ParticleID_* collection in the new schema.
// EDM4hep PID layout: each ParticleIDData has [parameters_begin,
// parameters_end] indices into a flat float vector
// `_ParticleID_<name>_parameters`. The link from PID -> RP is via the
// OneToOne `particle` relation stored as a parallel ObjectID collection
// `_ParticleID_<name>_particle.index` (per-event integer RVec).
//
// This is a single generic helper: pick `param_index` (offset within each
// PID's parameter slice) and the function returns the per-RP RVec of that
// parameter value, with 0.0 for RPs that have no PID attached. We use:
//   * dEdx (2 params/PID): param_index=0 -> dEdx value, 1 -> sigma
//   * HadronRich (18 params/PID): param_index 0..17 for individual tags
inline ROOT::VecOps::RVec<float>
getRPPIDParam(
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
    ROOT::VecOps::RVec<edm4hep::ParticleIDData> pids,
    ROOT::VecOps::RVec<float> pid_params_flat,
    ROOT::VecOps::RVec<int> pid_particle_idx,
    int param_index = 0) {
  ROOT::VecOps::RVec<float> out(reco.size(), 0.f);
  const int n_pids = static_cast<int>(pids.size());
  const int n_reco = static_cast<int>(reco.size());
  const int n_link = static_cast<int>(pid_particle_idx.size());
  const int n_flat = static_cast<int>(pid_params_flat.size());
  for (int i = 0; i < n_pids; ++i) {
    int rp_idx = (i < n_link) ? pid_particle_idx[i] : -1;
    if (rp_idx < 0 || rp_idx >= n_reco) continue;
    int b = pids[i].parameters_begin + param_index;
    int e = pids[i].parameters_end;
    if (b >= 0 && b < n_flat && b < e) {
      out[rp_idx] = pid_params_flat[b];
    }
  }
  return out;
}


// Per-hemisphere truth label assignment using gen-level B-hadrons.
//
// Given per-event RVecs of (px, py, pz) for each gen-B species
// (genBd, genBu, genBs, genBc, genLb) plus the event's thrust unit
// vector (computed from EVT_thrust_theta / EVT_thrust_phi), return a
// per-event integer "hemisphereLabelMask" of length 2 with:
//
//   index 0 = Emin side (thrust_angle > 0): which B species sits here
//   index 1 = Emax side (thrust_angle < 0): same
//
// Encoding: 0=none, 1=Bd, 2=Bu, 3=Bs, 4=has1Bc, 5=Lb
//
// Use case: prep_bhadron.py reads this and overrides the channel-based
// label with the per-hemisphere true content.
inline std::vector<int>
classifyHemispheresByGenB(
    float thrust_theta, float thrust_phi,
    ROOT::VecOps::RVec<float> Bd_px, ROOT::VecOps::RVec<float> Bd_py, ROOT::VecOps::RVec<float> Bd_pz,
    ROOT::VecOps::RVec<float> Bu_px, ROOT::VecOps::RVec<float> Bu_py, ROOT::VecOps::RVec<float> Bu_pz,
    ROOT::VecOps::RVec<float> Bs_px, ROOT::VecOps::RVec<float> Bs_py, ROOT::VecOps::RVec<float> Bs_pz,
    ROOT::VecOps::RVec<float> Bc_px, ROOT::VecOps::RVec<float> Bc_py, ROOT::VecOps::RVec<float> Bc_pz,
    ROOT::VecOps::RVec<float> Lb_px, ROOT::VecOps::RVec<float> Lb_py, ROOT::VecOps::RVec<float> Lb_pz) {
  // Thrust axis unit vector.
  const double tx = std::sin(thrust_theta) * std::cos(thrust_phi);
  const double ty = std::sin(thrust_theta) * std::sin(thrust_phi);
  const double tz = std::cos(thrust_theta);

  // For each side, find dominant B species (count of gen-particles in that
  // hemisphere); ties broken by Bd<Bu<Bs<Bc<Lb declaration order.
  int n_emin[6] = {0, 0, 0, 0, 0, 0};
  int n_emax[6] = {0, 0, 0, 0, 0, 0};
  auto fold = [&](const ROOT::VecOps::RVec<float> &px,
                  const ROOT::VecOps::RVec<float> &py,
                  const ROOT::VecOps::RVec<float> &pz,
                  int species) {
    for (size_t i = 0; i < px.size(); ++i) {
      double dot = tx * px[i] + ty * py[i] + tz * pz[i];
      if (dot > 0) n_emin[species]++; else n_emax[species]++;
    }
  };
  fold(Bd_px, Bd_py, Bd_pz, 1);
  fold(Bu_px, Bu_py, Bu_pz, 2);
  fold(Bs_px, Bs_py, Bs_pz, 3);
  fold(Bc_px, Bc_py, Bc_pz, 4);
  fold(Lb_px, Lb_py, Lb_pz, 5);

  auto pick = [](const int *cnt) {
    int best = 0, best_n = 0;
    for (int k = 1; k <= 5; ++k) {
      if (cnt[k] > best_n) {
        best_n = cnt[k];
        best = k;
      }
    }
    return best;
  };
  return std::vector<int>{pick(n_emin), pick(n_emax)};
}

// |PDG| of the first quark (|PDG| in 1..6) appearing in MCParticles.
// For inclusive Pythia 8 e+e- -> Z -> qq samples this is the Z-decay
// quark flavour (idx 5 after e-,e+,e-,e+,Z). Returns 0 if no quark is
// found (defensive, shouldn't happen for hadronic events).
inline int get_qqPDG(ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  for (size_t i = 0; i < in.size(); ++i) {
    int p = std::abs(in[i].PDG);
    if (p >= 1 && p <= 6) return p;
  }
  return 0;
}


// -----------------------------------------------------------------------
// arXiv:2510.18762v1 Table-1 selection (charged+neutral thrust).
//
//   Per-particle (charged):  20 <= theta <= 160 deg, pT > 0.4 GeV
//   Per-particle (neutral):  20 <= theta <= 160 deg, E  > 0.5 GeV
//   Per-event:               n_ch_sel >= 7,
//                            E_tot_sel >= 0.5 * E_cm = 45.6 GeV,
//                            30 <= theta_thrust <= 150 deg
//
// The helpers below compute the per-particle-cut quantities (n_ch_sel,
// n_neu_sel, E_tot_sel) so the analyzer can both *persist* them as
// per-event scalars in stage1 output and Filter on the combined pass
// flag (paperPassEvent).
// -----------------------------------------------------------------------
inline int paperNChSel(ROOT::VecOps::RVec<float> rp_theta,
                       ROOT::VecOps::RVec<float> rp_charge,
                       ROOT::VecOps::RVec<float> rp_px,
                       ROOT::VecOps::RVec<float> rp_py) {
  const float theta_lo = 20.0f  * static_cast<float>(M_PI) / 180.0f;
  const float theta_hi = 160.0f * static_cast<float>(M_PI) / 180.0f;
  const float pt_min   = 0.4f;
  int n = 0;
  for (size_t i = 0; i < rp_theta.size(); ++i) {
    if (rp_theta[i] < theta_lo || rp_theta[i] > theta_hi) continue;
    if (std::abs(rp_charge[i]) <= 0.1f) continue;
    const float pt = std::sqrt(rp_px[i] * rp_px[i] + rp_py[i] * rp_py[i]);
    if (pt > pt_min) ++n;
  }
  return n;
}

inline int paperNNeuSel(ROOT::VecOps::RVec<float> rp_theta,
                        ROOT::VecOps::RVec<float> rp_charge,
                        ROOT::VecOps::RVec<float> rp_e) {
  const float theta_lo = 20.0f  * static_cast<float>(M_PI) / 180.0f;
  const float theta_hi = 160.0f * static_cast<float>(M_PI) / 180.0f;
  const float e_min    = 0.5f;
  int n = 0;
  for (size_t i = 0; i < rp_theta.size(); ++i) {
    if (rp_theta[i] < theta_lo || rp_theta[i] > theta_hi) continue;
    if (std::abs(rp_charge[i]) > 0.1f) continue;
    if (rp_e[i] > e_min) ++n;
  }
  return n;
}

inline float paperETotSel(ROOT::VecOps::RVec<float> rp_theta,
                          ROOT::VecOps::RVec<float> rp_charge,
                          ROOT::VecOps::RVec<float> rp_px,
                          ROOT::VecOps::RVec<float> rp_py,
                          ROOT::VecOps::RVec<float> rp_e) {
  const float theta_lo = 20.0f  * static_cast<float>(M_PI) / 180.0f;
  const float theta_hi = 160.0f * static_cast<float>(M_PI) / 180.0f;
  const float pt_min   = 0.4f;
  const float e_min    = 0.5f;
  float E = 0.0f;
  for (size_t i = 0; i < rp_theta.size(); ++i) {
    if (rp_theta[i] < theta_lo || rp_theta[i] > theta_hi) continue;
    if (std::abs(rp_charge[i]) > 0.1f) {
      const float pt = std::sqrt(rp_px[i] * rp_px[i] + rp_py[i] * rp_py[i]);
      if (pt > pt_min) E += rp_e[i];
    } else {
      if (rp_e[i] > e_min) E += rp_e[i];
    }
  }
  return E;
}

// Final event pass flag, given the precomputed counters + scalar
// EVT_thrust_theta (radians).
inline int paperPassEvent(int n_ch_sel, float e_tot_sel, float theta_thrust) {
  constexpr float E_CM           = 91.2f;
  constexpr int   N_CH_MIN       = 7;
  constexpr float E_TOT_FRAC     = 0.5f;
  constexpr float THETA_TH_LO    = 30.0f  * static_cast<float>(M_PI) / 180.0f;
  constexpr float THETA_TH_HI    = 150.0f * static_cast<float>(M_PI) / 180.0f;
  if (n_ch_sel < N_CH_MIN) return 0;
  if (e_tot_sel < E_TOT_FRAC * E_CM) return 0;
  if (theta_thrust < THETA_TH_LO || theta_thrust > THETA_TH_HI) return 0;
  return 1;
}


} // namespace ZHfunctions
} // namespace FCCAnalyses

#endif
