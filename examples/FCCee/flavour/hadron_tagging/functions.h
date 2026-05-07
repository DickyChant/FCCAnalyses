#ifndef ZHfunctions_H
#define ZHfunctions_H

// Trimmed for the modern smoke-test stage1. Original had many helpers
// keyed off EDM4hep 0.x field names (vertex.primary, TrackData.dxQuantities_begin,
// TrackerHitData, RP::type). Restore them piecewise as needed when porting
// the full feature set (and updating to EDM4hep 1.0 schema).

#include <cmath>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"

namespace FCCAnalyses { namespace ZHfunctions {

using rp = edm4hep::ReconstructedParticleData;
using Vec_rp = ROOT::VecOps::RVec<rp>;

// missing energy (vs. an assumed CM-frame total energy)
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

}}  // namespace

#endif
