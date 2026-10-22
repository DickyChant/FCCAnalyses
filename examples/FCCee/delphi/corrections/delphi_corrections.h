// Corrections for DELPHI open data from delphi-correctionlib
// (https://github.com/delphi-fullDST-edm4hep/delphi-correctionlib), applied in
// RDataFrame to the output of the delphi-edm4hep converter.
//
// correctionlib has to be loaded into ROOT before this header is included;
// register() in delphi_corrections.py does both. See README.md.
#ifndef DELPHI_CORRECTIONS_H
#define DELPHI_CORRECTIONS_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "ROOT/RVec.hxx"
#include "correction.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/VertexData.h"
#include "podio/ObjectID.h"

namespace FCCAnalyses {
namespace DelphiCorrections {

// An int frame parameter (columns GPIntKeys, GPIntValues), or `fallback`.
inline int intParameter(const ROOT::RVec<std::string> &keys,
                        const ROOT::RVec<std::vector<int>> &values,
                        const std::string &key, int fallback = 0) {
  for (std::size_t i = 0; i < keys.size(); ++i) {
    if (keys[i] == key && !values[i].empty()) return values[i][0];
  }
  return fallback;
}

// A string frame parameter (columns GPStringKeys, GPStringValues), or "".
inline std::string
stringParameter(const ROOT::RVec<std::string> &keys,
                const ROOT::RVec<std::vector<std::string>> &values,
                const std::string &key) {
  for (std::size_t i = 0; i < keys.size(); ++i) {
    if (keys[i] == key && !values[i].empty()) return values[i][0];
  }
  return "";
}

// Production radius r0 in xy (mm) of each particle, as the momentum-scale
// correction defines it: the radius of the V0 vertex for daughters of a DELPHI
// V0 (K0S, Lambda or photon conversion) with tag >= 1 (parameters[0], KRV0(5)),
// as in the calibration; 0 for every other particle. A particle in more than
// one such V0 takes the first.
//
// v0Particles must index into `particles`: pass the V0 collection of the same
// pass (sDST_V0_V0Candidates with sDST_MAIN_Particles, fDST_V0_V0Candidates
// with fDST_MAIN_Particles).
inline ROOT::RVec<float>
v0DaughterRadius(const ROOT::RVec<edm4hep::ReconstructedParticleData> &particles,
                 const ROOT::RVec<edm4hep::VertexData> &v0s,
                 const ROOT::RVec<podio::ObjectID> &v0Particles,
                 const ROOT::RVec<float> &v0Parameters) {
  ROOT::RVec<float> r0(particles.size(), 0.f);
  for (const auto &v0 : v0s) {
    if (v0.parameters_begin >= v0.parameters_end ||
        v0Parameters[v0.parameters_begin] < 1.f)
      continue;
    const float radius = std::hypot(v0.position.x, v0.position.y);
    for (unsigned k = v0.particles_begin; k < v0.particles_end; ++k) {
      const int index = v0Particles[k].index;
      if (index >= 0 && index < static_cast<int>(particles.size()) &&
          r0[index] == 0.f)
        r0[index] = radius;
    }
  }
  return r0;
}

// The charged-track momentum scale, correction track_momentum_scale of
// delphi-correctionlib's tracking/momentumScale_9495.json.gz.
//
// Returns `particles` with the momentum of every charged particle multiplied
// by the factor s(dataset, variation, r0, theta, p, q) and the energy
// recomputed from the particle's own mass; neutral particles are unchanged.
// The order and the relations are kept, so indices into the collection stay
// valid.
//
// Data only: events with a negative run number (simulation) are returned
// unchanged. The dataset key is the processing tag of the file
// (sDST_EVT_dstProcessingTag, e.g. "94C2") unless one is given here; a
// processing the correction does not cover is an error, as the correction is
// only valid for the processings it was calibrated on.
//
// Only tracks inside the calibration's track selection are corrected: by
// default 0.4 < p < 50 GeV and 20 < theta < 160 degrees, the cuts of
// momentumScale_skim.py. Outside it the model is an extrapolation: its
// energy-loss term grows as 1/pT, and on 94C2 data it gives s < 0.5 for tracks
// of a few tens of MeV, and s < 0 for some. Tracks outside the window keep
// their momentum. To correct every charged track, pass pMin = 0, pMax = 0 (no
// upper limit), thetaMinDeg = 0 and thetaMaxDeg = 180.
class TrackMomentumScale {
public:
  TrackMomentumScale(const std::string &file,
                     const std::string &variation = "nominal",
                     const std::string &dataset = "", double pMin = 0.4,
                     double pMax = 50., double thetaMinDeg = 20.,
                     double thetaMaxDeg = 160.)
      : correction_(correction::CorrectionSet::from_file(file)->at(
            "track_momentum_scale")),
        variation_(variation), dataset_(dataset), pMin_(pMin),
        pMax_(pMax > 0. ? pMax : HUGE_VAL), thetaMin_(thetaMinDeg * M_PI / 180.),
        thetaMax_(thetaMaxDeg * M_PI / 180.) {}

  // The factor s per particle: 1 for neutral particles, for tracks outside the
  // window, and for simulation.
  ROOT::RVec<double>
  factors(const ROOT::RVec<edm4hep::ReconstructedParticleData> &particles,
          const ROOT::RVec<float> &r0, int runNumber,
          const std::string &processing) const {
    ROOT::RVec<double> s(particles.size(), 1.);
    if (runNumber < 0) return s;
    const std::string dataset = dataset_.empty() ? processing : dataset_;
    for (std::size_t i = 0; i < particles.size(); ++i) {
      const auto &rp = particles[i];
      if (rp.charge == 0.f) continue;
      const double px = rp.momentum.x, py = rp.momentum.y, pz = rp.momentum.z;
      const double p = std::sqrt(px * px + py * py + pz * pz);
      if (p <= 0.) continue;
      const double theta = std::atan2(std::hypot(px, py), pz);
      if (p <= pMin_ || p >= pMax_ || theta <= thetaMin_ || theta >= thetaMax_)
        continue;
      try {
        s[i] = correction_->evaluate({dataset, variation_,
                                      static_cast<double>(r0[i]), theta, p,
                                      static_cast<double>(rp.charge)});
      } catch (const std::exception &e) {
        throw std::runtime_error(
            "track_momentum_scale for processing '" + dataset +
            "', variation '" + variation_ + "': " + e.what() +
            ". The correction covers only the processings it was calibrated "
            "on; give the dataset explicitly if the file's processing tag is "
            "wrong.");
      }
    }
    return s;
  }

  ROOT::RVec<edm4hep::ReconstructedParticleData>
  operator()(const ROOT::RVec<edm4hep::ReconstructedParticleData> &particles,
             const ROOT::RVec<float> &r0, int runNumber,
             const std::string &processing) const {
    const auto s = factors(particles, r0, runNumber, processing);
    ROOT::RVec<edm4hep::ReconstructedParticleData> out(particles);
    for (std::size_t i = 0; i < out.size(); ++i) {
      if (s[i] == 1.) continue;
      auto &rp = out[i];
      rp.momentum.x *= s[i];
      rp.momentum.y *= s[i];
      rp.momentum.z *= s[i];
      const double p2 = double(rp.momentum.x) * rp.momentum.x +
                        double(rp.momentum.y) * rp.momentum.y +
                        double(rp.momentum.z) * rp.momentum.z;
      rp.energy = std::sqrt(p2 + double(rp.mass) * rp.mass);
    }
    return out;
  }

private:
  correction::Correction::Ref correction_;
  std::string variation_;
  std::string dataset_;
  double pMin_, pMax_, thetaMin_, thetaMax_; // GeV, rad
};

} // namespace DelphiCorrections
} // namespace FCCAnalyses

#endif
