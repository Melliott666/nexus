// ----------------------------------------------------------------------------
// nexus | CosmogenicMuonGenerator.h
//
// Bare-sky cosmic muons conditioned to cross a finite horizontal cylinder.
// ----------------------------------------------------------------------------

#ifndef COSMOGENIC_MUON_GENERATOR_H
#define COSMOGENIC_MUON_GENERATOR_H

#include <G4ThreeVector.hh>
#include <G4VPrimaryGenerator.hh>
#include <Randomize.hh>

class G4Event;
class G4GenericMessenger;

namespace nexus {

  class CosmogenicMuonGenerator: public G4VPrimaryGenerator
  {
  public:
    CosmogenicMuonGenerator();
    ~CosmogenicMuonGenerator();

    void GeneratePrimaryVertex(G4Event*) override;

  private:
    void InitializeSpectrum();
    void SampleEnergyAndDirection(G4double&, G4ThreeVector&);
    G4bool SampleVertex(const G4ThreeVector&, G4ThreeVector&) const;
    G4bool IntersectsTarget(const G4ThreeVector&, const G4ThreeVector&) const;
    G4double GuanIntensity(G4double, G4double) const;

    G4GenericMessenger* msg_;
    G4RandGeneral* spectrum_sampler_;
    G4bool initialized_;

    G4double energy_min_;
    G4double energy_max_;
    G4double positive_fraction_;
    G4double target_radius_;
    G4double target_length_;
    G4ThreeVector target_center_;
    G4double generation_distance_;
    G4bool paired_event_seeds_;
    G4int base_seed_;
    G4int theta_bins_;
    G4int energy_bins_;
  };

} // namespace nexus

#endif
