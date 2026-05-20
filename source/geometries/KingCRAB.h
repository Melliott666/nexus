//
// Geometry of the King CRAB detector
//

#ifndef NEXUS_KINGCRAB_H
#define NEXUS_KINGCRAB_H

#include "GeometryBase.h"

class G4GenericMessenger;

namespace nexus {

class BoxPointSampler;
class CylinderPointSampler;

class KingCRAB: public GeometryBase
{
  public:
    /// Constructor
    KingCRAB();

    /// Destructor
    ~KingCRAB();

    virtual void Construct();
    virtual G4ThreeVector GenerateVertex(const G4String& region) const;

  private:
    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;

    // Global geometry controls
    G4double Lab_size;
    G4double gas_pressure_;
    G4double max_step_size_;
    G4String gastype_;
    G4ThreeVector specific_vertex_;

    // Active/drift field controls
    G4bool drift_field_on_;
    G4double drift_field_int_;
    G4double drift_v_;
    G4double drift_e_lifetime_;

    // EL field controls
    G4bool el_field_on_;
    G4double el_field_int_;
    G4double EL_drift_v_;

    // Vertex generators
    CylinderPointSampler* anode_gen_;
    CylinderPointSampler* cathode_gen_;
    CylinderPointSampler* active_volume_gen_;
    BoxPointSampler* practice_track_;

    void ConstructLab();
    void PlaceVolumes();
    void AssignVisuals();
    void PrintParam();
};

} // namespace nexus

#endif // NEXUS_KINGCRAB_H
