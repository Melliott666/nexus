//
// Geometry of the King CRAB detector
//

#ifndef NEXUS_KINGCRAB_H
#define NEXUS_KINGCRAB_H
#include "GeometryBase.h"
class G4GenericMessenger;

namespace nexus { class BoxPointSampler; class CylinderPointSampler;}

namespace nexus {
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
            G4double Lab_size;
            G4double gas_pressure_;
            G4double max_step_size_;
            CylinderPointSampler* anode_gen_;
            CylinderPointSampler* cathode_gen_;
            CylinderPointSampler* active_volume_gen_;
            CylinderPointSampler* practice_track_;
            G4String gastype_;
            G4ThreeVector specific_vertex_;
            
            void ConstructLab();
            void PlaceVolumes();
            void AssignVisuals();
            void PrintParam();

    };

} // end namespace nexus


#endif //NEXUS_KingCRAB_H
