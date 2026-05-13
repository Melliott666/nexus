#ifndef NEXUS_KINGCRABPRACTICE_H
#define NEXUS_KINGCRABPRACTICE_H
#include "GeometryBase.h"
class G4GenericMessenger;

namespace nexus { class BoxPointSampler; class CylinderPointSampler;}

namespace nexus {
    
    class KingCRABPractice : public GeometryBase
{
        public:
            // Constructor / Destructor
            KingCRABPractice();
            ~KingCRABPractice();

    
    virtual void Construct();
    virtual G4ThreeVector GenerateVertex(const G4String& region) const;

        private:
            /// Messenger for the definition of control commands
            G4GenericMessenger* msg_;
            G4double Lab_size;
            G4double gas_pressure_;
            G4double max_step_size_;
            BoxPointSampler* practice_track_;
            G4String gastype_;
            G4ThreeVector specific_vertex_;
            
            void ConstructLab();
            void PlaceVolumes();
            void AssignVisuals();
            void PrintParam();

    };

} // end namespace nexus


#endif //NEXUS_KingCRAB_H