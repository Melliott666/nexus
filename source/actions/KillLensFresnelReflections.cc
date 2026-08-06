// ----------------------------------------------------------------------------
// nexus | KillLensFresnelReflections.cc
//
// Kill optical photons that undergo Fresnel reflection at either fused-silica
// lens while leaving normally refracted photons unchanged. This isolates lens
// reflection ghosts without modifying the lens geometry or material.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "FactoryBase.h"

#include <G4LogicalVolume.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4OpticalPhoton.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4UserSteppingAction.hh>
#include <G4VPhysicalVolume.hh>

#include <cstddef>


namespace nexus {

class KillLensFresnelReflections: public G4UserSteppingAction
{
public:
  KillLensFresnelReflections(): boundary_(0) {}

  void UserSteppingAction(const G4Step* step) override
  {
    G4Track* track = step->GetTrack();

    if (track->GetDefinition() != G4OpticalPhoton::Definition())
      return;

    G4StepPoint* post = step->GetPostStepPoint();
    if (post->GetStepStatus() != fGeomBoundary)
      return;

    FindBoundaryProcess(track);
    if (!boundary_ || boundary_->GetStatus() != FresnelReflection)
      return;

    const G4VPhysicalVolume* pre_volume =
      step->GetPreStepPoint()->GetPhysicalVolume();
    const G4VPhysicalVolume* post_volume = post->GetPhysicalVolume();
    const G4VPhysicalVolume* next_volume = track->GetNextVolume();

    if (IsLens(pre_volume) || IsLens(post_volume) || IsLens(next_volume))
      track->SetTrackStatus(fStopAndKill);
  }

private:
  static G4bool IsLens(const G4VPhysicalVolume* volume)
  {
    return volume && volume->GetLogicalVolume()->GetName() == "FS_LENS";
  }

  void FindBoundaryProcess(const G4Track* track)
  {
    if (boundary_)
      return;

    G4ProcessVector* processes =
      track->GetDefinition()->GetProcessManager()->GetProcessList();

    for (std::size_t i = 0; i < processes->size(); ++i) {
      boundary_ = dynamic_cast<G4OpBoundaryProcess*>((*processes)[i]);
      if (boundary_)
        return;
    }
  }

  G4OpBoundaryProcess* boundary_;
};

REGISTER_CLASS(KillLensFresnelReflections, G4UserSteppingAction)

} // namespace nexus
