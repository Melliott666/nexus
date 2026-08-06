// ----------------------------------------------------------------------------
// nexus | SaveAllSteppingAction.cc
//
// This class adds a new group and table to the output file, "/DEBUG/steps".
// This table contains information (position and volume of both the
// pre- and post-step points, average time, process name and other identifiers)
// of some steps of the simulation. By default all steps are stored. However,
// a subset of them can be selected by cherry-picking the volumes and particles
// involved in the step. This can be achieved with the commands
// /Actions/SaveAllSteppingAction/select_particle
// and
// /Actions/SaveAllSteppingAction/select_volume
// without the need for re-compilation.
// It must be noted that the files produced with this action become large
// very quickly. Therefore, strict filtering and small number of events are
// encouraged.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "SaveAllSteppingAction.h"
#include "PersistencyManager.h"
#include "FactoryBase.h"
#include "IonizationElectron.h"

#include <G4Step.hh>
#include <G4LogicalVolume.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4OpticalPhoton.hh>
#include <G4VPersistencyManager.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>
#include <G4ParticleTable.hh>
#include <G4Track.hh>
#include <G4VPhysicalVolume.hh>

using namespace nexus;

REGISTER_CLASS(SaveAllSteppingAction, G4UserSteppingAction)

SaveAllSteppingAction::SaveAllSteppingAction():
G4UserSteppingAction(),
msg_(0),
selected_volumes_(),
selected_particles_(),
initial_volumes_(),
final_volumes_(),
proc_names_(),
initial_poss_(),
final_poss_(),
times_(),
kill_after_selection_(false),
kill_lens_fresnel_reflections_(false),
record_selected_track_deaths_(false),
record_el_electron_entries_(false),
optical_boundary_(0),
selected_track_ids_(),
recorded_el_electron_ids_()
{
  msg_ = new G4GenericMessenger(this, "/Actions/SaveAllSteppingAction/");

  msg_->DeclareMethod("select_particle",
                      &SaveAllSteppingAction::AddSelectedParticle,
                      "Add a new particle to select");

  msg_->DeclareMethod("select_volume", &SaveAllSteppingAction::AddSelectedVolume,
                      "Add a new volume to select");

  msg_->DeclareProperty("kill_after_selection", kill_after_selection_,
                        "Whether to kill a particle after a step has been selected");

  msg_->DeclareProperty(
    "kill_lens_fresnel_reflections", kill_lens_fresnel_reflections_,
    "Kill optical photons Fresnel-reflected at either FS_LENS boundary");

  msg_->DeclareProperty(
    "record_selected_track_deaths", record_selected_track_deaths_,
    "Record the terminal step of tracks that previously touched a selected volume");

  msg_->DeclareProperty(
    "record_el_electron_entries", record_el_electron_entries_,
    "Record each ionization electron on its EL-photon-producing drift step");

  PersistencyManager* pm = dynamic_cast<PersistencyManager*>
        (G4VPersistencyManager::GetPersistencyManager());

  pm->StoreSteps(true);

}



SaveAllSteppingAction::~SaveAllSteppingAction()
{
}



void SaveAllSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track*              track         = step->GetTrack();
  G4ParticleDefinition* pdef          = track->GetDefinition();
  G4int                 track_id      = track->GetTrackID();
  G4String              particle_name = pdef->GetParticleName();

  G4bool is_ionization_electron =
    pdef == IonizationElectron::Definition();
  G4bool consider_el_entry =
    record_el_electron_entries_ && is_ionization_electron;

  if (!consider_el_entry && !KeepParticle(pdef)) return;

  // This optional filter lets KingCRAB test the lens-reflection ghost
  // hypothesis while retaining the tightly filtered /DEBUG/steps output.
  // Nexus accepts only one user stepping action, so this cannot be run as a
  // separate KillLensFresnelReflections action alongside this recorder.
  G4bool killed_lens_fresnel =
    kill_lens_fresnel_reflections_ && IsLensFresnelReflection(step);
  if (killed_lens_fresnel)
    track->SetTrackStatus(fStopAndKill);

  G4StepPoint* pre  = step->GetPreStepPoint();
  G4StepPoint* post = step->GetPostStepPoint();

  G4ThreeVector initial_pos = pre ->GetPosition();
  G4ThreeVector   final_pos = post->GetPosition();
  G4double        step_time = (pre->GetGlobalTime()  +
                              post->GetGlobalTime()) / 2.;

  const G4VPhysicalVolume* initial_physical = pre ->GetPhysicalVolume();
  const G4VPhysicalVolume*   final_physical = post->GetPhysicalVolume();
  G4String initial_volume = initial_physical ? initial_physical->GetName()
                                             : "OUT_OF_WORLD";
  G4String   final_volume = final_physical ? final_physical->GetName()
                                           : "OUT_OF_WORLD";

  // KingCRAB's fast drift transports an ionization electron directly to the
  // anode in one step and generates EL photons at sampled points along that
  // drift line. It therefore need not expose EL_GAP as either touchable. Use
  // the presence of optical secondaries to identify the EL-producing step;
  // its post-step x-y is the diffused charge-arrival coordinate.
  G4bool produced_el_photons = false;
  const auto* secondaries = step->GetSecondaryInCurrentStep();
  if (consider_el_entry && secondaries) {
    for (const G4Track* secondary : *secondaries) {
      if (secondary->GetDefinition() == G4OpticalPhoton::Definition()) {
        produced_el_photons = true;
        break;
      }
    }
  }
  G4bool first_el_entry = consider_el_entry && produced_el_photons &&
                          !recorded_el_electron_ids_.count(track_id);
  if (first_el_entry)
    recorded_el_electron_ids_.insert(track_id);

  G4TrackStatus status   = track->GetTrackStatus();
  G4bool        terminal = status == fStopAndKill ||
                           status == fKillTrackAndSecondaries;
  if (consider_el_entry && !first_el_entry) return;

  G4bool selected_step = !consider_el_entry &&
                         KeepVolume(initial_volume, final_volume);
  if (selected_step)
    selected_track_ids_.insert(track_id);

  G4bool selected_death = record_selected_track_deaths_ && terminal &&
                          selected_track_ids_.count(track_id);

  if (!selected_step && !selected_death && !first_el_entry) return;

  const G4VProcess* process = post->GetProcessDefinedStep();
  G4String proc_name = process ? process->GetProcessName() : "NoProcess";
  if (killed_lens_fresnel)
    proc_name = "KilledLensFresnelReflection";
  if (selected_death)
    proc_name = "DEATH:" + proc_name;
  if (first_el_entry)
    proc_name = "EL_ELECTRON_ENTRY:" + proc_name;

  std::pair<G4int, G4String> key = std::make_pair(track_id, particle_name);

  initial_volumes_[key].push_back(initial_volume);
    final_volumes_[key].push_back(  final_volume);
       proc_names_[key].push_back(     proc_name);

  initial_poss_   [key].push_back(initial_pos);
    final_poss_   [key].push_back(  final_pos);
         times_   [key].push_back(  step_time);

  if (kill_after_selection_)
    track->SetTrackStatus(fStopAndKill);
}


G4bool SaveAllSteppingAction::IsLensFresnelReflection(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  if (track->GetDefinition() != G4OpticalPhoton::Definition()) return false;

  G4StepPoint* post = step->GetPostStepPoint();
  if (post->GetStepStatus() != fGeomBoundary) return false;

  if (!optical_boundary_) {
    G4ProcessVector* processes =
      track->GetDefinition()->GetProcessManager()->GetProcessList();
    for (std::size_t i = 0; i < processes->size(); ++i) {
      optical_boundary_ =
        dynamic_cast<G4OpBoundaryProcess*>((*processes)[i]);
      if (optical_boundary_) break;
    }
  }

  if (!optical_boundary_ ||
      optical_boundary_->GetStatus() != FresnelReflection) return false;

  const auto is_lens = [](const G4VPhysicalVolume* volume) {
    return volume && volume->GetLogicalVolume()->GetName() == "FS_LENS";
  };

  const G4VPhysicalVolume*  pre_volume =
    step->GetPreStepPoint()->GetPhysicalVolume();
  const G4VPhysicalVolume* post_volume = post->GetPhysicalVolume();
  const G4VPhysicalVolume* next_volume = track->GetNextVolume();

  return is_lens(pre_volume) || is_lens(post_volume) || is_lens(next_volume);
}


void SaveAllSteppingAction::AddSelectedParticle(G4String particle_name)
{
  G4ParticleDefinition* pdef = G4ParticleTable::GetParticleTable()->FindParticle(particle_name);
  if (!pdef) {
    G4String msg = "No particle description was found for particle name " + particle_name;
    G4Exception("[SaveAllSteppingAction]", "AddSelectedParticle()", FatalException, msg);
  }
  selected_particles_.push_back(pdef);
}


void SaveAllSteppingAction::AddSelectedVolume(G4String volume_name)
{
  selected_volumes_.push_back(volume_name);
}


G4bool SaveAllSteppingAction::KeepParticle(G4ParticleDefinition* pdef)
{
  if (!selected_particles_.size()) return true;

  auto it = std::find(selected_particles_.begin(), selected_particles_.end(), pdef);
  return it != selected_particles_.end();
}


G4bool SaveAllSteppingAction::KeepVolume(G4String& initial_volume, G4String& final_volume)
{
  if (!selected_volumes_.size()) return true;

  for (auto volume=selected_volumes_.begin(); volume != selected_volumes_.end(); volume++)
  {
    if (G4StrUtil::contains(initial_volume, *volume)) return true;
    if (G4StrUtil::contains(  final_volume, *volume)) return true;
  }

  return false;
}



void SaveAllSteppingAction::Reset()
{
  initial_volumes_.clear();
    final_volumes_.clear();
       proc_names_.clear();

  initial_poss_   .clear();
    final_poss_   .clear();
         times_   .clear();

  selected_track_ids_.clear();
  recorded_el_electron_ids_.clear();
}
