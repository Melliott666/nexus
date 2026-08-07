// ----------------------------------------------------------------------------
// nexus | CosmogenicMuonGenerator.cc
//
// Generate sea-level muons from the modified Guan--Gaisser spectrum and
// condition them to intersect a finite cylinder whose axis is global z.
// Global -y is vertically downward, matching the Nexus muon convention.
// ----------------------------------------------------------------------------

#include "CosmogenicMuonGenerator.h"
#include "FactoryBase.h"

#include <G4Event.hh>
#include <G4GenericMessenger.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

using namespace CLHEP;

namespace nexus {

REGISTER_CLASS(CosmogenicMuonGenerator, G4VPrimaryGenerator)

CosmogenicMuonGenerator::CosmogenicMuonGenerator():
  msg_(nullptr), spectrum_sampler_(nullptr), initialized_(false),
  energy_min_(10.*MeV), energy_max_(100.*TeV),
  positive_fraction_(0.56), target_radius_(160.*mm),
  target_length_(840.*mm), target_center_(0., 0., 424.474*mm),
  generation_distance_(800.*mm), paired_event_seeds_(true), base_seed_(410000),
  theta_bins_(180), energy_bins_(240)
{
  msg_ = new G4GenericMessenger(
    this, "/Generator/CosmogenicMuon/",
    "Bare-sky Guan--Gaisser muons crossing a finite cylinder.");

  auto& min_energy = msg_->DeclareProperty(
    "min_energy", energy_min_, "Minimum surface-muon kinetic energy.");
  min_energy.SetUnitCategory("Energy");
  min_energy.SetParameterName("min_energy", false);
  min_energy.SetRange("min_energy>0.");

  auto& max_energy = msg_->DeclareProperty(
    "max_energy", energy_max_, "Maximum surface-muon kinetic energy.");
  max_energy.SetUnitCategory("Energy");
  max_energy.SetParameterName("max_energy", false);
  max_energy.SetRange("max_energy>0.");

  msg_->DeclareProperty(
    "positive_fraction", positive_fraction_,
    "Probability of generating mu+ rather than mu-.");

  auto& radius = msg_->DeclareProperty(
    "target_radius", target_radius_, "Fiducial-cylinder radius.");
  radius.SetUnitCategory("Length");
  radius.SetParameterName("target_radius", false);
  radius.SetRange("target_radius>0.");

  auto& length = msg_->DeclareProperty(
    "target_length", target_length_, "Fiducial-cylinder length along global z.");
  length.SetUnitCategory("Length");
  length.SetParameterName("target_length", false);
  length.SetRange("target_length>0.");

  msg_->DeclarePropertyWithUnit(
    "target_center", "mm", target_center_, "Fiducial-cylinder global center.");

  auto& distance = msg_->DeclareProperty(
    "generation_distance", generation_distance_,
    "Upstream distance of the generation disk from the target center.");
  distance.SetUnitCategory("Length");
  distance.SetParameterName("generation_distance", false);
  distance.SetRange("generation_distance>0.");

  msg_->DeclareProperty(
    "paired_event_seeds", paired_event_seeds_,
    "Reset the random stream by event ID for paired field-on/off primaries.");
  msg_->DeclareProperty(
    "base_seed", base_seed_, "Base seed used when paired_event_seeds is true.");
}

CosmogenicMuonGenerator::~CosmogenicMuonGenerator()
{
  delete spectrum_sampler_;
  delete msg_;
}


G4double CosmogenicMuonGenerator::GuanIntensity(G4double energy_gev,
                                                G4double theta) const
{
  constexpr G4double p1 = 0.102573;
  constexpr G4double p2 = -0.068287;
  constexpr G4double p3 = 0.958633;
  constexpr G4double p4 = 0.0407253;
  constexpr G4double p5 = 0.817285;

  G4double cos_theta = std::max(0., std::cos(theta));
  G4double numerator = cos_theta*cos_theta + p1*p1
                     + p2*std::pow(cos_theta, p3)
                     + p4*std::pow(cos_theta, p5);
  G4double denominator = 1. + p1*p1 + p2 + p4;
  G4double cos_star = std::sqrt(numerator / denominator);
  G4double corrected_energy = energy_gev *
    (1. + 3.64 / (energy_gev * std::pow(cos_star, 1.29)));
  G4double pion = 1. / (1. + 1.1*energy_gev*cos_star/115.);
  G4double kaon = 0.054 / (1. + 1.1*energy_gev*cos_star/850.);
  return 0.14 * std::pow(corrected_energy, -2.7) * (pion + kaon);
}


void CosmogenicMuonGenerator::InitializeSpectrum()
{
  if (energy_max_ <= energy_min_)
    G4Exception("[CosmogenicMuonGenerator]", "InitializeSpectrum()",
                FatalException, "max_energy must be greater than min_energy");
  if (positive_fraction_ < 0. || positive_fraction_ > 1.)
    G4Exception("[CosmogenicMuonGenerator]", "InitializeSpectrum()",
                FatalException, "positive_fraction must lie in [0, 1]");
  G4double bound = std::sqrt(target_radius_*target_radius_
                           + target_length_*target_length_/4.);
  if (generation_distance_ <= bound)
    G4Exception("[CosmogenicMuonGenerator]", "InitializeSpectrum()",
                FatalException,
                "generation_distance must exceed the target bounding radius");

  std::vector<G4double> weights(theta_bins_ * energy_bins_);
  G4double log_min = std::log(energy_min_ / GeV);
  G4double log_max = std::log(energy_max_ / GeV);
  for (G4int itheta = 0; itheta < theta_bins_; ++itheta) {
    G4double theta = (itheta + 0.5) * halfpi / theta_bins_;
    for (G4int ienergy = 0; ienergy < energy_bins_; ++ienergy) {
      G4double log_energy = log_min
        + (ienergy + 0.5) * (log_max - log_min) / energy_bins_;
      G4double energy_gev = std::exp(log_energy);
      // The table is uniform in theta and ln(E): dE = E dln(E), while
      // dOmega contributes sin(theta). Projected-area weighting is supplied
      // exactly by rejecting lines that miss the target below.
      weights[itheta*energy_bins_ + ienergy] =
        GuanIntensity(energy_gev, theta) * energy_gev * std::sin(theta);
    }
  }
  spectrum_sampler_ = new G4RandGeneral(weights.data(), weights.size());
  initialized_ = true;
}


void CosmogenicMuonGenerator::SampleEnergyAndDirection(
  G4double& kinetic_energy, G4ThreeVector& direction)
{
  G4double sampled = spectrum_sampler_->fire();
  G4int index = std::min(
    static_cast<G4int>(sampled * theta_bins_ * energy_bins_),
    theta_bins_ * energy_bins_ - 1);
  G4int itheta = index / energy_bins_;
  G4int ienergy = index % energy_bins_;

  G4double theta = (itheta + G4UniformRand()) * halfpi / theta_bins_;
  G4double log_min = std::log(energy_min_ / GeV);
  G4double log_max = std::log(energy_max_ / GeV);
  G4double log_energy = log_min + (ienergy + G4UniformRand())
    * (log_max - log_min) / energy_bins_;
  kinetic_energy = std::exp(log_energy) * GeV;

  G4double phi = twopi * G4UniformRand();
  direction = G4ThreeVector(std::sin(theta)*std::cos(phi),
                            -std::cos(theta),
                            std::sin(theta)*std::sin(phi));
}


G4bool CosmogenicMuonGenerator::IntersectsTarget(
  const G4ThreeVector& vertex, const G4ThreeVector& direction) const
{
  G4ThreeVector p = vertex - target_center_;
  G4double a = direction.x()*direction.x() + direction.y()*direction.y();
  G4double b = 2. * (p.x()*direction.x() + p.y()*direction.y());
  G4double c = p.x()*p.x() + p.y()*p.y()
             - target_radius_*target_radius_;

  G4double radial_min = 0.;
  G4double radial_max = DBL_MAX;
  if (a < 1.e-15) {
    if (c > 0.) return false;
  } else {
    G4double discriminant = b*b - 4.*a*c;
    if (discriminant < 0.) return false;
    G4double root = std::sqrt(discriminant);
    radial_min = (-b - root) / (2.*a);
    radial_max = (-b + root) / (2.*a);
  }

  G4double axial_min = 0.;
  G4double axial_max = DBL_MAX;
  G4double half_length = target_length_/2.;
  if (std::abs(direction.z()) < 1.e-15) {
    if (std::abs(p.z()) > half_length) return false;
  } else {
    G4double t1 = (-half_length - p.z()) / direction.z();
    G4double t2 = (+half_length - p.z()) / direction.z();
    axial_min = std::min(t1, t2);
    axial_max = std::max(t1, t2);
  }

  return std::max({0., radial_min, axial_min})
       <= std::min(radial_max, axial_max);
}


G4bool CosmogenicMuonGenerator::SampleVertex(
  const G4ThreeVector& direction, G4ThreeVector& vertex) const
{
  G4double bound = std::sqrt(target_radius_*target_radius_
                           + target_length_*target_length_/4.);
  G4ThreeVector helper = std::abs(direction.z()) < 0.9
    ? G4ThreeVector(0., 0., 1.) : G4ThreeVector(1., 0., 0.);
  G4ThreeVector basis_u = direction.cross(helper).unit();
  G4ThreeVector basis_v = direction.cross(basis_u).unit();

  G4double radius = bound * std::sqrt(G4UniformRand());
  G4double angle = twopi * G4UniformRand();
  G4ThreeVector offset = radius *
    (std::cos(angle)*basis_u + std::sin(angle)*basis_v);
  vertex = target_center_ - generation_distance_*direction + offset;
  return IntersectsTarget(vertex, direction);
}


void CosmogenicMuonGenerator::GeneratePrimaryVertex(G4Event* event)
{
  if (paired_event_seeds_)
    G4Random::setTheSeed(base_seed_ + event->GetEventID());
  if (!initialized_) InitializeSpectrum();

  G4double kinetic_energy = 0.;
  G4ThreeVector direction;
  G4ThreeVector vertex;
  G4bool accepted = false;
  do {
    SampleEnergyAndDirection(kinetic_energy, direction);
    accepted = SampleVertex(direction, vertex);
    // A direction is accepted in proportion to the cylinder silhouette it
    // presents. This supplies the A_perp(theta, phi) factor in the rate.
  } while (!accepted);

  const char* particle_name =
    G4UniformRand() < positive_fraction_ ? "mu+" : "mu-";
  G4ParticleDefinition* definition =
    G4ParticleTable::GetParticleTable()->FindParticle(particle_name);
  if (!definition)
    G4Exception("[CosmogenicMuonGenerator]", "GeneratePrimaryVertex()",
                FatalException, "muon particle definition is unavailable");

  G4double mass = definition->GetPDGMass();
  G4double total_energy = kinetic_energy + mass;
  G4double momentum = std::sqrt(total_energy*total_energy - mass*mass);
  G4ThreeVector p = momentum * direction;

  auto* primary_vertex = new G4PrimaryVertex(vertex, 0.);
  primary_vertex->SetPrimary(
    new G4PrimaryParticle(definition, p.x(), p.y(), p.z()));
  event->AddPrimaryVertex(primary_vertex);
}

} // namespace nexus
