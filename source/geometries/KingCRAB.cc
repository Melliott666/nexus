//
// Geometry of the King CRAB detector
//

#include "KingCRAB.h"

// Nexus geometry/utilities
#include "BoxPointSampler.h"
#include "CylinderPointSampler.h"
#include "FactoryBase.h"
#include "HexagonMeshTools.h"
#include "Visibilities.h"

// Nexus materials/optics/physics
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "UniformElectricDriftField.h"
#include "XenonProperties.h"

// Nexus sensitive detectors
#include "IonizationSD.h"
#include "SensorSD.h"

// Geant4 solids
#include <G4Box.hh>
#include <G4ExtrudedSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>

// Geant4 geometry/volumes
#include <G4LogicalVolume.hh>
#include <G4PVParameterised.hh>
#include <G4PVPlacement.hh>
#include <G4PVReplica.hh>
#include <G4Region.hh>

// Geant4 stores/managers
#include <G4LogicalVolumeStore.hh>
#include <G4NistManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4SDManager.hh>

// Geant4 UI/misc
#include <G4GenericMessenger.hh>
#include <G4UnitsTable.hh>
#include <G4UserLimits.hh>
#include <G4VisAttributes.hh>

// Geant4 optical surfaces
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>

// CLHEP
#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>


namespace nexus{
    using namespace CLHEP;
    REGISTER_CLASS(KingCRAB, GeometryBase)


    KingCRAB::KingCRAB():
        GeometryBase(),
        msg_(nullptr),
        Lab_size(550. * m),
        gas_pressure_(1. * bar),
        max_step_size_(0.1 * mm),
        gastype_("xenon"),
        specific_vertex_(0., 0., 0.),
        drift_field_on_(true),
        drift_v_(1.0 * mm/microsecond),
        drift_e_lifetime_(100. * ms),
        el_field_on_(true),
        el_field_int_(15.0 * kilovolt/cm),
        EL_drift_v_(2.5 * mm/microsecond),
        anode_gen_(nullptr),
        cathode_gen_(nullptr),
        active_volume_gen_(nullptr),
        practice_track_(nullptr)

    {
        msg_ = new G4GenericMessenger(this, "/Geometry/KingCRAB/","Control commands of geometry of KingCRAB TPC");
        G4GenericMessenger::Command&  Pressure_cmd =msg_->DeclarePropertyWithUnit("gas_pressure","bar",gas_pressure_,"Pressure of Gas");
        Pressure_cmd.SetParameterName("XenonPressure", false);

        G4GenericMessenger::Command&  step_size_cmd =msg_->DeclarePropertyWithUnit("max_step_size","mm",max_step_size_,"The maximum step size");
        step_size_cmd .SetParameterName("max_step_size", false);

        msg_->DeclareProperty("gastype", gastype_, "The GAS to use in the detector");

        msg_->DeclarePropertyWithUnit("specific_vertex_", "mm",  specific_vertex_, "Set generation vertex.");

        msg_->DeclareProperty("drift_field_on", drift_field_on_, "Turn drift field on/off.");

        G4GenericMessenger::Command& drift_v_cmd =msg_->DeclareProperty("drift_v", drift_v_,"The active region drift velocity.");
        drift_v_cmd.SetParameterName("drift_v", true);
        drift_v_cmd.SetRange("drift_v>=0.");

        G4GenericMessenger::Command& drift_lifetime_cmd =msg_->DeclareProperty("drift_e_lifetime", drift_e_lifetime_,"The electron lifetime in the drift region.");
        drift_lifetime_cmd.SetParameterName("drift_e_lifetime", true);
        drift_lifetime_cmd.SetUnitCategory("Time");

        msg_->DeclareProperty("el_field_on", el_field_on_, "Turn EL field on/off.");

        G4GenericMessenger::Command& el_field_cmd =msg_->DeclareProperty("EL_field_intensity", el_field_int_,"Electric field in the EL region.");
        el_field_cmd.SetUnitCategory("Electric field");
        el_field_cmd.SetParameterName("EL_field_intensity", true);

        G4GenericMessenger::Command& el_drift_cmd =msg_->DeclareProperty("EL_drift_v", EL_drift_v_,"The EL region drift velocity.");
        el_drift_cmd.SetParameterName("EL_drift_v", true);
        el_drift_cmd.SetRange("EL_drift_v>=0.");
    }


    KingCRAB::~KingCRAB()
    {
        delete msg_;
        delete anode_gen_;
        delete cathode_gen_;
        delete active_volume_gen_;
        delete practice_track_;
    }


    void KingCRAB::Construct(){

        // --------------------------
        // Materials
        // --------------------------
        G4Material *GAS;
        if (gastype_ == "xenon"){
            std::cout << "Using Xenon! "<< gas_pressure_/bar << " bar"  << std::endl;
            GAS = materials::GXeEnriched(gas_pressure_, 293. * kelvin);
            G4double sc_yield    = 25510. * 1/MeV;
            G4double e_lifetime  = 100*ms;
            GAS->SetMaterialPropertiesTable(opticalprops::GXe(gas_pressure_, 293*kelvin, sc_yield, e_lifetime));
        }
        else if (gastype_ == "argon"){
            std::cout << "Using Argon! "<< gas_pressure_/bar << " bar"  << std::endl;
            GAS = materials::GAr(gas_pressure_, 293. * kelvin);
            G4double sc_yield    = 25510. * 1/MeV; // NEEDS UPDATING Ws
            G4double e_lifetime  = 100*ms;
            GAS->SetMaterialPropertiesTable(opticalprops::GAr(sc_yield, e_lifetime));
        }
        else 
            std::cout << "Error in specified gas" << std::endl;
        
        G4Material *Steel=materials::Steel();
        Steel->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());

        G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
        G4Material* PTFE = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");


        // --------------------------
        // Optical Surfaces
        // --------------------------
        G4OpticalSurface* gas_steel_opsur = new G4OpticalSurface("GAS_STEEL_OPSURF", unified, ground, dielectric_metal);
        gas_steel_opsur->SetSigmaAlpha(0.0);
        gas_steel_opsur->SetMaterialPropertiesTable(opticalprops::NoRef()); // No Reflections
        // gas_steel_opsur->SetMaterialPropertiesTable(opticalprops::Steel()); // With Reflections

        G4OpticalSurface* gas_copper_opsur = new G4OpticalSurface("GAS_COPPER_OPSURF", unified, ground, dielectric_metal);
        gas_copper_opsur->SetSigmaAlpha(0.0);
        gas_copper_opsur->SetMaterialPropertiesTable(opticalprops::NoRef()); // No Reflections
        // gas_copper_opsur->SetMaterialPropertiesTable(opticalprops::Copper()); // With Reflections

        G4OpticalSurface* gas_poly_opsurf = new G4OpticalSurface("GAS_POLY_OPSURF", unified, ground, dielectric_dielectric, .01);
        gas_poly_opsurf->SetMaterialPropertiesTable(opticalprops::NoRef()); // No Reflections
        // gas_poly_opsurf->SetMaterialPropertiesTable(opticalprops::PTFE()); // With Reflections


        // --------------------------
        // Lab Space
        // --------------------------
        G4String lab_name="LAB";
        G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2.,Lab_size/2.,Lab_size/2.);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;
        auto labPhysical= new G4PVPlacement(0,G4ThreeVector(), lab_logic_volume, lab_logic_volume->GetName(), 0, false, 0, false);
    

        // --------------------------
        // Vessel, Flanges, Endcaps
        // --------------------------
        G4double vessel_OD      = 50.8*cm;
        G4double vessel_length  = 135.255*cm;
        G4double vessel_thickn  = 1.27*cm;

        G4double vessel_OR = vessel_OD/2.0;
        G4double vessel_IR = vessel_OR - vessel_thickn;

        G4double anode_pos = 13*cm; // Distance from anode-side vessel wall to EL center
        G4double z_shift   = vessel_length/2.0 - anode_pos; // Shifts so the center of the EL gap is z=0

        G4Tubs* vessel_solid = new G4Tubs("VESSEL", vessel_IR, vessel_OR, vessel_length/2.0, 0, twopi);
        G4LogicalVolume* vessel_logic = new G4LogicalVolume(vessel_solid, Steel, "VESSEL");
        G4VPhysicalVolume* vessel_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_shift), vessel_logic, "VESSEL", lab_logic_volume, false, 0, true);
        new G4LogicalSkinSurface("GAS_VESSEL_OPSURF", vessel_logic, gas_steel_opsur);

        G4double flange_diam  = 71.12*cm;
        G4double flange_thick = 3.81*cm;

        G4Tubs* flange_solid = new G4Tubs("FLANGE_RING", vessel_IR, flange_diam/2.0, flange_thick/2.0, 0, twopi);
        G4LogicalVolume* flange_logic = new G4LogicalVolume(flange_solid, Steel, "FLANGE_RING");

        G4double z_vessel_plus_end  = z_shift + vessel_length/2.0;
        G4double z_vessel_minus_end = z_shift - vessel_length/2.0;

        G4VPhysicalVolume* flange1_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_vessel_plus_end  + flange_thick/2.0), flange_logic, "FLANGE_RING_PLUS",  lab_logic_volume, false, 0, true);
        G4VPhysicalVolume* flange2_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_vessel_minus_end - flange_thick/2.0), flange_logic, "FLANGE_RING_MINUS", lab_logic_volume, false, 1, true);
        new G4LogicalSkinSurface("GAS_FLANGE_OPSURF", flange_logic, gas_steel_opsur);

        G4double endcap_thick = 3.81*cm;

        G4double detector_hole_diam = 2.0*2.54*cm;
        G4double detector_hole_rad  = detector_hole_diam/2.0;
        G4double detector_hole_xpos = -8.255*cm;
        G4double detector_hole_ypos = 14.3002*cm;

        G4Tubs* endcap_solid = new G4Tubs("ENDCAP", 0., flange_diam/2.0, endcap_thick/2.0, 0, twopi);
        G4Tubs* detector_hole_solid = new G4Tubs("DETECTOR_ENDCAP_HOLE", 0., detector_hole_rad, endcap_thick, 0, twopi);
        G4SubtractionSolid* endcap_plus_solid = new G4SubtractionSolid("ENDCAP_PLUS_WITH_HOLE", endcap_solid, detector_hole_solid, 0, G4ThreeVector(detector_hole_xpos, detector_hole_ypos, 0.));

        G4LogicalVolume* endcap_plus_logic  = new G4LogicalVolume(endcap_plus_solid, Steel, "ENDCAP_PLUS");
        G4LogicalVolume* endcap_minus_logic = new G4LogicalVolume(endcap_solid, Steel, "ENDCAP_MINUS");

        G4VPhysicalVolume* endcap1_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_vessel_plus_end + flange_thick + endcap_thick/2.0), endcap_plus_logic, "ENDCAP_PLUS", lab_logic_volume, false, 0, true);
        G4VPhysicalVolume* endcap2_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_vessel_minus_end - flange_thick - endcap_thick/2.0), endcap_minus_logic, "ENDCAP_MINUS", lab_logic_volume, false, 1, true);

        new G4LogicalSkinSurface("GAS_ENDCAP_PLUS_OPSURF",  endcap_plus_logic,  gas_steel_opsur);
        new G4LogicalSkinSurface("GAS_ENDCAP_MINUS_OPSURF", endcap_minus_logic, gas_steel_opsur);


        // --------------------------
        // Image Intensifier Region Parameters
        // --------------------------
        G4double II_OD     = 10.0*2.54*cm;
        G4double II_OR     = II_OD/2.0;
        G4double II_thickn = vessel_thickn;
        G4double II_IR     = II_OR - II_thickn;

        G4double II_length = 10.0*2.54*cm;
        G4double II_xpos   = detector_hole_xpos;
        G4double II_ypos   = detector_hole_ypos;

        G4double z_cathode_endcap_outer = z_vessel_plus_end + flange_thick + endcap_thick;
        G4double II_zpos = z_cathode_endcap_outer + II_length/2.0;


        // --------------------------
        // Gas Volume
        // --------------------------
        G4double gas_overlap = 0.1*mm;

        G4Tubs* gas_main_solid = new G4Tubs("GAS_MAIN", 0, vessel_IR, vessel_length/2.0 + flange_thick, 0, twopi);
        G4Tubs* detector_hole_gas_solid = new G4Tubs("DETECTOR_ENDCAP_HOLE_GAS_SOLID", 0., detector_hole_rad, endcap_thick/2.0 + gas_overlap, 0, twopi);

        G4double detector_hole_gas_zpos = z_vessel_plus_end + flange_thick + endcap_thick/2.0 - z_shift;

        G4UnionSolid* gas_plus_hole_solid = new G4UnionSolid("GAS_PLUS_DETECTOR_HOLE", gas_main_solid, detector_hole_gas_solid, 0, G4ThreeVector(detector_hole_xpos, detector_hole_ypos, detector_hole_gas_zpos));

        G4Tubs* II_gas_solid = new G4Tubs("II_REGION_GAS_SOLID", 0., II_IR, II_length/2.0 + gas_overlap, 0, twopi);
        G4double II_gas_zpos = II_zpos - z_shift;

        G4UnionSolid* gas_solid = new G4UnionSolid("GAS", gas_plus_hole_solid, II_gas_solid, 0, G4ThreeVector(II_xpos, II_ypos, II_gas_zpos));

        G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, GAS, "GAS");
        G4VPhysicalVolume * gas_phys= new G4PVPlacement(0, G4ThreeVector(0.,0., z_shift), gas_logic, gas_solid->GetName(), lab_logic_volume, false, 0, true);

        IonizationSD* gasSD = new IonizationSD("/KingCRAB/GAS");
        gas_logic->SetSensitiveDetector(gasSD);
        G4SDManager::GetSDMpointer()->AddNewDetector(gasSD);
        gas_logic->SetUserLimits(new G4UserLimits(max_step_size_));
        this->SetLogicalVolume(lab_logic_volume);


        // --------------------------
        // EL / Cathode / Drift Geometry Parameters
        // --------------------------
        G4double EL_ring_ID    = 382*mm;
        G4double EL_ring_OD    = 437*mm;
        G4double EL_ring_thick = 14*mm;
        G4double mesh_thick    = 130*um;
        G4double el_gap        = 7*mm;

        G4double cathode_ring_ID    = 321*mm;
        G4double cathode_ring_OD    = 373*mm;
        G4double cathode_ring_thick = 14*mm;

        G4double z_cathode_wall = +vessel_length/2.0;
        G4double cathode_inset_face = 39*cm;
        G4double z_cathode_mesh = z_cathode_wall - (cathode_inset_face + mesh_thick/2.0);
        G4double z_cathode_ring = z_cathode_mesh - (mesh_thick/2.0 + cathode_ring_thick/2.0);

        G4double z_anode_mesh = -el_gap/2. - mesh_thick/2. - z_shift;
        G4double z_gate_mesh  = +el_gap/2. + mesh_thick/2. - z_shift;

        G4double z_active_min = z_gate_mesh + mesh_thick/2.;
        G4double z_active_max = z_cathode_mesh - mesh_thick/2.;
        G4double active_length = z_active_max - z_active_min;
        G4double active_zpos   = 0.5 * (z_active_min + z_active_max);

        G4double el_gap_posZ = 0.5 * (z_anode_mesh + z_gate_mesh);


        // --------------------------
        // EL Mesh Logic
        // --------------------------
        G4double EL_mesh_diam = 2.5*mm;
        G4double hex_circumradius = EL_mesh_diam/std::sqrt(3)*mm;
        G4int n_hex = (G4int) ((EL_ring_ID/2.0) / hex_circumradius);

        G4Tubs* grid_solid = new G4Tubs("EL_GRID", 0., EL_ring_OD/2.0 , mesh_thick/2., 0., twopi);
        G4LogicalVolume* EL_grid_logic = new G4LogicalVolume(grid_solid, Steel, "EL_GRID");

        G4ExtrudedSolid* hex_prism = CreateHexagon(mesh_thick/2.0, hex_circumradius);
        G4LogicalVolume* EL_hex_logic  = new G4LogicalVolume(hex_prism, GAS, "MESH_HEX_GAS");

        PlaceHexagons(n_hex, EL_mesh_diam, mesh_thick, EL_grid_logic, EL_hex_logic, EL_ring_ID);
        new G4LogicalSkinSurface("GAS_EL_MESH_OPSURF", EL_grid_logic, gas_steel_opsur);

        CLHEP::HepRotationZ Roty(15*deg);
        G4RotationMatrix* pRot = new G4RotationMatrix();
        pRot->set(Roty);


        // --------------------------
        // Drift Field
        // CRAB0-style: attach the drift field to the full GAS logical volume.
        // The z positions still restrict the drift range, but the region root is GAS.
        // --------------------------
        if (drift_field_on_) {
            UniformElectricDriftField* drift_field = new UniformElectricDriftField();

            drift_field->SetCathodePosition(z_active_max);
            drift_field->SetAnodePosition(z_gate_mesh - mesh_thick/2.);
            drift_field->SetDriftVelocity(drift_v_);
            drift_field->SetLifetime(drift_e_lifetime_);

            G4Region* drift_region = new G4Region("DRIFT");
            drift_region->SetUserInformation(drift_field);
            drift_region->AddRootLogicalVolume(gas_logic);
        }


        // --------------------------
        // EL Gap Region
        // --------------------------
       G4Tubs* elgap_solid = new G4Tubs("EL_GAP", 0., EL_ring_ID/2., el_gap/2., 0., twopi);
        G4LogicalVolume* elgap_logic = new G4LogicalVolume(elgap_solid, GAS, "EL_GAP");

        new G4PVPlacement(0, G4ThreeVector(0., 0., el_gap_posZ), elgap_logic, "EL_GAP", gas_logic, false, 0, true);

        elgap_logic->SetUserLimits(new G4UserLimits(max_step_size_));

        if (el_field_on_) {
            G4double yield = XenonELLightYield(el_field_int_, gas_pressure_);

            UniformElectricDriftField* el_field = new UniformElectricDriftField();

            el_field->SetCathodePosition(z_gate_mesh - mesh_thick/2.);
            el_field->SetAnodePosition(z_anode_mesh + mesh_thick/2.);
            el_field->SetDriftVelocity(EL_drift_v_);
            el_field->SetLightYield(yield);

            G4Region* el_region = new G4Region("EL_GAP");
            el_region->SetUserInformation(el_field);
            el_region->AddRootLogicalVolume(elgap_logic);
        }


        // --------------------------
        // EL Mesh Placements
        // --------------------------
        new G4PVPlacement(0, G4ThreeVector(0., 0., z_anode_mesh), EL_grid_logic, "EL_MESH_ANODE", gas_logic, false, 0, false);
        new G4PVPlacement(pRot, G4ThreeVector(0., 0., z_gate_mesh), EL_grid_logic, "EL_MESH_GATE", gas_logic, false, 1, false);


        // --------------------------
        // Anode and EL Rings
        // --------------------------
        G4Tubs*          EL_solid    = new G4Tubs("RING_SOLID", EL_ring_ID/2., EL_ring_OD/2., EL_ring_thick/2. - mesh_thick/2., 0, twopi);
        G4LogicalVolume* anode_logic = new G4LogicalVolume(EL_solid, Steel, "ANODE_RING");
        G4LogicalVolume* EL_logic    = new G4LogicalVolume(EL_solid, Steel, "EL_RING");

        // new G4PVPlacement(0, G4ThreeVector(0., 0., -el_gap/2.0 - mesh_thick/2. - EL_ring_thick/2.0 - z_shift), anode_logic, "ANODE_RING", gas_logic, false, 0, true);
        // new G4PVPlacement(0, G4ThreeVector(0., 0., +el_gap/2.0 + mesh_thick/2. + EL_ring_thick/2.0 - z_shift), EL_logic, "EL_RING", gas_logic, false, 0, true);
        
        new G4LogicalSkinSurface("GAS_ANODE_OPSURF", anode_logic, gas_steel_opsur);
        new G4LogicalSkinSurface("GAS_EL_GATE_OPSURF", EL_logic, gas_steel_opsur);


        // --------------------------
        // Cathode Ring and Mesh
        // --------------------------
        G4Tubs* cathode_solid = new G4Tubs("CATHODE_RING_SOLID", cathode_ring_ID/2., cathode_ring_OD/2., cathode_ring_thick/2. - mesh_thick/2., 0, twopi);
        G4LogicalVolume* cathode_logic = new G4LogicalVolume(cathode_solid, Steel, "CATHODE_RING");

        // new G4PVPlacement(0, G4ThreeVector(0., 0., z_cathode_ring), cathode_logic, "CATHODE_RING", gas_logic, false, 0, true);

        G4int n_hex_cathode = (G4int)((cathode_ring_ID/2.0) / hex_circumradius);

        G4Tubs* cathode_grid_solid = new G4Tubs("CATHODE_GRID", 0., cathode_ring_OD/2.0, mesh_thick/2., 0, twopi);
        G4LogicalVolume* cathode_grid_logic = new G4LogicalVolume(cathode_grid_solid, Steel, "CATHODE_GRID");

        PlaceHexagons(n_hex_cathode, EL_mesh_diam, mesh_thick, cathode_grid_logic, EL_hex_logic, cathode_ring_ID);

        new G4LogicalSkinSurface("GAS_CATHODE_MESH_OPSURF", cathode_grid_logic, gas_steel_opsur);

        // new G4PVPlacement(0, G4ThreeVector(0., 0., z_cathode_mesh), cathode_grid_logic, "CATHODE_MESH", gas_logic, false, 0, false);


        // --------------------------
        // Poly Wrap, Field Cage, Staves
        // --------------------------
        G4double poly_inner_diam = 482.6*mm;
        G4double poly_length = 960*mm;
        G4double poly_thickn = 6.35*mm; 
        G4Tubs*  poly_solid = new G4Tubs("POLY_WRAP", poly_inner_diam/2.0, poly_inner_diam/2.0+poly_thickn, poly_length/2.0, 0, twopi);
        G4LogicalVolume* poly_logic = new G4LogicalVolume(poly_solid, Steel, "POLY_WRAP");

        // new G4PVPlacement(0, G4ThreeVector(0., 0., poly_length/2.0 - z_shift - el_gap/2. - mesh_thick/2. - EL_ring_thick), poly_logic, poly_solid->GetName(), gas_logic, false, 0, true);
        new G4LogicalSkinSurface("GAS_POLY_OPSURF", poly_logic, gas_poly_opsurf);

        G4double field_ring_OD = 430*mm;
        G4double field_ring_thickn = 5*mm;
        G4Tubs*  field_ring_solid = new G4Tubs("FIELD_RING", field_ring_OD/2.0 - field_ring_thickn, field_ring_OD/2.0, field_ring_thickn/2.0, 0, twopi);
        G4LogicalVolume* field_ring_logic = new G4LogicalVolume(field_ring_solid, Steel, "FIELD_RING");
        
        G4double field_ring_origin = 46*mm + EL_ring_thick + mesh_thick/2. + el_gap/2. - z_shift; 
        G4double active_ring_sep = 24*mm;
        G4double posz;
        G4double n_active_FR = 31;

        // for (G4int i=0; i<n_active_FR; i++) {posz = field_ring_origin + i*active_ring_sep; new G4PVPlacement(0, G4ThreeVector(0., 0., posz), field_ring_logic, field_ring_logic->GetName(), gas_logic, false, i, true);}

        G4double active_buffer_FR_dist = 72*mm;
        G4double buffer_ring_sep = 48*mm;

        // for (G4int i=0; i<4; i++) {posz = field_ring_origin + 30*active_ring_sep + active_buffer_FR_dist + i*buffer_ring_sep; new G4PVPlacement(0, G4ThreeVector(0., 0., posz), field_ring_logic, "FIELD_RING", gas_logic, false, i+31, false);}

        new G4LogicalSkinSurface("GAS_FIELDCAGE_OPSURF", field_ring_logic, gas_copper_opsur);

        G4double Stave_ID = 6*mm;
        G4double Stave_OD = 44.7*mm;
        G4double Stave_thick = (Stave_OD - Stave_ID)/2.0;
        G4int n_Staves = 11;

        G4double Stave_rpos = field_ring_OD/2.0 + Stave_OD/2.0 + 5*mm;
        G4double z_first_FR = field_ring_origin;
        G4double z_last_FR  = field_ring_origin + 30*active_ring_sep + active_buffer_FR_dist + 3*buffer_ring_sep;
        G4double Stave_length = (z_last_FR - z_first_FR); 
        G4double Stave_zpos   = (z_first_FR + z_last_FR)/2.0;

        G4Tubs*  Stave_solid = new G4Tubs("STAVES", Stave_ID/2.0, Stave_OD/2.0, Stave_length/2.0, 0, twopi);
        G4LogicalVolume* Stave_logic = new G4LogicalVolume(Stave_solid, PTFE, "STAVES");

        // for (G4int i=0; i<n_Staves; i++) {G4double phi = i*twopi/n_Staves; G4double xt = Stave_rpos*std::cos(phi); G4double yt = Stave_rpos*std::sin(phi); new G4PVPlacement(0, G4ThreeVector(xt, yt, Stave_zpos), Stave_logic, Stave_solid->GetName(), gas_logic, false, i, true);}


        // --------------------------
        // Fused Silica Plano-Convex Lens
        // --------------------------
        G4Material* fs_mat = materials::FusedSilica();
        fs_mat->SetMaterialPropertiesTable(opticalprops::FusedSilica());

        G4double Lens_D  = 50.8*mm;
        G4double Lens_R  = 34.39*mm;
        G4double Lens_tc = 12.5*mm;
        G4double Lens_a  = Lens_D/2.0;

        G4double Lens_zc_shift = Lens_tc/2.0 - Lens_R;

        G4Tubs*   Lens_blank  = new G4Tubs("FS_LENS_BLANK",  0., Lens_a, Lens_tc/2.0, 0., twopi);
        G4Sphere* Lens_sphere = new G4Sphere("FS_LENS_SPHERE", 0., Lens_R, 0., twopi, 0., pi);
        G4Tubs*   Lens_big    = new G4Tubs("FS_LENS_BIG",    0., Lens_a, Lens_R,      0., twopi);

        G4SubtractionSolid* Lens_outside = new G4SubtractionSolid("FS_LENS_OUTSIDE", Lens_big, Lens_sphere, 0, G4ThreeVector(0.,0.,Lens_zc_shift));
        G4SubtractionSolid* Lens_solid   = new G4SubtractionSolid("FS_LENS_SOLID",   Lens_blank, Lens_outside, 0, G4ThreeVector(0.,0.,0.));

        G4LogicalVolume* Lens_logic = new G4LogicalVolume(Lens_solid, fs_mat, "FS_LENS");

        G4double z_endcap_wall = vessel_length/2.0 + flange_thick;
        G4double lens_offset = 11.75*cm;
        G4double Lens_zpos = z_endcap_wall - (lens_offset + Lens_tc/2.0);

        G4RotationMatrix* Lens_rot = new G4RotationMatrix();
        Lens_rot->rotateY(180.0*deg);

        new G4PVPlacement(Lens_rot, G4ThreeVector(0., 0., Lens_zpos), Lens_logic, "FS_LENS", gas_logic, false, 0, true);


        // --------------------------
        // Mirrors
        // --------------------------
        G4double Mirror_D = 50.8*mm;
        G4double Mirror_T = .01*mm;
        G4double Mirror_zpos = Lens_zpos + (6.67*cm);

        G4Tubs* Mirror_solid = new G4Tubs("MIRROR", 0., Mirror_D/2.0, Mirror_T/2.0, 0., twopi);
        G4LogicalVolume* Mirror_logic = new G4LogicalVolume(Mirror_solid, Steel, "MIRROR");

        G4double theta = -29.9963208064*deg;

        G4RotationMatrix* Mirror_rot = new G4RotationMatrix();
        Mirror_rot->rotateZ(theta);
        Mirror_rot->rotateX(135.0*deg);

        new G4PVPlacement(Mirror_rot, G4ThreeVector(0., 0., Mirror_zpos), Mirror_logic, "MIRROR", gas_logic, false, 0, true);

        auto* mirror1_opsur = new G4OpticalSurface("MIRROR1_OPSURF", unified, polished, dielectric_metal);
        mirror1_opsur->SetMaterialPropertiesTable(opticalprops::MirrorReflectivity());
        new G4LogicalSkinSurface("GAS_MIRROR1_SKIN", Mirror_logic, mirror1_opsur);

        G4double Mirror2_xpos = detector_hole_xpos;
        G4double Mirror2_ypos = detector_hole_ypos;
        G4double Mirror2_zpos = Mirror_zpos;

        G4Tubs* Mirror2_solid = new G4Tubs("MIRROR2", 0., Mirror_D/2.0, Mirror_T/2.0, 0., twopi);
        G4LogicalVolume* Mirror2_logic = new G4LogicalVolume(Mirror2_solid, Steel, "MIRROR2");

        G4RotationMatrix* Mirror2_rot = new G4RotationMatrix();
        Mirror2_rot->rotateZ(theta);
        Mirror2_rot->rotateX(-45.0*deg);

        new G4PVPlacement(Mirror2_rot, G4ThreeVector(Mirror2_xpos, Mirror2_ypos, Mirror2_zpos), Mirror2_logic, "MIRROR2", gas_logic, false, 0, true);

        auto* mirror2_opsur = new G4OpticalSurface("MIRROR2_OPSURF", unified, polished, dielectric_metal);
        mirror2_opsur->SetMaterialPropertiesTable(opticalprops::MirrorReflectivity());
        new G4LogicalSkinSurface("GAS_MIRROR2_SKIN", Mirror2_logic, mirror2_opsur);


        // --------------------------
        // Detector and Image Intensifier Lens
        // --------------------------
        G4double lens_diam = 5.08*cm;
        G4double lens_thick = 2*mm;

        G4Tubs* lens_solid = new G4Tubs("LENS", 0, lens_diam/2.0, lens_thick/2.0, 0, twopi);
        G4LogicalVolume* lens_logic = new G4LogicalVolume(lens_solid, Steel, "LENS");
        
        G4OpticalSurface* lens_opsur = new G4OpticalSurface("LENS", unified, polished, dielectric_metal);
        lens_opsur->SetMaterialPropertiesTable(opticalprops::Absorber());
        new G4LogicalSkinSurface("LENS", lens_logic, lens_opsur);

        G4RotationMatrix* detRot = nullptr;
        G4double det_zpos = 1457.4*mm - z_shift;
        G4ThreeVector det_pos(II_xpos, II_ypos, det_zpos);

        new G4PVPlacement(detRot, det_pos, lens_logic, lens_solid->GetName(), gas_logic, false, 0, true);

        G4double II_lens_from_detector = 68.405974*mm;
        G4double II_lens_zpos = det_zpos - II_lens_from_detector;

        G4RotationMatrix* II_Lens_rot = new G4RotationMatrix();
        II_Lens_rot->rotateY(180.0*deg);

        new G4PVPlacement(II_Lens_rot, G4ThreeVector(II_xpos, II_ypos, II_lens_zpos), Lens_logic, "II_FS_LENS", gas_logic, false, 1, true);


        // ------------------------
        // Needle
        // ------------------------
        G4double needle_diam = 0.7*mm;
        G4double needle_length = 2.54*cm;
        G4double needle_xpos = 13.97*cm;
        G4Tubs* needle_solid = new G4Tubs("NEEDLE", 0, needle_diam/2.0, needle_length/2.0, 0, twopi);
        G4LogicalVolume* needle_logic  = new G4LogicalVolume(needle_solid, Steel, "NEEDLE_BODY");

        // new G4PVPlacement(0,G4ThreeVector(needle_xpos/2, 0., z_cathode_ring-needle_length/4), needle_logic, "NEEDLE_BODY", gas_logic, false, 0, true);


        // --------------------------
        // Image Intensifier Region Shell
        // --------------------------
        G4Tubs* II_shell_solid = new G4Tubs("II_REGION_SHELL", II_IR, II_OR, II_length/2.0, 0, twopi);
        G4LogicalVolume* II_shell_logic = new G4LogicalVolume(II_shell_solid, Steel, "II_REGION_SHELL");

        new G4PVPlacement(0, G4ThreeVector(II_xpos, II_ypos, II_zpos), II_shell_logic, "II_REGION_SHELL", lab_logic_volume, false, 0, true);
        new G4LogicalSkinSurface("GAS_II_REGION_SHELL_OPSURF", II_shell_logic, gas_steel_opsur);

        G4double II_endcap_thick = vessel_thickn;
        G4Tubs* II_endcap_solid = new G4Tubs("II_REGION_ENDCAP", 0., II_OR, II_endcap_thick/2.0, 0, twopi);
        G4LogicalVolume* II_endcap_logic = new G4LogicalVolume(II_endcap_solid, Steel, "II_REGION_ENDCAP");

        G4double II_endcap_zpos = II_zpos + II_length/2.0 + II_endcap_thick/2.0;

        new G4PVPlacement(0, G4ThreeVector(II_xpos, II_ypos, II_endcap_zpos), II_endcap_logic, "II_REGION_ENDCAP", lab_logic_volume, false, 0, true);
        new G4LogicalSkinSurface("GAS_II_REGION_ENDCAP_OPSURF", II_endcap_logic, gas_steel_opsur);


        // --------------------------
        // Vertex Generators 
        // --------------------------
        anode_gen_ = new CylinderPointSampler(4*mm, EL_ring_ID/2, 0.5*um, 0., twopi, nullptr, G4ThreeVector(0.,0.,el_gap_posZ));
        cathode_gen_ = new CylinderPointSampler(0.0, cathode_ring_ID/2.0, 0.5*um, 0.0, twopi, nullptr, G4ThreeVector(0., 0., z_cathode_mesh));
        active_volume_gen_ = new CylinderPointSampler(0.0, cathode_ring_ID/2.0, active_length/2.0, 0.0, twopi, nullptr, G4ThreeVector(0., 0., active_zpos));
        practice_track_ = new BoxPointSampler(0, 0, 0 , 0);

       
        // --------------------------
        // Visuals 
        // --------------------------
        AssignVisuals();

    }


    void KingCRAB::AssignVisuals() {
        G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

        G4LogicalVolume* Lab = lvStore->GetVolume("LAB");
        Lab->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4VisAttributes *VesselVa=new G4VisAttributes(nexus::TitaniumGreyAlpha());
        VesselVa->SetForceSolid(true);
        
        G4LogicalVolume* Vessel = lvStore->GetVolume("VESSEL");
        Vessel->SetVisAttributes(VesselVa);

        G4LogicalVolume* FlangeRing = lvStore->GetVolume("FLANGE_RING");
        if (FlangeRing) FlangeRing->SetVisAttributes(VesselVa);

        G4LogicalVolume* EndcapPlus = lvStore->GetVolume("ENDCAP_PLUS");
        if (EndcapPlus) EndcapPlus->SetVisAttributes(VesselVa);

        G4LogicalVolume* EndcapMinus = lvStore->GetVolume("ENDCAP_MINUS");
        if (EndcapMinus) EndcapMinus->SetVisAttributes(VesselVa);

        G4LogicalVolume* EL_Gate = lvStore->GetVolume("EL_RING");
        if (EL_Gate) EL_Gate->SetVisAttributes(VesselVa);
        
        G4LogicalVolume* Anode = lvStore->GetVolume("ANODE_RING");
        if (Anode) Anode->SetVisAttributes(VesselVa);

        G4VisAttributes *PolyWrapVa=new G4VisAttributes(nexus::WhiteAlpha());
        PolyWrapVa->SetForceSolid(true);
        G4LogicalVolume* PolyWrap = lvStore->GetVolume("POLY_WRAP");
        if (PolyWrap) PolyWrap->SetVisAttributes(PolyWrapVa);

        G4VisAttributes *LensVa=new G4VisAttributes(nexus::Red());
        LensVa->SetForceSolid(true);
        G4LogicalVolume* Lens = lvStore->GetVolume("LENS");
        if (Lens) Lens->SetVisAttributes(LensVa);

        G4VisAttributes *FieldRingVa=new G4VisAttributes(nexus::CopperBrownAlpha());
        FieldRingVa->SetForceSolid(true);
        G4LogicalVolume* FieldRing = lvStore->GetVolume("FIELD_RING");
        if (FieldRing) FieldRing->SetVisAttributes(FieldRingVa);

        G4LogicalVolume* Gas = lvStore->GetVolume("GAS");
        G4VisAttributes *GasVa=new G4VisAttributes(nexus::BlueAlpha());
        GasVa->SetForceCloud(true);
        if (Gas) Gas->SetVisAttributes(GasVa);

        G4LogicalVolume* Active = lvStore->GetVolume("ACTIVE");
        if (Active) Active->SetVisAttributes(GasVa);

        G4LogicalVolume* ELGap = lvStore->GetVolume("EL_GAP");
        if (ELGap) ELGap->SetVisAttributes(GasVa);

        G4VisAttributes *StaveVa =new G4VisAttributes(nexus::WhiteAlpha());
        StaveVa-> SetForceSolid(true);
        G4LogicalVolume* Tube = lvStore->GetVolume("STAVES");
        if (Tube) Tube->SetVisAttributes(StaveVa);

        G4VisAttributes* fsVa = new G4VisAttributes(nexus::Blue());
        fsVa->SetForceSolid(true);
        G4LogicalVolume* fsLV = lvStore->GetVolume("FS_LENS");
        if (fsLV) fsLV->SetVisAttributes(fsVa);

        G4VisAttributes* MirrorVa = new G4VisAttributes(nexus::Blue());
        MirrorVa->SetForceSolid(true);

        G4LogicalVolume* Mirror1 = lvStore->GetVolume("MIRROR");
        if (Mirror1) Mirror1->SetVisAttributes(MirrorVa);

        G4LogicalVolume* Mirror2 = lvStore->GetVolume("MIRROR2");
        if (Mirror2) Mirror2->SetVisAttributes(MirrorVa);

        auto* ELVa = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1.0));
        ELVa->SetForceSolid(true);
        if (EL_Gate) EL_Gate->SetVisAttributes(ELVa);

        auto* AnodeVa = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1.0));
        AnodeVa->SetForceSolid(true);
        if (Anode) Anode->SetVisAttributes(AnodeVa);

        G4VisAttributes* CathodeVa = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
        CathodeVa->SetForceSolid(true);

        G4LogicalVolume* CathodeRingLV = lvStore->GetVolume("CATHODE_RING");
        if (CathodeRingLV) CathodeRingLV->SetVisAttributes(CathodeVa);

        G4LogicalVolume* CathodeMeshLV = lvStore->GetVolume("CATHODE_GRID");
        if (CathodeMeshLV) CathodeMeshLV->SetVisAttributes(CathodeVa);

        G4VisAttributes* NeedleVa = new G4VisAttributes(G4Colour(0, 1.0, 2.0));
        NeedleVa->SetForceSolid(true);

        G4LogicalVolume* needleLV  = lvStore->GetVolume("NEEDLE_BODY");
        if (needleLV) needleLV->SetVisAttributes(NeedleVa);

        G4LogicalVolume* IIRegionShell = lvStore->GetVolume("II_REGION_SHELL");
        if (IIRegionShell) IIRegionShell->SetVisAttributes(VesselVa);

        G4LogicalVolume* IIRegionEndcap = lvStore->GetVolume("II_REGION_ENDCAP");
        if (IIRegionEndcap) IIRegionEndcap->SetVisAttributes(VesselVa);
    }


    void KingCRAB::PrintParam() {
    }


    G4ThreeVector KingCRAB::GenerateVertex(const G4String& region) const
    {
            G4ThreeVector pos;
            if (region == "AD_HOC") {
                pos = specific_vertex_;
            }
            else if((region == "ANODE")){
                pos = anode_gen_->GenerateVertex(VOLUME);
            }
            else if((region == "CATHODE")){
                pos = cathode_gen_->GenerateVertex(VOLUME);
            }
            else if((region == "ACTIVE")){
                pos = active_volume_gen_->GenerateVertex(VOLUME);
            }
            else if (region == "PRACTICE") {
                G4double t = G4UniformRand();

                G4double Lx = 10.0 * mm;
                G4double Ay = 2.0 * mm;
                G4double nwiggles = 3.0;

                G4double x = -0.5 * Lx + Lx * t;
                G4double y = Ay * std::sin(2.0 * pi * nwiggles * t);
                G4double z = 0.0;

                pos = G4ThreeVector(x, y, z);
            }
            else {
                G4Exception("[KingCRAB]", "GenerateVertex()", JustWarning, "Unknown vertex generation region. Setting default region as 0,0,0");
                pos=specific_vertex_;
            }
        return pos;
    }
}
