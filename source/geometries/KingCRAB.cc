//
// Geometry of the King CRAB detector
//

#include "KingCRAB.h"
#include "Visibilities.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "IonizationSD.h"
#include "FactoryBase.h"
#include <G4SubtractionSolid.hh>
#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4PVReplica.hh>
#include <G4SDManager.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "XenonProperties.h"
#include "G4UnitsTable.hh"
#include "G4ExtrudedSolid.hh"
#include "SensorSD.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <G4UserLimits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include "CylinderPointSampler.h"


namespace nexus{
    using namespace CLHEP;
    REGISTER_CLASS(KingCRAB, GeometryBase)


             KingCRAB::KingCRAB():
             GeometryBase(),
             msg_(nullptr),
             Lab_size(550. * m),
             vtx_(0,0,0), anode_gen_(nullptr),
             max_step_size_(0.1*mm), 
             gas_pressure_(1. * bar),
             gastype_("xenon")

    {
        msg_ = new G4GenericMessenger(this, "/Geometry/KingCRAB/","Control commands of geometry of KingCRAB TPC");
        G4GenericMessenger::Command&  Pressure_cmd =msg_->DeclarePropertyWithUnit("gas_pressure","bar",gas_pressure_,"Pressure of Gas");
        Pressure_cmd.SetParameterName("XenonPressure", false);

        G4GenericMessenger::Command&  step_size_cmd =msg_->DeclarePropertyWithUnit("max_step_size","mm",max_step_size_,"The maximum step size");
        step_size_cmd .SetParameterName("max_step_size", false);

        msg_->DeclareProperty("gastype", gastype_, "The GAS to use in the detector");

    }

    KingCRAB::~KingCRAB()
    {

        delete msg_;
        delete anode_gen_;
    }


    void KingCRAB::Construct(){

        // Materials
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

        // G4Material *Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
        // G4Material* PTFE = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");

        //Constructing Lab Space
        G4String lab_name="LAB";
        G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2.,Lab_size/2.,Lab_size/2.);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;
    
        // Vessel cylinder
        // Inner diameter = 19.25" = 48.9 cm, Vessel Length = 56.12" = 142.5 cm, vessel thick
        G4double vessel_diam   = 48.895*cm;
        G4double vessel_length = 142.5448*cm;
        G4double vessel_thickn = 2.54*cm;
        G4Tubs* vessel_solid = new G4Tubs("VESSEL", vessel_diam/2.0, vessel_diam/2.0+vessel_thickn, vessel_length/2.0, 0, twopi);
        G4LogicalVolume* vessel_logic = new G4LogicalVolume(vessel_solid, Steel, "VESSEL");

        // Flanges
        // Diameter = 28" = 71.12 cm, thickness 1"=2.54cm
        G4double flange_diam  = 71.12*cm;
        G4double flange_thick = 2.54*cm;
        G4Tubs*  flange_solid = new G4Tubs("FLANGE", 0, flange_diam/2.0, flange_thick/2.0, 0, twopi);
        G4LogicalVolume* flange_logic = new G4LogicalVolume(flange_solid, Steel, "FLANGE");

        // Field Cage -- Model as PTFE for now for simplicity
        G4double field_cage_diam = 34.3*cm;
        G4double field_cage_length = 142.5448*cm; // Use same as vessel length for now
        G4double field_cage_thickn = 1*cm; // Placeholder
        G4Tubs* field_cage_solid = new G4Tubs("FIELD_CAGE", field_cage_diam/2.0, field_cage_diam/2.0+field_cage_thickn, field_cage_length/2.0, 0, twopi);
        G4LogicalVolume* field_cage_logic = new G4LogicalVolume(field_cage_solid, Steel, "FIELD_CAGE"); // Set as steel for now so we can change reflectivity

        // Add optical surface to the field cage
        G4OpticalSurface* gas_FC_opsur = new G4OpticalSurface("GAS_FC_OPSURF");
        gas_FC_opsur->SetType(dielectric_metal);
        gas_FC_opsur->SetModel(unified);
        gas_FC_opsur->SetFinish(ground);
        gas_FC_opsur->SetSigmaAlpha(0.0);
        gas_FC_opsur->SetMaterialPropertiesTable(opticalprops::NoRef()); // No Reflections
        // gas_FC_opsur->SetMaterialPropertiesTable(opticalprops::Steel()); // With Reflections
        
        new G4LogicalSkinSurface("GAS_FIELDCAGE_OPSURF", field_cage_logic, gas_FC_opsur);
        new G4LogicalSkinSurface("GAS_FIELDCAGE_OPSURF", flange_logic, gas_FC_opsur);

        // Gas Volume
        G4Tubs* gas_solid = new G4Tubs("GAS", 0, vessel_diam/2., vessel_length/2., 0, twopi);
        G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, GAS, "GAS");

        // Lens
        G4double lens_diam = 2.54*cm; // 1 inch
        G4double lens_thick = 3*mm; // Choose arbitary for now
        G4Tubs* lens_solid = new G4Tubs("LENS", 0, lens_diam/2.0, lens_thick/2.0, 0, twopi);
        G4LogicalVolume* lens_logic = new G4LogicalVolume(lens_solid, Steel, "LENS"); // Set as steel for now so we can change reflectivity
        
        // Set optical properties
        G4OpticalSurface* lens_opsur = new G4OpticalSurface("LENS_OPSURF", unified, polished, dielectric_metal);
        lens_opsur->SetMaterialPropertiesTable(opticalprops::Absorber()); 

        // Sensitive detector
        SensorSD* lens_sd = new SensorSD("/LENS/LENS_OPSURF");
        lens_sd->SetDetectorVolumeDepth(2);
        lens_sd->SetTimeBinning(10000*ns);
        G4SDManager::GetSDMpointer()->AddNewDetector(lens_sd);
        lens_logic->SetSensitiveDetector(lens_sd);

        

        // Place the Volumes

        // LAB
        auto labPhysical= new G4PVPlacement(0,G4ThreeVector(),lab_logic_volume,lab_logic_volume->GetName(),0, false, 0, false);

        // Vessel
        G4VPhysicalVolume * vessel_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), vessel_logic, vessel_solid->GetName(), lab_logic_volume, false, 0, true);

        // Flanges
        G4VPhysicalVolume * flange1_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., (vessel_length+flange_thick)/2.0), flange_logic, "FLANGE1", lab_logic_volume, false, 0, true);
        G4VPhysicalVolume * flange2_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., -1*(vessel_length+flange_thick)/2.0), flange_logic, "FLANGE2", lab_logic_volume, false, 0, true);

        // Xenon Gas
        G4VPhysicalVolume * gas_phys= new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_logic, gas_solid->GetName(), lab_logic_volume, false, 0, true);

        // Field Cage -- Placed in the gas
        G4VPhysicalVolume * field_cage_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), field_cage_logic, field_cage_solid->GetName(), gas_logic, false, 0, true);

        // Lens -- Also placed in the gas
        // 2"(5.08 cm) shifted to right in x, 3.46" (8.7884 cm) down in y
        G4double lens_pos = 12.5*cm; // Position of the lens relative to cathode-side end-cap
        G4VPhysicalVolume * lens_phys = new G4PVPlacement(0, G4ThreeVector(-5.08*cm, -8.7884*cm, vessel_length/2.0-lens_thick/2.0-lens_pos), lens_logic, lens_solid->GetName(), gas_logic, false, 0, true);

        // VERTEX GENERATORS /////////////////////////////////////
        G4double anode_pos = 19*cm; // Distance from Anode side vessel wall -- Measured insitu
        anode_gen_          = new CylinderPointSampler(0, 10*mm, 0.5*um, 0., twopi, nullptr, G4ThreeVector (0,0,-vessel_length/2. + anode_pos));

        // Xenon Gas in Active Area and Non-Active Area
        IonizationSD* gasSD = new IonizationSD("/KingCRAB/GAS");
        gas_logic->SetSensitiveDetector(gasSD);
        G4SDManager::GetSDMpointer()->AddNewDetector(gasSD);

        /// Limit the step size in this volume for better tracking precision
        gas_logic->SetUserLimits(new G4UserLimits(max_step_size_));


        this->SetLogicalVolume(lab_logic_volume);

        AssignVisuals();

    }

    void KingCRAB::AssignVisuals() {
        // Chamber
        G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

        // Lab
        G4LogicalVolume* Lab = lvStore->GetVolume("LAB");
        Lab->SetVisAttributes(G4VisAttributes::GetInvisible());


        // Vessel
        G4VisAttributes *VesselVa=new G4VisAttributes(nexus::TitaniumGreyAlpha());
        VesselVa->SetForceSolid(true);
        
        G4LogicalVolume* Vessel = lvStore->GetVolume("VESSEL");
        Vessel->SetVisAttributes(VesselVa);
        // Chamber->SetVisAttributes(G4VisAttributes::GetInvisible());

        G4LogicalVolume* Flange1 = lvStore->GetVolume("FLANGE");
        Flange1->SetVisAttributes(VesselVa);

        G4LogicalVolume* Flange2 = lvStore->GetVolume("FLANGE");
        Flange2->SetVisAttributes(VesselVa);

        // Field Cage
        G4VisAttributes *FieldCageVa=new G4VisAttributes(nexus::WhiteAlpha());
        FieldCageVa->SetForceSolid(true);
        G4LogicalVolume* FieldCage = lvStore->GetVolume("FIELD_CAGE");
        FieldCage->SetVisAttributes(FieldCageVa);

        // Lens
        G4VisAttributes *LensVa=new G4VisAttributes(nexus::Red());
        LensVa->SetForceSolid(true);
        G4LogicalVolume* Lens = lvStore->GetVolume("LENS");
        Lens->SetVisAttributes(LensVa);


        // GAS
        G4LogicalVolume* Gas = lvStore->GetVolume("GAS");
        G4VisAttributes *GasVa=new G4VisAttributes(nexus::BlueAlpha());
        GasVa->SetForceCloud(true);
        Gas->SetVisAttributes(GasVa);

    }
    void KingCRAB::PrintParam() {
    }
    G4ThreeVector KingCRAB::GenerateVertex(const G4String& region) const
    {
            G4ThreeVector pos;
            if (region == "AD_HOC") {
                pos = vtx_;
            }
            else if((region == "ANODE")){
                pos = anode_gen_->GenerateVertex(VOLUME);
            }
            else {
                G4Exception("[KingCRAB]", "GenerateVertex()", JustWarning,
                            "Unknown vertex generation region. setting default region as 0,0,0");
                pos=vtx_;
            }
        return pos;
    }
}
