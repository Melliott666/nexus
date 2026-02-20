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
#include "HexagonMeshTools.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <G4UserLimits.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include "CylinderPointSampler.h"
#include <G4Sphere.hh>
#include <G4IntersectionSolid.hh>
#include <G4UnionSolid.hh>




namespace nexus{
    using namespace CLHEP;
    REGISTER_CLASS(KingCRAB, GeometryBase)


             KingCRAB::KingCRAB():
             GeometryBase(),
             msg_(nullptr),
             Lab_size(550. * m),
             anode_gen_(nullptr),
             max_step_size_(0.1*mm), 
             gas_pressure_(1. * bar),
             gastype_("xenon"),
             specific_vertex_(0,0,0)

    {
        msg_ = new G4GenericMessenger(this, "/Geometry/KingCRAB/","Control commands of geometry of KingCRAB TPC");
        G4GenericMessenger::Command&  Pressure_cmd =msg_->DeclarePropertyWithUnit("gas_pressure","bar",gas_pressure_,"Pressure of Gas");
        Pressure_cmd.SetParameterName("XenonPressure", false);

        G4GenericMessenger::Command&  step_size_cmd =msg_->DeclarePropertyWithUnit("max_step_size","mm",max_step_size_,"The maximum step size");
        step_size_cmd .SetParameterName("max_step_size", false);

        msg_->DeclareProperty("gastype", gastype_, "The GAS to use in the detector");



        msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_, "Set generation vertex.");

    }

    KingCRAB::~KingCRAB()
    {

        delete msg_;
        delete anode_gen_;


    }


    void KingCRAB::Construct(){

        // --------------------------
        // Materials
        // --------------------------
        
        // Setup the gas type
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
        
        // Other Materials
        G4Material *Steel=materials::Steel();
        Steel->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());

        G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu");
        G4Material* PTFE = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE");

        // --------------------------
        // Optical Surfaces
        // --------------------------
        // Add optical surface of gas to steel surfaces
        G4OpticalSurface* gas_steel_opsur = new G4OpticalSurface("GAS_STEEL_OPSURF", unified, ground, dielectric_metal);
        gas_steel_opsur->SetSigmaAlpha(0.0);
        // gas_steel_opsur->SetMaterialPropertiesTable(opticalprops::NoRef()); // No Reflections
        gas_steel_opsur->SetMaterialPropertiesTable(opticalprops::Steel()); // With Reflections

        // Add optical surface of gas to copper surfaces
        G4OpticalSurface* gas_copper_opsur = new G4OpticalSurface("GAS_COPPER_OPSURF", unified, ground, dielectric_metal);
        gas_copper_opsur->SetSigmaAlpha(0.0);
        gas_copper_opsur->SetMaterialPropertiesTable(opticalprops::Copper());

        // Optical surface gas to poly
        G4OpticalSurface* gas_poly_opsurf = new G4OpticalSurface("GAS_POLY_OPSURF", unified, ground, dielectric_metal, .01);
        gas_poly_opsurf->SetMaterialPropertiesTable(opticalprops::PTFE());

        // --------------------------
        // Constructing Lab Space
        // --------------------------
        G4String lab_name="LAB";
        G4Box * lab_solid_volume = new G4Box(lab_name,Lab_size/2.,Lab_size/2.,Lab_size/2.);
        G4LogicalVolume * lab_logic_volume= new G4LogicalVolume(lab_solid_volume,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),lab_name) ;
        auto labPhysical= new G4PVPlacement(0,G4ThreeVector(), lab_logic_volume, lab_logic_volume->GetName(), 0, false, 0, false);
    
        // --------------------------
        // Vessel cylinder
        // --------------------------
        // Inner diameter = 19.25" = 48.9 cm, Vessel Length = 56.12" = 142.5 cm, vessel thick
        G4double vessel_diam   = 48.895*cm;
        G4double vessel_length = 142.5448*cm;
        G4double vessel_thickn = 2.54*cm;
        G4double anode_pos = 19*cm; // Distance from Anode side vessel wall -- Measured insitu
        G4double z_shift = vessel_length/2.0 - anode_pos; // Shifts so the center of the EL gap is z=0
        
        G4Tubs* vessel_solid = new G4Tubs("VESSEL", vessel_diam/2.0, vessel_diam/2.0+vessel_thickn, vessel_length/2.0, 0, twopi);
        G4LogicalVolume* vessel_logic = new G4LogicalVolume(vessel_solid, Steel, "VESSEL");
        G4VPhysicalVolume * vessel_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., z_shift), vessel_logic, vessel_solid->GetName(), lab_logic_volume, false, 0, true);
        new G4LogicalSkinSurface("GAS_VESSEL_OPSURF", vessel_logic, gas_steel_opsur); // Reflections gas-steel vessel

        // --------------------------
        // Flanges
        // --------------------------
        // Diameter = 28" = 71.12 cm, thickness 1"=2.54cm
        G4double flange_diam  = 71.12*cm;
        G4double flange_thick = 2.54*cm;
        G4Tubs*  flange_solid = new G4Tubs("FLANGE", 0, flange_diam/2.0, flange_thick/2.0, 0, twopi);
        G4LogicalVolume* flange_logic = new G4LogicalVolume(flange_solid, Steel, "FLANGE");
        G4VPhysicalVolume * flange1_phys = new G4PVPlacement(0, G4ThreeVector(0., 0.,    (vessel_length+flange_thick)/2.0 + z_shift), flange_logic, "FLANGE1", lab_logic_volume, false, 0, true);
        G4VPhysicalVolume * flange2_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., -1*(vessel_length+flange_thick)/2.0 + z_shift), flange_logic, "FLANGE2", lab_logic_volume, false, 0, true);
        new G4LogicalSkinSurface("GAS_FLANGE_OPSURF", flange_logic, gas_steel_opsur); // Reflections gas-steel flange

        // --------------------------
        // Gas Volume
        // --------------------------
        G4Tubs* gas_solid = new G4Tubs("GAS", 0, vessel_diam/2., vessel_length/2., 0, twopi);
        G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, GAS, "GAS");
        G4VPhysicalVolume * gas_phys= new G4PVPlacement(0, G4ThreeVector(0.,0., z_shift), gas_logic, gas_solid->GetName(), lab_logic_volume, false, 0, true);

        // Xenon Gas sensitive volumne
        IonizationSD* gasSD = new IonizationSD("/KingCRAB/GAS");
        gas_logic->SetSensitiveDetector(gasSD);
        G4SDManager::GetSDMpointer()->AddNewDetector(gasSD);
        gas_logic->SetUserLimits(new G4UserLimits(max_step_size_)); // Limit the step size in this volume for better tracking precision
        this->SetLogicalVolume(lab_logic_volume);

        // --------------------------
        // Anode and EL rings.
        // --------------------------
        // Dim: thickness = 14mm, 382 mm inner diam, 437 mm outer diam
        G4double EL_ring_ID    = 382*mm;
        G4double EL_ring_OD    = 437*mm;
        G4double EL_ring_thick = 14*mm;
        G4double mesh_thick    = 130*um;
        G4double el_gap        = 7*mm;
        
        G4Tubs*          EL_solid    = new G4Tubs("RING_SOLID", EL_ring_ID/2., EL_ring_OD/2., EL_ring_thick/2. - mesh_thick/2., 0, twopi);
        G4LogicalVolume* anode_logic = new G4LogicalVolume(EL_solid, Steel, "ANODE_RING");
        G4LogicalVolume* EL_logic    = new G4LogicalVolume(EL_solid, Steel, "EL_RING");

        // Anode and EL Ring -- placement relative to gas volume centre
        // Shift in the -/+z direction by half-mesh thickness and reduce thickness
        // by the grid thickness. The grid thickness makes up the remaining ring thickness
        G4VPhysicalVolume * anode_ring = new G4PVPlacement(0, G4ThreeVector(0., 0., -el_gap/2.0 - mesh_thick/2. - EL_ring_thick/2.0 -z_shift), anode_logic, "ANODE_RING", gas_logic, false, 0, true);
        G4VPhysicalVolume * el_ring    = new G4PVPlacement(0, G4ThreeVector(0., 0., +el_gap/2.0 + mesh_thick/2. + EL_ring_thick/2.0 -z_shift), EL_logic,    "EL_RING",    gas_logic, false, 0, true);
        
        new G4LogicalSkinSurface("GAS_ANODE_OPSURF", anode_logic, gas_steel_opsur); // Reflections gas-steel anode/EL gate
        new G4LogicalSkinSurface("GAS_EL_GATE_OPSURF", EL_logic, gas_steel_opsur);
        
        // --------------------------
        // EL Mesh and Grids
        // --------------------------
        G4double EL_mesh_diam = 2.5*mm;

        // Dist from centre of hex to hex vertex, excluding the land width (circumradius)
        G4double hex_circumradius = EL_mesh_diam/std::sqrt(3)*mm;

        // Total number of hexagons that would fit side-by-side along the diameter
        G4int n_hex = (G4int) ((EL_ring_ID/2.0) / hex_circumradius);

        // Define the disk to punch hexagon holes through for the mesh
        G4Tubs* grid_solid = new G4Tubs("EL_GRID", 0., EL_ring_OD/2.0 , mesh_thick/2., 0., twopi);
        G4LogicalVolume* EL_grid_logic = new G4LogicalVolume(grid_solid, Steel, "EL_GRID");

        // Define a hexagonal prism
        G4ExtrudedSolid* hex_prism = CreateHexagon(mesh_thick/2.0, hex_circumradius);
        G4LogicalVolume* EL_hex_logic  = new G4LogicalVolume(hex_prism, GAS, "MESH_HEX_GAS");

        // Place GXe hexagons in the disk to make the mesh
        PlaceHexagons(n_hex, EL_mesh_diam, mesh_thick, EL_grid_logic, EL_hex_logic, EL_ring_ID);

        new G4LogicalSkinSurface("GAS_EL_MESH_OPSURF", EL_grid_logic, gas_steel_opsur); // Reflections gas-steel mesh

        // Create a rotation vector to change the orientation of the EL mesh
        CLHEP::HepRotationZ Roty(15*deg);
        G4RotationMatrix* pRot = new G4RotationMatrix();
        pRot->set(Roty);

        G4VPhysicalVolume * Anode_mesh = new G4PVPlacement(0,    G4ThreeVector(0., 0., -el_gap/2. - mesh_thick/2. - z_shift), EL_grid_logic, "EL_MESH_ANODE", gas_logic, false, 0, false);
        G4VPhysicalVolume * EL_mesh    = new G4PVPlacement(pRot, G4ThreeVector(0., 0., +el_gap/2. + mesh_thick/2. - z_shift), EL_grid_logic, "EL_MESH_GATE",  gas_logic, false, 1, false);

        // --------------------------
        // Poly Wrap
        // This is the plastic that wraps around the field cage
        // There is an error here. I think the diameter is wrongly extracted.
        // --------------------------
        G4double poly_inner_diam = 482.6*mm;
        G4double poly_length = 960*mm; // Use same as vessel length shifted by anode pos for now
        G4double poly_thickn = 6.35*mm; 
        G4Tubs*  poly_solid = new G4Tubs("POLY_WRAP", poly_inner_diam/2.0, poly_inner_diam/2.0+poly_thickn, poly_length/2.0, 0, twopi);
        G4LogicalVolume* poly_logic = new G4LogicalVolume(poly_solid, Steel, "POLY_WRAP"); // Set as steel for now so we can change reflectivity
        G4VPhysicalVolume * poly_phys = new G4PVPlacement(0, G4ThreeVector(0., 0.,  poly_length/2.0 - z_shift - el_gap/2. - mesh_thick/2. - EL_ring_thick), poly_logic, poly_solid->GetName(), gas_logic, false, 0, true);
        new G4LogicalSkinSurface("GAS_POLY_OPSURF", poly_logic, gas_poly_opsurf); // Reflections gas-poly

        // --------------------------
        // Field Cage Rings
        // --------------------------
        G4double field_ring_OD = 430*mm; // outer diameter
        G4double field_ring_thickn = 5*mm;
        G4Tubs*  field_ring_solid = new G4Tubs("FIELD_RING", field_ring_OD/2.0 - field_ring_thickn, field_ring_OD/2.0, field_ring_thickn/2.0, 0, twopi);
        G4LogicalVolume* field_ring_logic = new G4LogicalVolume(field_ring_solid, Steel, "FIELD_RING");
        
        // Positions FR closest to EL gate edge. 46mm is distance from EL gate edge to first FR centre
        G4double field_ring_origin = 46*mm + EL_ring_thick + mesh_thick/2. + el_gap/2. - z_shift; 

        G4double active_ring_sep = 24*mm; // separation between rings in active volume
        G4double posz;
        G4double n_active_FR = 31;
        for (G4int i=0; i<n_active_FR; i++) {
            posz = field_ring_origin + i*active_ring_sep;
            new G4PVPlacement(0, G4ThreeVector(0., 0., posz), field_ring_logic, field_ring_logic->GetName(), gas_logic, false, i, true);
        }

        // Add buffer rings
        G4double active_buffer_FR_dist = 72*mm; // Distance between last active vol field ring and buffer ring
        G4double buffer_ring_sep = 48*mm; // separation between rings in buffer volume

        for (G4int i=0; i<4; i++) {
            posz = field_ring_origin + 30*active_ring_sep + active_buffer_FR_dist + i*buffer_ring_sep;;
            new G4PVPlacement(0, G4ThreeVector(0., 0., posz), field_ring_logic, "FIELD_RING", gas_logic, false, i+31, false);
        }

        new G4LogicalSkinSurface("GAS_FIELDCAGE_OPSURF", field_ring_logic, gas_copper_opsur); // Reflections gas-copper field ring

        // --------------------------
        // Cathode Ring
        // --------------------------
        // Dim: thickness = 14mm, 321 mm inner diam, 373 mm outer diam
        G4double cathode_ring_ID    = 321*mm;
        G4double cathode_ring_OD    = 373*mm;
        G4double cathode_ring_thick = 14*mm;

        G4Tubs*          cathode_solid = new G4Tubs("CATHODE_RING_SOLID", cathode_ring_ID/2., cathode_ring_OD/2., cathode_ring_thick/2. - mesh_thick/2., 0, twopi);
        G4LogicalVolume* cathode_logic = new G4LogicalVolume(cathode_solid, Steel, "CATHODE_RING");

        // --------------------------
        // Cathode placement (ANCHOR TO END CAP WALL)
        // --------------------------
        // Cathode-side endcap wall (GAS frame)
        G4double z_cathode_wall = +vessel_length/2.0;

        // Measured from endcap wall to the OUTER mesh face (closest to wall)
        G4double cathode_inset_face = 39*cm;

        // Mesh center is half-thickness further into the detector
        G4double z_cathode_mesh = z_cathode_wall - (cathode_inset_face + mesh_thick/2.0);

        // Place ring so its +z face is flush with the mesh plane
        G4double z_cathode_ring = z_cathode_mesh - (mesh_thick/2.0 + cathode_ring_thick/2.0);

        // Cathode Ring -- placement relative to gas volume centre
        G4VPhysicalVolume* cathode_ring = new G4PVPlacement(0,G4ThreeVector(0., 0., z_cathode_ring), cathode_logic, "CATHODE_RING", gas_logic, false, 0, true);

        // --------------------------
        // Cathode Mesh and Grids
        // --------------------------

        // Total number of hexagons that would fit side-by-side along the diameter
        G4int n_hex_cathode = (G4int)((cathode_ring_ID/2.0) / hex_circumradius);

        // Define the disk to punch hexagon holes through for the mesh
        G4Tubs* cathode_grid_solid = new G4Tubs("CATHODE_GRID", 0., cathode_ring_OD/2.0, mesh_thick/2., 0, twopi);
        G4LogicalVolume* cathode_grid_logic = new G4LogicalVolume(cathode_grid_solid, Steel, "CATHODE_GRID");

        // Place GXe hexagons in the disk to make the mesh
        PlaceHexagons(n_hex_cathode, EL_mesh_diam, mesh_thick,
                      cathode_grid_logic, EL_hex_logic, cathode_ring_ID);

        new G4LogicalSkinSurface("GAS_CATHODE_MESH_OPSURF", cathode_grid_logic, gas_steel_opsur); // Reflections gas-steel mesh

        // Cathode Mesh -- placement relative to gas volume centre
        G4VPhysicalVolume* cathode_mesh = new G4PVPlacement(0, G4ThreeVector(0., 0., z_cathode_mesh), cathode_grid_logic, "CATHODE_MESH", gas_logic, false, 0, false);

        // --------------------------
        // Plastic Staves around Field Cage Rings
        // --------------------------
        // Stave dimensions 
        G4double Stave_ID = 6*mm;   // inner diameter (0 -> solid rod)
        G4double Stave_OD = 44.7*mm;  // outer diameter
        G4double Stave_thick = (Stave_OD - Stave_ID)/2.0;

        // How many staves around the cage 
        G4int n_Staves = 11;

        // Place Staves just outside the field rings but inside the poly wrap
        G4double Stave_rpos = field_ring_OD/2.0 + Stave_OD/2.0 + 5*mm;

        // Span tStaves over the full field-ring stack (active + buffer) with some margin
        G4double z_first_FR = field_ring_origin;
        G4double z_last_FR  = field_ring_origin + 30*active_ring_sep + active_buffer_FR_dist + 3*buffer_ring_sep;
        G4double Stave_length = (z_last_FR - z_first_FR); 
        G4double Stave_zpos   = (z_first_FR + z_last_FR)/2.0;

        G4Tubs*  Stave_solid = new G4Tubs("STAVES", Stave_ID/2.0, Stave_OD/2.0, Stave_length/2.0, 0, twopi);
        G4LogicalVolume* Stave_logic = new G4LogicalVolume(Stave_solid, PTFE, "STAVES");

        // Place Staves evenly spaced in phi
        for (G4int i=0; i<n_Staves; i++) {G4double phi = i*twopi/n_Staves; G4double xt = Stave_rpos*std::cos(phi); G4double yt = Stave_rpos*std::sin(phi); new G4PVPlacement(0, G4ThreeVector(xt, yt, Stave_zpos), Stave_logic, Stave_solid->GetName(), gas_logic, false, i, true);}

        // --------------------------
        // CaF2 Plano-Convex Lens (Thorlabs LA5210)
        // --------------------------
        G4Material* caf2_mat = materials::CaF2();
        caf2_mat->SetMaterialPropertiesTable(opticalprops::CaF2());

        // LA5210 (50.8 mm dia, f=100 mm) mechanical shape
        G4double Lens_D  = 50.8*mm;
        G4double Lens_R  = 43.4*mm;   // convex ROC
        G4double Lens_tc = 11.2*mm;   // center thickness
        G4double Lens_a  = Lens_D/2.0;

        // Build lens = (blank cylinder) - (outside-of-sphere within aperture)
        G4double Lens_zc_shift = Lens_tc/2.0 - Lens_R;

        G4Tubs*   Lens_blank  = new G4Tubs("CAF2_LENS_BLANK",  0., Lens_a,   Lens_tc/2.0, 0., twopi);
        G4Sphere* Lens_sphere = new G4Sphere("CAF2_LENS_SPHERE", 0., Lens_R, 0., twopi, 0., pi);
        G4Tubs*   Lens_big    = new G4Tubs("CAF2_LENS_BIG",    0., Lens_a,   Lens_R,      0., twopi);

        G4SubtractionSolid* Lens_outside = new G4SubtractionSolid("CAF2_LENS_OUTSIDE", Lens_big, Lens_sphere, 0, G4ThreeVector(0.,0.,Lens_zc_shift));

        G4SubtractionSolid* Lens_solid   = new G4SubtractionSolid("CAF2_LENS_SOLID", Lens_blank, Lens_outside, 0, G4ThreeVector(0.,0.,0.));

        G4LogicalVolume* Lens_logic = new G4LogicalVolume(Lens_solid, caf2_mat, "CAF2_LENS");

        // --------------------------
        // Lens placement (anchored to endcap) + ORIENTATION
        // --------------------------

        // Endcap wall position in GAS frame
        G4double z_endcap_wall = vessel_length/2.0;

        // 8 cm from endcap wall to LENS FACE (outer face nearest wall)
        G4double lens_offset = 8.0*cm;

        // Lens center position
        G4double Lens_zpos = z_endcap_wall - (lens_offset + Lens_tc/2.0);

        // Flip lens so convex side faces toward cathode (i.e. toward -z)
        G4RotationMatrix* Lens_rot = new G4RotationMatrix();
        Lens_rot->rotateY(180.0*deg);

        new G4PVPlacement(Lens_rot, G4ThreeVector(0., 0., Lens_zpos), Lens_logic, "CAF2_LENS", gas_logic, false, 0, true);

        // --------------------------
        // 45-deg Circular Mirror (Ø50.8 mm, T=10 mm)
        // --------------------------
        G4double Mirror_D = 50.8*mm;
        G4double Mirror_T = 10.0*mm;

        // 5 cm from lens toward endcap (same relative spacing as before)
        G4double Mirror_zpos = Lens_zpos + (5.0*cm);

        G4Tubs* Mirror_solid = new G4Tubs("MIRROR", 0., Mirror_D/2.0, Mirror_T/2.0, 0., twopi);
        G4LogicalVolume* Mirror_logic = new G4LogicalVolume(Mirror_solid, Steel, "MIRROR");

        // Yaw angle that rotates the old y-axis toward +x by theta.
        G4double theta = -29.9963208064*deg; // arctan(3.250/5.640)

        // Mirror #1 rotation
        G4RotationMatrix* Mirror_rot = new G4RotationMatrix();
        Mirror_rot->rotateZ(theta);
        Mirror_rot->rotateX(135.0*deg);

        new G4PVPlacement(Mirror_rot, G4ThreeVector(0., 0., Mirror_zpos),
                  Mirror_logic, "MIRROR", gas_logic, false, 0, true);

        // ---- Mirror #1 optical surface: PERFECT MIRROR
        auto* mirror1_opsur = new G4OpticalSurface("MIRROR1_OPSURF",
                                           unified, polished, dielectric_metal);
        mirror1_opsur->SetMaterialPropertiesTable(opticalprops::PerfectMirror());

        // Apply as skin surface (affects all faces of the solid)
        new G4LogicalSkinSurface("GAS_MIRROR1_SKIN", Mirror_logic, mirror1_opsur);


        // --------------------------
        // 45-deg Circular Mirror #2 (Ø50.8 mm, T=10 mm)
        // --------------------------
        G4double Mirror2_xpos = -8.255*cm;
        G4double Mirror2_ypos = 14.3002*cm;
        G4double Mirror2_zpos = Mirror_zpos;

        G4Tubs* Mirror2_solid = new G4Tubs("MIRROR2", 0., Mirror_D/2.0, Mirror_T/2.0, 0., twopi);
        G4LogicalVolume* Mirror2_logic = new G4LogicalVolume(Mirror2_solid, Steel, "MIRROR2");

        // Mirror #2 rotation
        G4RotationMatrix* Mirror2_rot = new G4RotationMatrix();
        Mirror2_rot->rotateZ(theta);
        Mirror2_rot->rotateX(-45.0*deg);

        new G4PVPlacement(Mirror2_rot, G4ThreeVector(Mirror2_xpos, Mirror2_ypos, Mirror2_zpos),
                  Mirror2_logic, "MIRROR2", gas_logic, false, 0, true);

        // ---- Mirror #2 optical surface: PERFECT MIRROR
        auto* mirror2_opsur = new G4OpticalSurface("MIRROR2_OPSURF",
                                           unified, polished, dielectric_metal);
        mirror2_opsur->SetMaterialPropertiesTable(opticalprops::PerfectMirror());

        new G4LogicalSkinSurface("GAS_MIRROR2_SKIN", Mirror2_logic, mirror2_opsur);
                
        // --------------------------
        // Window (Detector)
        // --------------------------
        G4double lens_diam = 5.08*cm; // 2 inch
        G4double lens_thick = 2*mm; // Choose arbitary for now
        G4Tubs* lens_solid = new G4Tubs("LENS", 0, lens_diam/2.0, lens_thick/2.0, 0, twopi);
        G4LogicalVolume* lens_logic = new G4LogicalVolume(lens_solid, Steel, "LENS"); // Set as steel for now, need to add CaF2 material properties
        
        // Set optical properties of lens
        G4OpticalSurface* lens_opsur = new G4OpticalSurface("LENS", unified, polished, dielectric_metal);
        lens_opsur->SetMaterialPropertiesTable(opticalprops::Absorber()); // 100% Quantum Efficiency
        new G4LogicalSkinSurface("LENS", lens_logic, lens_opsur);

        // Placement relative to gas volume centre
        G4VPhysicalVolume * lens_phys = new G4PVPlacement(0, G4ThreeVector(Mirror2_xpos, Mirror2_ypos, vessel_length/2 - lens_thick -2*mm), lens_logic, lens_solid->GetName(), gas_logic, false, 0, true);

        // --------------------------
        // VERTEX GENERATORS 
        // --------------------------
        anode_gen_          = new CylinderPointSampler(0, 10*mm, 0.5*um, 0., twopi, nullptr, G4ThreeVector (0,0,0)); // Generate in center of EL gap
        cathode_gen_        = new CylinderPointSampler(0.0, cathode_ring_ID/2.0, 0.5*um, 0.0, twopi, nullptr, G4ThreeVector(0., 0., z_shift + z_cathode_mesh)); // Generate in center of Cathode
        active_volume_gen_  = new CylinderPointSampler(0.0, vessel_diam/2.0, vessel_length/2.0, 0.0, twopi, nullptr, G4ThreeVector(0., 0., z_shift)); // Generate uniformly through the volume


        // --------------------------
        // Visuals 
        // --------------------------
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

        G4LogicalVolume* EL_Gate = lvStore->GetVolume("EL_RING");
        EL_Gate->SetVisAttributes(VesselVa);
        
        G4LogicalVolume* Anode = lvStore->GetVolume("ANODE_RING");
        Anode->SetVisAttributes(VesselVa);

        // Poly Wrap
        G4VisAttributes *PolyWrapVa=new G4VisAttributes(nexus::WhiteAlpha());
        PolyWrapVa->SetForceSolid(true);
        G4LogicalVolume* PolyWrap = lvStore->GetVolume("POLY_WRAP");
        PolyWrap->SetVisAttributes(PolyWrapVa);

        // Lens
        G4VisAttributes *LensVa=new G4VisAttributes(nexus::Red());
        LensVa->SetForceSolid(true);
        G4LogicalVolume* Lens = lvStore->GetVolume("LENS");
        Lens->SetVisAttributes(LensVa);

        // Field Rings
        G4VisAttributes *FieldRingVa=new G4VisAttributes(nexus::CopperBrownAlpha());
        FieldRingVa->SetForceSolid(true);
        G4LogicalVolume* FieldRing = lvStore->GetVolume("FIELD_RING");
        FieldRing->SetVisAttributes(FieldRingVa);


        // GAS
        G4LogicalVolume* Gas = lvStore->GetVolume("GAS");
        G4VisAttributes *GasVa=new G4VisAttributes(nexus::BlueAlpha());
        GasVa->SetForceCloud(true);
        Gas->SetVisAttributes(GasVa);

        // Staves
        G4VisAttributes *StaveVa =new G4VisAttributes(nexus::WhiteAlpha());
        StaveVa-> SetForceSolid(true);
        G4LogicalVolume* Tube = lvStore->GetVolume("STAVES");
        Tube->SetVisAttributes(StaveVa);

        // CaF2 Lens
        G4VisAttributes* caf2Va = new G4VisAttributes(nexus::Blue());
        caf2Va->SetForceSolid(true);
        G4LogicalVolume* caf2LV = lvStore->GetVolume("CAF2_LENS");
        caf2LV->SetVisAttributes(caf2Va);

        // Mirrors 
        G4VisAttributes* MirrorVa = new G4VisAttributes(nexus::Blue());
        MirrorVa->SetForceSolid(true);

        G4LogicalVolume* Mirror1 = lvStore->GetVolume("MIRROR");
        if (Mirror1) Mirror1->SetVisAttributes(MirrorVa);

        G4LogicalVolume* Mirror2 = lvStore->GetVolume("MIRROR2");
        if (Mirror2) Mirror2->SetVisAttributes(MirrorVa);

        // EL Gate (bright green)
        auto* ELVa = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1.0));
        ELVa->SetForceSolid(true);
        if (EL_Gate) EL_Gate->SetVisAttributes(ELVa);

        // Anode (bright green)
        auto* AnodeVa = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 1.0));
        AnodeVa->SetForceSolid(true);
        if (Anode) Anode->SetVisAttributes(AnodeVa);

        // --------------------------
        // Cathode (bright magenta)
        // --------------------------
        G4VisAttributes* CathodeVa =
            new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));  // bright magenta (RGB)

        CathodeVa->SetForceSolid(true);

        G4LogicalVolume* CathodeRingLV = lvStore->GetVolume("CATHODE_RING");
        if (CathodeRingLV) CathodeRingLV->SetVisAttributes(CathodeVa);

        G4LogicalVolume* CathodeMeshLV = lvStore->GetVolume("CATHODE_GRID");
        if (CathodeMeshLV) CathodeMeshLV->SetVisAttributes(CathodeVa);





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
            else {
                G4Exception("[KingCRAB]", "GenerateVertex()", JustWarning,
                            "Unknown vertex generation region. setting default region as 0,0,0");
                pos=specific_vertex_;
            }
        return pos;
    }
}
