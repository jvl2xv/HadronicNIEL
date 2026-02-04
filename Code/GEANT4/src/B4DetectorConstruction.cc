//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class 

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


#include "G4Region.hh"
#include "G4Proton.hh"
#include "G4UserLimits.hh"



// CADMESH //
#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4Tet.hh"
#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4EmCalculator.hh"

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction(int numMaterial, int numDimension)
 : G4VUserDetectorConstruction(),
   fnumMaterial(numMaterial),
   fnumDimension(numDimension),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();

G4GeometryManager::GetInstance()->OpenGeometry();
G4PhysicalVolumeStore::GetInstance()->Clean();
G4LogicalVolumeStore::GetInstance()->Clean();
G4SolidStore::GetInstance()->Clean();
// G4RunManager::GetRunManager()->DefineWorldVolume(nullptr); // Ensure world volume is redefined - this causes a seg fault

  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{


  // Geometry parameters

  G4NistManager* nistManager = G4NistManager::Instance();

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 

  // Get materials
  //G4Material* G4_water = nistManager->FindOrBuildMaterial("G4_WATER");
 
  G4Material* G4_vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");
  

G4Material* matH_1gcm3 = new G4Material("H_mat", 1.0*g/cm3, 1);
matH_1gcm3->AddElement(nistManager->FindOrBuildElement("H"), 1);

G4Material* matHe_1gcm3 = new G4Material("He_mat", 1.0*g/cm3, 1);
matHe_1gcm3->AddElement(nistManager->FindOrBuildElement("He"), 1);

G4Material* matLi_1gcm3 = new G4Material("Li_mat", 1.0*g/cm3, 1);
matLi_1gcm3->AddElement(nistManager->FindOrBuildElement("Li"), 1);

G4Material* matBe_1gcm3 = new G4Material("Be_mat", 1.0*g/cm3, 1);
matBe_1gcm3->AddElement(nistManager->FindOrBuildElement("Be"), 1);

G4Material* matB_1gcm3 = new G4Material("B_mat", 1.0*g/cm3, 1);
matB_1gcm3->AddElement(nistManager->FindOrBuildElement("B"), 1);

G4Material* matC_1gcm3 = new G4Material("C_mat", 1.0*g/cm3, 1);
matC_1gcm3->AddElement(nistManager->FindOrBuildElement("C"), 1);

G4Material* matN_1gcm3 = new G4Material("N_mat", 1.0*g/cm3, 1);
matN_1gcm3->AddElement(nistManager->FindOrBuildElement("N"), 1);

G4Material* matO_1gcm3 = new G4Material("O_mat", 1.0*g/cm3, 1);
matO_1gcm3->AddElement(nistManager->FindOrBuildElement("O"), 1);

G4Material* matF_1gcm3 = new G4Material("F_mat", 1.0*g/cm3, 1);
matF_1gcm3->AddElement(nistManager->FindOrBuildElement("F"), 1);

G4Material* matNe_1gcm3 = new G4Material("Ne_mat", 1.0*g/cm3, 1);
matNe_1gcm3->AddElement(nistManager->FindOrBuildElement("Ne"), 1);

G4Material* matNa_1gcm3 = new G4Material("Na_mat", 1.0*g/cm3, 1);
matNa_1gcm3->AddElement(nistManager->FindOrBuildElement("Na"), 1);

G4Material* matMg_1gcm3 = new G4Material("Mg_mat", 1.0*g/cm3, 1);
matMg_1gcm3->AddElement(nistManager->FindOrBuildElement("Mg"), 1);

G4Material* matAl_1gcm3 = new G4Material("Al_mat", 1.0*g/cm3, 1);
matAl_1gcm3->AddElement(nistManager->FindOrBuildElement("Al"), 1);

G4Material* matSi_1gcm3 = new G4Material("Si_mat", 1.0*g/cm3, 1);
matSi_1gcm3->AddElement(nistManager->FindOrBuildElement("Si"), 1);

G4Material* matP_1gcm3 = new G4Material("P_mat", 1.0*g/cm3, 1);
matP_1gcm3->AddElement(nistManager->FindOrBuildElement("P"), 1);

G4Material* matS_1gcm3 = new G4Material("S_mat", 1.0*g/cm3, 1);
matS_1gcm3->AddElement(nistManager->FindOrBuildElement("S"), 1);

G4Material* matCl_1gcm3 = new G4Material("Cl_mat", 1.0*g/cm3, 1);
matCl_1gcm3->AddElement(nistManager->FindOrBuildElement("Cl"), 1);

G4Material* matAr_1gcm3 = new G4Material("Ar_mat", 1.0*g/cm3, 1);
matAr_1gcm3->AddElement(nistManager->FindOrBuildElement("Ar"), 1);

G4Material* matK_1gcm3 = new G4Material("K_mat", 1.0*g/cm3, 1);
matK_1gcm3->AddElement(nistManager->FindOrBuildElement("K"), 1);

G4Material* matCa_1gcm3 = new G4Material("Ca_mat", 1.0*g/cm3, 1);
matCa_1gcm3->AddElement(nistManager->FindOrBuildElement("Ca"), 1);

G4Material* matSc_1gcm3 = new G4Material("Sc_mat", 1.0*g/cm3, 1);
matSc_1gcm3->AddElement(nistManager->FindOrBuildElement("Sc"), 1);

G4Material* matTi_1gcm3 = new G4Material("Ti_mat", 1.0*g/cm3, 1);
matTi_1gcm3->AddElement(nistManager->FindOrBuildElement("Ti"), 1);

G4Material* matV_1gcm3 = new G4Material("V_mat", 1.0*g/cm3, 1);
matV_1gcm3->AddElement(nistManager->FindOrBuildElement("V"), 1);

G4Material* matCr_1gcm3 = new G4Material("Cr_mat", 1.0*g/cm3, 1);
matCr_1gcm3->AddElement(nistManager->FindOrBuildElement("Cr"), 1);

G4Material* matMn_1gcm3 = new G4Material("Mn_mat", 1.0*g/cm3, 1);
matMn_1gcm3->AddElement(nistManager->FindOrBuildElement("Mn"), 1);

G4Material* matFe_1gcm3 = new G4Material("Fe_mat", 1.0*g/cm3, 1);
matFe_1gcm3->AddElement(nistManager->FindOrBuildElement("Fe"), 1);

G4Material* matCo_1gcm3 = new G4Material("Co_mat", 1.0*g/cm3, 1);
matCo_1gcm3->AddElement(nistManager->FindOrBuildElement("Co"), 1);

G4Material* matNi_1gcm3 = new G4Material("Ni_mat", 1.0*g/cm3, 1);
matNi_1gcm3->AddElement(nistManager->FindOrBuildElement("Ni"), 1);

G4Material* matCu_1gcm3 = new G4Material("Cu_mat", 1.0*g/cm3, 1);
matCu_1gcm3->AddElement(nistManager->FindOrBuildElement("Cu"), 1);

G4Material* matZn_1gcm3 = new G4Material("Zn_mat", 1.0*g/cm3, 1);
matZn_1gcm3->AddElement(nistManager->FindOrBuildElement("Zn"), 1);

G4Material* matGa_1gcm3 = new G4Material("Ga_mat", 1.0*g/cm3, 1);
matGa_1gcm3->AddElement(nistManager->FindOrBuildElement("Ga"), 1);

G4Material* matGe_1gcm3 = new G4Material("Ge_mat", 1.0*g/cm3, 1);
matGe_1gcm3->AddElement(nistManager->FindOrBuildElement("Ge"), 1);

G4Material* matAs_1gcm3 = new G4Material("As_mat", 1.0*g/cm3, 1);
matAs_1gcm3->AddElement(nistManager->FindOrBuildElement("As"), 1);

G4Material* matSe_1gcm3 = new G4Material("Se_mat", 1.0*g/cm3, 1);
matSe_1gcm3->AddElement(nistManager->FindOrBuildElement("Se"), 1);

G4Material* matBr_1gcm3 = new G4Material("Br_mat", 1.0*g/cm3, 1);
matBr_1gcm3->AddElement(nistManager->FindOrBuildElement("Br"), 1);

G4Material* matKr_1gcm3 = new G4Material("Kr_mat", 1.0*g/cm3, 1);
matKr_1gcm3->AddElement(nistManager->FindOrBuildElement("Kr"), 1);

G4Material* matRb_1gcm3 = new G4Material("Rb_mat", 1.0*g/cm3, 1);
matRb_1gcm3->AddElement(nistManager->FindOrBuildElement("Rb"), 1);

G4Material* matSr_1gcm3 = new G4Material("Sr_mat", 1.0*g/cm3, 1);
matSr_1gcm3->AddElement(nistManager->FindOrBuildElement("Sr"), 1);

G4Material* matY_1gcm3 = new G4Material("Y_mat", 1.0*g/cm3, 1);
matY_1gcm3->AddElement(nistManager->FindOrBuildElement("Y"), 1);

G4Material* matZr_1gcm3 = new G4Material("Zr_mat", 1.0*g/cm3, 1);
matZr_1gcm3->AddElement(nistManager->FindOrBuildElement("Zr"), 1);

G4Material* matNb_1gcm3 = new G4Material("Nb_mat", 1.0*g/cm3, 1);
matNb_1gcm3->AddElement(nistManager->FindOrBuildElement("Nb"), 1);

G4Material* matMo_1gcm3 = new G4Material("Mo_mat", 1.0*g/cm3, 1);
matMo_1gcm3->AddElement(nistManager->FindOrBuildElement("Mo"), 1);

G4Material* matRu_1gcm3 = new G4Material("Ru_mat", 1.0*g/cm3, 1);
matRu_1gcm3->AddElement(nistManager->FindOrBuildElement("Ru"), 1);

G4Material* matRh_1gcm3 = new G4Material("Rh_mat", 1.0*g/cm3, 1);
matRh_1gcm3->AddElement(nistManager->FindOrBuildElement("Rh"), 1);

G4Material* matPd_1gcm3 = new G4Material("Pd_mat", 1.0*g/cm3, 1);
matPd_1gcm3->AddElement(nistManager->FindOrBuildElement("Pd"), 1);

G4Material* matAg_1gcm3 = new G4Material("Ag_mat", 1.0*g/cm3, 1);
matAg_1gcm3->AddElement(nistManager->FindOrBuildElement("Ag"), 1);

G4Material* matCd_1gcm3 = new G4Material("Cd_mat", 1.0*g/cm3, 1);
matCd_1gcm3->AddElement(nistManager->FindOrBuildElement("Cd"), 1);

G4Material* matIn_1gcm3 = new G4Material("In_mat", 1.0*g/cm3, 1);
matIn_1gcm3->AddElement(nistManager->FindOrBuildElement("In"), 1);

G4Material* matSn_1gcm3 = new G4Material("Sn_mat", 1.0*g/cm3, 1);
matSn_1gcm3->AddElement(nistManager->FindOrBuildElement("Sn"), 1);

G4Material* matSb_1gcm3 = new G4Material("Sb_mat", 1.0*g/cm3, 1);
matSb_1gcm3->AddElement(nistManager->FindOrBuildElement("Sb"), 1);

G4Material* matTe_1gcm3 = new G4Material("Te_mat", 1.0*g/cm3, 1);
matTe_1gcm3->AddElement(nistManager->FindOrBuildElement("Te"), 1);

G4Material* matI_1gcm3 = new G4Material("I_mat", 1.0*g/cm3, 1);
matI_1gcm3->AddElement(nistManager->FindOrBuildElement("I"), 1);

G4Material* matXe_1gcm3 = new G4Material("Xe_mat", 1.0*g/cm3, 1);
matXe_1gcm3->AddElement(nistManager->FindOrBuildElement("Xe"), 1);

G4Material* matCs_1gcm3 = new G4Material("Cs_mat", 1.0*g/cm3, 1);
matCs_1gcm3->AddElement(nistManager->FindOrBuildElement("Cs"), 1);

G4Material* matBa_1gcm3 = new G4Material("Ba_mat", 1.0*g/cm3, 1);
matBa_1gcm3->AddElement(nistManager->FindOrBuildElement("Ba"), 1);

G4Material* matLa_1gcm3 = new G4Material("La_mat", 1.0*g/cm3, 1);
matLa_1gcm3->AddElement(nistManager->FindOrBuildElement("La"), 1);

G4Material* matCe_1gcm3 = new G4Material("Ce_mat", 1.0*g/cm3, 1);
matCe_1gcm3->AddElement(nistManager->FindOrBuildElement("Ce"), 1);

G4Material* matPr_1gcm3 = new G4Material("Pr_mat", 1.0*g/cm3, 1);
matPr_1gcm3->AddElement(nistManager->FindOrBuildElement("Pr"), 1);

G4Material* matNd_1gcm3 = new G4Material("Nd_mat", 1.0*g/cm3, 1);
matNd_1gcm3->AddElement(nistManager->FindOrBuildElement("Nd"), 1);

G4Material* matSm_1gcm3 = new G4Material("Sm_mat", 1.0*g/cm3, 1);
matSm_1gcm3->AddElement(nistManager->FindOrBuildElement("Sm"), 1);

G4Material* matEu_1gcm3 = new G4Material("Eu_mat", 1.0*g/cm3, 1);
matEu_1gcm3->AddElement(nistManager->FindOrBuildElement("Eu"), 1);

G4Material* matGd_1gcm3 = new G4Material("Gd_mat", 1.0*g/cm3, 1);
matGd_1gcm3->AddElement(nistManager->FindOrBuildElement("Gd"), 1);

G4Material* matTb_1gcm3 = new G4Material("Tb_mat", 1.0*g/cm3, 1);
matTb_1gcm3->AddElement(nistManager->FindOrBuildElement("Tb"), 1);

G4Material* matDy_1gcm3 = new G4Material("Dy_mat", 1.0*g/cm3, 1);
matDy_1gcm3->AddElement(nistManager->FindOrBuildElement("Dy"), 1);

G4Material* matHo_1gcm3 = new G4Material("Ho_mat", 1.0*g/cm3, 1);
matHo_1gcm3->AddElement(nistManager->FindOrBuildElement("Ho"), 1);

G4Material* matEr_1gcm3 = new G4Material("Er_mat", 1.0*g/cm3, 1);
matEr_1gcm3->AddElement(nistManager->FindOrBuildElement("Er"), 1);

G4Material* matTm_1gcm3 = new G4Material("Tm_mat", 1.0*g/cm3, 1);
matTm_1gcm3->AddElement(nistManager->FindOrBuildElement("Tm"), 1);

G4Material* matYb_1gcm3 = new G4Material("Yb_mat", 1.0*g/cm3, 1);
matYb_1gcm3->AddElement(nistManager->FindOrBuildElement("Yb"), 1);

G4Material* matLu_1gcm3 = new G4Material("Lu_mat", 1.0*g/cm3, 1);
matLu_1gcm3->AddElement(nistManager->FindOrBuildElement("Lu"), 1);

G4Material* matHf_1gcm3 = new G4Material("Hf_mat", 1.0*g/cm3, 1);
matHf_1gcm3->AddElement(nistManager->FindOrBuildElement("Hf"), 1);

G4Material* matTa_1gcm3 = new G4Material("Ta_mat", 1.0*g/cm3, 1);
matTa_1gcm3->AddElement(nistManager->FindOrBuildElement("Ta"), 1);

G4Material* matW_1gcm3 = new G4Material("W_mat", 1.0*g/cm3, 1);
matW_1gcm3->AddElement(nistManager->FindOrBuildElement("W"), 1);

G4Material* matRe_1gcm3 = new G4Material("Re_mat", 1.0*g/cm3, 1);
matRe_1gcm3->AddElement(nistManager->FindOrBuildElement("Re"), 1);

G4Material* matOs_1gcm3 = new G4Material("Os_mat", 1.0*g/cm3, 1);
matOs_1gcm3->AddElement(nistManager->FindOrBuildElement("Os"), 1);

G4Material* matIr_1gcm3 = new G4Material("Ir_mat", 1.0*g/cm3, 1);
matIr_1gcm3->AddElement(nistManager->FindOrBuildElement("Ir"), 1);

G4Material* matPt_1gcm3 = new G4Material("Pt_mat", 1.0*g/cm3, 1);
matPt_1gcm3->AddElement(nistManager->FindOrBuildElement("Pt"), 1);

G4Material* matAu_1gcm3 = new G4Material("Au_mat", 1.0*g/cm3, 1);
matAu_1gcm3->AddElement(nistManager->FindOrBuildElement("Au"), 1);

G4Material* matHg_1gcm3 = new G4Material("Hg_mat", 1.0*g/cm3, 1);
matHg_1gcm3->AddElement(nistManager->FindOrBuildElement("Hg"), 1);

G4Material* matTl_1gcm3 = new G4Material("Tl_mat", 1.0*g/cm3, 1);
matTl_1gcm3->AddElement(nistManager->FindOrBuildElement("Tl"), 1);

G4Material* matPb_1gcm3 = new G4Material("Pb_mat", 1.0*g/cm3, 1);
matPb_1gcm3->AddElement(nistManager->FindOrBuildElement("Pb"), 1);

G4Material* matBi_1gcm3 = new G4Material("Bi_mat", 1.0*g/cm3, 1);
matBi_1gcm3->AddElement(nistManager->FindOrBuildElement("Bi"), 1);


G4Material* materials[] = {matH_1gcm3, matHe_1gcm3, matLi_1gcm3, matBe_1gcm3, matB_1gcm3, matC_1gcm3, matN_1gcm3, matO_1gcm3, matF_1gcm3, matNe_1gcm3, matNa_1gcm3, matMg_1gcm3, matAl_1gcm3, matSi_1gcm3, matP_1gcm3, matS_1gcm3, matCl_1gcm3, matAr_1gcm3, matK_1gcm3, matCa_1gcm3, matSc_1gcm3, matTi_1gcm3, matV_1gcm3, matCr_1gcm3, matMn_1gcm3, matFe_1gcm3, matCo_1gcm3, matNi_1gcm3, matCu_1gcm3, matZn_1gcm3, matGa_1gcm3, matGe_1gcm3, matAs_1gcm3, matSe_1gcm3, matBr_1gcm3, matKr_1gcm3, matRb_1gcm3, matSr_1gcm3, matY_1gcm3, matZr_1gcm3, matNb_1gcm3, matMo_1gcm3, matRu_1gcm3, matRh_1gcm3, matPd_1gcm3, matAg_1gcm3, matCd_1gcm3, matIn_1gcm3, matSn_1gcm3, matSb_1gcm3, matTe_1gcm3, matI_1gcm3, matXe_1gcm3, matCs_1gcm3, matBa_1gcm3, matLa_1gcm3, matCe_1gcm3, matPr_1gcm3, matNd_1gcm3, matSm_1gcm3, matEu_1gcm3, matGd_1gcm3, matTb_1gcm3, matDy_1gcm3, matHo_1gcm3, matEr_1gcm3, matTm_1gcm3, matYb_1gcm3, matLu_1gcm3, matHf_1gcm3, matTa_1gcm3, matW_1gcm3, matRe_1gcm3, matOs_1gcm3, matIr_1gcm3, matPt_1gcm3, matAu_1gcm3, matHg_1gcm3, matTl_1gcm3, matPb_1gcm3, matBi_1gcm3};

G4Material* materialToUse = materials[fnumMaterial];





  G4Box* world_solid
    = new G4Box("World", 1.0*cm, 1.0*cm, 1.0*cm); // its size
                         
  G4LogicalVolume* world_logical
    = new G4LogicalVolume(
                 world_solid,           // its solid
                 G4_vacuum,  // its material
                 "World");         // its name
                                   
  G4VPhysicalVolume* world_physical
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 world_logical,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


// only lose 5% of it's energy through the depth (stopping power in MeV/mm and initial energy is in MeV so this is just units of mm)

double thicknesses_mm[] = {0.6e-5, 0.6e-5, 0.6e-5, 0.6e-5, 1e-5, 1e-5, 1e-5, 1e-5, 2e-5, 2e-5, 0.5e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 5e-3, 5e-3, 5e-3 , 5e-3 , 5e-3, 5e-3}; 

G4double thickness_mm = thicknesses_mm[fnumDimension];
cout << "thickness [mm] = " << thickness_mm << endl;


G4double sample_length_cm = 0.5;

 G4Box* Sample_solid
    = new G4Box("Sample", (sample_length_cm/2.0)*cm, (sample_length_cm/2.0)*cm, (thickness_mm/2.0)*mm); // its size
                         
  G4LogicalVolume* Sample_logical
    = new G4LogicalVolume(
                 Sample_solid,           // its solid
                 materialToUse,  // its material
                 "Sample");         // its name
                                   
  G4VPhysicalVolume* Sample_physical
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0.0*cm, 0.0*cm, 0.0*cm),  // at (0,0,0)
                 Sample_logical,          // its logical volume                         
                 "Sample",          // its name
                 world_logical,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  //G4double max_step = 0.01*nm;
  //Sample_logical->SetUserLimits(new G4UserLimits(max_step));
    





  // Create Target G4Region and add logical volume

  /* TURNS OUT THIS IS THE SAME AS NORMAL CUTS
  
   G4Region* sample_region = new G4Region("Sample_Region");
  
  G4ProductionCuts* cuts = new G4ProductionCuts();
  
  G4double defCut = 0.1*nm;
  cuts->SetProductionCut(defCut,"gamma");
  cuts->SetProductionCut(defCut,"e-");
  cuts->SetProductionCut(defCut,"e+");
  cuts->SetProductionCut(defCut,"proton");
  
  sample_region->SetProductionCuts(cuts);
  Sample_logical->SetRegion(sample_region);
  sample_region->AddRootLogicalVolume(Sample_logical, false); 
  
*/ 



  G4Colour colorRed(G4Colour::Red());
  G4Colour colorGreen(G4Colour::Green());
  G4Colour colorBlue(G4Colour::Blue());
  G4Colour colorYellow(G4Colour::Yellow());


  G4VisAttributes* redAtt= new G4VisAttributes(colorRed);
  G4VisAttributes* blueAtt= new G4VisAttributes(colorBlue);
  G4VisAttributes* yellowAtt= new G4VisAttributes(colorYellow);








  return world_physical;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}


 void B4DetectorConstruction::SetParameters(G4int p) 
    {
        param_num = p;
        G4RunManager::GetRunManager()->ReinitializeGeometry();
    }

















//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
