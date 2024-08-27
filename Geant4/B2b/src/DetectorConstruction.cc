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
//
/// \file B2/B2b/src/DetectorConstruction.cc
/// \brief Implementation of the B2b::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ChamberParameterisation.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


using namespace B2;

namespace B2b
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
{
  fMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fStepLimit;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Material definition

  G4NistManager* nistManager = G4NistManager::Instance();

  G4int ncomponents, natoms;

  G4double A;  // atomic mass
  G4double Z;  // atomic number
  G4double d;  // density

  /*
  A = 1.01*g/mole;
  G4Element* elH = new G4Element ("Hydrogen","H",Z = 1.,A);

  A = 12.011*g/mole;
  G4Element* elC = new G4Element("Carbon","C",Z = 6.,A);

  A = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",Z = 8.,A);

   // Perspex, plexiglass, lucite 
  d = 1.19*g/cm3;
  G4Material* plexiglass = new G4Material("Plexiglass",d,3);
  plexiglass->AddElement(elH,0.08);
  plexiglass->AddElement(elC,0.60);
  plexiglass->AddElement(elO,0.32);
  */

  //G4Isotope* Li6 = new G4Isotope("Li6", 3, 6);  
  //G4Isotope* Li7 = new G4Isotope("Li7", 3, 7);  
  //G4Element* Li = new G4Element("Lithium", "Li", ncomponents=2);
  //Li->AddIsotope(Li6, 4.85*perCent);
  //Li->AddIsotope(Li7, 95.15*perCent);
  //G4Material* natural_lithium = new G4Material("Lithium", 2.635*g/cm3, ncomponents=2, kStateSolid, 293.15*kelvin, 1*atmosphere);
  //natural_lithium->AddElement(Li, natoms=1);

  /*
  G4Isotope* Li7 = new G4Isotope("Li7", 3, 7);  
  G4Element* Li = new G4Element("Lithium", "Li", ncomponents=1);
  Li->AddIsotope(Li7, 100.*perCent);
  G4Material* natural_lithium = new G4Material("Lithium", 2.635*g/cm3, ncomponents=2, kStateSolid, 293.15*kelvin, 1*atmosphere);
  */

  //G4Element* Li = new G4Element("Lithium", "Li", 3., 6.941*g/mole);
  G4Element* Li = new G4Element("Lithium", "Li", 1);
  G4Isotope* Li7 = new G4Isotope("Li7", Z=3, A=7);  
  Li->AddIsotope(Li7, 100.*perCent);
  G4Material* natural_lithium = new G4Material("Lithium", 2.635*g/cm3, ncomponents=1, kStateSolid, 293.15*kelvin, 1*atmosphere);
  natural_lithium->AddElement(Li, natoms=1);

  // heavy water
  G4Element* O  = new G4Element("Oxygen", "O", 8., 16.00*g/mole);
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);  
  G4Material* D2O = new G4Material("HeavyWater", 1.11*g/cm3, ncomponents=2,
                        kStateLiquid, 293.15*kelvin, 1*atmosphere);
  D2O->AddElement(D, natoms=2);
  D2O->AddElement(O, natoms=1);

  // graphite
  G4Isotope* C12 = new G4Isotope("C12", 6, 12);  
  G4Element* C   = new G4Element("TS_C_of_Graphite","C", ncomponents=1);
  C->AddIsotope(C12, 100.*perCent);
  G4Material* graphite = 
  new G4Material("graphite", 2.27*g/cm3, ncomponents=1,
                         kStateSolid, 293*kelvin, 1*atmosphere);
  graphite->AddElement(C, natoms=1);    


  fTargetMaterial  = nistManager->FindOrBuildMaterial("Lithium");
  //fTargetMaterial  = graphite;
  //fTargetMaterial  = graphite;

  //fChamberMaterial = nistManager->FindOrBuildMaterial("Plexiglass");
  fChamberMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
  //fChamberMaterial = D2O;

  // Print materials
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  //G4Material* air  = G4Material::GetMaterial("G4_AIR");
  //G4Material* vacuum  = G4Material::GetMaterial("G4_Galactic");
  G4Material* vacuum = new G4Material("Galactic", 1, 1.01*g/mole, universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Sizes of the principal geometrical components (solids)

  G4int NbOfChambers = 4;
  G4double chamberSpacing = 10*cm; // from chamber center to center!

  G4double chamberWidth = 10.0*cm; // width of the chambers
  G4double targetLength =  3.0*mm; // full length of Target

  G4double trackerLength = (NbOfChambers+1)*chamberSpacing;

  G4double worldLength = 1.2 * (2*targetLength + trackerLength);

  G4double targetRadius  = 5.0*cm;   // Radius of Target
  targetLength = 0.5*targetLength;             // Half length of the Target
  G4double trackerSize   = 0.5*trackerLength;  // Half length of the Tracker

  // Definitions of Solids, Logical Volumes, Physical Volumes

  // World

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

  G4cout << "Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  auto worldS = new G4Box("world",                       // its name
    worldLength / 2, worldLength / 2, worldLength / 2 + 130 * cm);  // its size
  auto worldLV = new G4LogicalVolume(worldS,             // its solid
    vacuum,                                                 // its material
    "World");                                            // its name

  //  Must place the World Physical volume unrotated at (0,0,0).
  //
  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operations
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps

  // Target

  G4ThreeVector positionTarget = G4ThreeVector(0,0,0);

  auto targetS = new G4Tubs("target", 0., targetRadius, targetLength, 0. * deg, 360. * deg);
  fLogicTarget = new G4LogicalVolume(targetS, fTargetMaterial, "Target", nullptr, nullptr, nullptr);
  new G4PVPlacement(nullptr,  // no rotation
    positionTarget,           // at (x,y,z)
    fLogicTarget,             // its logical volume
    "Target",                 // its name
    worldLV,                  // its mother volume
    false,                    // no boolean operations
    0,                        // copy number
    fCheckOverlaps);          // checking overlaps

  G4cout << "Target is " << 2*targetLength/cm << " cm of "
         << fTargetMaterial->GetName() << G4endl;

  // Tracker

  G4ThreeVector positionTracker = G4ThreeVector(0,0,100*cm);

  auto trackerS = new G4Tubs("tracker", 0, trackerSize, trackerSize, 0. * deg, 360. * deg);
  auto trackerLV = new G4LogicalVolume(trackerS, vacuum, "Tracker", nullptr, nullptr, nullptr);
  new G4PVPlacement(nullptr,  // no rotation
    positionTracker,          // at (x,y,z)
    trackerLV,                // its logical volume
    "Tracker",                // its name
    worldLV,                  // its mother  volume
    false,                    // no boolean operations
    0,                        // copy number
    fCheckOverlaps);          // checking overlaps

  // Tracker segments

  // An example of Parameterised volumes
  // Dummy values for G4Tubs -- modified by parameterised volume

  auto chamberS = new G4Tubs("tracker", 0, 100 * cm, 100 * cm, 0. * deg, 360. * deg);
  fLogicChamber =
    new G4LogicalVolume(chamberS, fChamberMaterial, "Chamber", nullptr, nullptr, nullptr);

  G4double firstPosition = -trackerSize + chamberSpacing;
  G4double firstLength   = trackerLength;
  G4double lastLength    = trackerLength;

  G4VPVParameterisation* chamberParam =
    new ChamberParameterisation(
                                  NbOfChambers,   // NoChambers
                                  firstPosition,  // Z of center of first
                                  chamberSpacing, // Z spacing of centers
                                  chamberWidth,  // chamber width
                                  firstLength,    // initial length
                                  lastLength);    // final length

  // dummy value : kZAxis -- modified by parameterised volume

  new G4PVParameterised("Chamber",       // their name
                        fLogicChamber,   // their logical volume
                        trackerLV,       // Mother logical volume
                        kZAxis,          // Are placed along this axis
                        NbOfChambers,    // Number of chambers
                        chamberParam,    // The parametrisation
                        fCheckOverlaps); // checking overlaps

  G4cout << "There are " << NbOfChambers << " chambers in the tracker region. "
         << G4endl
         << "The chambers are " << chamberWidth/cm << " cm of "
         << fChamberMaterial->GetName() << G4endl
         << "The distance between chamber is " << chamberSpacing/cm << " cm"
         << G4endl;

  // Visualization attributes

  auto boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  worldLV   ->SetVisAttributes(boxVisAtt);
  fLogicTarget ->SetVisAttributes(boxVisAtt);
  trackerLV ->SetVisAttributes(boxVisAtt);

  auto chamberVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  fLogicChamber->SetVisAttributes(chamberVisAtt);

  // Example of User Limits
  //
  // Below is an example of how to set tracking constraints in a given
  // logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  G4double maxStep = 0.5*chamberWidth;
  fStepLimit = new G4UserLimits(maxStep);
  trackerLV->SetUserLimits(fStepLimit);

  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));

  // Always return the physical world

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B2/TrackerChamberSD";
  auto aTrackerSD = new TrackerSD(trackerChamberSDname, "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector( fLogicChamber,  aTrackerSD );

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fTargetMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fTargetMaterial = pttoMaterial;
        if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
        G4cout
          << G4endl
          << "----> The target is made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl
          << "-->  WARNING from SetTargetMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetChamberMaterial(G4String materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMaterial != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMaterial = pttoMaterial;
        if (fLogicChamber) fLogicChamber->SetMaterial(fChamberMaterial);
        G4cout
          << G4endl
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
