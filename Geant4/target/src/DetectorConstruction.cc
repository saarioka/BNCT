#include <string>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TargetSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace B2;

namespace B2b {

G4ThreadLocal G4GlobalMagFieldMessenger *DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction() { fMessenger = new DetectorMessenger(this); }

DetectorConstruction::~DetectorConstruction() {
    delete fStepLimit;
    delete fMessenger;
}

G4VPhysicalVolume *DetectorConstruction::Construct() {
    DefineMaterials();
    return DefineVolumes();
}

void DetectorConstruction::DefineMaterials() {
    G4NistManager *nistManager = G4NistManager::Instance();

    G4int ncomponents, natoms;

    G4double A; // atomic mass
    G4double Z; // atomic number
    G4double d; // density

    A = 1.01 * g / mole;
    G4Element *elH = new G4Element("Hydrogen", "H", Z = 1., A);

    A = 12.011 * g / mole;
    G4Element *elC = new G4Element("Carbon", "C", Z = 6., A);

    A = 16.00 * g / mole;
    G4Element *elO = new G4Element("Oxygen", "O", Z = 8., A);

    // Perspex, plexiglass, lucite
    d = 1.19 * g / cm3;
    G4Material *plexiglass = new G4Material("Plexiglass", d, 3);
    plexiglass->AddElement(elH, 0.08);
    plexiglass->AddElement(elC, 0.60);
    plexiglass->AddElement(elO, 0.32);

    // LiF
    G4Element *Li = new G4Element("Lithium", "Li", 2);
    G4Isotope *Li6 = new G4Isotope("Li6", Z = 3, A = 6);
    G4Isotope *Li7 = new G4Isotope("Li7", Z = 3, A = 7);
    Li->AddIsotope(Li6, 4.85 * perCent);
    Li->AddIsotope(Li7, 95.15 * perCent);

    G4Element *F = new G4Element("Fluorine", "F", 1);
    G4Isotope *F19 = new G4Isotope("F19", Z = 9, A = 19);
    F->AddIsotope(F19, 100 * perCent);

    G4Material *LiF = new G4Material("LiF", 2.635 * g / cm3, ncomponents = 2, kStateSolid, 293.15 * kelvin, 1 * atmosphere);
    LiF->AddElement(Li, natoms = 1);
    LiF->AddElement(F, natoms = 1);

    fTargetMaterial = nistManager->FindOrBuildMaterial("LiF");

    fWorldMaterial = nistManager->FindOrBuildMaterial("G4_Galactic");

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume *DetectorConstruction::DefineVolumes() {
    // Sizes of the principal geometrical components (solids)

    G4double targetLength = 1.5 * mm;      // half length of Target
    G4double targetRadius = 38.1 * mm / 2; // Radius of Target

    G4ThreeVector positionTarget = G4ThreeVector(0, 0, targetLength);

    // Definitions of Solids, Logical Volumes, Physical Volumes

    // World
    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(1.1 * m);

    G4cout << "Computed tolerance = " << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / mm << " mm" << G4endl;

    auto worldS = new G4Box("world", 60 * cm, 60 * cm, 130 * cm); // its size
    auto worldLV = new G4LogicalVolume(worldS,                    // its solid
                                       fWorldMaterial,                  // its material
                                       "World");                  // its name

    auto worldPV = new G4PVPlacement(nullptr,         // no rotation
                                     G4ThreeVector(), // at (0,0,0)
                                     worldLV,         // its logical volume
                                     "World",         // its name
                                     nullptr,         // its mother  volume
                                     false,           // no boolean operations
                                     0,               // copy number
                                     fCheckOverlaps); // checking overlaps

    // Target
    auto targetS = new G4Tubs("target", 0., targetRadius, targetLength, 0. * deg, 360. * deg);
    fLogicTarget = new G4LogicalVolume(targetS, fTargetMaterial, "Target", nullptr, nullptr, nullptr);
    new G4PVPlacement(nullptr,         // no rotation
                      positionTarget,  // at (x,y,z)
                      fLogicTarget,    // its logical volume
                      "Target",        // its name
                      worldLV,         // its mother volume
                      false,           // no boolean operations
                      0,               // copy number
                      fCheckOverlaps); // checking overlaps

    G4cout << "Target is " << fTargetMaterial->GetName() << ", " << 2 * targetLength / cm << " cm long and has radius of " << targetRadius / cm << " cm" << G4endl;

    // Visualization attributes
    auto boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    worldLV->SetVisAttributes(boxVisAtt);

    auto targetVisAtt = new G4VisAttributes(G4Colour(1, 1, 0));
    fLogicTarget->SetVisAttributes(targetVisAtt);

    // User Limits
    G4double maxStep = 0.1*cm;
    fStepLimit = new G4UserLimits(maxStep);

    // Set additional contraints on the track, with G4UserSpecialCuts
    // G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
    // trackerLV->SetUserLimits(new G4UserLimits(maxStep,
    //                                           maxLength,
    //                                           maxTime,
    //                                           minEkin));

    return worldPV;
}

void DetectorConstruction::ConstructSDandField() {
    auto targetSD = new TargetSD("TargetSD", "TargetHitCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(targetSD);
    SetSensitiveDetector(fLogicTarget, targetSD);
}

void DetectorConstruction::SetTargetMaterial(G4String materialName) {
    G4NistManager *nistManager = G4NistManager::Instance();

    G4Material *pttoMaterial = nistManager->FindOrBuildMaterial(materialName);

    if (fTargetMaterial != pttoMaterial) {
        if (pttoMaterial) {
            fTargetMaterial = pttoMaterial;
            if (fLogicTarget)
                fLogicTarget->SetMaterial(fTargetMaterial);
            G4cout << G4endl << "----> The target is made of " << materialName << G4endl;
        } else {
            G4cout << G4endl << "-->  WARNING from SetTargetMaterial : " << materialName << " not found" << G4endl;
        }
    }
}

void DetectorConstruction::SetMaxStep(G4double maxStep) {
    if ((fStepLimit) && (maxStep > 0.))
        fStepLimit->SetMaxAllowedStep(maxStep);
}

}
