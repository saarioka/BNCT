#ifndef B2bDetectorConstruction_h
#define B2bDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

namespace B2b
{

class DetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Set methods
    void SetTargetMaterial (G4String );
    void SetMaxStep (G4double );
    void SetCheckOverlaps(G4bool );

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    // static data members
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                         // magnetic field messenger
    // data members
    G4LogicalVolume*  fLogicTarget = nullptr;
    G4LogicalVolume*  fLogicTally = nullptr;

    G4Material*       fTargetMaterial = nullptr;
    G4Material*       fWorldMaterial = nullptr;
    G4Material*       fFlangeMaterial = nullptr;
    G4Material*       fTallyMaterial = nullptr;

    G4UserLimits* fStepLimit = nullptr; // pointer to user step limits

    DetectorMessenger* fMessenger = nullptr; // messenger

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

}

#endif
