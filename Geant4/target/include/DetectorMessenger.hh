#ifndef B2bDetectorMessenger_h
#define B2bDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

namespace B2b
{

class DetectorConstruction;

/// Messenger class that defines commands for B2b::DetectorConstruction.
///
/// It implements commands:
/// - /B2/det/setTargetMaterial name
/// - /B2/det/stepMax value unit

class DetectorMessenger: public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction* );
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    DetectorConstruction*  fDetectorConstruction = nullptr;

    G4UIdirectory*         fDirectory = nullptr;
    G4UIdirectory*         fDetDirectory = nullptr;

    G4UIcmdWithAString*    fTargMatCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fStepMaxCmd = nullptr;
};

}

#endif
