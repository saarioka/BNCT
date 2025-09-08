#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;

class PrimaryGeneratorMessenger : public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
    ~PrimaryGeneratorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    PrimaryGeneratorAction* fAction = nullptr;

    G4UIdirectory* fDir = nullptr;
    G4UIcmdWithAnInteger* fSelectActionCmd = nullptr;
    G4UIcmdWithAString* fSelectFolderCmd = nullptr;
};

#endif
