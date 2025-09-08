#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class PrimaryGeneratorAction0;
class PrimaryGeneratorAction1;
class PrimaryGeneratorMessenger;

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

  public:
    void GeneratePrimaries(G4Event*) override;

  public:
    G4ParticleGun* GetParticleGun() { return fParticleGun; };

    void SelectAction(G4int i) { fSelectedAction = i; };
    G4int GetSelectedAction() { return fSelectedAction; };

    void SelectFolder(G4String act) { fSelectedFolder = act; };
    G4String GetSelectedFolder() { return fSelectedFolder; };

    PrimaryGeneratorAction0* GetAction0() { return fAction0; };
    PrimaryGeneratorAction1* GetAction1() { return fAction1; };

  private:
    G4ParticleGun* fParticleGun = nullptr;

    PrimaryGeneratorAction0* fAction0 = nullptr;
    PrimaryGeneratorAction1* fAction1 = nullptr;
    G4int fSelectedAction = 0;
    G4String fSelectedFolder = "";

    PrimaryGeneratorMessenger* fGunMessenger = nullptr;
};

#endif
