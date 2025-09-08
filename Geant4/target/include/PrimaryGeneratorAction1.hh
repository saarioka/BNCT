#ifndef PrimaryGeneratorAction1_h
#define PrimaryGeneratorAction1_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4CsvAnalysisReader.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction1
{
  public:
    PrimaryGeneratorAction1(G4ParticleGun*);
    ~PrimaryGeneratorAction1() = default;

  public:
    void GeneratePrimaries(G4Event*);
    G4CsvAnalysisReader* GetAnalysisReader() { return fAnalysisReader; };

  private:
    G4ParticleGun* fParticleGun = nullptr;
    G4String fFilepattern = "";
    G4CsvAnalysisReader *fAnalysisReader = nullptr;
    G4double fE, fX, fY, fZ, fpX, fpY, fpZ;
    void InitializeAnalysisReader();
};

#endif
