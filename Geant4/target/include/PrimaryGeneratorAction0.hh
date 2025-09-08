#ifndef PrimaryGeneratorAction0_h
#define PrimaryGeneratorAction0_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction0
{
  public:
    PrimaryGeneratorAction0(G4ParticleGun*);
    ~PrimaryGeneratorAction0() = default;

  public:
    void GeneratePrimaries(G4Event*);

  private:
    G4ParticleGun* fParticleGun = nullptr;
};

#endif
