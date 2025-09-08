#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorAction0.hh"
#include "PrimaryGeneratorAction1.hh"
#include "PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // default particle kinematic
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("geantino");
  fParticleGun->SetParticleDefinition(particleDefinition);

  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  fAction0 = new PrimaryGeneratorAction0(fParticleGun);
  fAction1 = new PrimaryGeneratorAction1(fParticleGun);

  // create a messenger for this class
  fGunMessenger = new PrimaryGeneratorMessenger(this);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fAction0;
  delete fAction1;
  delete fParticleGun;
  delete fGunMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch (fSelectedAction) {
    case 0:
      fAction0->GeneratePrimaries(anEvent);
      break;
    case 1:
      fAction1->GeneratePrimaries(anEvent);
      break;
   default:
      G4cerr << "Invalid generator fAction" << G4endl;
  }
}
