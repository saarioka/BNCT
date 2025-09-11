#include "PrimaryGeneratorAction0.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ParticleTable.hh"

PrimaryGeneratorAction0::PrimaryGeneratorAction0(G4ParticleGun* gun) : fParticleGun(gun)
{
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  G4double energy = 2.31 * MeV;

  fParticleGun->SetParticleEnergy(energy);
}

void PrimaryGeneratorAction0::GeneratePrimaries(G4Event* anEvent)
{
  //G4cout << "Generating primary protons with energy "
  //       << fParticleGun->GetParticleEnergy()/keV << " keV"
  //       << G4endl;

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
