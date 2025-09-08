#include "PrimaryGeneratorAction1.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction1::PrimaryGeneratorAction1(G4ParticleGun* gun) : fParticleGun(gun) {
}

void PrimaryGeneratorAction1::InitializeAnalysisReader(){
  fAnalysisReader = G4CsvAnalysisReader::Instance();

  // For example ../run2_2310kev/Run0
  //fAnalysisReader->SetFileName(fFilepattern);

  // Read ntuple
  G4int ntupleId = fAnalysisReader->GetNtuple("Kinematics");
  if ( ntupleId >= 0 ) {
    fAnalysisReader->SetNtupleDColumn("ENeutron", fE);
    fAnalysisReader->SetNtupleDColumn("X", fX);
    fAnalysisReader->SetNtupleDColumn("Y", fY);
    fAnalysisReader->SetNtupleDColumn("Z", fZ);
    fAnalysisReader->SetNtupleDColumn("pX", fpX);
    fAnalysisReader->SetNtupleDColumn("pY", fpY);
    fAnalysisReader->SetNtupleDColumn("pZ", fpZ);
  }

  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particleDefinition);
  G4cout << "Generating primary neutrons from file " << fAnalysisReader->GetFileName() << G4endl;
}

void PrimaryGeneratorAction1::GeneratePrimaries(G4Event* anEvent)
{
  if ( fAnalysisReader == nullptr ) {
    InitializeAnalysisReader();
  }

  //G4cout << "PrimaryGeneratorAction1::GeneratePrimaries: Reading next row" << G4endl;
  //G4cout << fAnalysisReader->GetFileName() << G4endl;
  //G4cout << "PrimaryGeneratorAction1::GeneratePrimaries: Read row" << G4endl;

  if ( !fAnalysisReader->GetNtupleRow() ) {
    G4cerr << "No more entries in the input file." << G4endl;
    return;
  }

  G4cout << "PrimaryGeneratorAction1::GeneratePrimaries" << G4endl;

  G4cout << "Primary vertex with E: " << fE / eV << " keV"
         << ", X: " << fX / mm << " mm"
         << ", Y: " << fY / mm << " mm"
         << ", Z: " << fZ / mm << " mm"
         << ", pX: " << fpX
         << ", pY: " << fpY
         << ", pZ: " << fpZ
         << G4endl;

  fParticleGun->SetParticleEnergy(fE * eV);
  fParticleGun->SetParticlePosition(G4ThreeVector(fX * mm, fY * mm, fZ  * mm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fpX, fpY, fpZ));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
