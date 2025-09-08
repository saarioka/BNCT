#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "PrimaryGeneratorMessenger.hh"


namespace B2
{

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // No results yet, generate protons
  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  G4double energy = 2.31 * MeV;

  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(energy);

  G4cout << "Generating primary protons with energy " << energy << G4endl;

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore.

  G4double worldZHalfLength = 0;
  G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = nullptr;
  if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
  else  {
    G4cerr << "World volume of box not found." << G4endl;
    G4cerr << "Perhaps you have changed geometry." << G4endl;
    G4cerr << "The gun will be place in the center." << G4endl;
  }

  if ( fAnalysisReader->GetNtupleRow() ) {
    G4cout << "X: " << X << ", Y: " << Y << ", Z: " << Z
            << " pX: " << pX << ", pY: " << pY << ", pZ: " << pZ
            << " E: " << E << G4endl;
    fParticleGun->SetParticleEnergy(E * keV);
    fParticleGun->SetParticlePosition(G4ThreeVector(X * mm, Y * mm, Z  * mm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(pX, pY, pZ));
  } else {
    G4cerr << "No more entries in the input file." << G4endl;
    return;
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::InitializeAnalysisReader(G4String filepattern)
{
  fAnalysisReader = G4CsvAnalysisReader::Instance();

  // For example ../run2_2310kev/Run0
  fAnalysisReader->SetFileName(filepattern);

  // Read ntuple
  G4int ntupleId = fAnalysisReader->GetNtuple("Kinematics");
  if ( ntupleId >= 0 ) {
    G4double E, X, Y, Z, pX, pY, pZ;
    fAnalysisReader->SetNtupleDColumn("ENeutron", E);
    fAnalysisReader->SetNtupleDColumn("X", X);
    fAnalysisReader->SetNtupleDColumn("Y", Y);
    fAnalysisReader->SetNtupleDColumn("Z", Z);
    fAnalysisReader->SetNtupleDColumn("pX", pX);
    fAnalysisReader->SetNtupleDColumn("pY", pY);
    fAnalysisReader->SetNtupleDColumn("pZ", pZ);
  }

  G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particleDefinition);
  G4cout << "Generating primary neutrons from file " << filepattern << G4endl;
}

