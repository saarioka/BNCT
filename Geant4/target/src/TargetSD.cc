#include "TargetSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

namespace B2
{

TargetSD::TargetSD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

void TargetSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  // TODO what is SensitiveDetectorName?
  fHitsCollection = new TargetHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

G4bool TargetSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto track = step->GetTrack();

  G4double cutoffEnergy = 1800 * keV;
  G4ParticleDefinition* particleType = track->GetDefinition();

  if (track->GetKineticEnergy() < cutoffEnergy && particleType == G4Proton::Definition()) {
    // Too low energy proton, cannot produce neutrons anymore in LiF -> kill
    //G4cout << "Killing proton with E = " << track->GetKineticEnergy()/keV << " keV" << G4endl;
    track->SetTrackStatus(fStopAndKill);
    return false;
  }

  // neutron created in proton inelastic process
  if (particleType == G4Neutron::Definition() && track->GetCreatorProcess()->GetProcessName() == "protonInelastic") {
    //G4cout << "Got a neutron!" << G4endl; 
    G4double edep = step->GetTotalEnergyDeposit();
    G4double protonE = step->GetPreStepPoint()->GetKineticEnergy();
    G4double neutronE = step->GetPostStepPoint()->GetKineticEnergy();

    //G4cout << track->GetCreatorProcess()->GetProcessName() << " " 
    //       << step->GetPreStepPoint()->GetKineticEnergy()/keV << " keV -> "
    //       << step->GetPostStepPoint()->GetKineticEnergy()/keV << " keV" << G4endl;

    auto newHit = new TargetHit();

    newHit->SetEdep(edep);
    newHit->SetProtonE(protonE);
    newHit->SetNeutronE(neutronE);
    newHit->SetPos(step->GetPreStepPoint()->GetPosition());
    newHit->SetMom(step->GetPreStepPoint()->GetMomentum());

    //newHit->Print();

    fHitsCollection->insert( newHit );

    //G4cout << step->GetPreStepPoint()->GetKineticEnergy()/keV << " keV neutron" << G4endl;
    track->SetTrackStatus(fStopAndKill);
  }

  return true;
}

void TargetSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nofHits = fHitsCollection->entries();
  if ( nofHits > 0 ) {
    //G4cout << nofHits << " hits" << G4endl;
  }
  if ( verboseLevel>1 ) {
     G4cout << G4endl
            << "-------->Hits Collection: in this event there are " << nofHits << " hits: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ){
       auto hit = (*fHitsCollection)[i];
       if (hit->GetEdep() > 0.) {
          hit->Print();
       }
     }
  }
}

}

