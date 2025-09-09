#include "TallySD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4CsvAnalysisReader.hh"


namespace B2
{

TallySD::TallySD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

void TallySD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  // TODO what is SensitiveDetectorName?
  fHitsCollection = new TargetHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

G4bool TallySD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto track = step->GetTrack();

  G4ParticleDefinition* particleType = track->GetDefinition();

  if (particleType != G4Neutron::Definition()) {
    G4cout << "Killing particle with E = " << track->GetKineticEnergy()/keV << " keV" << G4endl;
    track->SetTrackStatus(fStopAndKill);
    return false;
  }

  G4double neutronE = track->GetKineticEnergy();

  auto newHit = new TargetHit();

  newHit->SetNeutronE(neutronE);
  newHit->SetPos(track->GetPosition());
  newHit->SetMom(track->GetMomentum());

  newHit->Print();

  fHitsCollection->insert( newHit );

  track->SetTrackStatus(fStopAndKill);
  return true;
}

void TallySD::EndOfEvent(G4HCofThisEvent*)
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
