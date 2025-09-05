#include "TargetSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

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

G4bool TargetSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  auto particleType = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  //G4cout << "Particle type: " << particleType << G4endl;
  //if (particleType != "neutron") return false;

  G4ThreeVector parentPos = aStep->GetPreStepPoint()->GetTouchableHandle()->GetTranslation();
  //4cout << "Parent pos: " << parentPos/cm << " cm" << G4endl;

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double e = aStep->GetPreStepPoint()->GetKineticEnergy();
  //if (edep==0.) return false;
  //G4cout << "Edep: " << edep/keV << " keV" << G4endl;

  auto newHit = new TargetHit();

  newHit->SetParticleName(particleType);
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetE(e);
  //newHit->SetPos(parentPos - aStep->GetPostStepPoint()->GetPosition());
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert( newHit );

  //newHit->Print();

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

