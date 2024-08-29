//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2/B2b/src/TrackerSD.cc
/// \brief Implementation of the B2::TrackerSD class

#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

namespace B2
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& name,
                     const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory*)
{
  G4int chamberNb = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
  if (chamberNb != 0) return false;

  auto particleType = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  //G4cout << "Particle type: " << particleType << G4endl;

  if (particleType != "neutron") return false;

  //auto postStep = aStep->GetPostStepPoint();
  //auto stepProcess = ((G4VProcess*)postStep->GetProcessDefinedStep())->GetProcessName();
  //G4cout << "Process: " << stepProcess << G4endl;

  G4ThreeVector parentPos = aStep->GetPreStepPoint()->GetTouchableHandle()->GetTranslation();

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double e = aStep->GetPreStepPoint()->GetKineticEnergy();
  //if (edep==0.) return false;

  auto newHit = new TrackerHit();

  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(chamberNb);
  newHit->SetEdep(edep);
  newHit->SetE(e);
  newHit->SetPos(parentPos - aStep->GetPostStepPoint()->GetPosition());

  /*
  G4cout
  << "TrackerSD:    Edep: " << newHit->GetEdep()
  << ", E: " << newHit->GetE()
  << ", x: " << newHit->GetPos().x()
  << ", y: " << newHit->GetPos().y()
  << ", z: " << newHit->GetPos().z()
  << G4endl;
  */

  fHitsCollection->insert( newHit );

  //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  G4int nofHits = fHitsCollection->entries();
  if ( nofHits > 0 ) {
    //G4cout << nofHits << " hits" << G4endl;
  }
  if ( verboseLevel>1 ) {
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ){
       auto hit = (*fHitsCollection)[i];
       if (hit->GetEdep() > 0.) {
          hit->Print();
       }
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

