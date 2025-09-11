#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include "TargetHit.hh"

namespace B2
{

void EventAction::BeginOfEventAction(const G4Event*)
{}

void EventAction::EndOfEventAction(const G4Event* event)
{
  // get number of stored trajectories
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  G4int eventID = event->GetEventID();

  auto analysisManager = G4AnalysisManager::Instance();

  for (G4int i_hc=0; i_hc<2; i_hc++)
  {
    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(i_hc);
    G4int nHit = hc->GetSize();
    if (nHit <= 0) {
      continue;
    }
    //G4cout << "Number of hits in this event: " << nHit << G4endl;

    for (G4int i=0; i<nHit; i++){
      auto hit = dynamic_cast<TargetHit*>(hc->GetHit(i));
      G4double neutronE = hit->GetNeutronE();
      G4double Edep = hit->GetEdep();
      G4ThreeVector pos = hit->GetPos();
      G4ThreeVector mom = hit->GetMom();

      analysisManager->FillNtupleIColumn(0, eventID);
      analysisManager->FillNtupleIColumn(1, hit->getHitCollection());
      analysisManager->FillNtupleDColumn(2, Edep / keV);
      analysisManager->FillNtupleDColumn(3, neutronE / keV);
      analysisManager->FillNtupleDColumn(4, pos.x() / cm);
      analysisManager->FillNtupleDColumn(5, pos.y() / cm);
      analysisManager->FillNtupleDColumn(6, pos.z() / cm);
      analysisManager->FillNtupleDColumn(7, mom.x() / keV);
      analysisManager->FillNtupleDColumn(8, mom.y() / keV);
      analysisManager->FillNtupleDColumn(9, mom.z() / keV);
      analysisManager->AddNtupleRow();
    }
  }
}
}