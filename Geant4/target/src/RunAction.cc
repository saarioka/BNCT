#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

namespace B2
{

RunAction::RunAction()
{
  G4RunManager::GetRunManager()->SetPrintProgress(1000000);
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  auto analysisManager = G4AnalysisManager::Instance();

  std::string runnumber = std::to_string( run->GetRunID() );

  const char* runID = std::getenv("RUN_ID");
  G4String identifier;
  if (runID != NULL) {
    identifier = runID;
    identifier = "_" + identifier;
  } else {
    identifier = "";
  }

  G4String fileName = "Run" + runnumber + identifier + ".csv";

  analysisManager->SetNtupleMerging(false);
  analysisManager->OpenFile(fileName);

  // Ntuples
  analysisManager->CreateNtuple("Kinematics", "Kinematics");
  analysisManager->CreateNtupleIColumn("Evt");
  analysisManager->CreateNtupleIColumn("HC");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("ENeutron");
  analysisManager->CreateNtupleDColumn("X");
  analysisManager->CreateNtupleDColumn("Y");
  analysisManager->CreateNtupleDColumn("Z");
  analysisManager->CreateNtupleDColumn("pX");
  analysisManager->CreateNtupleDColumn("pY");
  analysisManager->CreateNtupleDColumn("pZ");
  analysisManager->FinishNtuple();
}

void RunAction::EndOfRunAction(const G4Run* ){
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();
}

}

