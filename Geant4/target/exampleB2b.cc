#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "G4ScoringManager.hh"

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC_AllHP.hh"
#include "G4ParticleHPManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4GenericBiasingPhysics.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);

  //use G4SteppingVerboseWithUnits
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // Construct the default run manager
  //
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
  G4ScoringManager* scoringManager = G4ScoringManager::GetScoringManager();

  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new B2b::DetectorConstruction());

  //G4VModularPhysicsList* physicsList = new QGSP_BIC_HP;
  G4VModularPhysicsList* physicsList = new QGSP_BIC_AllHP;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());

  runManager->SetUserInitialization(physicsList);

  G4ParticleHPManager::GetInstance()->SetSkipMissingIsotopes(true);
  G4ParticleHPManager::GetInstance()->SetDoNotAdjustFinalState(true);
  //G4ParticleHPManager::GetInstance()->SetProduceFissionFragments(true);
  //G4ParticleHPManager::GetInstance()->SetUseOnlyPhotoEvaporation(true);

  // Set user action classes
  runManager->SetUserInitialization(new B2::ActionInitialization());

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // barch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}
