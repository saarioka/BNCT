#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4CsvAnalysisReader.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* Gun) : fAction(Gun)
{
  fDir = new G4UIdirectory("/BNCT/");
  //fDir->SetGuidance("this example");

  fSelectActionCmd = new G4UIcmdWithAnInteger("/BNCT/selectGunAction", this);
  fSelectActionCmd->SetGuidance("Select primary generator action");
  fSelectActionCmd->SetGuidance("0 proton at fixed energy");
  fSelectActionCmd->SetGuidance("1 neutrons from old output file used now as input");
  fSelectActionCmd->SetParameterName("id", false);
  fSelectActionCmd->SetRange("id>=0 && id<2");

  fSelectFolderCmd = new G4UIcmdWithAString("/BNCT/selectInputFolder", this);
  fSelectFolderCmd->SetGuidance("Select input folder for primary generator action 1");
  fSelectFolderCmd->SetGuidance("For example: ../run2_2310kev/Run0");
  fSelectFolderCmd->SetGuidance("Do not include _nt and so on");
  fSelectFolderCmd->SetParameterName("folder", false);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fSelectActionCmd;
  delete fSelectFolderCmd;
  delete fDir;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fSelectActionCmd){
    fAction->SelectAction(fSelectActionCmd->GetNewIntValue(newValue));
  }
  else if (command == fSelectFolderCmd){
    //fAction->SelectFolder(fSelectFolderCmd->GetCurrentValue());
    G4CsvAnalysisReader::Instance()->SetFileName(newValue);
    G4cout << "Input file set to " << G4CsvAnalysisReader::Instance()->GetFileName() << G4endl;
  }
}
