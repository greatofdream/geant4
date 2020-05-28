#include "B2SteppingAction.hh"
#include "B2EventAction.hh"
#include "B2bDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
// analysis
#include "g4root.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4Electron.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2SteppingAction::B2SteppingAction(B2EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B2SteppingAction::~B2SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B2SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const B2bDetectorConstruction* detectorConstruction
      = static_cast<const B2bDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume != fScoringVolume) return;
  G4Track* track = step->GetTrack();
  if(fEventAction->nu_eFlag && track->GetDefinition() == G4AntiNeutrinoE::AntiNeutrinoE()){
    fEventAction->nu_eEnergy = track->GetKineticEnergy();
    fEventAction->nu_eFlag = false;
  } 
  if (fEventAction->nu_muFlag && track->GetDefinition() == G4NeutrinoMu::NeutrinoMu())
  {
    fEventAction->nu_muEnergy = track->GetKineticEnergy();
    fEventAction->nu_muFlag = false;
  }
  if (fEventAction->eFlag && track->GetDefinition() == G4Electron::Electron())
  {
    fEventAction->eEnergy = track->GetKineticEnergy();
    fEventAction->eFlag = false;
  }
  
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  // fEventAction->AddEdep(edepStep);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

