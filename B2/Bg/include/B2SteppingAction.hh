#ifndef B2SteppingAction_h
#define B2SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class B2EventAction;

class G4LogicalVolume;

/// Stepping action class
/// 

class B2SteppingAction : public G4UserSteppingAction
{
  public:
    B2SteppingAction(B2EventAction* eventAction);
    virtual ~B2SteppingAction();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

  private:
    B2EventAction*  fEventAction;
    G4LogicalVolume* fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
