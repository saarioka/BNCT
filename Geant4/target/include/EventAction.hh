#ifndef B2EventAction_h
#define B2EventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"

namespace B2
{

/// Event action class

class EventAction : public G4UserEventAction
{
  public:
    EventAction() = default;
    ~EventAction() override = default;

    void  BeginOfEventAction(const G4Event* ) override;
    void    EndOfEventAction(const G4Event* ) override;
  
  private:
};

}

#endif
