#ifndef B2TallySD_h
#define B2TallySD_h 1

#include "G4VSensitiveDetector.hh"

#include "TargetHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

namespace B2
{

class TallySD : public G4VSensitiveDetector
{
  public:
    TallySD(const G4String& name,
                const G4String& hitsCollectionName);
    ~TallySD() override = default;

    // methods from base class
    void   Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void   EndOfEvent(G4HCofThisEvent* hitCollection) override;

  private:
    TargetHitsCollection* fHitsCollection = nullptr;
};

}

#endif
