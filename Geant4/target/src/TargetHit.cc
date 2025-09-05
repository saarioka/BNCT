#include "TargetHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

namespace B2
{

G4ThreadLocal G4Allocator<TargetHit>* TargetHitAllocator = nullptr;

G4bool TargetHit::operator==(const TargetHit& right) const
{
  return ( this == &right ) ? true : false;
}

void TargetHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void TargetHit::Print()
{
  G4cout
     << " Edep: " << std::setw(4) << G4BestUnit(fEdep, "Energy")
     << " ProtonE: " << std::setw(4) << G4BestUnit(fProtonE, "Energy")
     << " NeutronE: " << std::setw(4) << G4BestUnit(fNeutronE, "Energy")
     << " Position: " << std::setw(4) << G4BestUnit(fPos, "Length")
     << " Momentum: " << std::setw(4) << G4BestUnit(fMom, "Momentum")
     << G4endl;
}

}
