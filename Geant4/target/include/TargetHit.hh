#ifndef B2TargetHit_h
#define B2TargetHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "tls.hh"
#include <map>

namespace B2 {

class TargetHit : public G4VHit {
  public:
    TargetHit() = default;
    TargetHit(const TargetHit &) = default;
    ~TargetHit() override = default;

    // operators
    TargetHit &operator=(const TargetHit &) = default;
    G4bool operator==(const TargetHit &) const;

    inline void *operator new(size_t);
    inline void operator delete(void *);

    // methods from base class
    void Draw() override;
    void Print() override;

    // Set methods
    void setHitCollection(G4int hc) { fHitCollection = hc; };
    void SetEdep(G4double de) { fEdep = de; };
    void SetNeutronE(G4double e) { fNeutronE = e; };
    void SetPos(G4ThreeVector xyz) { fPos = xyz; };
    void SetMom(G4ThreeVector xyz) { fMom = xyz; };

    // Get methods
    G4int getHitCollection() const { return fHitCollection; };
    G4double GetEdep() const { return fEdep; };
    G4double GetNeutronE() const { return fNeutronE; };
    G4ThreeVector GetPos() const { return fPos; };
    G4ThreeVector GetMom() const { return fMom; };

  private:
    G4int fHitCollection = -1;  // 0 = target, 1 = tally
    G4double fEdep = 0.;
    G4double fNeutronE = 0.;
    G4ThreeVector fPos;
    G4ThreeVector fMom;
};

using TargetHitsCollection = G4THitsCollection<TargetHit>;

extern G4ThreadLocal G4Allocator<TargetHit> *TargetHitAllocator;

inline void *TargetHit::operator new(size_t) {
    if (!TargetHitAllocator)
        TargetHitAllocator = new G4Allocator<TargetHit>;
    return (void *)TargetHitAllocator->MallocSingle();
}

inline void TargetHit::operator delete(void *hit) { TargetHitAllocator->FreeSingle((TargetHit *)hit); }

} // namespace B2

#endif
