#ifndef B2TargetHit_h
#define B2TargetHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "tls.hh"
#include <map>

namespace B2 {

/// Tracker hit class
///
/// It defines data members to store the trackID, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fEdep, fPos

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
    void SetParticleName(G4String name) { fParticleName = name; };
    void SetTrackID(G4int track) { fTrackID = track; };
    void SetEdep(G4double de) { fEdep = de; };
    void SetE(G4double e) { fE = e; };
    void SetPos(G4ThreeVector xyz) { fPos = xyz; };

    // Get methods
    G4String GetParticleName() const { return fParticleName; };
    G4int GetTrackID() const { return fTrackID; };
    G4double GetEdep() const { return fEdep; };
    G4double GetE() const { return fE; };
    G4ThreeVector GetPos() const { return fPos; };

  private:
    G4String fParticleName = "";
    G4int fTrackID = -1;
    G4double fEdep = 0.;
    G4double fE = 0.;
    G4ThreeVector fPos;
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
