//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B2/B2b/include/TrackerHit.hh
/// \brief Definition of the B2::TrackerHit class

#ifndef B2TrackerHit_h
#define B2TrackerHit_h 1

#include "G4Allocator.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "tls.hh"
#include <map>

namespace B2 {

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

static std::map<G4int, G4String> fChamberNbToName = {{0, "Flange"}, {1, "Moderator"}, {2, "Panel"}, {3, "BertholdGas"}, {4, "Scorer1"}, {5, "Scorer2"}, {6, "Scorer3"}};
static std::map<G4String, G4int> fChamberNameToNb = {{"Flange", 0}, {"Moderator", 1}, {"Panel", 2}, {"BertholdGas", 3}, {"Scorer1", 4}, {"Scorer2", 5}, {"Scorer3", 6}};

class TrackerHit : public G4VHit {
  public:
    TrackerHit() = default;
    TrackerHit(const TrackerHit &) = default;
    ~TrackerHit() override = default;

    // operators
    TrackerHit &operator=(const TrackerHit &) = default;
    G4bool operator==(const TrackerHit &) const;

    inline void *operator new(size_t);
    inline void operator delete(void *);

    // methods from base class
    void Draw() override;
    void Print() override;

    // Set methods
    void SetParticleName(G4String name) { fParticleName = name; };
    void SetTrackID(G4int track) { fTrackID = track; };
    void SetChamberNb(G4int chamb) { fChamberNb = chamb; };
    void SetEdep(G4double de) { fEdep = de; };
    void SetE(G4double e) { fE = e; };
    void SetPos(G4ThreeVector xyz) { fPos = xyz; };

    // Get methods
    G4String GetParticleName() const { return fParticleName; };
    G4int GetTrackID() const { return fTrackID; };
    G4int GetChamberNb() const { return fChamberNb; };
    G4String GetChamberName() const { return fChamberNbToName.at(fChamberNb); };
    G4double GetEdep() const { return fEdep; };
    G4double GetE() const { return fE; };
    G4ThreeVector GetPos() const { return fPos; };
    static G4int NameToNb(G4String name) { return fChamberNameToNb.at(name); };

  private:
    G4String fParticleName = "";
    G4int fTrackID = -1;
    G4int fChamberNb = -1;
    G4double fEdep = 0.;
    G4double fE = 0.;
    G4ThreeVector fPos;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using TrackerHitsCollection = G4THitsCollection<TrackerHit>;

extern G4ThreadLocal G4Allocator<TrackerHit> *TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void *TrackerHit::operator new(size_t) {
    if (!TrackerHitAllocator)
        TrackerHitAllocator = new G4Allocator<TrackerHit>;
    return (void *)TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerHit::operator delete(void *hit) { TrackerHitAllocator->FreeSingle((TrackerHit *)hit); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

} // namespace B2

#endif
