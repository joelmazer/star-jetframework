#ifndef STJET_H
#define STJET_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// this class has been adapted from the AliEmcalJet class

#include <vector>
#include <algorithm>
#include <utility>
#include <iosfwd>

// ROOT includes
#include <TArrayI.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TString.h>

// my includes
#include "StVParticle.h"
#include "FJ_includes.h"
#include "StFJWrapper.h"

#include "StJetShapeProperties.h"

using std::vector;
namespace fastjet {
  class PseudoJet;
}

/**
 * @class StJet
 *
 * adapted from the AliEmcalJet class 
 * @author Salvatore Aiola <salvatore.aiola@yale.edu>, Yale University
 * @author Constantin Loizides <cloizides@lbl.gov>, Lawrence Berkeley National Laboratory
 *
 * This class encapsulates a jet reconstructed using the EMCal jet framework.
 * It can represent charged (tracks), neutral (EMCal clusters) or full (tracks+clusters) jet,
 * or a particle-level jet reconstructed from a Monte Carlo simulation.
 * Information contained in the class includes:
 * - reconstructed jet axis in cylindrical coordinates (eta, phi, pT);
 * - jet area (used for background subtraction);
 * - jet constituents (cluster, tracks, particles);
 * - flavor tagging;
 * - jet shape properties;
 * - matching with other reconstructed jets (e.g. detector level with particle level).
 * The class implements also a number of service function to calculate other observable,
 * such as fragmentation functions, subtracted jet momentum, etc.
 */
class StJet : public StVParticle
{
 public:
  
  StJet();
  StJet(Double_t px, Double_t py, Double_t pz);
  StJet(Double_t pt, Double_t eta, Double_t phi, Double_t m);
  StJet(const StJet &jet);
  StJet& operator=(const StJet &jet);
  virtual ~StJet();
  friend std::ostream &operator<<(std::ostream &in, const StJet &jet);
  Int_t Compare(const TObject *obj)  const;
  std::ostream &Print(std::ostream &in) const;
  TString toString() const;

  // Implementation of StVParticle interface
  Double_t          Px()                         const { return fPt*TMath::Cos(fPhi) ; }
  Double_t          Py()                         const { return fPt*TMath::Sin(fPhi) ; }
  Double_t          Pz()                         const { return fPt*TMath::SinH(fEta); }
  Double_t          Pt()                         const { return fPt                  ; }
  Double_t          P()                          const { return fPt*TMath::CosH(fEta); }
  Bool_t            PxPyPz(Double_t p[3])        const { p[0]=Px();p[1]=Py();p[2]=Pz(); return kTRUE; }
  Double_t          Xv()                         const { return 0.;      }
  Double_t          Yv()                         const { return 0.;      }
  Double_t          Zv()                         const { return 0.;      }
  Bool_t            XvYvZv(Double_t x[3])        const { x[0]=0;x[1]=0;x[2]=0         ; return kTRUE; }
  Double_t          OneOverPt()                  const { return 1./fPt ; }
  Double_t          Phi()                        const { return fPhi   ; }
  Double_t          Theta()                      const { return 2*TMath::ATan(TMath::Exp(-fEta))      ; }
  Double_t          E()                          const { Double_t p=P(); return TMath::Sqrt(fM*fM+p*p); }
  Double_t          M()                          const { return fM     ; }
  Double_t          Eta()                        const { return fEta   ; }
  Double_t          Y()                          const { Double_t e = E(); Double_t pz = Pz(); return 0.5*TMath::Log((e+pz)/(e-pz));    }
  Short_t           Charge()                     const { return 0      ; }
  Int_t             GetLabel()                   const { return fLabel ; }
  Int_t             PdgCode()                    const { return 0;       }
  const Double_t   *PID()                        const { return 0;       }

  // Other kinematic and jet properties
  Double_t          Phi_0_2pi()                  const { return TVector2::Phi_0_2pi(fPhi); }
  Double_t          Area()                       const { return fArea                    ; }
  Double_t          AreaPt()                     const { return fArea                    ; }
  Double_t          AreaEta()                    const { return fAreaEta                 ; }
  Double_t          AreaPhi()                    const { return fAreaPhi                 ; }
  Double_t          AreaE()                      const { return fAreaE                   ; }
  Int_t             ClusterAt(Int_t idx)         const { return fTowerIDs.At(idx)        ; }  // this stores ID of tower constituents - can deprecate
  Int_t             TowerAt(Int_t idx)           const { return fTowerIDs.At(idx)        ; }  // this stores ID of tower constituents
  Int_t             TrackAt(Int_t idx)           const { return fTrackIDs.At(idx)        ; }  // this stores ID of track constituents
  UShort_t          GetNumberOfClusters()        const { return fTowerIDs.GetSize()      ; }  // # of tower constituents - can deprecate
  UShort_t          GetNumberOfTowers()          const { return fTowerIDs.GetSize()      ; }  // # of tower constituents 
  UShort_t          GetNumberOfTracks()          const { return fTrackIDs.GetSize()      ; }  // # of track constituents
  UShort_t          GetNumberOfConstituents()    const { return GetNumberOfTowers()+GetNumberOfTracks(); } // towers + tracks
  Bool_t            IsMC()                       const { return (Bool_t)(MCPt() > 0)     ; }
  Bool_t            IsSortable()                 const { return kTRUE                    ; }
  Double_t          MaxNeutralPt()               const { return fMaxNPt                  ; }
  Double_t          MaxChargedPt()               const { return fMaxCPt                  ; }
  Double_t          NEF()                        const { return fNEF                     ; }
  UShort_t          Nn()                         const { return fNn                      ; }
  UShort_t          Nch()                        const { return fNch                     ; }
  UShort_t          N()                          const { return Nch()+Nn()               ; }
  Double_t          MCPt()                       const { return fMCPt                    ; }
  Double_t          MaxClusterPt()               const { return MaxNeutralPt()           ; }  // Use GetMaxClusterPt() - TODO double check this 
  Double_t          MaxTowerPt()                 const { return MaxNeutralPt()           ; }  // Use GetMaxTowerPt() - this fnc is fine
  Double_t          MaxTrackPt()                 const { return MaxChargedPt()           ; }  // Use GetMaxTrackPt()
  Double_t          MaxPartPt()                  const { return fMaxCPt < fMaxNPt ? fMaxNPt : fMaxCPt; }
  Double_t          PtSub()                      const { return fPtSub                   ; }
  Double_t          PtSubVect()                  const { return fPtSubVect               ; }

  // Background subtraction
  Double_t          PtSub(Double_t rho, Bool_t save = kFALSE)          ;
  Double_t          PtSubVect(Double_t rho, Bool_t save = kFALSE)      ;
  TLorentzVector    SubtractRhoVect(Double_t rho, Bool_t save = kFALSE);

  // Jet constituents
////  Int_t             ContainsCluster(AliVCluster *cluster, TClonesArray *clusters)  const;
  Int_t             ContainsCluster(Int_t ic)                                      const; // should deprecate this, but keep (ID - 1)
  Int_t             ContainsTower(Int_t it)                                        const; // use towers in STAR (ID - 1)
////  StVCluster      *GetLeadingCluster(TClonesArray *clusters)                      const;
////  StVTower      *GetLeadingTower(TClonesArray *towers)                      const;
  Int_t             ContainsTrack(StVParticle *track, TClonesArray *tracks)        const;
  Int_t             ContainsTrack(Int_t itr)                                       const; // (ID - 1)
////  StVParticle       *GetLeadingTrack(TClonesArray *tracks)                          const;

  // Fragmentation Function
  Double_t          GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz)  const;
  Double_t          GetZ(const StVParticle *trk )                                           const;
  Double_t          GetXi(const StVParticle *trk )                                          const;
  Double_t          GetXi(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz) const;

  // Other service methods // FIXME
  void              GetMomentum(TLorentzVector &vec)                                        const;
  Double_t          DeltaR(const StVParticle *part)                                         const;

  // Setters
  void              SetLabel(Int_t l)                  { fLabel   = l;                     }
  void              SetArea(Double_t a)                { fArea    = a;                     }
  void              SetAreaEta(Double_t a)             { fAreaEta = a;                     }
  void              SetAreaPhi(Double_t a)             { fAreaPhi = TVector2::Phi_0_2pi(a); }
  void              SetAreaE(Double_t a)               { fAreaE = a;                       }
  void              SetMaxNeutralPt(Double32_t t)      { fMaxNPt  = t;                     }
  void              SetMaxChargedPt(Double32_t t)      { fMaxCPt  = t;                     }
  void              SetNEF(Double_t nef)               { fNEF     = nef;                   }
  void              SetNumberOfClusters(Int_t n)       { fTowerIDs.Set(n);                 } // towers - can deprecate
  void              SetNumberOfTowers(Int_t n)         { fTowerIDs.Set(n);                 } // towers
  void              SetNumberOfTracks(Int_t n)         { fTrackIDs.Set(n);                 }
  void              SetNumberOfCharged(Int_t n)        { fNch = n;                         }
  void              SetNumberOfNeutrals(Int_t n)       { fNn = n;                          }
  void              SetMCPt(Double_t p)                { fMCPt = p;                        }
  void              SetPtSub(Double_t ps)              { fPtSub          = ps;             }
  void              SetPtSubVect(Double_t ps)          { fPtSubVect      = ps;             }
  void              AddClusterAt(Int_t clus, Int_t idx){ fTowerIDs.AddAt(clus, idx);       } // towers - can deprecate
  void              AddTowerAt(Int_t tow, Int_t idx)   { fTowerIDs.AddAt(tow, idx);        } // towers
  void              AddTrackAt(Int_t track, Int_t idx) { fTrackIDs.AddAt(track, idx);      }
  void              Clear(Option_t* /*option*/="");

  void              SetMaxTrackPt(Double32_t t)        { fMaxTrackPt = t;                  }
  Double_t          GetMaxTrackPt()                    const { return fMaxTrackPt;         } // ;
  void              SetMaxTowerEt(Double32_t t)          { fMaxTowerEt = t;                  }
  Double_t          GetMaxTowerEt()                    const { return fMaxTowerEt;         } // ;

  // Sorting methods
  void              SortConstituents();
  std::vector<int>  GetPtSortedTrackConstituentIndexes(TClonesArray *tracks) const;

  // Matching
  void              ResetMatching();
  void              SetClosestJet(StJet *j, Double_t d)             { fClosestJets[0] = j; fClosestJetsDist[0] = d    ; }
  void              SetSecondClosestJet(StJet *j, Double_t d)       { fClosestJets[1] = j; fClosestJetsDist[1] = d    ; }
  void              SetMatchedToClosest(UShort_t m)                 { fMatched        = 0; fMatchingType       = m    ; }
  void              SetMatchedToSecondClosest(UShort_t m)           { fMatched        = 1; fMatchingType       = m    ; }
  StJet            *ClosestJet()                              const { return fClosestJets[0]                          ; }
  Double_t          ClosestJetDistance()                      const { return fClosestJetsDist[0]                      ; }
  StJet            *SecondClosestJet()                        const { return fClosestJets[1]                          ; }
  Double_t          SecondClosestJetDistance()                const { return fClosestJetsDist[1]                      ; }
  StJet            *MatchedJet()                              const { return fMatched < 2 ? fClosestJets[fMatched] : 0; }
  UShort_t          GetMatchingType()                         const { return fMatchingType                            ; }

  // Ghosts
  void AddGhost(const Double_t dPx, const Double_t dPy, const Double_t dPz, const Double_t dE);
  Bool_t HasGhost() const                               { return fHasGhost; }
  const std::vector<TLorentzVector> GetGhosts()   const { return fGhosts  ; }

//TEST =========================
  void SetJetConstituents(std::vector<fastjet::PseudoJet> n)        { fJetConstit = n;                         }
  //const std::vector<fastjet::PseudoJet>&  GetInputVectors()    const { return fInputVectors;               }
  const std::vector<fastjet::PseudoJet>& GetJetConstituents()  const { return fJetConstit; }

  // Debug printouts
  void Print(Option_t* /*opt*/ = "") const;

  // Jet shape
  StJetShapeProperties *GetShapeProperties() const{ return fJetShapeProperties; }
  StJetShapeProperties *GetShapeProperties() { if (!fJetShapeProperties) CreateShapeProperties(); return fJetShapeProperties; }
  void CreateShapeProperties() { if (fJetShapeProperties) delete fJetShapeProperties; fJetShapeProperties = new StJetShapeProperties(); }
  
 protected:
  Bool_t            IsJetTrack(StJet *jet, Int_t itrack, Bool_t sorted = kFALSE)       const;
  Bool_t            IsJetCluster(StJet *jet, Int_t iclus, Bool_t sorted = kFALSE)      const; // should deprecate
  Bool_t            IsJetTower(StJet *jet, Int_t itow, Bool_t sorted = kFALSE)         const;

  /// Jet transverse momentum
  Double32_t        fPt;                  //[0,0,12]
  /// Jet pseudo-rapidity
  Double32_t        fEta;                 //[-1,1,12]
  /// Jet axis azimuthal angle
  Double32_t        fPhi;                 //[0,6.3,12]
  /// Jet mass
  Double32_t        fM;                   //[0,0,8]
  /// Jet Neutral Energy Fraction
  Double32_t        fNEF;                 //[0,1,8]
  /// Jet transverse area
  Double32_t        fArea;                //[0,0,12]
  /// Jet eta area
  Double32_t        fAreaEta;             //[0,0,12]
  /// Jet phi area
  Double32_t        fAreaPhi;             //[0,0,12]
  /// Jet temporal area component
  Double32_t        fAreaE;               //[0,0,12]
  /// Pt of maximum charged constituent
  Double32_t        fMaxCPt;              //[0,0,12]
  /// Pt of maximum neutral constituent
  Double32_t        fMaxNPt;              //[0,0,12]
  Double32_t        fMaxTrackPt;          // Max track pt
  Double32_t        fMaxClusterPt;        // Max cluster pt
  Double32_t        fMaxTowerEt;          // Max tower Et
  Double32_t        fMCPt;                ///<  Pt from MC particles contributing to the jet
  Int_t             fNn;                  ///<  Number of neutral constituents
  Int_t             fNch;                 ///<  Number of charged constituents
  TArrayI           fTowerIDs;            ///<  Array containing ids of tower constituents
  TArrayI           fTrackIDs;            ///<  Array containing ids of track constituents
  StJet             *fClosestJets[2];     //!<! If this is MC it contains the two closest detector level jets in order of distance and viceversa
  Double32_t        fClosestJetsDist[2];  //!<! Distance from the two closest jets
  UShort_t          fMatched;             //!<! 0 or 1 if it is matched with one of the closest jets; 2 if it is not matched
  UShort_t          fMatchingType;        //!<! Matching type
  Double_t          fPtSub;               //!<! Background subtracted pt (not stored set from outside)
  Double_t          fPtSubVect;           //!<! Background vector subtracted pt (not stored set from outside)
  UInt_t            fTriggers;            //!<! Triggers that the jet might have fired (AliVEvent::EOfflineTriggerTypes)
  Int_t             fLabel;               //!<! Label to inclusive jet for constituent subtracted jet

  Bool_t            fHasGhost;            //!<! Whether ghost particle are included within the constituents
  std::vector<TLorentzVector> fGhosts;    //!<! Vector containing the ghost particles

  std::vector<fastjet::PseudoJet> fJetConstit; //!<! Vector containing the jet constituents

  StJetShapeProperties *fJetShapeProperties; //!<! Pointer to the jet shape properties

 private:
  /**
   * @struct sort_descend
   * @brief Simple C structure to allow sorting in descending order
   */
  struct sort_descend {
    // first value of the pair is Pt and the second is entry index
    bool operator () (const std::pair<Double_t, Int_t>& p1, const std::pair<Double_t, Int_t>& p2)  { return p1.first > p2.first ; }
  };

  /// \cond CLASSIMP
  ClassDef(StJet,3);
  /// \endcond
};

std::ostream &operator<<(std::ostream &in, const StJet &jet);
#endif
