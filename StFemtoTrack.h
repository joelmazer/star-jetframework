#ifndef STFEMTOTRACK_H
#define STFEMTOTRACK_H

// C++ includes
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
#include "TObject.h"
#include "TVector3.h"

#include "StVParticle.h"
#include "StMaker.h"
class StVParticle;
class StPicoTrack;
class StMaker;

class StFemtoTrack : public TObject
{
 public:
  StFemtoTrack(const StPicoTrack*, Double_t Bfield, TVector3 mVertex, bool prim, Double_t weight, Int_t trig);
  StFemtoTrack(const StPicoTrack*, double Bfield, TVector3 mVertex, bool prim);
  StFemtoTrack();
  StFemtoTrack(Double_t px, Double_t py, Double_t pz);
  StFemtoTrack(Double_t pt, Double_t eta, Double_t phi, Double_t charge);
  StFemtoTrack(const StFemtoTrack &t);
  StFemtoTrack& operator=(const StFemtoTrack &t);
  virtual ~StFemtoTrack();
  friend std::ostream &operator<<(std::ostream &in, const StFemtoTrack &t);

  // Implementation of StVParticle interface //FIXME
  Double_t          Px()                         const { return fPt*TMath::Cos(fPhi) ; }
  Double_t          Py()                         const { return fPt*TMath::Sin(fPhi) ; }
  Double_t          Pz()                         const { return fPt*TMath::SinH(fEta); }
  Double_t          Pt()                         const { return fPt                  ; }
  Double_t          P()                          const { return fPt*TMath::CosH(fEta); }
  Double_t          Phi()                        const { return fPhi   ; }
  Double_t          Theta()                      const { return 2*TMath::ATan(TMath::Exp(-fEta))      ; }
  Double_t          Eta()                        const { return fEta   ; }
  Short_t           Charge()                     const { return fCharge     ; }
  Double_t          ReWeightCorr()               const { return fReWeight ; }
  Int_t             MBTrig()                     const { return fMBTrig   ; }

  // Phi from 0 to 2*pi
  Double_t          Phi_0_2pi()                  const { return TVector2::Phi_0_2pi(fPhi); }

  void              SetVals(Double_t pt, Double_t eta, Double_t phi, Short_t q, Double_t weights, Int_t trig)
			{ fPt = pt; fEta = eta; fPhi = phi; fCharge = q; fReWeight = weights; fMBTrig = trig; }
  void              SetQ(Double_t l)                  { fCharge   = l;                     }

 protected:
  /// track transverse momentum
  Double32_t        fPt;
  /// track pseudo-rapidity
  Double32_t        fEta;           
  /// track axis azimuthal angle
  Double32_t        fPhi;               
  /// track charge
  Short_t           fCharge;           
  /// MB event reweight (and scale factor for MB5 -> MB30)
  Double32_t        fReWeight;
  /// MB event trigger (used in mixed events, tagging track by the event)
  Int_t             fMBTrig;     
  /// reference multiplicity correction (refMultCorr) for mixed event tracks
  Double32_t        fRefMultCorrMix;

 private:
  /// \cond CLASSIMP
  ClassDef(StFemtoTrack, 3);
  /// \endcond
};

std::ostream &operator<<(std::ostream &in, const StFemtoTrack &t);
#endif
