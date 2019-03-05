#ifndef StPicoTrk_h
#define StPicoTrk_h

#include <cmath>
#include "TObject.h"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/SystemOfUnits.h"

class StMuTrack;
class StDcaGeometry;

class StPicoTrk : public TObject
{
public:

  StPicoTrk();
  /// ctor. Note: primary track should be associated with the StPicoEvent::mPrimaryVertex
  StPicoTrk(StMuTrack const* globalTrack, StMuTrack const* primaryTrack,
     double magField, TVector3 const& pVtx, StDcaGeometry const& dcaG);
  virtual ~StPicoTrk() {}

  /// track id, copied from StMuTrack, StTrack
  Int_t   id() const;
  /// primary momentum, only if track is primary with selected vertex StPicoEvent::mPrimaryVertex
  TVector3 const& pMom() const;
  /// global momentum at point of DCA to StPicoEvent::mPrimaryVertex
  TVector3 const& gMom() const;
  /// origin at DCA to StPicoEvent::mPrimaryVertex
  TVector3 const& origin() const;
  Float_t Phi() const;
  Float_t Eta() const;
  Float_t Pt() const;
  /// global momentum at point of DCA to pVtx, B should be in kilogauss
  TVector3 gMom(TVector3 const& pVtx, float B) const;
  Short_t Charge() const;

  /// helix at point of DCA to StPicoEvent::mPrimaryVertex
  StPhysicalHelixD helix(float B) const;

  /** Checks whether this track is associated with a primary vertex. */
  bool isPrimary() const;

  void setPrimaryMomentum(Double_t px, Double_t py, Double_t pz);
  void setPrimaryMomentum(Float_t px, Float_t py, Float_t pz);
  void setPrimaryMomentum(TVector3 mom);
  void setGlobalMomentum(Double_t px, Double_t py, Double_t pz);
  void setGlobalMomentum(Float_t px, Float_t py, Float_t pz);
  void setGlobalMomentum(TVector3 mom);
  void setOrigin(Double_t x, Double_t y, Double_t z);
  void setOrigin(Float_t x, Float_t y, Float_t z);
  void setOrigin(TVector3 origin);

protected:
  UShort_t mId;           // track Id, copied from StMuTrack, StTrack
  //TVector3 mPMomentum;    // primary momentum, (0.,0.,0.) if none
  //TVector3 mGMomentum;    // global momentum at point of DCA to StPicoEvent::mPrimaryVertex
  //TVector3 mOrigin;       // origin at dca to primary vertex
  Float_t mPMomentumX;
  Float_t mPMomentumY;
  Float_t mPMomentumZ;
  Float_t mGMomentumX;
  Float_t mGMomentumY;
  Float_t mGMomentumZ;
  Float_t mOriginX;
  Float_t mOriginY;
  Float_t mOriginZ;

  Char_t   mCharge;

  ClassDef(StPicoTrk, 1)
};

inline Int_t   StPicoTrk::id() const { return mId; }
inline Float_t StPicoTrk::Pt() const { return mGMomentum.Perp(); }
inline Float_t StPicoTrk::Eta() const { return mGMomentum.PseudoRapidity(); }
inline Float_t StPicoTrk::Phi() const { return mGMomentum.Phi(); }
inline TVector3 const& StPicoTrk::pMom() const { return TVector3(mPMomentumX, mPMomentumY, mPMomentumZ); }
inline TVector3 const& StPicoTrk::gMom() const { return TVector3(mGMomentumX, mGMomentumY, mGMomentumZ); }
inline TVector3 const& StPicoTrk::origin() const { return mOrigin; }
inline Short_t StPicoTrk::Charge() const { return static_cast<Short_t>(mCharge); }

inline void StPicoTrk::setPrimaryMomentum(Double_t px, Double_t py, Double_t pz) {
  mPMomentumX = (Float_t)px; mPMomentumY = (Float_t)py; mPMomentumZ = (Float_t)pz;
}
inline void StPicoTrk::setPrimaryMomentum(Float_t px, Float_t py, Float_t pz) {
  mPMomentumX = px; mPMomentumY = py; mPMomentumZ = pz;
}
inline void StPicoTrk::setPrimaryMomentum(TVector3 mom) {
  mPMomentumX = (Float_t)mom.X(); mPMomentumY = (Float_t)mom.Y(); mPMomentumZ = (Float_t)mom.Z();
}
inline void StPicoTrk::setGlobalMomentum(Double_t px, Double_t py, Double_t pz) {
  mGMomentumX = (Float_t)px; mGMomentumY = (Float_t)py; mGMomentumZ = (Float_t)pz;
}
inline void StPicoTrk::setGlobalMomentum(Float_t px, Float_t py, Float_t pz) {
  mGMomentumX = px; mGMomentumY = py; mGMomentumZ = pz;
}
inline void StPicoTrk::setGlobalMomentum(TVector3 mom) {
  mGMomentumX = (Float_t)mom.X(); mGMomentumY = (Float_t)mom.Y(); mGMomentumZ = (Float_t)mom.Z();
}
inline void StPicoTrk::setOrigin(Double_t x, Double_t y, Double_t z) {
  mOriginX = (Float_t)x; mOriginY = (Float_t)y; mOriginZ = (Float_t)z;
}
inline void StPicoTrk::setOrigin(Float_t x, Float_t y, Float_t z) {
  mOriginX = x; mOriginY = y; mOriginZ = z;
}
inline void StPicoTrk::setOrigin(TVector3 orig) {
  mOriginX = (Float_t)orig.X(); mOriginY = (Float_t)orig.Y(); mOriginZ = (Float_t)orig.Z();
}


/**
 * The default "primary" momentum is (0, 0, 0) but it is expected to have
 * a non-zero length when the track is associated with a primary vertex.
 */
inline bool StPicoTrk::isPrimary() const
{
  return mPMomentum.Magnitude() > 0;
}

/// Return the global momentum at the dca point to the pVtx (usually it is the primary vertex.   B - magnetic field from PicoEvent::bField()
inline StThreeVectorF StPicoTrk::gMom(StThreeVectorF const& pVtx, float const B) const
{
  StPhysicalHelixD gHelix = helix(B);
  return gHelix.momentumAt(gHelix.pathLength(pVtx), B * kilogauss);
}

inline StPhysicalHelixD StPicoTrk::helix(float const B) const
{
  return StPhysicalHelixD(mGMomentum, mOrigin, B * kilogauss, static_cast<float>(Charge()));
}

#endif
