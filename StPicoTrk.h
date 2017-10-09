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
     double magField, StThreeVectorD const& pVtx, StDcaGeometry const& dcaG);
  virtual ~StPicoTrk() {}

  /// track id, copied from StMuTrack, StTrack
  Int_t   id() const;
  /// primary momentum, only if track is primary with selected vertex StPicoEvent::mPrimaryVertex
  StThreeVectorF const& pMom() const;
  /// global momentum at point of DCA to StPicoEvent::mPrimaryVertex
  StThreeVectorF const& gMom() const;
  /// origin at DCA to StPicoEvent::mPrimaryVertex
  StThreeVectorF const& origin() const;
  Float_t Phi() const;
  Float_t Eta() const;
  Float_t Pt() const;
  /// global momentum at point of DCA to pVtx, B should be in kilogauss
  StThreeVectorF gMom(StThreeVectorF const& pVtx, float B) const;
  Short_t Charge() const;

  /// helix at point of DCA to StPicoEvent::mPrimaryVertex
  StPhysicalHelixD helix(float B) const;

  /** Checks whether this track is associated with a primary vertex. */
  bool isPrimary() const;

protected:
  UShort_t mId;               // track Id, copied from StMuTrack, StTrack
  StThreeVectorF mPMomentum;  // primary momentum, (0.,0.,0.) if none
  StThreeVectorF mGMomentum;  // global momentum at point of DCA to StPicoEvent::mPrimaryVertex
  StThreeVectorF mOrigin; // origin at dca to primary vertex
  Char_t   mCharge;

  ClassDef(StPicoTrk, 1)
};

inline Int_t   StPicoTrk::id() const { return mId; }
inline Float_t StPicoTrk::Pt() const { return mGMomentum.perp(); }
inline Float_t StPicoTrk::Eta() const { return mGMomentum.pseudoRapidity(); }
inline Float_t StPicoTrk::Phi() const { return mGMomentum.phi(); }
inline StThreeVectorF const& StPicoTrk::pMom() const { return mPMomentum; }
inline StThreeVectorF const& StPicoTrk::gMom() const { return mGMomentum; }
inline StThreeVectorF const& StPicoTrk::origin() const { return mOrigin; }
inline Short_t StPicoTrk::Charge() const { return static_cast<Short_t>(mCharge); }

/**
 * The default "primary" momentum is (0, 0, 0) but it is expected to have
 * a non-zero length when the track is associated with a primary vertex.
 */
inline bool StPicoTrk::isPrimary() const
{
  return mPMomentum.magnitude() > 0;
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
