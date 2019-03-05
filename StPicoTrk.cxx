#include <limits>
#include "TMath.h"
#include "StEvent/StDcaGeometry.h"
#include "St_base/StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StPicoTrk.h"

//----------------------------------------------------------------------------------
StPicoTrk::StPicoTrk() : TObject(),
  mId(0), mPMomentumX(0), mPMomentumY(0), mPMomentumZ(0), mGMomentumX(0), mGMomentumY(0), mGMomentumZ(0),
  mOriginX(0), mOriginY(0), mOriginZ(0), mCharge(0)
{
}

//----------------------------------------------------------------------------------
StPicoTrk::StPicoTrk(StMuTrack const* const gTrk, StMuTrack const* const pTrk, double const B, TVector3 const& pVtx, StDcaGeometry const& dcaG)
  : StPicoTrk()
{
  if (!gTrk || gTrk->type() != global || (pTrk && (pTrk->type() != primary || pTrk->id() != gTrk->id())))
  {
    LOG_WARN << "Invalid arguments passed to StPicoTrack constructor. Object is default initialized" << endm;
    return;
  }

  mId = (UShort_t)gTrk->id();
  if(pTrk) {
    mPMomentumX = track.mPMomentumX;
    mPMomentumY = track.mPMomentumY;
    mPMomentumZ = track.mPMomentumZ;
  } 
  mGMomentumX = track.mGMomentumX;
  mGMomentumY = track.mGMomentumY;
  mGMomentumZ = track.mGMomentumZ;
  mOriginX = track.mOriginX;
  mOriginY = track.mOriginY;
  mOriginZ = track.mOriginZ;

  // Calculate global momentum and position at point of DCA to the pVtx
  StPhysicalHelixD gHelix = dcaG.helix();
  gHelix.moveOrigin(gHelix.pathLength(pVtx));
  mGMomentum = gHelix.momentum(B * kilogauss);
  mCharge    = (Char_t)(gTrk->charge());
}
