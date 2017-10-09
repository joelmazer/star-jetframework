#include <limits>
#include "TMath.h"
#include "StEvent/StDcaGeometry.h"
#include "St_base/StMessMgr.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StPicoTrk.h"

//----------------------------------------------------------------------------------
StPicoTrk::StPicoTrk() : TObject(),
  mId(0), mPMomentum(0., 0., 0.), mGMomentum(0., 0., 0.), mOrigin(0., 0., 0.), mCharge(0)
{
}

//----------------------------------------------------------------------------------
StPicoTrk::StPicoTrk(StMuTrack const* const gTrk, StMuTrack const* const pTrk, double const B, StThreeVectorD const& pVtx, StDcaGeometry const& dcaG)
  : StPicoTrk()
{
  if (!gTrk || gTrk->type() != global || (pTrk && (pTrk->type() != primary || pTrk->id() != gTrk->id())))
  {
    LOG_WARN << "Invalid arguments passed to StPicoTrack constructor. Object is default initialized" << endm;
    return;
  }

  mId = (UShort_t)gTrk->id();
  if(pTrk) { mPMomentum = pTrk->p(); }

  // Calculate global momentum and position at point of DCA to the pVtx
  StPhysicalHelixD gHelix = dcaG.helix();
  gHelix.moveOrigin(gHelix.pathLength(pVtx));
  mGMomentum = gHelix.momentum(B * kilogauss);
  mCharge    = (Char_t)(gTrk->charge());
}
