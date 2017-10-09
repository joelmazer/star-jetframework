// this class is adapted from the class AliEmcalJet
//
#include "StFemtoTrack.h"
#include "StVParticle.h"
#include "StThreeVectorF.hh"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"

#include "Riostream.h"
/// \cond CLASSIMP
ClassImp(StFemtoTrack);
/// \endcond

/**
 * Default constructor
 */
StFemtoTrack::StFemtoTrack() :
  TObject(),
  //StVParticle(),//FIXME
  fPt(0),
  fEta(0),
  fPhi(0),
  fCharge(0)
{
}

/**
 * Constructor that uses the px, py, pz to define the track
 * @param px component of the track momentum
 * @param py component of the track momentum
 * @param pz component of the track momentum
 */
StFemtoTrack::StFemtoTrack(Double_t px, Double_t py, Double_t pz) :
  TObject(),
  //StVParticle(),//FIXME
  fPt(TMath::Sqrt(px * px + py* py)),
  fEta(TMath::ASinH(pz / fPt)),
  fPhi(0),
  fCharge(0)
{
  if (fPt != 0) {
    //fPhi = TVector2::Phi_0_2pi(TMath::ATan2(py, px));
  }
}

/**
 * Constructor that uses the pt, eta, phi, and q to define track: used for Mixed Events 
 * to reduce memory consumption
 * @param pt Transverse component of the track momentum
 * @param eta Pseudo-rapidity of the track
 * @param phi Azimuthal angle of the track
 * @param q Charge of the track
 */
StFemtoTrack::StFemtoTrack(Double_t pt, Double_t eta, Double_t phi, Double_t charge) :
  TObject(),
  //StVParticle(),//FIXME
  fPt(pt),
  fEta(eta),
  fPhi(phi),
  fCharge(charge)
{
  //fPhi = TVector2::Phi_0_2pi(fPhi);
}

StFemtoTrack::StFemtoTrack(const StPicoTrack *track, double Bfield, StThreeVectorF mVertex, bool prim)
{

  double pt, eta, phi;
  // primary track switch
  if(prim){
    // get primary track variables
    StThreeVectorF mPMomentum = track->pMom();
    phi = mPMomentum.phi();
    eta = mPMomentum.pseudoRapidity();
    pt = mPMomentum.perp();
  } else {
    // get global track variables
    StThreeVectorF mgMomentum = track->gMom(mVertex, Bfield);
    phi = mgMomentum.phi();
    eta = mgMomentum.pseudoRapidity();
    pt = mgMomentum.perp();
  }

  fPt = pt;
  fEta = eta;
  fPhi = phi;
  fCharge = track->charge();
}

/**
 * Copy constructor.
 * @param track Constant reference to copy from
 */
StFemtoTrack::StFemtoTrack(const StFemtoTrack& t) :
  //StVParticle(t),//FIXME
  fPt(t.fPt),
  fEta(t.fEta),
  fPhi(t.fPhi),
  fCharge(t.fCharge)
{
}

/**
 * Destructor.
 */
StFemtoTrack::~StFemtoTrack()
{
}

/**
 * Assignment operator
 * @param track Constant reference to copy from
 * @return A reference to this
 */
StFemtoTrack& StFemtoTrack::operator=(const StFemtoTrack& t)
{
  if (this != &t) {
    //StVParticle::operator=(t); //FIXME
    fPt                 = t.fPt;
    fEta                = t.fEta;
    fPhi                = t.fPhi;
    fCharge             = t.fCharge;
  }

  return *this;
}
