/*
 // class: StFemtoTrack
 // It is to reduce StPicoTrack down to a smaller object
 // for the purpose of reducing memory consumption, especially
 // when having many tracks in a buffer, such as Mixed Events
 //
 // Available members:
 // 1) px, py, pz, p
 // 2) pt,
 // 3) eta, theta,
 // 4) charge
 //
 // Author: Joel Mazer
 // Affiliation: Rutgers University
 // for: the STAR collaboration
 //
*/

#include "StFemtoTrack.h"

#include "StVParticle.h"
#include "TVector3.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

#include "Riostream.h"

/// \cond CLASSIMP
ClassImp(StFemtoTrack);
/// \endcond

/**
 * Default constructor
 */
//______________________________________________________________________________________________
StFemtoTrack::StFemtoTrack() :
  TObject(),
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
//______________________________________________________________________________________________
StFemtoTrack::StFemtoTrack(Double_t px, Double_t py, Double_t pz) :
  TObject(),
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
//______________________________________________________________________________________________
StFemtoTrack::StFemtoTrack(Double_t pt, Double_t eta, Double_t phi, Double_t charge) :
  TObject(),
  fPt(pt),
  fEta(eta),
  fPhi(phi),
  fCharge(charge)
{
  //fPhi = TVector2::Phi_0_2pi(fPhi);
}
//
//______________________________________________________________________________________________
StFemtoTrack::StFemtoTrack(const StPicoTrack *track, double Bfield, TVector3 mVertex, bool prim)
{
  // primary track switch: get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if(prim) {
    // get primary track vector
    mTrkMom = track->pMom();
  } else {
    // get global track vector
    mTrkMom = track->gMom(mVertex, Bfield);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  short charge = track->charge();

  fPt = pt;
  fEta = eta;
  fPhi = phi;
  fCharge = charge;
}

/**
 * Copy constructor.
 * @param track Constant reference to copy from
 */
//______________________________________________________________________________________________
StFemtoTrack::StFemtoTrack(const StFemtoTrack& t) :
  fPt(t.fPt),
  fEta(t.fEta),
  fPhi(t.fPhi),
  fCharge(t.fCharge)
{
}

/**
 * Destructor.
 */
//______________________________________________________________________________________________
StFemtoTrack::~StFemtoTrack()
{
}

/**
 * Assignment operator
 * @param track Constant reference to copy from
 * @return A reference to this
 */
//______________________________________________________________________________________________
StFemtoTrack& StFemtoTrack::operator=(const StFemtoTrack& t)
{
  if (this != &t) {
    fPt                 = t.fPt;
    fEta                = t.fEta;
    fPhi                = t.fPhi;
    fCharge             = t.fCharge;
  }

  return *this;
}
