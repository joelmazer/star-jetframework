// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// this class is a base-class for the jet framework used with 
// analysis over PicoDst's
// ################################################################

#include "StJetFrameworkPicoBase.h"

#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"

#include <sstream>
#include <fstream>

// STAR includes
#include "StThreeVectorF.hh"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"

// jet-framework STAR includes
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPlaneMaker.h" // new

// new includes
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  // NEW

// old file, kept for useful constants
#include "StPicoConstants.h"

// centrality
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"


#include "StEmcUtil/projection/StEmcPosition.h"
class StEmcPosition;

// classes
//class StJetMakerTask;

ClassImp(StJetFrameworkPicoBase)

//-----------------------------------------------------------------------------
StJetFrameworkPicoBase::StJetFrameworkPicoBase() :
  StMaker(),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),
  fCorrJetPt(kFALSE),
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection(kFALSE),
  doUseBBCCoincidenceRate(kTRUE), // kFALSE = use ZDC
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  fJetType(0),
  fMinPtJet(0.0),
  fTrackBias(0.2),
  fTowerBias(0.2),
  fJetRad(0.4),
  fEventZVtxMinCut(-40.0), fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fTrackPtMinCut(0.2), fTrackPtMaxCut(20.0),
  fTrackPhiMinCut(0.0), fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0), fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15), fTracknHitsRatio(0.52),
  fTowerEMinCut(0.2), fTowerEMaxCut(100.0),
  fTowerEtaMinCut(-1.0), fTowerEtaMaxCut(1.0),
  fTowerPhiMinCut(0.0), fTowerPhiMaxCut(2.0*TMath::Pi()),
  fLeadingJet(0), 
  fSubLeadingJet(0),
  fExcludeLeadingJetsFromFit(1.0), fTrackWeight(1),
  fTracksME(0x0),
  fJets(0x0),
  fBGJets(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  JetMaker(0x0),
  JetMakerBG(0x0),
  RhoMaker(0x0),
  EventPlaneMaker(0x0),
  grefmultCorr(0),
  refmultCorr(0),
  refmult2Corr(0),
  mOutName(""),
  mOutNameEP(""),
  mOutNameQA(""),
  fJetMakerName(""),
  fJetBGMakerName(""),
  fRhoMakerName(""),
  fRhoSparseMakerName(""),
  fEventPlaneMakerName(""),
  fRho(0x0),
  fRhoVal(0),
  fAddToHistogramsName(""),
  mEventCounter(0),
  mAllPVEventCounter(0),
  mInputEventCounter(0)
{

}

//-----------------------------------------------------------------------------
StJetFrameworkPicoBase::StJetFrameworkPicoBase(const char* name) :
  StMaker(name),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),
  fCorrJetPt(kFALSE),
  fCentralityDef(4), //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection(kFALSE),
  doUseBBCCoincidenceRate(kTRUE), // kFALSE = use ZDC
  fMinPtJet(0.0),
  fTrackBias(0.2),
  fTowerBias(0.2),
  fJetRad(0.4),
  fEventZVtxMinCut(-40.0), fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fTrackPtMinCut(0.2), fTrackPtMaxCut(20.0),
  fTrackPhiMinCut(0.0), fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0), fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15), fTracknHitsRatio(0.52), 
  fTowerEMinCut(0.2), fTowerEMaxCut(100.0),
  fTowerEtaMinCut(-1.0), fTowerEtaMaxCut(1.0),
  fTowerPhiMinCut(0.0), fTowerPhiMaxCut(2.0*TMath::Pi()),
  fLeadingJet(0), 
  fSubLeadingJet(0),
  fExcludeLeadingJetsFromFit(1.0), fTrackWeight(1),
  fTracksME(0x0),
  fJets(0x0),
  fBGJets(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  JetMaker(0x0),
  JetMakerBG(0x0),
  RhoMaker(0x0),
  EventPlaneMaker(0x0),
  grefmultCorr(0),
  refmultCorr(0),
  refmult2Corr(0),
  mOutName(""),
  mOutNameEP(""),
  mOutNameQA(""),
  fJetMakerName(""),
  fJetBGMakerName(""),
  fRhoMakerName(""),
  fRhoSparseMakerName(""),
  fEventPlaneMakerName(""),
  fRho(0x0),
  fRhoVal(0),
  fAddToHistogramsName(""),
  mEventCounter(0),
  mAllPVEventCounter(0),
  mInputEventCounter(0)
{

}

//----------------------------------------------------------------------------- 
StJetFrameworkPicoBase::~StJetFrameworkPicoBase()
{ /*  */
  // destructor
}

//-----------------------------------------------------------------------------
Int_t StJetFrameworkPicoBase::Init() {
  fAddToHistogramsName = "";

/*
  //AddBadTowers( TString( getenv("STARPICOPATH" )) + "/badTowerList_y11.txt");
  // Add dead + bad tower lists
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        //AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_BadTowers.txt");
        AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_AltBadTowers.txt");
        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2014_DeadTowers.txt");
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2016_BadTowers.txt");
        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Y2016_DeadTowers.txt");
        break;

    default :
      AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_BadTowers.txt");
      AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_DeadTowers.txt");
  }
*/

/*
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  // initialize centrality correction
  // switch on Run Flag to look for firing trigger specifically requested for given run period
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        // this is the default for Run14
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();        
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        switch(fCentralityDef) {      
          case StJetFrameworkPicoBase::kgrefmult :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P16id :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
              break;
          case StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx : 
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx();
              break;
          case StJetFrameworkPicoBase::kgrefmult_VpdMB30 : 
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMB30();
              break;
          default:
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
        }
        break;  // added May20

    default :
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  }
*/

  refmultCorr = CentralityMaker::instance()->getRefMultCorr(); // OLD
  refmult2Corr = CentralityMaker::instance()->getRefMult2Corr();  // OLD 

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StJetFrameworkPicoBase::Finish() { 
  //  Summarize the run.

  return kStOK;
}

//----------------------------------------------------------------------------- 
// OLD user code says: //  Called every event after Make(). 
void StJetFrameworkPicoBase::Clear(Option_t *opt) {

}

// 
//  Get the event-wise rho value
//_________________________________________________________________________________________
double StJetFrameworkPicoBase::GetRhoValue(TString fRhoMakerNametemp)
{
  // get RhoMaker from event: old names "StRho_JetsBG", "OutRho", "StMaker#0"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerNametemp));
  const char *fRhoMakerNameCh = fRhoMakerNametemp;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;
  }

  // get rho/area       fRho->ls("");
  fRhoVal = fRho->GetVal();

  return fRhoVal;
}

//________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetCentBin(Int_t cent, Int_t nBin) const
{  // Get centrality bin.
  Int_t centbin = -1;

  if(nBin == 16) { centbin = nBin - 1 - cent; }
  if(nBin == 9)  { centbin = nBin - 1 - cent; }

  return centbin;
}

// this function generate a jet name based on input
TString StJetFrameworkPicoBase::GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray* partCont, TClonesArray* clusCont, TString tag)
{
  TString algoString;
  switch (jetAlgo) {
      case kt_algorithm:
        algoString = "KT";
        break;
      case antikt_algorithm:
        algoString = "AKT";
        break;
      default:
        ::Warning("StJetFrameworkPicoBase::GenerateJetName", "Unknown jet finding algorithm '%d'!", jetAlgo);
        algoString = "";
  }

  TString typeString;
  switch (jetType) {
      case kFullJet:
        typeString = "Full";
        break;
      case kChargedJet:
        typeString = "Charged";
        break;
      case kNeutralJet:
        typeString = "Neutral";
        break;
  }

  TString radiusString = TString::Format("R%03.0f", radius*100.0);

  TString trackString;
  if (jetType != kNeutralJet && partCont) {
    trackString = "_" + TString(partCont->GetTitle());
  }

  TString clusterString;
  if (jetType != kChargedJet && clusCont) {
    clusterString = "_" + TString(clusCont->GetTitle());
  }

  TString recombSchemeString;
  switch (recoScheme) {
      case E_scheme:
        recombSchemeString = "E_scheme";
        break;
      case pt_scheme:
        recombSchemeString = "pt_scheme";
        break;
      case pt2_scheme:
        recombSchemeString = "pt2_scheme";
        break;
      case Et_scheme:
        recombSchemeString = "Et_scheme";
        break;
      case Et2_scheme:
        recombSchemeString = "Et2_scheme";
        break;
      case BIpt_scheme:
        recombSchemeString = "BIpt_scheme";
        break;
      case BIpt2_scheme:
        recombSchemeString = "BIpt2_scheme";
        break;
      case external_scheme:
        recombSchemeString = "ext_scheme";
        break;
      default:
        ::Error("StJetFrameworkPicoBase::GenerateJetName", "Recombination %d scheme not recognized.", recoScheme);
  }

  TString name = TString::Format("%s_%s%s%s%s%s_%s",
      tag.Data(), algoString.Data(), typeString.Data(), radiusString.Data(), trackString.Data(), clusterString.Data(), recombSchemeString.Data());

  return name;
}

//________________________________________________________________________
Double_t StJetFrameworkPicoBase::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi<-0.5*TMath::Pi())  dphi+=2.*TMath::Pi();
  if(dphi>3./2.*TMath::Pi()) dphi-=2.*TMath::Pi();

  // test
  if( dphi < -1.*TMath::Pi()/2 || dphi > 3.*TMath::Pi()/2 )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}

//_________________________________________________________________________
Double_t StJetFrameworkPicoBase::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{ // function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ){
    dphi = dphi + 1*TMath::Pi();
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test
  if( dphi < 0 || dphi > TMath::Pi()/2 ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dEP not in range [0, 0.5*Pi]!"<<" jetAng: "<<jetAng<<"  EPang: "<<EPAng<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

/*
//_________________________________________________
TClonesArray* StJetFrameworkPicoBase::CloneAndReduceTrackList(TClonesArray* tracksME)
{
  // clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
  TClonesArray* tracksClone = new TClonesArray("StPicoTrack");
  //tracksClone->SetName("tracksClone");
  //tracksClone->SetOwner(kTRUE);

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3 vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  Int_t nMixTracks = mPicoDst->numberOfTracks();
  Int_t iterTrk = 0;
  Double_t phi, eta, px, py, pt, pz, p, charge;
  const double pi = 1.0*TMath::Pi();
  for (Int_t i=0; i<nMixTracks; i++) { 
    // get tracks
    StPicoTrack* trk = mPicoDst->track(i);
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // add tracks passing cuts to tracksClone array
    ((*tracksClone)[iterTrk]) = trk;
    ++iterTrk;
  } // end of looping through tracks

  return tracksClone;
}
*/

//
// Track Quality Cuts
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert) {
  // constants: assume neutral pion mass
  //double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0*TMath::Pi();

  // primary track switch
  // get momentum vector of track - global or primary track
  StThreeVectorF mTrkMom;
  if(doUsePrimTracks) {
    if(!(trk->isPrimary())) return kFALSE; // check if primary

    // get primary track vector
    mTrkMom = trk->pMom();
  } else {
    // get global track vector
    mTrkMom = trk->gMom(Vert, B);
  }

  // track variables
  double pt = mTrkMom.perp();
  double phi = mTrkMom.phi();
  double eta = mTrkMom.pseudoRapidity();
  //double px = mTrkMom.x();
  //double py = mTrkMom.y();
  //double pz = mTrkMom.z();
  //double p = mTrkMom.mag();
  //double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
  //short charge = trk->charge();
  double dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // do pt cut here to accommadate either type
  if(doUsePrimTracks) { // primary  track
    if(pt < fTrackPtMinCut) return kFALSE;
  } else { // global track
    if(pt < fTrackPtMinCut) return kFALSE;
  }

  // jet track acceptance cuts now - after getting 3vector - hardcoded
  if(pt > fTrackPtMaxCut) return kFALSE; // 20.0 STAR, 100.0 ALICE
  if((eta < fTrackEtaMinCut) || (eta > fTrackEtaMaxCut)) return kFALSE;
  if(phi < 0)    phi += 2*pi;
  if(phi > 2*pi) phi -= 2*pi;
  if((phi < fTrackPhiMinCut) || (phi > fTrackPhiMaxCut)) return kFALSE;
    
  // additional quality cuts for tracks
  if(dca > fTrackDCAcut)            return kFALSE;
  if(nHitsFit < fTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track
  return kTRUE;
}

/*
//
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::AcceptTower(StPicoBTowHit *tower, StThreeVectorF Vertex) {
  // get EMCal position
  StEmcPosition *mPosition = new StEmcPosition();

  // constants:
  double pi = 1.0*TMath::Pi();

  // tower ID
  int towerID = tower->id();

  // make sure some of these aren't still in event array
  if(towerID < 0) return kFALSE; 

  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  StThreeVectorF towerPosition = mPosition->getPosFromVertex(mVertex, towerID);
  double phi = towerPosition.phi();
  if(phi < 0)    phi += 2.0*pi;
  if(phi > 2*pi) phi -= 2.0*pi;
  double eta = towerPosition.pseudoRapidity();
  int towerADC = tower->adc();
  double towerEunCorr = tower->energy();  // uncorrected energy

  // check for bad (and dead) towers
  bool TowerOK = IsTowerOK(towerID);      // kTRUE means GOOD
  bool TowerDead = IsTowerDead(towerID);  // kTRUE means BAD
  if(!TowerOK)  { return kFALSE; }
  if(TowerDead) { return kFALSE; }

  // jet track acceptance cuts now - after getting 3vector - hardcoded
  if((eta < fTowerEtaMinCut) || (eta > fTowerEtaMaxCut)) return kFALSE;
  if(phi < 0)    phi+= 2*pi;
  if(phi > 2*pi) phi-= 2*pi;
  if((phi < fTowerPhiMinCut) || (phi > fTowerPhiMaxCut)) return kFALSE;

  // passed all above cuts - keep tower and fill input vector to fastjet
  return kTRUE;
}
*/

//________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetReactionPlane() { 
  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3-vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();

  //if(mVerbose)cout << "----------- In GetReactionPlane() -----------------" << endl;
  TVector2 mQ;
  double mQx = 0., mQy = 0.;
  int order = 2;
  int n = mPicoDst->numberOfTracks();
  double pi = 1.0*TMath::Pi();

  // leading jet check and removal
  float excludeInEta = -999;
  ///fLeadingJet = GetLeadingJet();  // FIXME
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from ep estimate
    if(fLeadingJet) excludeInEta = fLeadingJet->Eta();
  }

  // loop over tracks
  for(int i = 0; i < n; i++) {
    StPicoTrack* track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // primary track switch
    // get momentum vector of track - global or primary track
    StThreeVectorF mTrkMom;
    if(doUsePrimTracks) {
      if(!(track->isPrimary())) return kFALSE; // check if primary
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.perp();
    double phi = mTrkMom.phi();
    double eta = mTrkMom.pseudoRapidity();

    // do pt cut here to accommadate either type
    if(doUsePrimTracks) { // primary  track
      if(pt < fTrackPtMinCut) continue;
    } else { // global track
      if(pt < fTrackPtMinCut) continue;
    }

    // more acceptance cuts now - after getting 3vector - hardcoded for now
    if(pt > 5.0) continue;   // 100.0
    if((1.0*TMath::Abs(eta)) > 1.0) continue;
    if(phi < 0)    phi+= 2*pi;
    if(phi > 2*pi) phi-= 2*pi;
    if((phi < 0) || (phi > 2*pi)) continue;

    // check for leading jet removal
    if(fExcludeLeadingJetsFromFit > 0 && ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) || (TMath::Abs(eta) - fJetRad - 1.0 ) > 0 )) continue;

    // configure track weight when performing Q-vector summation
    double trackweight;
    if(fTrackWeight == kNoWeight) {
      trackweight = 1.0;
    } else if(fTrackWeight == kPtLinearWeight) {
      trackweight = pt;
    } else if(fTrackWeight == kPtLinear2Const5Weight) {
      if(pt <= 2.0) trackweight = pt;
      if(pt > 2.0)  trackweight = 2.0;
    } else {
      // nothing choosen, so don't use weight
      trackweight = 1.0;
    }

    // sum up q-vectors
    mQx += trackweight * cos(phi * order);
    mQy += trackweight * sin(phi * order);
  }
 
  mQ.Set(mQx, mQy);
  double psi= mQ.Phi() / order;
  //return psi*180/pi;  // converted to degrees
  return psi;
}

// _____________________________________________________________________________________________
StJet* StJetFrameworkPicoBase::GetLeadingJet(TString fJetMakerNametemp, StRhoParameter* eventRho) {
  // return pointer to the highest pt jet (before background subtraction) within acceptance
  // only rudimentary cuts are applied on this level, hence the implementation outside of
  // the framework

  // ================= JetMaker ================ //
  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerNametemp));
  const char *fJetMakerNameCh = fJetMakerNametemp;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return 0x0;
  }

  // if we have JetMaker, get jet collection associated with it
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return 0x0;
  }

  if(fJets) {
    // get number of jets, initialize pt and leading jet pointer
    Int_t iJets(fJets->GetEntriesFast());
    Double_t pt(0);
    StJet* leadingJet(0x0);

    if(!eventRho) { // no rho parameter provided
      // loop over jets
      for(Int_t i(0); i < iJets; i++) {
        StJet* jet = static_cast<StJet*>(fJets->At(i));
        //if(!AcceptJet(jet)) continue;
        if(jet->Pt() > pt) {
          leadingJet = jet;
          pt = leadingJet->Pt();
        }
      }
      return leadingJet;
    } else { // rho parameter provided
      // return leading jet after background subtraction
      //Double_t rho(0);
      double fRhoValtemp = eventRho->GetVal(); // test

      // loop over jets
      for(Int_t i(0); i < iJets; i++) {
        StJet* jet = static_cast<StJet*>(fJets->At(i));
        //if(!AcceptJet(jet)) continue;
        if((jet->Pt() - jet->Area()*fRhoValtemp) > pt) {
          leadingJet = jet;
          pt = (leadingJet->Pt()-jet->Area()*fRhoValtemp);
        }
      }
      return leadingJet;
    }
  }

  return 0x0;
}

// _____________________________________________________________________________________________
StJet* StJetFrameworkPicoBase::GetSubLeadingJet(TString fJetMakerNametemp, StRhoParameter* eventRho) {
  // return pointer to the second highest pt jet (before background subtraction) within acceptance
  // only rudimentary cuts are applied on this level, hence the implementation outside of the framework

  // ================= JetMaker ================ //
  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerNametemp));
  const char *fJetMakerNameCh = fJetMakerNametemp;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return 0x0;
  }

  // if we have JetMaker, get jet collection associated with it
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return 0x0;
  }

  if(fJets) {
    // get number of jets 
    Int_t iJets(fJets->GetEntriesFast());

    // initialize pt's
    Double_t pt(0);
    Double_t pt2(0);

    // leading and subleading jet pointers
    StJet* leadingJet(0x0);
    StJet* subleadingJet(0x0);

    if(!eventRho) { // no rho parameter provided
      // loop over jets
      for(Int_t i(0); i < iJets; i++) {
        StJet* jet = static_cast<StJet*>(fJets->At(i));
        // leading jet
        if(jet->Pt() > pt) {
          leadingJet = jet;
          pt = leadingJet->Pt();
        }

        // subleading jet
        if((jet->Pt() > pt2) && (jet->Pt() < pt)) {
          subleadingJet = jet;
          pt2 = subleadingJet->Pt();
        }

      }
      return subleadingJet;

    } else { // rho parameter provided
      // return leading jet after background subtraction
      //Double_t rho(0);
      double fRhoValtemp = eventRho->GetVal(); // test

      // loop over jets
      for(Int_t i(0); i < iJets; i++) {
        StJet* jet = static_cast<StJet*>(fJets->At(i));
        // leading jet
        if((jet->Pt() - jet->Area()*fRhoValtemp) > pt) {
          leadingJet = jet;
          pt = (leadingJet->Pt()-jet->Area()*fRhoValtemp);
        }

        // subleading jet
        if((jet->Pt() > pt2) && (jet->Pt() < pt)) {
          subleadingJet = jet;
          pt2 = subleadingJet->Pt();
        }
      }
      return leadingJet;
    }
  }

  return 0x0;
}


//__________________________________________________________________________________________
Int_t StJetFrameworkPicoBase::EventCounter() {
  mEventCounter++;

  return mEventCounter;
}

//__________________________________________________________________________________________
Bool_t StJetFrameworkPicoBase::SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut) {
  // this function is written to cut on centrality in a task for a given range
  // STAR centrality is written in re-verse binning (0,15) where lowest bin is highest centrality
  // -- this is a very poor approach
  // -- Solution, Use:  Int_t StJetFrameworkPicoBase::GetCentBin(Int_t cent, Int_t nBin) const
  //    in order remap the centrality to a more appropriate binning scheme
  //    (0, 15) where 0 will correspond to the 0-5% bin
  // -- NOTE: Use the 16 bin scheme instead of 9, its more flexible
  // -- Example usage:
  //    Int_t cent16 = grefmultCorr->getCentralityBin16();
  //    Int_t centbin = GetCentBin(cent16, 16);
  //    if(!SelectAnalysisCentralityBin(centbin, StJetFrameworkPicoBase::kCent3050)) return kStWarn; (or StOk to suppress warnings)
  //
  // other bins can be added if needed...

  Bool_t doAnalysis = kFALSE; // set false by default, to make sure user chooses an available bin

  // switch on bin selection
  switch(fCentralitySelectionCut) {
    case kCent010 :  // 0-10%
      if((centbin>-1) && (centbin<2)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent020 :  // 0-20%
      if((centbin>-1) && (centbin<4)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }  
      break;

    case kCent1020 : // 10-20%
      if((centbin>1) && (centbin<4)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent1030 : // 10-30%
      if((centbin>1) && (centbin<6)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent1040 : // 10-40%
      if((centbin>1) && (centbin<8)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent2030 : // 20-30%
      if((centbin>3) && (centbin<6)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent2040 : // 20-40%
      if((centbin>3) && (centbin<8)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent2050 : // 20-50%
      if((centbin>3) && (centbin<10)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent2060 : // 20-60%
      if((centbin>3) && (centbin<12)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent3050 : // 30-50%
      if((centbin>5) && (centbin<10)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent3060 : // 30-60%
      if((centbin>5) && (centbin<12)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent4060 : // 40-60%
      if((centbin>7) && (centbin<12)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent4070 : // 40-70%
      if((centbin>7) && (centbin<14)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent4080 : // 40-80%
      if((centbin>7) && (centbin<16)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent5080 : // 50-80%
      if((centbin>9) && (centbin<16)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case kCent6080 : // 60-80%
      if((centbin>11) && (centbin<16)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    default : // wrong entry
      doAnalysis = kFALSE;

  }

  return doAnalysis;
}

// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i = 0; i < elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }
  //cout<<"elems: "<<elems<<"  match: "<<match<<endl;

  return match;
}

//_________________________________________________________________________
Bool_t StJetFrameworkPicoBase::CheckForMB(int RunFlag, int type) {
  // Run14 triggers:
  int arrMB_Run14[] = {450014};
  int arrMB30_Run14[] = {450010, 450020};
  int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
  // additional 30: 4, 5, 450201, 450202, 450211, 450212

  // Run16 triggers:
  int arrMB_Run16[] = {520021};
  int arrMB5_Run16[] = {520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
  int arrMB10_Run16[] = {520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};

  // run flag selection to check for MB firing
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        switch(type) { 
          case kRun14main :
              if((DoComparison(arrMB_Run14, sizeof(arrMB_Run14)/sizeof(*arrMB_Run14)))) { return kTRUE; }
              break;
          case kVPDMB5 :
              if((DoComparison(arrMB5_Run14, sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14)))) { return kTRUE; }
              break;
          case kVPDMB30 :
              if((DoComparison(arrMB30_Run14, sizeof(arrMB30_Run14)/sizeof(*arrMB30_Run14)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run14, sizeof(arrMB_Run14)/sizeof(*arrMB_Run14)))) { return kTRUE; }
        }
        break; // added May20

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        switch(type) {
          case kRun16main :
              if((DoComparison(arrMB_Run16, sizeof(arrMB_Run16)/sizeof(*arrMB_Run16)))) { return kTRUE; }
              break;
          case kVPDMB5 :
              if((DoComparison(arrMB5_Run16, sizeof(arrMB5_Run16)/sizeof(*arrMB5_Run16)))) { return kTRUE; }
              break;
          case kVPDMB10 :
              if((DoComparison(arrMB10_Run16, sizeof(arrMB10_Run16)/sizeof(*arrMB10_Run16)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run16, sizeof(arrMB_Run16)/sizeof(*arrMB_Run16)))) { return kTRUE; }
        }
        break; // added May20
  } // RunFlag switch

  // return status
  return kFALSE;
} // MB function

//
// check to see if the event was EMC triggered for High Towers
//____________________________________________________________________________
Bool_t StJetFrameworkPicoBase::CheckForHT(int RunFlag, int type) {

  // Run14 triggers:
  int arrHT1_Run14[] = {450201, 450211, 460201};
  int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
  int arrHT3_Run14[] = {450203, 450213, 460203};

  // Run16 triggers:
  int arrHT1_Run16[] = {520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
  int arrHT2_Run16[] = {530202, 540203};
  int arrHT3_Run16[] = {520203, 530213};

  // run flag selection to check for MB firing
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        switch(type) {
          case kIsHT1 :
              if((DoComparison(arrHT1_Run14, sizeof(arrHT1_Run14)/sizeof(*arrHT1_Run14)))) { return kTRUE; }
              break;
          case kIsHT2 :
              if((DoComparison(arrHT2_Run14, sizeof(arrHT2_Run14)/sizeof(*arrHT2_Run14)))) { return kTRUE; }
              break;
          case kIsHT3 :
              if((DoComparison(arrHT3_Run14, sizeof(arrHT3_Run14)/sizeof(*arrHT3_Run14)))) { return kTRUE; }
              break;
          default :  // default to HT2
              if((DoComparison(arrHT2_Run14, sizeof(arrHT2_Run14)/sizeof(*arrHT2_Run14)))) { return kTRUE; }
        }
        break; // added May20

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        switch(type) {
          case kIsHT1 :
              if((DoComparison(arrHT1_Run16, sizeof(arrHT1_Run16)/sizeof(*arrHT1_Run16)))) { return kTRUE; }
              break;
          case kIsHT2 :
              if((DoComparison(arrHT2_Run16, sizeof(arrHT2_Run16)/sizeof(*arrHT2_Run16)))) { return kTRUE; }
              break;
          case kIsHT3 :
              if((DoComparison(arrHT3_Run16, sizeof(arrHT3_Run16)/sizeof(*arrHT3_Run16)))) { return kTRUE; }
              break;
          default :  // Run16 only has HT1's
              if((DoComparison(arrHT1_Run16, sizeof(arrHT1_Run16)/sizeof(*arrHT1_Run16)))) { return kTRUE; }
        }
        break; // added May20

  } // RunFlag switch

  return kFALSE;
}

//________________________________________________________________________________________________________
Bool_t StJetFrameworkPicoBase::GetMomentum(StThreeVectorF &mom, const StPicoBTowHit* tower, Double_t mass, StPicoEvent *PicoEvent) const {
  // initialize Emc position objects
  StEmcPosition *Position = new StEmcPosition();

  // vertex components
  StThreeVectorF fVertex = PicoEvent->primaryVertex();
//  double xVtx = fVertex.x();
//  double yVtx = fVertex.y();
//  double zVtx = fVertex.z();

  // get mass, E, P, ID
  if(mass < 0) mass = 0;
  Double_t energy = tower->energy();
  Double_t p = TMath::Sqrt(energy*energy - mass*mass);
  int towerID = tower->id();

  // get tower position
  StThreeVectorF towerPosition = Position->getPosFromVertex(fVertex, towerID);
  double posX = towerPosition.x();
  double posY = towerPosition.y();
  double posZ = towerPosition.z();

  // shouldn't need correction with above method
//  posX-=xVtx;
//  posY-=yVtx;
//  posZ-=zVtx;

  // get r, set position components
  Double_t r = TMath::Sqrt(posX*posX + posY*posY + posZ*posZ) ;
  if(r > 1e-12) {
    mom.setX( p*posX/r );
    mom.setY( p*posY/r );
    mom.setZ( p*posZ/r );
    // energy) ;
  } else { return kFALSE; }
  return kTRUE;
}

/*
//____________________________________________________________________________________________
Bool_t StJetFrameworkPicoBase::IsTowerOK( Int_t mTowId ){
  //if( badTowers.size()==0 ){
  if( badTowers.empty() ){
    // maybe change class if calling FROM base
    __ERROR("StJetFrameworkPicoBase::IsTowerOK: WARNING: You're trying to run without a bad tower list. If you know what you're doing, deactivate this throw and recompile.");
    throw ( -1 );
  }
  if( badTowers.count( mTowId )>0 ){
    __DEBUG(9, Form("Reject. Tower ID: %d", mTowId));
    return kFALSE;
  } else {
    __DEBUG(9, Form("Accept. Tower ID: %d", mTowId));
    return kTRUE;
  }
}

//____________________________________________________________________________________________
Bool_t StJetFrameworkPicoBase::IsTowerDead( Int_t mTowId ){
  //if( deadTowers.size()==0 ){
  if( deadTowers.empty() ){
    // maybe change class if calling FROM base
    __ERROR("StJetFrameworkPicoBase::IsTowerDead: WARNING: You're trying to run without a dead tower list. If you know what you're doing, deactivate this throw and recompile.");
    throw ( -1 );
  }
  if( deadTowers.count( mTowId )>0 ){
    __DEBUG(9, Form("Reject. Tower ID: %d", mTowId));
    return kTRUE;
  } else {
    __DEBUG(9, Form("Accept. Tower ID: %d", mTowId));
    return kFALSE;
  }
}

//____________________________________________________________________________
void StJetFrameworkPicoBase::ResetBadTowerList( ){
  badTowers.clear();
}

// Add bad towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
Bool_t StJetFrameworkPicoBase::AddBadTowers(TString csvfile){
  // open infile
  std::string line;
  std::ifstream inFile ( csvfile );

  __DEBUG(2, Form("Loading bad towers from %s", csvfile.Data()) );

  if( !inFile.good() ) {
    __WARNING(Form("Can't open %s", csvfile.Data()) );
    return kFALSE;
  }
 
  while(std::getline (inFile, line) ){
    if( line.size()==0 ) continue; // skip empty lines
    if( line[0] == '#' ) continue; // skip comments

    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if(ientry) {
        badTowers.insert( ientry );
        __DEBUG(2, Form("Added bad tower # %d", ientry));
      }
    }
  }
 
  return kTRUE;
}

// Add dead towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
Bool_t StJetFrameworkPicoBase::AddDeadTowers(TString csvfile){
  // open infile
  std::string line;
  std::ifstream inFile ( csvfile );

  __DEBUG(2, Form("Loading bad towers from %s", csvfile.Data()) );

  if( !inFile.good() ) {
    __WARNING(Form("Can't open %s", csvfile.Data()) );
    return kFALSE;
  }

  while(std::getline (inFile, line) ){
    if( line.size()==0 ) continue; // skip empty lines
    if( line[0] == '#' ) continue; // skip comments

    std::istringstream ss( line );
    while( ss ){
      std::string entry;
      std::getline( ss, entry, ',' );
      int ientry = atoi(entry.c_str());
      if(ientry) {
        deadTowers.insert( ientry );
        __DEBUG(2, Form("Added bad tower # %d", ientry));
      }
    }
  }

  return kTRUE;
}

//____________________________________________________________________________
void StJetFrameworkPicoBase::ResetDeadTowerList( ){
  deadTowers.clear();
}
*/
