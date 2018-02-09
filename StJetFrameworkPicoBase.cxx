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

// STAR includes
#include "StThreeVectorF.hh"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"

// jet-framework STAR includes
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPoolManager.h"

// new includes
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  // NEW

// old file, kept for useful constants
#include "StPicoConstants.h"

// centrality
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// classes
class StJetMakerTask;

ClassImp(StJetFrameworkPicoBase)

//-----------------------------------------------------------------------------
StJetFrameworkPicoBase::StJetFrameworkPicoBase() :
  StMaker(),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),
  fCorrJetPt(kFALSE),
  fCentralityDef(4), //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection(kFALSE),
  fMinPtJet(0.0),
  fTrackBias(0.0),
  fJetRad(0.4),
  fEventZVtxMinCut(-40.0), fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fTrackPtMinCut(0.2), fTrackPtMaxCut(20.0),
  fTrackPhiMinCut(0.0), fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0), fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15), fTracknHitsRatio(0.52),
  fLeadingJet(0), fExcludeLeadingJetsFromFit(1.0), fTrackWeight(1),
  fTracksME(0x0),
  fJets(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  JetMaker(0),
  JetMakerBG(0),
  RhoMaker(0),
  grefmultCorr(0),
  refmultCorr(0),
  refmult2Corr(0),
  mOutName(""),
  fJetMakerName(""),
  fJetBGMakerName(""),
  fRhoMakerName(""),
  fRhoSparseMakerName(""),
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
  fMinPtJet(0.0),
  fTrackBias(0.0),
  fJetRad(0.4),
  fEventZVtxMinCut(-40.0), fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fTrackPtMinCut(0.2), fTrackPtMaxCut(20.0),
  fTrackPhiMinCut(0.0), fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0), fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15), fTracknHitsRatio(0.52), 
  fLeadingJet(0), fExcludeLeadingJetsFromFit(1.0), fTrackWeight(1),
  fTracksME(0x0),
  fJets(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  JetMaker(0),
  JetMakerBG(0),
  RhoMaker(0),
  grefmultCorr(0),
  refmultCorr(0),
  refmult2Corr(0),
  mOutName(""),
  fJetMakerName(""),
  fJetBGMakerName(""),
  fRhoMakerName(""),
  fRhoSparseMakerName(""),
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

  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  // initialize centrality correction
  if(fRunFlag == Run14_AuAu200) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr(); }
  if(fRunFlag == Run16_AuAu200) {
    if(fCentralityDef == kgrefmult) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr(); }
    if(fCentralityDef == kgrefmult_P16id) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id(); }
    if(fCentralityDef == kgrefmult_VpdMBnoVtx) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx(); }
    if(fCentralityDef == kgrefmult_VpdMB30) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMB30(); }
  } 

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

//----------------------------------------------------------------------------- 
//  This method is called every event.
Int_t StJetFrameworkPicoBase::Make() {
  const double pi = 1.0*TMath::Pi();
  bool printInfo = kFALSE;
  bool firstEvent = kFALSE;

  // update counter
  mEventCounter++;

  // get PicoDstMaker 
  mPicoDstMaker = (StPicoDstMaker*)GetMaker("picoDst");
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent 
  mPicoEvent = mPicoDst->event();
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField(); 

  // get vertex 3 vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();
  
  // Z-vertex cut 
  // per the Aj analysis (-40, 40)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk; //kStWarn;

  // let me know the Run #, fill, and event ID
  Int_t RunId = mPicoEvent->runId();
  Int_t fillId = mPicoEvent->fillId();
  Int_t eventId = mPicoEvent->eventId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Float_t fZDCCoincidenceRate = mPicoEvent->ZDCx();

  // get JetMaker
  JetMaker = (StJetMakerTask*)GetMaker(fJetMakerName);
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  fJets = JetMaker->GetJets();
  if(!fJets) {     
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn; 
  }

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // ========================= Trigger Info =============================== //
  // get trigger IDs from PicoEvent class and loop over them
  Int_t trId[20]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds(); 
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    // fill trigger array with trigger IDs
    trId[i] = mytriggers[i];
  }
  cout<<endl;

  // ============================ CENTRALITY ============================== //
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_refmult.txt
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_grefmult.txt
  int grefMult = mPicoEvent->grefMult();
  int refMult = mPicoEvent->refMult();
  grefmultCorr->init(RunId);
  grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate);
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  Int_t cent16 = grefmultCorr->getCentralityBin16();
  Int_t cent9 = grefmultCorr->getCentralityBin9();
  Int_t centbin = GetCentBin(cent16, 16);
  Double_t refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);

  // to limit filling unused entries in sparse, only fill for certain centrality ranges
  // ranges can be different than functional cent bin setter
  Int_t cbin = -1;
  // need to figure out centrality first in STAR: TODO
  if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
  else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
  else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
  else if (centbin>5 && centbin<10) cbin = 4; // 30-50%
  else if (centbin>9 && centbin<16) cbin = 5; // 50-80%
  else cbin = -99;

  // ============================ end of CENTRALITY ============================== //

  // reaction plane angle
  double rpAngle = GetReactionPlane();

  // get RhoMaker from event: old names "StRho_JetsBG", "OutRho", "StMaker#0"
  RhoMaker = (StRho*)GetMaker(fRhoMakerName);
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = (StRhoParameter*)RhoMaker->GetRho();
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;    
  } 
  
  // get rho/area       fRho->ls("");
  fRhoVal = fRho->GetVal();

  // event counter at end of maker
  mInputEventCounter++;

  return kStOK;
}

//________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetCentBin(Int_t cent, Int_t nBin) const
{  // Get centrality bin.
  Int_t centbin = -1;

  if(nBin == 16) {
    centbin = nBin - 1 - cent;
  }
  if(nBin == 9) {
    centbin = nBin - 1 - cent;
  }

  return centbin;
}

// this function generate a jet name based on input
TString StJetFrameworkPicoBase::GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray* partCont, TClonesArray* clusCont, TString tag)
{
  TString algoString;
  switch (jetAlgo)
  {
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
  if(dphi<-0.5*TMath::Pi()) dphi+=2.*TMath::Pi();
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
    cout<<"dPHI not in range [0, 0.5*Pi]!"<<endl;
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

//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert) {
  // declare kinematic variables
  double phi, eta, px, py, pz, pt, p, energy, charge, dca;
  int nHitsFit, nHitsMax;
  double nHitsRatio;

  // constants: assume neutral pion mass
  double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0*TMath::Pi();

  // primary track switch
  if(doUsePrimTracks) {
    if(!(trk->isPrimary())) return kFALSE; // check if primary

    // get primary track variables
    StThreeVectorF mPMomentum = trk->pMom();
    phi = mPMomentum.phi();
    eta = mPMomentum.pseudoRapidity();
    px = mPMomentum.x();
    py = mPMomentum.y();
    pt = mPMomentum.perp();
    pz = mPMomentum.z();
    p = mPMomentum.mag();
  } else {
    // get global track variables
    StThreeVectorF mgMomentum = trk->gMom(Vert, B);
    phi = mgMomentum.phi();
    eta = mgMomentum.pseudoRapidity();
    px = mgMomentum.x();
    py = mgMomentum.y();
    pt = mgMomentum.perp();
    pz = mgMomentum.z();
    p = mgMomentum.mag();
  }

  // additional calculations
  energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
  charge = trk->charge();
  dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  nHitsFit = trk->nHitsFit();
  nHitsMax = trk->nHitsMax();
  nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // do pt cut here to accommadate either type
  if(doUsePrimTracks) { // primary  track
    if(pt < fTrackPtMinCut) return kFALSE;
  } else { // global track
    if(pt < fTrackPtMinCut) return kFALSE;
  }

  // jet track acceptance cuts now - after getting 3vector - hardcoded
  if(pt > fTrackPtMaxCut) return kFALSE; // 20.0 STAR, 100.0 ALICE
  if((eta < fTrackEtaMinCut) || (eta > fTrackEtaMaxCut)) return kFALSE;
  if(phi < 0) phi+= 2*pi;
  if(phi > 2*pi) phi-= 2*pi;
  if((phi < fTrackPhiMinCut) || (phi > fTrackPhiMaxCut)) return kFALSE;
    
  // additional quality cuts for tracks
  if(dca > fTrackDCAcut) return kFALSE;
  if(nHitsFit < fTracknHitsFit) return kFALSE;
  if(nHitsRatio < fTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track
  return kTRUE;
}

//________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetReactionPlane() { 
 // get event B (magnetic) field
 Float_t Bfield = mPicoEvent->bField();

 // get vertex 3-vector and declare variables
 StThreeVectorF mVertex = mPicoEvent->primaryVertex();
 double zVtx = mVertex.z();

 //if(mVerbose)cout << "----------- In GetReactionPlane() -----------------" << endl;
 TVector2 mQ;
 double mQx = 0., mQy = 0.;
 int order = 2;
 int n,i;
 double phi, eta, pt;
 double pi = 1.0*TMath::Pi();

 // leading jet check and removal
 Float_t excludeInEta = -999;
 fLeadingJet = GetLeadingJet();
 if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from ep estimate
   if(fLeadingJet) excludeInEta = fLeadingJet->Eta();
 }

 n = mPicoDst->numberOfTracks();
 for (i=0; i<n; i++) {
   StPicoTrack* track = mPicoDst->track(i);
   if(!track) { continue; }

   // declare kinematic variables
   if(doUsePrimTracks) {
     if(!(track->isPrimary())) continue; // check if primary

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

   // do pt cut here to accommadate either type
   if(doUsePrimTracks) { // primary  track
     if(pt < fTrackPtMinCut) continue;
   } else { // global track
     if(pt < fTrackPtMinCut) continue;
   }

   // more acceptance cuts now - after getting 3vector - hardcoded for now
   if(pt > 5.0) continue;   // 100.0
   if((1.0*TMath::Abs(eta)) > 1.0) continue;
   if(phi < 0) phi+= 2*pi;
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
     if(pt > 2.0) trackweight = 2.0;
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
StJet* StJetFrameworkPicoBase::GetLeadingJet(StRhoParameter* eventRho) {
    // return pointer to the highest pt jet (before background subtraction) within acceptance
    // only rudimentary cuts are applied on this level, hence the implementation outside of
    // the framework
  if(fJets) {
    Int_t iJets(fJets->GetEntriesFast());
    Double_t pt(0);
    StJet* leadingJet(0x0);
    if(!eventRho) {
        for(Int_t i(0); i < iJets; i++) {
            StJet* jet = static_cast<StJet*>(fJets->At(i));
            if(jet->Pt() > pt) {
               leadingJet = jet;
               pt = leadingJet->Pt();
            }
        }
        return leadingJet;
    } else {
        // return leading jet after background subtraction
        Double_t rho(0);
        for(Int_t i(0); i < iJets; i++) {
            StJet* jet = static_cast<StJet*>(fJets->At(i));
            //if(!AcceptMyJet(jet)) continue;
            if((jet->Pt() - jet->Area()*fRhoVal) > pt) {
               leadingJet = jet;
               pt = (leadingJet->Pt()-jet->Area()*fRhoVal);
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

  Bool_t doAnalysis;
  doAnalysis = kFALSE; // set false by default, to make sure user chooses an available bin

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
