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
#include "TVector3.h"
#include <sstream>
#include <fstream>

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"

// jet-framework STAR includes
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPlaneMaker.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu

// old file, kept for useful constants
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

//#include "StEmcUtil/projection/StEmcPosition.h" // old
#include "StEmcPosition2.h"

ClassImp(StJetFrameworkPicoBase)

//
//________________________________________________________________________________________________
StJetFrameworkPicoBase::StJetFrameworkPicoBase() :
  StMaker(),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),
  doppAnalysis(kFALSE),
  doJetShapeAnalysis(kFALSE),
  fCorrJetPt(kFALSE),
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection(kFALSE),
  doUseBBCCoincidenceRate(kFALSE), // kFALSE = use ZDC
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
//  mVertex(0x0),
  zVtx(0.0),
  fMaxEventTrackPt(30.0),
  fMaxEventTowerE(30.0),
  doRejectBadRuns(kFALSE),
  fBadRunListVers(999),
  fBadTowerListVers(0),
  fJetType(0),
  fMinPtJet(0.0),
  fTrackBias(0.2),
  fTowerBias(0.2),
  fJetRad(0.4),
  fEventZVtxMinCut(-40.0), fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fTrackPtMinCut(0.2), fTrackPtMaxCut(30.0),
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
  fJets1(0x0),
  fJets2(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  JetMaker(0x0),
  JetMakerBG(0x0),
  JetMaker1(0x0),
  JetMaker2(0x0),
  RhoMaker(0x0),
  RhoMaker1(0x0),
  RhoMaker2(0x0),
  EventPlaneMaker(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0),
  refmultCorr(0x0),
  refmult2Corr(0x0),
  mOutName(""),
  mOutNameEP(""),
  mOutNameQA(""),
  fJetMakerName(""),
  fJetBGMakerName(""),
  fRhoMakerName(""),
  fRhoSparseMakerName(""),
  fEventPlaneMakerName(""),
  fRho(0x0),
  fRho1(0x0),
  fRho2(0x0),
  fRhoVal(0),
  fAddToHistogramsName(""),
  mEventCounter(0),
  mAllPVEventCounter(0),
  mInputEventCounter(0)
{

}
//
//______________________________________________________________________________________________
StJetFrameworkPicoBase::StJetFrameworkPicoBase(const char* name) :
  StMaker(name),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),
  doppAnalysis(kFALSE),
  doJetShapeAnalysis(kFALSE),
  fCorrJetPt(kFALSE),
  fCentralityDef(4), //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection(kFALSE),
  doUseBBCCoincidenceRate(kTRUE), // kFALSE = use ZDC
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
//  mVertex(0x0),
  zVtx(0.0),
  fMaxEventTrackPt(30.0),
  fMaxEventTowerE(30.0),
  doRejectBadRuns(kFALSE),
  fBadRunListVers(999),
  fBadTowerListVers(0),
  fJetType(0),
  fMinPtJet(0.0),
  fTrackBias(0.2),
  fTowerBias(0.2),
  fJetRad(0.4),
  fEventZVtxMinCut(-40.0), fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fTrackPtMinCut(0.2), fTrackPtMaxCut(30.0),
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
  fJets1(0x0),
  fJets2(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  JetMaker(0x0),
  JetMakerBG(0x0),
  JetMaker1(0x0),
  JetMaker2(0x0),
  RhoMaker(0x0),
  RhoMaker1(0x0),
  RhoMaker2(0x0),
  EventPlaneMaker(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0),
  refmultCorr(0x0),
  refmult2Corr(0x0),
  mOutName(""),
  mOutNameEP(""),
  mOutNameQA(""),
  fJetMakerName(""),
  fJetBGMakerName(""),
  fRhoMakerName(""),
  fRhoSparseMakerName(""),
  fEventPlaneMakerName(""),
  fRho(0x0),
  fRho1(0x0),
  fRho2(0x0),
  fRhoVal(0),
  fAddToHistogramsName(""),
  mEventCounter(0),
  mAllPVEventCounter(0),
  mInputEventCounter(0)
{

}
//
//___________________________________________________________________________________
StJetFrameworkPicoBase::~StJetFrameworkPicoBase()
{ /*  */
  // destructor
}
//
//___________________________________________________________________________________
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

    // TODO - not sure, may want to add other runs? this base class may not be best place for reading in lists

    default :
      AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_BadTowers.txt");
      AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_DeadTowers.txt");
  }
*/

/*
  // Add bad run lists
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp (200 GeV)
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_w_missing_HT)  AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2012_BadRuns_P12id_w_missing_HT.txt");
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_wo_missing_HT) AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2012_BadRuns_P12id_wo_missing_HT.txt");
        break;
  
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu (200 GeV)
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_w_missing_HT)  AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih_w_missing_HT.txt");
        if(fBadRunListVers == StJetFrameworkPicoBase::fBadRuns_wo_missing_HT) AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih_wo_missing_HT.txt");
        break; 
  
    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu (200 GeV)
        AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2016_BadRuns_P16ij.txt");
        break; 
  
    default :
      AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Empty_BadRuns.txt");
  }
*/

  ///fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  return kStOK;
}
//
//
//_________________________________________________________________________________________
Int_t StJetFrameworkPicoBase::Finish() { 
  return kStOK;
}
//
// OLD user code says: //  Called every event after Make(). 
//_________________________________________________________________________________________
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
//
// get centrality bin
// - input is raw cent bin & nBins to re-order from
//________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetCentBin(Int_t cent, Int_t nBin) const
{
  Int_t centbin = -1;

  if(nBin == 16) { centbin = nBin - 1 - cent; }
  if(nBin == 9)  { centbin = nBin - 1 - cent; }

  return centbin;
}
//
// get centrality bin in 4 different bins
// - input is the scaled centrality
//________________________________________________________________________
Int_t StJetFrameworkPicoBase::Get4CentBin(Double_t scaledCent) const
{
  // initialize centrality bin
  int centbin = -99;

  // get centrality bin number
  if(scaledCent >= 0 && scaledCent <  10.0)  { centbin = 0; }
  else if(scaledCent <  20.0)                { centbin = 1; }
  else if(scaledCent <  50.0)                { centbin = 2; }
  else if(scaledCent <= 80.0)                { centbin = 3; }

  return centbin;
}
//
// function to get annuli bin
//___________________________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetAnnuliBin(Double_t deltaR) const
{
  // initialize annuli bin
  int annuliBin = -99;

  // get annuli bin number
  if(deltaR >= 0.00 && deltaR <= 0.05)     { annuliBin = 0; }
  else if(deltaR > 0.05 && deltaR <= 0.10) { annuliBin = 1; }
  else if(deltaR > 0.10 && deltaR <= 0.15) { annuliBin = 2; }
  else if(deltaR > 0.15 && deltaR <= 0.20) { annuliBin = 3; } 
  else if(deltaR > 0.20 && deltaR <= 0.25) { annuliBin = 4; }
  else if(deltaR > 0.25 && deltaR <= 0.30) { annuliBin = 5; }
  else if(deltaR > 0.30 && deltaR <= 0.35) { annuliBin = 6; }
  else if(deltaR > 0.35 && deltaR <= 0.40) { annuliBin = 7; }
  else if(deltaR > 0.40 && deltaR <= 0.45) { annuliBin = 8; }
  else if(deltaR > 0.45 && deltaR <= 0.50) { annuliBin = 9; }

  return annuliBin;
}
//
// function to jet pt bin
//___________________________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetJetPtBin(Double_t jetpt) const
{
  // initialize jet pt bin
  int jetPtBin = -99;

  // get jet pt bin number
  if(jetpt >= 10.0 && jetpt < 15.0)      { jetPtBin = 0; } 
  else if(jetpt >= 15.0 && jetpt < 20.0) { jetPtBin = 1; }
  else if(jetpt >= 20.0 && jetpt < 40.0) { jetPtBin = 2; } 
  else if(jetpt >= 40.0 && jetpt < 60.0) { jetPtBin = 3; }

  return jetPtBin;
}
//
// function to jet event plane bin
//___________________________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetJetEPBin(Double_t dEP) const
{
  // constants
  double pi = 1.0*TMath::Pi();

  // initialize jet event plane bin
  int jetEPBin = -99;

  // get jet event plane bin number
  if(dEP >= 0.0*pi/6.0 && dEP <= 1.0*pi/6.0)     { jetEPBin = 0; }
  else if(dEP > 1.0*pi/6.0 && dEP <= 2.0*pi/6.0) { jetEPBin = 1; }
  else if(dEP > 2.0*pi/6.0 && dEP <= 3.0*pi/6.0) { jetEPBin = 2; }

  return jetEPBin;
}
//
// this function generates a jet name based on input
//___________________________________________________________________________________________
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
//
// Function: to calculate relative Phi
//________________________________________________________________________
Double_t StJetFrameworkPicoBase::RelativePhi(Double_t jphi, Double_t tphi) const
{
  // calculate relative phi
  double dphi = jphi - tphi;

  // set dphi to operate on adjusted scale
  if(dphi < -0.5*TMath::Pi()) dphi += 2.*TMath::Pi();
  if(dphi >  1.5*TMath::Pi()) dphi -= 2.*TMath::Pi();

  // test check
  if( dphi < -0.5*TMath::Pi() || dphi > 1.5*TMath::Pi() )
    Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName());

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}
//
// Function: calculate angle between jet and EP in the 1st quadrant (0,Pi/2) 
//_________________________________________________________________________
Double_t StJetFrameworkPicoBase::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{
  Double_t pi = 1.0*TMath::Pi();
  Double_t dphi = 1.0*TMath::Abs(EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi < -pi ){
    dphi = dphi + pi;
  } // this assumes we are doing full jets currently 
 
  if(dphi > 1.5*pi) dphi -= 2.0*pi;
  if((dphi > 1.0*pi) && (dphi < 1.5*pi)) dphi -= 1.0*pi;
  if((dphi > 0.5*pi) && (dphi < 1.0*pi)) dphi -= 1.0*pi;
  dphi = 1.0*TMath::Abs(dphi);

  // test check
  if( dphi < 0 || dphi > 0.5*TMath::Pi() ) {
    //Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());
    cout<<"dEP not in range [0, 0.5*Pi]!"<<" jetAng: "<<jetAng<<"  EPang: "<<EPAng<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}
//
/*
//
// Function: clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
//_________________________________________________
TClonesArray* StJetFrameworkPicoBase::CloneAndReduceTrackList(TClonesArray* tracksME)
{
  TClonesArray *tracksClone = new TClonesArray("StPicoTrack");
  //tracksClone->SetName("tracksClone");
  //tracksClone->SetOwner(kTRUE);

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3 vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  Int_t nMixTracks = mPicoDst->numberOfTracks();
  Int_t iterTrk = 0;
  Double_t phi, eta, px, py, pt, pz, p, charge;
  const double pi = 1.0*TMath::Pi();
  for (Int_t i = 0; i < nMixTracks; i++) { 
    // get track pointer
    StPicoTrack *trk = mPicoDst->track(i);
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
// - can call from a class that inherits from this base class and that classes' set 
// members will be used here
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert) {
  // constants: assume neutral pion mass
  //double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0*TMath::Pi();

  // primary track switch: get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    if(!(trk->isPrimary())) return kFALSE; // check if primary
    // get primary track vector
    mTrkMom = trk->pMom();
  } else {
    // get global track vector
    mTrkMom = trk->gMom(Vert, B);
  }

  // track variables
  double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  //double px = mTrkMom.x();
  //double py = mTrkMom.y();
  //double pz = mTrkMom.z();
  //double p = mTrkMom.Mag();
  //double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
  //short charge = trk->charge();
  ////double dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // track acceptance cuts now - after getting 3vector - hardcoded
  if(pt < fTrackPtMinCut) return kFALSE;
  if(pt > fTrackPtMaxCut) return kFALSE; // 20.0 STAR, (increased to 30.0) 100.0 ALICE
  if((eta < fTrackEtaMinCut) || (eta > fTrackEtaMaxCut)) return kFALSE;
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  if((phi < fTrackPhiMinCut) || (phi > fTrackPhiMaxCut)) return kFALSE;
    
  // additional quality cuts for tracks
  if(dca > fTrackDCAcut)            return kFALSE;
  if(nHitsFit < fTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fTracknHitsRatio) return kFALSE;

  //cout<<"doUsePrimTracks: "<<doUsePrimTracks<<"  fTrackPtMinCut: "<<fTrackPtMinCut<<"  fTrackPtMaxCut: "<<fTrackPtMaxCut<<"  fTrackEtaMinCut: "<<fTrackEtaMinCut<<"  fTrackEtaMaxCut: "<<fTrackEtaMaxCut<<"  fTrackPhiMinCut: "<<fTrackPhiMinCut<<"  fTrackPhiMaxCut: "<<fTrackPhiMaxCut<<"  fTrackDCAcut: "<<fTrackDCAcut<<"  fTracknHitsFit: "<<fTracknHitsFit<<"  fTracksnHitsRatio: "<<fTracknHitsRatio<<endl;

  // passed all above cuts - keep track
  return kTRUE;
}
/*
//
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::AcceptTower(StPicoBTowHit *tower, TVector3 Vertex, Int_t towerID) {
  // get EMCal position - FIXME if using
  StEmcPosition2 *mPosition = new StEmcPosition2();

  // constants:
  double pi = 1.0*TMath::Pi();

  // tower ID - passed into function: make sure some of these aren't still in event array
  if(towerID < 0) return kFALSE; 

  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  TVector3 towerPosition = mPosition->getPosFromVertex(mVertex, towerID);

  // free up memory
  delete mPosition;

  double phi = towerPosition.Phi();
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  double eta = towerPosition.PseudoRapidity();
  int towerADC = tower->adc();
  double towerEunCorr = tower->energy();  // uncorrected energy

  // check for bad (and dead) towers
  bool TowerOK = IsTowerOK(towerID);      // kTRUE means GOOD
  bool TowerDead = IsTowerDead(towerID);  // kTRUE means BAD
  if(!TowerOK)  { return kFALSE; }
  if(TowerDead) { return kFALSE; }

  // tower acceptance cuts now - after getting 3vector - hardcoded
  if((eta < fTowerEtaMinCut) || (eta > fTowerEtaMaxCut)) return kFALSE;
  if((phi < fTowerPhiMinCut) || (phi > fTowerPhiMaxCut)) return kFALSE;

  // passed all above cuts - keep tower and fill input vector to fastjet
  return kTRUE;
}
*/
//
// general reaction plane function
//________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetReactionPlane() { 
  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3-vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();

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
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      if(!(track->isPrimary())) return kFALSE; // check if primary
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

    // more acceptance cuts now - after getting 3vector - hardcoded for now
    if(pt < fTrackPtMinCut) continue; 
    if(pt > 5.0) continue;   // 100.0
    if((1.0*TMath::Abs(eta)) > 1.0) continue;
    if(phi < 0.0)    phi+= 2.0*pi;
    if(phi > 2.0*pi) phi-= 2.0*pi;
    if((phi < 0) || (phi > 2.0*pi)) continue;

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
      if(pt >  2.0) trackweight = 2.0;
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
//
// get dijet asymmetry
// _____________________________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetDiJetAj(StJet *jet1, StJet *jet2, StRhoParameter *eventRho, Bool_t doCorrJetPt) {
  // returns dijet asymmetry Aj of 2 jets
  // - meant for leading/subleading jets
  // - jet1 should have pt > jet2 pt

  // get jet1 and jet2 pt's
  double jetPt1, jetPt2;
  if(doCorrJetPt) {
    jetPt1 = jet1->Pt() - jet1->Area()*eventRho->GetVal();
    jetPt2 = jet2->Pt() - jet2->Area()*eventRho->GetVal();
  } else {
    jetPt1 = jet1->Pt();
    jetPt2 = jet2->Pt();
  }

  // calculate dijet imbalance Aj
  double Aj = (jetPt1 - jetPt2) / (jetPt1 + jetPt2);

  // return absolute value in case user input wrong
  return 1.0*TMath::Abs(Aj);
}
//
// get leading jet pointer
// _____________________________________________________________________________________________
StJet* StJetFrameworkPicoBase::GetLeadingJet(TString fJetMakerNametemp, StRhoParameter *eventRho) {
  // return pointer to the highest pt jet (before/after background subtraction) within acceptance
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
    // get number of jets, initialize pt and leading jet pointer
    Int_t iJets(fJets->GetEntriesFast());
    Double_t pt(0);
    StJet *leadingJet(0x0);

    if(!eventRho) { // no rho parameter provided
      // loop over jets
      for(Int_t i(0); i < iJets; i++) {
        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJets->At(i));
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

      // loop over corrected jets
      for(Int_t i(0); i < iJets; i++) {
        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJets->At(i));
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
//
// get subledaing jet pointer
// _____________________________________________________________________________________________
StJet* StJetFrameworkPicoBase::GetSubLeadingJet(TString fJetMakerNametemp, StRhoParameter *eventRho) {
  // return pointer to the second highest pt jet (before/after background subtraction) within acceptance
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
    StJet *leadingJet(0x0);
    StJet *subleadingJet(0x0);

    if(!eventRho) { // no rho parameter provided
      // loop over jets
      for(Int_t i(0); i < iJets; i++) {
        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJets->At(i));
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
      // return subleading jet after background subtraction
      // BUG: Fixed October 26, 2018: returned leadingjet, and didn't subtract bg from subleading
      //Double_t rho(0);
      double fRhoValtemp = eventRho->GetVal(); // test

      // loop over corrected jets
      for(Int_t i(0); i < iJets; i++) {
        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJets->At(i));
        // leading jet
        if((jet->Pt() - jet->Area()*fRhoValtemp) > pt) {
          leadingJet = jet;
          pt = (leadingJet->Pt()-jet->Area()*fRhoValtemp);
        }

        // subleading jet
        if((jet->Pt() - jet->Area()*fRhoValtemp > pt2) && (jet->Pt() - jet->Area()*fRhoValtemp < pt)) {
          subleadingJet = jet;
          pt2 = subleadingJet->Pt() - jet->Area()*fRhoValtemp;
        }
      }
      return subleadingJet;
    }
  }

  return 0x0;
}
//
// Function: Event counter
//__________________________________________________________________________________________
Int_t StJetFrameworkPicoBase::EventCounter() {
  mEventCounter++;

  return mEventCounter;
}
//
// Function: select analysis centrality bin, if performing class wide cut 
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
  //    if(!SelectAnalysisCentralityBin(centbin, StJetFrameworkPicoBase::kCent3050)) return kStOK; (or StOk to suppress warnings)
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
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::DoComparison(int myarr[], int elems) {
  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) { LOG_WARN << " No PicoDstMaker! Skip! " << endm;  return kStWarn; }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) { LOG_WARN << " No PicoDst! Skip! " << endm; return kStWarn; }

  // create pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) { LOG_WARN << " No PicoEvent! Skip! " << endm; return kStWarn; }

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
//
// Function: check if event fired min-bias (MB) trigger
//_________________________________________________________________________
Bool_t StJetFrameworkPicoBase::CheckForMB(int RunFlag, int type) {
  // Run11 triggers: pp
  int arrMB_Run11[] = {320000, 320001, 320011, 320021, 330021};

  // Run12 (200 GeV pp) triggers: 1) VPDMB, 370011-main
  int arrMB_Run12[] = {370011};
  int arrMB_Run12extra[] = {370001, 370011};

  // Run13 triggers: pp
  int arrMB_Run13[] = {430001, 430011, 430021, 430031};

  // Run14 triggers: 200 GeV AuAu
  int arrMB_Run14[] = {450014};
  int arrMB30_Run14[] = {450010, 450020};
  int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
  // additional 30: 450201, 450202, 450211, 450212
  // 1: VPDMB-5   Run:            15075055 - 15076099
  // 1: VPDMB-5-p-nobsmd-hlt Run: 15081020 - 15090048
  // 4: VPDMB-5-p-nobsmd-hlt Run: 15090049 - 15167007

  // Run16 triggers: 200 GeV AuAu
  int arrMB_Run16[] = {520021};
  int arrMB5_Run16[] = {520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
  int arrMB10_Run16[] = {520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};

  // Run17 triggers: 510 GeV pp
  int arrMB30_Run17[] = {570001, 590001};
  int arrMB100_Run17[] = {590002};
  int arrMBnovtx_Run17[] = {570004};

  // run flag selection to check for MB firing
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run11_pp500 : // Run11 pp
        switch(type) {
          case StJetFrameworkPicoBase::kVPDMB :
              if((DoComparison(arrMB_Run11, sizeof(arrMB_Run11)/sizeof(*arrMB_Run11)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run11, sizeof(arrMB_Run11)/sizeof(*arrMB_Run11)))) { return kTRUE; }
        }
        break;

    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp (200 GeV)
        switch(type) {
          case StJetFrameworkPicoBase::kRun12main :  // update if needed
              if((DoComparison(arrMB_Run12, sizeof(arrMB_Run12)/sizeof(*arrMB_Run12)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMB :
              if((DoComparison(arrMB_Run12extra, sizeof(arrMB_Run12extra)/sizeof(*arrMB_Run12extra)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run12, sizeof(arrMB_Run12)/sizeof(*arrMB_Run12)))) { return kTRUE; }
        }
        break;

    case StJetFrameworkPicoBase::Run13_pp510 : // Run13 pp
        switch(type) {
          case StJetFrameworkPicoBase::kVPDMB :
              if((DoComparison(arrMB_Run13, sizeof(arrMB_Run13)/sizeof(*arrMB_Run13)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run13, sizeof(arrMB_Run13)/sizeof(*arrMB_Run13)))) { return kTRUE; }
        }
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu (200 GeV)
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
          default : // kRun14Main or kVPDMB
              if((DoComparison(arrMB_Run14, sizeof(arrMB_Run14)/sizeof(*arrMB_Run14)))) { return kTRUE; }
        }
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu (200 GeV)
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
        break;

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17 pp (510 GeV)
        switch(type) {
          case StJetFrameworkPicoBase::kVPDMB30 :
              if((DoComparison(arrMB30_Run17, sizeof(arrMB30_Run17)/sizeof(*arrMB30_Run17)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMB100 :
              if((DoComparison(arrMB100_Run17, sizeof(arrMB100_Run17)/sizeof(*arrMB100_Run17)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMBnovtx :
              if((DoComparison(arrMBnovtx_Run17, sizeof(arrMBnovtx_Run17)/sizeof(*arrMBnovtx_Run17)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB30_Run17, sizeof(arrMB30_Run17)/sizeof(*arrMB30_Run17)))) { return kTRUE; }
        }
        break;

  } // RunFlag switch

  // return status
  return kFALSE;
} // MB function
//
// check to see if the event was EMC triggered for High Towers
//____________________________________________________________________________
Bool_t StJetFrameworkPicoBase::CheckForHT(int RunFlag, int type) {
  // Run12 (200 GeV pp) triggers:
  int arrHT1_Run12[] = {370511, 370546};
  int arrHT2_Run12[] = {370521, 370522, 370531, 370980};
  int arrHT3_Run12[] = {380206, 380216}; // NO HT3 triggers in this dataset

  // Run14 triggers: 200 GeV AuAu
  int arrHT1_Run14[] = {450201, 450211, 460201};
  int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
  int arrHT3_Run14[] = {450203, 450213, 460203};

  // Run16 triggers: 200 GeV AuAu
  int arrHT1_Run16[] = {520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
  int arrHT2_Run16[] = {530202, 540203};
  int arrHT3_Run16[] = {520203, 530213};

  // Run17 triggers: (HT1 and HT2 not exclusive) 510 GeV pp
  int arrHT1_Run17[] = {570204, 570214};
  int arrHT2_Run17[] = {570205, 570215};
  int arrHT3_Run17[] = {570201, 590201};

  // run flag selection to check for MB firing
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp (200 GeV)
        switch(type) {
          case StJetFrameworkPicoBase::kIsHT1 :
              if((DoComparison(arrHT1_Run12, sizeof(arrHT1_Run12)/sizeof(*arrHT1_Run12)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT2 :
              if((DoComparison(arrHT2_Run12, sizeof(arrHT2_Run12)/sizeof(*arrHT2_Run12)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT3 :
              if((DoComparison(arrHT3_Run12, sizeof(arrHT3_Run12)/sizeof(*arrHT3_Run12)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrHT2_Run12, sizeof(arrHT2_Run12)/sizeof(*arrHT2_Run12)))) { return kTRUE; }
        }
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu (200 GeV)
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
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu (200 GeV)
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
        break;

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17 pp (510 GeV)
        switch(type) {
          case StJetFrameworkPicoBase::kIsHT1 :
              if((DoComparison(arrHT1_Run17, sizeof(arrHT1_Run17)/sizeof(*arrHT1_Run17)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT2 :
              if((DoComparison(arrHT2_Run17, sizeof(arrHT2_Run17)/sizeof(*arrHT2_Run17)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT3 :
              if((DoComparison(arrHT3_Run16, sizeof(arrHT3_Run17)/sizeof(*arrHT3_Run17)))) { return kTRUE; }
              break;
          default : // HT3
              if((DoComparison(arrHT3_Run17, sizeof(arrHT3_Run17)/sizeof(*arrHT3_Run17)))) { return kTRUE; }
        }
        break; 

  } // RunFlag switch

  return kFALSE;
}
//
// Function: calculate momentum of a tower
//________________________________________________________________________________________________________
Bool_t StJetFrameworkPicoBase::GetMomentum(TVector3 &mom, const StPicoBTowHit* tower, Double_t mass, StPicoEvent *PicoEvent, Int_t towerID) const {
  // initialize Emc position objects - FIXME if using
  StEmcPosition2 *Position = new StEmcPosition2();

  // vertex components
  TVector3 fVertex = PicoEvent->primaryVertex();
//  double xVtx = fVertex.x();
//  double yVtx = fVertex.y();
//  double zVtx = fVertex.z();

  // get mass, E, P, ID
  if(mass < 0) mass = 0;
  Double_t energy = tower->energy();
  Double_t p = TMath::Sqrt(energy*energy - mass*mass);

  // tower ID - passed into function
  // get tower position
  TVector3 towerPosition = Position->getPosFromVertex(fVertex, towerID);
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
    mom.SetX( p*posX/r );
    mom.SetY( p*posY/r );
    mom.SetZ( p*posZ/r );
    // energy) ;
  } else { return kFALSE; }

  // free up memory
  delete Position;

  return kTRUE;
}

//
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
//
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
//
//____________________________________________________________________________
void StJetFrameworkPicoBase::ResetBadTowerList( ){
  badTowers.clear();
}
//
// Add bad towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
//____________________________________________________________________________
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
//
// Add dead towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
//_________________________________________________________________________
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
//
//____________________________________________________________________________
void StJetFrameworkPicoBase::ResetDeadTowerList( ){
  deadTowers.clear();
}

//
// function to convert 5% centrality to 10% bins
// must already be 'properly calculated' i.e. increasing bin# -> increasing centrality
//__________________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetCentBin10(Int_t cbin) const {
  int cbin10;
  if(cbin== 0 || cbin== 1) cbin10 = 0; //  0-10%
  if(cbin== 2 || cbin== 3) cbin10 = 1; // 10-20%
  if(cbin== 4 || cbin== 5) cbin10 = 2; // 20-30%
  if(cbin== 6 || cbin== 7) cbin10 = 3; // 30-40%
  if(cbin== 8 || cbin== 9) cbin10 = 4; // 40-50%
  if(cbin==10 || cbin==11) cbin10 = 5; // 50-60%
  if(cbin==12 || cbin==13) cbin10 = 6; // 60-70%
  if(cbin==14 || cbin==15) cbin10 = 7; // 70-80%

  return cbin10;
}
//
// Returns pt of hardest track in the event
//__________________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetMaxTrackPt()
{
  // get # of tracks
  int nTrack = mPicoDst->numberOfTracks();
  double fMaxTrackPt = -99;

  // loop over all tracks
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();

    // get max track
    if(pt > fMaxTrackPt) { fMaxTrackPt = pt; }
  }

  return fMaxTrackPt;
}
//
// Function: Returns E of most energetic tower in the event
// TODO this function needs to be re-thought, as select 'bad towers' have static Energy reading which is meaningless
//      and sometimes over the requested threshold, thus excluding event.  Set default value to 1000 for now.. July 11, 2019
//_________________________________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetMaxTowerE()
{
  // get # of towers
  int nTowers = mPicoDst->numberOfBTowHits();
  double fMaxTowerE = -99;

  // loop over all towers
  for(int i = 0; i < nTowers; i++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(i));
    if(!tower) continue;

    // tower position - from vertex and ID: shouldn't need additional eta correction
    //TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    //double towerPhi = towerPosition.Phi();
    //if(towerPhi < 0.0)    towerPhi += 2.0*pi;
    //if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;
    //double towerEta = towerPosition.PseudoRapidity();
    double towerEuncorr = tower->energy();               // uncorrected energy
    //double towEt = towerE / (1.0*TMath::CosH(towerEta)); // should this be used instead?

    // get max tower
    if(towerEuncorr > fMaxTowerE) { fMaxTowerE = towerEuncorr; }
  }

  return fMaxTowerE;
}
//
// Returns correction for tracking efficiency
//
//Double_t StJetFrameworkPicoBase::ApplyTrackingEffpp(StPicoTrack *trk)
//Double_t StJetFrameworkPicoBase::ApplyTrackingEffAuAu(StPicoTrack *trk)
//____________________________________________________________________________________________
Double_t StJetFrameworkPicoBase::ApplyTrackingEff(StPicoTrack *trk, Bool_t applyEff)
{
  // initialize effieciency
  double trEff = 1.0;
  if(!applyEff) return trEff;

  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    // get primary track vector
    mTrkMom = trk->pMom();
  } else {
    // get global track vector
    mTrkMom = trk->gMom(mVertex, Bfield);
  }

  // track variables
  double tpt = mTrkMom.Perp();
  double teta = mTrkMom.PseudoRapidity();
  double tphi = mTrkMom.Phi();

  // track efficiency has pt, eta, and centrality dependence 
  //
  // RunFlag switch
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 :   // Run12 pp
        break;
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        break;
    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        break; 
  } // RunFlag switch

  // return the single track reconstruction efficiency for the corresponding dataset
  return trEff;
}
//
// Function: check for and reject bad runs
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::RejectRun(int RunFlag, int nRun) const {
  bool fRejectEvent = kFALSE;

  int fBadRuns[] = {0,1};
  //int *fBadRuns;

  // RunFlag switch
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp
        //fBadRuns = new int [3]; // FIXME 
        break;
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        //fBadRuns = new int{0, 1, 4, 6};
        break;
    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        //fBadRuns = new int{1, 2, 3}; // FIXME
        break;

  } // RunFlag switch

  size_t nBadRuns = sizeof(fBadRuns)/sizeof(fBadRuns[0]);

  // loop over bad runs
  for(int i = 0; i < nBadRuns; i++) {
    // compare run to bad runs
    if(nRun == fBadRuns[i]) { 
      fRejectEvent = kTRUE;
      break;
    }
  }

  return fRejectEvent;
}
//
// Trigger QA histogram, label bins 
// check and fill a Event Selection QA histogram for different trigger selections after cuts
//_____________________________________________________________________________
TH1* StJetFrameworkPicoBase::FillEventTriggerQA(TH1* h) {
  // Run12 pp 200 GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    // Run12 (200 GeV pp) triggers:
    int arrHT1[] = {370511, 370546};
    int arrHT2[] = {370521, 370522, 370531, 370980};
    //int arrHT3[] = {380206, 380216}; // NO HT3 triggered events
    int arrMB[] = {370001, 370011, 370983};

    // fill for kAny
    h->Fill(1);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { h->Fill(2); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { h->Fill(3); } // HT2
    //if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { h->Fill(4); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { h->Fill(10); } // VPDMB

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(2, "BHT1");
    h->GetXaxis()->SetBinLabel(3, "BHT2");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, ""); //"VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, ""); //"Central-5");
    h->GetXaxis()->SetBinLabel(8, ""); //"Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB");
  }

  // Run14 AuAu 200 GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    int arrBHT1[] = {450201, 450211, 460201};
    int arrBHT2[] = {450202, 450212, 460202, 460212};
    int arrBHT3[] = {460203, 450213, 460203};
    int arrMB[] = {450014};
    int arrMB5[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
    int arrMB30[] = {450010, 450020};
    int arrCentral5[] = {450010, 450020};
    int arrCentral[] = {460101, 460111};

    // fill for kAny
    h->Fill(1);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1))) { h->Fill(2); } // HT1
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2))) { h->Fill(3); } // HT2
    if(DoComparison(arrBHT3, sizeof(arrBHT3)/sizeof(*arrBHT3))) { h->Fill(4); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { h->Fill(5); } // MB 
    if(DoComparison(arrCentral5, sizeof(arrCentral5)/sizeof(*arrCentral5))) { h->Fill(7); } // Central-5
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { h->Fill(8); }    // Central & Central-mon
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { h->Fill(10); }    // VPDMB-5 
    if(DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { h->Fill(11); } // VPDMB-30

    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2)) && DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { h->Fill(13); } // HT2 && MB
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2)) && DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { h->Fill(14); } // HT2 && MB30
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1)) && DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { h->Fill(15); } // HT1 && MB
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1)) && DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { h->Fill(16); } // HT1 && MB30

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(2, "BHT1*VPDMB-30");
    h->GetXaxis()->SetBinLabel(3, "BHT2*VPDMB-30");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, "Central-5");
    h->GetXaxis()->SetBinLabel(8, "Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-30");
    h->GetXaxis()->SetBinLabel(13, "HT2*VPDMB30 && MB");
    h->GetXaxis()->SetBinLabel(14, "HT2*VPDMB30 && MB30");
    h->GetXaxis()->SetBinLabel(15, "HT1*VPDMB30 && MB");
    h->GetXaxis()->SetBinLabel(16, "HT1*VPDMB30 && MB30");
  }

  // Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    // hard-coded trigger Ids for run16
    //int arrBHT0[] = {520606, 520616, 520626, 520636, 520646, 520656};
    int arrBHT1[] = {520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
    int arrBHT2[] = {530202, 540203};
    int arrBHT3[] = {520203, 530213};
    int arrMB[] = {520021};
    int arrMB5[] = {520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
    int arrMB10[] = {520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};
    int arrCentral[] = {520101, 520111, 520121, 520131, 520141, 520103, 520113, 520123};

    // fill for kAny
    h->Fill(1);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1))) { h->Fill(2); } // HT1
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2))) { h->Fill(3); } // HT2
    if(DoComparison(arrBHT3, sizeof(arrBHT3)/sizeof(*arrBHT3))) { h->Fill(4); } // HT3
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { h->Fill(5); }  // MB
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { h->Fill(7); } // Central-5 & Central-novtx
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { h->Fill(10); }    // VPDMB-5 
    if(DoComparison(arrMB10, sizeof(arrMB10)/sizeof(*arrMB10))) { h->Fill(11); } // VPDMB-10

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(2, "BHT1");
    h->GetXaxis()->SetBinLabel(3, "BHT2");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-p-sst");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, "Central");
    h->GetXaxis()->SetBinLabel(8, "");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-10");
  }

  // set general label
  h->GetXaxis()->SetBinLabel(1, "Any trigger");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}
//
// this function checks for the bin number of the run from a runlist header 
// in order to apply various corrections and fill run-dependent histograms
// 1287 - Liang
// _________________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetRunNo(int runid){
  // Run12 pp (200 GeV)
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    for(int i = 0; i < 857; i++) {
      if(Run12pp_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run14 AuAu
  // Run14AuAu_IdNo: SL17id
  // 1654 for Run14 AuAu, new picoDst production is 830
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    for(int i = 0; i < 830; i++) {
      if(Run14AuAu_P18ih_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run16 AuAu
  // 1359 for Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    for(int i = 0; i < 1359; i++){
      if(Run16AuAu_IdNo[i] == runid) {
        return i;
      }
    }
  }

  cout<<" *********** RunID not matched with list ************!!!! "<<endl;
  return -999;
}
/*
//
// additional jet quality / acceptance cuts - done in jet finder
//________________________________________________________________________
Bool_t StJetFrameworkPicoBase::AcceptJet(StJet *jet) { // for jets
  // applies all jet cuts except pt
  if ((jet->Phi() < fPhimin) || (jet->Phi() > fPhimax)) return kFALSE;
  if ((jet->Eta() < fEtamin) || (jet->Eta() > fEtamax)) return kFALSE;
  if (jet->Area() < fAreacut) return kFALSE;

  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area() == 0) return kFALSE;

  // exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt() > 30) return kFALSE;

  // jet passed all above cuts
  return kTRUE;
}
*/
//
// function to calcuate delta R between a jet centroid and a track
//___________________________________________________________________________________________
Double_t StJetFrameworkPicoBase::GetDeltaR(StJet *jet, StPicoTrack *trk) {
  // constants
  double deltaR = -99.;
  double pi = 1.0*TMath::Pi();

  // get track momentum vector 
  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    if(!(trk->isPrimary())) return -99.;
    // get primary track vector
    mTrkMom = trk->pMom();
  } else {
    // get global track vector
    mTrkMom = trk->gMom(mVertex, Bfield);
  }

  // track variables
  double tphi = mTrkMom.Phi();
  if(tphi > 2.0*pi) tphi -= 2.0*pi;
  if(tphi < 0.0   ) tphi += 2.0*pi;
  double teta = mTrkMom.PseudoRapidity();

  // jet variables
  double jphi = jet->Phi();
  double jeta = jet->Eta();

  // calculate radial distance
  double deltaEta = 1.0*TMath::Abs(jeta - teta);
  double deltaPhi = 1.0*TMath::Abs(jphi - tphi);
  if(deltaPhi > 1.0*pi) deltaPhi = 2.0*pi - deltaPhi;
  deltaR = 1.0*TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

  return deltaR;
}
//
//
// this is code from Liang to get the Vz region for event plane corrections
// __________________________________________________________________________________
Int_t StJetFrameworkPicoBase::GetVzRegion(double Vz) // 0-14, 15          0-19, 20
{
  /////////if(fabs(Vz >= 30.)) return 999;  // THIS CUT IS ALREADY DONE
  //int regionvz=int((Vz+30)/6.);  // bin width is equal to 6 centi-meters
  //int regionvz = int((Vz+30.)/4.); // bin width is equal to 4 centi-meters
  //if(regionvz >= 15 || regionvz <= -1) return 999;   // FIXME! don't need this

  // 0-19, 20 bins (-40 to 40)
  int regionvz = int((Vz+40.)/4.); // bin width is equal to 4 centi-meters
  if(regionvz >= 20 || regionvz <= -1) return 999;

  return regionvz;
}
