// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This code is set as an AnalysisMaker task, where it can perform:
// - event mixing setup
//
//
// ################################################################

#include "StEventPoolMaker.h"
#include "StRoot/StarRoot/StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TParameter.h"
#include <TProfile.h>
#include "TVector3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h" // NEW name
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  // NEW (OLD: StPicoEmcPidTraits.h)

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StEventPoolManager.h"
#include "StFemtoTrack.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

ClassImp(StEventPoolMaker)

//
//__________________________________________________________________________________________
StEventPoolMaker::StEventPoolMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName = "", bool mDoComments = kFALSE)
  : StJetFrameworkPicoBase(name)  //StMaker(name):
{
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  doRejectBadRuns = kFALSE;
  fHistCentBinMin = 0;
  fHistCentBinMax = 9;               // 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
  fHistZvertBinMin = 0;
  fHistZvertBinMax = 20;             // (-40, 40) 4cm bins
  fPoolMgr = 0x0;
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doPrintEventCounter = kFALSE;
  fDoEffCorr = kFALSE;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0; fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0; fTowerPhiMaxCut = 2.0*TMath::Pi();
  fDoEventMixing = 0; fMixingTracks = 50000; fNMIXtracks = 5000; fNMIXevents = 5;
  fCentBinSize = 5; fReduceStatsCent = -1;
  fDoUseMultBins = kFALSE;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2; fMixingEventType = 0;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4800; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }
  doComments = mDoComments;
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
  fEventPlaneMakerName = "";
}
//
//__________________________________________________________________________________________
StEventPoolMaker::~StEventPoolMaker()
{ /*  */
  // destructor
  if(hEventZVertex)       delete hEventZVertex;
  if(hCentrality)         delete hCentrality;
  if(hMultiplicity)       delete hMultiplicity;
  if(hTrackEtavsPhi)      delete hTrackEtavsPhi;
  if(fHistEventSelectionQA)           delete fHistEventSelectionQA;
  if(fHistEventSelectionQAafterCuts)  delete fHistEventSelectionQAafterCuts;
  if(hTriggerIds)         delete hTriggerIds;
  if(hEmcTriggers)        delete hEmcTriggers;
  if(hMixEvtStatZVtx)     delete hMixEvtStatZVtx;
  if(hMixEvtStatCent)     delete hMixEvtStatCent;
  if(hMixEvtStatZvsCent)  delete hMixEvtStatZvsCent;

  fPoolMgr->Clear(); delete fPoolMgr;
}
//
//__________________________________________________________________________________________
Int_t StEventPoolMaker::Init() {
  //StJetFrameworkPicoBase::Init();

  // initialize the histograms
  DeclareHistograms();

  return kStOK;
}
//
// Function: write to and close output file
//__________________________________________________________________________________________
Int_t StEventPoolMaker::Finish() { 
  cout << "StEventPoolMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
   
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StEventPoolMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}
//
// Function: declare histograms and set up some global objects
//__________________________________________________________________________________________
void StEventPoolMaker::DeclareHistograms() {
  // constants 
  double pi = 1.0*TMath::Pi();

  // cent histo bins
  int nHistCentBins;
  if(fCentBinSize == 20) nHistCentBins = 5;
  if(fCentBinSize == 10) nHistCentBins = 10;
  if(fCentBinSize ==  5) nHistCentBins = 20;

  // QA histos
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 100, -50, 50);
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", 160, 0, 800);
  hTrackEtavsPhi = new TH2F(Form("hTrackEtavsPhi"), Form("track distribution: #eta vs #phi"), 144, 0, 2*pi, 40, -1.0, 1.0);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);
  hMixEvtStatZVtx = new TH1F("hMixEvtStatZVtx", "no of events in pool vs zvtx", 20, -40.0, 40.0);
  hMixEvtStatCent = new TH1F("hMixEvtStatCent", "no of events in pool vs Centrality", nHistCentBins, 0, 100);
  hMixEvtStatZvsCent = new TH2F("hMixEvtStatZvsCent", "no of events: zvtx vs Centality", nHistCentBins, 0, 100, 20, -40.0, 40.0);

  // Setup for Au-Au collisions: cent bin size can only be 5 or 10% bins
  int nCentralityBins = 100;
  double mult = 1.0;
  if(fCentBinSize==1) { 
    nCentralityBins = 100;
    mult = 1.0;  
  } else if(fCentBinSize==2){
    nCentralityBins = 50;
    mult = 2.0;
  } else if(fCentBinSize==5){ // will be most commonly used
    nCentralityBins = 20;
    mult = 5.0;
  } else if(fCentBinSize==10){
    nCentralityBins = 10;
    mult = 10.0;
  } else if(fCentBinSize==20){
    nCentralityBins = 5;
    mult = 20.0;
  }

  // set bin edges
  Double_t cBins[nCentralityBins + 1]; // nCentralityBins
  for(Int_t ic = 0; ic < nCentralityBins + 1; ic++){
    //cBins[ic] = mult*ic;
    cBins[ic] = 1.0*ic;
  }

  Double_t *centralityBins = cBins;

  // multiplicity bins
  Int_t nMultBins = 25;  // Alt-2
  Double_t multBins[] = {10,15,21,31,42,53,66,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469};
  Double_t *multiplicityBins = multBins;

  // z-vertex bins for mixed events
  Int_t nZvBins  = 20;
  Double_t vBins[] = {-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40};
  Double_t *zvbins = vBins;
  // =================================================================================================
  //Int_t nbinsjetMIX = sizeof(vBinsJS)/sizeof(Double_t) - 1;

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implementation of AliEventPoolManager
  //fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinspp, centralityBinspp, nZvtxBins, zvtxbin);
  // binning schemes
  if(fDoUseMultBins) {
    fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBins, (Double_t*)multiplicityBins, nZvBins, (Double_t*)zvbins);
  } else {
    fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBins, (Double_t*)centralityBins, nZvBins, (Double_t*)zvbins);
  }

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}
//
// write histograms
//_____________________________________________________________________________
void StEventPoolMaker::WriteHistograms() {
  // default histos
  hEventZVertex->Write();
  hCentrality->Write();
  hMultiplicity->Write();
  hTrackEtavsPhi->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
  hMixEvtStatZVtx->Write();
  hMixEvtStatCent->Write();
  hMixEvtStatZvsCent->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StEventPoolMaker::Clear(Option_t *opt) {

}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StEventPoolMaker::Make() {
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent 
  mPicoEvent = static_cast<StPicoEvent*>(mPicoDst->event());
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // get base class pointer
  mBaseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!mBaseMaker) {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }

  // get bad run, dead & bad tower lists
  badRuns = mBaseMaker->GetBadRuns();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  // cut event on max track pt > 35.0 GeV (30 Oct25, 2018)
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  hEventZVertex->Fill(zVtx);

  // ============================ CENTRALITY ============================== //
  // get CentMaker pointer
  mCentMaker = static_cast<StCentMaker*>(GetMaker("CentMaker"));
  if(!mCentMaker) {
    LOG_WARN << " No CenttMaker! Skip! " << endm;
    return kStWarn;
  }

  // centrality variables
  int grefMult = mCentMaker->GetgrefMult(); // see StPicoEvent
  int refMult =  mCentMaker->GetrefMult();  // see StPicoEvent
  ref9 = mCentMaker->GetRef9();   // binning from central -> peripheral
  ref16 = mCentMaker->GetRef16(); // binning from central -> peripheral
  int cent16 = mCentMaker->GetCent16(); // centrality bin from StRefMultCorr (increasing bin corresponds to decreasing cent %) - Don't use except for cut below
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // bin-age to use for mixed event and sparses
  Int_t centbin10 = GetCentBin10(centbin);
  double centBinToUse;
  if(fCentBinSize==10) { centBinToUse = (double)centbin10 * 10.0;
  } else if(fCentBinSize==5) { centBinToUse = (double)centbin * 5.0; }

  // centrality / multiplicity histograms
  hCentrality->Fill(fCentralityScaled);
  hMultiplicity->Fill(refCorr2);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }
  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggersHist(hEmcTriggers);

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds(); 
  if(fDebugLevel == kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", "; 
  }
  if(fDebugLevel == kDebugEmcTrigger) cout<<endl;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30); 
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu
  if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

  // fill arrays for towers that fired trigger
  FillTowerTriggersArr();

  // ========================== Event pool setup ===================================== //
  // create pool pointer
  StEventPool *pool = 0x0;

  // require event mixing
  if(fDoEventMixing > 0) {
    // convert back to integer bins for mixed event pool - 10% bins (0, 7), 5% bins (0, 15)
    Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSize);

    // initialize event pools
    if(fDoUseMultBins) { pool = fPoolMgr->GetEventPool(refCorr2, zVtx);
    } else { pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); }
    if(!pool) {
      Form("No pool found for centrality = %i, zVtx = %f", mixcentbin, zVtx); // FIXME if cent changes to double
      return kTRUE;
    }

    // use only tracks from MB events
    if(fRunForMB) { // kMB5 or kMB30 (don't exclude HT)
    //if(fRunForMB && (!fHaveEmcTrigger)) { // kMB5 or kMB30 excluding HT
      // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());

    } // MB 
  }

  // ============================================================================================= //

  // get number of jets, tracks, and global tracks in events
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();

  // ***************************************************************************************************************
  // ******************************** Event MIXING *****************************************************************
  // ***************************************************************************************************************
  //Prepare to do event mixing
  if(fDoEventMixing>0){
    // event mixing

    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    // mix jets from triggered events with tracks from MB events
    // get the trigger bit, need to change trigger bits between different runs

    // mixed event centbin
    Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSize);
    //cout<<"mixcentbin = "<<mixcentbin<<"  centbin = "<<centbin<<" centbin10 = "<<centbin10<<"  zvtx = "<<zVtx<<endl;

    // initialize event pools
    StEventPool *pool = 0x0;
    if(fDoUseMultBins) { pool = fPoolMgr->GetEventPool(refCorr2, zVtx);
    } else { pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); } 
    if(!pool) {
      Form("No pool found for centrality = %i, zVtx = %f", mixcentbin, zVtx); // FIXME if cent changes to double
      return kTRUE;
    }

    ///if(fMixingEventType) { //kMB) {
    if(fRunForMB) { // kMB or kMB30 (don't exclude HT)
    //if(fRunForMB && (!fHaveEmcTrigger)) { // kMB or kMB30 (excluding HT)
      // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());

      // fill QA histo's
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB
  } // end of event mixing

  return kStOK;
}
//
//________________________________________________________________________________
TClonesArray *StEventPoolMaker::CloneAndReduceTrackList()
{
  // clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
  //  TClonesArray *tracksClone = new TClonesArray("StPicoTrack"); // original way
  TClonesArray *tracksClone = new TClonesArray("StFemtoTrack");
//  tracksClone->SetName("tracksClone");
//  tracksClone->SetOwner(kTRUE);

  // construct variables, get # of tracks
  int nMixTracks = mPicoDst->numberOfTracks();
  int iterTrk = 0;
  //const double pi = 1.0*TMath::Pi();

  // loop over tracks
  for(int i = 0; i < nMixTracks; i++) { 
    // get track pointer
    StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      if(!(trk->isPrimary())) continue; // check if primary
      mTrkMom = trk->pMom();                 // get primary track vector
    } else {
      mTrkMom = trk->gMom(mVertex, Bfield);  // get global track vector
    }

/*
    // track variables - used with alt method below
    double pt = mTrkMom.Perp();

    // when doing event plane calculation via pt assoc bin
    // this is TEMP, it will filter track by the pt bin used for analysis - double check syntax if using
    if(doTPCptassocBin && fDoFilterPtMixEvents) {
      if(fTPCptAssocBin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(fTPCptAssocBin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }
*/

    // create StFemtoTracks out of accepted tracks - light-weight object for mixing
    //  StFemtoTrack *t = new StFemtoTrack(pt, eta, phi, charge);
    StFemtoTrack *t = new StFemtoTrack(trk, Bfield, mVertex, doUsePrimTracks);
    if(!t) continue;

    // add light-weight tracks passing cuts to TClonesArray
    ((*tracksClone)[iterTrk]) =  t;

    //delete t;
    ++iterTrk;
  } // end of looping through tracks

  return tracksClone;
}
//
//
//_________________________________________________________________________
TH1 *StEventPoolMaker::FillEmcTriggersHist(TH1 *h) {
  // number of Emcal Triggers
  for(int i = 0; i < 8; i++) { fEmcTriggerArr[i] = 0; }
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT0 = emcTrig->isHT0();
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    bool isJP0 = emcTrig->isJP0();
    bool isJP1 = emcTrig->isJP1();
    bool isJP2 = emcTrig->isJP2();

    // fill for valid triggers
    if(isHT0) { h->Fill(1); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(isHT1) { h->Fill(2); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(isHT2) { h->Fill(3); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(isHT3) { h->Fill(4); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(isJP0) { h->Fill(5); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(isJP1) { h->Fill(6); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(isJP2) { h->Fill(7); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
  // kAny trigger - filled once per event
  h->Fill(10); 

  return h;
}
//
// Set the bin errors on histograms
// __________________________________________________________________________________
void StEventPoolMaker::SetSumw2() {
  // set sum weights
  hEventZVertex->Sumw2();
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hTrackEtavsPhi->Sumw2();
  fHistEventSelectionQA->Sumw2();
  fHistEventSelectionQAafterCuts->Sumw2();
  hTriggerIds->Sumw2();
  hEmcTriggers->Sumw2();
  hMixEvtStatZVtx->Sumw2();
  hMixEvtStatCent->Sumw2();
  hMixEvtStatZvsCent->Sumw2();
}
//
//
//_________________________________________________________________________
void StEventPoolMaker::FillTowerTriggersArr() {
  // tower - HT trigger types array: zero these out - so they are refreshed for each event
  for(int i = 0; i < 4800; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }

  // get number of Emc triggers
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // emc trigger parameters
    int emcTrigID = emcTrig->id();
    int emcTrigIDindex = emcTrigID - 1;

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigIDindex] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigIDindex] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigIDindex] = kTRUE;
  }
}
