// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This code is set as an AnalysisMaker task, where it can perform:
// 1) jet analysis
// 	- tagging
// 	- jet-hadron correlations
// 	- mixed events: use of an event pool to mix triggers with
//      - Rho (underlying event) subtraction to jets
//      - leading jet tag
//      - access to jet constituents
//      - general QA
//      
// can get a pointer to:
// 1) collection of jets  	
// 2) event wise rho parameter
// 3) jet constituents (4 vectors)
//
// ################################################################

#include "StMyAnalysisMaker.h"

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
//#include "StRoot/StPicoDstMaker/StPicoV0.h"

// my STAR includes
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPoolManager.h"

// new includes
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
//#include "StRoot/StPicoEvent/StPicoBTOWHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
//#include "StRoot/StPicoEvent/StPicoEmcPidTraits.h" // OLD
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  // NEW
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

#include "StPicoConstants.h"

// centrality
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// classes
class StJetMakerTask;

namespace fastjet {
  class PseudoJet;
}

ClassImp(StMyAnalysisMaker)

//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", bool mDoComments = kFALSE, double minJetPt = 1.0, double trkbias = 0.15, const char* jetMakerName = "", const char* rhoMakerName = "")
  : StMaker(name)
{
  fLeadingJet = 0; fExcludeLeadingJetsFromFit = 1.0; 
  fTrackWeight = 1; fEvenPlaneMaxTrackPtCut = 5.0;
  //mPicoDstMaker = picoMaker; // don't need this because will be using multiple tasks
  fPoolMgr = 0x0;
  fTracksME = 0x0;
  fJets = 0x0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0;
  refmultCorr = 0;
  refmult2Corr = 0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;  // see StMyAnalysisMaker::fRunFlagEnum
  fDoEffCorr = kFALSE;
  fCorrJetPt = kFALSE;
  fMinPtJet = minJetPt;
  fTrackBias = trkbias;
  fJetRad = 0.4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 20.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15; fTracknHitsRatio = 0.52; 
  fDoEventMixing = 0; fMixingTracks = 50000; fNMIXtracks = 5000; fNMIXevents = 5;
  fCentBinSize = 5; fReduceStatsCent = -1;
  fTriggerEventType = 0; fMixingEventType = 0;
  for(int i=0; i<7; i++) { fEmcTriggerArr[i] = 0; }
  doComments = mDoComments;
  fhnJH = 0x0;
  fhnMixedEvents = 0x0;
  fhnCorr = 0x0;
  fRho = 0x0;
  fRhoVal = 0;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;

  mEventCounter = 0;
  mAllPVEventCounter = 0;
  mInputEventCounter = 0;
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */
  // destructor
  delete hEventPlane;
  delete hEventZVertex;
  delete hTriggerPt;
  delete hJetPt;
  delete hJetCorrPt;
  delete hJetPt2;
  delete hJetE;
  delete hJetEta;
  delete hJetPhi;
  delete hJetNEF;
  delete hJetArea;
  delete hJetTracksPt;
  delete hJetTracksPhi;
  delete hJetTracksEta;
  delete hJetPtvsArea;
  delete fHistJetHEtaPhi;
  delete fHistEventSelectionQA;
  delete fHistEventSelectionQAafterCuts;
  delete hTriggerIds;
  delete hEmcTriggers;

  delete fhnJH;
  delete fhnMixedEvents;
  delete fhnCorr;

//  fJets->Clear(); delete fJets;
//  fRho->Clear(); delete fRho; 
  fPoolMgr->Clear(); delete fPoolMgr;
}

//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();

  //Create user objects.
  //  fRho = (StRhoParameter*)

/*
  // clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
  //TClonesArray* tracksClone = new TClonesArray("StPicoTrack");
  fTracksME = new TClonesArray("StPicoTrack");
  fTracksME->SetName("fTracksME");
  fTracksME->SetOwner(kTRUE);
*/

  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);

  // may not need, used for old RUNS
  // StRefMultCorr* getgRefMultCorr() ; // For grefmult //Run14 AuAu200GeV
  //refmultCorr = CentralityMaker::instance()->getRefMultCorr(); // OLD
  //refmult2Corr = CentralityMaker::instance()->getRefMult2Corr();  // OLD 
  if(fRunFlag == Run14_AuAu200) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr(); }
  if(fRunFlag == Run16_AuAu200) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id(); }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() { 
  //  Summarize the run.
  cout << "StMyAnalysisMaker::Finish()\n";
  //cout << "\tProcessed " << mEventCounter << " events." << endl;
  //cout << "\tWithout PV cuts: "<< mAllPVEventCounter << " events"<<endl;
  //cout << "\tInput events: "<<mInputEventCounter<<endl;

  //  Write histos to file and close it.
/*
    if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Close();
  }
*/

  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    //fout->ls();
    fout->cd();
    fout->mkdir(fAnalysisMakerName);
    fout->cd(fAnalysisMakerName);
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StMyAnalysisMaker::Finish"<<endl;

  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {
  double pi = 1.0*TMath::Pi();

  // some default pid histograms left
  hTriggerPt = new TH1F("hTriggerPt","Trigger Pt",100,1,10);

  // QA histos
  hEventPlane = new TH1F("hEventPlane", "Event Plane distribution", 144, 0.0, 1.0*pi);
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 100, -50, 50);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
  hJetPt2 = new TH1F("hJetPt2", "Jet p_{T}", 100, 0, 100);
  hJetE = new TH1F("hJetE", "Jet energy distribution", 100, 0, 100);
  hJetEta = new TH1F("hJetEta", "Jet #eta distribution", 24, -1.2, 1.2);
  hJetPhi = new TH1F("hJetPhi", "Jet #phi distribution", 72, 0.0, 2*pi);
  hJetNEF = new TH1F("hJetNEF", "Jet NEF", 100, 0, 1);
  hJetArea = new TH1F("hJetArea", "Jet Area", 100, 0, 1);
  hJetTracksPt = new TH1F("hJetTracksPt", "Jet track constituent p_{T}", 80, 0, 20.0);
  hJetTracksPhi = new TH1F("hJetTracksPhi", "Jet track constituent #phi", 72, 0, 2*pi);
  hJetTracksEta = new TH1F("hJetTracksEta", "Jet track constituent #eta", 56, -1.4, 1.4);
  hJetPtvsArea = new TH2F("hJetPtvsArea", "Jet p_{T} vs Jet area", 100, 0, 100, 100, 0, 1);

  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi", "Jet-Hadron #Delta#eta-#Delta#phi", 56, -1.4, 1.4, 72, -0.5*pi, 1.5*pi);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);

  // set up jet-hadron sparse
  UInt_t bitcodeMESE = 0; // bit coded, see GetDimParams() below
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; // | 1<<8 | 1<<9 | 1<<10;
  //if(fDoEventMixing) {
    fhnJH = NewTHnSparseF("fhnJH", bitcodeMESE);
  //}

  // set up centrality bins for mixed events
  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it
  Int_t nCentralityBinspp = 8;
  Double_t centralityBinspp[9] = {0.0, 4., 9, 15, 25, 35, 55, 100.0, 500.0};  

  // Setup for Au-Au collisions: cent bin size can only be 5 or 10% bins
  Int_t nCentralityBinsAuAu = 100;
  Double_t mult = 1.0;
  if(fCentBinSize==1) { 
    nCentralityBinsAuAu = 100;
    mult = 1.0;  
  } else if(fCentBinSize==2){
    nCentralityBinsAuAu = 50;
    mult = 2.0;
  } else if(fCentBinSize==5){ // will be most commonly used
    nCentralityBinsAuAu = 20;
    mult = 5.0;
  } else if(fCentBinSize==10){
    nCentralityBinsAuAu = 10;
    mult = 10.0;
  }

  // not used right now
  Double_t centralityBinsAuAu[nCentralityBinsAuAu]; // nCentralityBinsAuAu
  for(Int_t ic=0; ic<nCentralityBinsAuAu; ic++){
   centralityBinsAuAu[ic]=mult*ic;
  }

  // temp FIXME:  Currently 5% centrality bins 0-80%, 4cm z-vtx bins
  Int_t nCentBins = 16;
  Double_t cenBins[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  Double_t* centralityBin = cenBins;
  Int_t nZvBins  = 10+1+10;
  Double_t vBins[] = {-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40};
  //Int_t nZvBins  = 2+1+2;
  //Double_t vBins[] = {-40, -20, -0, 20, 40};
  Double_t* zvbin = vBins;

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
  //fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinspp, centralityBinspp, nZvtxBins, zvtxbin);
  //fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinsAuAu, centralityBinsAuAu, nZvtxBins, zvtxbin);
  fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBins, centralityBin, nZvBins, zvbin);

  // set up event mixing sparse
  //if(fDoEventMixing){
    bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; // | 1<<8 | 1<<9;
    fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", bitcodeMESE);
  //} // end of do-eventmixing

  UInt_t bitcodeCorr = 0; // bit coded, see GetDimparamsCorr() below
  bitcodeCorr = 1<<0 | 1<<1 | 1<<2 | 1<<3; // | 1<<4;
  fhnCorr = NewTHnSparseFCorr("fhnCorr", bitcodeCorr);

/*
  // REWRITE THIS FOR STAR // TODO
  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    TH2 *h2 = dynamic_cast<TH2*>(fOutput->At(i));
    if (h2){
      h2->Sumw2();
      continue;
    }
    TH3 *h3 = dynamic_cast<TH3*>(fOutput->At(i));
    if (h3){
      h3->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }
*/

  fhnJH->Sumw2();  
  fhnMixedEvents->Sumw2();
  fhnCorr->Sumw2(); 
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
  // default histos
  hTriggerPt->Write();
  hEventPlane->Write();
  hEventZVertex->Write();

  // jet histos
  hJetPt->Write();
  hJetCorrPt->Write();
  hJetPt2->Write();
  hJetE->Write();
  hJetEta->Write();
  hJetPhi->Write();
  hJetNEF->Write();
  hJetArea->Write();
  hJetTracksPt->Write();
  hJetTracksPhi->Write();
  hJetTracksEta->Write();
  hJetPtvsArea->Write();

  fHistJetHEtaPhi->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();

  // jet sparse
  fhnJH->Write();
  fhnMixedEvents->Write();
  fhnCorr->Write();
}

//----------------------------------------------------------------------------- 
// OLD user code says: //  Called every event after Make(). 
void StMyAnalysisMaker::Clear(Option_t *opt) {
  fJets->Clear();
  //fRho->Clear();

/*
  delete [] fTracksME; fTracksME=0;
  delete fhnJH; fhnJH = 0;
*/

}

//----------------------------------------------------------------------------- 
//  This method is called every event.
Int_t StMyAnalysisMaker::Make() {
  const double pi = 1.0*TMath::Pi();
  bool printInfo = kFALSE;
  bool firstEvent = kFALSE;
  bool fHaveEmcTrigger = kFALSE;
  bool fHaveMBevent = kFALSE;

  //StMemStat::PrintMem("MyAnalysisMaker at beginning of make");

  // clear mixed event tracks array: double check this!! TODO
  //fTracksME->Clear();

  // update counter
  mEventCounter++;
  //cout<<"StMyANMaker event# = "<<mEventCounter<<endl;

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

  // trigger information: the below is all seemingly incomplete, coming from StPicoEvent class
  // trigger information from the PicoEvent object: they seem to not be set... same results
  //cout<<"istrigger = "<<mPicoEvent->isTrigger(450011)<<endl; //NEW
  //cout<<"istrigger = "<<mPicoEvent->isTrigger(450021)<<endl; // NEW
  //cout<<"trigger ids = "<<mPicoEvent->triggerIds()<<endl;
/*
  cout<<"MB: "<<mPicoEvent->isMinBias()
  <<"  MBlow: "<<mPicoEvent->isMBSlow()
  <<"  Central: "<<mPicoEvent->isCentral()
  <<"  HT: "<<mPicoEvent->isHT()
  <<"  HT11: "<<mPicoEvent->isHT11()
  <<"  HT15: "<<mPicoEvent->isHT15()
  <<"  HT18: "<<mPicoEvent->isHT18()
  <<"  MtdTrig: "<<mPicoEvent->isMtdTrig()
  <<"  DiMuon: "<<mPicoEvent->isDiMuon()
  <<"  DiMuonHFT: "<<mPicoEvent->isDiMuonHFT()
  <<"  SingleMuon: "<<mPicoEvent->isSingleMuon()
  <<"  EMuon: "<<mPicoEvent->isEMuon()<<endl;
*/

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();
  
  // zvertex cut - per the Aj analysis
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStWarn;
  hEventZVertex->Fill(zVtx);

  // let me know the Run #, fill, and event ID
  Int_t RunId = mPicoEvent->runId();
  Int_t fillId = mPicoEvent->fillId();
  Int_t eventId = mPicoEvent->eventId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Float_t fZDCCoincidenceRate = mPicoEvent->ZDCx();
  //cout<<"RunID = "<<RunId<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

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

  // create TClonesArray for use with mixed events
  //TClonesArray* fTracksME = mPicoDst->picoArray(picoTrack);  // CORRECT, but OLD
  //TClonesArray* fTracksME = mPicoDst->picoArray(StPicoArrays::Tracks);  // Track ->Tracks Aug17: should be the NEW way
/*
  fTracksME = mPicoDst->picoArray(StPicoArrays::Tracks);  // Track ->Tracks Aug17: should be the NEW way
  if(!fTracksME) { return kStWarn; }
*/

  //static StPicoTrack* track(int i) { return (StPicoTrack*)picoArrays[picoTrack]->UncheckedAt(i); }
/*
  // printing available information from the PicoDst objects
  StPicoTrack* t = mPicoDst->track(1);
  if(t) t->Print();
  StPicoEmcTrigger *emcTrig = mPicoDst->emcTrigger(0);
  if(emcTrig) emcTrig->Print();
  StPicoMtdTrigger *mtdTrig = mPicoDst->mtdTrigger(0);
  if(mtdTrig) mtdTrig->Print();
  StPicoBTOWHit *btowhit = mPicoDst->btowHit(0); // OLD NAME
  if(btowhit) btowhit->Print();
  StPicoBTofHit *btofhit = mPicoDst->btofHit(0);
  if(btofhit) btofhit->Print();
  //StPicoMtdHit *mtdhit = mPicoDst->mtdHit(0);
  //mtdhit->Print();
  StPicoEmcPidTraits *emcpid = mPicoDst->emcPidTraits(0); // OLD NAME now its StPicoBEmcPidTraits
  if(emcpid) emcpid->Print();
  StPicoBTofPidTraits *tofpid = mPicoDst->btofPidTraits(0);
  if(tofpid) tofpid->Print();
  StPicoMtdPidTraits *mtdpid = mPicoDst->mtdPidTraits(0);
  if(mtdpid) mtdpid->Print();
*/

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggersHist(hEmcTriggers);

  // get trigger IDs from PicoEvent class and loop over them
  Int_t trId[20]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds(); 
  if(fDebugLevel == 7) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == 7) cout<<"i = "<<i<<": "<<mytriggers[i] << ", "; 
    if(mytriggers[i] == 450014) { fHaveMBevent = kTRUE; }   
    //if(mytriggers[i] == -9999) { fHaveEmcTrigger = kTRUE; } // TODO - may not use this

    // fill trigger array with trigger IDs
    trId[i] = mytriggers[i];
  }
  if(fDebugLevel == 7) cout<<endl;

  // ======================== end of Triggers ============================= //

  // ============================ CENTRALITY ============================== //
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_refmult.txt
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_grefmult.txt
  int grefMult = mPicoEvent->grefMult();
  int refMult = mPicoEvent->refMult();
  grefmultCorr->init(RunId);
  grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate);
//  if(grefmultCorr->isBadRun(RunId)) cout << "Run is bad" << endl; 
//  if(grefmultCorr->isIndexOk()) cout << "Index Ok" << endl;
//  if(grefmultCorr->isZvertexOk()) cout << "Zvertex Ok" << endl;
//  if(grefmultCorr->isRefMultOk()) cout << "RefMult Ok" << endl;
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  Int_t cent16 = grefmultCorr->getCentralityBin16();
  Int_t cent9 = grefmultCorr->getCentralityBin9();
  Int_t centbin = GetCentBin(cent16, 16);
  Double_t refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
  Double_t refCorr1 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 1);
  Double_t refCorr0 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 0);
//  cout<<"grefMult = "<<grefMult<<"  refMult = "<<refMult<<"  "; //endl;
//  cout<<"refCorr2 = "<<refCorr2<<"  refCorr1 = "<<refCorr1<<"  refCorr0 = "<<refCorr0;
//  cout<<"  cent16 = "<<cent16<<"   cent9 = "<<cent9<<"  centbin = "<<centbin<<endl;
//  cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<"  refCorr2 = "<<refCorr2<<"  grefMult = "<<grefMult<<"  centbin = "<<centbin<<endl;

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
  // cache the leading jet within acceptance
  fLeadingJet = GetLeadingJet();
  double rpAngle = GetReactionPlane();
  hEventPlane->Fill(rpAngle);

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
  //cout<<"   fRhoVal = "<<fRhoVal<<"   Correction = "<<1.0*TMath::Pi()*0.2*0.2*fRhoVal<<endl;

  // loop over Jets in the event: initialize some parameter variables
  double jetarea, jetpt, corrjetpt, jetptselected, jetE, jetEta, jetPhi, jetNEF, maxtrackpt, NtrackConstit;
  double gpt, gp, phi, eta, px, py, pt, pz, charge;
  double gphi, geta, gpt2, gpx, gpy, gpz;
  double deta, dphijh, dEP;
  Int_t ijethi = -1;
  Double_t highestjetpt = 0.0;

  //if((fDebugLevel == 7) && (fEmcTriggerArr[kIsHT2])) cout<<"Have HT2!! "<<fEmcTriggerArr[kIsHT2]<<endl;
  if((fDebugLevel == 7) && (fEmcTriggerArr[fTriggerEventType])) cout<<"Have selected trigger type for signal jet!!"<<endl;
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    if(!fEmcTriggerArr[fTriggerEventType]) continue;

    // get pointer to jets
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    // get some jet parameters
    jetarea = jet->Area();
    jetpt = jet->Pt();
    corrjetpt = jet->Pt() - jetarea*fRhoVal;
    jetE = jet->E();
    jetEta = jet->Eta();
    jetPhi = jet->Phi();    
    jetNEF = jet->NEF();
    dEP = RelativeEPJET(jet->Phi(), rpAngle);         // difference between jet and EP

    // some threshold cuts for tests
    if(fCorrJetPt) {
      if(corrjetpt < fMinPtJet) continue;
    } else { if(jetpt < fMinPtJet) continue; }
    //if(jet->MaxTrackPt() < fTrackBias) continue; // INCOMPLETE MEMBER //FIXME
    if(jet->GetMaxTrackPt() < fTrackBias) continue;

    // test QA stuff...
    //cout<<"ijet = "<<ijet<<"  dEP = "<<dEP<<"  jetpt = "<<jetpt<<"  corrjetpt = "<<corrjetpt<<"  maxtrackpt = "<<jet->GetMaxTrackPt()<<endl;
    //cout<<"Fragfunc Z = "<<jet->GetZ(px, py, pz)<<"  trk pt = "<<pt<<"  jet pt = "<<jetpt<<endl;

    // fill some histos
    hJetPt->Fill(jetpt);
    hJetCorrPt->Fill(corrjetpt);
    hJetE->Fill(jetE);
    hJetEta->Fill(jetEta);
    hJetPhi->Fill(jetPhi);
    hJetNEF->Fill(jetNEF);
    hJetArea->Fill(jetarea);
    hJetPtvsArea->Fill(jetpt, jetarea);

    // get nTracks and maxTrackPt
    maxtrackpt = jet->GetMaxTrackPt();
    NtrackConstit = jet->GetNumberOfTracks();
    if(doComments) cout<<"Jet# = "<<ijet<<"  JetPt = "<<jetpt<<"  JetE = "<<jetE<<endl;
    if(doComments) cout<<"MaxtrackPt = "<<maxtrackpt<<"  NtrackConstit = "<<NtrackConstit<<endl;

    // get highest Pt jet in event
    if(highestjetpt<jetpt){
      ijethi=ijet;
      highestjetpt=jetpt;
    }

    // the below track and cluster cut is already done (FOR TRACK only as still looking at CHARGED only)
    //if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
    // set up and fill Jet-Hadron trigger jets THnSparse
    // ====================================================================================
    if(fCorrJetPt) { jetptselected = corrjetpt; 
    } else { jetptselected = jetpt; }

    Double_t CorrEntries[4] = {centbin, jetptselected, dEP, zVtx};
    if(fReduceStatsCent > 0) {
      if(cbin == fReduceStatsCent) fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    } else fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    // ====================================================================================
    //} // check on max track and cluster pt/Et

    // get jet constituents
    //std::vector<fastjet::PseudoJet>  GetJetConstituents(ijet)
    //const std::vector<TLorentzVector> aGhosts = jet->GetGhosts();
    const std::vector<fastjet::PseudoJet> myConstituents = jet->GetMyJets();
    //cout<<"Jet = "<<ijet<<"  pt = "<<jetpt<<endl;
    for (UInt_t i=0; i<myConstituents.size(); i++) {
      // cut on ghost particles
      if(myConstituents[i].perp() < 0.1) continue; 
      //cout<<"  pt trk "<<i<<" = "<<myConstituents[i].perp()<<endl;

      hJetTracksPt->Fill(myConstituents[i].perp());
      hJetTracksPhi->Fill(myConstituents[i].phi());
      hJetTracksEta->Fill(myConstituents[i].eta());
    }
    //cout<<endl;

    // track loop inside jet loop - loop over ALL tracks in PicoDst
    for(unsigned short itrack=0;itrack<ntracks;itrack++){
      // get tracks
      StPicoTrack* trk = mPicoDst->track(itrack);
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

      // declare kinematic variables
      //gpt = trk->gPt();
      //gp = trk->gPtot();
      if(doUsePrimTracks) {
        // get primary track variables
        StThreeVectorF mPMomentum = trk->pMom();
        phi = mPMomentum.phi();
        eta = mPMomentum.pseudoRapidity();
        px = mPMomentum.x();
        py = mPMomentum.y();
        pt = mPMomentum.perp();
        pz = mPMomentum.z();
      } else {
        // get global track variables
        StThreeVectorF mgMomentum = trk->gMom(mVertex, Bfield);
        phi = mgMomentum.phi();
        eta = mgMomentum.pseudoRapidity();
        px = mgMomentum.x();
        py = mgMomentum.y();
        pt = mgMomentum.perp();
        pz = mgMomentum.z();
      }

      charge = trk->charge();

      //cout<<"gphi = "<<gphi<<"  phi = "<<phi<<"  geta = "<<geta<<"  eta = "<<eta<<"  gpt = "<<gpt<<"  pt = "<<pt<<"  gpt2 = "<<gpt2<<endl;
//      if(jet->ContainsTrack(itrack)) cout<<" track part of jet, pt = "<<pt<<endl;
//      const std::vector<TLorentzVector> constit =  jet->GetConstits();
//      cout<<" constit pt = "<<constit[0].perp()<<endl;

      // get jet - track relations
      //deta = eta - jetEta;               // eta betweeen hadron and jet
      deta = jetEta - eta;               // eta betweeen jet and hadron
      dphijh = RelativePhi(jetPhi, phi); // angle between jet and hadron

      // fill jet sparse 
      Double_t triggerEntries[8] = {centbin, jetptselected, pt, deta, dphijh, dEP, zVtx, charge};
      Double_t trefficiency = 1.0;
      //if(fDoEventMixing) {
        if(fReduceStatsCent > 0) {
          if(cbin == fReduceStatsCent) fhnJH->Fill(triggerEntries, 1.0/trefficiency);    // fill Sparse Histo with trigger entries
        } else fhnJH->Fill(triggerEntries, 1.0/trefficiency);
      //}

      fHistJetHEtaPhi->Fill(deta,dphijh);                          // fill jet-hadron  eta--phi distributio
    // =====================================================================================
    } // track loop

  } // jet loop

// ***************************************************************************************************************
// ******************************** Event MIXING *****************************************************************
// ***************************************************************************************************************

// =======================================================================================================================

  // initialize object array of cloned picotracks
  //TObjArray* tracksClone = 0x0;
//  TClonesArray* tracksClone = 0x0;

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

    // need to see if this event was selected by MIXED event trigger
    //FIXME UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    // if event was not selected (triggered) for any reseason (should never happen) then return
    //if (trigger==0)  return kTRUE; //FIXME

    // initialize event pools
    StEventPool* pool = 0x0;
    pool = fPoolMgr->GetEventPool(centbin, zVtx); //FIXME AuAu fcent: cent bin? cent16
    if (!pool) {
      Form("No pool found for centrality = %i, zVtx = %f", centbin, zVtx); // FIXME if cent changes to double
      return kTRUE; //FIXME
    }

    // initialize background tracks array
    TObjArray* bgTracks;

  ///FIXME if(trigger && fTriggerEventType) { //kEMCEJE)) {     
  bool runthisnow = kTRUE;
  //if(runthisnow) {
  // do event mixing when Signal Jet is part of event with a HT1 or HT2 trigger firing
  if(fEmcTriggerArr[fTriggerEventType]) {
    if (pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {

      // loop over Jets in the event
      double Mixjetarea, Mixjetpt, Mixcorrjetpt, Mixjetptselected, MixjetE, MixjetEta, MixjetPhi, MixjetNEF, Mixmaxtrackpt, MixNtrackConstit;
      double gMixpt, gMixp, Mixphi, Mixeta, Mixpx, Mixpy, Mixpt, Mixpz, Mixcharge;
      double dMixeta, dMixphijh, dEP;
      // loop over jets (passing cuts?)
      for (Int_t ijet = 0; ijet < njets; ijet++) {
        Double_t leadjet=0;
        if (ijet==ijethi) leadjet=1; //FIXME for leading jet

        // get jet object
        StJet *jet = static_cast<StJet*>(fJets->At(ijet));
        if (!jet) continue;

        // get some get parameters of jets for mixing
        Mixjetarea = jet->Area();
        Mixjetpt = jet->Pt();
        Mixcorrjetpt = Mixjetpt - Mixjetarea*fRhoVal;
        MixjetE = jet->E();
        MixjetEta = jet->Eta();
        MixjetPhi = jet->Phi();
        MixjetNEF = jet->NEF();
        dEP = RelativeEPJET(jet->Phi(), rpAngle);         // difference between jet and EP

        // some threshold cuts - do mixing only if we have a jet meeting out pt threshold and bias
        if(fCorrJetPt) {
          if(Mixcorrjetpt < fMinPtJet) continue;
        } else { if(Mixjetpt < fMinPtJet) continue; }

        //if(jet->GetMaxTrackPt() < fTrackBias) continue; 
   	//TODO if (!AcceptMyJet(jet)) continue;  // acceptance cuts done to jet in JetMaker

        // get number of current events in pool
        Int_t nMix = pool->GetCurrentNEvents();
        //cout<<"nMix = "<<nMix<<endl;

        // Fill for biased jet triggers only
        //FIXME if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)) {  // && jet->Pt() > fJetPtcut) {
        if(jet->GetMaxTrackPt() > fTrackBias) {
          // Fill mixed-event histos here: loop over nMix events
          for (Int_t jMix=0; jMix<nMix; jMix++) {
 
            // get jMix'th event
            bgTracks = pool->GetEvent(jMix);
            const Int_t Nbgtrks = bgTracks->GetEntries();

            for(Int_t ibg=0; ibg<Nbgtrks; ibg++) {
              // get tracks
              //StPicoTrack* trk = mPicoDst->track(ibg);
              StPicoTrack* trk = static_cast<StPicoTrack*>(bgTracks->At(ibg));
              if(!trk){ continue; }

              // acceptance and kinematic quality cuts
              if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }
              
              // declare kinematic variables
              if(doUsePrimTracks) {
                // get primary track variables
                StThreeVectorF mPMomentum = trk->pMom();
                Mixphi = mPMomentum.phi();
                Mixeta = mPMomentum.pseudoRapidity();
                Mixpt = mPMomentum.perp();
              } else {
                // get global track variables
                StThreeVectorF mgMomentum = trk->gMom(mVertex, Bfield);
                Mixphi = mgMomentum.phi();
                Mixeta = mgMomentum.pseudoRapidity();
                Mixpt = mgMomentum.perp();
              }

              Mixcharge = trk->charge();

              // get jet - track relations
              //deta = eta - jetEta;               // eta betweeen hadron and jet
              dMixeta = MixjetEta - Mixeta;               // eta betweeen jet and hadron
              dMixphijh = RelativePhi(MixjetPhi, Mixphi); // angle between jet and hadron

              // calculate single particle tracking efficiency of mixed events for correlations
              Double_t mixefficiency = -999;
              //FIXME mixefficiency = EffCorrection(part->Eta(), part->Pt(), fDoEffCorr);                           
              mixefficiency = 1.0;

              if(fCorrJetPt) { Mixjetptselected = corrjetpt;
              } else { Mixjetptselected = jetpt; }

              // create / fill mixed event sparse
              Double_t triggerEntries[8] = {centbin, Mixjetptselected, Mixpt, dMixeta, dMixphijh, dEP, zVtx, Mixcharge}; //array for ME sparse
              if(fReduceStatsCent > 0) {
                if(cbin == fReduceStatsCent) fhnMixedEvents->Fill(triggerEntries,1./(nMix*mixefficiency));   // fill Sparse histo of mixed events
              } else fhnMixedEvents->Fill(triggerEntries,1./(nMix*mixefficiency));   // fill Sparse histo of mixed events

            } // end of background track loop
          } // end of filling mixed-event histo's:  jth mix event loop
        } // end of check for biased jet triggers
      } // end of jet loop
    } // end of check for pool being ready
  } // end EMC triggered loop

    // use only tracks from MB and Semi-Central events
    ///FIXME if(trigger && fMixingEventType) { //kMB) {
    if(fHaveMBevent) { // kMB
      if(fDebugLevel == 6) cout<<"...MB event... update event pool"<<endl;

      TClonesArray* tracksClone2 = 0x0;
      //cout<<"tracks in clone before reduce: "<<tracksClone2->GetEntriesFast();
      // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      tracksClone2 = CloneAndReduceTrackList(fTracksME);
      //cout<<"   entries after reduced = "<<tracksClone2->GetEntriesFast()<<endl; //<<"  ME = "<<fTracksME->GetEntriesFast()<<endl;

      // update pool if jet in event or not
      pool->UpdatePool(tracksClone2);
    } ///FIXME } // MB and Central and Semi-Central events
    //pool->Clear();    delete pool;

  } // end of event mixing

// =======================================================================================================================

  // event counter at end of maker
  mInputEventCounter++;
  //cout<<"end of event counter = "<<mInputEventCounter<<endl;

  // fill Event Trigger QA
  //FillEventTriggerQA(fHistEventSelectionQAafterCuts, trig);

  //StMemStat::PrintMem("MyAnalysisMaker at end of make");

  return kStOK; //for tests = don't need rest of this functions code since passing collection of jets

  // ======================================================================
  // ============================ Do some jet stuff =======================
  // Fastjet analysis - select algorithm and parameters
  // recombination schemes: 
  // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme

/*
  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  //fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme; // was set as the default
  fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  fastjet::JetAlgorithm          algorithm;
  int    power   = -1;     // -1 = anti-kT; 0 = C/A; 1 = kT
  if (power == -1)      algorithm = fastjet::antikt_algorithm;
  if (power ==  0)      algorithm = fastjet::cambridge_algorithm;
  if (power ==  1)      algorithm = fastjet::kt_algorithm; // was set as the default

  // set jet definition
  jetDef = new fastjet::JetDefinition(algorithm, Rparam, recombScheme, strategy);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;

  // Reset Fastjet input
  fjInputs.resize(0);

  // loop over tracks
  for(unsigned short itrack=0;itrack<ntracks;itrack++){
    StPicoTrack* trk = mPicoDst->track(itrack);
    if(!trk){ continue; }
    if(trk->gPt() < 0.2){ continue; }

    // check if primary track
    if((!trk->isPrimary())) continue;

    StThreeVectorF mPMomentum = trk->pMom(); 

    // declare kinematic variables
    double pt = trk->gPt();
    double p = trk->gPtot();
    double phi = mPMomentum.phi();
    double eta = mPMomentum.pseudoRapidity();
    double px = mPMomentum.x();
    double py = mPMomentum.y();
    double pz = mPMomentum.z();

    // Store as input to Fastjet:    (px, py, pz, E) => use p for E for now
    fjInputs.push_back( fastjet::PseudoJet(px, py, pz, p ) );
  }

  // no jets (particle) in event
  if (fjInputs.size() == 0) {
    cout << "Error: event with no final state particles" << endl;
    return kFALSE;
  }

  // Run Fastjet algorithm
  vector <fastjet::PseudoJet> inclusiveJets, sortedJets, selJets;
  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

  // Extract inclusive jets sorted by pT (note minimum pT of 10.0 GeV)
  inclusiveJets = clustSeq.inclusive_jets(10.0);
  sortedJets    = sorted_by_pt(inclusiveJets);

  // apply some cuts for a subset of jets
  fastjet::Selector select_eta = fastjet::SelectorEtaRange(-0.5, 0.5); // selects |eta| < 0.5
  fastjet::Selector select_phi = fastjet::SelectorPhiRange(0, 3.1); // selecta 1.6 < phi < 2.94
  fastjet::Selector select_pt = fastjet::SelectorPtRange(10, 40); // selecta 20 < pt < 40
  fastjet::Selector select_Nhard = fastjet::SelectorNHardest(2);
  fastjet::Selector select_eta_phi = select_eta && select_phi; //&& select_Nhard;

  // fill subset
  selJets = select_eta_phi(sortedJets);

  // For the first event, print the FastJet details
  if (firstEvent) {
    cout << "Ran " << jetDef->description() << endl;
    cout << "Strategy adopted by FastJet was "
         << clustSeq.strategy_string() << endl << endl;
    firstEvent = false;
  }

  // loop over reconstructed jets
  ///if(printInfo) cout<<"       pt       eta       phi"<<endl;
  //cout<<"       pt       eta       phi"<<endl;
  for(unsigned i = 0; i< sortedJets.size(); i++) {
    ///if(printInfo) cout<<"jet " <<i<<": "<<sortedJets[i].pt()<<" "<<sortedJets[i].eta()<<" "<<sortedJets[i].phi()<<endl;
    if(i==0) { 
      cout<<"       pt       eta       phi"<<endl;
      cout<<"jet " <<i<<": "<<sortedJets[i].pt()<<" "<<sortedJets[i].eta()<<" "<<sortedJets[i].phi()<<endl;
    }

    // fill some jet histos
    hJetPt2->Fill(sortedJets[i].pt());

    // Loop over jet constituents
    vector <fastjet::PseudoJet> constituents = sortedJets[i].constituents();
    vector <fastjet::PseudoJet> sortedconstituents = sorted_by_pt(constituents);
    for(unsigned j = 0; j < sortedconstituents.size(); j++) {
      if(printInfo) cout<<"    constituent "<<j<<"'s pt: "<<sortedconstituents[j].pt()<<endl;

      Double_t z = sortedconstituents[j].pt() / sortedJets[i].Et();
    }

    //Double_t var[4] = {sortedJets[i].pt(), sortedJets[i].phi(), sortedJets[i].eta(), sortedconstituents.size()}; //,sortedJets[i].area()};
  }

  delete jetDef;

*/
  // ======================================================================

  return kStOK;
}

//________________________________________________________________________
Int_t StMyAnalysisMaker::GetCentBin(Int_t cent, Int_t nBin) const
{  // Get centrality bin.
  Int_t centbin = -1;

  if(nBin == 16) {
    if(cent == 0)  centbin = 15;
    if(cent == 1)  centbin = 14;
    if(cent == 2)  centbin = 13;
    if(cent == 3)  centbin = 12;
    if(cent == 4)  centbin = 11;
    if(cent == 5)  centbin = 10;
    if(cent == 6)  centbin = 9;
    if(cent == 7)  centbin = 8;
    if(cent == 8)  centbin = 7;
    if(cent == 9)  centbin = 6;
    if(cent == 10) centbin = 5;
    if(cent == 11) centbin = 4;
    if(cent == 12) centbin = 3;
    if(cent == 13) centbin = 2;
    if(cent == 14) centbin = 1;
    if(cent == 15) centbin = 0;
  }
  if(nBin == 9) {
    if(cent == 0)  centbin = 8;
    if(cent == 1)  centbin = 7;
    if(cent == 2)  centbin = 6;
    if(cent == 3)  centbin = 5;
    if(cent == 4)  centbin = 4;
    if(cent == 5)  centbin = 3;
    if(cent == 6)  centbin = 2;
    if(cent == 7)  centbin = 1;
    if(cent == 8)  centbin = 0;
  }

  return centbin;
}

// this function generate a jet name based on input
TString StMyAnalysisMaker::GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray* partCont, TClonesArray* clusCont, TString tag)
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
    ::Warning("StMyAnalysisMaker::GenerateJetName", "Unknown jet finding algorithm '%d'!", jetAlgo);
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
    ::Error("StMyAnalysisMaker::GenerateJetName", "Recombination %d scheme not recognized.", recoScheme);
  }

  TString name = TString::Format("%s_%s%s%s%s%s_%s",
      tag.Data(), algoString.Data(), typeString.Data(), radiusString.Data(), trackString.Data(), clusterString.Data(), recombSchemeString.Data());

  return name;
}

//________________________________________________________________________
Double_t StMyAnalysisMaker::RelativePhi(Double_t mphi,Double_t vphi) const
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
Double_t StMyAnalysisMaker::RelativeEPJET(Double_t jetAng, Double_t EPAng) const
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
    cout<<"dPhi not in range [0, 0.5*Pi]!"<<endl;
  }

  return dphi;   // dphi in [0, Pi/2]
}

//______________________________________________________________________
THnSparse* StMyAnalysisMaker::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()
   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){
         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }

      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF

void StMyAnalysisMaker::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "centrality 5% bin";
    /*
      nbins = 10;
      xmin = 0.;
      xmax = 100.;
    */ 
      // think about how I want to do this here TODO
      nbins = 16;
      xmin = 0.;
      xmax = 16.;     
      break;

   case 1:
      if(fCorrJetPt) { // correct jet pt
        label = "Jet Corrected p_{T}";
        nbins = 30;
        xmin = -50.;
        xmax = 100.;
      } else { // don't correct jet pt
        label = "Jet p_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      }
      break;

   case 2:
      label = "Track p_{T}";
      nbins = 80; 
      xmin = 0.;
      xmax = 20.;
      break;

   case 3:
      label = "Relative Eta";
      nbins = 72; // 48
      xmin = -1.8;
      xmax = 1.8;
      break;

   case 4: 
      label = "Relative Phi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;

   case 5:
      label = "Relative angle of Jet and Reaction Plane";
      nbins = 3; // (12) 72
      xmin = 0;
      xmax = 0.5*pi;
      break;

   case 6:
      label = "z-vertex";
      nbins = 20; // 10
      xmin = -40; //-10
      xmax =  40; //+10
      break;

   case 7:
      label = "track charge";
      nbins = 3;
      xmin = -1.5;
      xmax = 1.5;
      break;

   case 8:
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;

   case 9: // need to update
      label = "leading track";
      nbins = 10;
      xmin = 0;
      xmax = 50;
      break; 

   } // end of switch
} // end of getting dim-params

//______________________________________________________________________
THnSparse* StMyAnalysisMaker::NewTHnSparseFCorr(const char* name, UInt_t entries) {
  // generate new THnSparseD, axes are defined in GetDimParamsD()
  Int_t count = 0;
  UInt_t tmp = entries;
  while(tmp!=0){
    count++;
    tmp = tmp &~ -tmp;  // clear lowest bit
  }

  TString hnTitle(name);
  const Int_t dim = count;
  Int_t nbins[dim];
  Double_t xmin[dim];
  Double_t xmax[dim];

  Int_t i=0;
  Int_t c=0;
  while(c<dim && i<32){
    if(entries&(1<<i)){
      TString label("");
      GetDimParamsCorr(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }

    i++;
  }
  hnTitle += ";";

  return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF

//______________________________________________________________________________________________
void StMyAnalysisMaker::GetDimParamsCorr(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   //stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

    case 0:
      label = "centrality 5% bin";
    /*
      nbins = 10;
      xmin = 0.;
      xmax = 100.;
    */
      // think about how I want to do this here TODO
      nbins = 16;
      xmin = 0.;
      xmax = 16.;
      break;

    case 1:
      if(fCorrJetPt) { // correct jet pt
        label = "Jet Corrected p_{T}";
        nbins = 30;
        xmin = -50.;
        xmax = 100.;
      } else { // don't correct jet pt
        label = "Jet p_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      }
      break;

    case 2:
      label = "Relative angle: Jet and Reaction Plane";
      nbins = 3; // (12) 48
      xmin = 0.;
      xmax = 0.5*pi;
      break;

    case 3:
      label = "Z-vertex";
      nbins = 20;
      xmin = -40.;
      xmax = 40.;
      break;

    case 4: // may delete this case
      label = "Jet p_{T} corrected with Rho";
      nbins = 50; // 250
      xmin = -50.;
      xmax = 200.;  
      break;

   }// end of switch
} // end of Correction (ME) sparse


//_________________________________________________
// From CF event mixing code PhiCorrelations
TClonesArray* StMyAnalysisMaker::CloneAndReduceTrackList(TClonesArray* tracksME)
{
  // clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
  //TClonesArray* tracksClone = new TClonesArray;
  TClonesArray* tracksClone = new TClonesArray("StPicoTrack");
  //tracksClone->SetName("tracksClone");
  //tracksClone->SetOwner(kTRUE);

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3 vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  //Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nMixTracks = mPicoDst->numberOfTracks();
  //Int_t nMixTracks = tracksME->GetEntriesFast();
  Int_t iterTrk = 0;
  Double_t phi, eta, px, py, pt, pz, p, charge;
  const double pi = 1.0*TMath::Pi();
  for (Int_t i=0; i<nMixTracks; i++) { 
    // get tracks
    StPicoTrack* trk = mPicoDst->track(i);
    //StPicoTrack* trk = (StPicoTrack*)tracksME->At(i);
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // add tracks passing cuts to tracksClone array
    ((*tracksClone)[iterTrk]) = trk;
    ++iterTrk;
  } // end of looping through tracks

  //cout<<"tracksclone Ntracks = "<<tracksClone->GetEntriesFast()<<endl;
  return tracksClone;
}

//________________________________________________________________________
Bool_t StMyAnalysisMaker::AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert) {
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

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}

/*
//________________________________________________________________________
Int_t StMyAnalysisMaker::AcceptMyJet(StJet *jet) { // for jets
  //applies all jet cuts except pt
  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax)) return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax)) return 0;
  if (jet->Area()<fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0) return 0;
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100) return 0;

  //passed all above cuts
  return 1;
}
*/

//_________________________________________________________________________
TH1* StMyAnalysisMaker::FillEmcTriggersHist(TH1* h) {
  // number of Emcal Triggers
  for(int i=0; i<7; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  if(fDebugLevel == 7) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

  //static StPicoEmcTrigger* emcTrigger(int i) { return (StPicoEmcTrigger*)picoArrays[picoEmcTrigger]->UncheckedAt(i); }
  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    StPicoEmcTrigger *emcTrig = mPicoDst->emcTrigger(i);
    if(!emcTrig) continue;

    // print some EMCal Trigger info
    if(fDebugLevel == 7) {
      cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
      cout<<"  isHT0: "<<emcTrig->isHT0()<<"  isHT1: "<<emcTrig->isHT1();
      cout<<"  isHT2: "<<emcTrig->isHT2()<<"  isHT3: "<<emcTrig->isHT3();
      cout<<"  isJP0: "<<emcTrig->isJP0()<<"  isJP1: "<<emcTrig->isJP1()<<"  isJP2: "<<emcTrig->isJP2()<<endl;
    }

    // fill for valid triggers
    if(emcTrig->isHT0()) { h->Fill(1); fEmcTriggerArr[kIsHT0] = 1; }
    if(emcTrig->isHT1()) { h->Fill(2); fEmcTriggerArr[kIsHT1] = 1; }
    if(emcTrig->isHT2()) { h->Fill(3); fEmcTriggerArr[kIsHT2] = 1; }
    if(emcTrig->isHT3()) { h->Fill(4); fEmcTriggerArr[kIsHT3] = 1; }
    if(emcTrig->isJP0()) { h->Fill(5); fEmcTriggerArr[kIsJP0] = 1; }
    if(emcTrig->isJP1()) { h->Fill(6); fEmcTriggerArr[kIsJP1] = 1; }
    if(emcTrig->isJP2()) { h->Fill(7); fEmcTriggerArr[kIsJP2] = 1; }
  }
 
  h->GetXaxis()->SetBinLabel(1, "HT0");
  h->GetXaxis()->SetBinLabel(2, "HT1");
  h->GetXaxis()->SetBinLabel(3, "HT2");
  h->GetXaxis()->SetBinLabel(4, "HT3");
  h->GetXaxis()->SetBinLabel(5, "JP0");
  h->GetXaxis()->SetBinLabel(6, "JP1");
  h->GetXaxis()->SetBinLabel(7, "JP2");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}

//_____________________________________________________________________________
// Trigger QA histogram, label bins 
TH1* StMyAnalysisMaker::FillEventTriggerQA(TH1* h) {
  // check and fill a Event Selection QA histogram for different trigger selections after cuts

  // Run14 AuAu
  if(fRunFlag == Run14_AuAu200) {

    //Int_t trId[20]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
    int bin = 0;
    for(unsigned int i=0; i<mytriggers.size(); i++) {
      //cout<<"MyTriggers, i = "<<i<<": "<<mytriggers[i] << " "; 
      bin = 1; 

      if((mytriggers[i] == 450201) || (mytriggers[i] == 450211) || (mytriggers[i] == 460201)) { bin = 2; h->Fill(bin); } // HT1
      if((mytriggers[i] == 450202) || (mytriggers[i] == 450212)) { bin = 3; h->Fill(bin); }// HT2
      if((mytriggers[i] == 460203) || (mytriggers[i] == 6) || (mytriggers[i] == 10) || (mytriggers[i] == 14) || (mytriggers[i] == 31) || (mytriggers[i] == 450213)) { bin = 4; h->Fill(bin); } // HT3
      if(mytriggers[i] == 450014) { bin = 5; h->Fill(bin); }// MB
      if((mytriggers[i] == 20) || (mytriggers[i] == 450010) || (mytriggers[i] == 450020)) { bin = 6; h->Fill(bin); } // VPDMB-30
      if((mytriggers[i] == 26) || (mytriggers[i] == 450103)) { bin = 7; h->Fill(bin); }// Central-5
      if((mytriggers[i] == 15) || (mytriggers[i] == 460101) || (mytriggers[i] == 460111)) { bin = 8; h->Fill(bin); } // Central & Central-mon
      if((mytriggers[i] == 1) || (mytriggers[i] == 4) || (mytriggers[i] == 16) || (mytriggers[i] == 32) || (mytriggers[i] == 450005) || (mytriggers[i] == 450008) || (mytriggers[i] == 450009) || (mytriggers[i] == 450014) || (mytriggers[i] == 450015) || (mytriggers[i] == 450018) || (mytriggers[i] == 450024) || (mytriggers[i] == 450025) || (mytriggers[i] == 450050) || (mytriggers[i] == 450060)) { bin = 10; h->Fill(bin); }// VPDMB-5 
      if((mytriggers[i] == 20) || (mytriggers[i] == 20) || (mytriggers[i] == 450010) || (mytriggers[i] == 450020)) { bin = 11; h->Fill(bin); } // VPDMB-30
 
      // fill trigger selection histo
      //h->Fill(bin); // - do per trigger in case an above classification matches 2 simultaneously
    }

    // TODO need to do this once and not for every event!!!!!
    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
    h->GetXaxis()->SetBinLabel(2, "BHT1*VPDMB-30");
    h->GetXaxis()->SetBinLabel(3, "BHT2*VPDMB-30");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "VPDMB-30");
    h->GetXaxis()->SetBinLabel(7, "Central-5");
    h->GetXaxis()->SetBinLabel(8, "Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-30");
  }

  // Run16 AuAu
  if(fRunFlag == Run16_AuAu200) {
    //Int_t trId[20]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999,-999,-999,-999,-999,-999,-999,-999,-999,-999};
    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
    int bin = 0;
    for(unsigned int i=0; i<mytriggers.size(); i++) {
      //cout<<"MyTriggers, i = "<<i<<": "<<mytriggers[i] << " "; 
      bin = 1;

      if((mytriggers[i] == 450201) || (mytriggers[i] == 450211) || (mytriggers[i] == 460201)) { bin = 2; h->Fill(bin); } // HT1
      if((mytriggers[i] == 450202) || (mytriggers[i] == 450212)) { bin = 3; h->Fill(bin); }// HT2
      if((mytriggers[i] == 460203) || (mytriggers[i] == 6) || (mytriggers[i] == 10) || (mytriggers[i] == 14) || (mytriggers[i] == 31) || (mytriggers[i] == 450213)) { bin = 4; h->Fill(bin); } // HT3
      if(mytriggers[i] == 450014) { bin = 5; h->Fill(bin); }// MB
      if((mytriggers[i] == 20) || (mytriggers[i] == 450010) || (mytriggers[i] == 450020)) { bin = 6; h->Fill(bin); } // VPDMB-30
      if((mytriggers[i] == 26) || (mytriggers[i] == 450103)) { bin = 7; h->Fill(bin); }// Central-5
      if((mytriggers[i] == 15) || (mytriggers[i] == 460101) || (mytriggers[i] == 460111)) { bin = 8; h->Fill(bin); } // Central & Central-mon
      if((mytriggers[i] == 1) || (mytriggers[i] == 4) || (mytriggers[i] == 16) || (mytriggers[i] == 32) || (mytriggers[i] == 450005) || (mytriggers[i] == 450008) || (mytriggers[i] == 450009) || (mytriggers[i] == 450014) || (mytriggers[i] == 450015) || (mytriggers[i] == 450018) || (mytriggers[i] == 450024) || (mytriggers[i] == 450025) || (mytriggers[i] == 450050) || (mytriggers[i] == 450060)) { bin = 10; h->Fill(bin); }// VPDMB-5 
      if((mytriggers[i] == 20) || (mytriggers[i] == 20) || (mytriggers[i] == 450010) || (mytriggers[i] == 450020)) { bin = 11; h->Fill(bin); } // VPDMB-30

      // fill trigger selection histo
      //h->Fill(bin); // - do per trigger in case an above classification matches 2 simultaneously
    }

    // TODO need to do this once and not for every event!!!!!
    // label bins of the analysis trigger selection summary
    // TODO may need to re-label for run16
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
    h->GetXaxis()->SetBinLabel(2, "BHT1*VPDMB-30");
    h->GetXaxis()->SetBinLabel(3, "BHT2*VPDMB-30");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "VPDMB-30");
    h->GetXaxis()->SetBinLabel(7, "Central-5");
    h->GetXaxis()->SetBinLabel(8, "Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-30");
  }

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}

//________________________________________________________________________
Double_t StMyAnalysisMaker::GetReactionPlane() { 
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

  // loop over tracks
  n = mPicoDst->numberOfTracks();
  for (i=0; i<n; i++) {
    StPicoTrack* track = mPicoDst->track(i);
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below
    // may change this back in future
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // declare kinematic variables
    if(doUsePrimTracks) {
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

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3vector - hardcoded for now
    if(pt > fEvenPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    if(phi < 0) phi += 2*pi;
    if(phi > 2*pi) phi -= 2*pi;

    // check for leading jet removal - taken from Redmers approach (CHECK! TODO!)
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
  double psi = mQ.Phi() / order;
  //return psi*180/pi;  // converted to degrees
  return psi;
}

// _____________________________________________________________________________________________
StJet* StMyAnalysisMaker::GetLeadingJet(StRhoParameter* eventRho) {
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
            //if(!AcceptMyJet(jet)) continue;
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

/*
Trigger information for Physics stream:

Run14: AuAu
physics	
1	BBC_mb_fast_prepost	15062057	15064011
1	MB_sel1	15046070	15046094
1	MB_sel9	15045068	15045069
1	VPDMB-5	15075055	15076099
1	VPDMB-5-p-nobsmd-hlt	15081020	15090048
1	VPD_mb-emcped	15052044	15052047
4	BBC_mb_fast_prepost	15083023	15083023
4	BHT2*VPDMB-30	15076087	15076099
4	VPDMB-5-p-nobsmd-hlt	15090049	15167007
4	VPD_mb-emcped	15052058	15056025
5	BHT1*VPDMB-30	15076087	15076099
6	BHT3	15076073	15076099
6	ZDC-mb	15046070	15059070
7	ZDC-mb	15056026	15057020
8	BBC_mb	15046070	15046094
9	EHT	15077045	15167014
9	VPD_mb	15046070	15047100
9	VPD_mb-emcped	15056026	15070021
10	BBC_mb_fast_prepost	15064061	15070021
10	BHT3-L2Gamma	15076089	15077035
13	BHT2*BBC	15171058	15174040
14	BHT3*BBC	15171058	15174040
15	BHT1*BBC	15171058	15174040
15	Central-mon	15079048	15167014
16	VPDMB-5-p-nobsmd	15077042	15081025
19	ZDC-mb	15061031	15061032
20	VPDMB-30	15074070	15076099
21	MB-mon	15076072	15076099
26	Central-5	15075071	15079046
27	VPDMB-5-p-nobsmd-hltheavyfrag	15083024	15083024
28	Central-5-hltheavyfrag	15083024	15083024
30	Central-5-hltDiElectron	15083024	15083024
31	BHT3-hltDiElectron	15076089	15076092
32	VPDMB-5-ssd	15080029	15083062
88	VPDMB-5-p-nobsmd-hltDiElectron	15083024	15083024
440001	VPD_mb	15047102	15070021
440004	ZDC-mb	15057021	15070021
440005	BBC_mb	15046095	15049033
440015	BBC_mb	15049037	15070021
440050	MB_sel1	15046095	15047093
450005	VPDMB-5-p-nobsmd	15081026	15086018
450008	VPDMB-5	15076101	15167014
450009	VPDMB-5-p-nobsmd-ssd-hlt	15143032	15167014
450010	VPDMB-30	15076101	15167014
450011	MB-mon	15076101	15167014
450014	VPDMB-5-nobsmd	15077056	15167014
450015	VPDMB-5-p-nobsmd	15086042	15097029
450018	VPDMB-5	15153036	15167007
450020	VPDMB-30	15153036	15167007
450021	MB-mon	15153036	15167007
450024	VPDMB-5-nobsmd	15153036	15167007
450025	VPDMB-5-p-nobsmd	15097030	15167014
450050	VPDMB-5-p-nobsmd-hlt	15097030	15112023
450060	VPDMB-5-p-nobsmd-hlt	15112024	15167014
450103	Central-5	15079047	15167014
450201	BHT1*VPDMB-30	15076101	15167014
450202	BHT2*VPDMB-30	15076101	15167014
450203	BHT3	15076101	15167014
450211	BHT1*VPDMB-30	15153036	15167007
450212	BHT2*VPDMB-30	15153036	15167007
450213	BHT3	15153036	15167007
460101	Central	15170104	15174040
460111	Central	15174041	15187006
460201	BHT1*VPD	15174041	15174041
460202	BHT2*BBC	15174041	15175041
460203	BHT3	15174041	15187006
460212	BHT2*ZDCE	15175042	15187006


Run16: AuAu
physics	
15	BHT1-VPD-10	17132061	17142053
15	BHT1-VPD-30	17142054	17142058
16	BHT2-VPD-30	17134025	17142053
17	BHT2	17142054	17142058
17	BHT3	17132061	17169117
43	VPDMB-5-hlt	17169022	17170041
45	VPDMB-5-p-hlt	17038047	17041016
56	vpdmb-10	17038047	17038047
520001	VPDMB-5-p-sst	17039044	17056036
520002	VPDMB-5-p-nosst	17043050	17056017
520003	VPDMB-5	17038047	17039043
520007	vpdmb-10	17038050	17057056
520011	VPDMB-5-p-sst	17056051	17057056
520012	VPDMB-5-p-nosst	17056050	17056050
520013	VPDMB-5	17039044	17057056
520017	vpdmb-10	17058002	17076012
520021	VPDMB-5-p-sst	17058002	17076012
520022	VPDMB-5-p-nosst	17058022	17076010
520023	VPDMB-5	17058002	17076012
520027	vpdmb-10	17076040	17100024
520031	VPDMB-5-p-sst	17076043	17100024
520032	VPDMB-5-p-nosst	17076042	17099014
520033	VPDMB-5	17076040	17100024
520037	vpdmb-10	17100030	17130013
520041	VPDMB-5-p-sst	17100030	17130013
520042	VPDMB-5-p-nosst	17102032	17128029
520043	VPDMB-5	17100030	17130013
520051	VPDMB-5-p-sst	17109050	17111019
520051	VPDMB-5-sst	17169020	17169021
520052	VPDMB-5-nosst	17169020	17169021
520101	central-5	17040052	17056036
520103	central-novtx	17041033	17056036
520111	central-5	17056050	17057056
520113	central-novtx	17056050	17057056
520121	central-5	17058002	17076012
520123	central-novtx	17058002	17130013
520131	central-5	17076040	17100024
520141	central-5	17100030	17130013
520201	BHT1*VPDMB-10	17038050	17039043
520203	BHT3	17038047	17169021
520211	BHT1*VPDMB-10	17039044	17042044
520221	BHT1*VPDMB-10	17042074	17057056
520231	BHT1*VPDMB-10	17058002	17060034
520241	BHT1*VPDMB-10	17060035	17076012
520251	BHT1*VPDMB-10	17076042	17100024
520261	BHT1*VPDMB-10	17100030	17130013
520601	dimuon*VPDMB-10	17041033	17045018
520605	e-muon-BHT1	17041033	17042044
520606	e-muon-BHT0	17041033	17045011
520611	dimuon*VPDMB-10	17045039	17057056
520615	e-muon-BHT1	17042074	17057056
520616	e-muon-BHT0	17045039	17056036
520621	dimuon*VPDMB-10	17058002	17076012
520625	e-muon-BHT1	17058002	17060034
520626	e-muon-BHT0	17056050	17057056
520631	dimuon*VPDMB-10	17076040	17100024
520635	e-muon-BHT1	17060035	17076012
520636	e-muon-BHT0	17058002	17076012
520641	dimuon*VPDMB-10	17100030	17110037
520645	e-muon-BHT1	17076042	17100024
520646	e-muon-BHT0	17076042	17100024
520655	e-muon-BHT1	17100030	17130013
520656	e-muon-BHT0	17100030	17130013
520802	VPDMB-5-p-hlt	17041033	17057056
520812	VPDMB-5-p-hlt	17058002	17076012
520822	VPDMB-5-p-hlt	17076042	17100024
520832	VPDMB-5-hlt	17169020	17169021
520832	VPDMB-5-p-hlt	17100030	17130013
520842	VPDMB-5-p-hlt	17109050	17111019
530201	BHT1-VPD-10	17133038	17141005
530202	BHT2-VPD-30	17136037	17141005
530213	BHT3	17133038	17141005
540201	BHT1-VPD-30	17142059	17148003
540203	BHT2	17142059	17148003
550201	BHT1	17149053	17160009
560201	BHT1	17161024	17169018
560202	BHT1_rhic_feedback	17161024	17169018
570203	BHT3	17169118	17170012
570802	VPDMB-5-hlt	17169118	17172018
}
*/
