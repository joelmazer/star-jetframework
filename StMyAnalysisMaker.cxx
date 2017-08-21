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
//      - general QA
//      
// it gets a pointer to a collection of jets  	
//
// ################################################################

#include "StMyAnalysisMaker.h"

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
///#include "StRoot/StPicoDstMaker/StPicoEvent.h"
///#include "StRoot/StPicoDstMaker/StPicoTrack.h"
//#include "StRoot/StPicoDstMaker/StPicoBTOWHit.h"
///#include "StRoot/StPicoDstMaker/StPicoBTofHit.h"
///#include "StRoot/StPicoDstMaker/StPicoEmcTrigger.h"
///#include "StRoot/StPicoDstMaker/StPicoMtdTrigger.h"
///#include "StRoot/StPicoDstMaker/StPicoEmcPidTraits.h" // OLD
///#include "StRoot/StPicoDstMaker/StPicoBEmcPidTraits.h" // NEW
///#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
///#include "StRoot/StPicoDstMaker/StPicoMtdPidTraits.h"
///#include "StRoot/StPicoDstMaker/StPicoConstants.h"


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

#include "StRoot/StPicoEvent/StPicoConstants.h"

// centrality
//#include "eventStructuredAu/TStarJetPicoEventHeader.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// This is the minimal interface needed to access FastJet. 
// A more sophisticated interface is demonstrated in main72.cc.
//#include "FastJet/fastjet-install/include/fastjet/PseudoJet.hh"
//#include "FastJet/fastjet-install/include/fastjet/ClusterSequence.hh"
//#include "FastKet/fastjet-install/include/fastjet/Selector.hh"
//gSystem->AddIncludePath("$FASTJET/include");
/*
gSystem->Load("$FASTJET/lib/libfastjet.so");
gSystem->AddIncludePath("$FASTJET/include");
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

using namespace fastjet;
*/

#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/SISConePlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#include <fastjet/Selector.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/Subtractor.hh>
//from contrib package
/*
#include <fastjet/contrib/GenericSubtractor.hh>
#include <fastjet/contrib/ShapeWithComponents.hh>
#include <fastjet/contrib/ConstituentSubtractor.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/Njettiness.hh>
#include <fastjet/contrib/NjettinessPlugin.hh>
#include <fastjet/contrib/SoftDrop.hh>
*/

class StJetMakerTask;

namespace fastjet {
  class PseudoJet;
}

ClassImp(StMyAnalysisMaker)

//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName, bool mDoComments, double minJetPt = 1.0, double trkbias = 0.15, const char* jetMakerName = "", const char* rhoMakerName = "")
  : StMaker(name)
{
  //mPicoDstMaker = picoMaker; // don't need this because will be using multiple tasks
  fPoolMgr = 0x0;
  //tracksClone = 0x0;
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
  fDoEffCorr = kFALSE;
  fCorrJetPt = kFALSE;
  fMinPtJet = minJetPt;
  fTrackBias = trkbias;
  fTrackPtCut = 0.15;
  mCentrality = 0;
  fDoEventMixing = 0; fMixingTracks = 50000; fNMIXtracks = 5000; fNMIXevents = 5;
  fCentBinSize = 5; fReduceStatsCent = -1;
  //fTriggerEventType = StVEvent::kAny; fMixingEventType = StVEvent::kAny;
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


}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();

  //Create user objects.
  //  fRho = (StRhoParameter*)

  // clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
  //TClonesArray* tracksClone = new TClonesArray("StPicoTrack");
  //tracksClone = new TClonesArray("StPicoTrack");
  //tracksClone->SetName("tracksClone");
  //tracksClone->SetOwner(kTRUE);

  fTracksME = new TClonesArray("StPicoTrack");
  fTracksME->SetName("fTracksME");
  fTracksME->SetOwner(kTRUE);

  fJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);

  // may not need, used for old RUNS
  // StRefMultCorr* getgRefMultCorr() ; // For grefmult //Run14 AuAu200GeV
  grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  //refmultCorr = CentralityMaker::instance()->getRefMultCorr(); // OLD
  //refmult2Corr = CentralityMaker::instance()->getRefMult2Corr();  // OLD 

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

  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {
  double pi = 1.0*TMath::Pi();

  // some default pid histograms 
  mKsM     = new TH1F("KsM","",300,0.35,0.65);
  mLambdaM = new TH1F("LambdaM","",400,1.05,1.25);
  mLbarM   = new TH1F("LbarM","",400,1.05,1.25);
  mKsRM     = new TH1F("KsRM","",300,0.35,0.65);
  mLambdaRM = new TH1F("LambdaRM","",400,1.05,1.25);
  mLbarRM   = new TH1F("LbarRM","",400,1.05,1.25);
  mdedxvspt     = new TH2F("dedxvspt","",500,0.,5.0,500,0.,50.);
  mKsDecayL     = new TH2F("KsDecayL","",500,0.,5.0,500,0.,50.);
  mLambdaDecayL = new TH2F("LambdaDecayL","",500,0.,5.0,500,0.,50.);
  mLbarDecayL   = new TH2F("LbarDecayL","",500,0.,5.0,500,0.,50.);
  hTriggerPt = new TH1F("hTriggerPt","Trigger Pt",100,1,10);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
  hJetPt2 = new TH1F("hJetPt2", "Jet p_{T}", 100, 0, 100);
  hJetE = new TH1F("hJetE", "Jet energy distribution", 100, 0, 100);
  hJetEta = new TH1F("hJetEta", "Jet #eta distribution", 24, -1.2, 1.2);
  hJetPhi = new TH1F("hJetPhi", "Jet #phi distribution", 72, 0.0, 2*pi);
  hJetNEF = new TH1F("hJetNEF", "Jet NEF", 100, 0, 1);
  hJetArea = new TH1F("hJetArea", "Jet Area", 100, 0, 1);

  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi", "Jet-Hadron #Delta#eta-#Delta#phi", 56, -1.4, 1.4, 72, -0.5*pi, 1.5*pi);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);

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

  // temp FIXME
  Int_t nCentBins = 16;
  Double_t cenBins[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  Double_t* centralityBin = cenBins;
  //Int_t nZvBins  = 8+1+8;
  //Double_t vBins[] = {-40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40};
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
  mKsM->Write();
  mLambdaM->Write();
  mLbarM->Write();
  mKsRM->Write();
  mLambdaRM->Write();
  mLbarRM->Write();
  mdedxvspt->Write();
  mKsDecayL->Write();
  mLambdaDecayL->Write();
  mLbarDecayL->Write();
  hTriggerPt->Write();

  // jet histos
  hJetPt->Write();
  hJetCorrPt->Write();
  hJetPt2->Write();
  hJetE->Write();
  hJetEta->Write();
  hJetPhi->Write();
  hJetNEF->Write();
  hJetArea->Write();

  fHistJetHEtaPhi->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();

  // jet sparse
  fhnJH->Write();
  fhnMixedEvents->Write();
  fhnCorr->Write();
}

//----------------------------------------------------------------------------- 
// OLD user code says: //  Called every event after Make(). 
void StMyAnalysisMaker::Clear(Option_t *opt) {
/*
  delete [] fJets; fJets=0;
  delete [] fRho; fRho=0;
  delete [] fTracksME; fTracksME=0;
*/

  //delete fhnJH; fhnJH = 0;

/*
  delete mKsM;
  delete mLambdaM;
  delete mLbarM;
  delete mKsRM;
  delete mLambdaRM;
  delete mLbarRM;
  delete mdedxvspt;
  delete mKsDecayL;
  delete mLambdaDecayL;
  delete mLbarDecayL;
  delete hTriggerPt;
  delete hJetPt;
  delete hJetCorrPt;
  delete hJetPt2;
  delete hJetE;
  delete hJetEta;
  delete hJetPhi;
  delete hJetNEF;
  delete hJetArea;
  delete fHistJetHEtaPhi;
  delete fHistEventQA; 
  delete fHistEventSelectionQA; 
  delete fHistEventSelectionQAafterCuts;

  delete fhnJH;
  delete fhnMixedEvents;
  delete fhnCorr;
*/
}

//----------------------------------------------------------------------------- 
//  This method is called every event.
Int_t StMyAnalysisMaker::Make() {
  const double pi = 1.0*TMath::Pi();

  // clear mixed event tracks array: double check this!! TODO
  //tracksClone->Clear();
  //fTracksME->Clear();

  // update counter
  mEventCounter++;
  //cout<<"StMyANMaker event# = "<<mEventCounter<<endl;

  bool printInfo = kFALSE;
  bool firstEvent = kFALSE;

  // get PicoDstMaker 
  mPicoDstMaker = (StPicoDstMaker*)GetMaker("picoDst");
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    //return kStWarn;
    return kStFatal;
  }

  // construct PicoDst object from maker
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    //return kStWarn;
    return kStFatal;
  }

  // create pointer to PicoEvent 
  mPicoEvent = mPicoDst->event();
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
    //return kStFatal;
  }

  // fill Event Trigger QA
  //FillEventTriggerQA(fHistEventSelectionQA, trig);

  // trigger information: the below is all seemingly incomplete, coming from StPicoEvent class
  // trigger information from the PicoEvent object: they seem to not be set... same results
  //cout<<"istrigger = "<<mPicoEvent->isTrigger(450011)<<endl; //NEW
  //cout<<"istrigger = "<<mPicoEvent->isTrigger(450021)<<endl; // NEW
  //cout<<"trigger ids = "<<mPicoEvent->triggerIds()<<endl;
  //cout<<"trigger word = "<<mPicoEvent->triggerWord()<<"  trigword mtd = "<<mPicoEvent->triggerWordMtd()<<endl; // no longer used
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

  // get vertex 3 vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();
  
  // zvertex cut - per the Aj analysis
  if((1.0*TMath::Abs(zVtx)) > 40) return kStWarn; //kStFatal;

  // let me know the Run #, fill, and event ID
  Int_t RunId = mPicoEvent->runId();
  Int_t fillId = mPicoEvent->fillId();
  Int_t eventId = mPicoEvent->eventId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Float_t fZDCCoincidenceRate = mPicoEvent->ZDCx();
  //cout<<"RunID = "<<RunId<<endl;
  //cout<<"fillID = "<<fillId<<endl; 
  //cout<<"eventID = "<<eventId<<endl; // what is eventID?

  // get JetMaker
  JetMaker = (StJetMakerTask*)GetMaker(fJetMakerName);
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn; //kStFatal;
  }

  // if we have JetMaker, get jet collection associated with it
  fJets = JetMaker->GetJets();
  if(!fJets) {     
    cout<<"no fJets array"<<endl;
    return kStWarn; 
  }

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();
  //if(doComments) cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<endl;

  // create TClonesArray for use with mixed events
  //(StPicoTrack*)picoArray[picoTrack]
  //TClonesArray* fTracksME = mPicoDst->picoArray(picoTrack);  // CORRECT, but OLD
  TClonesArray* fTracksME = mPicoDst->picoArray(StPicoArrays::Tracks);  // Track ->Tracks Aug17: should be the NEW way
  //TClonesArray* fTracksME = mPicoDst->picoArray(1);
  if(!fTracksME) {
    cout<<"don't have mixed event TClonesArray.. thats weird, should be same as tracks."<<endl;
    return kStWarn;
  }

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

/*
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  cout<<"nEmcTrigger = "<<nEmcTrigger<<endl;
  //static StPicoEmcTrigger* emcTrigger(int i) { return (StPicoEmcTrigger*)picoArrays[picoEmcTrigger]->UncheckedAt(i); }
  for(int i = 0; i < nEmcTrigger; i++) {
    StPicoEmcTrigger *emcTrig = mPicoDst->emcTrigger(i);
    if(emcTrig) {
      //emcTrig->Print();
      cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc()<<endl;
      cout<<"isHT0: "<<emcTrig->isHT0()<<endl;
      cout<<"isHT1: "<<emcTrig->isHT1()<<endl;
      cout<<"isHT2: "<<emcTrig->isHT2()<<endl;
      cout<<"isHT3: "<<emcTrig->isHT3()<<endl;
      cout<<"isJP0: "<<emcTrig->isJP0()<<endl;
      cout<<"isJP1: "<<emcTrig->isJP1()<<endl;
      cout<<"isJP2: "<<emcTrig->isJP2()<<endl;

    } else {
      cout<<"don't have triggers"<<endl;
      return kStWarn;
    }
  }
*/

  // get centrality test::
//  double cent = TStarJetPicoEventHeader::GetCorrectedReferenceMultiplicity();
  //Int_t Pico::mCent_Year11_200GeV[nCen] = {10,21,41,72,118,182,266,375,441}; // Run11 200 GeV (Copy from Run10 200 GeV)
  int grefMult = mPicoEvent->grefMult();
  int refMult = mPicoEvent->refMult();
  //int mCentrality = centrality(refMult); // Run11: OLD
  //int mgCentrality = centrality(grefMult); // Run11: OLD
  //refmultCorrUtil->init(RunId); // Run11: OLD
  //refmultCorrUtil->initEvent(refmult, zVtx, zdcCoincidenceRate);  // Run11: OLD  
  // need to find out why grefMult and refMult are different //TODO

  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_refmult.txt
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_grefmult.txt

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
//  cout<<"grefMult = "<<grefMult<<"  refMult = "<<refMult<<"  "; //endl; //"  mCentrality = "<<mCentrality<<endl;
//  cout<<"refCorr2 = "<<refCorr2<<"  refCorr1 = "<<refCorr1<<"  refCorr0 = "<<refCorr0;
//  cout<<"  cent16 = "<<cent16<<"   cent9 = "<<cent9<<"  centbin = "<<centbin<<endl;
  cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<"  refCorr2 = "<<refCorr2<<"  grefMult = "<<grefMult<<"  centbin = "<<centbin<<endl;

  // to limit filling unused entries in sparse, only fill for certain centrality ranges
  // ranges can be different than functional cent bin setter
  Int_t cbin = -1;
  // need to figure out centrality first in STAR: TODO
  if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
  else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
  else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
  else if (centbin>5 && centbin<10) cbin = 4; // 30-50%
  else if (centbin>9 && centbin<16) cbin = 5; // 50-90%
  else cbin = -99;

  // get rho from event and fill relative histo's: "StRho_JetsBG", "OutRho", "StMaker#0"
  RhoMaker = (StRho*)GetMaker(fRhoMakerName);
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn; //kStFatal;
  }

  //fRho = GetRhoFromEvent(fRhoName);
  fRho = (StRhoParameter*)RhoMaker->GetRho();
  if(!fRho) {
    //Form("Couldn't get fRho named %s\n", fRhoName.Data());
    cout<<"Couldn't get fRho"<<endl;
    return kStWarn;    
  } 
  
  // get rho/area for standard area with R=0.2; fRho->ls("");
  fRhoVal = fRho->GetVal();
  //cout<<"   fRhoVal = "<<fRhoVal<<"   Correction = "<<1.0*TMath::Pi()*0.2*0.2*fRhoVal<<endl;

  // loop over Jets in the event: initialize some parameter variables
  double jetarea, jetpt, corrjetpt, jetptselected, jetE, jetEta, jetPhi, jetNEF, maxtrackpt, NtrackConstit;
  double gpt, gp, phi, eta, px, py, pt, pz, charge;
  double gphi, geta, gpt2, gpx, gpy, gpz;
  double deta, dphijh;
  Int_t ijethi = -1;
  Double_t highestjetpt = 0.0;
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // get our jets
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    // get some get parameters
    jetarea = jet->Area();
    jetpt = jet->Pt();
    corrjetpt = jet->Pt() - jetarea*fRhoVal;
    jetE = jet->E();
    jetEta = jet->Eta();
    jetPhi = jet->Phi();    
    jetNEF = jet->NEF();

    // some threshold cuts for tests
    if(fCorrJetPt) {
      if(corrjetpt < fMinPtJet) continue;
    } else { if(jetpt < fMinPtJet) continue; }
    //if(jet->MaxTrackPt() < fTrackBias) continue; // INCOMPLETE MEMBER //FIXME
    if(jet->GetMaxTrackPt() < fTrackBias) continue;

    // fill some histos
    hJetPt->Fill(jetpt);
    hJetCorrPt->Fill(corrjetpt);
    hJetE->Fill(jetE);
    hJetEta->Fill(jetEta);
    hJetPhi->Fill(jetPhi);
    hJetNEF->Fill(jetNEF);
    hJetArea->Fill(jetarea);

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

    Double_t CorrEntries[4] = {centbin, jetptselected, 0, zVtx};
    if(fReduceStatsCent > 0) {
      if(cbin == fReduceStatsCent) fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    } else fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    // ====================================================================================
    //} // check on max track and cluster pt/Et


    // track loop inside jet loop - loop over ALL tracks in PicoDst
    for(unsigned short itrack=0;itrack<ntracks;itrack++){
      // get tracks
      StPicoTrack* trk = mPicoDst->track(itrack);
      if(!trk){ continue; }

      // MAKE SURE TO figure out whether to select primary or global tracks
      //if(trk->gPt() < 0.15){ continue; }
 
      // declare kinematic variables
      //gpt = trk->gPt();
      //gp = trk->gPtot();
      if(doUsePrimTracks) {
        if(!(trk->isPrimary())) continue; // check if primary

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

      // do pt cut here to accommadate either type
      if(doUsePrimTracks) { // primary  track
        if(pt < fTrackPtCut) continue;
      } else { // global track
        if(pt < fTrackPtCut) continue;
      }

      // more acceptance cuts now - after getting 3vector - hardcoded for now
      //if(pt<0.15) continue;
      if((1.0*TMath::Abs(eta)) > 1.0) continue;      
      if(phi < 0) phi+= 2*pi;
      if(phi > 2*pi) phi-= 2*pi;
      if((phi < 0) || (phi > 2*pi)) continue;

      // get jet - track relations
      //deta = eta - jetEta;               // eta betweeen hadron and jet
      deta = jetEta - eta;               // eta betweeen jet and hadron
      dphijh = RelativePhi(jetPhi, phi); // angle between jet and hadron

      // fill jet sparse 
      Double_t triggerEntries[8] = {centbin, jetptselected, pt, deta, dphijh, 0, zVtx, charge};
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
    //cout<<"bgtracks"<<endl;

  ///FIXME if(trigger && fTriggerEventType) { //kEMCEJE)) {     
  bool runthisnow = kTRUE;
  if(runthisnow) {
    if (pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {

      // loop over Jets in the event
      double Mixjetarea, Mixjetpt, Mixcorrjetpt, Mixjetptselected, MixjetE, MixjetEta, MixjetPhi, MixjetNEF, Mixmaxtrackpt, MixNtrackConstit;
      double gMixpt, gMixp, Mixphi, Mixeta, Mixpx, Mixpy, Mixpt, Mixpz, Mixcharge;
      double dMixeta, dMixphijh;
      // loop over jets (passing cuts?)
      for (Int_t ijet = 0; ijet < njets; ijet++) {
        Double_t leadjet=0;
        if (ijet==ijethi) leadjet=1; //FIXME

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

        // some threshold cuts - do mixing only if we have a jet meeting out pt threshold and bias
        if(fCorrJetPt) {
          if(Mixcorrjetpt < fMinPtJet) continue;
        } else { if(Mixjetpt < fMinPtJet) continue; }

        //if(jet->GetMaxTrackPt() < fTrackBias) continue; 
   	//TODO if (!AcceptMyJet(jet)) continue;

        // get number of current events in pool
        Int_t nMix = pool->GetCurrentNEvents();
        //cout<<"nMix = "<<nMix<<endl;

        // Fill for biased jet triggers only
        //FIXME if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)) {  // && jet->Pt() > fJetPtcut) {
        if(jet->GetMaxTrackPt() > fTrackBias) {
          //cout<<"have jet with pt and tracks > bias"<<endl;

          // Fill mixed-event histos here: loop over nMix events
          for (Int_t jMix=0; jMix<nMix; jMix++) {
 
            // get jMix'th event
            bgTracks = pool->GetEvent(jMix);
            const Int_t Nbgtrks = bgTracks->GetEntries();
            //cout<<"nbgtrks = "<<Nbgtrks<<endl;         

            for(Int_t ibg=0; ibg<Nbgtrks; ibg++) {
              // get tracks
              //StPicoTrack* trk = mPicoDst->track(ibg);
              StPicoTrack* trk = static_cast<StPicoTrack*>(bgTracks->At(ibg));
              if(!trk){ continue; }

              // MAKE SURE TO figure out whether to select primary or global tracks TODO
              //if(trk->gPt() < 0.15){ continue; }

              // declare kinematic variables
              //gpt = trk->gPt();
              //gp = trk->gPtot();
              if(doUsePrimTracks) {
                if(!(trk->isPrimary())) continue; // check if primary

                // get primary track variables
                StThreeVectorF mPMomentum = trk->pMom();
                Mixphi = mPMomentum.phi();
                Mixeta = mPMomentum.pseudoRapidity();
                Mixpx = mPMomentum.x();
                Mixpy = mPMomentum.y();
                Mixpt = mPMomentum.perp();
                Mixpz = mPMomentum.z();
              } else {
                // get global track variables
                StThreeVectorF mgMomentum = trk->gMom(mVertex, Bfield);
                Mixphi = mgMomentum.phi();
                Mixeta = mgMomentum.pseudoRapidity();
                Mixpx = mgMomentum.x();
                Mixpy = mgMomentum.y();
                Mixpt = mgMomentum.perp();
                Mixpz = mgMomentum.z();
              }

              Mixcharge = trk->charge();

              // do pt cut here to accommadate either type
              if(doUsePrimTracks) { // primary  track
                if(pt < fTrackPtCut) continue;
              } else { // global track
                if(pt < fTrackPtCut) continue;
              }

              // more acceptance cuts now - after getting 3vector - hardcoded for now
              //if(Mixpt<0.15) continue;
              if((1.0*TMath::Abs(Mixeta)) > 1.0) continue;
              if(Mixphi < 0) Mixphi+= 2*pi;
              if(Mixphi > 2*pi) Mixphi-= 2*pi;
              if((Mixphi < 0) || (Mixphi > 2*pi)) continue;

              // get jet - track relations
              //deta = eta - jetEta;               // eta betweeen hadron and jet
              dMixeta = MixjetEta - Mixeta;               // eta betweeen jet and hadron
              dMixphijh = RelativePhi(MixjetPhi, Mixphi); // angle between jet and hadron
              //Double_t dEP = RelativeEPJET(jet->Phi(),fEPV0);	     // difference between jet and EP

              // calculate single particle tracking efficiency of mixed events for correlations
              Double_t mixefficiency = -999;
              //FIXME mixefficiency = EffCorrection(part->Eta(), part->Pt(), fDoEffCorr);                           
              mixefficiency = 1.0;

              if(fCorrJetPt) { Mixjetptselected = corrjetpt;
              } else { Mixjetptselected = jetpt; }

              // create / fill mixed event sparse
              Double_t triggerEntries[8] = {centbin, Mixjetptselected, Mixpt, dMixeta, dMixphijh, 0, zVtx, Mixcharge}; //array for ME sparse
              if(fReduceStatsCent > 0) {
                if(cbin == fReduceStatsCent) fhnMixedEvents->Fill(triggerEntries,1./(nMix*mixefficiency));   // fill Sparse histo of mixed events
              } else fhnMixedEvents->Fill(triggerEntries,1./(nMix*mixefficiency));   // fill Sparse histo of mixed events

            } // end of background track loop
          } // end of filling mixed-event histo's:  jth mix event loop
        } // end of check for biased jet triggers
      } // end of jet loop
    } // end of check for pool being ready
  } // end EMC triggered loop

    // use only tracks from MB and Central and Semi-Central events
    ///FIXME if(trigger && fMixingEventType) { //kMB) {

      TClonesArray* tracksClone2 = 0x0;
      //cout<<"tracks in clone before reduce: "<<tracksClone->GetEntriesFast();
      // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      tracksClone2 = CloneAndReduceTrackList(fTracksME);
      //cout<<"   entries after reduced = "<<tracksClone2->GetEntriesFast()<<"  ME = "<<fTracksME->GetEntriesFast()<<endl;

      // update pool if jet in event or not
      pool->UpdatePool(tracksClone2);
    ///FIXME } // MB and Central and Semi-Central events
  } // end of event mixing

// =======================================================================================================================

  // event counter at end of maker
  mInputEventCounter++;
  //cout<<"end of event counter = "<<mInputEventCounter<<endl;

  // fill Event Trigger QA
  //FillEventTriggerQA(fHistEventSelectionQAafterCuts, trig);

  return kStOK; //for tests = don't need rest of this functions code since passing collection of jets

  int ntracks2 = mPicoDst->numberOfTracks();
  if(doComments) cout << "ntracks = " << ntracks2 << endl;
  // loop over ALL tracks in PicoDst
  for(unsigned short itrack=0;itrack<ntracks2;itrack++){
    // will need to check if primary later
    StPicoTrack* trk = mPicoDst->track(itrack);
    if(!trk){ continue; }
    if(trk->gPt() < 0.2){ continue; }
             
    // fill some new histos
    hTriggerPt->Fill(trk->gPt());
    mdedxvspt->Fill(trk->gPt(),trk->dEdx());

    // some tests on tracks
    StThreeVectorF mPMomentum = trk->pMom();
    //cout<<"mass hypothesis = "<<mPMomentum.massHypothesis(1.0)<<endl;
  }

  // ======================================================================
  // ============================ Do some jet stuff =======================
  // Fastjet analysis - select algorithm and parameters
  // recombination schemes: 
  // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme

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
    //if(mPMomentum.x() > 3.0) cout<<"track# = "<<itrack<<"  Xmom = "<<mPMomentum.x()<<endl;

    // declare kinematic variables
    double pt = trk->gPt();
    double p = trk->gPtot();
    double phi = mPMomentum.phi();
    double eta = mPMomentum.pseudoRapidity();
    double px = mPMomentum.x();
    double py = mPMomentum.y();
    double pz = mPMomentum.z();

    // constituent cut of particles
    // Store as input to Fastjet
    // (px, py, pz, E) => use p for E for now
    fjInputs.push_back( fastjet::PseudoJet(px, py, pz, p ) );
  }

  // no jets (particle) in event
  if (fjInputs.size() == 0) {
    cout << "Error: event with no final state particles" << endl;
    //continue;
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

/*
    // Loop over jet constituents
    vector <fastjet::PseudoJet> constituents = sortedJets[i].constituents();
    vector <fastjet::PseudoJet> sortedconstituents = sorted_by_pt(constituents);
    for(unsigned j = 0; j < sortedconstituents.size(); j++) {
      if(printInfo) cout<<"    constituent "<<j<<"'s pt: "<<sortedconstituents[j].pt()<<endl;

      Double_t z = sortedconstituents[j].pt() / sortedJets[i].Et();
    }
*/
    //Double_t var[4] = {sortedJets[i].pt(), sortedJets[i].phi(), sortedJets[i].eta(), sortedconstituents.size()}; //,sortedJets[i].area()};
  }

  delete jetDef;
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
  Double_t dphi = (EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ){
    dphi = dphi + 1*TMath::Pi();
  } // this assumes we are doing full jets currently 
  
  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){
    dphi = dphi + 1*TMath::Pi();
  } 
  
  // test
  if( dphi < 0 || dphi > TMath::Pi()/2 )
    Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName());

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
  Int_t nMixTracks = tracksME->GetEntriesFast();
  Int_t iterTrk = 0;
  Double_t phi, eta, px, py, pt, pz, p, charge;
  const double pi = 1.0*TMath::Pi();
  for (Int_t i=0; i<nMixTracks; i++) { 
    // get tracks
    //StPicoTrack* trk = mPicoDst->track(i);
    StPicoTrack* trk = (StPicoTrack*)tracksME->At(i);
    if(!trk){ continue; }

    // MAKE SURE TO figure out whether to select primary or global tracks: HARD coded right now TODO
    //if(trk->gPt() < 0.15){ continue; }

    // declare kinematic variables
    //gpt = trk->gPt();
    //gp = trk->gPtot();
    if(doUsePrimTracks) {
      if(!(trk->isPrimary())) continue; // check if primary

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
      StThreeVectorF mgMomentum = trk->gMom(mVertex, Bfield);
      phi = mgMomentum.phi();
      eta = mgMomentum.pseudoRapidity();
      px = mgMomentum.x();
      py = mgMomentum.y();
      pt = mgMomentum.perp();
      pz = mgMomentum.z();
      p = mgMomentum.mag();
    }

    charge = trk->charge();

    // do pt cut here to accommadate either type
    if(doUsePrimTracks) { // primary  track
      if(pt < fTrackPtCut) continue;
    } else { // global track
      if(pt < fTrackPtCut) continue;
    }

    // more acceptance cuts now - after getting 3vector - hardcoded for now
    //if(pt < 0.15) continue; // if primary
    if((1.0*TMath::Abs(eta)) > 1.0) continue;
    if(phi < 0) phi+= 2*pi;
    if(phi > 2*pi) phi-= 2*pi;
    if((phi < 0) || (phi > 2*pi)) continue;

    // add tracks passing cuts to tracksClone array
    ((*tracksClone)[iterTrk]) = trk;
    ++iterTrk;
  } // end of looping through tracks

  //cout<<"tracksclone Ntracks = "<<tracksClone->GetEntriesFast()<<endl;
  return tracksClone;
}


/*
//________________________________________________________________________
Int_t StMyAnalysisMaker::AcceptMyTrack(StPicoTrack *trk) { //FIXME for jets
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

//-----------------------------------------------------------------------
// copied from the StPicoDstMaker class - set up a another parameter for datasets //TODO
Int_t StMyAnalysisMaker::centrality(int refMult) {
  //Int_t Pico::mCent_Year11_200GeV[nCen] = {10,21,41,72,118,182,266,375,441}; // Run11 200 GeV (Copy from Run10 200 GeV)
  for(int i=0;i<nCen;i++) {
//    if(refMult <= Pico::mCent_Year10_39GeV[i]) {
//    if(refMult <= Pico::mCent_Year10_7_7GeV[i]) {
//    if(refMult <= Pico::mCent_Year11_19_6GeV[i]) {
    if(refMult <= Pico::mCent_Year11_200GeV[i]) {
      return i;
    }
  }
  return nCen;
}

// Trigger QA histogram, label bins 
TH1* StMyAnalysisMaker::FillEventTriggerQA(TH1* h, UInt_t trig) {
  // check and fill a Event Selection QA histogram for different trigger selections after cuts
/*
  if(trig == 0) h->Fill(1);
  if(trig & AliVEvent::kAny) h->Fill(2);
  if(trig & AliVEvent::kAnyINT) h->Fill(3);
  if(trig & AliVEvent::kMB) h->Fill(4);
  if(trig & AliVEvent::kINT7) h->Fill(5);
  if(trig & AliVEvent::kEMC1) h->Fill(6);
  if(trig & AliVEvent::kEMC7) h->Fill(7);
  if(trig & AliVEvent::kEMC8) h->Fill(8);
  if(trig & AliVEvent::kEMCEJE) h->Fill(9);
  if(trig & AliVEvent::kEMCEGA) h->Fill(10);
  if(trig & AliVEvent::kCentral) h->Fill(11);
  if(trig & AliVEvent::kSemiCentral) h->Fill(12);
  if(trig & AliVEvent::kINT8) h->Fill(13);

  if(trig & (AliVEvent::kEMCEJE | AliVEvent::kMB)) h->Fill(14);
  if(trig & (AliVEvent::kEMCEGA | AliVEvent::kMB)) h->Fill(15);
  if(trig & (AliVEvent::kAnyINT | AliVEvent::kMB)) h->Fill(16);

  if(trig & (AliVEvent::kEMCEJE & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral))) h->Fill(17);
  if(trig & (AliVEvent::kEMCEGA & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral))) h->Fill(18);
  if(trig & (AliVEvent::kAnyINT & (AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral))) h->Fill(19);
*/
  // label bins of the analysis trigger selection summary
  h->GetXaxis()->SetBinLabel(1, "no trigger");
  h->GetXaxis()->SetBinLabel(2, "kAny");
  h->GetXaxis()->SetBinLabel(3, "kAnyINT");
  h->GetXaxis()->SetBinLabel(4, "kMB");
  h->GetXaxis()->SetBinLabel(5, "kINT7");
  h->GetXaxis()->SetBinLabel(6, "kEMC1");
  h->GetXaxis()->SetBinLabel(7, "kEMC7");
  h->GetXaxis()->SetBinLabel(8, "kEMC8");
  h->GetXaxis()->SetBinLabel(9, "kEMCEJE");
  h->GetXaxis()->SetBinLabel(10, "kEMCEGA");
  h->GetXaxis()->SetBinLabel(11, "kCentral");
  h->GetXaxis()->SetBinLabel(12, "kSemiCentral");
  h->GetXaxis()->SetBinLabel(13, "kINT8");
  h->GetXaxis()->SetBinLabel(14, "kEMCEJE or kMB");
  h->GetXaxis()->SetBinLabel(15, "kEMCEGA or kMB");
  h->GetXaxis()->SetBinLabel(16, "kAnyINT or kMB");
  h->GetXaxis()->SetBinLabel(17, "kEMCEJE & (kMB or kCentral or kSemiCentral)");
  h->GetXaxis()->SetBinLabel(18, "kEMCEGA & (kMB or kCentral or kSemiCentral)");
  h->GetXaxis()->SetBinLabel(19, "kAnyINT & (kMB or kCentral or kSemiCentral)");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}
