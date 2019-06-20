// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// Centrality QA
//
// ################################################################

#include "StCentralityQA.h"
#include "StRoot/StarRoot/StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h" 
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  

// jet-framework includes
#include "StJetFrameworkPicoBase.h"
#include "StFemtoTrack.h"

// old file kept
#include "StPicoConstants.h"

// extra includes
#include "StJetPicoDefinitions.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StCentralityQA)

//______________________________________________________________________________
StCentralityQA::StCentralityQA(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", bool mDoComments = kFALSE)
  : StJetFrameworkPicoBase(name)
{
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  fCentralityDef = 4; // see StJetFrameworkPicoBase::fCentralityDefEnum //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doUseBBCCoincidenceRate = kFALSE; // kFALSE = use ZDC
  fMaxEventTrackPt = 30.0;
  fHistCentBinMin = 0;
  fHistCentBinMax = 9;               // 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
  fHistZvertBinMin = 0;
  fHistZvertBinMax = 20;             // (-40, 40) 4cm bins
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  grefmultCorrNEW = 0x0;
  mOutName = outName;
  fDoEffCorr = kFALSE;
  doRejectBadRuns = kFALSE;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0; fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0; fTowerPhiMaxCut = 2.0*TMath::Pi();
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4801; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }
  doComments = mDoComments;
  fAnalysisMakerName = name;
}

//_____________________________________________________________________________
StCentralityQA::~StCentralityQA()
{ /*  */
  // destructor
  delete hEventZVertex;
  delete hCentrality;
  delete hMultiplicity;
  delete hCentralityNEW;
  delete hMultiplicityNEW;
  delete hCentralityDiff;
  delete hMultiplicityDiff;

  delete fHistEventSelectionQA;
  delete fHistEventSelectionQAafterCuts;
  delete hTriggerIds;
  delete hEmcTriggers;
}

//_____________________________________________________________________________
Int_t StCentralityQA::Init() {
  //StJetFrameworkPicoBase::Init();

  // initialize the histograms
  DeclareHistograms();

  // Add bad run lists
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp (200 GeV)
        AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2012_BadRuns_P12id.txt");
        break;
  
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu (200 GeV)
        //AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P17id.txt");
        AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2014_BadRuns_P18ih.txt");
        break; 
  
    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu (200 GeV)
        AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Y2016_BadRuns_P16ij.txt");
        break; 
  
    default :
      AddBadRuns("StRoot/StMyAnalysisMaker/runLists/Empty_BadRuns.txt");
  }

  // switch on Run Flag to look for firing trigger specifically requested for given run period
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run11_pp500 : // Run11: 500 GeV pp
        break;

    case StJetFrameworkPicoBase::Run12_pp200 : // Run12: 200 GeV pp
        break;

    case StJetFrameworkPicoBase::Run12_pp500 : // Run12: 500 GeV pp
        break;

    case StJetFrameworkPicoBase::Run13_pp510 : // Run13: 510 (500) GeV pp
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
              //grefmultCorrNEW = CentralityMaker::instance()->getgRefMultCorr_P17id_VpdMB30();
              grefmultCorrNEW = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30();
 
/*
        switch(fCentralityDef) {
          case StJetFrameworkPicoBase::kgrefmult :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P17id_VpdMB30();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30 :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P16id :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
              break;
          default: // this is the default for Run14
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
        }
        break;
*/
             break;

    case StJetFrameworkPicoBase::Run15_pp200 : // Run15: 200 GeV pp
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
        break;

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17: 510 (500) GeV pp
        // this is the default for Run17 pp - don't set anything for pp
        break;

    default :
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  }

  return kStOK;
}

//____________________________________________________________________________
Int_t StCentralityQA::Finish() { 
  cout << "StCentralityQA::Finish()\n";

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

  cout<<"End of StCentralityQA::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//_____________________________________________________________________________
void StCentralityQA::DeclareHistograms() {
  int nHistCentBins;
  int fCentBinSize = 5;
  if(fCentBinSize == 10) nHistCentBins = 10;
  if(fCentBinSize ==  5) nHistCentBins = 20;

  // QA histos
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 100, -50, 50);
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hCentralityNEW = new TH1F("hCentralityNEW", "No. events vs centrality - NEW", nHistCentBins, 0, 100);
  hCentralityDiff = new TH1F("hCentralityDiff", "No. events vs centrality - Diff", 2*nHistCentBins, -100, 100);

  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", 160, 0, 800);
  hMultiplicityNEW = new TH1F("hMultiplicityNEW", "No. events vs multiplicity - NEW", 160, 0, 800);
  hMultiplicityDiff = new TH1F("hMultiplicityDiff", "No. events vs multiplicity - Diff", 200, -100, 100);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}

//
// write histograms
//_____________________________________________________________________________
void StCentralityQA::WriteHistograms() {
  // default histos
  hEventZVertex->Write();
  hCentrality->Write();
  hCentralityNEW->Write();
  hCentralityDiff->Write();
  hMultiplicity->Write();
  hMultiplicityNEW->Write();
  hMultiplicityDiff->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StCentralityQA::Clear(Option_t *opt) {

}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StCentralityQA::Make() {
  //StMemStat::PrintMem("MyAnalysisMaker at beginning of make");

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

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !IsRunOK(fRunNumber) ) return kStOK;
  }

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower E > 30.0 GeV
  //if(GetMaxTowerE() > fMaxEventTowerE) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  hEventZVertex->Fill(zVtx);

  // let me know the Run #, fill, and event ID
  int RunId = mPicoEvent->runId();
  int fillId = mPicoEvent->fillId();
  int eventId = mPicoEvent->eventId();
  double fBBCCoincidenceRate = mPicoEvent->BBCx();
  double fZDCCoincidenceRate = mPicoEvent->ZDCx();
  if(fDebugLevel == kDebugGeneralEvt) cout<<"RunID = "<<RunId<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

  // ============================ CENTRALITY ============================== //
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_refmult.txt
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_grefmult.txt
  int grefMult = mPicoEvent->grefMult();
  //int refMult = mPicoEvent->refMult();
  Int_t centbin, cent9, cent16;
  Double_t refCorr2;

  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorr->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }
//    if(grefmultCorr->isBadRun(RunId)) cout << "Run is bad" << endl; 
//    if(grefmultCorr->isIndexOk()) cout << "Index Ok" << endl;
//    if(grefmultCorr->isZvertexOk()) cout << "Zvertex Ok" << endl;
//    if(grefmultCorr->isRefMultOk()) cout << "RefMult Ok" << endl;
    // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning

    // get centrality bin: either 0-7 or 0-15
    cent16 = grefmultCorr->getCentralityBin16();
    cent9 = grefmultCorr->getCentralityBin9();

    // re-order binning to be from central -> peripheral
    ref9 = GetCentBin(cent9, 9);
    ref16 = GetCentBin(cent16, 16);
    centbin = GetCentBin(cent16, 16);  // 0-16

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fZDCCoincidenceRate, 2); }

    //Double_t refCorr1 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 1);
    //Double_t refCorr0 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 0);
    //grefmultCorr->isCentralityOk(cent16)
  } else { // for pp
    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;
  }

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them

  // centrality / multiplicity histograms
  ///hMultiplicity->Fill(refCorr2);
  if(fDebugLevel == kDebugCentrality) { if(centbin > 15) cout<<"centbin = "<<centbin<<"  mult = "<<refCorr2<<"  Centbin*5.0 = "<<centbin*5.0<<"  cent16 = "<<cent16<<endl; }
  fCentralityScaled = centbin*5.0;
  ///hCentrality->Fill(fCentralityScaled);

  // =======================================================================
  // new centrality determination

  Int_t centbinNEW, cent9NEW, cent16NEW;
  Int_t ref9NEW, ref16NEW;
  Double_t refCorr2NEW;
  Double_t fCentralityScaledNEW;

  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorrNEW->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorrNEW->initEvent(grefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorrNEW->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }
//    if(grefmultCorrNEW->isBadRun(RunId)) cout << "Run is bad" << endl; 
//    if(grefmultCorrNEW->isIndexOk()) cout << "Index Ok" << endl;
//    if(grefmultCorrNEW->isZvertexOk()) cout << "Zvertex Ok" << endl;
//    if(grefmultCorrNEW->isRefMultOk()) cout << "RefMult Ok" << endl;
    // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning

    // get centrality bin: either 0-7 or 0-15
    cent16NEW = grefmultCorrNEW->getCentralityBin16();
    cent9NEW = grefmultCorrNEW->getCentralityBin9();

    // re-order binning to be from central -> peripheral
    ref9NEW = GetCentBin(cent9NEW, 9);
    ref16NEW = GetCentBin(cent16NEW, 16);
    centbinNEW = GetCentBin(cent16NEW, 16);  // 0-16

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { refCorr2NEW = grefmultCorrNEW->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ refCorr2NEW = grefmultCorrNEW->getRefMultCorr(grefMult, zVtx, fZDCCoincidenceRate, 2); }

    //Double_t refCorr1NEW = grefmultCorrNEW->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 1);
    //Double_t refCorr0NEW = grefmultCorrNEW->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 0);
    //grefmultCorrNEW->isCentralityOk(cent16NEW)
  } else { // for pp
    centbinNEW = 0, cent9NEW = 0, cent16NEW = 0, refCorr2NEW = 0.0, ref9NEW = 0, ref16NEW = 0;
  }

  // cut on unset centrality, > 80%
  if(cent16NEW == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them

  // centrality / multiplicity histograms
  ///hMultiplicityNEW->Fill(refCorr2NEW);
  if(fDebugLevel == kDebugCentrality) { if(centbinNEW > 15) cout<<"centbin = "<<centbinNEW<<"  mult = "<<refCorr2NEW<<"  CentbinNEW*5.0 = "<<centbinNEW*5.0<<"  cent16NEW = "<<cent16NEW<<endl; }
  fCentralityScaledNEW = centbinNEW*5.0;
  ///hCentralityNEW->Fill(fCentralityScaledNEW);

  // ====================================================================================
  // differences
  double refCorr2Diff = refCorr2 - refCorr2NEW;
  double fCentralityScaledDiff = fCentralityScaled - fCentralityScaledNEW;
  ///hMultiplicityDiff->Fill(refCorr2Diff);
  ///hCentralityDiff->Fill(fCentralityScaledDiff);

  // to limit filling unused entries in sparse, only fill for certain centrality ranges
  // ranges can be different than functional cent bin setter
  Int_t cbin = -1;
  // this is actually not used since the below line does the cut:  if(fRequireCentSelection)
  if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
  else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
  else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
  else if (centbin>5 && centbin<10) cbin = 4; // 30-50%
  else if (centbin>9 && centbin<16) cbin = 5; // 50-80%
  else cbin = -99;

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
  if(fDebugLevel == kDebugEmcTrigger) 
  cout<<endl;

  // check for MB/HT event
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30); 
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

  // fill arrays for towers that fired trigger
  FillTowerTriggersArr();
  // ======================== end of Triggers ============================= //

  if(fHaveMB30event || fHaveEmcTrigger) {
    // centrality / multiplicity histograms - OLD
    hMultiplicity->Fill(refCorr2);
    hCentrality->Fill(fCentralityScaled);

    // centrality / multiplicity histograms - NEW
    hMultiplicityNEW->Fill(refCorr2NEW);
    hCentralityNEW->Fill(fCentralityScaledNEW);

    // diff
    hMultiplicityDiff->Fill(refCorr2Diff);
    hCentralityDiff->Fill(fCentralityScaledDiff);

  }

  return kStOK;
}

//_________________________________________________________________________
TH1* StCentralityQA::FillEmcTriggersHist(TH1* h) {
  // number of Emcal Triggers
  for(int i = 0; i < 8; i++) { fEmcTriggerArr[i] = 0; }
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  //if(fDebugLevel == kDebugEmcTrigger) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

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

    // print some EMCal Trigger info
    if(fDebugLevel == kDebugEmcTrigger) {
      cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
      cout<<"  isHT0: "<<isHT0<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3;
      cout<<"  isJP0: "<<isJP0<<"  isJP1: "<<isJP1<<"  isJP2: "<<isJP2<<endl;
    }

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

  // set bin labels
  h->GetXaxis()->SetBinLabel(1, "HT0");
  h->GetXaxis()->SetBinLabel(2, "HT1");
  h->GetXaxis()->SetBinLabel(3, "HT2");
  h->GetXaxis()->SetBinLabel(4, "HT3");
  h->GetXaxis()->SetBinLabel(5, "JP0");
  h->GetXaxis()->SetBinLabel(6, "JP1");
  h->GetXaxis()->SetBinLabel(7, "JP2");
  h->GetXaxis()->SetBinLabel(10, "Any");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}
//
// Set the bin errors on histograms, set sum weights
// __________________________________________________________________________________
void StCentralityQA::SetSumw2() {
  //hEventZVertex->Sumw2();
  //hCentrality->Sumw2();
  //hMultiplicity->Sumw2();
  //hCentralityNEW->Sumw2();
  //hMultiplicityNEW->Sumw2();
  //hCentralityDiff->Sumw2();
  //hMultiplicityDiff->Sumw2();
  
  //fHistEventSelectionQA->Sumw2();
  //fHistEventSelectionQAafterCuts->Sumw2();
  //hTriggerIds->Sumw2();
  //hEmcTriggers->Sumw2();
}

//_________________________________________________________________________
void StCentralityQA::FillTowerTriggersArr() {
  // tower - HT trigger types array: zero these out - so they are refreshed for each event
  for(int i = 0; i < 4801; i++) {
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

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigID] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigID] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigID] = kTRUE;

    //cout<<"i = "<<i<<"  EmcTrigID = "<<emcTrigID<<"  adc = "<<emcTrig->adc()<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3<<endl;
  }

/*
  // loop over towers and add input vectors to fastjet
  int nTowers = mPicoDst->numberOfBTOWHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: get from index of array shifted by +1
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    //cout<<"itow = "<<itow<<"  towerID = "<<towerID<<"  HT1: "<<fTowerToTriggerTypeHT1[towerID]<<"  adc = "<<tower->adc()<<endl;
    //cout<<"itow = "<<itow<<"  towerID = "<<towerID<<"  HT2: "<<fTowerToTriggerTypeHT2[towerID]<<"  adc = "<<tower->adc()<<endl;
    //cout<<"itow = "<<itow<<"  towerID = "<<towerID<<"  HT3: "<<fTowerToTriggerTypeHT3[towerID]<<"  adc = "<<tower->adc()<<endl;
  }
*/

}
//
//
//____________________________________________________________________________
void StCentralityQA::ResetBadRunList( ){
  badRuns.clear();
}
//
// Add bad runs from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
//_________________________________________________________________________________
Bool_t StCentralityQA::AddBadRuns(TString csvfile){
  // open infile
  std::string line;
  std::ifstream inFile ( csvfile );

  __DEBUG(2, Form("Loading bad runs from %s", csvfile.Data()) );

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
        badRuns.insert( ientry );
        __DEBUG(2, Form("Added bad run # %d", ientry));
      }
    }
  }

  return kTRUE;
}
//
// Function: check on if Run is OK or not
//____________________________________________________________________________________________
Bool_t StCentralityQA::IsRunOK( Int_t mRunId ){
  //if( badRuns.size()==0 ){
  if( badRuns.empty() ){
    __ERROR("StCentralityQA::IsRunOK: WARNING: You're trying to run without a bad run list. If you know what you're doing, deactivate this throw and recompile.");
    throw ( -1 );
  }
  if( badRuns.count( mRunId )>0 ){
    __DEBUG(9, Form("Reject. Run ID: %d", mRunId));
    return kFALSE;
  } else {
    __DEBUG(9, Form("Accept. Run ID: %d", mRunId));
    return kTRUE;
  }
}
