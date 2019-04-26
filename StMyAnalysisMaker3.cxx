// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This code is set as an AnalysisMaker task, where it can perform:
// 1) jet analysis
// 	- tagging
// 	- jet-hadron correlations
//      - jet shape analysis
// 	- mixed events: use of an event pool to mix triggers with
//      - Rho (underlying event) subtraction to jets
//      - leading jet and subleading tag
//      - event plane calculation with BBC, ZDC, TPC
//      - event plane corrections with BBC, ZDC, TPC
//      - access to jet constituents
//      - general QA
//
// can get a pointer to:
// 1) collection of jets  	
// 2) event wise rho parameter
// 3) jet constituents (4 vectors)
// 4) leading + subleading jets
// 5) event plane calculation for various pt bins
// 6) externally read in Event Pool for event mixing
// ################################################################

#include "StMyAnalysisMaker3.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include <sstream>
#include <fstream>

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"

// my STAR includes
#include "StEventPlaneMaker.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPoolManager.h"
#include "StFemtoTrack.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu

// new includes
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h" // NEW name
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"  // NEW (OLD: StPicoEmcPidTraits.h)
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StMyAnalysisMaker3)

//______________________________________________________________________________
StMyAnalysisMaker3::StMyAnalysisMaker3(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", bool mDoComments = kFALSE, double minJetPt = 1.0, double trkbias = 0.15, const char* jetMakerName = "", const char* rhoMakerName = "")
  : StJetFrameworkPicoBase(name)  //StMaker(name): Oct3
{
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  doJetShapeAnalysis = kFALSE;
  fJetShapeJetType = kLeadingJets;  // see header, LJ, SubLJ, inclusive
  doRequireAjSelection = kFALSE;
  fCorrJetPt = kFALSE;
  fCentralityDef = 4; // see StJetFrameworkPicoBase::fCentralityDefEnum //(kgrefmult_P16id, default for Run16AuAu200)
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doWriteTrackQAHist = kTRUE;
  doWriteJetQAHist = kTRUE;
  doUseBBCCoincidenceRate = kFALSE; // kFALSE = use ZDC
  fMaxEventTrackPt = 30.0;
  fLeadingJet = 0x0; fSubLeadingJet = 0x0; fExcludeLeadingJetsFromFit = 1.0;
  fTrackWeight = 1; //StJetFrameworkPicoBase::kPtLinear2Const5Weight; // see StJetFrameworkPicoBase::EPtrackWeightType 
  fEventPlaneMaxTrackPtCut = 5.0;
  fTPCEPmethod = 1;
  fHistCentBinMin = 0;
  fHistCentBinMax = 9;               // 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
  fHistZvertBinMin = 0;
  fHistZvertBinMax = 20;             // (-40, 40) 4cm bins
  TPC_PSI2 = -999;
  TPCA_PSI2 = -999; // subevent A
  TPCB_PSI2 = -999; // subevent B
  BBC_PSI2 = -999; ZDC_PSI2 = -999; 
  BBC_PSI1 = -999; ZDC_PSI1 = -999;
  PSI2 = -999;
  RES = -999;
  TPC_raw_comb = 0.; TPC_raw_neg = 0.; TPC_raw_pos = 0.;
  BBC_raw_comb = 0.; BBC_raw_east = 0.; BBC_raw_west = 0.;
  ZDC_raw_comb = 0.; ZDC_raw_east = 0.; ZDC_raw_west = 0.;
  fPoolMgr = 0x0;
  fJets = 0x0;
  fRunNumber = 0;
  fEPTPCn = 0.; fEPTPCp = 0.;
  fEPTPC = 0.; fEPBBC = 0.; fEPZDC = 0.;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  refmultCorr = 0x0;  // FIXME TEST
  refmult2Corr = 0x0; // FIXME TEST
  mOutName = outName;
  mOutNameEP = "";
  mOutNameQA = "";
  mOutNameME = "";
  doPrintEventCounter = kFALSE;
  fDoEffCorr = kFALSE;
  doEventPlaneRes = kFALSE;
  doTPCptassocBin = kFALSE;
  fTPCptAssocBin = -99;
  doUseMainEPAngle = kFALSE;  // kTRUE: use 0.2-2.0 GeV charged tracks for event plane
  fMinPtJet = minJetPt;
  fJetConstituentCut = 2.0;
  fTrackBias = trkbias;
  fTowerBias = 0.2;
  fJetRad = 0.4;
  fJetShapeTrackPtMin = 0.2;  fJetShapeTrackPtMax = 30.0;
  fJetShapePtAssocBin = 4;
  fLeadJetPtMin = 20.0; fSubLeadJetPtMin = 10.0;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0; fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2; fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0; fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0; fTowerPhiMaxCut = 2.0*TMath::Pi();
  fDoEventMixing = 0; fMixingTracks = 50000; fNMIXtracks = 5000; fNMIXevents = 5;
  fCentBinSize = 5; fCentBinSizeJS = 10; fReduceStatsCent = -1;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
  //mVertex = 0x0;
  zVtx = 0.0;
  fDoFilterPtMixEvents = kFALSE;
  fDoUseMultBins = kFALSE;
  doUseEPBins = kFALSE;
  doIgnoreExternalME = kTRUE;
  fEmcTriggerEventType = 0; fMBEventType = 2; fMixingEventType = 0;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4801; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }
  doComments = mDoComments;
  fhnJH = 0x0;
  fhnMixedEvents = 0x0;
  fhnCorr = 0x0;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  fEventPlaneMakerName = "";
  //fEventPoolOutputList = {{}}; // FIXME - does this make sense as is?
  fUsePtBinnedEventPool = kFALSE;
  fCheckEventNumberInMixedEvent = kFALSE;
  fListOfPools = 0x0;
}
//
//
//__________________________________________________________________________________
StMyAnalysisMaker3::~StMyAnalysisMaker3()
{ /*  */
  // destructor
  delete hdEPReactionPlaneFnc;
  delete hdEPEventPlaneFncN2;
  delete hdEPEventPlaneFncP2;
  delete hdEPEventPlaneFnc2;
  delete hdEPEventPlaneClass;
  delete hReactionPlaneFnc;
  delete hEventPlaneFncN2;
  delete hEventPlaneFncP2;
  delete hEventPlaneFnc2;
  delete hEventPlaneClass;

  delete hEventPlane;
  delete fHistEPTPCn;
  delete fHistEPTPCp;
  delete fHistEPBBC;
  delete fHistEPZDC;
  delete hEventZVertex;
  delete hCentrality;
  delete hMultiplicity;
  delete hRhovsCent;
  for(int i=0; i<5; i++) { delete hdEPtrk[i]; }
  for(int i=0; i<9; i++){ // centrality
    delete hTrackPhi[i];
    delete hTrackEta[i];
    delete hTrackPt[i];
  }
  delete hTrackEtavsPhi;

  delete hJetPt;
  delete hJetCorrPt;
  delete hJetLeadingPt;
  delete hJetSubLeadingPt;
  delete hJetLeadingPtAj;
  delete hJetSubLeadingPtAj;
  delete hJetDiJetAj;
  delete hJetE;
  delete hJetEta;
  delete hJetPhi;
  delete hJetNEF;
  delete hJetArea;
  delete hJetTracksPt;
  delete hJetTracksPhi;
  delete hJetTracksEta;
  delete hJetTracksZ;
  delete hJetPtvsArea;
  delete hJetEventEP;
  delete hJetPhivsEP;

  delete hJetPtIn;
  delete hJetPhiIn;
  delete hJetEtaIn;
  delete hJetEventEPIn;
  delete hJetPhivsEPIn;
  delete hJetPtMid;
  delete hJetPhiMid;
  delete hJetEtaMid;
  delete hJetEventEPMid;
  delete hJetPhivsEPMid;
  delete hJetPtOut;
  delete hJetPhiOut;
  delete hJetEtaOut;
  delete hJetEventEPOut;
  delete hJetPhivsEPOut;

  delete fHistJetHEtaPhi;
  delete fHistEventSelectionQA;
  delete fHistEventSelectionQAafterCuts;
  delete hTriggerIds;
  delete hEmcTriggers;
  delete hMixEvtStatZVtx;
  delete hMixEvtStatCent;
  delete hMixEvtStatZvsCent;
  delete hTriggerEvtStatZVtx;
  delete hTriggerEvtStatCent;
  delete hTriggerEvtStatZvsCent;
  delete hMBvsMult;
  delete hMB5vsMult;
  delete hMB30vsMult;
  delete hHTvsMult;
  delete hNMixEvents;

  if(hTPCvsBBCep) delete hTPCvsBBCep;
  if(hTPCvsZDCep) delete hTPCvsZDCep;
  if(hBBCvsZDCep) delete hBBCvsZDCep;

  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        delete fProfJetV2[k][j][i];
      }
    }
  }

  if(doJetShapeAnalysis) {
    for(int k=0; k<4; k++) {
      for(int j=0; j<4; j++) {
        for(int i=0; i<4; i++) {
          delete hJetCounter[k][j][i];
          delete hJetCounterCase1[k][j][i];
          delete hJetCounterCase2[k][j][i];
          delete hJetCounterCase3BG[k][j][i];

          for(int p=0; p<9; p++) {
            delete hJetShape[k][j][i][p];
            delete hJetShapeCase1[k][j][i][p];
            delete hJetShapeCase2[k][j][i][p];
            delete hJetShapeBG[k][j][i][p];
            delete hJetShapeBGCase1[k][j][i][p];
            delete hJetShapeBGCase2[k][j][i][p];
            delete hJetShapeBGCase3[k][j][i][p];

            delete hJetPtProfile[k][j][i][p];
            delete hJetPtProfileCase1[k][j][i][p];
            delete hJetPtProfileCase2[k][j][i][p];
            delete hJetPtProfileBG[k][j][i][p];
            delete hJetPtProfileBGCase1[k][j][i][p];
            delete hJetPtProfileBGCase2[k][j][i][p];
            delete hJetPtProfileBGCase3[k][j][i][p];
          }
        }
      }
    }
  }

  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      if(fProfV2Resolution[i]) delete fProfV2Resolution[i];
      if(fProfV3Resolution[i]) delete fProfV3Resolution[i];
      if(fProfV4Resolution[i]) delete fProfV4Resolution[i];
      if(fProfV5Resolution[i]) delete fProfV5Resolution[i];
    }
  }

  // delete sparses
  delete fhnJH;
  delete fhnMixedEvents;
  delete fhnCorr;

  // clear and delete objects
//  fJets->Clear(); delete fJets;
//  fRho->Clear(); delete fRho; 
//  fPoolMgr->Clear(); delete fPoolMgr;

  // Clear unnecessary pools before saving - FIXME
  fPoolMgr->ClearPools();
  delete fListOfPools;
}
//
// initialize objects & set up
//_________________________________________________________________________________________
Int_t StMyAnalysisMaker3::Init() {
  //StJetFrameworkPicoBase::Init();

  // Initialize output list of event pools - TEST
  if (fListOfPools != NULL){
        delete fListOfPools;
        fListOfPools = NULL;
  }
  if (!fListOfPools){
        fListOfPools = new TList();
        fListOfPools->SetOwner(kTRUE);
  }

  // initialize the histograms
  DeclareHistograms();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  // switch on Run Flag to look for firing trigger specifically requested for given run period
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        switch(fCentralityDef) {
          case StJetFrameworkPicoBase::kgrefmult :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P17id_VpdMB30();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P16id :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
              break;
          default: // this is the default for Run14
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
              refmult2Corr = CentralityMaker::instance()->getgRefMultCorr();
              refmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
        }
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
        break; // added May20

    case StJetFrameworkPicoBase::Run11_pp500 : // Run11: 500 GeV pp
        break;

    case StJetFrameworkPicoBase::Run12_pp200 : // Run12: 200 GeV pp
        break;

    case StJetFrameworkPicoBase::Run12_pp500 : // Run12: 500 GeV pp
        break;

    case StJetFrameworkPicoBase::Run13_pp510 : // Run13: 510 (500) GeV pp
        break;

    case StJetFrameworkPicoBase::Run15_pp200 : // Run15: 200 GeV pp
        break;

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17: 510 (500) GeV pp
        // this is the default for Run17 pp - don't set anything for pp
        break;

    default :
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  }

  return kStOK;
}
//
//
//______________________________________________________________________________________
Int_t StMyAnalysisMaker3::Finish() { 
  //  Summarize the run.
  cout << "StMyAnalysisMaker3::Finish()\n";
  //cout << "\tProcessed " << mEventCounter << " events." << endl;
  //cout << "\tWithout PV cuts: "<< mAllPVEventCounter << " events"<<endl;
  //cout << "\tInput events: "<<mInputEventCounter<<endl;

  // Write event pool manager object to file and close it
  if(mOutNameME!="") {
    TFile *fOutME = new TFile(mOutNameME.Data(), "RECREATE");
    fOutME->cd();
    fListOfPools->Write();  // write pools to file if we have them and want to
    fOutME->Write();
    fOutME->Close();
  }

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    //fout->ls();
    fout->cd();
    fout->mkdir(fAnalysisMakerName);
    fout->cd(fAnalysisMakerName);
    WriteHistograms();
   
    // jet shape analysis
    if(doJetShapeAnalysis) {
      for(int j=0; j<9; j++) {
        fout->cd();
        //fout->mkdir(Form("JetShapeAnalysis_bin%i", fTPCptAssocBin));
        //fout->cd(Form("JetShapeAnalysis_bin%i", fTPCptAssocBin));
        //fout->mkdir(Form("JetShapeAnalysis%i_bin%i", fJetShapeJetType, fJetShapePtAssocBin));
        //fout->cd(Form("JetShapeAnalysis%i_bin%i", fJetShapeJetType, fJetShapePtAssocBin));
        fout->mkdir(Form("JetShapeAnalysis_bin%i", j));
        fout->cd(Form("JetShapeAnalysis_bin%i", j));
        WriteJetShapeHistograms(j);
      }
    }

    fout->cd();
    fout->Write();
    fout->Close();
  }

  //  Write QA histos to file and close it.
  if(mOutNameQA!=""  && fJetShapePtAssocBin < 5) {
    TFile *fQAout = new TFile(mOutNameQA.Data(), "UPDATE");
    fQAout->cd();

    // track QA
    if(doWriteTrackQAHist && (fJetShapePtAssocBin < 5)) {
      fQAout->mkdir(Form("TrackQA_bin%i", fTPCptAssocBin));
      fQAout->cd(Form("TrackQA_bin%i", fTPCptAssocBin));
      WriteTrackQAHistograms();
      fQAout->cd();
    }

    // jet QA
    if(doWriteJetQAHist && (fJetShapePtAssocBin < 5)) {
      fQAout->mkdir(Form("JetEPQA_bin%i", fTPCptAssocBin));
      fQAout->cd(Form("JetEPQA_bin%i", fTPCptAssocBin));
      WriteJetEPQAHistograms();
      fQAout->cd();
    }

    fQAout->Write();
    fQAout->Close();
  }

  cout<<"End of StMyAnalysisMaker3::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}
//
// function to declare histograms and set up global objects
//__________________________________________________________________________________________
void StMyAnalysisMaker3::DeclareHistograms() {
  // constants
  double pi = 1.0*TMath::Pi();

  // binning for cent histograms
  int nHistCentBins = 0;
  if(fCentBinSize == 10) nHistCentBins = 10;
  if(fCentBinSize ==  5) nHistCentBins = 20;

  // binning for mult histograms
  double kHistMultMax = 800.;
  int kHistMultBins = 400;

  // pp specific settings
  if(doppAnalysis) {
    kHistMultMax = 100.;
    kHistMultBins = 100.;
  }

  // QA histos
  hdEPReactionPlaneFnc = new TH1F("hdEPReactionPlaneFnc", "jets relative to EP from reaction plane function", 3, 0.0, 0.5*pi);
  hdEPEventPlaneFncN2 = new TH1F("hdEPEventPlaneFncN2", "jets relative to EP from event plane function SUB1 meth2", 3, 0.0, 0.5*pi);
  hdEPEventPlaneFncP2 = new TH1F("hdEPEventPlaneFncP2", "jets relative to EP from event plane function SUB2 meth2", 3, 0.0, 0.5*pi);
  hdEPEventPlaneFnc2 = new TH1F("hdEPEventPlaneFnc2", "jets relative to EP from event plane function meth2", 3, 0.0, 0.5*pi);
  hdEPEventPlaneClass = new TH1F("hdEPEventPlaneClass", "jets relative to EP from event plane class", 3, 0.0, 0.5*pi);
  hReactionPlaneFnc = new TH1F("hReactionPlaneFnc", "Event plane distribution from reaction plane function", 72, 0.0, 1.0*pi);
  hEventPlaneFncN2 = new TH1F("hEventPlaneFncN2", "Event plane distribution from event plane function SUB1 meth2", 72, 0.0, 1.0*pi);
  hEventPlaneFncP2 = new TH1F("hEventPlaneFncP2", "Event plane distribution from event plane function SUB2 meth2", 72, 0.0, 1.0*pi);
  hEventPlaneFnc2 = new TH1F("hEventPlaneFnc2", "Event plane distribution from event plane function meth2", 72, 0.0, 1.0*pi);
  hEventPlaneClass = new TH1F("hEventPlaneClass", "Event plane distribution from event plane class", 72, 0.0, 1.0*pi);

  hEventPlane = new TH1F("hEventPlane", "Event plane distribution", 72, 0.0, 1.0*pi);
  fHistEPTPCn = new TH2F("fHistEPTPCn", "", 20, 0., 100., 72, -pi, pi);
  fHistEPTPCp = new TH2F("fHistEPTPCp", "", 20, 0., 100., 72, -pi, pi);
  fHistEPBBC = new TH2F("fHistEPBBC", "", 20, 0., 100., 72, -pi, pi);
  fHistEPZDC = new TH2F("fHistEPZDC", "", 20, 0., 100., 72, -pi, pi);
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 100, -50, 50);
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hRhovsCent = new TH2F("hRhovsCent", "#rho vs centrality", 20, 0, 100, 200, 0, 200);

  for(int i=0; i<5; i++) { // pt bins
    hdEPtrk[i] = new TH1F(Form("hdEPtrk%i", i), Form("tracks relative to event plane, p_{T} bin=%i", i), 72, 0, 0.5*pi);
  }

  // track phi distribution for centrality
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i] = new TH1F(Form("hTrackPhi%d", i), Form("track distribution vs #phi, centr%d", i), 144, 0, 2.0*pi);
    hTrackEta[i] = new TH1F(Form("hTrackEta%d", i), Form("track distribution vs #eta, centr%d", i), 40, -1.0, 1.0);
    hTrackPt[i] = new TH1F(Form("hTrackPt%d", i), Form("track distribution vs p_{T}, centr%d", i), 120, 0., 30.0);
  }
  hTrackEtavsPhi = new TH2F(Form("hTrackEtavsPhi"), Form("track distribution: #eta vs #phi"), 144, 0, 2.0*pi, 40, -1.0, 1.0);

  // jet QA histos
  hJetPt = new TH1F("hJetPt", "Jet p_{T}", 100, 0, 100);
  hJetCorrPt = new TH1F("hJetCorrPt", "Corrected Jet p_{T}", 125, -25, 100);
  hJetLeadingPt = new TH1F("hJetLeadingPt", "Leading Jet p_{T}", 100, 0, 100);
  hJetSubLeadingPt = new TH1F("hJetSubLeadingPt", "SubLeading Jet p_{T}", 100, 0, 100);
  hJetLeadingPtAj = new TH1F("hJetLeadingPtAj", "Leading Jet p_{T} with Aj cut", 100, 0, 100);
  hJetSubLeadingPtAj = new TH1F("hJetSubLeadingPtAj", "SubLeading Jet p_{T} with Aj cut", 100, 0, 100);
  hJetDiJetAj = new TH1F("hJetDiJetAj", "DiJet imbalance: Aj", 100, 0., 1.);
  hJetE = new TH1F("hJetE", "Jet energy distribution", 100, 0, 100);
  hJetEta = new TH1F("hJetEta", "Jet #eta distribution", 24, -1.2, 1.2);
  hJetPhi = new TH1F("hJetPhi", "Jet #phi distribution", 72, 0.0, 2.0*pi);
  hJetNEF = new TH1F("hJetNEF", "Jet NEF", 100, 0, 1);
  hJetArea = new TH1F("hJetArea", "Jet Area", 100, 0, 1);
  hJetTracksPt = new TH1F("hJetTracksPt", "Jet track constituent p_{T}", 120, 0, 30.0);
  hJetTracksPhi = new TH1F("hJetTracksPhi", "Jet track constituent #phi", 72, 0, 2.0*pi);
  hJetTracksEta = new TH1F("hJetTracksEta", "Jet track constituent #eta", 56, -1.4, 1.4);
  hJetTracksZ = new TH1F("hJetTracksZ", "Jet track fragmentation function", 144, 0, 1.44);
  hJetPtvsArea = new TH2F("hJetPtvsArea", "Jet p_{T} vs Jet area", 100, 0, 100, 100, 0, 1);
  hJetEventEP = new TH1F("hJetEventEP", "no of jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEP = new TH2F("hJetPhivsEP", "Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);

  // jet EP QA histos
  hJetPtIn = new TH1F("hJetPtIn", "no of jets in-plane vs p_{T}", 100, 0.0, 100);
  hJetPhiIn = new TH1F("hJetPhiIn", "no of jets in-plane vs #phi", 72, 0.0, 2.0*pi);
  hJetEtaIn = new TH1F("hJetEtaIn", "no of jets in-plane vs #eta", 40, -1.0, 1.0);
  hJetEventEPIn = new TH1F("hJetEventEPIn", "no of in-plane jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEPIn = new TH2F("hJetPhivsEPIn", "in-plane Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
  hJetPtMid = new TH1F("hJetPtMid", "no of jets mid-plane vs p_{T}", 100, 0.0, 100);
  hJetPhiMid = new TH1F("hJetPhiMid", "no of jets mid-plane vs #phi", 72, 0.0, 2.0*pi);
  hJetEtaMid = new TH1F("hJetEtaMid", "no of jets mid-plane vs #eta", 40, -1.0, 1.0);
  hJetEventEPMid = new TH1F("hJetEventEPMid", "no of mid-plane jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEPMid = new TH2F("hJetPhivsEPMid", "mid-plane Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
  hJetPtOut = new TH1F("hJetPtOut", "no of jets out-of-plane vs p_{T}", 100, 0.0, 100);
  hJetPhiOut = new TH1F("hJetPhiOut", "no of jets out-of-plane vs #phi", 72, 0.0, 2.0*pi);
  hJetEtaOut = new TH1F("hJetEtaOut", "no of jets out-of-plane vs #eta", 40, -1.0, 1.0);
  hJetEventEPOut = new TH1F("hJetEventEPOut", "no of out-of-plane jet events vs event plane", 72, 0.0, 1.0*pi);
  hJetPhivsEPOut = new TH2F("hJetPhivsEPOut", "out-of-plane Jet #phi vs event plane", 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);

  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi", "Jet-hadron #Delta#eta-#Delta#phi", 72, -1.8, 1.8, 72, -0.5*pi, 1.5*pi);

  // Event Selection QA histo + event QA
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);
  hMixEvtStatZVtx = new TH1F("hMixEvtStatZVtx", "no of events in pool vs zvtx", 20, -40.0, 40.0);
  hMixEvtStatCent = new TH1F("hMixEvtStatCent", "no of events in pool vs Centrality", nHistCentBins, 0, 100);
  hMixEvtStatZvsCent = new TH2F("hMixEvtStatZvsCent", "no of events: zvtx vs Centality", nHistCentBins, 0, 100, 20, -40.0, 40.0);
  hTriggerEvtStatZVtx = new TH1F("hTriggerEvtStatZVtx", "no of trigger events used vs zvtx", 20, -40.0, 40.0);
  hTriggerEvtStatCent = new TH1F("hTriggerEvtStatCent", "no of trigger events used vs Centrality", nHistCentBins, 0, 100);
  hTriggerEvtStatZvsCent = new TH2F("hTriggerEvtStatZvsCent", "no of trigger events used: zvtx vs Centality", nHistCentBins, 0, 100, 20, -40.0, 40.0);
  hMBvsMult = new TH1F("hMBvsMult", "# MB events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hMB5vsMult = new TH1F("hMB5vsMult", "# MB5 events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hMB30vsMult = new TH1F("hMB30vsMult", "# MB30 events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hHTvsMult = new TH1F("hHTvsMult", "# HT events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  hNMixEvents = new TH1F("hNMixEvents", "number of mixing events", 200, 0, 200);

  //// res_cen=new TProfile("res_cen","res vs. cen",10,0,10,-2,2);
  // set binning for run based corrections - run dependent
  //Int_t nRunBins = 1; // - just a default
  //if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) nRunBins = 830; //1654;
  //if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) nRunBins = 1359;

  // 2-D event plane differences - (updated ranges on February 7th)
  hTPCvsBBCep = new TH2F("hTPCvsBBCep", "TPC vs BBC 2nd order event plane", 144, 0.*pi, 1.0*pi, 144, 0.*pi, 1.*pi);
  hTPCvsZDCep = new TH2F("hTPCvsZDCep", "TPC vs ZDC 2nd order event plane", 144, 0.*pi, 1.0*pi, 144, 0.*pi, 1.*pi);
  hBBCvsZDCep = new TH2F("hBBCvsZDCep", "BBC vs ZDC 2nd order event plane", 144, 0.*pi, 1.0*pi, 144, 0.*pi, 1.*pi);

  // jet v2 histogram - primarily used for EP resolution unfolding
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        fProfJetV2[k][j][i] = new TProfile(Form("fProfJetV2_%i_%i_%i", k, j, i), Form("fProfJetV2_%i_%i_%i", k, j, i), 11, -0.5, 10.5);
      }
    }
  }

  // jet shape analysis histograms
  if(doJetShapeAnalysis) {
    for(int k=0; k<4; k++) {
      for(int j=0; j<4; j++) {
        for(int i=0; i<4; i++) {
          hJetCounter[k][j][i] = new TH1F(Form("hJetCounter_%i_%i_%i", k, j, i), Form("Jet counter - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);
          hJetCounterCase1[k][j][i] = new TH1F(Form("hJetCounterCase1_%i_%i_%i", k, j, i), Form("Jet counter case2 - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);
          hJetCounterCase2[k][j][i] = new TH1F(Form("hJetCounterCase2_%i_%i_%i", k, j, i), Form("Jet counter case2 - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);
          hJetCounterCase3BG[k][j][i] = new TH1F(Form("hJetCounterCase3BG_%i_%i_%i", k, j, i), Form("Jet counter case3 BG only - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);

          for(int p=0; p<9; p++) {
            hJetShape[k][j][i][p] = new TH1F(Form("hJetShape_%i_%i_%i_%i", k, j, i, p), Form("Jet shape #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetShapeCase1[k][j][i][p] = new TH1F(Form("hJetShapeCase1_%i_%i_%i_%i", k, j, i, p), Form("Jet shape case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetShapeCase2[k][j][i][p] = new TH1F(Form("hJetShapeCase2_%i_%i_%i_%i", k, j, i, p), Form("Jet shape case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetShapeBG[k][j][i][p] = new TH1F(Form("hJetShapeBG_%i_%i_%i_%i", k, j, i, p), Form("Jet shape BG #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetShapeBGCase1[k][j][i][p] = new TH1F(Form("hJetShapeBGCase1_%i_%i_%i_%i", k, j, i, p), Form("Jet shape BG case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetShapeBGCase2[k][j][i][p] = new TH1F(Form("hJetShapeBGCase2_%i_%i_%i_%i", k, j, i, p), Form("Jet shape BG case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetShapeBGCase3[k][j][i][p] = new TH1F(Form("hJetShapeBGCase3_%i_%i_%i_%i", k, j, i, p), Form("Jet shape BG case3 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);

            hJetPtProfile[k][j][i][p] = new TH1F(Form("hJetPtProfile_%i_%i_%i_%i", k, j, i, p), Form("Jet pt profile #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetPtProfileCase1[k][j][i][p] = new TH1F(Form("hJetPtProfileCase1_%i_%i_%i_%i", k, j, i, p), Form("Jet pt profile case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetPtProfileCase2[k][j][i][p] = new TH1F(Form("hJetPtProfileCase2_%i_%i_%i_%i", k, j, i, p), Form("Jet pt profile case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetPtProfileBG[k][j][i][p] = new TH1F(Form("hJetPtProfileBG_%i_%i_%i_%i", k, j, i, p), Form("Jet pt profile BG #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetPtProfileBGCase1[k][j][i][p] = new TH1F(Form("hJetPtProfileBGCase1_%i_%i_%i_%i", k, j, i, p), Form("Jet profile BG case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetPtProfileBGCase2[k][j][i][p] = new TH1F(Form("hJetPtProfileBGCase2_%i_%i_%i_%i", k, j, i, p), Form("Jet profile BG case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
            hJetPtProfileBGCase3[k][j][i][p] = new TH1F(Form("hJetPtProfileBGCase3_%i_%i_%i_%i", k, j, i, p), Form("Jet profile BG case3 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i, associated p_{T} bin %i", k, j, i, p), 10, 0.0, 0.50);
          }
        }
      }
    }
  }

  if(doEventPlaneRes){
    // Reaction Plane resolution as function of centrality - corrected for 2nd order event plane
    for (Int_t i=0; i<9; i++){
      // 2nd order correction
      fProfV2Resolution[i] = new TProfile(Form("fProfV2Resolution_%i", i), Form("fProfV2Resolution_%i", i), 25, 0.5, 25.5);
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(2(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(2(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(2(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(2(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(2(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(2(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(2(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(2(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(2(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(2(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(2(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(2(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(2(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(2(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(2(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(2(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 3rd order correction
      fProfV3Resolution[i] = new TProfile(Form("fProfV3Resolution_%i", i), Form("fProfV3Resolution_%i", i), 25, 0.5, 25.5);
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(3(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(3(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(3(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(3(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(3(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(3(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(3(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(3(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(3(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(3(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(3(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(3(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(3(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(3(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(3(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(3(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(3(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(3(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 4th order correction
      fProfV4Resolution[i] = new TProfile(Form("fProfV4Resolution_%i", i), Form("fProfV4Resolution_%i", i), 25, 0.5, 25.5);
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(4(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(4(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(4(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(4(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(4(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(4(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(4(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(4(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(4(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(4(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(4(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(4(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(4(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(4(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(4(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(4(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(4(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(4(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 5th order correction
      fProfV5Resolution[i] = new TProfile(Form("fProfV5Resolution_%i", i), Form("fProfV5Resolution_%i", i), 25, 0.5, 25.5);
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(5(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(5(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(5(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(5(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(5(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(5(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(5(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(5(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(5(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(5(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(5(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(5(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(5(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(5(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(5(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(5(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(5(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(5(#Psi_{ZDC1} - #Psi_{TPCp}))>");
    }
  }

  // set up jet-hadron sparse
  UInt_t bitcodeMESE = 0; // bit coded, see GetDimParams() below
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; // | 1<<8 | 1<<9 | 1<<10;
  fhnJH = NewTHnSparseF("fhnJH", bitcodeMESE);

  // =================================================================================================
  // set up centrality bins for mixed events
  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it

  // Setup for Au-Au collisions: cent bin size can only be 5 or 10% bins
  // centrality bins for mixed events
  int nCentralityBinsJS = 100 + 1;
  double multJS = 1.0;
  if(fCentBinSizeJS==5){ // will be most commonly used
    nCentralityBinsJS = 20;
    multJS = 5.0;
  } else if(fCentBinSizeJS==10){
    nCentralityBinsJS = 10;
    multJS = 10.0;
  } else if(fCentBinSizeJS==20){
    nCentralityBinsJS = 5;
    multJS = 20.0;
  }

  // set bin edges for centrality mixed events
  Double_t centralityBinsJS[nCentralityBinsJS +1];
  for(Int_t ic = 0; ic < nCentralityBinsJS + 1; ic++){
    //centralityBinsJS[ic] = multJS*ic; 
    centralityBinsJS[ic] = 1.0*ic;
  }

  // this is temp as the above and various other implementation attempts would not work for both cases
  // need to look into this, but a few hours already is too much.  Don't want to have to have this hard-coded
  Int_t nCentBins = 8;
  Double_t* centralityBins = GenerateFixedBinArray(nCentBins, 0., 8.);

/*
Centrality_def_grefmult.txt               - OLD (16)
 10 14 21 29 40 54 | 71 92 116 145 179 218 | 263 315 | 373 441  
Centrality_def_grefmult_P17id_VpdMB30.txt - NEW (16)
 10 15 22 31 43 58 | 76 97 123 154 189 230 | 276 329 | 390 459  
*/

  // for AuAu data
//  Int_t nMultBinsJS = 29; 
//  Double_t multBinsJS[] = {10, 14, 19, 25, 31, 37, 44, 52, 61, 71, 82, 95, 109, 124, 140, 157, 175, 194, 214, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469};
//  Int_t nMultBinsJS = 24;  // Alt-1: Best Yet
//  Double_t multBinsJS[] = {10,16,24,34,46,61,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469};
  ///Int_t nMultBinsJS = 25;  // Alt-2
  ///Double_t multBinsJS[] = {10,15,21,31,42,53,66,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469};

  // for AuAu data
  // +1 to accomodate the fact that we define bins rather than array entries
  const int nMultBinsJS = 26;  // Alt-2 - 27 values, 26 ranges
  Double_t multBinsJS[nMultBinsJS + 1] = {0, 10,15,21,31,42,53,66,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469, 800};
  Double_t *multiplicityBinsJS = multBinsJS;

  // cent bins for AuAu data
  Int_t nCentBinsJS = 20;
  Double_t* centralityBinsJSnew = GenerateFixedBinArray(nCentBinsJS, 0., 20.); 

  // for pp data
  const int nMultBinsJSpp = 7;
  //Double_t multBinsJSpp[nMultBinsJSpp + 1] = {0.0, 4., 9, 15, 25, 35, 55, 100.0, 500.0};  // 8 (9)
  Double_t multBinsJSpp[nMultBinsJSpp + 1] = {0.0, 4.0, 6.0, 8.0, 10.0, 13.0, 30., 100.};   // 7 (8)
  Double_t *multiplicityBinsJSpp = multBinsJSpp;

  // z-vertex bins for mixed events
  Int_t nZvBins  = 20; // 4 cm wide, 40 for 2 cm wide
  Double_t* zvBins = GenerateFixedBinArray(nZvBins, -40., 40.); // min/max doesn't matter as data is cut zmin/zmax

  // event plane bins for mixed events
  Int_t nEPBins = 6; // 6 from 0-180 degrees, (0 - pi)
  Double_t* epBins = GenerateFixedBinArray(nEPBins, 0., 1.0*pi);

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implementation of AliEventPoolManager
  if(doIgnoreExternalME) {
    if(doJetShapeAnalysis) {
      if(fDoUseMultBins) {
        if(!doppAnalysis && !doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins); // not pp
        if(!doppAnalysis &&  doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins); // not pp
        if( doppAnalysis) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJSpp, (Double_t*)multiplicityBinsJSpp, nZvBins, (Double_t*)zvBins); // is pp

      } else { // centrality binning  
        //if(!doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinsJS, (Double_t*)centralityBinsJS, nZvBins, (Double_t*)zvBins);
        //if( doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinsJS, (Double_t*)centralityBinsJS, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins);
        if(!doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsJS, (Double_t*)centralityBinsJSnew, nZvBins, (Double_t*)zvBins);
        if( doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsJS, (Double_t*)centralityBinsJSnew, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins);

      } // correlation analysis setup

    } else { fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBins, (Double_t*)centralityBins, nZvBins, (Double_t*)zvBins); }
  }

  // ================================================================================================
  // TEMP testing 
  if(!doIgnoreExternalME) {
    Int_t nPtBins = 1;
    Double_t defaultPtBins[2] = {-9999., 9999.};
    Double_t* ptBins = defaultPtBins;

    if(!fPoolMgr) {
      cout<<"Don't have pre-existing fPoolMgr object, We will now create it!"<<endl;
      fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBins, (Double_t*)zvBins, nEPBins, epBins, nPtBins, ptBins);
    } else {
      cout<<"User has provided fPoolMgr object.. goodluck!"<<endl;
    }

    // If some bins of the pool should be saved, fEventPoolOutputList must be given using AddEventPoolsToOutput()
    // Note that this is in principle also possible, if an external poolmanager was given
    for(UInt_t i = 0; i < fEventPoolOutputList.size(); i++) {
      Double_t minCent = fEventPoolOutputList[i][0];
      Double_t maxCent = fEventPoolOutputList[i][1];
      Double_t minZvtx = fEventPoolOutputList[i][2];
      Double_t maxZvtx = fEventPoolOutputList[i][3];
      Double_t minPsi2 = fEventPoolOutputList[i][4];
      Double_t maxPsi2 = fEventPoolOutputList[i][5];
      Double_t minPt   = fEventPoolOutputList[i][6];
      Double_t maxPt   = fEventPoolOutputList[i][7];

      fPoolMgr->SetSaveFlag(minCent, maxCent, minZvtx, maxZvtx, minPsi2, maxPsi2, minPt, maxPt);
    }

    // Basic checks and printing of pool properties
    fPoolMgr->Validate();

    // save to output if requested
    if(fEventPoolOutputList.size()) fListOfPools->Add(fPoolMgr);

  } // external use switch
  // ================================================================================================

  // set up event mixing sparse
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; // | 1<<8 | 1<<9;
  fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", bitcodeMESE);

  // jet counter for normalizations in correlation analysis
  UInt_t bitcodeCorr = 0; // bit coded, see GetDimparamsCorr() below
  bitcodeCorr = 1<<0 | 1<<1 | 1<<2 | 1<<3; // | 1<<4;
  fhnCorr = NewTHnSparseFCorr("fhnCorr", bitcodeCorr);

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}
//
// write track QA histograms
//_____________________________________________________________________________
void StMyAnalysisMaker3::WriteTrackQAHistograms() {
  // track phi distribution for centrality
  for(int i=0; i<5; i++) { hdEPtrk[i]->Write(); }
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i]->Write();
    hTrackEta[i]->Write();
    hTrackPt[i]->Write();
  }
  hTrackEtavsPhi->Write();
}
//
// write Jet event plane QA histograms
//______________________________________________________________________________
void StMyAnalysisMaker3::WriteJetEPQAHistograms() {
  hdEPReactionPlaneFnc->Write();
  hdEPEventPlaneFncN2->Write();
  hdEPEventPlaneFncP2->Write();
  hdEPEventPlaneFnc2->Write();
  hdEPEventPlaneClass->Write();
  hReactionPlaneFnc->Write();
  hEventPlaneFncN2->Write();
  hEventPlaneFncP2->Write();
  hEventPlaneFnc2->Write();
  hEventPlaneClass->Write();

  hJetPtIn->Write();
  hJetPhiIn->Write();
  hJetEtaIn->Write();
  hJetEventEPIn->Write();
  hJetPhivsEPIn->Write();
  hJetPtMid->Write();
  hJetPhiMid->Write();
  hJetEtaMid->Write();
  hJetEventEPMid->Write();
  hJetPhivsEPMid->Write();
  hJetPtOut->Write();
  hJetPhiOut->Write();
  hJetEtaOut->Write();
  hJetEventEPOut->Write();
  hJetPhivsEPOut->Write();
}
//
// write jet shape histograms
//________________________________________________________________________________
void StMyAnalysisMaker3::WriteJetShapeHistograms(Int_t option) {
  if(option >= 0) {
    // jet shape histos
    if(doJetShapeAnalysis) {
      for(int k=0; k<4; k++) {
        for(int j=0; j<4; j++) {
          for(int i=0; i<4; i++) {
            hJetCounter[k][j][i]->Write();
            hJetCounterCase1[k][j][i]->Write();
            hJetCounterCase2[k][j][i]->Write();
            hJetCounterCase3BG[k][j][i]->Write();

            hJetShape[k][j][i][option]->Write();
            hJetShapeCase1[k][j][i][option]->Write();
            hJetShapeCase2[k][j][i][option]->Write();
            hJetShapeBG[k][j][i][option]->Write();
            hJetShapeBGCase1[k][j][i][option]->Write();
            hJetShapeBGCase2[k][j][i][option]->Write();
            hJetShapeBGCase3[k][j][i][option]->Write();

            hJetPtProfile[k][j][i][option]->Write();
            hJetPtProfileCase1[k][j][i][option]->Write();
            hJetPtProfileCase2[k][j][i][option]->Write();
            hJetPtProfileBG[k][j][i][option]->Write();
            hJetPtProfileBGCase1[k][j][i][option]->Write();
            hJetPtProfileBGCase2[k][j][i][option]->Write();
            hJetPtProfileBGCase3[k][j][i][option]->Write();
          }
        }
      }
    }

  } else {
    // jet shape histos
    if(doJetShapeAnalysis) {
      for(int k=0; k<4; k++) {
        for(int j=0; j<4; j++) {
          for(int i=0; i<4; i++) {
            hJetCounter[k][j][i]->Write();
            hJetCounterCase1[k][j][i]->Write();
            hJetCounterCase2[k][j][i]->Write();
            hJetCounterCase3BG[k][j][i]->Write();

            for(int p=0; p<9; p++) {
              hJetShape[k][j][i][p]->Write();
              hJetShapeCase1[k][j][i][p]->Write();
              hJetShapeCase2[k][j][i][p]->Write();
              hJetShapeBG[k][j][i][p]->Write();
              hJetShapeBGCase1[k][j][i][p]->Write();
              hJetShapeBGCase2[k][j][i][p]->Write();
              hJetShapeBGCase3[k][j][i][p]->Write();

              hJetPtProfile[k][j][i][p]->Write();
              hJetPtProfileCase1[k][j][i][p]->Write();
              hJetPtProfileCase2[k][j][i][p]->Write();
              hJetPtProfileBG[k][j][i][p]->Write();
              hJetPtProfileBGCase1[k][j][i][p]->Write();
              hJetPtProfileBGCase2[k][j][i][p]->Write();
              hJetPtProfileBGCase3[k][j][i][p]->Write();
            }
          }
        }
      }
    }
  }

}
//
// write histograms
//_____________________________________________________________________________
void StMyAnalysisMaker3::WriteHistograms() {
  // default histos
  hEventPlane->Write();
  fHistEPTPCn->Write();
  fHistEPTPCp->Write();
  fHistEPBBC->Write();
  fHistEPZDC->Write();
  hEventZVertex->Write();
  hCentrality->Write();
  hMultiplicity->Write();
  hRhovsCent->Write();

  // jet histos
  hJetPt->Write();
  hJetCorrPt->Write();
  hJetLeadingPt->Write();
  hJetSubLeadingPt->Write();
  hJetLeadingPtAj->Write();
  hJetSubLeadingPtAj->Write();
  hJetDiJetAj->Write();
  hJetE->Write();
  hJetEta->Write();
  hJetPhi->Write();
  hJetNEF->Write();
  hJetArea->Write();
  hJetTracksPt->Write();
  hJetTracksPhi->Write();
  hJetTracksEta->Write();
  hJetTracksZ->Write();
  hJetPtvsArea->Write();
  hJetEventEP->Write();
  hJetPhivsEP->Write();

  fHistJetHEtaPhi->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
  hMixEvtStatZVtx->Write();
  hMixEvtStatCent->Write();
  hMixEvtStatZvsCent->Write();
  hTriggerEvtStatZVtx->Write();
  hTriggerEvtStatCent->Write();
  hTriggerEvtStatZvsCent->Write();
  hMBvsMult->Write();
  hMB5vsMult->Write();
  hMB30vsMult->Write();
  hHTvsMult->Write();
  hNMixEvents->Write();

  // jet sparse
  if(!doJetShapeAnalysis) { // don't write these when doing jet shape analysis
    fhnJH->Write();
    fhnMixedEvents->Write();
    fhnCorr->Write();
  }

  // (perhaps temp - save resolution hists to main output file instead of event plane calibration file)
  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      fProfV2Resolution[i]->Write();
      fProfV3Resolution[i]->Write();
      fProfV4Resolution[i]->Write();
      fProfV5Resolution[i]->Write();
    }
  }

  // TODO - may decide to move elsewhere
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        fProfJetV2[k][j][i]->Write();
      }
    }
  }

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StMyAnalysisMaker3::Clear(Option_t *opt) {
  fJets->Clear();
}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StMyAnalysisMaker3::Make() {
  const double pi = 1.0*TMath::Pi();
  //StMemStat::PrintMem("MyAnalysisMaker at beginning of make");

  // update counter
//  mEventCounter++;
//  cout<<"StMyANMaker event# = "<<mEventCounter<<endl;
  if(doPrintEventCounter) cout<<"StMyAnMaker event# = "<<EventCounter()<<endl;

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

  // cut event on max track pt > 35.0 GeV (30 Oct25, 2018)
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  // cut on (-30, 30) when using NEW centrality definitions - perhaps cut on (-28, 28)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  hEventZVertex->Fill(zVtx);

  // get the Run #, fill, and event ID
  int RunId = mPicoEvent->runId();
  fRunNumber = mPicoEvent->runId();
  int fillId = mPicoEvent->fillId();
  int eventId = mPicoEvent->eventId();
  double fBBCCoincidenceRate = mPicoEvent->BBCx();
  double fZDCCoincidenceRate = mPicoEvent->ZDCx();
  if(fDebugLevel == kDebugGeneralEvt) cout<<"RunID = "<<RunId<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

  // ============================ CENTRALITY ============================== //
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_refmult.txt
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_grefmult.txt
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  int grefMult = mPicoEvent->grefMult();
  //int refMult = mPicoEvent->refMult();
  Int_t centbin, cent9, cent16;
  Double_t refCorr2;

  // check for AuAu analyses
  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorr->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }
//    if(grefmultCorr->isBadRun(RunId)) cout << "Run is bad" << endl; 
//    if(grefmultCorr->isIndexOk()) cout << "Index Ok" << endl;
//    if(grefmultCorr->isZvertexOk()) cout << "Zvertex Ok" << endl;
//    if(grefmultCorr->isRefMultOk()) cout << "RefMult Ok" << endl;

    // get centrality bin: either 0-7 or 0-15
    cent16 = grefmultCorr->getCentralityBin16();
    cent9 = grefmultCorr->getCentralityBin9();

    // re-order binning to be from central -> peripheral
    ref9 = GetCentBin(cent9, 9);
    ref16 = GetCentBin(cent16, 16);  
    centbin = GetCentBin(cent16, 16);  // 0-16

    // calculate corrected multiplicity
    // refCorr2 = grefmultCorr->getRefMultCorr();  // already initialized

    // calculate corrected multiplicity: 
    // Double_t getRefMultCorr(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag=2) const ;
    // flag=0:  Luminosity only
    // flag=1:  z-vertex only
    // flag=2:  full correction (default)
    if(doUseBBCCoincidenceRate) { refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fZDCCoincidenceRate, 2); }

    //grefmultCorr->isCentralityOk(cent16)
  } else { // for pp
    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;
  }

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStWarn; // maybe kStOk; - this is for lowest multiplicity events 80%+ centrality, cut on them
 
  // bin-age to use for mixed event and sparses
  Int_t centbin10 = GetCentBin10(centbin);
  double centBinToUse;
  if(fCentBinSize==10)       { centBinToUse = (double)centbin10 * 10.0;
  } else if(fCentBinSize==5) { centBinToUse = (double)centbin * 5.0; }

  // centrality / multiplicity histograms
  // event activity - compensate for pp or AuAu
  double kEventActivity = (doppAnalysis) ? (double)grefMult : refCorr2;
  hMultiplicity->Fill(kEventActivity);
  if(fDebugLevel == kDebugCentrality) { if(centbin > 15) cout<<"centbin = "<<centbin<<"  mult = "<<refCorr2<<"  Centbin*5.0 = "<<centbin*5.0<<"  cent16 = "<<cent16<<endl; }
  fCentralityScaled = centbin*5.0;
  hCentrality->Fill(fCentralityScaled);

  // to limit filling unused entries in sparse, only fill for certain centrality ranges
  // ranges can be different than functional cent bin setter
  Int_t cbin = -1;
  // this is actually not used since the below line does the cut:  if(fRequireCentSelection)
  if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
  else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
  else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
  else if (centbin>9 && centbin<16) cbin = 5; // 50-80%
  else cbin = -99;

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

  // switches for Jet and Event Plane analysis
  Bool_t doJetAnalysis = kFALSE; // set false by default
  Bool_t doEPAnalysis = kFALSE;  // set false by default

  // if we have trigger: perform jet analysis
  if(fHaveEmcTrigger) { doJetAnalysis = kTRUE; }

  // if we have trigger && AuAu dataset: run event plane analysis
  if(fHaveEmcTrigger && (fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200 || fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200)) {
    doEPAnalysis = kTRUE;
  }

  // ======================== end of Triggers ============================= //

  // ================= JetMaker ================ //
  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get JetMaker collection of jets
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

/*
  // February 7th, 2019 - FIXME, fix this chunk!
  // check if bad/dead towers fired trigger and kill event if true
  badTowers = JetMaker->GetBadTowers();
  deadTowers = JetMaker->GetDeadTowers();
  if(DidBadTowerFireTrigger()) { 
    cout<<"FOUND PROBLEM.. a bad tower fired the trigger: "<<GetName()<<endl;  
    return kStWarn;
  }
*/

/*
  // TEST: the below is snippet of code for getting jets and their cosntituents using
  // constituent subtractor method in StJetMakerTask
  //

  // get JetMaker collection of jets with background subtraction
  TClonesArray *fJetsBGsub = static_cast<TClonesArray*>(JetMaker->GetJetsBGsub());
  if(!fJetsBGsub) {
    LOG_WARN << Form(" No fJetsBGsub object! Skip! ") << endm;
    return kStWarn;
  }

  // loop over background subtracted jets
  Int_t njetsbgsub = fJetsBGsub->GetEntries();
  for(int ijet = 0; ijet < njetsbgsub; ijet++) {  // JET LOOP
    // get pointer to jets
    StJet *jet = static_cast<StJet*>(fJetsBGsub->At(ijet));
    if(!jet) continue;
  
    vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();
    for(UInt_t ic = 0; ic < fConstituents.size(); ++ic) {
      // get user defined index
      Int_t uid = fConstituents[ic].user_index();
      double cpt = fConstituents[ic].perp();
      double ceta = fConstituents[ic].eta();
      double cphi = fConstituents[ic].phi();
      cout<<"ic = "<<ic<<", uid = "<<uid<<", cpt = "<<cpt<<", ceta = "<<ceta<<", cphi = "<<cphi<<endl;

    } // loop over constituents

  } // loop over jets

*/

  // ================================================================================

  // ============== RhoMaker =============== //
  // get RhoMaker from event: old names "StRho_JetsBG", "OutRho", "StMaker#0"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
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

  // get rho/area value from rho object     fRho->ls("");
  //double value = GetRhoValue(fRhoMakerName);
  fRhoVal = fRho->GetVal();
  hRhovsCent->Fill(centbin*5.0, fRhoVal);
  if(fDebugLevel == kDebugEmcTrigger) cout<<"   fRhoVal = "<<fRhoVal<<"   Correction = "<<1.0*TMath::Pi()*fJetRad*fJetRad*fRhoVal<<endl;

  // =========== Leading and Subleading Jets ============= //
  // cache the leading + subleading jets within acceptance
  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  // fill leading and subleading pt spectra
  if(fLeadingJet)    hJetLeadingPt->Fill(fLeadingJet->Pt());
  if(fSubLeadingJet) hJetSubLeadingPt->Fill(fSubLeadingJet->Pt());   

  // get dijet imbalance Aj -  z = if(condition) then(?) <do this> else(:) <do this>  
  double fDiJetAj = (fLeadingJet && fSubLeadingJet) ? GetDiJetAj(fLeadingJet, fSubLeadingJet, fRho, fCorrJetPt) : -999.;
  if(fDiJetAj > 0)  hJetDiJetAj->Fill(fDiJetAj);

  // =========== Event Plane Angle ============= //
  // basic method - not used for anything..
  double rpAngle = GetReactionPlane();
  hEventPlane->Fill(rpAngle);

  // ============== EventPlaneMaker =============== //
  // get StEventPlaneMaker from event
  double angle2, angle2n, angle2p, angle2comb, angle1, angle4;
  StEventPlaneMaker *EventPlaneMaker0, *EventPlaneMaker1, *EventPlaneMaker2, *EventPlaneMaker3, *EventPlaneMaker4;
  double tpc2EP_bin0, tpc2EP_bin1, tpc2EP_bin2, tpc2EP_bin3, tpc2EP_bin4;
  if(doEPAnalysis) {
    if(!doTPCptassocBin) {
      EventPlaneMaker = static_cast<StEventPlaneMaker*>(GetMaker(fEventPlaneMakerName));
      const char *fEventPlaneMakerNameCh = fEventPlaneMakerName;
      if(!EventPlaneMaker) {
        LOG_WARN << Form(" No %s! Skip! ", fEventPlaneMakerNameCh) << endm;
        return kStWarn;
      }

      // set event plane
      double tpc2EP = (double)EventPlaneMaker->GetTPCEP();
      TPC_PSI2 = tpc2EP;

    } else { // pt-dependent bin mode
      const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
      EventPlaneMaker0 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 0)));
      EventPlaneMaker1 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 1)));
      EventPlaneMaker2 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 2)));
      EventPlaneMaker3 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 3)));
      EventPlaneMaker4 = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 4)));

      // check for requested EventPlaneMaker pointer
      if((fTPCptAssocBin == 0) && (!EventPlaneMaker0)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 1) && (!EventPlaneMaker1)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 2) && (!EventPlaneMaker2)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 3) && (!EventPlaneMaker3)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}
      if((fTPCptAssocBin == 4) && (!EventPlaneMaker4)) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}

      // get event plane angle for different pt bins
      // could also write this as:  tpc2EP_bin = (EventPlaneMaker) ? (double)EventPlaneMaker->GetTPCEP() : -999;
      if(EventPlaneMaker0) { tpc2EP_bin0 = (double)EventPlaneMaker0->GetTPCEP(); } else { tpc2EP_bin0 = -999; }
      if(EventPlaneMaker1) { tpc2EP_bin1 = (double)EventPlaneMaker1->GetTPCEP(); } else { tpc2EP_bin1 = -999; }
      if(EventPlaneMaker2) { tpc2EP_bin2 = (double)EventPlaneMaker2->GetTPCEP(); } else { tpc2EP_bin2 = -999; }
      if(EventPlaneMaker3) { tpc2EP_bin3 = (double)EventPlaneMaker3->GetTPCEP(); } else { tpc2EP_bin3 = -999; }
      if(EventPlaneMaker4) { tpc2EP_bin4 = (double)EventPlaneMaker4->GetTPCEP(); } else { tpc2EP_bin4 = -999; }

      // assign global event plane to selected pt-dependent bin
      if(fTPCptAssocBin == 0) TPC_PSI2 = tpc2EP_bin0;
      if(fTPCptAssocBin == 1) TPC_PSI2 = tpc2EP_bin1;
      if(fTPCptAssocBin == 2) TPC_PSI2 = tpc2EP_bin2;
      if(fTPCptAssocBin == 3) TPC_PSI2 = tpc2EP_bin3;
      if(fTPCptAssocBin == 4) TPC_PSI2 = tpc2EP_bin4;
      //cout<<"EP angles:  "<<tpc2EP_bin0<<"   "<<tpc2EP_bin1<<"   "<<tpc2EP_bin2<<"   "<<tpc2EP_bin3<<"   "<<tpc2EP_bin4<<endl;
    }

    // fill histogram with angle from GetReactionPlane() function
    angle1 = GetReactionPlane();
    hReactionPlaneFnc->Fill(angle1);

    // fill histogram with angle from GetEventPlane() function
    GetEventPlane(kFALSE, 2, fTPCEPmethod, 2.0, fTPCptAssocBin); // 2nd to last param not used (ptcut)
    angle2     = fEPTPC;
    angle2n    = fEPTPCn;
    angle2p    = fEPTPCp;
    hEventPlaneFncN2->Fill(angle2n);
    hEventPlaneFncP2->Fill(angle2p);
    hEventPlaneFnc2->Fill(angle2);

    // fill histogram with angle from StEventPlaneMaker class
    angle4 = TPC_PSI2;
    hEventPlaneClass->Fill(angle4);

    //cout<<"1: "<<angle1<<"  2: "<<angle2<<"  2n: "<<angle2n<<"  2p: "<<angle2p<<"  4: "<<angle4<<endl;
  }

  // FIXME FIXME - March 3, 2019 - this is outdated at this point, need to recode below chunk if doing ptbin approach
  // switch to require specific trigger (for Event Plane corrections + Resolution)
  if(doEPAnalysis) {
    // compare BBC, ZDC, TPC event planes
    // only truely relevant for STEP3 - when both recentering and shifting corrections are read in
    hTPCvsBBCep->Fill(BBC_PSI2, TPC_PSI2);
    hTPCvsZDCep->Fill(ZDC_PSI2, TPC_PSI2);
    hBBCvsZDCep->Fill(ZDC_PSI2, BBC_PSI2);
 
    // calculate / fill event plane resolution histograms
    if(doEventPlaneRes){
      // corrected event planes used for Resolution
      // MAKE sure for meaningful results to be in mode: STEP3
      CalculateEventPlaneResolution(BBC_PSI2, ZDC_PSI2, TPC_PSI2, TPCA_PSI2, TPCB_PSI2, BBC_PSI1, ZDC_PSI1);
    }
  } // have Emc HT trigger - process event plane calculation/ corrections / resolutions

  // ===================================================================================

  // ========================== Jet Shape Analysis ===================================== //
  int jsret = -99;
  if(doJetShapeAnalysis) {
    // declare pool pointer
    StEventPool *pool = 0x0;

    // require event mixing
    if(fDoEventMixing > 0) {
      // convert back to integer bins for mixed event pool - 10% bins (0, 7), 5% bins (0, 15)
      Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSizeJS);
      //cout<<"fCentralityScaled: "<<fCentralityScaled<<"  fCentBinSizeJS: "<<fCentBinSizeJS<<"  mixcentbin: "<<mixcentbin<<"  zVtx: "<<zVtx<<endl;
  
      // for an angle (0-pi)    
      const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
      StEventPlaneMaker *EPMaker = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 4)));
      double psi2 = (EPMaker) ? (double)EPMaker->GetTPCEP() : -999;

      double psibin = -1.;
      if(psi2 >=       0.0 && psi2 <  1.*pi/6.0) psibin = 0.5; // psibin = 1;
      if(psi2 >= 1.*pi/6.0 && psi2 <  2.*pi/6.0) psibin = 1.5; // psibin = 2;
      if(psi2 >= 2.*pi/6.0 && psi2 <  4.*pi/6.0) psibin = 2.5; // psibin = 3;
      if(psi2 >= 4.*pi/6.0 && psi2 <  6.*pi/6.0) psibin = 1.5; // psibin = 2;
      if(psi2 >= 5.*pi/6.0 && psi2 <= 6.*pi/6.0) psibin = 0.5; // psibin = 1;
      //cout<<"psi2 = "<<psi2<<"  psibin = "<<psibin<<endl;
    
      // initialize event pools - different cases for each dataset
      if(fDoUseMultBins) {
        if(!doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx);
        //if( doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx, psibin);
        if( doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx, psi2);
      } else {
        if(!doUseEPBins) pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); // FIXME AuAu fcent: cent bin? cent16
        if( doUseEPBins) pool = fPoolMgr->GetEventPool(mixcentbin, zVtx, psi2); // psibin
      }

      // check if pool exists
      if(!pool) {
        Form("No pool found for centrality = %i, zVtx = %f", mixcentbin, zVtx); // FIXME if cent changes to double
        return kTRUE;
      }

      //cout<<"nPsiIndex: "<<pool->PsiBinIndex()<<"   zvertexbinIndex: "<<pool->ZvtxBinIndex()<<"   multbinindex: "<<pool->MultBinIndex()<<endl;
    }

    // initialize some variables
    double ljpttemp = 0., sljpttemp = 0.;
    bool isBackToBack = kFALSE;

    // check for back-to-back jets
    if(fLeadingJet && fSubLeadingJet) { 
      double BackToBackPhi = fLeadingJet->Phi() - fSubLeadingJet->Phi() - pi;
      if(BackToBackPhi < 0.0) BackToBackPhi += 2.0*pi;
      if(BackToBackPhi < 0.4) isBackToBack = kTRUE;

      // get lj and sublj pt
      if(fCorrJetPt) {
        ljpttemp = fLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
        sljpttemp = fSubLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
      } else {
        ljpttemp = fLeadingJet->Pt();
        sljpttemp = fSubLeadingJet->Pt(); 
      }

    } // lj + sublj

    // Aj selection for Jet Shape Analysis: LJ pt > 20, SubLJ pt > 10 defaults
    bool doAjSelection = (isBackToBack && ljpttemp > fLeadJetPtMin && sljpttemp > fSubLeadJetPtMin) ? kTRUE : kFALSE;
    if(doAjSelection) hJetLeadingPtAj->Fill(ljpttemp);      // leading jet w/ Aj cut
    if(doAjSelection) hJetSubLeadingPtAj->Fill(sljpttemp);  // sub-leading jet w/ Aj cut

    // Triggered events and leading/subleading jets - do Jet Shape Analysis
    // check for back to back jets: must have leading + subleading jet, subleading jet must be > 10 GeV, subleading jet must be within 0.4 of pi opposite of leading jet
    if(doRequireAjSelection) {
      for(int ptbin=0; ptbin<9; ptbin++) {
        if(doAjSelection && fHaveEmcTrigger && fJetShapeJetType == kLeadingJets && fLeadingJet)       jsret = JetShapeAnalysis(fLeadingJet, pool, kEventActivity, ptbin);
        if(doAjSelection && fHaveEmcTrigger && fJetShapeJetType == kSubLeadingJets && fSubLeadingJet) jsret = JetShapeAnalysis(fSubLeadingJet, pool, kEventActivity, ptbin);
      }
    } else { // don't require back-to-back jets meeting Aj criteria
      for(int ptbin=0; ptbin<9; ptbin++) {
        if(fHaveEmcTrigger && fJetShapeJetType == kLeadingJets && fLeadingJet)       jsret = JetShapeAnalysis(fLeadingJet, pool, kEventActivity, ptbin);
        if(fHaveEmcTrigger && fJetShapeJetType == kSubLeadingJets && fSubLeadingJet) jsret = JetShapeAnalysis(fSubLeadingJet, pool, kEventActivity, ptbin);
      }
    }

    // jet shape - case for: inclusive jets
    if(fHaveEmcTrigger && fJetShapeJetType == kInclusiveJets) {

      // loop over Jets in the event
      for(int ijet = 0; ijet < fJets->GetEntries(); ijet++) {  // JET LOOP
        // get pointer to jet
        StJet *jet = static_cast<StJet*>(fJets->At(ijet));
        if(!jet) continue;

        // loop over pt associated bins
        for(int ptbin = 0; ptbin < 9; ptbin++) {
          if(doRequireAjSelection && doAjSelection) { jsret = JetShapeAnalysis(jet, pool, kEventActivity, ptbin);
          } else { jsret = JetShapeAnalysis(jet, pool, kEventActivity, ptbin); }

        }

      } // loop over jets
    }   // inclusive jet case

    // use only tracks from MB events
    //if(fDoEventMixing > 0 && fRunForMB && (!fHaveEmcTrigger)) { // kMB5 or kMB30 - AuAu, kMB - pp (excluding HT)
    if(fDoEventMixing > 0 && fRunForMB) { // kMB5 or kMB30 - AuAu, kMB - pp (don't exclude HT)
      // update pool: create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());

      // fill QA histo's
      hMBvsMult->Fill(kEventActivity);    // MB5 || MB30 AuAu and MB pp
      if(fHaveMB5event)  hMB5vsMult->Fill(refCorr2);   // MB5
      if(fHaveMB30event) hMB30vsMult->Fill(refCorr2);  // MB30
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB 

    // only return if doing jet shape analysis
    FillEventTriggerQA(fHistEventSelectionQAafterCuts);

    // only return if doing jet shape analysis
    return kStOK;  
  }

  // ======================================================================================================
  // ======================================================================================================
  // the code is cut off here when doing a jet-shape analysis
  // ======================================================================================================
  // ======================================================================================================
  // ======================================================================================================

  // run Track QA and fill histograms
  if((doWriteTrackQAHist) && (doJetAnalysis)) TrackQA();

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();
  if(fDebugLevel == kDebugGeneralEvt) {
    //cout<<"grefMult = "<<grefMult<<"  refMult = "<<refMult<<"  refCorr2 = "<<refCorr2<<"  cent16 = "<<cent16<<"   cent9 = "<<cent9<<"  centbin = "<<centbin<<endl;
    cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<"  refCorr2 = "<<refCorr2<<"  grefMult = "<<grefMult<<"  centbin = "<<centbin<<endl;
  }

  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables
  Int_t ijethi = -1;
  Double_t highestjetpt = 0.0;
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // Run - Trigger Selection to process jets from
    if(!doJetAnalysis) continue;

    // get pointer to jets
    StJet *jet = static_cast<StJet*>(fJets->At(ijet));
    if(!jet) continue;

    // get some jet parameters
    double jetarea = jet->Area();
    double jetpt = jet->Pt();
    double corrjetpt = jet->Pt() - jetarea*fRhoVal;
    double jetE = jet->E();
    double jetEta = jet->Eta();
    double jetPhi = jet->Phi();    
    double jetNEF = jet->NEF();
    //double dEP = RelativeEPJET(jet->Phi(), rpAngle);         // difference between jet and EP
    if(!doppAnalysis) {
      double dEPmethod1      = RelativeEPJET(jetPhi, angle1);       // from Reaction Plane function
      double dEPmethod2p     = RelativeEPJET(jetPhi, angle2p);      // from Event Plane function: sub event 1 - method 2
      double dEPmethod2n     = RelativeEPJET(jetPhi, angle2n);      // from Event Plane function: sub event 2 - method 2
      double dEPmethod2      = RelativeEPJET(jetPhi, angle2);       // from Event Plane function - method 2
      double dEPmethod4      = RelativeEPJET(jetPhi, angle4);       // from Event Plane class
      hdEPReactionPlaneFnc->Fill(dEPmethod1);
      hdEPEventPlaneFncN2->Fill(dEPmethod2n);
      hdEPEventPlaneFncP2->Fill(dEPmethod2p);
      hdEPEventPlaneFnc2->Fill(dEPmethod2);
      hdEPEventPlaneClass->Fill(dEPmethod4);
    }

    // get pt dependent event plane calculation
    // FIXME, double check, might want to code this nicer
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, TPC_PSI2) : 0; // CORRECTED event plane angle - STEP3
    if(doTPCptassocBin && !doppAnalysis) {
      // z = if(condition) then(?) <do this> else(:) <do this>  
      double dEP0 = (EventPlaneMaker0) ? RelativeEPJET(jetPhi, tpc2EP_bin0) : -999;
      double dEP1 = (EventPlaneMaker1) ? RelativeEPJET(jetPhi, tpc2EP_bin1) : -999;
      double dEP2 = (EventPlaneMaker2) ? RelativeEPJET(jetPhi, tpc2EP_bin2) : -999;
      double dEP3 = (EventPlaneMaker3) ? RelativeEPJET(jetPhi, tpc2EP_bin3) : -999;
      double dEP4 = (EventPlaneMaker4) ? RelativeEPJET(jetPhi, tpc2EP_bin4) : -999;
      if(fTPCptAssocBin == 0) dEP = dEP0;
      if(fTPCptAssocBin == 1) dEP = dEP1;
      if(fTPCptAssocBin == 2) dEP = dEP2;
      if(fTPCptAssocBin == 3) dEP = dEP3;
      if(fTPCptAssocBin == 4) dEP = dEP4;

      if(fDebugLevel == kDebugJetvsEPtype) cout<<"jetPhi = "<<jetPhi<<"  RP = "<<angle1<<"  bin0 = "<<tpc2EP_bin0<<"  bin1 = "<<tpc2EP_bin1<<"  bin2 = "<<tpc2EP_bin2<<"  bin3 = "<<tpc2EP_bin3<<"  bin4 = "<<tpc2EP_bin4<<endl;
      if(fDebugLevel == kDebugJetvsEPtype) cout<<"dRP = "<<dEP1<<"  dEP0: "<<dEP0<<"  dEP1: "<<dEP1<<"  dEP2: "<<dEP2<<"  dEP3: "<<dEP3<<"  dEP4: "<<dEP4<<endl;
    }

    // this is a double check - should not happen, but if it does -> kill the job so it can be fixed - most likely in the runPico macro
    if(dEP < -900) return kStFatal;

    // some threshold cuts
    if(fCorrJetPt) {  // background subtracted jet pt
      if(corrjetpt < fMinPtJet) continue;
    } else { if(jetpt < fMinPtJet) continue; }
    if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerE() < fTowerBias)) continue;

    // TODO check that jet contains a tower that fired the trigger
    if(!DidTowerConstituentFireTrigger(jet)) { continue; }

    // loop over constituent tracks
    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);      

      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
      if(!trk){ continue; }

      // track momentum vector
      TVector3 mTrkMom;
      if(doUsePrimTracks) {
        if(!(trk->isPrimary())) return kFALSE; // check if primary
        // get primary track vector
        mTrkMom = trk->pMom();
      } else {
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield);
      }

      // track variables
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      double pt = mTrkMom.Perp();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      double jetZ = jet->GetZ(px, py, pz);

      // shift angle (0, 2*pi) 
      if(phi < 0.0)    phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;

      // fill jet track constituent histograms
      hJetTracksPt->Fill(pt);
      hJetTracksPhi->Fill(phi);
      hJetTracksEta->Fill(eta);
      hJetTracksZ->Fill(jetZ);

    } // constituent track loop

/*
    // loop over constituent towers
    for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
      int ArrayIndex = jet->ClusterAt(itow);
      StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(ArrayIndex));
      if(!tow){ continue; }

      //int towID = tow->id(); // ArrayIndex = towID - 1 because of array element numbering different than ids which start at 1
      int towID = -1;
      if( gROOT->GetClass("StPicoBTowHit")->GetClassVersion() < 3) {
        //towID = tow->id();
      } else {
        //towID = tow->numericIndex2SoftId(ArrayIndex);
      }

      // tower ID: get from index of array shifted by +1
      towID = ArrayIndex + 1;
      if(towID < 0) continue;

      int containsTower = jet->ContainsTower(ArrayIndex);
      cout<<">= 0: "<<containsTower<<"  itow = "<<itow<<"  id = "<<towID<<"  ArrIndex = "<<ArrayIndex<<"  towE = "<<tow->energy()<<endl;
    } // constituent tower loop
*/

    // fill some jet histograms
    hJetPt->Fill(jetpt);
    hJetCorrPt->Fill(corrjetpt);
    hJetE->Fill(jetE);
    hJetEta->Fill(jetEta);
    hJetPhi->Fill(jetPhi);
    hJetNEF->Fill(jetNEF);
    hJetArea->Fill(jetarea);
    hJetPtvsArea->Fill(jetpt, jetarea);
    hJetEventEP->Fill(TPC_PSI2);
    hJetPhivsEP->Fill(jetPhi, TPC_PSI2);

    // fill some jet QA plots for each orientation
    if(dEP >= 0.0*pi/6.0 && dEP < 1.0*pi/6.0) {
      hJetPtIn->Fill(jetpt);
      hJetPhiIn->Fill(jetPhi);
      hJetEtaIn->Fill(jetEta);
      hJetEventEPIn->Fill(TPC_PSI2);
      hJetPhivsEPIn->Fill(jetPhi, TPC_PSI2);
    } else if(dEP >= 1.0*pi/6.0 && dEP < 2.0*pi/6.0) {
      hJetPtMid->Fill(jetpt);
      hJetPhiMid->Fill(jetPhi);
      hJetEtaMid->Fill(jetEta);
      hJetEventEPMid->Fill(TPC_PSI2);
      hJetPhivsEPMid->Fill(jetPhi, TPC_PSI2);
    } else if(dEP >= 2.0*pi/6.0 && dEP <= 3.0*pi/6.0) {
      hJetPtOut->Fill(jetpt);
      hJetPhiOut->Fill(jetPhi);
      hJetEtaOut->Fill(jetEta);
      hJetEventEPOut->Fill(TPC_PSI2);
      hJetPhivsEPOut->Fill(jetPhi, TPC_PSI2);
    }

    // get nTracks and maxTrackPt
    double maxtrackpt = jet->GetMaxTrackPt();
    int NtrackConstit = jet->GetNumberOfTracks();
    if(doComments) cout<<"Jet# = "<<ijet<<"  JetPt = "<<jetpt<<"  JetE = "<<jetE<<endl;
    if(doComments) cout<<"MaxtrackPt = "<<maxtrackpt<<"  NtrackConstit = "<<NtrackConstit<<endl;

    // get highest Pt jet in event (leading jet)
    if(highestjetpt < jetpt){
      ijethi = ijet;
      highestjetpt = jetpt;
    }

    // the below track and cluster cut is already done
    //if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fTowerBias)){
    // set up and fill Jet-Hadron trigger jets THnSparse
    // ====================================================================================
    double jetptselected;
    if(fCorrJetPt) { jetptselected = corrjetpt; 
    } else { jetptselected = jetpt; }

    // fill jet sparse for trigger normalization
    Double_t CorrEntries[4] = {centBinToUse, jetptselected, dEP, zVtx};
    if(fReduceStatsCent > 0) {
      if(cbin == fReduceStatsCent) fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    } else fhnCorr->Fill(CorrEntries);    // fill Sparse Histo with trigger Jets entries
    // ====================================================================================
    //} // check on max track and tower pt/E

    // ======================================================================================

/*
    // this is for jet constituents
    vector<fastjet::PseudoJet> fConstituents = jet->GetJetConstituents();
    for(UInt_t ic = 0; ic < fConstituents.size(); ++ic) {
      // get user defined index
      Int_t uid = fConstituents[ic].user_index();
      double cpt = fConstituents[ic].perp();
      double ceta = fConstituents[ic].eta();
      double cphi = fConstituents[ic].phi();
      cout<<"ic = "<<ic<<", uid = "<<uid<<", cpt = "<<cpt<<", ceta = "<<ceta<<", cphi = "<<cphi<<endl;
    }
*/

    // =============== jet shape analysis ================
    for(int ptbin = 0; ptbin < 9; ptbin++) {
    ////// if(doJetShapeAnalysis) JetShapeAnalysis(jet, 0x0, refCorr2, ptbin); // FIXME
    }

    // track loop inside jet loop - loop over ALL tracks in PicoDst
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

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
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      short charge = trk->charge();

      // get jet - track relations
      //deta = eta - jetEta;               // eta betweeen hadron and jet
      double deta = jetEta - eta;               // eta betweeen jet and hadron
      double dphijh = RelativePhi(jetPhi, phi); // angle between jet and hadron

      // === June 1, 2018 - these cuts should not be here!
/*
      // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
      // when doing event plane calculation via pt assoc bin
      if(doTPCptassocBin) {
        if(fTPCptAssocBin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 4) { if((pt > 2.00) && (pt <= 30.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
        if(fTPCptAssocBin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(fTPCptAssocBin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
      }
*/

      // fill jet sparse 
      double triggerEntries[8] = {centBinToUse, jetptselected, pt, deta, dphijh, dEP, zVtx, (double)charge};
      double trefficiency = 1.0;
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
  if(fDebugLevel == kDebugMixedEvents) cout<<"StMyAnMaker event# = "<<EventCounter()<<"  Centbin = "<<centbin<<"  zVtx = "<<zVtx<<endl;

  // Prepare to do event mixing
  if(fDoEventMixing > 0){
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

    //Double_t mycentbin = (Double_t)centbin + 0.001;
    // mixed event centbin - FIXME
    Int_t mixcentbin;
    if(fCentBinSize==10)       { mixcentbin = centbin10;    // 10% bins (0,  7)
    } else if(fCentBinSize==5) { mixcentbin = centbin; }    //  5% bins (0, 15)
    //cout<<"mixcentbin = "<<mixcentbin<<"  centbin = "<<centbin<<" centbin10 = "<<centbin10<<"  zvtx = "<<zVtx<<endl;

    // initialize event pools
    StEventPool *pool = 0x0;
    pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); //FIXME AuAu fcent: cent bin? cent16
    if(!pool) {
      Form("No pool found for centrality = %i, zVtx = %f", mixcentbin, zVtx); // FIXME if cent changes to double
      return kTRUE; //FIXME
    }

    if(fDebugLevel == kDebugMixedEvents) cout<<"NtracksInPool = "<<pool->NTracksInPool()<<"  CurrentNEvents = "<<pool->GetCurrentNEvents()<<endl;

    // initialize background tracks array
    TObjArray *bgTracks;

  // do event mixing when Signal Jet is part of event with a HT1 or HT2 or HT3 trigger firing
  if(doJetAnalysis) { // trigger type requested was fired for this event - do mixing
    if(pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {

      // loop over jets (passing cuts - set by jet maker)
      for(int ijet = 0; ijet < njets; ijet++) {
        // leading jet - why was this a double?
        Double_t leadjet = 0;
        if(ijet == ijethi) leadjet = 1; //FIXME for leading jet

        // get jet pointer
        StJet *jet = static_cast<StJet*>(fJets->At(ijet));
        if(!jet) continue;

        // get some get parameters of jets for mixing
        double Mixjetarea = jet->Area();
        double Mixjetpt = jet->Pt();
        double Mixcorrjetpt = Mixjetpt - Mixjetarea*fRhoVal;
        double MixjetEta = jet->Eta();
        double MixjetPhi = jet->Phi();
        //double dMixEP = RelativeEPJET(jet->Phi(), rpAngle);         // difference between jet and EP
        double dMixEP = (!doppAnalysis) ? RelativeEPJET(jet->Phi(), TPC_PSI2) : 0; // CORRECTED event plane angle - STEP3
        if(doTPCptassocBin && !doppAnalysis) {
          // z = if(condition) then(?) <do this> else(:) <do this>  
          double dMixEP0 = (EventPlaneMaker0) ? RelativeEPJET(MixjetPhi, tpc2EP_bin0) : -999;
          double dMixEP1 = (EventPlaneMaker1) ? RelativeEPJET(MixjetPhi, tpc2EP_bin1) : -999;
          double dMixEP2 = (EventPlaneMaker2) ? RelativeEPJET(MixjetPhi, tpc2EP_bin2) : -999;
          double dMixEP3 = (EventPlaneMaker3) ? RelativeEPJET(MixjetPhi, tpc2EP_bin3) : -999;
          double dMixEP4 = (EventPlaneMaker4) ? RelativeEPJET(MixjetPhi, tpc2EP_bin4) : -999;
          if(fTPCptAssocBin == 0) dMixEP = dMixEP0;
          if(fTPCptAssocBin == 1) dMixEP = dMixEP1;
          if(fTPCptAssocBin == 2) dMixEP = dMixEP2;
          if(fTPCptAssocBin == 3) dMixEP = dMixEP3;
          if(fTPCptAssocBin == 4) dMixEP = dMixEP4;
        }

        // this is a double check - should not happen, but if it does -> kill the job so it can be fixed - most likely in the runPico macro
        if(dMixEP < -900) return kStFatal;

        // some threshold cuts - do mixing only if we have a jet meeting our pt threshold and bias
        if(fCorrJetPt) {
          if(Mixcorrjetpt < fMinPtJet) continue;
        } else { if(Mixjetpt < fMinPtJet) continue; }
        //if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerE() < fTowerBias)) continue;
   	//TODO if (!AcceptJet(jet)) continue;  // acceptance cuts done to jet in JetMaker

        // get number of current events in pool
        int nMix = pool->GetCurrentNEvents();

        // check that jet contains a tower that fired the trigger
        if( !DidTowerConstituentFireTrigger(jet) ) continue;

        // Fill for biased jet triggers only
        if((jet->GetMaxTrackPt() > fTrackBias) || (jet->GetMaxTowerE() > fTowerBias)) {
          // Fill mixed-event histos here: loop over nMix events
          for(int jMix = 0; jMix < nMix; jMix++) {
 
            // get jMix'th event
            bgTracks = pool->GetEvent(jMix);
            //TObjArray *bgTracks = pool->GetEvent(jMix);
            const Int_t Nbgtrks = bgTracks->GetEntries();

            // loop over background (mixed event) tracks
            for(int ibg = 0; ibg < Nbgtrks; ibg++) {
              // get Femto track pointer
              StFemtoTrack *trk = static_cast<StFemtoTrack*>(bgTracks->At(ibg));
              if(!trk) continue;
              double Mixphi = trk->Phi();
              double Mixeta = trk->Eta();
              double Mixpt = trk->Pt();
              short Mixcharge = trk->Charge();

              // shift angle (0, 2*pi) 
              if(Mixphi < 0.0)    Mixphi += 2.0*pi;
              if(Mixphi > 2.0*pi) Mixphi -= 2.0*pi;
              //cout<<"itrack = "<<ibg<<"  phi = "<<Mixphi<<"  eta = "<<Mixeta<<"  pt = "<<Mixpt<<"  q = "<<Mixcharge<<endl;

              // get jet - track relations
              //double deta = eta - jetEta;               // eta betweeen hadron and jet
              double dMixeta = MixjetEta - Mixeta;               // eta betweeen jet and hadron
              double dMixphijh = RelativePhi(MixjetPhi, Mixphi); // angle between jet and hadron

              // print tracks outside of acceptance somehow
              if(fDebugLevel == kDebugMixedEvents) if(TMath::Abs(dMixeta > 1.6)) cout<<"DELTA ETA out of bounds... deta = "<<dMixeta<<"   iTrack = "<<ibg<<"  jetEta = "<<MixjetEta<<"  trk eta = "<<Mixeta<<endl;

              // calculate single particle tracking efficiency of mixed events for correlations (-999)
              double mixefficiency = 1.0;
              //FIXME mixefficiency = EffCorrection(part->Eta(), part->Pt(), fDoEffCorr);                           

              // select which jet pt to use for filling
              double Mixjetptselected;
              if(fCorrJetPt) { Mixjetptselected = Mixcorrjetpt;
              } else { Mixjetptselected = Mixjetpt; }

              // === June 1, 2018 - these cuts should not be here!
/*
              // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
              // when doing event plane calculation via pt assoc bin
              if(doTPCptassocBin) {
                if(fTPCptAssocBin == 0) { if((Mixpt > 0.20) && (Mixpt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 1) { if((Mixpt > 0.50) && (Mixpt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 2) { if((Mixpt > 1.00) && (Mixpt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 3) { if((Mixpt > 1.50) && (Mixpt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 4) { if((Mixpt > 2.00) && (Mixpt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
                if(fTPCptAssocBin == 5) { if((Mixpt > 2.00) && (Mixpt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 6) { if((Mixpt > 3.00) && (Mixpt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
                if(fTPCptAssocBin == 7) { if((Mixpt > 4.00) && (Mixpt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
              }
*/

              // create / fill mixed event sparse  (centbin*5.0)
              double triggerEntries[8] = {centBinToUse, Mixjetptselected, Mixpt, dMixeta, dMixphijh, dMixEP, zVtx, (double)Mixcharge};

              if(fReduceStatsCent > 0) {
                if(cbin == fReduceStatsCent) fhnMixedEvents->Fill(triggerEntries,1./(nMix*mixefficiency));
              } else fhnMixedEvents->Fill(triggerEntries,1./(nMix*mixefficiency));   // fill Sparse histo of mixed events

            } // end of background track loop
          } // end of filling mixed-event histo's:  jth mix event loop
        } // end of check for biased jet triggers
      } // end of jet loop
    } // end of check for pool being ready
  } // end EMC triggered loop

    // use only tracks from MB (and Semi-Central) events
    ///if(fMixingEventType) { //kMB) {
    //if(fRunForMB && (!fHaveEmcTrigger)) { // kMB or kMB30 - AuAu, kMB - pp, TODO probably want to use this line in future, may not matter
    //if(fRunForMB) { // kMB or kMB30 (don't exclude HT)
    if(fHaveMBevent) { // kMB - FIXME, this is how it was set up for correlations
      // update pool - create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());

      // fill QA histo's
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB event - updating pool

  } // end of event mixing

  // event counter at end of maker
  mInputEventCounter++;
  //cout<<"end of event counter = "<<mInputEventCounter<<endl;

  // fill Event Trigger QA
  //FillEventTriggerQA(fHistEventSelectionQAafterCuts);
  //StMemStat::PrintMem("MyAnalysisMaker at end of make");

  return kStOK;
}
//
//______________________________________________________________________
THnSparse* StMyAnalysisMaker3::NewTHnSparseF(const char* name, UInt_t entries)
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
//
//________________________________________________________________________
void StMyAnalysisMaker3::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "centrality 5% bin";
      if(fCentBinSize==10) {
        nbins = 10;
      } else {
        nbins = 20; //16;
      } 
      xmin = 0.;
      xmax = 100.; //16.;     
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
//
//______________________________________________________________________
THnSparse* StMyAnalysisMaker3::NewTHnSparseFCorr(const char* name, UInt_t entries) {
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
//
//______________________________________________________________________________________________
void StMyAnalysisMaker3::GetDimParamsCorr(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   //stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

    case 0:
      label = "centrality 5% bin";
      if(fCentBinSize==10) {
        nbins = 10;
      } else {
        nbins = 20; //16;
      }
      xmin = 0.;
      xmax = 100.; //16.;
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
//
// From CF event mixing code PhiCorrelations
// clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
//____________________________________________________________________________
TClonesArray* StMyAnalysisMaker3::CloneAndReduceTrackList()
{
  TClonesArray *tracksClone = new TClonesArray("StFemtoTrack");
//  tracksClone->SetName("tracksClone");
//  tracksClone->SetOwner(kTRUE);

  // construct variables, get # of tracks
  int nMixTracks = mPicoDst->numberOfTracks();
  int iterTrk = 0;

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
      // get primary track vector
      mTrkMom = trk->pMom();
    } else {
      // get global track vector
      mTrkMom = trk->gMom(mVertex, Bfield);
    }

    // track variables - used with alt method below
    double pt = mTrkMom.Perp();
    //double phi = mTrkMom.Phi();
    //double eta = mTrkMom.PseudoRapidity();
    //short charge = trk->charge();

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    // this is TEMP, it will filter track by the pt bin used for analysis
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

    // create StFemtoTracks out of accepted tracks - light-weight object for mixing
    //  StFemtoTrack *t = new StFemtoTrack(pt, eta, phi, charge);
    StFemtoTrack* t = new StFemtoTrack(trk, Bfield, mVertex, doUsePrimTracks);
    if(!t) continue;

    // add light-weight tracks passing cuts to TClonesArray
    ((*tracksClone)[iterTrk]) =  t;

    //delete t;
    ++iterTrk;
  } // end of looping through tracks

  return tracksClone;
}
//
/*
//________________________________________________________________________
Bool_t StMyAnalysisMaker3::AcceptJet(StJet *jet) { // for jets
  // applies all jet cuts except pt
  if ((jet->Phi() < fPhimin) || (jet->Phi() > fPhimax)) return kFALSE;
  if ((jet->Eta() < fEtamin) || (jet->Eta() > fEtamax)) return kFALSE;
  if (jet->Area() < fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area() == 0) return kFALSE;
  // exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt() > 20) return kFALSE;

  // jet passed all above cuts
  return kTRUE;
}
*/
//
//
//_________________________________________________________________________
TH1* StMyAnalysisMaker3::FillEmcTriggersHist(TH1* h) {
  // number of Emcal Triggers
  for(int i = 0; i < 8; i++) { fEmcTriggerArr[i] = 0; }
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  if(fDebugLevel == kDebugEmcTrigger) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

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
// Trigger QA histogram, label bins 
// check and fill a Event Selection QA histogram for different trigger selections after cuts
//_____________________________________________________________________________
TH1* StMyAnalysisMaker3::FillEventTriggerQA(TH1* h) {
  // Run12 pp 200 GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    // Run12 (200 GeV pp) triggers:
    int arrHT1[] = {370511, 370546};
    int arrHT2[] = {370521, 370522, 370531, 370980};
    //int arrHT3[] = {380206, 380216}; // NO HT3 triggered events
    int arrMB[] = {370001, 370011, 370983};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { bin = 3; h->Fill(bin); } // HT2
    //if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { bin = 4; h->Fill(bin); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 10; h->Fill(bin); } // VPDMB

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
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
    int arrMB30[] = {450010, 450020};
    int arrCentral5[] = {450010, 450020};
    int arrCentral[] = {460101, 460111};
    int arrMB5[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrBHT3, sizeof(arrBHT3)/sizeof(*arrBHT3))) { bin = 4; h->Fill(bin); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); } // MB 
    if(DoComparison(arrCentral5, sizeof(arrCentral5)/sizeof(*arrCentral5))) { bin = 7; h->Fill(bin); }// Central-5
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 8; h->Fill(bin); } // Central & Central-mon
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { bin = 10; h->Fill(bin); }// VPDMB-5 
    if(DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { bin = 11; h->Fill(bin); } // VPDMB-30

    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2)) && DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 13; h->Fill(bin); } // HT2 && MB
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2)) && DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { bin = 14; h->Fill(bin); } // HT2 && MB30
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1)) && DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 15; h->Fill(bin); } // HT1 && MB
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1)) && DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { bin = 16; h->Fill(bin); } // HT1 && MB30

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
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrBHT3, sizeof(arrBHT3)/sizeof(*arrBHT3))) { bin = 4; h->Fill(bin); } // HT3
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); }  // MB
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 7; h->Fill(bin); }// Central-5 & Central-novtx
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { bin = 10; h->Fill(bin); } // VPDMB-5 
    if(DoComparison(arrMB10, sizeof(arrMB10)/sizeof(*arrMB10))) { bin = 11; h->Fill(bin); }// VPDMB-10

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
  h->GetXaxis()->SetBinLabel(1, "un-identified trigger");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");
  
  return h;
}
//
// general function to get local reaction plane calculation
//________________________________________________________________________
Double_t StMyAnalysisMaker3::GetReactionPlane() { 
  // initliaze some variables and constants
  TVector2 mQ;
  double mQx = 0., mQy = 0.;
  int order = 2;
  double pi = 1.0*TMath::Pi();

  // leading jet check and removal
  Double_t excludeInEta = -999;
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from EP estimate
    if(fLeadingJet) excludeInEta = fLeadingJet->Eta();
  }

  // loop over tracks
  int nTrack = mPicoDst->numberOfTracks();
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
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // check for leading jet removal - taken from Redmers approach (CHECK! TODO!)
    if((fLeadingJet) && 
       (fExcludeLeadingJetsFromFit > 0) && 
       ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit) || 
       (TMath::Abs(eta) - fJetRad - 1.0 ) > 0 )) continue;

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

  // set q-vector components 
  mQ.Set(mQx, mQy);
  double psi = mQ.Phi() / order;

  return psi;
}
//
// Set the bin errors on histograms
// __________________________________________________________________________________
void StMyAnalysisMaker3::SetSumw2() {
  // set sum weights
  hdEPReactionPlaneFnc->Sumw2();
  hdEPEventPlaneFncN2->Sumw2();
  hdEPEventPlaneFncP2->Sumw2();
  hdEPEventPlaneFnc2->Sumw2();
  hdEPEventPlaneClass->Sumw2();
  hReactionPlaneFnc->Sumw2();
  hEventPlaneFncN2->Sumw2();
  hEventPlaneFncP2->Sumw2();
  hEventPlaneFnc2->Sumw2();
  hEventPlaneClass->Sumw2();

  hEventPlane->Sumw2();
  fHistEPTPCn->Sumw2();
  fHistEPTPCp->Sumw2();
  fHistEPBBC->Sumw2();
  fHistEPZDC->Sumw2();
  //hEventZVertex->Sumw2();
  //hCentrality->Sumw2();
  //hMultiplicity->Sumw2();
  //hRhovsCent->Sumw2();
  for(int i=0; i<5; i++) { hdEPtrk[i]->Sumw2(); }
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i]->Sumw2();
    hTrackEta[i]->Sumw2();
    hTrackPt[i]->Sumw2();
  }
  hTrackEtavsPhi->Sumw2();

  hJetPt->Sumw2();
  hJetCorrPt->Sumw2();
  hJetLeadingPt->Sumw2();
  hJetSubLeadingPt->Sumw2();
  hJetLeadingPtAj->Sumw2();
  hJetSubLeadingPtAj->Sumw2();
  hJetDiJetAj->Sumw2();
  hJetE->Sumw2();
  hJetEta->Sumw2();
  hJetPhi->Sumw2();
  hJetNEF->Sumw2();
  hJetArea->Sumw2();
  hJetTracksPt->Sumw2();
  hJetTracksPhi->Sumw2();
  hJetTracksEta->Sumw2();
  hJetTracksZ->Sumw2();
  hJetPtvsArea->Sumw2();
  hJetEventEP->Sumw2();
  hJetPhivsEP->Sumw2();

  hJetPtIn->Sumw2();
  hJetPhiIn->Sumw2();
  hJetEtaIn->Sumw2();
  hJetEventEPIn->Sumw2();
  hJetPhivsEPIn->Sumw2();
  hJetPtMid->Sumw2();
  hJetPhiMid->Sumw2();
  hJetEtaMid->Sumw2();
  hJetEventEPMid->Sumw2();
  hJetPhivsEPMid->Sumw2();
  hJetPtOut->Sumw2();
  hJetPhiOut->Sumw2();
  hJetEtaOut->Sumw2();
  hJetEventEPOut->Sumw2();
  hJetPhivsEPOut->Sumw2();

  fHistJetHEtaPhi->Sumw2();
  //fHistEventSelectionQA->Sumw2();
  //fHistEventSelectionQAafterCuts->Sumw2();
  //hTriggerIds->Sumw2();
  //hEmcTriggers->Sumw2();
  hMixEvtStatZVtx->Sumw2();
  hMixEvtStatCent->Sumw2();
  hMixEvtStatZvsCent->Sumw2();
  hTriggerEvtStatZVtx->Sumw2();
  hTriggerEvtStatCent->Sumw2();
  hTriggerEvtStatZvsCent->Sumw2();
  hMBvsMult->Sumw2();
  hMB5vsMult->Sumw2();
  hMB30vsMult->Sumw2();
  hHTvsMult->Sumw2();
  hNMixEvents->Sumw2();

  hTPCvsBBCep->Sumw2();
  hTPCvsZDCep->Sumw2();
  hBBCvsZDCep->Sumw2();

/*
  // don't set Sumw2 for TProfile plots
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        fProfJetV2[k][j][i]->Sumw2();
      }
    }
  }
*/

  if(doJetShapeAnalysis) {
    for(int k=0; k<4; k++) {
      for(int j=0; j<4; j++) {
        for(int i=0; i<4; i++) {
          hJetCounter[k][j][i]->Sumw2();
          hJetCounterCase1[k][j][i]->Sumw2();
          hJetCounterCase2[k][j][i]->Sumw2();
          hJetCounterCase3BG[k][j][i]->Sumw2();

          for(int p=0; p<9; p++) {
            hJetShape[k][j][i][p]->Sumw2();
            hJetShapeCase1[k][j][i][p]->Sumw2();
            hJetShapeCase2[k][j][i][p]->Sumw2();
            hJetShapeBG[k][j][i][p]->Sumw2();
            hJetShapeBGCase1[k][j][i][p]->Sumw2();
            hJetShapeBGCase2[k][j][i][p]->Sumw2();
            hJetShapeBGCase3[k][j][i][p]->Sumw2();

            hJetPtProfile[k][j][i][p]->Sumw2();
            hJetPtProfileCase1[k][j][i][p]->Sumw2();
            hJetPtProfileCase2[k][j][i][p]->Sumw2();
            hJetPtProfileBG[k][j][i][p]->Sumw2();
            hJetPtProfileBGCase1[k][j][i][p]->Sumw2();
            hJetPtProfileBGCase2[k][j][i][p]->Sumw2();
            hJetPtProfileBGCase3[k][j][i][p]->Sumw2();
          }
        }
      }
    }
  }

/*
  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      fProfV2Resolution[i]->Sumw2();
      fProfV3Resolution[i]->Sumw2();
      fProfV4Resolution[i]->Sumw2();
      fProfV5Resolution[i]->Sumw2();
    }
  }
*/

  fhnJH->Sumw2();
  fhnMixedEvents->Sumw2();
  fhnCorr->Sumw2();
}
//
//
// ________________________________________________________________________________________
void StMyAnalysisMaker3::GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin)
{ 
  // local variables
  TVector2 mQtpcn, mQtpcp, mQtpc;
  double mQtpcnx = 0., mQtpcny = 0., mQtpcpx = 0., mQtpcpy = 0., mQtpcX = 0., mQtpcY = 0.;
  int order = n;
  double pi = 1.0*TMath::Pi();
  int ntracksNEG = 0, ntracksPOS = 0;

  // leading jet check and removal
  double excludeInEta = -999., excludeInPhi = -999.;
  double excludeInEtaSub = -999., excludeInPhiSub = -999.;
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from EP estimate
    // check for leading jet
    if(fLeadingJet) {
      excludeInEta = fLeadingJet->Eta();
      excludeInPhi = fLeadingJet->Phi();
      //cout<<"leading: pt = "<<fLeadingJet->Pt()<<"  eta = "<<fLeadingJet->Eta()<<"  phi = "<<fLeadingJet->Phi()<<endl;
    }

    // check for subleading jet
    if(fSubLeadingJet) {
      excludeInEtaSub = fSubLeadingJet->Eta();
      excludeInPhiSub = fSubLeadingJet->Phi();
      //cout<<"subleading: pt = "<<fSubLeadingJet->Pt()<<"  eta = "<<fSubLeadingJet->Eta()<<"  phi = "<<fSubLeadingJet->Phi()<<endl;
    }
  } // leading jets

  // loop over tracks
  TRandom3 *rand = new TRandom3();
  int nTOT = 0, nA = 0, nB = 0;
  int nTrack = mPicoDst->numberOfTracks();
  for(int i=0; i<nTrack; i++) {
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
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    if(doTPCptassocBin) {
      if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations - FIXME why is this 0.25??
      if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(ptbin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(ptbin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(ptbin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    // remove strip only when we have a leading jet
    // Method1: kRemoveEtaStrip
    if(fTPCEPmethod == 1){
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && 
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    } else if(fTPCEPmethod == 2){
      // remove cone (in eta and phi) around leading jet
      // Method2: kRemoveEtaPhiCone - FIXME found bug May25
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    } else if(fTPCEPmethod == 3){
      // remove tracks above 2 GeV in cone around leading jet
      // Method3: kRemoveLeadingJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
    } else if(fTPCEPmethod == 4){
      // remove strip only when we have a leading + subleading jet
      // Method4: kRemoveEtaStripLeadSubLead
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEtaSub) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    } else if(fTPCEPmethod == 5){
      // remove cone (in eta and phi) around leading + subleading jet
      // Method5: kRemoveEtaPhiConeLeadSubLead
      double deltaR    = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi-excludeInPhiSub)*(phi-excludeInPhiSub));
      if((fLeadingJet)    && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR    < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaRSub < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    } else if(fTPCEPmethod == 6){
      // remove tracks above 2 GeV in cone around leading + subleading jet
      // Method6: kRemoveLeadingSubJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi - excludeInPhiSub)*(phi - excludeInPhiSub));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaRSub < fJetRad)) continue;
    } else {
      // DO NOTHING! nothing is removed...
    }

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

    // test - Jan15 for random subevents
    // generate random distribution from 0 -> 1: and split subevents for [0,0.5] and [0.5, 1]
    double randomNum = rand->Rndm();
    ////double randomNum = gRandom->Rndm();  // > 0.5?

    // split up Q-vectors into 2 random TPC sub-events (A and B)
    if(doTPCptassocBin) {
      if(randomNum >= 0.5) {
        // sum up q-vectors for random event A
        mQtpcpx += trackweight * cos(phi * order);
        mQtpcpy += trackweight * sin(phi * order);
        nA++;
      }
      if(randomNum < 0.5) {
        // sum up q-vectors for random event B
        mQtpcnx += trackweight * cos(phi * order);
        mQtpcny += trackweight * sin(phi * order);
        nB++;
      }
    }  // pt-dependent mode

    // non-pt dependent mode
    // split up Q-vectors to +/- eta regions
    if(!doTPCptassocBin) {
      if(eta < 0) {
        // sum up q-vectors on negative eta side
        mQtpcnx += trackweight * cos(phi * order);
        mQtpcny += trackweight * sin(phi * order);
        ntracksNEG++;
      }
      if(eta >= 0) {
        // sum up q-vectors on positive eta side
        mQtpcpx += trackweight * cos(phi * order);
        mQtpcpy += trackweight * sin(phi * order);
        ntracksPOS++;
      }
    }  

    // combined TPC event plane q-vectors
    mQtpcX += trackweight * cos(phi * order);
    mQtpcY += trackweight * sin(phi * order);
    nTOT++;

  } // end of track loop

  // test statements
  //cout<<"ntracksNEG = "<<ntracksNEG<<"  ntracksPOS = "<<ntracksPOS<<endl;
  //cout<<"mQtpcpx = "<<mQtpcpx<<"  mQtpcpy = "<<mQtpcpy<<"  mQtpcnx = "<<mQtpcnx<<"  mQtpcny = "<<mQtpcny<<endl;
  //cout<<"nTOT = "<<nTOT<<"  nA = "<<nA<<"  nB = "<<nB<<endl;

  // set q-vector components 
  mQtpcn.Set(mQtpcnx, mQtpcny);
  mQtpcp.Set(mQtpcpx, mQtpcpy);
  mQtpc.Set(mQtpcX, mQtpcY);

  // test..
  TVector2 mQtpcComb;
  mQtpcComb.Set(mQtpcnx + mQtpcpx, mQtpcny + mQtpcpy);
  fEPTPCcomb = mQtpcComb.Phi() / order;

  // Calculate the Event Plane
  fEPTPCn = mQtpcn.Phi() / order;
  fEPTPCp = mQtpcp.Phi() / order;
  fEPTPC = mQtpc.Phi() / order;

  // make sure event plane is contrained from [0, pi]
  if((fEPTPCn > pi) || (fEPTPCn < 0)) cout<<"fEPTPCn out of range..... "<<fEPTPCn<<endl;
  if((fEPTPCp > pi) || (fEPTPCp < 0)) cout<<"fEPTPCp out of range..... "<<fEPTPCp<<endl;
  if((fEPTPC  > pi) || (fEPTPC  < 0)) cout<<"fEPTPC out of range..... "<<fEPTPC<<endl;
  if(fEPTPCn > pi) fEPTPCn -= pi;
  if(fEPTPCn <  0) fEPTPCn += pi;
  if(fEPTPCp > pi) fEPTPCp -= pi;
  if(fEPTPCp <  0) fEPTPCp += pi;
  if(fEPTPC > pi)  fEPTPC  -= pi;
  if(fEPTPC <  0)  fEPTPC  += pi;

  // combine x-y vectors for neg and pos Eta ranges (cross-check)
  //double tpc2 = (0.5*TMath::ATan2(mQtpcY, mQtpcX));

  // standard event plane distributions as function of centrality
  fHistEPTPCn->Fill(fCentralityScaled, fEPTPCn);
  fHistEPTPCp->Fill(fCentralityScaled, fEPTPCp);
}
//
// 1) get the binning for: ref9 and region_vz
// 2) get function: GetRunNo( );
// 3) 
//
//
// this function checks for the bin number of the run from a runlist header 
// in order to apply various corrections and fill run-dependent histograms
// _________________________________________________________________________________
Int_t StMyAnalysisMaker3::GetRunNo(int runid){ 
  //1287 - Liang

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
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    // 1654 for Run14 AuAu, new picoDst production is 830
    for(int i = 0; i < 830; i++) {
      if(Run14AuAu_P18ih_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    // 1359 for Run16 AuAu
    for(int i = 0; i < 1359; i++){
      if(Run16AuAu_IdNo[i] == runid) {
        return i;
      }
    }
  }

  cout<<" *********** RunID not matched with list ************!!!! "<<endl;
  return -999;
}
//
// Fill event plane resolution histograms
//_____________________________________________________________________________
void StMyAnalysisMaker3::CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1)
{ 
    // fill the profiles for the resolution parameters
    // R2 resolution for 2nd order event plane
    fProfV2Resolution[ref9]->Fill(2., TMath::Cos(2.*(bbc - tpc)));
    fProfV2Resolution[ref9]->Fill(3., TMath::Cos(2.*(bbc - zdc)));
    fProfV2Resolution[ref9]->Fill(4., TMath::Cos(2.*(tpc - bbc)));  // bin2
    fProfV2Resolution[ref9]->Fill(5., TMath::Cos(2.*(tpc - zdc)));  
    fProfV2Resolution[ref9]->Fill(6., TMath::Cos(2.*(zdc - tpc)));  // bin6
    fProfV2Resolution[ref9]->Fill(7., TMath::Cos(2.*(zdc - bbc)));  // bin3
    fProfV2Resolution[ref9]->Fill(8., TMath::Cos(2.*(bbc - tpcN)));
    fProfV2Resolution[ref9]->Fill(9., TMath::Cos(2.*(bbc - tpcP)));
    fProfV2Resolution[ref9]->Fill(10., TMath::Cos(2.*(zdc - tpcN)));
    fProfV2Resolution[ref9]->Fill(11., TMath::Cos(2.*(zdc - tpcP)));
    fProfV2Resolution[ref9]->Fill(12., TMath::Cos(2.*(tpcP - tpcN)));
    fProfV2Resolution[ref9]->Fill(17., TMath::Cos(2.*(bbc1 - tpc)));
    fProfV2Resolution[ref9]->Fill(18., TMath::Cos(2.*(bbc1 - tpcN)));
    fProfV2Resolution[ref9]->Fill(19., TMath::Cos(2.*(bbc1 - tpcP)));
    fProfV2Resolution[ref9]->Fill(20., TMath::Cos(2.*(bbc1 - zdc1)));
    fProfV2Resolution[ref9]->Fill(21., TMath::Cos(2.*(zdc1 - tpc)));
    fProfV2Resolution[ref9]->Fill(22., TMath::Cos(2.*(zdc1 - tpcN)));
    fProfV2Resolution[ref9]->Fill(23., TMath::Cos(2.*(zdc1 - tpcP)));

    // R3 resolution for 2nd order event plane
    fProfV3Resolution[ref9]->Fill(2., TMath::Cos(3.*(bbc - tpc)));
    fProfV3Resolution[ref9]->Fill(3., TMath::Cos(3.*(bbc - zdc)));
    fProfV3Resolution[ref9]->Fill(4., TMath::Cos(3.*(tpc - bbc)));
    fProfV3Resolution[ref9]->Fill(5., TMath::Cos(3.*(tpc - zdc)));
    fProfV3Resolution[ref9]->Fill(6., TMath::Cos(3.*(zdc - tpc)));
    fProfV3Resolution[ref9]->Fill(7., TMath::Cos(3.*(zdc - bbc)));
    fProfV3Resolution[ref9]->Fill(8., TMath::Cos(3.*(bbc - tpcN)));
    fProfV3Resolution[ref9]->Fill(9., TMath::Cos(3.*(bbc - tpcP)));
    fProfV3Resolution[ref9]->Fill(10., TMath::Cos(3.*(zdc - tpcN)));
    fProfV3Resolution[ref9]->Fill(11., TMath::Cos(3.*(zdc - tpcP)));
    fProfV3Resolution[ref9]->Fill(12., TMath::Cos(3.*(tpcP - tpcN)));
    fProfV3Resolution[ref9]->Fill(17., TMath::Cos(3.*(bbc1 - tpc)));
    fProfV3Resolution[ref9]->Fill(18., TMath::Cos(3.*(bbc1 - tpcN)));
    fProfV3Resolution[ref9]->Fill(19., TMath::Cos(3.*(bbc1 - tpcP)));
    fProfV3Resolution[ref9]->Fill(20., TMath::Cos(3.*(bbc1 - zdc1)));
    fProfV3Resolution[ref9]->Fill(21., TMath::Cos(3.*(zdc1 - tpc)));
    fProfV3Resolution[ref9]->Fill(22., TMath::Cos(3.*(zdc1 - tpcN)));
    fProfV3Resolution[ref9]->Fill(23., TMath::Cos(3.*(zdc1 - tpcP)));

    // R4 resolution for 2nd order event plane
    fProfV4Resolution[ref9]->Fill(2., TMath::Cos(4.*(bbc - tpc)));
    fProfV4Resolution[ref9]->Fill(3., TMath::Cos(4.*(bbc - zdc)));
    fProfV4Resolution[ref9]->Fill(4., TMath::Cos(4.*(tpc - bbc)));
    fProfV4Resolution[ref9]->Fill(5., TMath::Cos(4.*(tpc - zdc)));
    fProfV4Resolution[ref9]->Fill(6., TMath::Cos(4.*(zdc - tpc)));
    fProfV4Resolution[ref9]->Fill(7., TMath::Cos(4.*(zdc - bbc)));
    fProfV4Resolution[ref9]->Fill(8., TMath::Cos(4.*(bbc - tpcN)));
    fProfV4Resolution[ref9]->Fill(9., TMath::Cos(4.*(bbc - tpcP)));
    fProfV4Resolution[ref9]->Fill(10., TMath::Cos(4.*(zdc - tpcN)));
    fProfV4Resolution[ref9]->Fill(11., TMath::Cos(4.*(zdc - tpcP)));
    fProfV4Resolution[ref9]->Fill(12., TMath::Cos(4.*(tpcP - tpcN)));
    fProfV4Resolution[ref9]->Fill(17., TMath::Cos(4.*(bbc1 - tpc)));
    fProfV4Resolution[ref9]->Fill(18., TMath::Cos(4.*(bbc1 - tpcN)));
    fProfV4Resolution[ref9]->Fill(19., TMath::Cos(4.*(bbc1 - tpcP)));
    fProfV4Resolution[ref9]->Fill(20., TMath::Cos(4.*(bbc1 - zdc1)));
    fProfV4Resolution[ref9]->Fill(21., TMath::Cos(4.*(zdc1 - tpc)));
    fProfV4Resolution[ref9]->Fill(22., TMath::Cos(4.*(zdc1 - tpcN)));
    fProfV4Resolution[ref9]->Fill(23., TMath::Cos(4.*(zdc1 - tpcP)));

    // R5 resolution for 2nd order event plane
    fProfV5Resolution[ref9]->Fill(2., TMath::Cos(5.*(bbc - tpc)));
    fProfV5Resolution[ref9]->Fill(3., TMath::Cos(5.*(bbc - zdc)));
    fProfV5Resolution[ref9]->Fill(4., TMath::Cos(5.*(tpc - bbc)));
    fProfV5Resolution[ref9]->Fill(5., TMath::Cos(5.*(tpc - zdc)));
    fProfV5Resolution[ref9]->Fill(6., TMath::Cos(5.*(zdc - tpc)));
    fProfV5Resolution[ref9]->Fill(7., TMath::Cos(5.*(zdc - bbc)));
    fProfV5Resolution[ref9]->Fill(8., TMath::Cos(5.*(bbc - tpcN)));
    fProfV5Resolution[ref9]->Fill(9., TMath::Cos(5.*(bbc - tpcP)));
    fProfV5Resolution[ref9]->Fill(10., TMath::Cos(5.*(zdc - tpcN)));
    fProfV5Resolution[ref9]->Fill(11., TMath::Cos(5.*(zdc - tpcP)));
    fProfV5Resolution[ref9]->Fill(12., TMath::Cos(5.*(tpcP - tpcN)));
    fProfV5Resolution[ref9]->Fill(17., TMath::Cos(5.*(bbc1 - tpc)));
    fProfV5Resolution[ref9]->Fill(18., TMath::Cos(5.*(bbc1 - tpcN)));
    fProfV5Resolution[ref9]->Fill(19., TMath::Cos(5.*(bbc1 - tpcP)));
    fProfV5Resolution[ref9]->Fill(20., TMath::Cos(5.*(bbc1 - zdc1)));
    fProfV5Resolution[ref9]->Fill(21., TMath::Cos(5.*(zdc1 - tpc)));
    fProfV5Resolution[ref9]->Fill(22., TMath::Cos(5.*(zdc1 - tpcN)));
    fProfV5Resolution[ref9]->Fill(23., TMath::Cos(5.*(zdc1 - tpcP)));

    // for the resolution of the combined vzero event plane, use two tpc halves as uncorrelated subdetectors
} 
//
// function to calculate: event plane chi
//_____________________________________________________________________________
Double_t StMyAnalysisMaker3::CalculateEventPlaneChi(Double_t res) {
  // return chi for given resolution to combine event plane estimates from two subevents
  // see Phys. Rev. C no. CS6346 (http://arxiv.org/abs/nucl-ex/9805001)
  Double_t chi(2.), delta(1.), con((TMath::Sqrt(TMath::Pi()))/(2.*TMath::Sqrt(2)));
  for (Int_t i(0); i < 15; i++) {
    chi = ((con*chi*TMath::Exp(-chi*chi/4.)*(TMath::BesselI0(chi*chi/4.)+TMath::BesselI1(chi*chi/4.))) < res) ? chi + delta : chi - delta;
    delta = delta / 2.;
  }

  return chi;
}
//
// track QA function to fill some histograms with track information
//_____________________________________________________________________________
void StMyAnalysisMaker3::TrackQA()
{
  // get # of tracks
  int nTrack = mPicoDst->numberOfTracks();
  double pi = 1.0*TMath::Pi();

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
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // for heavy ion collisions, get EP and fill QA plots
    // get angle between tracks and event plane
    double dEPtrk = (!doppAnalysis) ? RelativeEPJET(phi, TPC_PSI2) : 0;
    //cout<<"dEPtrk = "<<dEPtrk<<"  phi = "<<phi<<"  EP2 = "<<TPC_PSI2<<endl;
    if((pt > 0.20) && (pt <= 0.5)) hdEPtrk[0]->Fill(dEPtrk);
    if((pt > 0.50) && (pt <= 1.0)) hdEPtrk[1]->Fill(dEPtrk);
    if((pt > 1.00) && (pt <= 1.5)) hdEPtrk[2]->Fill(dEPtrk);
    if((pt > 1.50) && (pt <= 2.0)) hdEPtrk[3]->Fill(dEPtrk);
    if((pt > 2.00) && (pt <= 20.)) hdEPtrk[4]->Fill(dEPtrk);

    // fill other track QA plots
    hTrackPhi[ref9]->Fill(phi);
    hTrackEta[ref9]->Fill(eta);
    hTrackPt[ref9]->Fill(pt);
    hTrackEtavsPhi->Fill(phi, eta);

  }
}
//
// function filling trigger arrays with booleans set to kTRUE for fired triggers
//_________________________________________________________________________
void StMyAnalysisMaker3::FillTowerTriggersArr() {
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
    // get pointer to trigger
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // emc trigger parameters
    int emcTrigID = emcTrig->id();

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();

    // set trigger level type
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigID] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigID] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigID] = kTRUE;

    //cout<<"i = "<<i<<"  EmcTrigID = "<<emcTrigID<<"  adc = "<<emcTrig->adc()<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3<<endl;
  }

}
//
// function to require that a jet constituent tower fired a HT trigger
//___________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::DidTowerConstituentFireTrigger(StJet *jet) {  
  // tower constituent fired trigger
  Bool_t mFiredTrigger = kFALSE;

  // loop over constituents towers
  for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
    int towerIndex = jet->ClusterAt(itow);
    // get tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towerIndex));
    if(!tow){ continue; }
    
    // tower ID: get from index of array shifted by +1
    int towID = -1;
    towID = towerIndex + 1;
    if(towID < 0) continue;

    // change flag to true if jet tower fired trigger
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towID]) mFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towID]) mFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towID]) mFiredTrigger = kTRUE;
  
  } // tower constituent loop

  return mFiredTrigger;
}  
//
// function to require that a jet constituent tower fired a HT trigger
//___________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::DidBadTowerFireTrigger() {
  // bad/dead tower fired trigger
  Bool_t mBadTowerFiredTrigger = kFALSE;

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // get tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tow) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: get from index of array shifted by +1
    int towID = itow + 1;
    if(towID < 0) continue; // double check these aren't still in the event list

    // check if tower is bad or dead - functions return kTRUE if ok and kFALSE if NOT dead 
    // isTowerOK(towID) = kTRUE if tower is OK and not dead
    // isTowerDead(towID) = kTRUE if tower is dead, = kFALSE if tower is NOT dead
    bool isTowOk = (IsTowerOK(towID) && !IsTowerDead(towID)); // FIXME

    // change flag to true if jet tower fired trigger
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towID] && !isTowOk) mBadTowerFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towID] && !isTowOk) mBadTowerFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towID] && !isTowOk) mBadTowerFiredTrigger = kTRUE;

    if(mBadTowerFiredTrigger == kTRUE) cout<<"Bad tower fired trigger, towID: "<<towID<<endl;
  } // tower constituent loop

  return mBadTowerFiredTrigger;
}
//
// function to calcuate delta R between a jet centroid and a track
//___________________________________________________________________________________________
Double_t StMyAnalysisMaker3::GetDeltaR(StJet *jet, StPicoTrack *trk) {
  // constants
  double deltaR = -99.;
  double pi = 1.0*TMath::Pi();

  // get track momentum vector 
  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    if(!(trk->isPrimary())) return -99.; // check if primary
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
// ===========================================================================================
// ===========================================================================================
//
// function that does jet shape analysis
// 	refCorr2 is refCorr2 for AuAu, but is grefMult for pp
//___________________________________________________________________________________________
Int_t StMyAnalysisMaker3::JetShapeAnalysis(StJet *jet, StEventPool *pool, Double_t refCorr2, Int_t assocPtBin) {
    // constants
    double pi = 1.0*TMath::Pi();
    double rbinSize = 0.05;

    // get centrality bin
    int centBin = (!doppAnalysis) ? Get4CentBin(fCentralityScaled) : 0; // this secures it will still run for pp, and use cent bin = 0
    if(centBin < 0) return kStOK;

    // get jet info
    double jetPhi = jet->Phi();
    double jetEta = jet->Eta();

    // get jet pt bin and value
    double jetPt = -99.;
    if(fCorrJetPt) { jetPt = jet->Pt() - jet->Area()*fRhoVal;
    } else { jetPt = jet->Pt(); }
    int jetPtBin = GetJetPtBin(jetPt);
    if(jetPtBin < 0) return kStOK;

    // require tower and or track bias for jet
    //if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerE() < fTowerBias)) return kStOK;

    // check that jet contains a tower that fired the trigger
    //if(!DidTowerConstituentFireTrigger(jet)) return kStOK;

    // ==========================================================================================================
    // get StEventPlaneMaker from event
    StEventPlaneMaker *EventPlaneMaker[4];

    // pt-dependent bin mode
    const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
    for(int i=0; i<5; i++) {
      EventPlaneMaker[i] = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, i)));
    }

    // event plane bin to use
    int ptAssocBins[9] = {0, 1, 2, 3, 4, 4,4,4,4};
    int EPBinToUse = ptAssocBins[assocPtBin];

    // check for requested EventPlaneMaker pointer
    if(!EventPlaneMaker[EPBinToUse]) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn;}

    // get event plane angle for different pt bins
    // assign global event plane to selected pt-dependent bin
    // could also write this as:  tpc2EP_bin = (EventPlaneMaker) ? (double)EventPlaneMaker->GetTPCEP() : -999;
    double tpc2EP = (EventPlaneMaker[EPBinToUse]) ? (double)EventPlaneMaker[EPBinToUse]->GetTPCEP() : -999;
    double jetV2EP = (EventPlaneMaker[EPBinToUse]) ? (double)EventPlaneMaker[EPBinToUse]->GetTPCEP() : -999;

    // if requiring a single event plane angle (non-pt dependent): use charged tracks 0.2-2.0 GeV to calculate EP
    if(doUseMainEPAngle) {
      tpc2EP = (EventPlaneMaker[4]) ? (double)EventPlaneMaker[4]->GetTPCEP() : -999;
      jetV2EP = (EventPlaneMaker[4]) ? (double)EventPlaneMaker[4]->GetTPCEP() : -999;
    }

    //cout<<"assocPtBin: "<<assocPtBin<<"  tpc2EP: "<<tpc2EP<<"  jetV2EP: "<<jetV2EP<<endl;

    // get relative angle between jet and event plane
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, tpc2EP) : -99.; // CORRECTED event plane angle - STEP3

    // get relative jet-event plane bin
    int EPBin = (!doppAnalysis) ? GetJetEPBin(dEP) : 0; // this secures it will still run for pp, and use EP bin = 0
    if(EPBin < 0) return kStOK;

    // calculate jet v2 here, can move around for other analysis, but here for now
    // jetPtr, EPangle, ptBin
    // only calculate and fill jet v2 for assocPtBin == 0 for now - TODO
    //    if(assocPtBin == 0) GetJetV2(jet, tpc2EP, 1);
    //if(assocPtBin == 0 ) GetJetV2(jet, jetV2EP, 1); 
    GetJetV2(jet, jetV2EP, assocPtBin);
    // ==========================================================================================================================

    // annuli sum - initialize
    double rsum[10] = {0.0};
    double rsumBG[10] = {0.0};
    //double rsumBG3[10] = {0.0};

    // track loop inside jet loop - loop over ALL tracks in PicoDst
    Int_t ntracks = mPicoDst->numberOfTracks();
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

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
      double tphi = mTrkMom.Phi();
      if(tphi < 0.0)    tphi += 2.0*pi;   // require 0,2pi interval
      if(tphi > 2.0*pi) tphi -= 2.0*pi;
      double teta = mTrkMom.PseudoRapidity();

      // cut on track pt
      if(tpt < fJetShapeTrackPtMin) { continue; }
      if(tpt > fJetShapeTrackPtMax) { continue; }

      // additional pt selection when doing pt associated bin method
      if(doTPCptassocBin) {
        if(assocPtBin == 0) { if((tpt < 0.20) || (tpt  > 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(assocPtBin == 1) { if((tpt < 0.50) || (tpt  > 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(assocPtBin == 2) { if((tpt < 1.00) || (tpt  > 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(assocPtBin == 3) { if((tpt < 1.50) || (tpt  > 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(assocPtBin == 4) { if((tpt < 2.00) || (tpt  > 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(assocPtBin == 5) { if((tpt < 3.00) || (tpt  > 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(assocPtBin == 6) { if((tpt < 4.00) || (tpt  > 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations, 8->6 Oct23
        if(assocPtBin == 7) { if((tpt <  6.0))                 continue; }  //       6.0+ GeV assoc bin used for correlations, 8->6 Oct23
        if(assocPtBin == 8) { if((tpt <  0.5))                 continue; }  //       0.5+ GeV assoc bin used for correlations
      }

      // get radial distance between track and jet axis
      double deltaR = GetDeltaR(jet, trk);

      // get annuli bin
      int annuliBin = GetAnnuliBin(deltaR);
      if(annuliBin < 0) continue;

      // calculate radial pt sum
      rsum[annuliBin] += tpt;

    } // track loop

    // QA histograms
    // fill only for case of kJetShapePtAssocBin == 0 (assocPtBin)
    if(assocPtBin == 0) {
      hHTvsMult->Fill(refCorr2);
      hTriggerEvtStatZVtx->Fill(zVtx);
      hTriggerEvtStatCent->Fill(fCentralityScaled); // bad for pp
      hTriggerEvtStatZvsCent->Fill(fCentralityScaled, zVtx); // bad for pp
    }

    // fill jet shape histograms
    for(int i = 0; i < 10; i++) {
      hJetShape[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);
      hJetShape[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);

      hJetPtProfile[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
      hJetPtProfile[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
    }
    if(assocPtBin == 0) {
      hJetCounter[jetPtBin][centBin][EPBin]->Fill(0.5);
      hJetCounter[jetPtBin][centBin][3]->Fill(0.5); // ALL angles
    }

    // background calculation variables
    bool case1 = kFALSE, case2 = kFALSE, case3 = kFALSE;
    double jetPhiBG = -99., jetEtaBG = -99.;
    double jetPhiBG3 = -99., jetEtaBG3 = -99.;

    // fiducial cuts
    double etaMin = fJetRad;
    double etaMax = 1.0 - fJetRad;

    // CASE 1: Eta reflection
    if(TMath::Abs(jetEta) > etaMin && TMath::Abs(jetEta) < etaMax) {
      for(int i = 0; i < 10; i++) {
        hJetShapeCase1[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);
        hJetShapeCase1[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);

        hJetPtProfileCase1[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
        hJetPtProfileCase1[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
      }
      if(assocPtBin == 0){
        hJetCounterCase1[jetPtBin][centBin][EPBin]->Fill(0.5);
        hJetCounterCase1[jetPtBin][centBin][3]->Fill(0.5); // ALL angles 
      }

      // background
      case1 = kTRUE;
      jetEtaBG = - jetEta;
      jetPhiBG =   jetPhi;
    }

    // CASE 2: Phi shifted
    if(TMath::Abs(jetEta) < etaMin) {
      for(int i = 0; i < 10; i++) {
        hJetShapeCase2[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);
        hJetShapeCase2[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);

        hJetPtProfileCase2[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
        hJetPtProfileCase2[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
      }
      if(assocPtBin == 0) {
        hJetCounterCase2[jetPtBin][centBin][EPBin]->Fill(0.5);
        hJetCounterCase2[jetPtBin][centBin][3]->Fill(0.5); // ALL angles
      }

      // background
      case2 = kTRUE;
      jetEtaBG = jetEta;
      jetPhiBG = jetPhi + 0.5*pi; // 90 degree shift
      if(jetPhiBG > 2.0*pi) jetPhiBG -= 2.0*pi;
    }

    // CASE 3: Background eta is calculated based on where jet is
    if(TMath::Abs(jetEta) < etaMax) {
      // background
      if(fJetRad < 0.4) { // for R=0.3 and R=0.2 jets
        case3 = kTRUE;

        // what about when jetEta is > etaMax? - FIXME
        if(jetEta > -0.3 && jetEta < -0.2) jetEtaBG3 = jetEta + 2.0*fJetRad + 0.20;
        if(jetEta > -0.2 && jetEta < -0.1) jetEtaBG3 = jetEta + 2.0*fJetRad + 0.20;
        if(jetEta > -0.1 && jetEta <  0.0) jetEtaBG3 = jetEta + 2.0*fJetRad + 0.10;
        if(jetEta >  0.0 && jetEta <  0.1) jetEtaBG3 = jetEta - 2.0*fJetRad - 0.10;
        if(jetEta >  0.1 && jetEta <  0.2) jetEtaBG3 = jetEta - 2.0*fJetRad - 0.20;
        if(jetEta >  0.2 && jetEta <  0.3) jetEtaBG3 = jetEta - 2.0*fJetRad - 0.20;

        jetPhiBG3 = jetPhi;
      }
    }

    // BACKGROUND tracks
    // track loop inside jet loop - loop over ALL tracks in PicoDst - for BG
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack* trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

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
      double tphi = mTrkMom.Phi();
      if(tphi < 0.0)    tphi += 2.0*pi; // require 0,2pi interval
      if(tphi > 2.0*pi) tphi -= 2.0*pi;
      double teta = mTrkMom.PseudoRapidity();

      // cut on track pt
      if(tpt < fJetShapeTrackPtMin) { continue; }
      if(tpt > fJetShapeTrackPtMax) { continue; }

      // additional pt selection when doing pt associated bin method
      if(doTPCptassocBin) {
        if(assocPtBin == 0) { if((tpt < 0.20) || (tpt  > 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(assocPtBin == 1) { if((tpt < 0.50) || (tpt  > 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(assocPtBin == 2) { if((tpt < 1.00) || (tpt  > 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(assocPtBin == 3) { if((tpt < 1.50) || (tpt  > 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(assocPtBin == 4) { if((tpt < 2.00) || (tpt  > 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(assocPtBin == 5) { if((tpt < 3.00) || (tpt  > 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(assocPtBin == 6) { if((tpt < 4.00) || (tpt  > 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations, 8->6 Oct23
        if(assocPtBin == 7) { if((tpt <  6.0))                 continue; }  //       6.0+ GeV assoc bin used for correlations, 8->6 Oct23
        if(assocPtBin == 8) { if((tpt <  0.5))                 continue; }  //       0.5+ GeV assoc bin used for correlations
      }

      // get radial distance between track and jet axis
      double deltaEtaBG = 1.0*TMath::Abs(jetEtaBG - teta);
      double deltaPhiBG = 1.0*TMath::Abs(jetPhiBG - tphi);
      if(deltaPhiBG > 1.0*pi) deltaPhiBG = 2.0*pi - deltaPhiBG;
      double deltaR = 1.0*TMath::Sqrt(deltaEtaBG*deltaEtaBG + deltaPhiBG*deltaPhiBG);

      // get annuli bin
      int annuliBin = GetAnnuliBin(deltaR);
      //if(annuliBin < 0) continue;

      // calculate radial pt sum
      //rsumBG[annuliBin] += tpt;
      if(annuliBin >= 0) rsumBG[annuliBin] += tpt; // use this instance to avoid 2 separate track loops

/*
      // =============================================================================
      // this is for CASE 3, to keep separate
      // get radial distance between track and jet axis
      double deltaEtaBG3 = 1.0*TMath::Abs(jetEtaBG3 - teta);
      double deltaPhiBG3 = 1.0*TMath::Abs(jetPhiBG3 - tphi);
      if(deltaPhiBG3 > 1.0*pi) deltaPhiBG3 = 2.0*pi - deltaPhiBG3;
      double deltaR3 = 1.0*TMath::Sqrt(deltaEtaBG3*deltaEtaBG3 + deltaPhiBG3*deltaPhiBG3);

      // get annuli bin
      int annuliBin3 = GetAnnuliBin(deltaR3);
      //if(annuliBin3 < 0) continue;

      // calculate radial pt sum
      //rsumBG3[annuliBin3] += tpt;
      if(annuliBin3 >= 0) rsumBG3[annuliBin3] += tpt; // use this instance to avoid 2 separate track loops
*/
      // ==============================================================================
    } // track loop

    // inclusive case: Background
    for(int i = 0; i < 10; i++) {
      hJetShapeBG[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);
      hJetShapeBG[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);

      hJetPtProfileBG[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
      hJetPtProfileBG[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
    }

    // Case 1: background
    if(case1) {
      for(int i = 0; i < 10; i++) {
        hJetShapeBGCase1[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);
        hJetShapeBGCase1[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);

        hJetPtProfileBGCase1[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
        hJetPtProfileBGCase1[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
      }
    }

    // Case 2: background
    if(case2) {
      for(int i = 0; i < 10; i++) {
        hJetShapeBGCase2[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);
        hJetShapeBGCase2[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);

        hJetPtProfileBGCase2[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
        hJetPtProfileBGCase2[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
      }
    }

/*
    // Case 3: background
    if(case3) {
      for(int i = 0; i < 10; i++) {
        hJetShapeBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]/jetPt);
        hJetShapeBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]/jetPt);

        hJetPtProfileBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]);
        hJetPtProfileBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]);
      }
    }
*/

    // event mixing for background jet cones
    if(fDoEventMixing > 0){
      // initialize background tracks array
      TObjArray *bgTracks;

      // do event mixing when Signal Jet is part of event with a HT1 or HT2 or HT3 trigger firing
      if(pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {
        // get number of current events in pool
        int nMix = pool->GetCurrentNEvents();

        // QA histogram
        if(assocPtBin == 0) hNMixEvents->Fill(nMix);

        // reset annuli sums here - do this inside loop over mixed events
        //double rsumBG3[10] = {0.0};

        // Fill mixed-event histos here: loop over nMix events
        for(int jMix = 0; jMix < nMix; jMix++) {
          // get jMix'th event
          bgTracks = pool->GetEvent(jMix);
          const Int_t Nbgtrks = bgTracks->GetEntries();

          // reset annuli sums here - when NOT normalizing by nMix
          double rsumBG3[10] = {0.0};

          // loop over background (mixed event) tracks
          for(int ibg = 0; ibg < Nbgtrks; ibg++) {
            // trying new slimmed PicoTrack class: StFemtoTrack
            StFemtoTrack* trk = static_cast<StFemtoTrack*>(bgTracks->At(ibg));
            if(!trk) continue;

            double Mphi = trk->Phi();
            double Meta = trk->Eta();
            double Mpt = trk->Pt();

            // shift angle (0, 2*pi) 
            if(Mphi < 0.0)    Mphi += 2.0*pi;
            if(Mphi > 2.0*pi) Mphi -= 2.0*pi;

            // cut on track pt
            if(Mpt < fJetShapeTrackPtMin) { continue; }
            if(Mpt > fJetShapeTrackPtMax) { continue; }

            // additional pt selection when doing pt associated bin method
            if(doTPCptassocBin) {
              if(assocPtBin == 0) { if((Mpt < 0.20) || (Mpt  > 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
              if(assocPtBin == 1) { if((Mpt < 0.50) || (Mpt  > 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
              if(assocPtBin == 2) { if((Mpt < 1.00) || (Mpt  > 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
              if(assocPtBin == 3) { if((Mpt < 1.50) || (Mpt  > 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
              if(assocPtBin == 4) { if((Mpt < 2.00) || (Mpt  > 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
              if(assocPtBin == 5) { if((Mpt < 3.00) || (Mpt  > 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
              if(assocPtBin == 6) { if((Mpt < 4.00) || (Mpt  > 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations, 8->6 Oct23
              if(assocPtBin == 7) { if((Mpt <  6.0))                 continue; }  //       6.0+ GeV assoc bin used for correlations, 8->6 Oct23
              if(assocPtBin == 8) { if((Mpt <  0.5))                 continue; }  //       0.5+ GeV assoc bin used for correlations
            }

            // get radial distance between track and jet axis
            double deltaEtaBG = 1.0*TMath::Abs(jetEta - Meta);
            double deltaPhiBG = 1.0*TMath::Abs(jetPhi - Mphi);
            if(deltaPhiBG > 1.0*pi) deltaPhiBG = 2.0*pi - deltaPhiBG;
            double deltaR = 1.0*TMath::Sqrt(deltaEtaBG*deltaEtaBG + deltaPhiBG*deltaPhiBG);

            // get annuli bin
            int annuliBin = GetAnnuliBin(deltaR);
            if(annuliBin < 0) continue;

            // calculate radial pt sum
            rsumBG3[annuliBin] += Mpt;

            // calculate single particle tracking efficiency of mixed events for correlations (-999)
            //double mixefficiency = 1.0;
            //FIXME mixefficiency = EffCorrection(part->Eta(), part->Pt(), fDoEffCorr);                           
          } // end of background track loop

          // fill BG histos here
          for(int i = 0; i < 10; i++) {
            hJetShapeBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
            hJetShapeBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
            hJetPtProfileBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
            hJetPtProfileBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
          } // loop over annuli bins
          if(assocPtBin == 0) { 
            //hJetCounterCase3BG[jetPtBin][centBin][EPBin]->Fill(0.5);
            //hJetCounterCase3BG[jetPtBin][centBin][3]->Fill(0.5); // ALL angles
          }

        }   // end of filling mixed-event histo's:  jth mix event loop

/*  ------ October 29, 2018
        // fill histo here... FIXME
        for(int i = 0; i < 10; i++) {
          hJetShapeBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
          hJetShapeBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
          hJetPtProfileBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
          hJetPtProfileBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
        } // loop over annuli bins
        if(assocPtBin == 0) {
          hJetCounterCase3BG[jetPtBin][centBin][EPBin]->Fill(0.5);
          hJetCounterCase3BG[jetPtBin][centBin][3]->Fill(0.5); // ALL angles
        }
*/

        if(assocPtBin == 0) {
          hJetCounterCase3BG[jetPtBin][centBin][EPBin]->Fill(0.5);
          hJetCounterCase3BG[jetPtBin][centBin][3]->Fill(0.5); // ALL angles
        }
      }     // end of check for pool being ready
    }       // end of event mixing

    return kStOK;
}
//
// This is strictly set up for a test
//___________________________________________________________________________________________
Double_t StMyAnalysisMaker3::TestBool() {
  return kFALSE;
}
//
// Function to check if Tower is OK or NOT
//____________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::IsTowerOK( Int_t mTowId ){
  //if( badTowers.size()==0 ){
  if( badTowers.empty() ){
    __ERROR("StMyAnalysisMaker3::IsTowerOK: WARNING: You're trying to run without a bad tower list. If you know what you're doing, deactivate this throw and recompile.");
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
// Function to check if Tower is DEAD or NOT
//____________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::IsTowerDead( Int_t mTowId ){
  //if( deadTowers.size()==0 ){
  if( deadTowers.empty() ){
    __ERROR("StMyAnalysisMaker3::IsTowerDead: WARNING: You're trying to run without a dead tower list. If you know what you're doing, deactivate this throw and recompile.");
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
// function to calculate Jet v2 using the event plane method
//____________________________________________________________________________________________
void StMyAnalysisMaker3::GetJetV2(StJet *jet, Double_t EPangle, Int_t ptAssocBin)
{
  // get centrality bin
  int centBin = (!doppAnalysis) ? Get4CentBin(fCentralityScaled) : 0; // this secures it will still run for pp
  if(centBin < 0) return; // kStOK;

  // get jet pt bin and value
  double jetPt = -99.;
  if(fCorrJetPt) { jetPt = jet->Pt() - jet->Area()*fRhoVal;
  } else { jetPt = jet->Pt(); }
  int jetPtBin = GetJetPtBin(jetPt);
  if(jetPtBin < 0) return; // kStOK;

  // require tower and or track bias for jet
  //if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerE() < fTowerBias)) return kStOK;

  // check that jet contains a tower that fired the trigger
  //if(!DidTowerConstituentFireTrigger(jet)) return kStOK;

  // jet parameters
  double jetPhi = jet->Phi();

  // get event plane bin
  double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, EPangle) : -99.; // CORRECTED event plane angle - STEP3
  int EPBin = (!doppAnalysis) ? GetJetEPBin(dEP) : 0; // this secures it will still run for pp, and use EP bin = 0
  if(EPBin < 0) return; // kStOK;

  // fill histogram
  fProfJetV2[jetPtBin][centBin][EPBin]->Fill(ptAssocBin, TMath::Cos(2.*(jetPhi - EPangle)));
  fProfJetV2[jetPtBin][centBin][3]->Fill(ptAssocBin, TMath::Cos(2.*(jetPhi - EPangle)));

}
//
// Function to add event pools to output file
//____________________________________________________________________
void StMyAnalysisMaker3::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPsi2, Double_t maxPsi2, Double_t minPt, Double_t maxPt)
{
  // create vector of vectors
  std::vector<Double_t> binVec;
  binVec.push_back(minCent);
  binVec.push_back(maxCent);
  binVec.push_back(minZvtx);
  binVec.push_back(maxZvtx);
  binVec.push_back(minPsi2);
  binVec.push_back(maxPsi2);
  binVec.push_back(minPt);
  binVec.push_back(maxPt);

  // add vectors to output list
  fEventPoolOutputList.push_back(binVec);
}
