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
//      - check on tower constituents firing trigger threshold (HT, JP)
//      - check for bad runs, bad towers, dead towers,
//      - 
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
#include "TF2.h"
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
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"

// jet-framework includes
#include "StEventPlaneMaker.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPoolManager.h"
#include "StFemtoTrack.h"
#include "StCentMaker.h"
//#include "trackingEfficiency_Run14.h"

// old file kept
#include "StPicoConstants.h"

ClassImp(StMyAnalysisMaker3)

//______________________________________________________________________________
StMyAnalysisMaker3::StMyAnalysisMaker3(const char *name, StPicoDstMaker *picoMaker, const char *outName = "", bool mDoComments = kFALSE, double minJetPt = 1.0, const char *jetMakerName = "", const char *rhoMakerName = "")
  : StJetFrameworkPicoBase(name)
{
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  doJetShapeAnalysis = kFALSE;
  fJetAnalysisJetType = kLeadingJets;  // see header, LJ, SubLJ, inclusive
  doRequireAjSelection = kFALSE;
  fCorrJetPt = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doWriteTrackQAHist = kTRUE;
  doWriteJetQAHist = kTRUE;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  fLeadingJet = 0x0; fSubLeadingJet = 0x0; fExcludeLeadingJetsFromFit = 1.0;
  fTrackWeight = 2; //StJetFrameworkPicoBase::kPtLinear2Const5Weight; // see StJetFrameworkPicoBase::EPtrackWeightType 
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
  mOutName = outName;
  mOutNameEP = "";
  mOutNameQA = "";
  mOutNameME = "";
  doPrintEventCounter = kFALSE;
  fDoEffCorr = kFALSE;
  fTrackEfficiencyType = StJetFrameworkPicoBase::kNormalPtEtaBased;
  doEventPlaneRes = kFALSE;
  doTPCptassocBin = kFALSE;
  fTPCptAssocBin = -99;
  doUseMainEPAngle = kFALSE;  // kTRUE: use 0.2-2.0 GeV charged tracks for event plane
  doRejectBadRuns = kFALSE;
  fMinPtJet = minJetPt;
  fJetConstituentCut = 2.0;
  fTrackBias = 0.2;
  fTowerBias = 0.2;
  fJetRad = 0.4;
  fJetShapeTrackPtMin = 0.2;  fJetShapeTrackPtMax = 30.0;
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
  fCentBinSize = 5; fReduceStatsCent = -1;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
  //mVertex = 0x0;
  zVtx = 0.0;
  fDoFilterPtMixEvents = kFALSE;
  fDoUseMultBins = kFALSE;
  doUseEPBins = kFALSE;
  fnEPBins = 6;
  doIgnoreExternalME = kTRUE;
  fBackgroundConeFractionCut = 5.0; // high number by default to not use
  doGenerateBadMixEventBGcone = kFALSE;
  doRunAnalysis = kTRUE;
  doJetHadronCorrelationAnalysis = kFALSE;
  doSkip1ParticleJets = kFALSE;
  doBiasJetLeadConstituent = kFALSE;
  doRequireJetTowFireTrig = kFALSE;
  fSysUncType = kDoNothing;
  fEmcTriggerEventType = 0; fMBEventType = 2;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4800; i++) {
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
  mCentMaker = 0x0;
  mBaseMaker = 0x0;
  //fEventPoolOutputList = {{}}; // FIXME - does this make sense as is?
  fUsePtBinnedEventPool = kFALSE;
  fCheckEventNumberInMixedEvent = kFALSE;
  fListOfPools = 0x0;
  fEfficiencyInputFile = 0x0;
}
//
//
//__________________________________________________________________________________
StMyAnalysisMaker3::~StMyAnalysisMaker3()
{ /*  */
  // destructor
  if(hdEPReactionPlaneFnc) delete hdEPReactionPlaneFnc;
  if(hEventPlaneFncN2)     delete hEventPlaneFncN2;
  if(hEventPlaneFncP2)     delete hEventPlaneFncP2;
  if(hEventPlaneFnc2)      delete hEventPlaneFnc2;
  if(hEventPlaneClass)     delete hEventPlaneClass;

  if(hEventPlane)   delete hEventPlane;
  if(fHistEPTPCn)   delete fHistEPTPCn;
  if(fHistEPTPCp)   delete fHistEPTPCp;
  if(fHistEPBBC)    delete fHistEPBBC;
  if(fHistEPZDC)    delete fHistEPZDC;
  if(hEventZVertex) delete hEventZVertex;
  if(hCentrality)   delete hCentrality;
  if(hCentralityPostCut) delete hCentralityPostCut;
  if(hMultiplicity) delete hMultiplicity;
  if(hStats)        delete hStats;
  if(hRhovsCent)    delete hRhovsCent;
  for(int i=0; i<5; i++) { 
    if(hdEPtrk[i])  delete hdEPtrk[i]; 
  }
  for(int i=0; i<9; i++){ // centrality
    if(hTrackPhi[i]) delete hTrackPhi[i];
    if(hTrackEta[i]) delete hTrackEta[i];
    if(hTrackPt[i])  delete hTrackPt[i];
  }
  if(hTrackEtavsPhi) delete hTrackEtavsPhi;

  if(hJetPt)             delete hJetPt;
  if(hJetCorrPt)         delete hJetCorrPt;
  if(hJetLeadingPt)      delete hJetLeadingPt;
  if(hJetSubLeadingPt)   delete hJetSubLeadingPt;
  if(hJetLeadingPtAj)    delete hJetLeadingPtAj;
  if(hJetSubLeadingPtAj) delete hJetSubLeadingPtAj;
  if(hJetDiJetAj)        delete hJetDiJetAj;
  if(hJetE)              delete hJetE;
  if(hJetEta)            delete hJetEta;
  if(hJetPhi)            delete hJetPhi;
  if(hJetNEF)            delete hJetNEF;
  if(hJetArea)           delete hJetArea;
  if(hJetMass)           delete hJetMass;
  if(hJetTracksPt)       delete hJetTracksPt;
  if(hJetTracksPhi)      delete hJetTracksPhi;
  if(hJetTracksEta)      delete hJetTracksEta;
  if(hJetTracksZ)        delete hJetTracksZ;
  if(hJetPtvsArea)       delete hJetPtvsArea;

  for(int i = 0; i < 5; i++) {
    if(hJetEventEP[i])    delete hJetEventEP[i];
    if(hJetPhivsEP[i])    delete hJetPhivsEP[i];

    if(hJetPtIn[i])       delete hJetPtIn[i];
    if(hJetPhiIn[i])      delete hJetPhiIn[i];
    if(hJetEtaIn[i])      delete hJetEtaIn[i];
    if(hJetEventEPIn[i])  delete hJetEventEPIn[i];
    if(hJetPhivsEPIn[i])  delete hJetPhivsEPIn[i];
    if(hJetPtMid[i])      delete hJetPtMid[i];
    if(hJetPhiMid[i])     delete hJetPhiMid[i];
    if(hJetEtaMid[i])     delete hJetEtaMid[i];
    if(hJetEventEPMid[i]) delete hJetEventEPMid[i];
    if(hJetPhivsEPMid[i]) delete hJetPhivsEPMid[i];
    if(hJetPtOut[i])      delete hJetPtOut[i];
    if(hJetPhiOut[i])     delete hJetPhiOut[i];
    if(hJetEtaOut[i])     delete hJetEtaOut[i];
    if(hJetEventEPOut[i]) delete hJetEventEPOut[i];
    if(hJetPhivsEPOut[i]) delete hJetPhivsEPOut[i];
  }

  if(hJetHTrigMaxTowEt)      delete hJetHTrigMaxTowEt;
  if(hJetHTrigMaxTrkPt)      delete hJetHTrigMaxTrkPt;
  if(fHistJetHEtaPhi)        delete fHistJetHEtaPhi;
  if(fHistEventSelectionQA)  delete fHistEventSelectionQA;
  if(fHistEventSelectionQAafterCuts) delete fHistEventSelectionQAafterCuts;
  if(hEmcTriggers)           delete hEmcTriggers;
  if(hEventTriggerIDs)       delete hEventTriggerIDs;
  if(hBadTowerFiredTrigger)  delete hBadTowerFiredTrigger;
  if(hNGoodTowersFiringTrigger)delete hNGoodTowersFiringTrigger;
  if(hMixEvtStatZVtx)        delete hMixEvtStatZVtx;
  if(hMixEvtStatCent)        delete hMixEvtStatCent;
  if(hMixEvtStatZvsCent)     delete hMixEvtStatZvsCent;
  if(hTriggerEvtStatZVtx)    delete hTriggerEvtStatZVtx;
  if(hTriggerEvtStatCent)    delete hTriggerEvtStatCent;
  if(hTriggerEvtStatZvsCent) delete hTriggerEvtStatZvsCent;
  if(hMBvsMult)              delete hMBvsMult;
  if(hMB5vsMult)             delete hMB5vsMult;
  if(hMB30vsMult)            delete hMB30vsMult;
  if(hHTvsMult)              delete hHTvsMult;
  if(hNMixEvents)            delete hNMixEvents;

  if(hNEventsvsZvtxMB5)      delete hNEventsvsZvtxMB5;
  if(hNEventsvsCentMB5)      delete hNEventsvsCentMB5;
  if(hNEventsvsMultMB5)      delete hNEventsvsMultMB5;
  if(hNEventsvsZvsCentMB5)   delete hNEventsvsZvsCentMB5;
  if(hNEventsvsZvtxMB30)     delete hNEventsvsZvtxMB30;
  if(hNEventsvsCentMB30)     delete hNEventsvsCentMB30;
  if(hNEventsvsMultMB30)     delete hNEventsvsMultMB30;
  if(hNEventsvsZvsCentMB30)  delete hNEventsvsZvsCentMB30;
  if(hNEventsvsZvtxMB5Wt)      delete hNEventsvsZvtxMB5Wt;    
  if(hNEventsvsCentMB5Wt)      delete hNEventsvsCentMB5Wt;
  if(hNEventsvsMultMB5Wt)      delete hNEventsvsMultMB5Wt;
  if(hNEventsvsZvsCentMB5Wt)   delete hNEventsvsZvsCentMB5Wt;
  if(hNEventsvsZvtxMB30Wt)     delete hNEventsvsZvtxMB30Wt;
  if(hNEventsvsCentMB30Wt)     delete hNEventsvsCentMB30Wt;
  if(hNEventsvsMultMB30Wt)     delete hNEventsvsMultMB30Wt;
  if(hNEventsvsZvsCentMB30Wt)  delete hNEventsvsZvsCentMB30Wt;
  if(hNEventsvsZvtxHT2)        delete hNEventsvsZvtxHT2;
  if(hNEventsvsCentHT2)        delete hNEventsvsCentHT2;
  if(hNEventsvsMultHT2)        delete hNEventsvsMultHT2;
  if(hNEventsvsZvsCentHT2)     delete hNEventsvsZvsCentHT2;
  if(hNPairsvsZvtxMB5)         delete hNPairsvsZvtxMB5;
  if(hNPairsvsZvtxMB30)        delete hNPairsvsZvtxMB30;
  if(hNPairsvsZvtxHT2)         delete hNPairsvsZvtxHT2;
  if(hNPairsvsZvtxMB5Wt)       delete hNPairsvsZvtxMB5Wt;
  if(hNPairsvsZvtxMB30Wt)      delete hNPairsvsZvtxMB30Wt;

  for(int k = 0; k<9; k++) {
    if(hNMixNormBefore[k])   delete hNMixNormBefore[k];
    if(hNMixNormAfter[k])    delete hNMixNormAfter[k];
  }

  if(hBGconeFractionOfJetPt) delete hBGconeFractionOfJetPt;
  if(hJetPtvsBGconeFraction) delete hJetPtvsBGconeFraction;
  if(hJetPtvsBGconePt)       delete hJetPtvsBGconePt;

  if(hMB5TrkPtRaw)           delete hMB5TrkPtRaw;
  if(hMB5TrkPtReWeight)      delete hMB5TrkPtReWeight;
  if(hMB30TrkPtRaw)          delete hMB30TrkPtRaw;
  if(hMB30TrkPtReWeight)     delete hMB30TrkPtReWeight;

  if(hTPCvsBBCep) delete hTPCvsBBCep;
  if(hTPCvsZDCep) delete hTPCvsZDCep;
  if(hBBCvsZDCep) delete hBBCvsZDCep;

  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        if(fProfJetV2[k][j][i]) delete fProfJetV2[k][j][i];
      }
    }
  }

  if(doJetShapeAnalysis) {
    for(int k=0; k<4; k++) {
      for(int j=0; j<4; j++) {
        for(int i=0; i<4; i++) {
          if(hJetCounter[k][j][i])        delete hJetCounter[k][j][i];
          if(hJetCounterCase1[k][j][i])   delete hJetCounterCase1[k][j][i];
          if(hJetCounterCase2[k][j][i])   delete hJetCounterCase2[k][j][i];
          if(hJetCounterCase3BG[k][j][i]) delete hJetCounterCase3BG[k][j][i];

          for(int p=0; p<9; p++) {
            if(hJetShape[k][j][i][p])        delete hJetShape[k][j][i][p];
            if(hJetShapeCase1[k][j][i][p])   delete hJetShapeCase1[k][j][i][p];
            if(hJetShapeCase2[k][j][i][p])   delete hJetShapeCase2[k][j][i][p];
            if(hJetShapeBG[k][j][i][p])      delete hJetShapeBG[k][j][i][p];
            if(hJetShapeBGCase1[k][j][i][p]) delete hJetShapeBGCase1[k][j][i][p];
            if(hJetShapeBGCase2[k][j][i][p]) delete hJetShapeBGCase2[k][j][i][p];
            if(hJetShapeBGCase3[k][j][i][p]) delete hJetShapeBGCase3[k][j][i][p];

            if(hJetPtProfile[k][j][i][p])        delete hJetPtProfile[k][j][i][p];
            if(hJetPtProfileCase1[k][j][i][p])   delete hJetPtProfileCase1[k][j][i][p];
            if(hJetPtProfileCase2[k][j][i][p])   delete hJetPtProfileCase2[k][j][i][p];
            if(hJetPtProfileBG[k][j][i][p])      delete hJetPtProfileBG[k][j][i][p];
            if(hJetPtProfileBGCase1[k][j][i][p]) delete hJetPtProfileBGCase1[k][j][i][p];
            if(hJetPtProfileBGCase2[k][j][i][p]) delete hJetPtProfileBGCase2[k][j][i][p];
            if(hJetPtProfileBGCase3[k][j][i][p]) delete hJetPtProfileBGCase3[k][j][i][p];
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
  if(fhnJH)          delete fhnJH;
  if(fhnMixedEvents) delete fhnMixedEvents;
  if(fhnCorr)        delete fhnCorr;

  // clear and delete objects
//  fJets->Clear();    delete fJets;
//  fRho->Clear();     delete fRho; 
//  fPoolMgr->Clear(); delete fPoolMgr;

  // Clear unnecessary pools before saving - FIXME
  fPoolMgr->ClearPools();
  if(fListOfPools) delete fListOfPools;

  // track reconstruction efficiency input file
  if(fEfficiencyInputFile) {
    fEfficiencyInputFile->Close();
    delete fEfficiencyInputFile;
  }
}
//
// initialize objects & set up
//_________________________________________________________________________________________
Int_t StMyAnalysisMaker3::Init() {
  //StJetFrameworkPicoBase::Init();

  // input file - for tracking efficiency: Run14 AuAu
  //const char *input = Form("./StRoot/StMyAnalysisMaker/Run14_efficiency.root");
  ////const char *input = Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller.root");
  const char *input = "";
//if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) input=Form("./StRoot/StMyAnalysisMaker/Run14_AuAu_200_tracking_efficiency_and_momentum_smearing_dca_3p0_nhit_15_nhitfrac_0p52.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200)    input = Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) input = Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)      input = Form("./StRoot/StMyAnalysisMaker/Run12_efficiency_New.root"); // Oct17, 2019 added
  if(fDoEffCorr) {
    ///fEfficiencyInputFile = TFile::Open(input);
    //fEfficiencyInputFile = new TFile(input);
    fEfficiencyInputFile = new TFile(input, "READ");
    if(!fEfficiencyInputFile) cout<<Form("do not have input file: %s", input);
  }

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

  return kStOK;
}
//
//
//______________________________________________________________________________________
Int_t StMyAnalysisMaker3::Finish() { 
  //  Summarize the run.
  cout << "StMyAnalysisMaker3::Finish()\n";

  // Write event pool manager object to file and close it
  if(mOutNameME != "") {
    TFile *fOutME = new TFile(mOutNameME.Data(), "RECREATE");
    fOutME->cd();
    fListOfPools->Write();  // write pools to file if we have them and want to
    fOutME->Write();
    fOutME->Close();
  }

  //  Write histos to file and close it.
  if(mOutName != "") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(fAnalysisMakerName);
    fout->cd(fAnalysisMakerName);
    WriteHistograms();
   
    // jet shape analysis
    if(doJetShapeAnalysis) {
      for(int j=0; j<9; j++) {
        fout->cd();
        fout->mkdir(Form("JetShapeAnalysis_bin%i", j));
        fout->cd(Form("JetShapeAnalysis_bin%i", j));
        WriteJetShapeHistograms(j);
      }
    }

    fout->cd();
    fout->Write();
    fout->Close();
  }

  //  Write QA histos to file and close it - FIXME, check if statement requirements now that some revamp has been done
  if(mOutNameQA != ""  && (doWriteTrackQAHist || doWriteJetQAHist)) {
    TFile *fQAout = new TFile(mOutNameQA.Data(), "UPDATE");
    fQAout->cd();

    cout<<"print out to see how many time this happens"<<endl;

    // track QA
    if(doWriteTrackQAHist) {
      fQAout->mkdir(Form("TrackQA_bin%i", fTPCptAssocBin));
      fQAout->cd(Form("TrackQA_bin%i", fTPCptAssocBin));
      WriteTrackQAHistograms();
      fQAout->cd();
    }

    // jet QA
    if(doWriteJetQAHist) {
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
  int nHistCentBins = 20;
  if(fCentBinSize == 10) nHistCentBins = 10;
  if(fCentBinSize ==  5) nHistCentBins = 20;

  // binning for mult histograms:       (dopp) ? then A else B
  double kHistMultMax = (doppAnalysis) ? 100. : 800.;
  int kHistMultBins = (doppAnalysis) ? 100 : 800;

  // QA histos
  hdEPReactionPlaneFnc = new TH1F("hdEPReactionPlaneFnc", "jets relative to EP from reaction plane function", 3, 0.0, 0.5*pi);
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
  hCentralityPostCut = new TH1F("hCentralityPostCut", "No. events vs centrality, after cut", nHistCentBins, 0, 100);
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);

  hStats = new TH1F("hStats", "QA stats for cuts and objects", 30, 0.5, 30.5);
  hStats->GetXaxis()->SetBinLabel(1, "total events");
  hStats->GetXaxis()->SetBinLabel(2, "good run event");
  hStats->GetXaxis()->SetBinLabel(3, "passed max trk pt cut");
  hStats->GetXaxis()->SetBinLabel(4, "passed max tow E cut");
  hStats->GetXaxis()->SetBinLabel(5, "passed z-vtx cut");
  hStats->GetXaxis()->SetBinLabel(6, "passed cent cut");
  hStats->GetXaxis()->SetBinLabel(7, "post cent selection");
  hStats->GetXaxis()->SetBinLabel(8, "have jet ptr");
  hStats->GetXaxis()->SetBinLabel(9, "have jets (n>0)");
  hStats->GetXaxis()->SetBinLabel(10, "EmcTriggered HT#"); 
  hStats->GetXaxis()->SetBinLabel(11, "bad tow fired trg");
  hStats->GetXaxis()->SetBinLabel(12, "NO rho object");
  hStats->GetXaxis()->SetBinLabel(13, "past getting rho");
  hStats->GetXaxis()->SetBinLabel(14, "");
  hStats->GetXaxis()->SetBinLabel(15, "do jet analysis");
  hStats->GetXaxis()->SetBinLabel(16, "no pool");
  hStats->GetXaxis()->SetBinLabel(17, "leading jets");
  hStats->GetXaxis()->SetBinLabel(18, "subleading jets");
  hStats->GetXaxis()->SetBinLabel(19, "aj criteria met");
  hStats->GetXaxis()->SetBinLabel(20, "bad cent bin");
  hStats->GetXaxis()->SetBinLabel(21, "1-constit jets");
  hStats->GetXaxis()->SetBinLabel(22, "bad jet(pt) bin");
  hStats->GetXaxis()->SetBinLabel(23, "jets w/o lead bias");
  hStats->GetXaxis()->SetBinLabel(24, "jets w/o firing tow");
  hStats->GetXaxis()->SetBinLabel(25, "no EP(pt) bin");
  hStats->GetXaxis()->SetBinLabel(26, "bad ep bin");
  hStats->GetXaxis()->SetBinLabel(27, "");

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
  hJetE = new TH1F("hJetE", "Jet energy distribution", 100, 0., 100.);
  hJetEta = new TH1F("hJetEta", "Jet #eta distribution", 24, -1.2, 1.2);
  hJetPhi = new TH1F("hJetPhi", "Jet #phi distribution", 72, 0.0, 2.0*pi);
  hJetNEF = new TH1F("hJetNEF", "Jet NEF", 100, 0., 1.);
  hJetArea = new TH1F("hJetArea", "Jet Area", 100, 0., 1.);
  hJetMass = new TH1F("hJetMass", "Jet mass", 65, -5., 60.);
  hJetTracksPt = new TH1F("hJetTracksPt", "Jet track constituent p_{T}", 120, 0, 30.0);
  hJetTracksPhi = new TH1F("hJetTracksPhi", "Jet track constituent #phi", 72, 0, 2.0*pi);
  hJetTracksEta = new TH1F("hJetTracksEta", "Jet track constituent #eta", 56, -1.4, 1.4);
  hJetTracksZ = new TH1F("hJetTracksZ", "Jet track fragmentation function", 144, 0, 1.44);
  hJetPtvsArea = new TH2F("hJetPtvsArea", "Jet p_{T} vs Jet area", 100, 0, 100, 100, 0, 1);

  // jet EP QA histos
  for(int i=0; i<5; i++) {
    hJetEventEP[i] = new TH1F(Form("hJetEventEP%i", i), Form("no of jet events vs event plane, epbin=%i", i), 72, 0.0, 1.0*pi);
    hJetPhivsEP[i] = new TH2F(Form("hJetPhivsEP%i", i), Form("Jet #phi vs event plane, epbin=%i", i), 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);

    hJetPtIn[i] = new TH1F(Form("hJetPtIn%i", i), Form("no of jets in-plane vs p_{T}, epbin=%i", i), 100, 0.0, 100);
    hJetPhiIn[i] = new TH1F(Form("hJetPhiIn%i", i), Form("no of jets in-plane vs #phi, epbin=%i", i), 72, 0.0, 2.0*pi);
    hJetEtaIn[i] = new TH1F(Form("hJetEtaIn%i", i), Form("no of jets in-plane vs #eta, epbin=%i", i), 40, -1.0, 1.0);
    hJetEventEPIn[i] = new TH1F(Form("hJetEventEPIn%i", i), Form("no of in-plane jet events vs event plane, epbin=%i", i), 72, 0.0, 1.0*pi);
    hJetPhivsEPIn[i] = new TH2F(Form("hJetPhivsEPIn%i", i), Form("in-plane Jet #phi vs event plane, epbin=%i", i), 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
    hJetPtMid[i] = new TH1F(Form("hJetPtMid%i", i), Form("no of jets mid-plane vs p_{T}, epbin=%i", i), 100, 0.0, 100);
    hJetPhiMid[i] = new TH1F(Form("hJetPhiMid%i", i), Form("no of jets mid-plane vs #phi, epbin=%i", i), 72, 0.0, 2.0*pi);
    hJetEtaMid[i] = new TH1F(Form("hJetEtaMid%i", i), Form("no of jets mid-plane vs #eta, epbin=%i", i), 40, -1.0, 1.0);
    hJetEventEPMid[i] = new TH1F(Form("hJetEventEPMid%i", i), Form("no of mid-plane jet events vs event plane, epbin=%i", i), 72, 0.0, 1.0*pi);
    hJetPhivsEPMid[i] = new TH2F(Form("hJetPhivsEPMid%i", i), Form("mid-plane Jet #phi vs event plane, epbin=%i", i), 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
    hJetPtOut[i] = new TH1F(Form("hJetPtOut%i", i), Form("no of jets out-of-plane vs p_{T}, epbin=%i", i), 100, 0.0, 100);
    hJetPhiOut[i] = new TH1F(Form("hJetPhiOut%i", i), Form("no of jets out-of-plane vs #phi, epbin=%i", i), 72, 0.0, 2.0*pi);
    hJetEtaOut[i] = new TH1F(Form("hJetEtaOut%i", i), Form("no of jets out-of-plane vs #eta, epbin=%i", i), 40, -1.0, 1.0);
    hJetEventEPOut[i] = new TH1F(Form("hJetEventEPOut%i", i), Form("no of out-of-plane jet events vs event plane, epbin=%i", i), 72, 0.0, 1.0*pi);
    hJetPhivsEPOut[i] = new TH2F(Form("hJetPhivsEPOut%i", i), Form("out-of-plane Jet #phi vs event plane, epbin=%i", i), 72, 0.0, 2.0*pi, 72, 0.0, 1.0*pi);
  }

  hJetHTrigMaxTowEt = new TH1F("hJetHTrigMaxTowEt", "NJet-hadron trigger vs max tower E_{T}", 120, 0., 30.);
  hJetHTrigMaxTrkPt = new TH1F("hJetHTrigMaxTrkPt", "NJet-hadron trigger vs max track p_{T}", 120, 0., 30.);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi", "Jet-hadron #Delta#eta-#Delta#phi", 72, -1.8, 1.8, 72, -0.5*pi, 1.5*pi);

  // Event Selection QA histo + event QA
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);
  hEventTriggerIDs = new TH1F("hEventTriggerIDs", "NTriggers vs trigger IDs", 60, 0.5, 60.5);
  hBadTowerFiredTrigger = new TH1F("hBadTowerFiredTrigger", "# bad tower fired trigger vs tower Id", 4800, 0.5, 4800.5);
  hNGoodTowersFiringTrigger = new TH1F("hNGoodTowersFiringTrigger", "# good towers firing trigger per event", 16, -0.5, 15.5);
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

  // QA for trigged and MB sets
  hNEventsvsZvtxMB5 = new TH1F("hNEventsvsZvtxMB5", "no of event vs zvtx - MB5", 20, -40.0, 40.0);
  hNEventsvsCentMB5 = new TH1F("hNEventsvsCentMB5", "no of event vs centrality - MB5", 20, 0., 100.); 
  hNEventsvsMultMB5 = new TH1F("hNEventsvsMultMB5", "no of event vs multiplicity - MB5", 200, 0., 800.);
  hNEventsvsZvsCentMB5 = new TH2F("hNEventsvsZvsCentMB5", "no of events: zvtx vs Centality - MB5", 20, 0., 100, 20, -40.0, 40.0);
  hNEventsvsZvtxMB30 = new TH1F("hNEventsvsZvtxMB30", "no of events vs zvtx - MB30", 20, -40.0, 40.0);
  hNEventsvsCentMB30 = new TH1F("hNEventsvsCentMB30", "no of event vs centrality - MB30", 20, 0., 100.);
  hNEventsvsMultMB30 = new TH1F("hNEventsvsMultMB30", "no of event vs multiplicity - MB30", 200, 0., 800.);
  hNEventsvsZvsCentMB30 = new TH2F("hNEventsvsZvsCentMB30", "no of events: zvtx vs Centality - MB30", 20, 0., 100, 20, -40.0, 40.0);
  hNEventsvsZvtxMB5Wt = new TH1F("hNEventsvsZvtxMB5Wt", "no of event vs zvtx - MB5, weighted", 20, -40.0, 40.0);
  hNEventsvsCentMB5Wt = new TH1F("hNEventsvsCentMB5Wt", "no of event vs centrality - MB5, weighted", 20, 0., 100.);
  hNEventsvsMultMB5Wt = new TH1F("hNEventsvsMultMB5Wt", "no of event vs multiplicity - MB5, weighted", 200, 0., 800.);
  hNEventsvsZvsCentMB5Wt = new TH2F("hNEventsvsZvsCentMB5Wt", "no of events: zvtx vs Centality - MB5, weighted", 20, 0., 100, 20, -40.0, 40.0);
  hNEventsvsZvtxMB30Wt = new TH1F("hNEventsvsZvtxMB30Wt", "no of events vs zvtx - MB30, weighted", 20, -40.0, 40.0);
  hNEventsvsCentMB30Wt = new TH1F("hNEventsvsCentMB30Wt", "no of event vs centrality - MB30, weighted", 20, 0., 100.);
  hNEventsvsMultMB30Wt = new TH1F("hNEventsvsMultMB30Wt", "no of event vs multiplicity - MB30, weighted", 200, 0., 800.);
  hNEventsvsZvsCentMB30Wt = new TH2F("hNEventsvsZvsCentMB30Wt", "no of events: zvtx vs Centality - MB30, weighted", 20, 0., 100, 20, -40.0, 40.0);
  hNEventsvsZvtxHT2 = new TH1F("hNEventsvsZvtxHT2", "no of events vs zvtx - HT2", 20, -40.0, 40.0);
  hNEventsvsCentHT2 = new TH1F("hNEventsvsCentHT2", "no of event vs centrality - HT2", 20, 0., 100.);
  hNEventsvsMultHT2 = new TH1F("hNEventsvsMultHT2", "no of event vs multiplicity - HT2", 200, 0., 800.);
  hNEventsvsZvsCentHT2 = new TH2F("hNEventsvsZvsCentHT2", "no of events: zvtx vs Centality - HT2", 20, 0., 100, 20, -40.0, 40.0);
  hNPairsvsZvtxMB5 = new TH1F("hNPairsvsZvtxMB5", "no of jet-track pairs vs zvtx - MB5", 20, -40.0, 40.0);
  hNPairsvsZvtxMB30 = new TH1F("hNPairsvsZvtxMB30", "no of jet-track pairs vs zvtx - MB30", 20, -40.0, 40.0);
  hNPairsvsZvtxHT2 = new TH1F("hNPairsvsZvtxHT2", "no of jet-track pairs vs zvtx - HT2", 20, -40.0, 40.0);
  hNPairsvsZvtxMB5Wt = new TH1F("hNPairsvsZvtxMB5Wt", "no of jet-track pairs vs zvtx - MB5, weighted", 20, -40.0, 40.0);
  hNPairsvsZvtxMB30Wt = new TH1F("hNPairsvsZvtxMB30Wt", "no of jet-track pairs vs zvtx - MB30, weighted", 20, -40.0, 40.0);

  for(int k = 0; k<9; k++) {
    hNMixNormBefore[k] = new TH1F(Form("hNMixNormBefore%i", k), Form("number of mixing events - same time filling hnMixEvents: pt bin %i", k), 200, 0, 200);
    hNMixNormAfter[k] = new TH1F(Form("hNMixNormAfter%i", k), Form("number of mixing events after exclusion: pt bin %i", k), 220, -20, 200);
  }

  //// res_cen=new TProfile("res_cen","res vs. cen",10,0,10,-2,2);
  // set binning for run based corrections - run dependent
  //Int_t nRunBins = 1; // - just a default
  //if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) nRunBins = 830;
  //if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) nRunBins = 1378 + 22;
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

  // jet shape analysis background QA histograms
  hBGconeFractionOfJetPt = new TH1F("hBGconeFractionOfJetPt", "Background cone fraction of jet p_{T}", 100, 0.0, 2.0);
  hJetPtvsBGconeFraction = new TH2F("hJetPtvsBGconeFraction", "Jet p_{T} vs Background cone fraction of jet p_{T}", 100, 0.0, 100.0, 100, 0.0, 2.0);
  hJetPtvsBGconePt = new TH2F("hJetPtvsBGconePt", "Jet p_{T} vs Background cone p_{T}", 100, 0.0, 100.0, 100, 0.0, 100.0);

  // QA to check on pt distributions before and after weighting mixed events
  hMB5TrkPtRaw = new TH1F("hMB5TrkPtRaw", "track distribution vs p_{T}, raw", 100, 0., 20.0);
  hMB5TrkPtReWeight = new TH1F("hMB5TrkPtReWeight", "track distribution vs p_{T}, reweight", 100, 0., 20.0);
  hMB30TrkPtRaw = new TH1F("hMB30TrkPtRaw", "track distribution vs p_{T}, raw", 100, 0., 20.0);
  hMB30TrkPtReWeight = new TH1F("hMB30TrkPtReWeight", "track distribution vs p_{T}, reweight", 100, 0., 20.0);

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

  //===============================================
  // Set up mixed event pool settings: binning, etc 
  SetupMixEvtPool();

  // set up jet-hadron sparse
  UInt_t bitcodeMESE = 0; // bit coded, see GetDimParams() below
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8; // | 1<<9 | 1<<10;
  fhnJH = NewTHnSparseF("fhnJH", bitcodeMESE);

  // set up event mixing sparse
  bitcodeMESE = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8; // | 1<<9;
  fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", bitcodeMESE);

  // jet counter for normalizations in correlation analysis
  UInt_t bitcodeCorr = 0; // bit coded, see GetDimparamsCorr() below
  bitcodeCorr = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4;
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
  hEventPlaneFncN2->Write();
  hEventPlaneFncP2->Write();
  hEventPlaneFnc2->Write();
  hEventPlaneClass->Write();

  for(int i=0; i<5; i++) {
    hJetEventEP[i]->Write();
    hJetPhivsEP[i]->Write();

    hJetPtIn[i]->Write();
    hJetPhiIn[i]->Write();
    hJetEtaIn[i]->Write();
    hJetEventEPIn[i]->Write();
    hJetPhivsEPIn[i]->Write();
    hJetPtMid[i]->Write();
    hJetPhiMid[i]->Write();
    hJetEtaMid[i]->Write();
    hJetEventEPMid[i]->Write();
    hJetPhivsEPMid[i]->Write();
    hJetPtOut[i]->Write();
    hJetPhiOut[i]->Write();
    hJetEtaOut[i]->Write();
    hJetEventEPOut[i]->Write();
    hJetPhivsEPOut[i]->Write();
  }
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
  hCentralityPostCut->Write();
  hMultiplicity->Write();
  hStats->Write();
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
  hJetMass->Write();
  hJetTracksPt->Write();
  hJetTracksPhi->Write();
  hJetTracksEta->Write();
  hJetTracksZ->Write();
  hJetPtvsArea->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
  fHistEventSelectionQAafterCuts->Write();
  hEmcTriggers->Write();
  hEventTriggerIDs->Write();
  hBadTowerFiredTrigger->Write();
  hNGoodTowersFiringTrigger->Write();
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

  hNEventsvsZvtxMB5->Write();
  hNEventsvsCentMB5->Write();
  hNEventsvsMultMB5->Write();
  hNEventsvsZvsCentMB5->Write();
  hNEventsvsZvtxMB30->Write();
  hNEventsvsCentMB30->Write();
  hNEventsvsMultMB30->Write();
  hNEventsvsZvsCentMB30->Write();
  hNEventsvsZvtxMB5Wt->Write();
  hNEventsvsCentMB5Wt->Write();
  hNEventsvsMultMB5Wt->Write();
  hNEventsvsZvsCentMB5Wt->Write();
  hNEventsvsZvtxMB30Wt->Write();
  hNEventsvsCentMB30Wt->Write();
  hNEventsvsMultMB30Wt->Write();
  hNEventsvsZvsCentMB30Wt->Write();
  hNEventsvsZvtxHT2->Write();
  hNEventsvsCentHT2->Write();
  hNEventsvsMultHT2->Write();
  hNEventsvsZvsCentHT2->Write();
  hNPairsvsZvtxMB5->Write();
  hNPairsvsZvtxMB30->Write();
  hNPairsvsZvtxHT2->Write();
  hNPairsvsZvtxMB5Wt->Write();
  hNPairsvsZvtxMB30Wt->Write();

  for(int k = 0; k<9; k++) {
    hNMixNormBefore[k]->Write();
    hNMixNormAfter[k]->Write();
  }

  hBGconeFractionOfJetPt->Write();
  hJetPtvsBGconeFraction->Write();
  hJetPtvsBGconePt->Write();

  hMB5TrkPtRaw->Write();
  hMB5TrkPtReWeight->Write();
  hMB30TrkPtRaw->Write();
  hMB30TrkPtReWeight->Write();

  // jet sparses for jet hadron correlations
  if(doJetHadronCorrelationAnalysis) {
    hJetHTrigMaxTowEt->Write();
    hJetHTrigMaxTrkPt->Write();
    fHistJetHEtaPhi->Write();
    //-------------------------------
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

  // event plane method: jet v2 histograms
  if(doJetShapeAnalysis) {
    for(int k=0; k<4; k++) {
      for(int j=0; j<4; j++) {
        for(int i=0; i<4; i++) {
          fProfJetV2[k][j][i]->Write();
        }
      }
    } 
  } // JS
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
  //StMemStat::PrintMem("MyAnalysisMaker at beginning of make");
  hStats->Fill(1);

  // zero out these global variables - may want to initialize these to negative obscure values
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // constants
  const double pi = 1.0*TMath::Pi();

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
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();

  // get run number, check bad runs list if desired (kFALSE if bad)
  fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
    hStats->Fill(2);
  }

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;
  hStats->Fill(3);

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;
  hStats->Fill(4);

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  // cut on (-30, 30) when using NEW centrality definitions - perhaps cut on (-28, 28) when using 4 cm z-vtx bins
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  hEventZVertex->Fill(zVtx);
  hStats->Fill(5);

  // get the fill, and event ID
  int fillId = mPicoEvent->fillId();
  int eventId = mPicoEvent->eventId();
  if(fDebugLevel == kDebugGeneralEvt) cout<<"RunID = "<<fRunNumber<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

  // ============================ CENTRALITY ============================== //
  // see StCentMaker and https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
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
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // calculate corrected multiplicity: 
  // Double_t getRefMultCorr(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag=2) const ;
  // flag=0:  Luminosity only
  // flag=1:  z-vertex only
  // flag=2:  full correction (default)
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 
  hStats->Fill(6);

  // bin-age to use for mixed event and sparses
  Int_t centbin10 = GetCentBin10(centbin);
  double centBinToUse;
  if(fCentBinSize==10)       { centBinToUse = (double)centbin10 * 10.0;
  } else if(fCentBinSize==5) { centBinToUse = (double)centbin * 5.0; }

  // centrality / multiplicity histograms:  event activity - compensate for pp or AuAu
  double kEventActivity = (doppAnalysis) ? (double)grefMult : refCorr2;
  hMultiplicity->Fill(kEventActivity);
  hCentrality->Fill(fCentralityScaled);

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

  // cut on centrality for analysis before doing anything - allow analysis to only run for a specific centrality range
  // see StJetFrameworkPicoBase::fCentralityBinEnum
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }
  hCentralityPostCut->Fill(fCentralityScaled);
  hStats->Fill(7);

  // ========================= Trigger Info =============================== //
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggersHist(hEmcTriggers);

  // fill trigger IDs
  FillTriggerIDs(hEventTriggerIDs);

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds(); 
  if(fDebugLevel == kDebugEmcTrigger) cout<<"Event Trigger-IDs: ";
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
  if( doppAnalysis) fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

  // fill arrays for towers that fired trigger
  FillTowerTriggersArr();

  // switches for Jet and Event Plane analysis
  Bool_t doJetAnalysis = kFALSE; // set false by default
  Bool_t doEPAnalysis = kFALSE;  // set false by default

  // if we have trigger: perform jet analysis
  if(fHaveEmcTrigger) { doJetAnalysis = kTRUE; }

  // if we have trigger && AuAu dataset: run event plane analysis (ADD new runs as needed)
  if(fHaveEmcTrigger && (fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200 || fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB || fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200)) {
    doEPAnalysis = kTRUE;
  }

  // ============================== JetMaker ============================== //
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

  // fill stats histo to do QA
  hStats->Fill(8);
  if(fJets->GetEntries() > 0) hStats->Fill(9);

  // check if bad/dead towers fired trigger and kill event if true
  // 1) check that event is an HT event (specifics: HT1, HT2, HT3 set by user)
  // 2) once we have an HT event, check to see how many triggers/towers fired
  //	- how many were good/bad?
  // 3) did a non-bad tower fire a trigger in this event? 
  if(fHaveEmcTrigger) hStats->Fill(10);
  if(fHaveEmcTrigger && DidBadTowerFireHTTrigger()) { 
    hStats->Fill(11); 
    return kStOK;
  }

  // ======================================================

/*
  // TEST: the below is snippet of code for getting jets and their constituents using constituent subtractor method in StJetMakerTask
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

  // =========================== RhoMaker ============================ //
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
    hStats->Fill(12);
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;
  }

  // get rho/area value from rho object     fRho->ls("");
  //double value = GetRhoValue(fRhoMakerName);
  fRhoVal = fRho->GetVal();
  hRhovsCent->Fill(centbin*5.0, fRhoVal);
  hStats->Fill(13);
  if(fDebugLevel == kDebugRhoEstimate) cout<<"   fRhoVal = "<<fRhoVal<<"   Correction = "<<1.0*TMath::Pi()*fJetRad*fJetRad*fRhoVal<<" GeV/A"<<endl;

  // ==================== Leading and Subleading Jets ===================== //
  // cache the leading + subleading jets within acceptance
  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  // fill leading and subleading pt spectra
  if(fLeadingJet) { 
    hJetLeadingPt->Fill(fLeadingJet->Pt()); 
    hStats->Fill(17);
  }
  if(fSubLeadingJet) { 
    hJetSubLeadingPt->Fill(fSubLeadingJet->Pt());   
    hStats->Fill(18);
  }

  // get dijet imbalance Aj -  z = if(condition) then(?) <do this> else(:) <do this>  
  double fDiJetAj = (fLeadingJet && fSubLeadingJet) ? GetDiJetAj(fLeadingJet, fSubLeadingJet, fRho, fCorrJetPt) : -999.;
  if(fDiJetAj > 0)  hJetDiJetAj->Fill(fDiJetAj);

  // ======================== Event Plane Angle ========================= //
  // basic method - not used for anything..
  double rpAngle = GetReactionPlane();
  hEventPlane->Fill(rpAngle);

  // fill histogram with angle from GetEventPlane() function
  GetEventPlane(kFALSE, 2, fTPCEPmethod, 2.0, 4); // 2nd to last param not used (ptcut) FIXME August 2019 - made other deletions
  hEventPlaneFncN2->Fill(fEPTPCn);
  hEventPlaneFncP2->Fill(fEPTPCp);
  hEventPlaneFnc2->Fill(fEPTPC);
  hEventPlaneClass->Fill(TPC_PSI2);

  // FIXME FIXME - March 3, 2019 - this is outdated at this point, need to recode below chunk if doing ptbin event plane approach
  // 	and wanting to fill histograms additionally here, the event plane maker output
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

  // ========================== Jet Analysis ===================================== //
  if(doRunAnalysis) {
    hStats->Fill(15);

    if(fDebugLevel == kDebugMixedEvents) cout<<"event # = "<<EventCounter()<<"  Centbin = "<<centbin<<"  zVtx = "<<zVtx<<endl;

    // require event mixing
    // 1. First get an event pool corresponding in mult (cent) and zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into FillCorrelations() (diff. now). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    // mix jets from triggered events with tracks from MB events
    // get the trigger bit, need to change trigger bits between different runs

    // declare pool pointer
    StEventPool *pool = 0x0;

    // require event mixing
    if(fDoEventMixing > 0) {
      //===============================================================================================================================
      // convert back to integer bins for mixed event pool - 10% bins (0, 7), 5% bins (0, 15)
      Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSize);
      //cout<<"fCentralityScaled: "<<fCentralityScaled<<"  fCentBinSize: "<<fCentBinSize<<"  mixcentbin: "<<mixcentbin<<"  zVtx: "<<zVtx<<endl;
 
/*
      // mixed event centbin - part of old test
      Int_t mixcentbin;
      if(fCentBinSize==10)       { mixcentbin = centbin10;    // 10% bins (0,  7)
      } else if(fCentBinSize==5) { mixcentbin = centbin; }    //  5% bins (0, 15)
      //cout<<"mixcentbin = "<<mixcentbin<<"  centbin = "<<centbin<<" centbin10 = "<<centbin10<<"  zvtx = "<<zVtx<<endl;
*/
      //===============================================================================================================================

      // get the generic event plane angle of the event: using tracks 0.2-2.0 GeV for calculation (option 4)
      // for an angle (0, pi)    
      const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
      StEventPlaneMaker *EPMaker = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, 4)));
      double psi2 = (EPMaker) ? (double)EPMaker->GetTPCEP() : -999;

      // initialize event pools - different cases for each dataset
      if(fDoUseMultBins) {
        if(!doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx);
        if( doUseEPBins) pool = fPoolMgr->GetEventPool(kEventActivity, zVtx, psi2);
      } else {
        if(!doUseEPBins) pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); // AuAu cent16
        if( doUseEPBins) pool = fPoolMgr->GetEventPool(mixcentbin, zVtx, psi2); 
      }

      // check if pool exists
      if(!pool) {
        Form("No pool found for centrality = %.1f, zVtx = %f", (double)mixcentbin, zVtx);
        hStats->Fill(16);
        return kStOK;
      }

      if(fDebugLevel == kDebugMixedEvents) cout<<"NtracksInPool = "<<pool->NTracksInPool()<<"  CurrentNEvents = "<<pool->GetCurrentNEvents()<<endl;
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

      // get leading jet and subleading jet pt
      if(fCorrJetPt) {
        ljpttemp = fLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
        sljpttemp = fSubLeadingJet->Pt() - fSubLeadingJet->Area()*fRhoVal;
      } else {
        ljpttemp = fLeadingJet->Pt();
        sljpttemp = fSubLeadingJet->Pt(); 
      }

    } // lj + sublj

    // Aj selection for Jet Analysis: LJ pt > 20, SubLJ pt > 10 defaults
    bool doAjSelection = (isBackToBack && ljpttemp > fLeadJetPtMin && sljpttemp > fSubLeadJetPtMin) ? kTRUE : kFALSE;
    if(doAjSelection) {
      hStats->Fill(19);
      hJetLeadingPtAj->Fill(ljpttemp);      // leading jet w/ Aj cut
      hJetSubLeadingPtAj->Fill(sljpttemp);  // sub-leading jet w/ Aj cut
    }

    // ==================================================================================================
    // ==================================================================================================
    // run jet shape analysis
    if(doJetShapeAnalysis) {
      // Triggered events and leading/subleading jets - do Jet Shape Analysis
      // check for back to back jets: must have leading + subleading jet, subleading jet must be > 10 GeV, subleading jet must be within 0.4 of pi opposite of leading jet
      if(doRequireAjSelection) {
        for(int ptbin=0; ptbin<9; ptbin++) {
          if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet)       JetShapeAnalysis(fLeadingJet, pool, kEventActivity, ptbin);
          if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets && fSubLeadingJet) JetShapeAnalysis(fSubLeadingJet, pool, kEventActivity, ptbin);
        }
      } else { // don't require back-to-back jets meeting Aj criteria
        for(int ptbin=0; ptbin<9; ptbin++) {
          if(fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet)       JetShapeAnalysis(fLeadingJet, pool, kEventActivity, ptbin);
          if(fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets && fSubLeadingJet) JetShapeAnalysis(fSubLeadingJet, pool, kEventActivity, ptbin);
        }
      }

      // jet shape - case for: inclusive jets
      if(fHaveEmcTrigger && fJetAnalysisJetType == kInclusiveJets) {
        // loop over Jets in the event
        for(int ijet = 0; ijet < fJets->GetEntries(); ijet++) {  // JET LOOP
          // get pointer to jet
          StJet *jet = static_cast<StJet*>(fJets->At(ijet));
          if(!jet) continue;

          // loop over pt associated bins
          for(int ptbin = 0; ptbin < 9; ptbin++) {
            if(doRequireAjSelection && doAjSelection) { JetShapeAnalysis(jet, pool, kEventActivity, ptbin);
            } else { JetShapeAnalysis(jet, pool, kEventActivity, ptbin); }
          } // loop over pt bins
        }   // loop over jets
      }     // inclusive jet case

    }       // JET SHAPE analysis

    // =====================================================================================================
    // =====================================================================================================
    // run jet hadron correlation analysis
    if(doJetHadronCorrelationAnalysis) {
      // Triggered events and leading/subleading jets - do Jet Hadron Correlation Analysis
      // check for back to back jets: must have leading + subleading jet, subleading jet must be > 10 GeV, subleading jet must be within 0.4 of pi opposite of leading jet
      if(doRequireAjSelection) {
        for(int ptbin=0; ptbin<5; ptbin++) {
          if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet)       JetHadronCorrelationAnalysis(fLeadingJet, pool, centbin, ptbin);
          if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets && fSubLeadingJet) JetHadronCorrelationAnalysis(fSubLeadingJet, pool, centbin, ptbin);
        }
      } else { // don't require back-to-back jets meeting Aj criteria
        for(int ptbin=0; ptbin<5; ptbin++) {
          if(fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet)       JetHadronCorrelationAnalysis(fLeadingJet, pool, centbin, ptbin);
          if(fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets && fSubLeadingJet) JetHadronCorrelationAnalysis(fSubLeadingJet, pool, centbin, ptbin);
        }
      }

      // jet-hadron correlation - case for: inclusive jets
      if(fHaveEmcTrigger && fJetAnalysisJetType == kInclusiveJets) {
        // loop over jets
        for(int ijet = 0; ijet < fJets->GetEntries(); ijet++) {  // JET LOOP
          // get pointer to jet
          StJet *jet = static_cast<StJet*>(fJets->At(ijet));
          if(!jet) continue;

          // loop over pt associated bins for analysis
          for(int ptbin = 0; ptbin < 5; ptbin++) {
            if(doRequireAjSelection && doAjSelection) { JetHadronCorrelationAnalysis(jet, pool, centbin, ptbin);
            } else { JetHadronCorrelationAnalysis(jet, pool, centbin, ptbin); }
          } // loop over pt bins
        }   // loop over jets
      }     // inclusive jet case

    }       // JET-HADRON CORRELATION analysis switch
    // ==================================================================================================

    // use only tracks from MB events
    //if(fDoEventMixing > 0 && fRunForMB && (!fHaveEmcTrigger)) { // kMB5 or kMB30 - AuAu, kMB - pp (excluding HT)
    if(fDoEventMixing > 0 && fRunForMB) { // kMB5 or kMB30 - AuAu, kMB - pp (don't exclude HT)
      ///==///    if(fDoEventMixing > 0 && fHaveMB5event && !fHaveMB30event) { // kMB5 and NOT kMB30
      ///==///    if(fDoEventMixing > 0 && !fHaveMB5event && fHaveMB30event) { // kMB30 and NOT kMB5

      // kill mixing when both MB5 and MB30 trigger for the event - FIXME - hardcoded cutoff
      // update pool: create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      if((fHaveMB5event && !fHaveMB30event && zVtx <= 16.0) || 
         (fHaveMB30event && !fHaveMB5event)) {
        ///==///      if(!fHaveMB5event && fHaveMB30event) {
        ///==///      if(fHaveMB5event && !fHaveMB30event) {
        pool->UpdatePool(CloneAndReduceTrackList());
      }

      // fill QA histo's
      hMBvsMult->Fill(kEventActivity);    // MB5 || MB30 AuAu and MB pp
      if(fHaveMB5event)  hMB5vsMult->Fill(refCorr2);   // MB5
      if(fHaveMB30event) hMB30vsMult->Fill(refCorr2);  // MB30
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB + event mixing 

  } // end of jet analysis

  // QA plots - December2020
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200 || fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) {
    // do stuff here...
    //	1. zvertex (cut at +/- 28cm)
    //	2. multiplicity (limited from cut)
    //	3. centrality (limited from cut)
    //  4. trigger counts - from just the fills

    //----------------------==---------------------==
    // correction factors generated for mixed event triggered events
    // bin 4-17 (-28 cm < z < +28 cm)
    // === Event weight ===
    double mCorrEvtMB5[] = {0.,0.,0., 187.641, 186.714, 54.9565, 9.34451, 2.15293, 0.509046, 0.279472, 0.375596, 1.2418, 8.87955, 55.0963, 219.312, 222.645, 100.253, 0.,0.,0. };
    double mCorrEvtMB30[]= {0.,0.,0., 0.965911, 0.944671, 0.961855, 0.969395, 0.992681, 1.00764, 1.01291, 1.01672, 1.02394, 1.01706, 1.01427, 1.02102, 0.995457, 0.933781, 0.,0.,0., };
    // === Pair weight ===
    //double mCorrPairMB5[] = {0.,0.,0.,   136.554, 57.9349, 11.4968, 2.20195, 0.873841, 0.502317, 0.477275, 0.486535, 0.604795, 2.26412, 11.5256, 62.9154, 64.5642, 316.419,   0.,0.,0.};
    //double mCorrPairMB30[]= {0.,0.,0.,   1.73895, 1.3147, 1.11781, 1.00823, 0.967721, 0.906712, 0.900274, 0.894385, 0.900898, 0.921108, 1.01382, 1.13631, 1.44455, 2.22927,   0.,0.,0.};
    double mMBTrigZWeight = 1.;
    int zbin = GetZVertex4cmBin(zVtx);
    // FIXME //if(assocPtBin == 0 && ibg < 2) cout<<"mMBTrig: "<<mMBTrig<<"   mMBTrigZWeight: "<<mMBTrigZWeight<<endl;  // FIXME  - here for test
    //----------------------==---------------------==

    if(fHaveMB5event)   {
      hNEventsvsZvtxMB5->Fill(zVtx);
      hNEventsvsCentMB5->Fill(fCentralityScaled);
      hNEventsvsMultMB5->Fill(kEventActivity);
      hNEventsvsZvsCentMB5->Fill(fCentralityScaled, zVtx);
      //---------------------------------------------------- Weighted distributions
      mMBTrigZWeight = mCorrEvtMB5[zbin];
      hNEventsvsZvtxMB5Wt->Fill(zVtx, mMBTrigZWeight);
      hNEventsvsCentMB5Wt->Fill(fCentralityScaled, mMBTrigZWeight);
      hNEventsvsMultMB5Wt->Fill(kEventActivity, mMBTrigZWeight);
      hNEventsvsZvsCentMB5Wt->Fill(fCentralityScaled, zVtx, mMBTrigZWeight);
    }
    if(fHaveMB30event)  { 
      hNEventsvsZvtxMB30->Fill(zVtx);
      hNEventsvsCentMB30->Fill(fCentralityScaled);
      hNEventsvsMultMB30->Fill(kEventActivity);
      hNEventsvsZvsCentMB30->Fill(fCentralityScaled, zVtx);
      //----------------------------------------------------- Weighted distributions
      mMBTrigZWeight = mCorrEvtMB30[zbin];
      hNEventsvsZvtxMB30Wt->Fill(zVtx, mMBTrigZWeight);
      hNEventsvsCentMB30Wt->Fill(fCentralityScaled, mMBTrigZWeight);
      hNEventsvsMultMB30Wt->Fill(kEventActivity, mMBTrigZWeight);
      hNEventsvsZvsCentMB30Wt->Fill(fCentralityScaled, zVtx, mMBTrigZWeight);
    }    
    if(fHaveEmcTrigger) {  
      hNEventsvsZvtxHT2->Fill(zVtx);
      hNEventsvsCentHT2->Fill(fCentralityScaled);
      hNEventsvsMultHT2->Fill(kEventActivity);
      hNEventsvsZvsCentHT2->Fill(fCentralityScaled, zVtx);
    }
  }  // do QA plots for Run14

  // run Track QA and fill histograms
  ///if((doWriteTrackQAHist) && (doJetAnalysis)) TrackQA();
  ///RunJetQA();

  // event counter at end of maker
  mInputEventCounter++;

  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQAafterCuts);
  //StMemStat::PrintMem("MyAnalysisMaker at end of make");

  //===========================================================================
  if(fRunForMB) {
    // loop over tracks
    int nTracks = mPicoDst->numberOfTracks();
    for(int i = 0; i < nTracks; i++) {
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

      // track variables - used with alt method below
      double pt = mTrkMom.Perp();
      double eta = mTrkMom.PseudoRapidity();
      double phi = mTrkMom.Phi();      

      // shift angle (0, 2*pi) 
      if(phi < 0.0)    phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;

      // calculate single particle tracking efficiency
      int effCent   = mCentMaker->GetRef16();
      double fZDCx  = mPicoEvent->ZDCx();
      double trkEff = ApplyTrackingEff(fDoEffCorr, pt, eta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile);

      //  - MB30: weight is just 're-weight'
      //  - MB5: weight is 're-weight' scaled by additional weight for MB5 -> MB30 
      double refMultReWeightCorr = 1.0;
      if( fHaveMB5event && !fHaveMB30event) {
        refMultReWeightCorr = mCentMaker->GetMB5toMB30ReWeight();
        hMB5TrkPtRaw->Fill(pt, 1.0/trkEff);
        hMB5TrkPtReWeight->Fill(pt, refMultReWeightCorr/trkEff);
      }
      if(fHaveMB30event &&  !fHaveMB5event) {
        refMultReWeightCorr = mCentMaker->GetReWeight();
        hMB30TrkPtRaw->Fill(pt, 1.0/trkEff);
        hMB30TrkPtReWeight->Fill(pt, refMultReWeightCorr/trkEff);
      }

    } // end of trk loop
  }  // MB event

  return kStOK;
}
//
//______________________________________________________________________
THnSparse *StMyAnalysisMaker3::NewTHnSparseF(const char *name, UInt_t entries)
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
         hnTitle += Form(";%s", label.Data());
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
      label = "centrality";
      if(fCentBinSize==10) {
        nbins = 10;
      } else {
        nbins = 20;
      }
      xmin = 0.;
      xmax = 100.;
      break;

   case 1:
      if(fCorrJetPt) { // correct jet pt
        label = "Jet corrected #it{p}_{T}";
        nbins = 30;
        xmin = -50.;
        xmax = 100.;
      } else { // don't correct jet pt
        label = "Jet #it{p}_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      }
      break;

   case 2:
      label = "Track #it{p}_{T}";
      nbins = 80; 
      xmin =  0.;
      xmax = 20.;
      break;

   case 3:
      label = "Relative eta";
      nbins = 72; // 48
      xmin = -1.8;
      xmax =  1.8;
      break;

   case 4: 
      label = "Relative phi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax =  1.5*pi;
      break;

   case 5:
      label = "Relative angle of jet and event plane";
      nbins = 3; // (12) 72
      xmin = 0.;
      xmax = 0.5*pi;
      break;

   case 6:
      label = "z-vertex";
      nbins = 20; // 10
      xmin = -40.; //-10
      xmax =  40.; //+10
      break;

   case 7:
      label = "track charge";
      nbins = 3;
      xmin = -1.5;
      xmax =  1.5;
      break;

   case 8:
      label = "TPC associated #it{p}_{T} EPbin";
      nbins = 5;
      xmin = -0.5;
      xmax =  4.5;
      break;

   case 9:
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax =  2.5;
      break;

   } // end of switch
} // end of getting dim-params
//
//______________________________________________________________________
THnSparse *StMyAnalysisMaker3::NewTHnSparseFCorr(const char *name, UInt_t entries) {
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
      hnTitle += Form(";%s", label.Data());
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
      label = "centrality";
      if(fCentBinSize==10) {
        nbins = 10;
      } else {
        nbins = 20;
      }
      xmin = 0.;
      xmax = 100.;
      break;

    case 1:
      if(fCorrJetPt) { // correct jet pt
        label = "Jet corrected #it{p}_{T}";
        nbins = 30;
        xmin = -50.;
        xmax = 100.;
      } else { // don't correct jet pt
        label = "Jet #it{p}_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      }
      break;

    case 2:
      label = "Relative angle: jet and event plane";
      nbins = 3; // (12)
      xmin = 0.;
      xmax = 0.5*pi;
      break;

    case 3:
      label = "z-vertex";
      nbins = 20;
      xmin = -40.;
      xmax =  40.;
      break;

    case 4:
      label = "TPC associated #it{p}_{T} EPbin";
      nbins = 5;
      xmin = -0.5;
      xmax =  4.5;  
      break;

   }// end of switch
} // end of Correction (ME) sparse
//
// From CF event mixing code PhiCorrelations
// clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
//____________________________________________________________________________
TClonesArray *StMyAnalysisMaker3::CloneAndReduceTrackList()
{
  // create array for Femto tracks
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
      mTrkMom = trk->pMom();                 // get primary track vector
    } else {
      mTrkMom = trk->gMom(mVertex, Bfield);  // get global track vector
    }

    // track variables - used with alt method below
    double pt = mTrkMom.Perp();

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    // this is TEMP, it will filter track by the pt bin used for analysis - double check this syntax
    if(doTPCptassocBin && fDoFilterPtMixEvents) {
      if(fTPCptAssocBin == 0) { if((pt < 0.20) || (pt >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 1) { if((pt < 0.50) || (pt >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 2) { if((pt < 1.00) || (pt >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 3) { if((pt < 1.50) || (pt >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 4) { if((pt < 2.00) || (pt >= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(fTPCptAssocBin == 5) { if((pt < 2.00) || (pt >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 6) { if((pt < 3.00) || (pt >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(fTPCptAssocBin == 7) { if((pt < 4.00) || (pt >= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    //==============================================================================================
    // get trigger to separate correction weight for the min bias events (kVPDMB5 and kVPDMB30)
    bool fHaveMB5event  = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
    bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
    int fMixMBTrig = 0;
    if( fHaveMB5event && !fHaveMB30event) fMixMBTrig = 5;
    if(fHaveMB30event &&  !fHaveMB5event) fMixMBTrig = 30;

    //  - MB30: weight is just 're-weight'
    //  - MB5: weight is 're-weight' scaled by additional weight for MB5 -> MB30 
    double refMultReWeightCorr = 1.0;
    if( fHaveMB5event && !fHaveMB30event) refMultReWeightCorr = mCentMaker->GetMB5toMB30ReWeight();
    if(fHaveMB30event && !fHaveMB5event)  refMultReWeightCorr = mCentMaker->GetReWeight();
    //==============================================================================================

    // create StFemtoTracks out of accepted tracks - light-weight object for mixing
    //  StFemtoTrack *t = new StFemtoTrack(pt, eta, phi, charge);
    //StFemtoTrack *t = new StFemtoTrack(trk, Bfield, mVertex, doUsePrimTracks);
    StFemtoTrack *t = new StFemtoTrack(trk, Bfield, mVertex, doUsePrimTracks, refMultReWeightCorr, fMixMBTrig);
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
TH1 *StMyAnalysisMaker3::FillEmcTriggersHist(TH1 *h) {
  // set bin labels
  h->GetXaxis()->SetBinLabel(1, "HT0");
  h->GetXaxis()->SetBinLabel(2, "HT1");
  h->GetXaxis()->SetBinLabel(3, "HT2");
  h->GetXaxis()->SetBinLabel(4, "HT3");
  h->GetXaxis()->SetBinLabel(5, "JP0");
  h->GetXaxis()->SetBinLabel(6, "JP1");
  h->GetXaxis()->SetBinLabel(7, "JP2");
  h->GetXaxis()->SetBinLabel(10, "Any");
  h->LabelsOption("v");  // set x-axis labels vertically
  //h->LabelsDeflate("X");

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

    // print some EMCal Trigger info: exclude JP2, JP1
    if(emcTrig->flag() != 112 && emcTrig->flag() != 96) { // FIXME - here for test May 30, 2019
      if(fDebugLevel == kDebugEmcTrigger) {
        cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
        cout<<"  isHT0: "<<isHT0<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3;
        cout<<"  isJP0: "<<isJP0<<"  isJP1: "<<isJP1<<"  isJP2: "<<isJP2<<endl;
      }
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

  return h;
}
//
// general function to get local reaction plane calculation
//________________________________________________________________________
Double_t StMyAnalysisMaker3::GetReactionPlane() { 
  // initialize some variables and constants
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
    if(doUsePrimTracks) { // get primary track vector
      mTrkMom = track->pMom();
    } else {              // get global track vector
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

    // check for leading jet removal - taken from Redmers approach (CHECK!)
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
  hEventPlaneFncN2->Sumw2();
  hEventPlaneFncP2->Sumw2();
  hEventPlaneFnc2->Sumw2();
  hEventPlaneClass->Sumw2();

  hEventPlane->Sumw2();
  fHistEPTPCn->Sumw2();
  fHistEPTPCp->Sumw2();
  fHistEPBBC->Sumw2();
  fHistEPZDC->Sumw2();
  hEventZVertex->Sumw2();
  hCentrality->Sumw2();
  hCentralityPostCut->Sumw2();
  hMultiplicity->Sumw2();
  //hStats->Sumw2();
  hRhovsCent->Sumw2();
  for(int i=0; i<5; i++) { hdEPtrk[i]->Sumw2(); }
  for(int i=0; i<9; i++) { // centrality
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
  hJetMass->Sumw2();
  hJetTracksPt->Sumw2();
  hJetTracksPhi->Sumw2();
  hJetTracksEta->Sumw2();
  hJetTracksZ->Sumw2();
  hJetPtvsArea->Sumw2();

  for(int i=0; i<5; i++) {
    hJetEventEP[i]->Sumw2();
    hJetPhivsEP[i]->Sumw2();

    hJetPtIn[i]->Sumw2();
    hJetPhiIn[i]->Sumw2();
    hJetEtaIn[i]->Sumw2();
    hJetEventEPIn[i]->Sumw2();
    hJetPhivsEPIn[i]->Sumw2();
    hJetPtMid[i]->Sumw2();
    hJetPhiMid[i]->Sumw2();
    hJetEtaMid[i]->Sumw2();
    hJetEventEPMid[i]->Sumw2();
    hJetPhivsEPMid[i]->Sumw2();
    hJetPtOut[i]->Sumw2();
    hJetPhiOut[i]->Sumw2();
    hJetEtaOut[i]->Sumw2();
    hJetEventEPOut[i]->Sumw2();
    hJetPhivsEPOut[i]->Sumw2();
  }

  hJetHTrigMaxTowEt->Sumw2();
  hJetHTrigMaxTrkPt->Sumw2();
  fHistJetHEtaPhi->Sumw2();
  fHistEventSelectionQA->Sumw2();
  fHistEventSelectionQAafterCuts->Sumw2();
  hEmcTriggers->Sumw2();
  hEventTriggerIDs->Sumw2();
  hBadTowerFiredTrigger->Sumw2();
  hNGoodTowersFiringTrigger->Sumw2();
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

  hNEventsvsZvtxMB5->Sumw2();
  hNEventsvsCentMB5->Sumw2();
  hNEventsvsMultMB5->Sumw2();
  hNEventsvsZvsCentMB5->Sumw2();
  hNEventsvsZvtxMB30->Sumw2();
  hNEventsvsCentMB30->Sumw2();
  hNEventsvsMultMB30->Sumw2();
  hNEventsvsZvsCentMB30->Sumw2();
  hNEventsvsZvtxMB5Wt->Sumw2();
  hNEventsvsCentMB5Wt->Sumw2();
  hNEventsvsMultMB5Wt->Sumw2();
  hNEventsvsZvsCentMB5Wt->Sumw2();
  hNEventsvsZvtxMB30Wt->Sumw2();
  hNEventsvsCentMB30Wt->Sumw2();
  hNEventsvsMultMB30Wt->Sumw2();
  hNEventsvsZvsCentMB30Wt->Sumw2();
  hNEventsvsZvtxHT2->Sumw2();
  hNEventsvsCentHT2->Sumw2();
  hNEventsvsMultHT2->Sumw2();
  hNEventsvsZvsCentHT2->Sumw2();
  hNPairsvsZvtxMB5->Sumw2();
  hNPairsvsZvtxMB30->Sumw2();
  hNPairsvsZvtxHT2->Sumw2();
  hNPairsvsZvtxMB5Wt->Sumw2();
  hNPairsvsZvtxMB30Wt->Sumw2();

  for(int k = 0; k<9; k++) {
    hNMixNormBefore[k]->Sumw2();
    hNMixNormAfter[k]->Sumw2();
  }

  hBGconeFractionOfJetPt->Sumw2();
  hJetPtvsBGconeFraction->Sumw2(); 
  hJetPtvsBGconePt->Sumw2();

  hMB5TrkPtRaw->Sumw2();
  hMB5TrkPtReWeight->Sumw2();
  hMB30TrkPtRaw->Sumw2();
  hMB30TrkPtReWeight->Sumw2();

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
      excludeInPhi = fLeadingJet->Phi();
      if(fDebugLevel == kDebugLeadSubLeadJets) cout<<"leading: pt = "<<fLeadingJet->Pt()<<"  eta = "<<fLeadingJet->Eta()<<"  phi = "<<fLeadingJet->Phi()<<endl;
    }

    // check for subleading jet
    if(fSubLeadingJet) {
      excludeInEtaSub = fSubLeadingJet->Eta();
      excludeInPhiSub = fSubLeadingJet->Phi();
      if(fDebugLevel == kDebugLeadSubLeadJets) cout<<"subleading: pt = "<<fSubLeadingJet->Pt()<<"  eta = "<<fSubLeadingJet->Eta()<<"  phi = "<<fSubLeadingJet->Phi()<<endl;
    }
  } // leading jets

  // loop over tracks
  TRandom3 *rand = new TRandom3();
  int nTOT = 0, nA = 0, nB = 0;
  int nTrack = mPicoDst->numberOfTracks();
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) { // get primary track vector
      mTrkMom = track->pMom();
    } else {              // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double eta = mTrkMom.PseudoRapidity();
    double phi = mTrkMom.Phi();

    // shift phi to be [0, 2pi]
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    if(doTPCptassocBin) {
      if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(ptbin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(ptbin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(ptbin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    // Method1: kRemoveEtaStrip - remove strip only when we have a leading jet
    if(fTPCEPmethod == kRemoveEtaStrip){
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && 
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;

    } else if(fTPCEPmethod == kRemoveEtaPhiCone){
      // Method2: kRemoveEtaPhiCone - remove cone (in eta and phi) around leading jet
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && (deltaR < fJetRad )) continue;

    } else if(fTPCEPmethod == kRemoveLeadingJetConstituents){
      // Method3: kRemoveLeadingJetConstituents - remove tracks above 2 GeV in cone around leading jet
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;

    } else if(fTPCEPmethod == kRemoveEtaStripLeadSub){
      // Method4: kRemoveEtaStripLeadSub - remove strip only when we have a leading + subleading jet
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && ((TMath::Abs(eta - excludeInEtaSub) < fJetRad*fExcludeLeadingJetsFromFit ) )) continue;

    } else if(fTPCEPmethod == kRemoveEtaPhiConeLeadSub){
      // Method5: kRemoveEtaPhiConeLeadSub - remove cone (in eta and phi) around leading + subleading jet
      double deltaR    = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi - excludeInPhiSub)*(phi - excludeInPhiSub));
      if((fLeadingJet)    && (fExcludeLeadingJetsFromFit > 0) && (deltaR    < fJetRad )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && (deltaRSub < fJetRad )) continue;

    } else if(fTPCEPmethod == kRemoveLeadingSubJetConstituents){
      // Method6: kRemoveLeadingSubJetConstituents - remove tracks above 2 GeV in cone around leading + subleading jet
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi - excludeInPhiSub)*(phi - excludeInPhiSub));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && (pt > fJetConstituentCut) && (deltaRSub < fJetRad)) continue;

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
    if(doUsePrimTracks) { // get primary track vector
      mTrkMom = track->pMom();
    } else {              // get global track vector
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
    // FIXME - the angle TPC_PSI2 is no longer meaningful - August8, 2019
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
  for(int i = 0; i < 4800; i++) {
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
    int emcTrigIDindex = emcTrigID - 1;

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();

    // set trigger level type
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigIDindex] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigIDindex] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigIDindex] = kTRUE;

    //cout<<"i = "<<i<<"  EmcTrigID = "<<emcTrigID<<"  adc = "<<emcTrig->adc()<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3<<endl;
  }

}
//
// function to require that a jet constituent tower fired a HT trigger
//___________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::DidTowerConstituentFireTrigger(StJet *jet) {  
  // tower constituent fired trigger
  Bool_t mFiredTrigger = kFALSE;

  // loop over constituent towers
  for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
    int towerIndex = jet->ClusterAt(itow);

    // get tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towerIndex));
    if(!tow) continue;
    
    // tower ID: get from index of array shifted by +1
    int towID = towerIndex + 1;
    int towIDindex = towID - 1;
    if(towID < 0) continue;

    // change flag to true if jet tower fired trigger
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towIDindex]) mFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towIDindex]) mFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towIDindex]) mFiredTrigger = kTRUE;
  
  } // tower constituent loop

  return mFiredTrigger;
}  
//
// FIXME TODO - still working on how this function will perform.. May 31, 2019
// function to check if a bad tower fired the events HT trigger
//	- this function is called if the event is flagged by HT trigger threshold (set by USER)
//___________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::DidBadTowerFireTrigger() {
  // bad/dead tower fired trigger
  bool mBadTowerFiredTrigger = kFALSE;
  bool mBadTowerFiring = kFALSE;
  int nGood = 0;

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    // update status
    mBadTowerFiredTrigger = kFALSE;

    // get tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tow) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // tower ID: get from index of array shifted by +1
    int towID = itow + 1;
    int towIDindex = towID - 1;
    if(towID < 0) continue; // double check these aren't still in the event list

    // check if tower is bad or dead - functions return kTRUE if ok and kFALSE if NOT dead 
    // isTowerOK(towID) = kTRUE if tower is OK and not dead
    // isTowerDead(towID) = kTRUE if tower is dead, = kFALSE if tower is NOT dead
    //bool isTowOk = (mBaseMaker->IsTowerOK(towID) && !mBaseMaker->IsTowerDead(towID)); // FIXME
    bool isTowOk = (mBaseMaker->IsTowerOK(towID));

    // change flag to true if bad tower fired trigger
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towIDindex] && !isTowOk) mBadTowerFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towIDindex] && !isTowOk) mBadTowerFiredTrigger = kTRUE;
    if((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towIDindex] && !isTowOk) mBadTowerFiredTrigger = kTRUE;

    // check if good tower fired trigger
    if( ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && fTowerToTriggerTypeHT1[towIDindex] && isTowOk ) ||
        ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && fTowerToTriggerTypeHT2[towIDindex] && isTowOk ) ||
        ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && fTowerToTriggerTypeHT3[towIDindex] && isTowOk )) {
      nGood++;
      //cout<<"TriggerType: "<<fEmcTriggerEventType<<"    Legit tower firing, towID: "<<towID<<endl;
    }

    // for bad tower firings: fill histo, inform user, return kTRUE
    if(mBadTowerFiredTrigger == kTRUE) {
      hBadTowerFiredTrigger->Fill(towID);
      mBadTowerFiring = kTRUE;
      cout<<"Bad tower fired trigger, towID: "<<towID<<endl;
      //break;
    }

  } // tower loop

  // total good firings:
  //cout<<"nGood: "<<nGood<<endl;
  //hNGoodTowersFiringTrigger->Fill(nGood);

  // ===========================================
  // TODO - TODO
  // problem is: a few towers always fire when BAD
  // - 1) reject event if any bad tower fires a HT of interest?
  // - 2) keep event if there are good tower(s) that fire an event that a bad tower also fired in?
  // 	> this seems the way to go, else we would be eliminating nearly all events
  //
  // - Note: tower 3405!
  // ===========================================

  // kTRUE: bad tower fired trigger
  //return mBadTowerFiredTrigger;
  return mBadTowerFiring;
}
//
// FIXME TODO - still working on how this function will perform.. May 31, 2019
// function to check if a bad tower fired the events HT trigger
//      - this function is called if the event is flagged by HT trigger threshold (set by USER)
//___________________________________________________________________________________________
Bool_t StMyAnalysisMaker3::DidBadTowerFireHTTrigger() {
  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers2 = mPicoEvent->triggerIds();
  if(fDebugLevel == kDebugTowersFiringTriggers) cout<<endl<<"Event Trigger-IDs: ";
  for(unsigned int i=0; i<mytriggers2.size(); i++) {
    if(fDebugLevel == kDebugTowersFiringTriggers) cout<<"i = "<<i<<": "<<mytriggers2[i] << ", ";
  }
  if(fDebugLevel == kDebugTowersFiringTriggers) cout<<endl;

  // bad/dead tower fired trigger
  bool mBadTowerFiredTrigger = kTRUE; //kFALSE;
  int nGood = 0;

  // get number of Emc triggers
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get pointer to trigger
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // emc trigger parameters
    int emcTrigID = emcTrig->id();
    int towerID = emcTrigID;        // tower ID (1-4800) - tower ID is directly related to referenced emcTrigID
    if(towerID < 0) { cout<<"tower ID < 0, tower ID = "<<towerID<<endl; continue; } // double check these aren't still in the event list

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT0 = emcTrig->isHT0();
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    bool isJP0 = emcTrig->isJP0();
    bool isJP1 = emcTrig->isJP1();
    bool isJP2 = emcTrig->isJP2();

    // get associated tower pointer - minus 1 to get array element
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towerID - 1));
    if(!tower) { cout<<"No tower pointer... iTow = "<<i<<"  towerID = "<<towerID<<endl; continue; }

    // check if tower is bad or dead - functions return kTRUE if ok and kFALSE if NOT dead 
    // isTowerOK(towerID) = kTRUE if tower is OK and not dead
    // isTowerDead(towerID) = kTRUE if tower is dead, = kFALSE if tower is NOT dead
    bool isTowOk = (mBaseMaker->IsTowerOK(towerID));

    // eliminate jet-patch JP triggers
    if(emcTrig->flag() == 112 || emcTrig->flag() == 96) continue;

    // check if good tower fired trigger
    if( ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && isHT1 && isTowOk ) ||
        ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && isHT2 && isTowOk ) ||
        ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && isHT3 && isTowOk )) {
      nGood++;
    } 

    // bad towers firing the trigger
    if( ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT1) && isHT1 && !isTowOk ) ||
        ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT2) && isHT2 && !isTowOk ) ||
        ((fEmcTriggerEventType == StJetFrameworkPicoBase::kIsHT3) && isHT3 && !isTowOk )) {
      hBadTowerFiredTrigger->Fill(towerID);
    }

    // print some EMCal Trigger info: exclude JP2, JP1
    if((fDebugLevel == kDebugTowersFiringTriggers) && (emcTrig->flag() != 112) && (emcTrig->flag() != 96)) {
      cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
      cout<<"  isHT0: "<<isHT0<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3;
      cout<<"  isJP0: "<<isJP0<<"  isJP1: "<<isJP1<<"  isJP2: "<<isJP2<<endl;
    }

  }

  // determine if we should keep event
  if(nGood > 0) mBadTowerFiredTrigger = kFALSE;
  //cout<<"nGood: "<<nGood<<"  get rid of event? - "<<mBadTowerFiredTrigger<<endl;

  // total triggers which had good towers that fired them in this event:
  hNGoodTowersFiringTrigger->Fill(nGood);
  
  return mBadTowerFiredTrigger;
}
// ===========================================================================================
// ===========================================================================================
//
// function thats runs jet shape analysis
// 	refCorr2 is refCorr2 for AuAu, but is grefMult for pp
//___________________________________________________________________________________________
void StMyAnalysisMaker3::JetShapeAnalysis(StJet *jet, StEventPool *pool, Double_t refCorr2, Int_t assocPtBin) {
    // constants
    double pi = 1.0*TMath::Pi();
    double rbinSize = 0.05;

    // get centrality bin
    int centBin = (!doppAnalysis) ? Get4CentBin(fCentralityScaled) : 0; // this secures it will still run for pp, and use cent bin = 0
    if(centBin < 0) {
      if(assocPtBin==0) hStats->Fill(20);
      return;
    }

    // check for jets with only 1 constituent
    int nJetConstituents = jet->GetNumberOfConstituents();
    if(nJetConstituents == 1) {
      if(assocPtBin==0) hStats->Fill(21);
      if(doSkip1ParticleJets)  return;
    }

    // get jet info
    double jetPhi = jet->Phi();
    double jetEta = jet->Eta();
    double jetArea = jet->Area();
    double jetMass = jet->M();
    double uncorrjetPt = jet->Pt();
    double corrjetPt = jet->Pt() - jetArea*fRhoVal;
    double jetE = jet->E();
    double jetNEF = jet->NEF();

    // get jet pt bin and value
    double jetPt = -99.;
    if(fCorrJetPt) { jetPt = jet->Pt() - jet->Area()*fRhoVal; // remove underlying event background is desired/needed
    } else { jetPt = jet->Pt(); }
    int jetPtBin = GetJetPtBin(jetPt);
    if(jetPtBin < 0) { // cut on jets outside analysis pt range
      if(assocPtBin==0) hStats->Fill(22);
      return;
    }

    // require tower and or track bias for jet
    // apply leading bias to jet
    if(doBiasJetLeadConstituent && (jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) {
      if(assocPtBin==0) hStats->Fill(23);
      return;
    }

    // check that jet contains a tower that fired the trigger
    if(doRequireJetTowFireTrig && !DidTowerConstituentFireTrigger(jet)) {
      if(assocPtBin==0) hStats->Fill(24);
      return;
    }

    // ==========================================================================================================
    // get StEventPlaneMaker from event
    StEventPlaneMaker *EventPlaneMaker[5];

    // pt-dependent bin mode
    const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
    for(int i = 0; i < 5; i++) {
      EventPlaneMaker[i] = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, i)));
    }

    // event plane bin to use: pt dependent ranges
    int ptAssocBins[9] = {0, 1, 2, 3, 4, 4,4,4,4};
    int EPBinToUse = ptAssocBins[assocPtBin];

    // check for requested EventPlaneMaker pointer
    if(!EventPlaneMaker[EPBinToUse] && !doppAnalysis) { 
      LOG_WARN<<Form("No EventPlaneMaker bin: %i!", EPBinToUse)<<endm; 
      if(assocPtBin==0) hStats->Fill(25);
      return; 
    }

    // get event plane angle for different pt bins - assign global event plane to selected pt-dependent bin
    // could also write this as:  tpc2EP_bin = (EventPlaneMaker) ? (double)EventPlaneMaker->GetTPCEP() : -999;
    double tpc2EP  = (EventPlaneMaker[EPBinToUse]) ? (double)EventPlaneMaker[EPBinToUse]->GetTPCEP() : -999;
    double jetV2EP = (EventPlaneMaker[EPBinToUse]) ? (double)EventPlaneMaker[EPBinToUse]->GetTPCEP() : -999;

    // if requiring a single event plane angle (non-pt dependent): use charged tracks 0.2-2.0 GeV to calculate EP
    if(doUseMainEPAngle) {
      tpc2EP  = (EventPlaneMaker[4]) ? (double)EventPlaneMaker[4]->GetTPCEP() : -999;
      jetV2EP = (EventPlaneMaker[4]) ? (double)EventPlaneMaker[4]->GetTPCEP() : -999;
    }
    //cout<<"assocPtBin: "<<assocPtBin<<"  tpc2EP: "<<tpc2EP<<"  jetV2EP: "<<jetV2EP<<endl;

    // get relative angle between jet and event plane
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, tpc2EP) : -99.; // CORRECTED event plane angle - STEP3

    // get relative jet-event plane bin
    int EPBin = (!doppAnalysis) ? GetJetEPBin(dEP) : 0; // this secures it will still run for pp, and use EP bin = 0
    if(EPBin < 0) {
      if(assocPtBin==0) hStats->Fill(26);
      return;
    }

    // calculate jet v2 here, can move around for other analysis, but here for now
    // jetPtr, EPangle, ptBin
    GetJetV2(jet, jetV2EP, assocPtBin);

    // fill some jet histograms
    if(assocPtBin == 0) { // fill only once, so do it for lowest pt assoc bin
      hJetPt->Fill(uncorrjetPt);
      hJetCorrPt->Fill(corrjetPt);
      hJetE->Fill(jetE);
      hJetEta->Fill(jetEta);
      hJetPhi->Fill(jetPhi);
      hJetNEF->Fill(jetNEF);
      hJetArea->Fill(jetArea);
      hJetMass->Fill(jetMass);
      hJetPtvsArea->Fill(uncorrjetPt, jetArea);
    }
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
      if(doUsePrimTracks) { // get primary track vector
        mTrkMom = trk->pMom();
      } else {              // get global track vector
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
        if(assocPtBin == 0) { if((tpt < 0.20) || (tpt  >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for jet shapes
        if(assocPtBin == 1) { if((tpt < 0.50) || (tpt  >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for jet shapes
        if(assocPtBin == 2) { if((tpt < 1.00) || (tpt  >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for jet shapes
        if(assocPtBin == 3) { if((tpt < 1.50) || (tpt  >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for jet shapes
        if(assocPtBin == 4) { if((tpt < 2.00) || (tpt  >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for jet shapes
        if(assocPtBin == 5) { if((tpt < 3.00) || (tpt  >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for jet shapes
        if(assocPtBin == 6) { if((tpt < 4.00) || (tpt  >= 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for jet shapes
        if(assocPtBin == 7) { if((tpt <  6.0))                  continue; }  //       6.0+ GeV assoc bin used for jet shapes
        if(assocPtBin == 8) { if((tpt <  0.5))                  continue; }  //       0.5+ GeV assoc bin used for jet shapes
        //if(thisbin == 4) { if((tpt < 2.00) || (tpt >= 4.0)) continue; } // 2.00 - 4.0 GeV assoc bin used for jet shapes
        //if(thisbin == 4) { if((tpt < 4.00) || (tpt >= 8.0)) continue; } // 4.00 - 8.0 GeV assoc bin used for jet shapes
      }

      // get radial distance between track and jet axis
      double deltaR = GetDeltaR(jet, trk);

      // get annuli bin
      int annuliBin = GetAnnuliBin(deltaR);
      if(annuliBin < 0) { continue; }

      // calculate single particle tracking efficiency
      int effCent   = mCentMaker->GetRef16();
      double fZDCx  = mPicoEvent->ZDCx();
      double trkEff = ApplyTrackingEff(fDoEffCorr, tpt, teta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile);

      // calculate radial pt sum
      rsum[annuliBin] += tpt*(1.0/trkEff);

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

    // BACKGROUND tracks
    // track loop inside jet loop - loop over ALL tracks in PicoDst - for BG
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

      // get momentum vector of track - global or primary track
      TVector3 mTrkMom;
      if(doUsePrimTracks) { // get primary track vector
        mTrkMom = trk->pMom();
      } else {              // get global track vector
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
        if(assocPtBin == 0) { if((tpt < 0.20) || (tpt  >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(assocPtBin == 1) { if((tpt < 0.50) || (tpt  >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(assocPtBin == 2) { if((tpt < 1.00) || (tpt  >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(assocPtBin == 3) { if((tpt < 1.50) || (tpt  >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(assocPtBin == 4) { if((tpt < 2.00) || (tpt  >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(assocPtBin == 5) { if((tpt < 3.00) || (tpt  >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(assocPtBin == 6) { if((tpt < 4.00) || (tpt  >= 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations
        if(assocPtBin == 7) { if((tpt <  6.0))                  continue; }  //       6.0+ GeV assoc bin used for correlations
        if(assocPtBin == 8) { if((tpt <  0.5))                  continue; }  //       0.5+ GeV assoc bin used for correlations
      }

      // get radial distance between track and jet axis
      double deltaEtaBG = 1.0*TMath::Abs(jetEtaBG - teta);
      double deltaPhiBG = 1.0*TMath::Abs(jetPhiBG - tphi);
      if(deltaPhiBG > 1.0*pi) deltaPhiBG = 2.0*pi - deltaPhiBG;
      double deltaR = 1.0*TMath::Sqrt(deltaEtaBG*deltaEtaBG + deltaPhiBG*deltaPhiBG);
      int annuliBin = GetAnnuliBin(deltaR);
      //if(annuliBin < 0) continue;

      // ApplyTrackingEff
      // calculate single particle tracking efficiency
      int effCent   = mCentMaker->GetRef16();
      double fZDCx  = mPicoEvent->ZDCx();
      double trkEff = ApplyTrackingEff(fDoEffCorr, tpt, teta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile);

      // calculate radial pt sum
      //rsumBG[annuliBin] += tpt;
      if(annuliBin >= 0) rsumBG[annuliBin] += tpt*(1.0/trkEff); // use this instance to avoid 2 separate track loops

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

    // Case 3: mixed events!
    // event mixing for background jet cones
    if(fDoEventMixing > 0){
      // initialize background tracks array
      TObjArray *bgTracks;
      TObjArray *bgTracksCheck; // FIXME probably don't need another array to hold the *same* tracks

      // do event mixing when Signal Jet is part of event with a HT1 or HT2 or HT3 trigger firing
      if(pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {
        // get number of current events in pool
        int nMix = pool->GetCurrentNEvents();
        int nMixNorm = pool->GetCurrentNEvents();

        // QA histogram
        if(assocPtBin == 0) hNMixEvents->Fill(nMix);
        hNMixNormBefore[assocPtBin]->Fill(nMixNorm);

        // reset annuli sums here - do this inside loop over mixed events
        //double rsumBG3[10] = {0.0};

        // =============================================================================================================
        if(doGenerateBadMixEventBGcone) {
          // loop over nMix events: find 'bad' mixed events to drop - dropping events where background jet cone reaches threshold fraction
          for(int jMix = 0; jMix < nMix; jMix++) {
            // get jMix'th event
            bgTracksCheck = pool->GetEvent(jMix);
            const Int_t Nbgtrks = bgTracksCheck->GetEntries();

            // sum everything in background cone to compare
            double backgroundPtConeSum = 0.0;

            // loop over background (mixed event) tracks
            for(int ibg = 0; ibg < Nbgtrks; ibg++) {
              // slimmed PicoTrack class: StFemtoTrack
              StFemtoTrack *trk = static_cast<StFemtoTrack*>(bgTracksCheck->At(ibg));
              if(!trk) continue;

              // mixed track variables
              double Mphi = trk->Phi();
              double Meta = trk->Eta();
              double Mpt = trk->Pt();
              //int mMBTrig = trk->MBTrig(); // FIXME don't need to worry much about this for jet shape
            
              // shift angle (0, 2*pi) 
              if(Mphi < 0.0)    Mphi += 2.0*pi;
              if(Mphi > 2.0*pi) Mphi -= 2.0*pi;

              // sum up background cone: (constitCut, MaxTrack) range
              // assocPtBin == 0 && 
              if(Mpt >= fJetConstituentCut && Mpt < fJetShapeTrackPtMax) {
                double dEtaBG = 1.0*TMath::Abs(jetEta - Meta);
                double dPhiBG = 1.0*TMath::Abs(jetPhi - Mphi);
                if(dPhiBG > 1.0*pi) dPhiBG = 2.0*pi - dPhiBG;

                // delta R
                double dRbg = 1.0*TMath::Sqrt(dEtaBG*dEtaBG + dPhiBG*dPhiBG);
                if(dRbg <= fJetRad) backgroundPtConeSum += Mpt;
              }
            } // end of background track loop

            // calculate fraction of jetpt that background cone is - FIXME
            double backgroundFraction = 1.0*(backgroundPtConeSum / jetPt);

            // cut on background fraction level for mixed events
            // change the normalization to use the parameter below:
            if(backgroundFraction > fBackgroundConeFractionCut) {
              // THIS IS WHERE A CUT would need to be done - INSIDE the event loop
              // a *continue* here would skip to the next mixed event
              nMixNorm--;  // this reduces an event from the mixevent normalization
              //  continue;
            }

          }   // end of jth mix event loop
        }  //doGenerateBadMixEventBGcone check
        // =============================================================================================================

        // get z-vtx and corresponding bins
        int zbin = GetZVertex4cmBin(zVtx);
        int ZvtxBin = GetZvtxBin(zVtx);
        if(ZvtxBin < 0) return;

        // loop over nMix events: fill mixed event histos here
        for(int jMix = 0; jMix < nMix; jMix++) {
          // get jMix'th event
          bgTracks = pool->GetEvent(jMix);
          const Int_t Nbgtrks = bgTracks->GetEntries();

          // reset annuli sums here - when NOT normalizing by nMix
          double rsumBG3[10] = {0.0};

          // sum everything in background cone to compare
          double backgroundPtConeSum = 0.0;
          hNMixNormAfter[assocPtBin]->Fill(nMixNorm);

          //============================================================================================
          // centrality & z-vtx correction to event level:
          //	- initialize here
          //	- since all tracks in same mixed event will have same MB trigger from same event (zvtx + RefMultCorr)
          //      get the trigger and therefore weight of first track, keep it until end of track loop when filling histos 
          // mixed track z-vtx weights based on MB5/MB30
          // combined correction factor - zvtx and RefMultCorr dependency
          double fMixEvtWeightCorrFactor = 1.0;
          //=============================================================================================

          // loop over background (mixed event) tracks
          for(int ibg = 0; ibg < Nbgtrks; ibg++) {
            // slimmed PicoTrack class: StFemtoTrack
            StFemtoTrack *trk = static_cast<StFemtoTrack*>(bgTracks->At(ibg));
            if(!trk) continue;

            // mixed track variables
            double Mphi = trk->Phi();
            double Meta = trk->Eta();
            double Mpt = trk->Pt();
            //int mMBTrig = trk->MBTrig(); // FIXME don't need to worry much about this for jet shape

            // shift angle (0, 2*pi) 
            if(Mphi < 0.0)    Mphi += 2.0*pi;
            if(Mphi > 2.0*pi) Mphi -= 2.0*pi;

            // cut on track pt
            if(Mpt < fJetShapeTrackPtMin) { continue; } // min track pt cut
            if(Mpt > fJetShapeTrackPtMax) { continue; } // max track pt cut

            // sum up background cone: (constitCut, MaxTrack) range
            // fJetShapeTrackPtMax would equal the max track used in jet reconstruction, which is why its also used here
            //if(assocPtBin == 0 && 
            if(Mpt >= fJetConstituentCut && Mpt < fJetShapeTrackPtMax) {
              double dEtaBG = 1.0*TMath::Abs(jetEta - Meta);
              double dPhiBG = 1.0*TMath::Abs(jetPhi - Mphi);
              if(dPhiBG > 1.0*pi) dPhiBG = 2.0*pi - dPhiBG;

              // delta R
              double dRbg = 1.0*TMath::Sqrt(dEtaBG*dEtaBG + dPhiBG*dPhiBG);
              if(dRbg <= fJetRad) backgroundPtConeSum += Mpt;
            }

            // additional pt selection when doing pt associated bin method
            if(doTPCptassocBin) {
              if(assocPtBin == 0) { if((Mpt < 0.20) || (Mpt  >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for jet shapes
              if(assocPtBin == 1) { if((Mpt < 0.50) || (Mpt  >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for jet shapes
              if(assocPtBin == 2) { if((Mpt < 1.00) || (Mpt  >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for jet shapes
              if(assocPtBin == 3) { if((Mpt < 1.50) || (Mpt  >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for jet shapes
              if(assocPtBin == 4) { if((Mpt < 2.00) || (Mpt  >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for jet shapes
              if(assocPtBin == 5) { if((Mpt < 3.00) || (Mpt  >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for jet shapes
              if(assocPtBin == 6) { if((Mpt < 4.00) || (Mpt  >= 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for jet shapes
              if(assocPtBin == 7) { if((Mpt <  6.0))                  continue; }  //       6.0+ GeV assoc bin used for jet shapes
              if(assocPtBin == 8) { if((Mpt <  0.5))                  continue; }  //       0.5+ GeV assoc bin used for jet shapes
              //if(thisbin == 4) { if((tpt < 2.00) || (tpt >= 4.0)) continue; } // 2.00 - 4.0 GeV assoc bin used for jet shapes
              //if(thisbin == 4) { if((tpt < 4.00) || (tpt >= 8.0)) continue; } // 4.00 - 8.0 GeV assoc bin used for jet shapes
            }

            // get radial distance between track and jet axis
            double deltaEtaBG = 1.0*TMath::Abs(jetEta - Meta);
            double deltaPhiBG = 1.0*TMath::Abs(jetPhi - Mphi);
            if(deltaPhiBG > 1.0*pi) deltaPhiBG = 2.0*pi - deltaPhiBG;
            double deltaR = 1.0*TMath::Sqrt(deltaEtaBG*deltaEtaBG + deltaPhiBG*deltaPhiBG);

            // get annuli bin
            int annuliBin = GetAnnuliBin(deltaR);
            if(annuliBin < 0) continue;

            //============================================================================================
            // centrality & z-vtx correction to event level: tracks
            // mixed track z-vtx weights based on MB5/MB30
            ////int mMBTrig = trk->MBTrig();      // mixed event trigger: where track came from

            // combined correction factor - zvtx and RefMultCorr dependency
            ////double fMixWeightCorrFactor = trk->ReWeightCorr();

            // get info from first track for scaling/correction of MB5->MB30
            //	- entire mixed event will have same MB trigger (zvtx + RefMultCorr)
            if(ibg == 0)  fMixEvtWeightCorrFactor = trk->ReWeightCorr();

            //=============================================================================================

            // calculate single particle tracking efficiency of mixed events for correlations (-999)
            int effCent   = mCentMaker->GetRef16();
            double fZDCx  = mPicoEvent->ZDCx();
            double mixefficiency = ApplyTrackingEff(fDoEffCorr, Mpt, Meta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile);

            // calculate radial pt sum - 
            //rsumBG3[annuliBin] += Mpt*(1.0/mixefficiency) * fMixWeightCorrFactor;
            rsumBG3[annuliBin] += Mpt*(1.0/mixefficiency);

          } // end of background track loop

          // calculate fraction of jetpt that background cone is
          double backgroundFraction = 1.0*(backgroundPtConeSum / jetPt);
          if(assocPtBin == 0) { // only fill this for one pt associated bin, as it is the EXACT same, because they are mixed *events*
            //cout<<"jMix: "<<jMix<<"   backgroundFraction: "<<backgroundFraction<<"   bgSum: "<<backgroundPtConeSum<<"   jetPt: "<<jetPt<<endl;
            hBGconeFractionOfJetPt->Fill(backgroundFraction);
            hJetPtvsBGconeFraction->Fill(jetPt, backgroundFraction);
            hJetPtvsBGconePt->Fill(jetPt, backgroundPtConeSum);
          }
        
          // cut on background fraction level for mixed events
          // change the normalization to use the parameter below:
          if(backgroundFraction > fBackgroundConeFractionCut && doGenerateBadMixEventBGcone) { 
            // THIS IS WHERE A CUT would need to be done
            // a *continue* here would skip to the next mixed event
            //  nMixNorm--; // implemented in above loop
            //cout<<"Skipping this mixed event because backgroundFraction: "<<backgroundFraction<<endl;
            continue;
          }
 
          // fill BG histos here
          for(int i = 0; i < 10; i++) {
            //hJetShapeBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
            //hJetShapeBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
            //hJetPtProfileBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
            //hJetPtProfileBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));

            hJetShapeBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0 * fMixEvtWeightCorrFactor * rsumBG3[i] / (nMixNorm*jetPt));
            hJetShapeBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0 * fMixEvtWeightCorrFactor * rsumBG3[i] / (nMixNorm*jetPt));
            hJetPtProfileBGCase3[jetPtBin][centBin][EPBin][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0 * fMixEvtWeightCorrFactor * rsumBG3[i] / (nMixNorm));
            hJetPtProfileBGCase3[jetPtBin][centBin][3][assocPtBin]->Fill(i*rbinSize + 1e-3, 1.0 * fMixEvtWeightCorrFactor * rsumBG3[i] / (nMixNorm));
          } // loop over annuli bins

        }   // end of filling mixed-event histo's:  jth mix event loop

        // counter - only need to fill once - think about filling this elsewhere?
        if(assocPtBin == 0) {
          hJetCounterCase3BG[jetPtBin][centBin][EPBin]->Fill(0.5);
          hJetCounterCase3BG[jetPtBin][centBin][3]->Fill(0.5); // ALL angles
        }
      }     // end of check for pool being ready
    }       // end of event mixing

}
//
// function to calculate Jet v2 using the event plane method
//____________________________________________________________________________________________
void StMyAnalysisMaker3::GetJetV2(StJet *jet, Double_t EPangle, Int_t ptAssocBin)
{
  // get centrality bin
  int centBin = (!doppAnalysis) ? Get4CentBin(fCentralityScaled) : 0; // this secures it will still run for pp
  if(centBin < 0) return;

  // get jet pt bin and value
  double jetPt = -99.;
  if(fCorrJetPt) { jetPt = jet->Pt() - jet->Area()*fRhoVal;
  } else { jetPt = jet->Pt(); }
  int jetPtBin = GetJetPtBin(jetPt);
  if(jetPtBin < 0) return;

  // the main function GetJetV2() is called in the jet shape analysis, and the below cuts are already performed
  // require tower and or track bias for jet
  if(doBiasJetLeadConstituent && (jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) {
    return;
  }

  // check that jet contains a tower that fired the trigger
  if(doRequireJetTowFireTrig && !DidTowerConstituentFireTrigger(jet)) {
    return;
  }

  // jet parameters
  double jetPhi = jet->Phi();

  // get event plane bin
  double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, EPangle) : -99.; // CORRECTED event plane angle - STEP3
  int EPBin = (!doppAnalysis) ? GetJetEPBin(dEP) : 0; // this secures it will still run for pp, and use EP bin = 0
  if(EPBin < 0) return;

  // fill histogram
  fProfJetV2[jetPtBin][centBin][EPBin]->Fill(ptAssocBin, TMath::Cos(2.*(jetPhi - EPangle)));
  fProfJetV2[jetPtBin][centBin][3]->Fill(ptAssocBin, TMath::Cos(2.*(jetPhi - EPangle)));

}
//
// Function to add event pools to output file
//____________________________________________________________________
void StMyAnalysisMaker3::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent, Double_t minZvtx, Double_t maxZvtx, Double_t minPsi2, Double_t maxPsi2, Double_t minPt, Double_t maxPt)
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
//
// Jet hadron correlation analysis function
//______________________________________________________________________________________________________________________________
void StMyAnalysisMaker3::JetHadronCorrelationAnalysis(StJet *jet, StEventPool *pool, Int_t centbin, Int_t assocPtBin) {
    // get number of jets, tracks, and global tracks in events
    Int_t njets = fJets->GetEntries();
    const Int_t ntracks = mPicoDst->numberOfTracks();
    Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();
    if(fDebugLevel == kDebugGeneralEvt) {
      //cout<<"grefMult = "<<grefMult<<"  refMult = "<<refMult<<"  refCorr2 = "<<refCorr2<<"  cent16 = "<<cent16<<"   cent9 = "<<cent9<<"  centbin = "<<centbin<<endl;
      /////cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<"  refCorr2 = "<<refCorr2<<"  grefMult = "<<grefMult<<"  centbin = "<<centbin<<endl;
    }

    // bin-age to use for mixed event and sparses - FIXME not set up right yet
    Int_t centbin10 = GetCentBin10(centbin);
    double centBinToUse;
    if(fCentBinSize==10)       { centBinToUse = (double)centbin10 * 10.0;
    } else if(fCentBinSize==5) { centBinToUse = (double)centbin * 5.0; }
    //cout<<"fCentBinSize: "<<fCentBinSize<<"   centBinToUse: "<<centBinToUse<<"   centbin = "<<centbin<<"   centbin10 = "<<centbin10<<endl;

    // to limit filling unused entries in sparse, only fill for certain centrality ranges
    // ranges can be different than functional cent bin setter
    Int_t cbin = -1;
    if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
    else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
    else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
    else if (centbin>5 && centbin<10) cbin = 4; // 30-50%
    else if (centbin>9 && centbin<16) cbin = 5; // 50-80%
    else cbin = -99;
    // ============================================

    // get some jet parameters
    double jetArea = jet->Area();
    double jetMass = jet->M();
    double jetPt = jet->Pt();
    double corrjetPt = jet->Pt() - jetArea*fRhoVal;
    double jetE = jet->E();
    double jetEta = jet->Eta();
    double jetPhi = jet->Phi();
    double jetNEF = jet->NEF();

    // select which jet pt to use
    double jetPtselected;
    if(fCorrJetPt) { jetPtselected = corrjetPt;
    } else { jetPtselected = jetPt; }

    // Cuts are done below, only once as there are no jet loops in this function, the function is passed a jet pointer
    // therefore, same jet used for signal (same event) as for mixed event and thus same cuts! - some threshold cuts
    if(fCorrJetPt) {  // background subtracted jet pt
      if(corrjetPt < fMinPtJet)    return;
    } else { if(jetPt < fMinPtJet) return; }

    // this is a test
/*
    cout<<"MaxJetTrack = "<<jet->GetMaxTrackPt()<<"   MaxJetTow = "<<jet->GetMaxTowerEt()<<endl;
    if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) { cout<<"1. tMax < TBias AND TowMax < TowBias = TRUE, thus RETURN"<<endl; }
    //if((jet->GetMaxTrackPt() < fTrackBias) || (jet->GetMaxTowerEt() < fTowerBias)) { cout<<"2.  tMax < TBias   OR TowMax < TowBias    = TRUE, thus RETURN"<<endl; }
    if(jet->GetMaxTowerEt() > fTowerBias) { cout<<"3. TowMax > TowBias = TRUE, thus I want this!!"<<endl; }
    if(jet->GetMaxTowerEt() < fTowerBias) { cout<<"4. TowMax < TowBias = TRUE, thus RETURN"<<endl; }
    if(! (jet->GetMaxTowerEt() > fTowerBias) ) { cout<<"5. tower did not quality"<<endl; }
    if(DidTowerConstituentFireTrigger(jet)) { cout<<"6. jet tower fired trigger!!"<<endl; }
//    if( ((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias))
    if( !((jet->GetMaxTrackPt() > fTrackBias) || (jet->GetMaxTowerEt() > fTowerBias)) ) { cout<<"7. I don't have anything!! "<<endl; }
*/

    // check for jets with only 1 constituent
    int nJetConstituents = jet->GetNumberOfConstituents();
    if(nJetConstituents == 1) {
      if(assocPtBin==0) hStats->Fill(21);
      if(doSkip1ParticleJets)  return;
    }

    // require tower and or track bias for jet
    if(doBiasJetLeadConstituent && (jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) {
      if(assocPtBin==0) hStats->Fill(23);
      return;
    }

    // check that jet contains a tower that fired the trigger
    if(doRequireJetTowFireTrig && !DidTowerConstituentFireTrigger(jet)) {
      if(assocPtBin==0) hStats->Fill(24);
      return;
    }

    // ==========================================================================================================
    // get StEventPlaneMaker from event
    StEventPlaneMaker *EventPlaneMaker[5];

    // pt-dependent bin mode
    const char *fEventPlaneMakerNameChTemp = fEventPlaneMakerName;
    for(int i = 0; i < 5; i++) {
      EventPlaneMaker[i] = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, i)));
    }

    // event plane bin to use: pt dependent ranges
    int ptAssocBins[9] = {0, 1, 2, 3, 4, 4,4,4,4};
    int EPBinToUse = ptAssocBins[assocPtBin];

    // check for requested EventPlaneMaker pointer
    if(!EventPlaneMaker[EPBinToUse]) { LOG_WARN<<Form("No EventPlaneMaker bin: %i!", EPBinToUse)<<endm; return; }

    // get event plane angle for different pt bins - assign global event plane to selected pt-dependent bin
    // could also write this as:  tpc2EP_bin = (EventPlaneMaker) ? (double)EventPlaneMaker->GetTPCEP() : -999;
    double tpc2EP  = (EventPlaneMaker[EPBinToUse]) ? (double)EventPlaneMaker[EPBinToUse]->GetTPCEP() : -999;

    // if requiring a single event plane angle (non-pt dependent): use charged tracks 0.2-2.0 GeV to calculate EP
    if(doUseMainEPAngle) {
      tpc2EP  = (EventPlaneMaker[4]) ? (double)EventPlaneMaker[4]->GetTPCEP() : -999;
    }

    // get relative angle between jet and event plane
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, tpc2EP) : -99.; // CORRECTED event plane angle - STEP3

    // get relative jet-event plane bin
    int EPBin = (!doppAnalysis) ? GetJetEPBin(dEP) : 0; // this secures it will still run for pp, and use EP bin = 0
    if(EPBin < 0) return;

    // =======================================================================================================================================
    // =======================================================================================================================================

    // fill some jet histograms
    if(assocPtBin == 0) { // fill only once, so do it for lowest pt assoc bin
      hJetPt->Fill(jetPt);
      hJetCorrPt->Fill(corrjetPt);
      hJetE->Fill(jetE);
      hJetEta->Fill(jetEta);
      hJetPhi->Fill(jetPhi);
      hJetNEF->Fill(jetNEF);
      hJetArea->Fill(jetArea);
      hJetMass->Fill(jetMass);
      hJetPtvsArea->Fill(jetPt, jetArea);

      hJetHTrigMaxTowEt->Fill(jet->GetMaxTowerEt());
      hJetHTrigMaxTrkPt->Fill(jet->GetMaxTrackPt());
    }
    hJetEventEP[assocPtBin]->Fill(tpc2EP);
    hJetPhivsEP[assocPtBin]->Fill(jetPhi, tpc2EP);

    // fill some jet QA plots for each orientation
    if(dEP >= 0.0*pi/6.0 && dEP < 1.0*pi/6.0) {
      hJetPtIn[assocPtBin]->Fill(jetPt);
      hJetPhiIn[assocPtBin]->Fill(jetPhi);
      hJetEtaIn[assocPtBin]->Fill(jetEta);
      hJetEventEPIn[assocPtBin]->Fill(tpc2EP);
      hJetPhivsEPIn[assocPtBin]->Fill(jetPhi, tpc2EP);
    } else if(dEP >= 1.0*pi/6.0 && dEP < 2.0*pi/6.0) {
      hJetPtMid[assocPtBin]->Fill(jetPt);
      hJetPhiMid[assocPtBin]->Fill(jetPhi);
      hJetEtaMid[assocPtBin]->Fill(jetEta);
      hJetEventEPMid[assocPtBin]->Fill(tpc2EP);
      hJetPhivsEPMid[assocPtBin]->Fill(jetPhi, tpc2EP);
    } else if(dEP >= 2.0*pi/6.0 && dEP <= 3.0*pi/6.0) {
      hJetPtOut[assocPtBin]->Fill(jetPt);
      hJetPhiOut[assocPtBin]->Fill(jetPhi);
      hJetEtaOut[assocPtBin]->Fill(jetEta);
      hJetEventEPOut[assocPtBin]->Fill(tpc2EP);
      hJetPhivsEPOut[assocPtBin]->Fill(jetPhi, tpc2EP);
    }

    // get nTracks and maxTrackPt
    double maxtrackPt = jet->GetMaxTrackPt();
    int NtrackConstit = jet->GetNumberOfTracks();
    if(fDebugLevel == kDebugJetConstituents) cout<<"JetPt = "<<jetPt<<"  JetE = "<<jetE<<"   MaxtrackPt = "<<maxtrackPt<<"  NtrackConstit = "<<NtrackConstit<<endl;
    
    // ====================================================================================
    // set up and fill jet THnSparse for trigger jet normalization
    Double_t CorrEntries[5] = {centBinToUse, jetPtselected, dEP, zVtx, (double)assocPtBin};
    if(fReduceStatsCent > 0) {
      if(cbin == fReduceStatsCent) fhnCorr->Fill(CorrEntries); // fill Sparse Histo with trigger Jets entries
    } else fhnCorr->Fill(CorrEntries);                         // fill Sparse Histo with trigger Jets entries
    // ======================================================================================

    // track loop inside jet loop - loop over ALL tracks in PicoDst
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk) continue; 

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) continue;

      // get momentum vector of track - global or primary track
      TVector3 mTrkMom;
      if(doUsePrimTracks) { mTrkMom = trk->pMom();                  // primary track vector
      } else {              mTrkMom = trk->gMom(mVertex, Bfield); } // global track vector

      // track variables
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      short charge = trk->charge();

      // shift angle (0, 2*pi)
      if(phi < 0.0)    phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;

      // get jet - track variables
      //deta = eta - jetEta;                    // eta betweeen hadron and jet
      double deta = jetEta - eta;               // eta betweeen jet and hadron
      double dphijh = RelativePhi(jetPhi, phi); // angle between jet and hadron

      // exclusion of track pt ranges from filling - need to use specific EP value
      // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0, 2.0-20.0 GeV - when doing event plane calculation via pt assoc bin
      // when doing event plane calculation via pt assoc bin
      // 	if(assocPtBin == 0) { if(pt < 0.20 || pt >= 0.5) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(doTPCptassocBin) {
        if(assocPtBin == 0) { if(pt < 0.20 || pt >= 0.5) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(assocPtBin == 1) { if(pt < 0.50 || pt >= 1.0) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(assocPtBin == 2) { if(pt < 1.00 || pt >= 1.5) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(assocPtBin == 3) { if(pt < 1.50 || pt >= 2.0) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(assocPtBin == 4) { if(pt < 2.00 || pt >= 20.) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      }

      // calculate single particle tracking efficiency of mixed events for correlations
      int effCent   = mCentMaker->GetRef16();
      double fZDCx  = mPicoEvent->ZDCx();
      double trkEfficiency = ApplyTrackingEff(fDoEffCorr, pt, eta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile);
 
      // fill jet sparse (signal jets correlated with tracks from the same event)
      double triggerEntries[9] = {centBinToUse, jetPtselected, pt, deta, dphijh, dEP, zVtx, (double)charge, (double)assocPtBin};
      if(fReduceStatsCent > 0) {
        if(cbin == fReduceStatsCent) fhnJH->Fill(triggerEntries, 1.0/trkEfficiency);  // fill Sparse Histo with trigger entries
      } else fhnJH->Fill(triggerEntries, 1.0/trkEfficiency);

      fHistJetHEtaPhi->Fill(deta, dphijh); // fill jet-hadron  eta--phi distribution
    } // track loop

    // ***************************************************************************************************************
    // ******************************** Event MIXING *****************************************************************
    // ***************************************************************************************************************
    if(fDoEventMixing) {
      // initialize background tracks array
      TObjArray *bgTracks;

      // check for readiness of using the event pool
      if(pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {
        // get some get parameters of jets for mixing - cuts on jets ALREADY done above
        // same jet as used in signal (same) event
        //double dMixEP = (!doppAnalysis) ? RelativeEPJET(jet->Phi(), TPC_PSI2) : 0; // CORRECTED event plane angle - STEP3

        // get number of current mixing events in pool
        int nMix = pool->GetCurrentNEvents();

        // QA histogram
        if(assocPtBin == 0) hNMixEvents->Fill(nMix);
        //double sumMixTrkPt = 0.0;

        //----------------------==---------------------==
        // correction factors generated for MB5 and MB30 mixed event triggered events
        // bins 4-17 (-28 cm < z < +28 cm)
        // === z-vtx event weight ===
        double mCorrEvtMB5[] = {0.,0.,0., 187.641, 186.714, 54.9565, 9.34451, 2.15293, 0.509046, 0.279472, 0.375596, 1.2418, 8.87955, 55.0963, 219.312, 222.645, 100.253, 0.,0.,0.};
        double mCorrEvtMB30[]= {0.,0.,0., 0.965911, 0.944671, 0.961855, 0.969395, 0.992681, 1.00764, 1.01291, 1.01672, 1.02394, 1.01706, 1.01427, 1.02102, 0.995457, 0.933781, 0.,0.,0. };
        // === z-vtx pair weight ===
        double mCorrPairMB5[]  = {0.,0.,0., 77.7527, 118.323, 19.5824, 3.22157, 1.10644, 0.618831, 0.57329, 0.596431, 0.787707, 4.03378, 22.5516, 83.528, 60.2196, 61.6985, 0.,0.,0.};
        double mCorrPairMB30[] = {0.,0.,0., 0.487022, 0.485908, 0.496887, 0.5761, 0.906906, 2.91743, 4.85915, 3.59698, 1.40364, 0.554785, 0.495122, 0.486797, 0.487971, 0.48787, 0.,0.,0.};

        // refMultCorr weight
        double corrWeightMB5[400] = { 1, 1, 1, 1, 0.463132, 0.447399, 0.456894, 0.457888, 0.467025, 0.467753, 0.477598, 0.481366, 0.484178, 0.493661, 0.496542, 0.49247, 0.509866, 0.510392, 0.521268, 0.518482, 0.518401, 0.526676, 0.532358, 0.533281, 0.530995, 0.550151, 0.549585, 0.554638, 0.567733, 0.56153, 0.573564, 0.572213, 0.585449, 0.577954, 0.585248, 0.59389, 0.597588, 0.611146, 0.609299, 0.617343, 0.628944, 0.615209, 0.617127, 0.625614, 0.633712, 0.639474, 0.647667, 0.648651, 0.656923, 0.656513, 0.673715, 0.671305, 0.677296, 0.673676, 0.695857, 0.687348, 0.684638, 0.685504, 0.706412, 0.703819, 0.71537, 0.710097, 0.725844, 0.723751, 0.72968, 0.737968, 0.730032, 0.755147, 0.734857, 0.747815, 0.761912, 0.774854, 0.783057, 0.776158, 0.78032, 0.791162, 0.809159, 0.7881, 0.803188, 0.81818, 0.801941, 0.8227, 0.820565, 0.827572, 0.842145, 0.853173, 0.843049, 0.866117, 0.850376, 0.85609, 0.885141, 0.89255, 0.916774, 0.899683, 0.89318, 0.895233, 0.907137, 0.905279, 0.912297, 0.915815, 0.948379, 0.962052, 0.957042, 0.941978, 0.959425, 0.981227, 0.969597, 0.990494, 0.977178, 0.970842, 0.977962, 0.989116, 0.984259, 1.00604, 1.00933, 1.03903, 1.02118, 1.01738, 1.04185, 1.06259, 1.06095, 1.06253, 1.08362, 1.08336, 1.07505, 1.10287, 1.10887, 1.13512, 1.11754, 1.1389, 1.12494, 1.14724, 1.12168, 1.13811, 1.18041, 1.14112, 1.16631, 1.16014, 1.1974, 1.18266, 1.21974, 1.20836, 1.20947, 1.24196, 1.2366, 1.30608, 1.25214, 1.27501, 1.29812, 1.29725, 1.32099, 1.30637, 1.31007, 1.35369, 1.32903, 1.34519, 1.35968, 1.38831, 1.37858, 1.40441, 1.37679, 1.42124, 1.41042, 1.44202, 1.43458, 1.51292, 1.42363, 1.46818, 1.47146, 1.4959, 1.51551, 1.54003, 1.54928, 1.56103, 1.56671, 1.55361, 1.61416, 1.63908, 1.6406, 1.62072, 1.63858, 1.64112, 1.66592, 1.6472, 1.65476, 1.69166, 1.72988, 1.70032, 1.78666, 1.73831, 1.77935, 1.77891, 1.82334, 1.84147, 1.89095, 1.81703, 1.81755, 1.84787, 1.88709, 1.98063, 1.8915, 1.94406, 1.95225, 2.01096, 1.9701, 1.96114, 2.05686, 2.04201, 2.07128, 2.05617, 2.0701, 2.13762, 2.1401, 2.16123, 2.19201, 2.17974, 2.17145, 2.21869, 2.22246, 2.19042, 2.20791, 2.23285, 2.26786, 2.29109, 2.31364, 2.32266, 2.29967, 2.30896, 2.46165, 2.37093, 2.41638, 2.4304, 2.43803, 2.46856, 2.47553, 2.48636, 2.53278, 2.44217, 2.55204, 2.53229, 2.57215, 2.60888, 2.61448, 2.63351, 2.63772, 2.64363, 2.69693, 2.60885, 2.74891, 2.66265, 2.73691, 2.81716, 2.73062, 2.74717, 2.81193, 2.89897, 2.76424, 2.92959, 2.88237, 2.86829, 2.92863, 3.07455, 3.12348, 2.9811, 2.99424, 2.88269, 3.03216, 3.13267, 3.09151, 3.00951, 3.15617, 3.14428, 3.17672, 3.10598, 3.28992, 3.30476, 3.20371, 3.25868, 3.15695, 3.32213, 3.38476, 3.24909, 3.24376, 3.21991, 3.18219, 3.31691, 3.4702, 3.42831, 3.68671, 3.69256, 3.35984, 3.53007, 3.05662, 3.50486, 3.23201, 3.62708, 3.49475, 3.70496, 3.22615, 3.1339, 3.82289, 3.82499, 3.49086, 4.21653, 3.7773, 3.94624, 3.45887, 4.11111, 3.96686, 3.73737, 3.47863, 2.50356, 2.14591, 2.89886, 4.51771, 2.24003, 3.68946, 6.32479, 3.86515, 2.63533, 6.85185, 6.32479, 4.74359, 3.16239, 1, 4.21653, 1.5812, 1, 1, 1, 1, 0.263533, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

        // get z-vtx and corresponding bins
        int zbin = GetZVertex4cmBin(zVtx);
        int ZvtxBin = GetZvtxBin(zVtx);
        if(ZvtxBin < 0) return;
        //----------------------==---------------------==

        // fill mixed-event histos here: loop over nMix events
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

            // mixed track variables
            double Mixphi = trk->Phi();       // mixed event track phi
            double Mixeta = trk->Eta();     
            double Mixpt = trk->Pt();
            short Mixcharge = trk->Charge();

            // shift angle (0, 2*pi)
            if(Mixphi < 0.0)    Mixphi += 2.0*pi;
            if(Mixphi > 2.0*pi) Mixphi -= 2.0*pi;

            // get jet - track relations
            //double deta = eta - jetEta;                   // eta betweeen hadron and jet
            double dMixeta = jetEta - Mixeta;               // eta betweeen jet and hadron
            double dMixphijh = RelativePhi(jetPhi, Mixphi); // angle between jet and hadron

            // print tracks outside of acceptance somehow and track characteristics
            if(fDebugLevel == kDebugMixedEvents) {
              cout<<"itrack = "<<ibg<<"  phi = "<<Mixphi<<"  eta = "<<Mixeta<<"  pt = "<<Mixpt<<"  q = "<<Mixcharge<<endl;
              if(TMath::Abs(dMixeta > 1.6)) cout<<"DELTA ETA out of bounds... deta = "<<dMixeta<<"   iTrack = "<<ibg<<"  jetEta = "<<jetEta<<"  trk eta = "<<Mixeta<<endl;
            }

            // exclusion of track pt ranges from filling - need to use specific EP value
            // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0, 2.0-20.0 - when doing event plane calculation via pt assoc bin
            //        if(assocPtBin == 0) { if(pt < 0.20 || pt >= 0.5) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
            if(doTPCptassocBin) {
              if(assocPtBin == 0) { if(Mixpt < 0.20 || Mixpt >= 0.5) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
              if(assocPtBin == 1) { if(Mixpt < 0.50 || Mixpt >= 1.0) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
              if(assocPtBin == 2) { if(Mixpt < 1.00 || Mixpt >= 1.5) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
              if(assocPtBin == 3) { if(Mixpt < 1.50 || Mixpt >= 2.0) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
              if(assocPtBin == 4) { if(Mixpt < 2.00 || Mixpt >= 20.) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
            }

            //============================================================================================
            // centrality & z-vtx correction to event level: tracks

            // mixed track z-vtx weights based on MB5/MB30
            int mMBTrig = trk->MBTrig();      // mixed event trigger: where track came from
            double mMBTrigZWeight = 1.;
            if(mMBTrig ==  5) mMBTrigZWeight = mCorrEvtMB5[zbin];  // mCorrPairMB5[zbin];   // weight for MB5:  corr(z)
            if(mMBTrig == 30) mMBTrigZWeight = mCorrEvtMB30[zbin]; // mCorrPairMB30[zbin];  // weight for MB30: corr(z)

            // combined correction factor
            double fMixWeightCorrFactor = trk->ReWeightCorr();

            //=============================================================================================

            // calculate single particle tracking efficiency of mixed events for correlations (-999)
            int effCent   = mCentMaker->GetRef16();
            double fZDCx  = mPicoEvent->ZDCx();
            double mixtrkEfficiency = ApplyTrackingEff(fDoEffCorr, Mixpt, Mixeta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile);

            // create / fill mixed event sparse  (centbin*5.0) - (signal jets correlated with tracks from mixed events)
            double triggerEntries[9] = {centBinToUse, jetPtselected, Mixpt, dMixeta, dMixphijh, dEP, zVtx, (double)Mixcharge, (double)assocPtBin};
            if(fReduceStatsCent > 0) {
              if(cbin == fReduceStatsCent) fhnMixedEvents->Fill(triggerEntries, 1./(nMix*mixtrkEfficiency) * fMixWeightCorrFactor);
            } else fhnMixedEvents->Fill(triggerEntries, 1./(nMix*mixtrkEfficiency) * fMixWeightCorrFactor);   // fill Sparse histo of mixed events

            // testing QA for trigger weights - Jan 2021
            //if(assocPtBin == 0) {
              bool fEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);              

              if(mMBTrig ==  5) hNPairsvsZvtxMB5->Fill(zVtx,  1./(nMix*mixtrkEfficiency));
              if(mMBTrig == 30) hNPairsvsZvtxMB30->Fill(zVtx, 1./(nMix*mixtrkEfficiency));
              if(fEmcTrigger)   hNPairsvsZvtxHT2->Fill(zVtx,  1./(nMix*mixtrkEfficiency));
              if(mMBTrig ==  5) hNPairsvsZvtxMB5Wt->Fill(zVtx,  mMBTrigZWeight * 1./(nMix*mixtrkEfficiency));
              if(mMBTrig == 30) hNPairsvsZvtxMB30Wt->Fill(zVtx, mMBTrigZWeight * 1./(nMix*mixtrkEfficiency));
            //}

          } // end of background track loop
        }   // end of filling mixed-event histo's:  jth mix event loop
      }     // end of check for pool being ready
    }       // end of fDoEventMixing
 
}
//
// fill trigger ids of the dataset into a histogram
// only currently set up for Run12 pp (200 GeV) and Run14 AuAu (200 GeV) 
// - set up for more runs
//________________________________________________________________________
void StMyAnalysisMaker3::FillTriggerIDs(TH1 *h) {
  // All non-test triggers for Run12 Run14

  // Run14 AuAu (200 GeV) - 51, 0-50
  unsigned int triggersRun14[] = {440001, 440004, 440005, 440006, 440007, 440015, 440016, 440017, 440050, 440061, 440064, 450005, 450008, 450009, 450010, 450011, 450012, 450013, 450014, 450015, 450018, 450020, 450021, 450023, 450024, 450025, 450050, 450060, 450103, 450201, 450202, 450203, 450211, 450212, 450213, 450600, 450601, 460001, 460002, 460003, 460005, 460007, 460012, 460101, 460102, 460111, 460201, 460202, 460203, 460212, 490016};

  // Run12 pp (200 GeV) - 27, 0-26
  unsigned int triggersRun12[] = {370001, 370011, 370021, 370022, 370031, 370032, 370301, 370341, 370361, 370501, 370511, 370521, 370522, 370531, 370541, 370542, 370546, 370601, 370611, 370621, 370641, 370701, 370801, 370980, 370981, 370982, 370983};

  // get size of trigger ID arrays:
  size_t nRun12IDs = sizeof(triggersRun12)/sizeof(triggersRun12[0]);
  size_t nRun14IDs = sizeof(triggersRun14)/sizeof(triggersRun14[0]);
  int nLoopMax = 0;
  if(StJetFrameworkPicoBase::Run12_pp200)      nLoopMax = nRun12IDs;
  if(StJetFrameworkPicoBase::Run14_AuAu200)    nLoopMax = nRun14IDs;
  if(StJetFrameworkPicoBase::Run14_AuAu200_MB) nLoopMax = nRun14IDs;

  // label bins of the analysis trigger selection summary
  for(int i = 0; i < nLoopMax; i++) {
    if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)      h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggersRun12[i]));
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200)    h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggersRun14[i]));
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggersRun14[i]));
  }

  // set x-axis labels vertically
  h->LabelsOption("v");

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  for(unsigned int i = 0; i < mytriggers.size(); i++) {
    // check for valid, non-test trigger ID
    if(mytriggers[i] > 1000) {
      for(int j = 0; j < nLoopMax; j++) {
        if(mytriggers[i] == triggersRun12[j] && fRunFlag == StJetFrameworkPicoBase::Run12_pp200)      h->Fill(j + 1);
        if(mytriggers[i] == triggersRun14[j] && fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200)    h->Fill(j + 1);
        if(mytriggers[i] == triggersRun14[j] && fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) h->Fill(j + 1);
      } // loops over ID's
    }   // non-test trigger
  }     // loop over triggers

}
// Returns correction for tracking efficiency
//
//Double_t StJetFrameworkPicoBase::ApplyTrackingEffpp(StPicoTrack *trk)
//Double_t StJetFrameworkPicoBase::ApplyTrackingEffAuAu(StPicoTrack *trk)
//____________________________________________________________________________________________
Double_t StMyAnalysisMaker3::ApplyTrackingEff(Bool_t applyEff, Double_t tpt, Double_t teta, Int_t cbin, Double_t ZDCx, Int_t effType, TFile *infile)
{
  // initialize effieciency - check we want to apply it
  double trkEff = 1.0;
  if(!applyEff) return trkEff;

  // x-variable = track pt, y-variable = track eta
  double x = tpt;
  double y = teta;
  //y = TMath::Abs(teta);
  double pi = 1.0*TMath::Pi();
  double effBinContent = -99; // value extracted from histogram
  const char *species =  "pion"; // FIXME
  int lumiBin = GetLuminosityBin(ZDCx);

  //const int nspecies = 3;
  //std::string species[nspecies] = {"pion", "kaon", "proton"};
  //std::string leg_species[nspecies] = {"#pion", "K", "p"};
  //std::string particles[nparticles] = {"piplus", "piminus", "kaonplus", "kaonminus", "proton", "antiproton"};
  //std::string leg_particles[nparticles] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"};
  //hTrack_Efficiency_pTEta[is][ic][il] = new TH2D(Form("hTrack_%s_Efficiency_pTEta_%s_centbin%d", species[is].c_str(), lumi[il].c_str(), ic), "", 100, 0, 10, 10, -1, 1);

  // ========== AuAu Run14 ===========
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200 || fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB) {
    // cut on pt: 0.2 - 30.0 GeV
    // cut on eta: -1.0 < eta < 1.0
    // pt is flat for AuAu above 5.0 (4.5) GeV
    if(x > 4.5) x = 4.5;  // for pt > 4.5 use value at 4.5
    int effBin = -99;

    // 2-D pt/eta dependent efficiency
    if(effType == StJetFrameworkPicoBase::kNormalPtEtaBased) {
      // get 2D pTEta efficiency histo
      //const char *histName = Form("hTrack_%s_Efficiency_pTEta_centbin%d_lumibin%d", species, cbin, lumiBin);
      const char *histName = Form("hTrack_%s_Efficiency_pTEta_final_centbin%d_lumibin%d", species, cbin, lumiBin);
      TH2F *h = static_cast<TH2F*>(infile->Get(Form("%s", histName)));
      if(!h) cout<<"don't have requested histogram! "<<Form("%s", histName)<<endl;
      h->SetName(Form("%s", histName));

      // get efficiency 
      effBin = h->FindBin(x, y); // pt, eta
      effBinContent = h->GetBinContent(effBin);
      //double effBinContentErr = h->GetBinError(effBin);
      //cout<<"effBinContent: "<<effBinContent<<"   effBinErr: "<<effBinContentErr<<endl;

      // delete histo and close input file
      delete h;
    }

    // test statement
    //cout<<"efficiency: "<<effBinContent<<"  pt: "<<x<<"  eta: "<<y<<"   cbin: "<<cbin<<"  lumi: "<<ZDCx<<endl;
  } // Run14 AuAu

  // =========================================================================================================
  // =========================================================================================================
  // =========  Run12 pp tracking efficiency: ========
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) { // the below is probably true of all pp datasets 
    // cut on pt: 0.2 - 30.0 GeV
    // cut on eta: -1.0 < eta < 1.0
    // pt is flat for pp above 2.0 (1.8) GeV 
    if(x > 1.8) x = 1.8;  // for pt > 1.8 use value at 1.8

    // some of the below is not used for the current pp parametrization
    int effBin = -99;

    // 2-D pt/eta dependent efficiency
    if(effType == StJetFrameworkPicoBase::kNormalPtEtaBased) {
      // get 2D pTEta efficiency histo
      const char *histName = Form("hppRun12_PtEtaEfficiency_data_aacuts");

      // changed from double to float
      TH2F *h = static_cast<TH2F*>(infile->Get(Form("%s", histName)));
      if(!h) cout<<"don't have requested histogram! "<<Form("%s", histName)<<endl;
      h->SetName(Form("%s", histName));

      // get efficiency 
      effBin = h->FindBin(x, y); // pt, eta
      effBinContent = h->GetBinContent(effBin);
      double effBinContentErr = h->GetBinError(effBin);
      //cout<<"effBinContent: "<<effBinContent<<"   effBinErr: "<<effBinContentErr<<endl;

      // delete histo and close input file
      delete h;
    }

    // test statement
    //cout<<"efficiency: "<<effBinContent<<"  pt: "<<x<<"  eta: "<<y<<"  lumi: "<<ZDCx<<endl;
  } // Run12 AuAu

  // track efficiency has pt, eta, and centrality dependence
  // RunFlag switch
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 :   // Run12 pp        
        // pp function is flat (and reliable per Raghav) from 2.0 - 30.0 GeV
        trkEff = effBinContent;
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        trkEff = effBinContent;
        break;

    case StJetFrameworkPicoBase::Run14_AuAu200_MB : // Run14 AuAu - HFT dataset (MB triggers)
        trkEff = effBinContent;
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        // don't have parametrization for Run16
        trkEff = 1.0;
        break;

  } // RunFlag switch

  // return the single track reconstruction efficiency for the corresponding dataset
  if(fSysUncType == kDoNothing) { trkEff = 1.00*trkEff; } 
  if(fSysUncType == kTrkEffMin) { trkEff = 0.95*trkEff; }
  if(fSysUncType == kTrkEffMax) { trkEff = 1.05*trkEff; }
  return trkEff;
}
// 
// Jet QA function
//____________________________________________________________________________________________________________________________________________
void StMyAnalysisMaker3::RunJetQA()
{
  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();
  const Int_t ntracks = mPicoDst->numberOfTracks();
  Int_t nglobaltracks = mPicoEvent->numberOfGlobalTracks();
  if(fDebugLevel == kDebugGeneralEvt) {
    cout<<"njets = "<<njets<<"  ntracks = "<<ntracks<<"  nglobaltracks = "<<nglobaltracks<<endl;
  }

  // ====================== Jet loop below ============================
  // loop over Jets in the event: initialize some parameter variables
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
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
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, TPC_PSI2) : 0; // CORRECTED event plane angle - STEP3 - FIXME angle needs updating

    // this is a double check - should not happen, but if it does -> kill the job so it can be fixed - most likely in the runPico macro
    if(dEP < -900) return;

    // some threshold cuts
    if(fCorrJetPt) {  // background subtracted jet pt
      if(corrjetpt < fMinPtJet) continue;
    } else { if(jetpt < fMinPtJet) continue; }

    // apply leading bias to jet
    if(jet->GetMaxTrackPt() < fTrackBias && 
       jet->GetMaxTowerEt() < fTowerBias) 
       continue;

    // TODO check that jet contains a tower that fired the trigger
    if(!DidTowerConstituentFireTrigger(jet)) { continue; }

    // =======================================================================================================================================

    // loop over constituent tracks
    for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++) {
      int trackid = jet->TrackAt(itrk);

      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackid));
      if(!trk){ continue; }

      // track momentum vector
      TVector3 mTrkMom;
      if(doUsePrimTracks) {
        if(!(trk->isPrimary())) return; // check if primary
        mTrkMom = trk->pMom();                // primary track vector
      } else {
        mTrkMom = trk->gMom(mVertex, Bfield); // global track vector
      }

      // track variables
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();
      double pt = mTrkMom.Perp();
      double jetZ = jet->GetZ(mTrkMom.x(), mTrkMom.y(), mTrkMom.z());

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

      // get tower pointer
      StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(ArrayIndex));
      if(!tow){ continue; }

      //int towID = tow->id(); // ArrayIndex = towID - 1 because of array element numbering different than ids which start at 1
      // tower ID: get from index of array shifted by +1
      int towID = ArrayIndex + 1;
      if(towID < 0) continue;

      int containsTower = jet->ContainsTower(ArrayIndex);
      cout<<">= 0: "<<containsTower<<"  itow = "<<itow<<"  id = "<<towID<<"  ArrIndex = "<<ArrayIndex<<"  towE = "<<tow->energy()<<endl;
    } // constituent tower loop
*/

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

    // fill some jet histograms
    //hJetPt->Fill(jetpt);
    //hJetCorrPt->Fill(corrjetpt);
    //hJetE->Fill(jetE);
    //hJetEta->Fill(jetEta);
    //hJetPhi->Fill(jetPhi);
    //hJetNEF->Fill(jetNEF);
    //hJetArea->Fill(jetarea);
    //hJetMass->Fill(jetMass);
    //hJetPtvsArea->Fill(jetpt, jetarea);
  } // jet loop

}

//
// write Jet event plane QA histograms
//______________________________________________________________________________
void StMyAnalysisMaker3::SetupMixEvtPool() {
  // set up centrality bins for mixed events, for pp we need mult bins for event mixing. 
  // Create binning here.
  Int_t nCentBinsTest = 8;
  Double_t *centralityBinsTest = GenerateFixedBinArray(nCentBinsTest, 0., 8.);

/*
Centrality_def_grefmult.txt               - OLD (16)
 10 14 21 29 40 54 | 71 92 116 145 179 218 | 263 315 | 373 441  
Centrality_def_grefmult_P17id_VpdMB30.txt - (16) - P17id
 10 15 22 31 43 58 | 76 97 123 154 189 230 | 276 329 | 390 459  
Centrality_def_grefmult_P18ih_VpdMB30.txt - Aug2019 (16) - P18ih low/mid/high lumi
 10 14 21 29 40 54 | 71 91 115 143 176 214 | 257 307 | 364 430
*/

  // for AuAu data
  // +1 to accomodate the fact that we define bins rather than array entries
  const int nMultBinsAuAu = 26;  // Alt-2 - 27 values, 26 ranges
  Double_t multBinsAuAu[nMultBinsAuAu + 1] = {0, 10,15,21,31,42,53,66,   80, 95, 112, 130, 149, 169, 190, 212, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469, 800};
  // TEST August 23, 2019
//  const int nMultBinsAuAu = 17; // 18 values, 17 ranges
//  Double_t multBinsJS[nMultBinsAuAU + 1] = {0, 10, 14, 21, 29, 40, 54, 71, 91, 115, 143, 176, 214, 257, 307, 364, 430, 800};
  Double_t *multiplicityBinsAuAu = multBinsAuAu;

  // cent bins for AuAu data
  //Int_t nCentBinsAuAu = 20;
  Int_t nCentBinsAuAu;
  if(fCentBinSize == 5) {          nCentBinsAuAu = 20;
  } else if(fCentBinSize == 10) {  nCentBinsAuAu = 10;
  }

  //Double_t *centralityBinsAuAu = GenerateFixedBinArray(nCentBinsAuAu, 0., 20.); 
  Double_t *centralityBinsAuAu = GenerateFixedBinArray(nCentBinsAuAu, 0., (double)nCentBinsAuAu);
  //-------------------------------

  // for pp data
  const int nMultBinspp = 7;
  //Double_t multBinspp[nMultBinspp + 1] = {0.0, 4., 9, 15, 25, 35, 55, 100.0, 500.0};  // 8 (9)
  Double_t multBinspp[nMultBinspp + 1] = {0.0, 4.0, 6.0, 8.0, 10.0, 13.0, 30., 100.};   // 7 (8)
  Double_t *multiplicityBinspp = multBinspp;

  // z-vertex bins for mixed events
  Int_t nZvBins  = 20; // 4 cm wide, 40 for 2 cm wide
  Double_t *zvBins = GenerateFixedBinArray(nZvBins, -40., 40.); // min/max doesn't matter as data is cut zmin/zmax

  // event plane bins for mixed events
  Int_t nEPBins = 6;   // default 6 from 0-180 degrees, (0 - pi) range
  if(fnEPBins != 6) nEPBins = fnEPBins;
  Double_t *epBins = GenerateFixedBinArray(nEPBins, 0., 1.0*pi);

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implementation of AliEventPoolManager
  if(doIgnoreExternalME) {
    if(doJetShapeAnalysis) { // jet shape analysis setup
      if(fDoUseMultBins) {
        if(!doppAnalysis && !doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsAuAu, (Double_t*)multiplicityBinsAuAu, nZvBins, (Double_t*)zvBins); // not pp
        if(!doppAnalysis &&  doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsAuAu, (Double_t*)multiplicityBinsAuAu, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins); // not pp
        if( doppAnalysis) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinspp, (Double_t*)multiplicityBinspp, nZvBins, (Double_t*)zvBins); // is pp

      } else { // centrality binning  
        if(!doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsAuAu, (Double_t*)centralityBinsAuAu, nZvBins, (Double_t*)zvBins);
        if( doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsAuAu, (Double_t*)centralityBinsAuAu, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins);
      }

    } else { // correlation analysis setup
      if(fDoUseMultBins) {
        if(!doppAnalysis && !doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsAuAu, (Double_t*)multiplicityBinsAuAu, nZvBins, (Double_t*)zvBins); // not pp
        if(!doppAnalysis &&  doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsAuAu, (Double_t*)multiplicityBinsAuAu, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins); // not pp
        if( doppAnalysis) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinspp, (Double_t*)multiplicityBinspp, nZvBins, (Double_t*)zvBins); // is pp

      } else { // centrality binning - updated cent pars June13 2019
        if(!doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsAuAu, (Double_t*)centralityBinsAuAu, nZvBins, (Double_t*)zvBins);
        if( doUseEPBins) fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBinsAuAu, (Double_t*)centralityBinsAuAu, nZvBins, (Double_t*)zvBins, nEPBins, (Double_t*)epBins);
      }
    }   // correlations
  }     // mixed-event setup

  // ================================================================================================
  // External Mixed Event object
  // TEMP testing 
  if(!doIgnoreExternalME) {
    Int_t nPtBins = 1;
    Double_t defaultPtBins[2] = {-9999., 9999.};
    Double_t *ptBins = defaultPtBins;

    // check if fPoolMgr is set, else construct object
    if(!fPoolMgr) {
      cout<<"Don't have pre-existing fPoolMgr object, We will now create it!"<<endl;
      fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsAuAu, (Double_t*)multiplicityBinsAuAu, nZvBins, (Double_t*)zvBins, nEPBins, epBins, nPtBins, ptBins);
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

}
//
// function to converte "zvertex" to a bin
//__________________________________________________________________________________
Int_t StMyAnalysisMaker3::GetZvtxBin(Double_t zvertex) const
{
  // cut on +/- 30 cm
  if(TMath::Abs(zvertex) >= 30.0) return -99;

  // initialize z-vtx bin
  Int_t zBin = int((zvertex + 30.) / 2.);  // bin width is equal to 2 centimeters 

  return zBin;
}
