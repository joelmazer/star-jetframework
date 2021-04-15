// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This code is set as an AnalysisMaker task, where it can perform:
// 1) jet analysis
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
// 4) leading + subleading jets
//
// ################################################################

#include "StJetShapeAnalysis.h"
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
#include "StEventPlaneMaker.h"
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StEventPoolManager.h"
#include "StFemtoTrack.h"
#include "StCentMaker.h"

// old file kept
#include "StPicoConstants.h"

ClassImp(StJetShapeAnalysis)

//
//______________________________________________________________________________________
StJetShapeAnalysis::StJetShapeAnalysis(const char *name, StPicoDstMaker *picoMaker, const char *outName = "", bool mDoComments = kFALSE, double minJetPt = 1.0, double trkbias = 0.15, const char *jetMakerName = "", const char *rhoMakerName = "")
  : StJetFrameworkPicoBase(name)
{
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  doJetShapeAnalysis = kFALSE;
  fJetAnalysisJetType = kLeadingJets; // see header, LJ, SubLJ, inclusive
  doRequireAjSelection = kFALSE;
  fCorrJetPt = kFALSE;
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doWriteTrackQAHist = kTRUE;
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerEt = 1000.0; // 30.0
  doRejectBadRuns = kFALSE;
  fLeadingJet = 0x0; fSubLeadingJet = 0x0; fExcludeLeadingJetsFromFit = 1.0;
  fTrackWeight = 1; //StJetFrameworkPicoBase::kPtLinear2Const5Weight; // see StJetFrameworkPicoBase::EPtrackWeightType 
  fEventPlaneMaxTrackPtCut = 2.0;
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
  fJets = 0x0;
  fRunNumber = 0;
  fEPTPCn = 0.; fEPTPCp = 0.; fEPTPC = 0.; fEPBBC = 0.; fEPZDC = 0.;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutName = outName;
  mOutNameEP = "";
  mOutNameQA = "";
  fDoEffCorr = kFALSE;
  fTrackEfficiencyType = StJetFrameworkPicoBase::kNormalPtEtaBased;
  doTPCptassocBin = kFALSE;
  fTPCptAssocBin = -99;
  fMinPtJet = minJetPt;
  fJetConstituentCut = 2.0;
  fTrackBias = trkbias;
  fTowerBias = 0.2;
  fJetRad = 0.4;
  fJetShapeTrackPtMin = 0.2;  fJetShapeTrackPtMax = 30.0;
  fJetShapePtAssocBin = 4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2;     fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0;    fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0;   fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;       fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fTowerEMinCut = 0.2;      fTowerEMaxCut = 100.0;
  fTowerEtaMinCut = -1.0;   fTowerEtaMaxCut = 1.0;
  fTowerPhiMinCut = 0.0;    fTowerPhiMaxCut = 2.0*TMath::Pi();
  fDoEventMixing = 0;  fMixingTracks = 50000;  fNMIXtracks = 5000;  fNMIXevents = 5;
  fCentBinSizeJS = 5;  fReduceStatsCent = -1;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
  //mVertex = 0x0;
  zVtx = 0.0;
  fDoFilterPtMixEvents = kFALSE;
  fEmcTriggerEventType = 0; fMBEventType = 2; fMixingEventType = 0;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4800; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }
  doComments = mDoComments;
  mBaseMaker = 0x0;
  mCentMaker = 0x0;
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;
  fEventPlaneMakerName = "";
}
//
//__________________________________________________________________________________________
StJetShapeAnalysis::~StJetShapeAnalysis()
{ /*  */
  // destructor
  if(hEventPlane)    delete hEventPlane;
  if(hEventZVertex)  delete hEventZVertex;
  if(hCentrality)    delete hCentrality;
  if(hMultiplicity)  delete hMultiplicity;
  if(hRhovsCent)     delete hRhovsCent;
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
  if(hJetTracksPt)       delete hJetTracksPt;
  if(hJetTracksPhi)      delete hJetTracksPhi;
  if(hJetTracksEta)      delete hJetTracksEta;
  if(hJetTracksZ)        delete hJetTracksZ;
  if(hJetPtvsArea)       delete hJetPtvsArea;
  if(hJetEventEP)        delete hJetEventEP;
  if(hJetPhivsEP)        delete hJetPhivsEP;

  if(fHistEventSelectionQA)          delete fHistEventSelectionQA;
  if(fHistEventSelectionQAafterCuts) delete fHistEventSelectionQAafterCuts;
  if(hTriggerIds)            delete hTriggerIds;
  if(hEmcTriggers)           delete hEmcTriggers;
  if(hMixEvtStatZVtx)        delete hMixEvtStatZVtx;
  if(hMixEvtStatCent)        delete hMixEvtStatCent;
  if(hMixEvtStatZvsCent)     delete hMixEvtStatZvsCent;
  if(hTriggerEvtStatZVtx)    delete hTriggerEvtStatZVtx;
  if(hTriggerEvtStatCent)    delete hTriggerEvtStatCent;
  if(hTriggerEvtStatZvsCent) delete hTriggerEvtStatZvsCent;
  if(hMBvsMult)   delete hMBvsMult;
  if(hMB5vsMult)  delete hMB5vsMult;
  if(hMB30vsMult) delete hMB30vsMult;
  if(hHTvsMult)   delete hHTvsMult;
  if(hNMixEvents) delete hNMixEvents;

  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        if(hJetShape[k][j][i])            delete hJetShape[k][j][i];
        if(hJetShapeCase1[k][j][i])       delete hJetShapeCase1[k][j][i];
        if(hJetShapeCase2[k][j][i])       delete hJetShapeCase2[k][j][i];
        if(hJetShapeBG[k][j][i])          delete hJetShapeBG[k][j][i];
        if(hJetShapeBGCase1[k][j][i])     delete hJetShapeBGCase1[k][j][i];
        if(hJetShapeBGCase2[k][j][i])     delete hJetShapeBGCase2[k][j][i];
        if(hJetShapeBGCase3[k][j][i])     delete hJetShapeBGCase3[k][j][i];
        if(hJetCounter[k][j][i])          delete hJetCounter[k][j][i];
        if(hJetCounterCase1[k][j][i])     delete hJetCounterCase1[k][j][i];
        if(hJetCounterCase2[k][j][i])     delete hJetCounterCase2[k][j][i];
        if(hJetCounterCase3BG[k][j][i])   delete hJetCounterCase3BG[k][j][i];

        if(hJetPtProfile[k][j][i])        delete hJetPtProfile[k][j][i];
        if(hJetPtProfileCase1[k][j][i])   delete hJetPtProfileCase1[k][j][i];
        if(hJetPtProfileCase2[k][j][i])   delete hJetPtProfileCase2[k][j][i];
        if(hJetPtProfileBG[k][j][i])      delete hJetPtProfileBG[k][j][i];
        if(hJetPtProfileBGCase1[k][j][i]) delete hJetPtProfileBGCase1[k][j][i];
        if(hJetPtProfileBGCase2[k][j][i]) delete hJetPtProfileBGCase2[k][j][i];
        if(hJetPtProfileBGCase3[k][j][i]) delete hJetPtProfileBGCase3[k][j][i];
      }
    }
  }
}
//
//_________________________________________________________________________________________
Int_t StJetShapeAnalysis::Init() {
  //StJetFrameworkPicoBase::Init();

  // initialize the histograms
  DeclareHistograms();

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it
  //fJets->SetName(fJetsName);
  //fJets->SetOwner(kTRUE);

  return kStOK;
}
//
//______________________________________________________________________________________________
Int_t StJetShapeAnalysis::Finish() { 
  //  Summarize the run.
  cout << "StJetShapeAnalysis::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(fAnalysisMakerName);
    fout->cd(fAnalysisMakerName);
    WriteHistograms();
   
    // jet shape analysis
    if(doJetShapeAnalysis) {
      fout->cd();
      //fout->mkdir(Form("JetShapeAnalysis_bin%i", fTPCptAssocBin));
      //fout->cd(Form("JetShapeAnalysis_bin%i", fTPCptAssocBin));
      //fout->mkdir(Form("JetShapeAnalysis%i_bin%i", fJetAnalysisJetType, fJetShapePtAssocBin));
      //fout->cd(Form("JetShapeAnalysis%i_bin%i", fJetAnalysisJetType, fJetShapePtAssocBin));
      fout->mkdir(Form("JetShapeAnalysis_bin%i", fJetShapePtAssocBin));
      fout->cd(Form("JetShapeAnalysis_bin%i", fJetShapePtAssocBin));
      WriteJetShapeHistograms();
    }

    fout->cd();
    fout->Write();
    fout->Close();
  }

  //  Write QA histos to file and close it.
  if(mOutNameQA!="" && fJetShapePtAssocBin < 5) {
    TFile *fQAout = new TFile(mOutNameQA.Data(), "UPDATE");
    fQAout->cd();

    // track QA
    if(doWriteTrackQAHist && (fJetShapePtAssocBin < 5)) {
      fQAout->mkdir(Form("TrackQA_bin%i", fTPCptAssocBin));
      fQAout->cd(Form("TrackQA_bin%i", fTPCptAssocBin));
      WriteTrackQAHistograms();
      fQAout->cd();
    }

    fQAout->Write();
    fQAout->Close();
  }

  cout<<"End of StJetShapeAnalysis::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}
//
// Function: declare histograms and output objects
//________________________________________________________________________________________________
void StJetShapeAnalysis::DeclareHistograms() {
  // constants
  double pi = 1.0*TMath::Pi();

  int nHistCentBins;
  if(fCentBinSizeJS == 10) nHistCentBins = 10;
  if(fCentBinSizeJS ==  5) nHistCentBins = 20;

  hEventPlane = new TH1F("hEventPlane", "Event plane distribution", 72, 0.0, 1.0*pi);
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 100, -50, 50);
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", 160, 0, 800);
  hRhovsCent = new TH2F("hRhovsCent", "#rho vs centrality", 20, 0, 100, 200, 0, 200);

  // track phi distribution for centrality
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i] = new TH1F(Form("hTrackPhi%d", i), Form("track distribution vs #phi, centr%d", i), 144, 0, 2*pi);
    hTrackEta[i] = new TH1F(Form("hTrackEta%d", i), Form("track distribution vs #eta, centr%d", i), 40, -1.0, 1.0);
    hTrackPt[i] = new TH1F(Form("hTrackPt%d", i), Form("track distribution vs p_{T}, centr%d", i), 150, 0., 30.0);
  }
  hTrackEtavsPhi = new TH2F(Form("hTrackEtavsPhi"), Form("track distribution: #eta vs #phi"), 144, 0, 2*pi, 40, -1.0, 1.0);

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
  hJetPhivsEP = new TH2F("hJetPhivsEP", "Jet #phi vs event plane", 72, 0.0, 2*pi, 72, 0.0, 1.0*pi);

  // Event Selection QA histo
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
  hMBvsMult = new TH1F("hMBvsMult", "# MB events vs multiplicity", 350, 0, 700);
  hMB5vsMult = new TH1F("hMB5vsMult", "# MB5 events vs multiplicity", 350, 0, 700);
  hMB30vsMult = new TH1F("hMB30vsMult", "# MB30 events vs multiplicity", 350, 0, 700);
  hHTvsMult = new TH1F("hHTvsMult", "# HT events vs multiplicity", 350, 0, 700);
  hNMixEvents = new TH1F("hNMixEvents", "number of mixing events", 200, 0, 200);

  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        hJetShape[k][j][i] = new TH1F(Form("hJetShape_%i_%i_%i", k, j, i), Form("Jet shape #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetShapeCase1[k][j][i] = new TH1F(Form("hJetShapeCase1_%i_%i_%i", k, j, i), Form("Jet shape case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetShapeCase2[k][j][i] = new TH1F(Form("hJetShapeCase2_%i_%i_%i", k, j, i), Form("Jet shape case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetShapeBG[k][j][i] = new TH1F(Form("hJetShapeBG_%i_%i_%i", k, j, i), Form("Jet shape BG #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetShapeBGCase1[k][j][i] = new TH1F(Form("hJetShapeBGCase1_%i_%i_%i", k, j, i), Form("Jet shape BG case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetShapeBGCase2[k][j][i] = new TH1F(Form("hJetShapeBGCase2_%i_%i_%i", k, j, i), Form("Jet shape BG case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetShapeBGCase3[k][j][i] = new TH1F(Form("hJetShapeBGCase3_%i_%i_%i", k, j, i), Form("Jet shape BG case3 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetCounter[k][j][i] = new TH1F(Form("hJetCounter_%i_%i_%i", k, j, i), Form("Jet counter - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);
        hJetCounterCase1[k][j][i] = new TH1F(Form("hJetCounterCase1_%i_%i_%i", k, j, i), Form("Jet counter case2 - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);
        hJetCounterCase2[k][j][i] = new TH1F(Form("hJetCounterCase2_%i_%i_%i", k, j, i), Form("Jet counter case2 - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);
        hJetCounterCase3BG[k][j][i] = new TH1F(Form("hJetCounterCase3BG_%i_%i_%i", k, j, i), Form("Jet counter case3 BG only - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 1, 0.0, 1.0);

        hJetPtProfile[k][j][i] = new TH1F(Form("hJetPtProfile_%i_%i_%i", k, j, i), Form("Jet pt profile #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetPtProfileCase1[k][j][i] = new TH1F(Form("hJetPtProfileCase1_%i_%i_%i", k, j, i), Form("Jet pt profile case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetPtProfileCase2[k][j][i] = new TH1F(Form("hJetPtProfileCase2_%i_%i_%i", k, j, i), Form("Jet pt profile case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetPtProfileBG[k][j][i] = new TH1F(Form("hJetPtProfileBG_%i_%i_%i", k, j, i), Form("Jet pt profile BG #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetPtProfileBGCase1[k][j][i] = new TH1F(Form("hJetPtProfileBGCase1_%i_%i_%i", k, j, i), Form("Jet profile BG case1 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetPtProfileBGCase2[k][j][i] = new TH1F(Form("hJetPtProfileBGCase2_%i_%i_%i", k, j, i), Form("Jet profile BG case2 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
        hJetPtProfileBGCase3[k][j][i] = new TH1F(Form("hJetPtProfileBGCase3_%i_%i_%i", k, j, i), Form("Jet profile BG case3 #rho(r) - p_{T} bin %i, centrality bin %i, EP bin %i", k, j, i), 10, 0.0, 0.50);
      }
    }
  }

  // set up centrality bins for mixed events
  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it
  //int nCentralityBinspp = 8;
  //double centralityBinspp[9] = {0.0, 4., 9, 15, 25, 35, 55, 100.0, 500.0};  

  // FIXME
  // this is temp as the above and various other implementation attempts would not work for both cases
  Int_t nCentBins = 8;
  Double_t cenBins[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  Double_t *centralityBin = cenBins;

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

  // centrality bins for mixed events
  //Int_t nCentBinsJS = 8;
  //Double_t cenBinsJS[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  //Double_t *centralityBinJS = cenBinsJS;
  Int_t nMultBinsJS = 29;
  Double_t multBinsJS[] = {10, 14, 19, 25, 31, 37, 44, 52, 61, 71, 82, 95, 109, 124, 140, 157, 175, 194, 214, 235, 257, 280, 304, 329, 355, 382, 410, 439, 469};
  Double_t *multiplicityBinsJS = multBinsJS;

  // z-vertex bins for mixed events
  Int_t nZvBinsJS = 20;
  Double_t vBinsJS[] = {-40,-36,-32,-28,-24,-20,-16,-12,-8,-4,0,4,8,12,16,20,24,28,32,36,40};
  Double_t *zvbinJS = vBinsJS;
  //Int_t nbinsjetMIX = sizeof(vBinsJS)/sizeof(Double_t) - 1;

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implementation of AliEventPoolManager
  //fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinspp, centralityBinspp, nZvtxBins, zvtxbin);
  if(fDoUseMultBins) {
    fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nMultBinsJS, (Double_t*)multiplicityBinsJS, nZvBinsJS, (Double_t*)zvbinJS);
  } else {
    fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentralityBinsJS, (Double_t*)centralityBinsJS, nZvBinsJS, (Double_t*)zvbinJS);
  }
  //fPoolMgr = new StEventPoolManager(poolsize, trackDepth, nCentBins, (Double_t*)centralityBin, nZvBins, (Double_t*)zvbin);

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}
//
// write track QA histograms
//_____________________________________________________________________________
void StJetShapeAnalysis::WriteTrackQAHistograms() {
  // track phi distribution for centrality
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i]->Write();
    hTrackEta[i]->Write();
    hTrackPt[i]->Write();
  }
  hTrackEtavsPhi->Write();
}
//
// write jet shape histograms
//________________________________________________________________________________
void StJetShapeAnalysis::WriteJetShapeHistograms() {
  // jet shape histos
  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        hJetShape[k][j][i]->Write();
        hJetShapeCase1[k][j][i]->Write();
        hJetShapeCase2[k][j][i]->Write();
        hJetShapeBG[k][j][i]->Write();
        hJetShapeBGCase1[k][j][i]->Write();
        hJetShapeBGCase2[k][j][i]->Write();
        hJetShapeBGCase3[k][j][i]->Write();
        hJetCounter[k][j][i]->Write();
        hJetCounterCase1[k][j][i]->Write();
        hJetCounterCase2[k][j][i]->Write();
        hJetCounterCase3BG[k][j][i]->Write();

        hJetPtProfile[k][j][i]->Write();
        hJetPtProfileCase1[k][j][i]->Write();
        hJetPtProfileCase2[k][j][i]->Write();
        hJetPtProfileBG[k][j][i]->Write();
        hJetPtProfileBGCase1[k][j][i]->Write();
        hJetPtProfileBGCase2[k][j][i]->Write();
        hJetPtProfileBGCase3[k][j][i]->Write();
      }
    }
  }

}
//
// write histograms
//_____________________________________________________________________________
void StJetShapeAnalysis::WriteHistograms() {
  // default histos
  hEventPlane->Write();
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
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StJetShapeAnalysis::Clear(Option_t *opt) {
  fJets->Clear();
}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StJetShapeAnalysis::Make() {
  // zero out these global variables
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
  }

  // cut event on max track pt > 30.0 GeV
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

  // let me know the fill, and event ID
  int fillId = mPicoEvent->fillId();
  int eventId = mPicoEvent->eventId();
  if(fDebugLevel == kDebugGeneralEvt) cout<<"RunID = "<<fRunNumber<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

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
  if(fCentBinSizeJS==10) { centBinToUse = (double)centbin10 * 10.0;
  } else if(fCentBinSizeJS==5) { centBinToUse = (double)centbin * 5.0; }

  // centrality / multiplicity histograms
  hMultiplicity->Fill(refCorr2);
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
  if(fHaveEmcTrigger && (fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200 || fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200_MB || fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200)) {
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
  fRhoVal = fRho->GetVal();
  hRhovsCent->Fill(centbin*5.0, fRhoVal);
  if(fDebugLevel == kDebugRhoEstimate) cout<<"   fRhoVal = "<<fRhoVal<<"   Correction = "<<1.0*TMath::Pi()*fJetRad*fJetRad*fRhoVal<<endl;

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

  // basic method - not used for anything..
  double rpAngle = GetReactionPlane();
  hEventPlane->Fill(rpAngle);

  // ============== EventPlaneMaker =============== //
  // get StEventPlaneMaker from event
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
      EventPlaneMaker = static_cast<StEventPlaneMaker*>(GetMaker(Form("%s%i", fEventPlaneMakerNameChTemp, fTPCptAssocBin)));

      // check for requested EventPlaneMaker pointer
      if(!EventPlaneMaker) {LOG_WARN<<Form("No EventPlaneMaker bin: %i!", fTPCptAssocBin)<<endm; return kStWarn; }

      // get event plane angle for different pt bins
      double tpc2EP = (EventPlaneMaker) ? (double)EventPlaneMaker->GetTPCEP() : -999;

      // assign global event plane to selected pt-dependent bin
      TPC_PSI2 = tpc2EP;
    }
  }
  // ===================================================================================

  // ========================== Jet Shape Analysis ===================================== //
  if(doJetShapeAnalysis) {
    StEventPool *pool = 0x0;

    // require event mixing
    if(fDoEventMixing > 0) {
      // convert back to integer bins for mixed event pool - 10% bins (0, 7), 5% bins (0, 15)
      Int_t mixcentbin = TMath::Floor(fCentralityScaled / fCentBinSizeJS);
      //cout<<"fCentralityScaled: "<<fCentralityScaled<<"  fCentBinSizeJS: "<<fCentBinSizeJS<<"  mixcentbin: "<<mixcentbin<<"  zVtx: "<<zVtx<<endl;

      // initialize event pools
      if(fDoUseMultBins) { pool = fPoolMgr->GetEventPool(refCorr2, zVtx);
      } else { pool = fPoolMgr->GetEventPool(mixcentbin, zVtx); } 
      if(!pool) {
        // FIXME update this error line, to display multiplicity bin
        Form("No pool found for centrality = %i, zVtx = %f", mixcentbin, zVtx); // FIXME if cent changes to double
        return kTRUE;
      }
    }

    // check for back-to-back jets
    double ljpttemp = 0., sljpttemp = 0.;
    bool isBackToBack = kFALSE;
    if(fLeadingJet && fSubLeadingJet) {
      double BackToBackPhi = fLeadingJet->Phi() - fSubLeadingJet->Phi() - pi;
      if(BackToBackPhi < 0) BackToBackPhi += 2*pi;
      if(BackToBackPhi < 0.4) isBackToBack = kTRUE;

      // get lj and sublj pt
      if(fCorrJetPt) {
        ljpttemp = fLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
        sljpttemp = fSubLeadingJet->Pt() - fLeadingJet->Area()*fRhoVal;
      } else {
        ljpttemp = fLeadingJet->Pt();
        sljpttemp = fSubLeadingJet->Pt();
      }
    }

    // Aj selection for Jet Shape Analysis
    bool doAjSelection = (isBackToBack && ljpttemp > 20. && sljpttemp > 10.) ? kTRUE : kFALSE;
    if(doAjSelection) hJetLeadingPtAj->Fill(ljpttemp);      // leading jets
    if(doAjSelection) hJetSubLeadingPtAj->Fill(sljpttemp);  // sub-leading jets

    // Triggered events and leading/subleading jets - do Jet Shape Analysis
    // check for back to back jets: must have leading + subleading jet, subleading jet must be > 10 GeV, subleading jet must be within 0.4 of pi opposite of leading jet
    if(doRequireAjSelection) {
      if(doAjSelection && fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet) JetShapeAnalysis(fLeadingJet, pool, refCorr2);
    } else {
      if(fHaveEmcTrigger && fJetAnalysisJetType == kLeadingJets && fLeadingJet) JetShapeAnalysis(fLeadingJet, pool, refCorr2);
    }
    // subleading jets
    if(fHaveEmcTrigger && fJetAnalysisJetType == kSubLeadingJets  && fSubLeadingJet) JetShapeAnalysis(fSubLeadingJet, pool, refCorr2);

    // use only tracks from MB events
    //if(fDoEventMixing > 0 && fRunForMB && (!fHaveEmcTrigger)) { // kMB5 or kMB30 - AuAu, kMB - pp (excluding HT)
    if(fDoEventMixing > 0 && fRunForMB) { // kMB5 or kMB30 - AuAu, kMB - pp (don't exclude HT)
      // update pool: create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      pool->UpdatePool(CloneAndReduceTrackList());
      hMBvsMult->Fill(refCorr2);                       // MB5 || MB30
      if(fHaveMB5event)  hMB5vsMult->Fill(refCorr2);   // MB5
      if(fHaveMB30event) hMB30vsMult->Fill(refCorr2);  // MB30

      // fill QA histo's
      hMixEvtStatZVtx->Fill(zVtx);
      hMixEvtStatCent->Fill(centBinToUse);
      hMixEvtStatZvsCent->Fill(centBinToUse, zVtx);
    } // MB 

    // only return if doing jet shape analysis
    return kStOK;
  }
  // =================================================================================================================

  // run Track QA and fill histograms
  if((doWriteTrackQAHist) && (doJetAnalysis)) TrackQA();
  return kStOK;

  // get number of jets, tracks, and global tracks in events
  Int_t njets = fJets->GetEntries();

  // ====================== Jet loop below ============================
  for(int ijet = 0; ijet < njets; ijet++) {  // JET LOOP
    // Run - Trigger Selection to process jets from
    if(!doJetAnalysis) continue;

    // get jet pointer
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

    // FIXME, double check, might want to code this nicer
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, TPC_PSI2) : -999; // CORRECTED event plane angle - STEP3

    // this is a double check - should not happen, but if it does -> kill the job so it can be fixed - most likely in the runPico macro
    if(dEP < -900) return kStFatal;

    // some threshold cuts
    if(fCorrJetPt) {  // background subtracted jet pt
      if(corrjetpt < fMinPtJet) continue;
    } else { if(jetpt < fMinPtJet) continue; }
    if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) continue;

    // TODO check that jet contains a tower that fired the trigger
    if(!DidTowerConstituentFireTrigger(jet)) { continue; }

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
  } // jet loop

  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQAafterCuts);

  return kStOK;
}
//
//
//_________________________________________________________________________
TH1 *StJetShapeAnalysis::FillEmcTriggersHist(TH1 *h) {
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

  return h;
}
//
// From CF event mixing code PhiCorrelations
//_________________________________________________
TClonesArray *StJetShapeAnalysis::CloneAndReduceTrackList()
{
  // clones a track list by using StPicoTrack which uses much less memory (used for event mixing)
//  TClonesArray *tracksClone = new TClonesArray("StPicoTrack");// original way
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

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    // this is TEMP, it will filter track by the pt bin used for analysis - double check syntax is set up correct 
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
// basic function to get event plane angle
//________________________________________________________________________
Double_t StJetShapeAnalysis::GetReactionPlane() { 
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
    // may change this back in future
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
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 2.0 GeV (up to 5.0 GeV)
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // check for leading jet removal
    if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit) )) continue;

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
void StJetShapeAnalysis::SetSumw2() {
  // set sum weights
  hEventPlane->Sumw2();
  hEventZVertex->Sumw2();
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  hRhovsCent->Sumw2();
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

  fHistEventSelectionQA->Sumw2();
  fHistEventSelectionQAafterCuts->Sumw2();
  hTriggerIds->Sumw2();
  hEmcTriggers->Sumw2();
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

  for(int k=0; k<4; k++) {
    for(int j=0; j<4; j++) {
      for(int i=0; i<4; i++) {
        hJetShape[k][j][i]->Sumw2();
        hJetShapeCase1[k][j][i]->Sumw2();
        hJetShapeCase2[k][j][i]->Sumw2();
        hJetShapeBG[k][j][i]->Sumw2();
        hJetShapeBGCase1[k][j][i]->Sumw2();
        hJetShapeBGCase2[k][j][i]->Sumw2();
        hJetShapeBGCase3[k][j][i]->Sumw2();
        hJetCounter[k][j][i]->Sumw2();
        hJetCounterCase1[k][j][i]->Sumw2();
        hJetCounterCase2[k][j][i]->Sumw2();
        hJetCounterCase3BG[k][j][i]->Sumw2();

        hJetPtProfile[k][j][i]->Sumw2();
        hJetPtProfileCase1[k][j][i]->Sumw2();
        hJetPtProfileCase2[k][j][i]->Sumw2();
        hJetPtProfileBG[k][j][i]->Sumw2();
        hJetPtProfileBGCase1[k][j][i]->Sumw2();
        hJetPtProfileBGCase2[k][j][i]->Sumw2();
        hJetPtProfileBGCase3[k][j][i]->Sumw2();
      }
    }
  }

}
//
// track QA function to fill some histograms with track information
//
//_____________________________________________________________________________
void StJetShapeAnalysis::TrackQA()
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

    // fill other track QA plots
    hTrackPhi[ref9]->Fill(phi);
    hTrackEta[ref9]->Fill(eta);
    hTrackPt[ref9]->Fill(pt);
    hTrackEtavsPhi->Fill(phi, eta);

  }
}
//
//_________________________________________________________________________
void StJetShapeAnalysis::FillTowerTriggersArr() {
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
    int emcTrigIDindex = emcTrigIDindex - 1;

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigIDindex] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigIDindex] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigIDindex] = kTRUE;

    //cout<<"i = "<<i<<"  EmcTrigID = "<<emcTrigID<<"  adc = "<<emcTrig->adc()<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3<<endl;
  }

}
//
//
// function to require that a jet constituent tower fired a HT trigger
//___________________________________________________________________________________________
Bool_t StJetShapeAnalysis::DidTowerConstituentFireTrigger(StJet *jet) {  
  // tower constituent fired trigger
  Bool_t mFiredTrigger = kFALSE;

  // loop over constituent towers
  for(int itow = 0; itow < jet->GetNumberOfClusters(); itow++) {
    int ArrayIndex = jet->ClusterAt(itow);

    // get jet tower pointer
    StPicoBTowHit *tow = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(ArrayIndex));
    if(!tow){ continue; }
    
    // tower ID: get from index of array shifted by +1
    int towID = ArrayIndex + 1;
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
// function that does jet shape analysis
//___________________________________________________________________________________________
void StJetShapeAnalysis::JetShapeAnalysis(StJet *jet, StEventPool *pool, Double_t refCorr2) {
    // constants
    double pi = 1.0*TMath::Pi();
    double rbinSize = 0.05;

    // get centrality bin
    int centbin = Get4CentBin(fCentralityScaled);
    if(centbin < 0) return;

    // get jet info
    double jetPhi = jet->Phi();
    double jetEta = jet->Eta();
    
    // get jet pt bin and value
    double jetPt = -99.;
    if(fCorrJetPt) { jetPt = jet->Pt() - jet->Area()*fRhoVal;
    } else { jetPt = jet->Pt(); }
    int jetPtBin = GetJetPtBin(jetPt);
    if(jetPtBin < 0) return;

    // require tower and or track bias for jet
    //if((jet->GetMaxTrackPt() < fTrackBias) && (jet->GetMaxTowerEt() < fTowerBias)) return kStOK;

    // check that jet contains a tower that fired the trigger
    //if(!DidTowerConstituentFireTrigger(jet)) return kStOK;

    // get event plane bin - FIXME, double check, might want to code this nicer
    double dEP = (!doppAnalysis) ? RelativeEPJET(jetPhi, TPC_PSI2) : -99.; // CORRECTED event plane angle - STEP3
    if(doTPCptassocBin && !doppAnalysis) {
      // z = if(condition) then(?) <do this> else(:) <do this>  
      double dEP0 = RelativeEPJET(jetPhi, TPC_PSI2);
      double dEP1 = RelativeEPJET(jetPhi, TPC_PSI2);
      double dEP2 = RelativeEPJET(jetPhi, TPC_PSI2);
      double dEP3 = RelativeEPJET(jetPhi, TPC_PSI2);
      double dEP4 = RelativeEPJET(jetPhi, TPC_PSI2);
      if(fJetShapePtAssocBin == 0) dEP = dEP0; // 0.2-0.5 GeV
      if(fJetShapePtAssocBin == 1) dEP = dEP1; // 0.5-1.0 GeV
      if(fJetShapePtAssocBin == 2) dEP = dEP2; // 1.0-1.5 GeV
      if(fJetShapePtAssocBin == 3) dEP = dEP3; // 1.5-2.0 GeV
      if(fJetShapePtAssocBin == 4) dEP = dEP4; // 2.0-3.0 GeV
      if(fJetShapePtAssocBin == 5) dEP = dEP4; // 3.0-4.0 GeV
      if(fJetShapePtAssocBin == 6) dEP = dEP4; // 4.0-6.0 GeV
      if(fJetShapePtAssocBin == 7) dEP = dEP4; // 6.0+    GeV
      if(fJetShapePtAssocBin == 8) dEP = dEP4; // 0.5+    GeV (inclusive)
    }
    int EPBin = GetJetEPBin(dEP);
    if(EPBin < 0) return;

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

      // primary track switch: get momentum vector of track - global or primary track
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
      if(tphi < 0.0   ) tphi += 2.0*pi; // require 0,2pi interval
      if(tphi > 2.0*pi) tphi -= 2.0*pi; 
      double teta = mTrkMom.PseudoRapidity();

      // cut on track pt range 
      if(tpt < fJetShapeTrackPtMin) { continue; }
      if(tpt > fJetShapeTrackPtMax) { continue; }

      // additional pt selection when doing pt associated bin method
      if(doTPCptassocBin) {
        if(fJetShapePtAssocBin == 0) { if((tpt < 0.20) || (tpt  >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 1) { if((tpt < 0.50) || (tpt  >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 2) { if((tpt < 1.00) || (tpt  >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 3) { if((tpt < 1.50) || (tpt  >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 4) { if((tpt < 2.00) || (tpt  >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 5) { if((tpt < 3.00) || (tpt  >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 6) { if((tpt < 4.00) || (tpt  >= 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations, 8->6 Oct23
        if(fJetShapePtAssocBin == 7) { if((tpt <  6.0))                 continue; }  //       6.0+ GeV assoc bin used for correlations, 8->6 Oct23
        if(fJetShapePtAssocBin == 8) { if((tpt <  0.5))                 continue; }  //       0.5+ GeV assoc bin used for correlations
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
    hHTvsMult->Fill(refCorr2);
    hTriggerEvtStatZVtx->Fill(zVtx);
    hTriggerEvtStatCent->Fill(fCentralityScaled);
    hTriggerEvtStatZvsCent->Fill(fCentralityScaled, zVtx);

    // fill jet shape histograms
    for(int i = 0; i < 10; i++) { 
      hJetShape[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt); 
      hJetShape[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt); 

      hJetPtProfile[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
      hJetPtProfile[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
    }
    hJetCounter[jetPtBin][centbin][EPBin]->Fill(0.5);
    hJetCounter[jetPtBin][centbin][3]->Fill(0.5); // ALL angles

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
        hJetShapeCase1[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt); 
        hJetShapeCase1[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);    

        hJetPtProfileCase1[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
        hJetPtProfileCase1[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
      }
      hJetCounterCase1[jetPtBin][centbin][EPBin]->Fill(0.5);
      hJetCounterCase1[jetPtBin][centbin][3]->Fill(0.5); // ALL angles 

      // background
      case1 = kTRUE;
      jetEtaBG = - jetEta;
      jetPhiBG =   jetPhi;
    }

    // CASE 2: Phi shifted
    if(TMath::Abs(jetEta) < etaMin) {
      for(int i = 0; i < 10; i++) { 
        hJetShapeCase2[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);
        hJetShapeCase2[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]/jetPt);    

        hJetPtProfileCase2[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
        hJetPtProfileCase2[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsum[i]);
      }
      hJetCounterCase2[jetPtBin][centbin][EPBin]->Fill(0.5);
      hJetCounterCase2[jetPtBin][centbin][3]->Fill(0.5); // ALL angles

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

    // track loop inside jet loop - loop over ALL tracks in PicoDst - for BG
    for(int itrack = 0; itrack < ntracks; itrack++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(itrack));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

      // primary track switch: get momentum vector of track - global or primary track
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
      if(tphi < 0.0) tphi += 2.0*pi; // require 0,2pi interval
      double teta = mTrkMom.PseudoRapidity();

      // cut on track pt
      if(tpt < fJetShapeTrackPtMin) { continue; }
      if(tpt > fJetShapeTrackPtMax) { continue; }

      // additional pt selection when doing pt associated bin method
      if(doTPCptassocBin) {
        if(fJetShapePtAssocBin == 0) { if((tpt < 0.20) || (tpt  >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 1) { if((tpt < 0.50) || (tpt  >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 2) { if((tpt < 1.00) || (tpt  >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 3) { if((tpt < 1.50) || (tpt  >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 4) { if((tpt < 2.00) || (tpt  >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 5) { if((tpt < 3.00) || (tpt  >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
        if(fJetShapePtAssocBin == 6) { if((tpt < 4.00) || (tpt  >= 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations, 8->6 Oct23
        if(fJetShapePtAssocBin == 7) { if((tpt <  6.0))                 continue; }  //       6.0+ GeV assoc bin used for correlations, 8->6 Oct23
        if(fJetShapePtAssocBin == 8) { if((tpt <  0.5))                 continue; }  //       0.5+ GeV assoc bin used for correlations
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
      hJetShapeBG[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);
      hJetShapeBG[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);    

      hJetPtProfileBG[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
      hJetPtProfileBG[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
    }

    // Case 1: background
    if(case1) { 
      for(int i = 0; i < 10; i++) {
        hJetShapeBGCase1[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt); 
        hJetShapeBGCase1[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);    

        hJetPtProfileBGCase1[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
        hJetPtProfileBGCase1[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
      }
    }

    // Case 2: background
    if(case2) {
      for(int i = 0; i < 10; i++) {
        hJetShapeBGCase2[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);
        hJetShapeBGCase2[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]/jetPt);    

        hJetPtProfileBGCase2[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
        hJetPtProfileBGCase2[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG[i]);
      }
    }

/*
    // Case 3: background
    if(case3) {
      for(int i = 0; i < 10; i++) {
        hJetShapeBGCase3[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]/jetPt);
        hJetShapeBGCase3[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]/jetPt);
        hJetPtProfileBGCase3[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]);
        hJetPtProfileBGCase3[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i]);
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
        hNMixEvents->Fill(nMix);

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
            // get Femto track pointer
            StFemtoTrack *trk = static_cast<StFemtoTrack*>(bgTracks->At(ibg));
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
              if(fJetShapePtAssocBin == 0) { if((Mpt < 0.20) || (Mpt  >= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
              if(fJetShapePtAssocBin == 1) { if((Mpt < 0.50) || (Mpt  >= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
              if(fJetShapePtAssocBin == 2) { if((Mpt < 1.00) || (Mpt  >= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
              if(fJetShapePtAssocBin == 3) { if((Mpt < 1.50) || (Mpt  >= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
              if(fJetShapePtAssocBin == 4) { if((Mpt < 2.00) || (Mpt  >= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
              if(fJetShapePtAssocBin == 5) { if((Mpt < 3.00) || (Mpt  >= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
              if(fJetShapePtAssocBin == 6) { if((Mpt < 4.00) || (Mpt  >= 6.0)) continue; }  // 4.00 - 6.0 GeV assoc bin used for correlations, 8->6 Oct23
              if(fJetShapePtAssocBin == 7) { if((Mpt <  6.0))                 continue; }  //       6.0+ GeV assoc bin used for correlations, 8->6 Oct23
              if(fJetShapePtAssocBin == 8) { if((Mpt <  0.5))                 continue; }  //       0.5+ GeV assoc bin used for correlations
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
            //int effCent   = mCentMaker->GetRef16();
            //double fZDCx  = mPicoEvent->ZDCx();
            //double mixEfficiency = ApplyTrackingEff(fDoEffCorr, Mpt, Meta, effCent, fZDCx, fTrackEfficiencyType);

          } // end of background track loop

          // fill BG histos here
          for(int i = 0; i < 10; i++) {
            hJetShapeBGCase3[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
            hJetShapeBGCase3[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix*jetPt));
            hJetPtProfileBGCase3[jetPtBin][centbin][EPBin]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
            hJetPtProfileBGCase3[jetPtBin][centbin][3]->Fill(i*rbinSize + 1e-3, 1.0*rsumBG3[i] / (nMix));
          } // loop over annuli bins
          //hJetCounterCase3BG[jetPtBin][centbin][EPBin]->Fill(0.5);
          //hJetCounterCase3BG[jetPtBin][centbin][3]->Fill(0.5); // ALL angles

        }   // end of filling mixed-event histo's:  jth mix event loop

        // fill background counters for Case 3
        hJetCounterCase3BG[jetPtBin][centbin][EPBin]->Fill(0.5);
        hJetCounterCase3BG[jetPtBin][centbin][3]->Fill(0.5); // ALL angles
      }     // end of check for pool being ready
    }       // end of event mixing

    return;
}
