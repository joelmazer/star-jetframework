// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################
// $Id$

#include "StPicoTrackClusterQA.h"
#include "StMemStat.h"
#include <sstream>
#include <fstream>

// root classes
#include <TChain.h>
#include <TClonesArray.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "TFile.h"
#include <THnSparse.h>
#include "TVector3.h"

// StRoot classes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StPicoConstants.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu

// test for clusters:
#include "StMuDSTMaker/COMMON/StMuTrack.h"
// StEmc:
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
///#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StEmcClusterCollection.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
//#include "StEmcUtil/others/emcDetectorName.h"
//#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcRawHit.h"
#include "StEmcModule.h"
#include "StEmcDetector.h"

// extra includes
#include "StEmcPosition2.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StJetFrameworkPicoBase.h"

// centrality
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// towers
#include "StJetPicoTower.h"
#include "StJetPicoDefinitions.h"

class StJetFrameworkPicoBase;

ClassImp(StPicoTrackClusterQA)
//
//________________________________________________________________________
StPicoTrackClusterQA::StPicoTrackClusterQA() : 
//  StJetFrameworkPicoBase(),
  StMaker(),
  doWriteHistos(kFALSE),
  doUsePrimTracks(kFALSE), 
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE),
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum
  fDoEffCorr(kFALSE),
  fDoTowerQAforHT(kFALSE),
  fMaxEventTrackPt(30.0),
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fRequireCentSelection(kFALSE),
  doUseBBCCoincidenceRate(kFALSE), // kFALSE = use ZDC
  mOutName(""),
  fAnalysisMakerName(""),
  fTracksName(""),
  fCaloName(""),
  fTrackPtMinCut(0.2),
  fTrackPtMaxCut(30.0),
  fClusterPtMinCut(0.2),
  fClusterPtMaxCut(100.0),
  fTrackPhiMinCut(0.0),
  fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0),
  fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15),
  fTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fGoodTrackCounter(0),
  fTowerEMinCut(0.2),
  fTowerEMaxCut(100.0),
  fTowerEtaMinCut(-1.0),
  fTowerEtaMaxCut(1.0),
  fTowerPhiMinCut(0.0),
  fTowerPhiMaxCut(2.0*TMath::Pi()),
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  //mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fEmcTriggerEventType(0),
  fMBEventType(2),  // kVPDMB5
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
  mBemcTables(0x0),
  mBemcMatchedTracks(),
  mTowerStatusMode(AcceptAllTowers),
  mTowerEnergyMin(0.2),
  mHadronicCorrFrac(1.0),
  mMuDstMaker(0x0),
  mMuDst(0x0),
  mMuInputEvent(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0),
  fhnTrackQA(0x0),
  fhnTowerQA(0x0)
{
  // Default constructor.
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4801; i++) { 
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE; 
  }
}
//
//________________________________________________________________________
StPicoTrackClusterQA::StPicoTrackClusterQA(const char *name, bool doHistos = kFALSE, const char* outName = "") : 
//  StJetFrameworkPicoBase(name),
  StMaker(name),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE),
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum
  fDoEffCorr(kFALSE),
  fDoTowerQAforHT(kFALSE),
  fMaxEventTrackPt(30.0),
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fRequireCentSelection(kFALSE),
  doUseBBCCoincidenceRate(kFALSE), // kFALSE = use ZDC
  mOutName(outName),
  fAnalysisMakerName(name),
  fTracksName("Tracks"),
  fCaloName("Clusters"),
  fTrackPtMinCut(0.2), //0.20
  fTrackPtMaxCut(30.0), 
  fClusterPtMinCut(0.2),
  fClusterPtMaxCut(100.0),
  fTrackPhiMinCut(0.0),
  fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0), 
  fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15),
  fTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fGoodTrackCounter(0),
  fTowerEMinCut(0.2),
  fTowerEMaxCut(100.0),
  fTowerEtaMinCut(-1.0),
  fTowerEtaMaxCut(1.0),
  fTowerPhiMinCut(0.0),
  fTowerPhiMaxCut(2.0*TMath::Pi()),
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  //mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fEmcTriggerEventType(0),
  fMBEventType(2),   // kVPDMB5
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
  mBemcTables(0x0),
  mBemcMatchedTracks(),
  mTowerStatusMode(AcceptAllTowers),
  mTowerEnergyMin(0.2),
  mHadronicCorrFrac(1.0),
  mMuDstMaker(0x0),
  mMuDst(0x0),
  mMuInputEvent(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0),
  fhnTrackQA(0x0),
  fhnTowerQA(0x0)
{
  // Standard constructor.
  if (!name) return;
  SetName(name);

  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  for(int i=0; i<4801; i++) { 
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }
}
//
//________________________________________________________________________
StPicoTrackClusterQA::~StPicoTrackClusterQA()
{
  // free up histogram objects if they exist

  // Destructor
  if(fHistNTrackvsPt)   delete fHistNTrackvsPt;
  if(fHistNTrackvsPhi)  delete fHistNTrackvsPhi;
  if(fHistNTrackvsEta)  delete fHistNTrackvsEta;
  if(fHistNTrackvsPhivsEta) delete fHistNTrackvsPhivsEta;
  if(fHistNHadCorrTowervsE)    delete fHistNHadCorrTowervsE;
  if(fHistNHadCorrTowervsEt)   delete fHistNHadCorrTowervsEt;
  if(fHistNHadCorrTowervsPhi)  delete fHistNHadCorrTowervsPhi;
  if(fHistNHadCorrTowervsEta)  delete fHistNHadCorrTowervsEta;
  if(fHistNHadCorrTowervsPhivsEta) delete fHistNHadCorrTowervsPhivsEta;
  if(fHistNHadCorrTowerHOTvsTowID) delete fHistNHadCorrTowerHOTvsTowID;
  if(fHistNTowervsE)    delete fHistNTowervsE;
  if(fHistNTowervsEt)   delete fHistNTowervsEt;
  if(fHistNTowervsPhi)  delete fHistNTowervsPhi;
  if(fHistNTowervsEta)  delete fHistNTowervsEta;
  if(fHistNTowervsPhivsEta) delete fHistNTowervsPhivsEta;
  if(fHistNTowerHOTvsTowID) delete fHistNTowerHOTvsTowID;

  delete fHistEventSelectionQA;
  delete fHistEventSelectionQAafterCuts;
  delete hTriggerIds;
  delete hEmcTriggers;
  delete fHistTriggerIDs;

  if(fProfEventRefMult) delete fProfEventRefMult;
  if(fProfEventRanking) delete fProfEventRanking;
  if(fProfEventZvtx) delete fProfEventZvtx;
  if(fProfEventYvtx) delete fProfEventYvtx;
  if(fProfEventXvtx) delete fProfEventXvtx;
  if(fProfEventVzVPD)delete fProfEventVzVPD;
  if(fProfEventBBCx) delete fProfEventBBCx;
  if(fProfEventZDCx) delete fProfEventZDCx;

  if(fHistNZeroEHT1vsID) delete fHistNZeroEHT1vsID;
  if(fHistNZeroEHT2vsID) delete fHistNZeroEHT2vsID;
  if(fHistNZeroEHT3vsID) delete fHistNZeroEHT3vsID;
  if(fHistNNegEHT1vsID)  delete fHistNNegEHT1vsID;
  if(fHistNNegEHT2vsID)  delete fHistNNegEHT2vsID;
  if(fHistNNegEHT3vsID)  delete fHistNNegEHT3vsID;
  if(fHistNFiredHT0vsID) delete fHistNFiredHT0vsID;
  if(fHistNFiredHT1vsID) delete fHistNFiredHT1vsID;
  if(fHistNFiredHT2vsID) delete fHistNFiredHT2vsID;
  if(fHistNFiredHT3vsID) delete fHistNFiredHT3vsID;
  if(fHistHT0FiredEtvsID) delete fHistHT0FiredEtvsID;
  if(fHistHT1FiredEtvsID) delete fHistHT1FiredEtvsID;
  if(fHistHT2FiredEtvsID) delete fHistHT2FiredEtvsID;
  if(fHistHT3FiredEtvsID) delete fHistHT3FiredEtvsID;
  if(fHistHT0IDvsFiredEt) delete fHistHT0IDvsFiredEt;
  if(fHistHT1IDvsFiredEt) delete fHistHT1IDvsFiredEt;
  if(fHistHT2IDvsFiredEt) delete fHistHT2IDvsFiredEt;
  if(fHistHT3IDvsFiredEt) delete fHistHT3IDvsFiredEt;

  if(fHistNFiredHT0vsFlag) delete fHistNFiredHT0vsFlag;
  if(fHistNFiredHT1vsFlag) delete fHistNFiredHT1vsFlag;
  if(fHistNFiredHT2vsFlag) delete fHistNFiredHT2vsFlag;
  if(fHistNFiredHT3vsFlag) delete fHistNFiredHT3vsFlag;
  if(fHistNFiredJP0vsFlag) delete fHistNFiredJP0vsFlag;
  if(fHistNFiredJP1vsFlag) delete fHistNFiredJP1vsFlag;
  if(fHistNFiredJP2vsFlag) delete fHistNFiredJP2vsFlag;

  if(fHistNFiredHT0vsADC) delete fHistNFiredHT0vsADC;
  if(fHistNFiredHT1vsADC) delete fHistNFiredHT1vsADC;
  if(fHistNFiredHT2vsADC) delete fHistNFiredHT2vsADC;
  if(fHistNFiredHT3vsADC) delete fHistNFiredHT3vsADC;
  if(fHistNFiredJP0vsADC) delete fHistNFiredJP0vsADC;
  if(fHistNFiredJP1vsADC) delete fHistNFiredJP1vsADC;
  if(fHistNFiredJP2vsADC) delete fHistNFiredJP2vsADC;

  delete fhnTrackQA;
  delete fhnTowerQA;
  delete mEmcPosition;
}
//
//_____________________________________________________________________________
Int_t StPicoTrackClusterQA::Init() {
  DeclareHistograms();

  // test placement
  mBemcTables = new StBemcTables();

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  //AddBadTowers( TString( getenv("STARPICOPATH" )) + "/badTowerList_y11.txt");
  // Add dead + bad tower lists
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp 200 GeV
        AddBadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_BadTowers.txt");
        AddDeadTowers("StRoot/StMyAnalysisMaker/towerLists/Empty_DeadTowers.txt");
        break;

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

  // Centrality object setup
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
        break;

    default :
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  }

  return kStOK;
}
//
//_____________________________________________________________________________
Int_t StPicoTrackClusterQA::Finish() {
  //  Summarize the run.
  cout << "StPicoTrackClusterQA::Finish()\n";

  // open output file
  if(doWriteHistos && mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE"); //"RECREATE");
    fout->cd();
    // GetName()
    fout->mkdir(fAnalysisMakerName);
    fout->cd(fAnalysisMakerName);

    cout<<fAnalysisMakerName<<endl;

    // write histograms and output file before closing
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StPicoTrackClusterQA::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}
//
//________________________________________________________________________
void StPicoTrackClusterQA::DeclareHistograms() {
    // declare histograms
    double pi = 1.0*TMath::Pi();

    // set binning for run based corrections - run dependent
    Int_t nRunBins = 1; // - just a default
    if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   nRunBins = 857 + 43;
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) nRunBins = 830; //1654;
    if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) nRunBins = 1359;
    Double_t nRunBinsMax = (Double_t)nRunBins + 0.5;

    // track histograms
    fHistNTrackvsPt = new TH1F("fHistNTrackvsPt", "Ntracks vs p_{T}", 150, 0., 30.);
    fHistNTrackvsPhi = new TH1F("fHistNTrackvsPhi", "Ntracks vs #phi", 144, 0., 2.0*pi);
    fHistNTrackvsEta = new TH1F("fHistNTrackvsEta", "Ntracks vs #eta", 40, -1.0, 1.0);
    fHistNTrackvsPhivsEta = new TH2F("fHistNTrackvsPhivsEta", "Ntrack vs #phi vs #eta", 144, 0, 2.0*pi, 40, -1.0, 1.0);

    // Tower histograms
    fHistNHadCorrTowervsE = new TH1F("fHistNHadCorrTowervsE", "NHadCorrTowers vs energy", 100, 0., 20.0);
    fHistNHadCorrTowervsEt = new TH1F("fHistNHadCorrTowervsEt", "NHadCorrTowers vs transverse energy", 100, 0., 20.0);
    fHistNHadCorrTowervsPhi = new TH1F("fHistNHadCorrTowervsPhi", "NHadCorrTowers vs #phi", 144, 0., 2.0*pi);
    fHistNHadCorrTowervsEta = new TH1F("fHistNHadCorrTowervsEta", "NHadCorrTowers vs #eta", 40, -1.0, 1.0);
    fHistNHadCorrTowervsPhivsEta = new TH2F("fHistNHadCorrTowervsPhivsEta", "NHadCorrTowers vs #phi vs #eta", 144, 0, 2.0*pi, 40, -1.0, 1.0);
    fHistNHadCorrTowerHOTvsTowID = new TH1F("fHistNHadCorrTowerHOTvsTowID", "NHadCorrTowers HOT vs tower ID", 4800, 0.5, 4800.5);
    fHistNTowervsE = new TH1F("fHistNTowervsE", "Ntowers vs energy", 100, 0., 20.0);
    fHistNTowervsEt = new TH1F("fHistNTowervsEt", "Ntowers vs transverse energy", 100, 0., 20.0);
    fHistNTowervsPhi = new TH1F("fHistNTowervsPhi", "Ntowers vs #phi", 144, 0., 2.0*pi);
    fHistNTowervsEta = new TH1F("fHistNTowervsEta", "Ntowers vs #eta", 40, -1.0, 1.0);
    fHistNTowervsPhivsEta = new TH2F("fHistNTowervsPhivsEta", "Ntowers vs #phi vs #eta", 144, 0, 2.0*pi, 40, -1.0, 1.0);
    fHistNTowerHOTvsTowID = new TH1F("fHistNTowerHOTvsTowID", "NTowerHOT vs tower ID", 4800, 0.5, 4800.5);

    // Event Selection QA histograms
    fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
    fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
    hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
    hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);
    fHistTriggerIDs = new TH1F("fHistTriggerIDs", "NTriggers vs trigger IDs", 30, 0.5, 30.5);

    // event QA
    fProfEventRefMult = new TProfile("fProfEventRefMult", "Event averaged refMult", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventRanking = new TProfile("fProfEventRanking", "Event averaged vertex ranking", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventZvtx = new TProfile("fProfEventZvtx", "Event averaged primary z-Vertex", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventYvtx = new TProfile("fProfEventYvtx", "Event averaged primary y-Vertex", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventXvtx = new TProfile("fProfEventXvtx", "Event averaged primary x-Vertex", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventVzVPD = new TProfile("fProfEventVzVPD", "Event averaged VzVPD", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventBBCx = new TProfile("fProfEventBBCx", "Event averaged BBC coincidence rate", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
    fProfEventZDCx = new TProfile("fProfEventZDCx", "Event averaged ZDC coincidence rate", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);

    // trigger histograms, zero and negative entries QA
    fHistNZeroEHT1vsID = new TH1F("fHistNZeroEHT1vsID", "NTowers fired HT1 with zero E vs tower ID", 4800, 0.5, 4800.5);
    fHistNZeroEHT2vsID = new TH1F("fHistNZeroEHT2vsID", "NTowers fired HT2 with zero E vs tower ID", 4800, 0.5, 4800.5);
    fHistNZeroEHT3vsID = new TH1F("fHistNZeroEHT3vsID", "NTowers fired HT3 with zero E vs tower ID", 4800, 0.5, 4800.5);
    fHistNNegEHT1vsID = new TH1F("fHistNNegEHT1vsID", "NTowers fired HT1 with negative E vs tower ID", 4800, 0.5, 4800.5);
    fHistNNegEHT2vsID = new TH1F("fHistNNegEHT2vsID", "NTowers fired HT2 with negative E vs tower ID", 4800, 0.5, 4800.5);
    fHistNNegEHT3vsID = new TH1F("fHistNNegEHT3vsID", "NTowers fired HT3 with negative E vs tower ID", 4800, 0.5, 4800.5);

    // trigger histograms: firing towers QA
    fHistNFiredHT0vsID = new TH1F("fHistNFiredHT0vsID", "NTrig fired HT0 vs tower ID", 4800, 0.5, 4800.5);
    fHistNFiredHT1vsID = new TH1F("fHistNFiredHT1vsID", "NTrig fired HT1 vs tower ID", 4800, 0.5, 4800.5);
    fHistNFiredHT2vsID = new TH1F("fHistNFiredHT2vsID", "NTrig fired HT2 vs tower ID", 4800, 0.5, 4800.5);
    fHistNFiredHT3vsID = new TH1F("fHistNFiredHT3vsID", "NTrig fired HT3 vs tower ID", 4800, 0.5, 4800.5);
    fHistHT0FiredEtvsID = new TH1F("fHistHT0FiredEtvsID", "HT0 fired transverse energy vs tower ID", 4800, 0.5, 4800.5);
    fHistHT1FiredEtvsID = new TH1F("fHistHT1FiredEtvsID", "HT1 fired transverse energy vs tower ID", 4800, 0.5, 4800.5);
    fHistHT2FiredEtvsID = new TH1F("fHistHT2FiredEtvsID", "HT2 fired transverse energy vs tower ID", 4800, 0.5, 4800.5);
    fHistHT3FiredEtvsID = new TH1F("fHistHT3FiredEtvsID", "HT3 fired transverse energy vs tower ID", 4800, 0.5, 4800.5);
    fHistHT0IDvsFiredEt = new TH2F("fHistHT0IDvsFiredEt", "HT0 tower ID vs fired transverse energy", 200, 0.0, 20.0, 4800, 0.5, 4800.5);
    fHistHT1IDvsFiredEt = new TH2F("fHistHT1IDvsFiredEt", "HT1 tower ID vs fired transverse energy", 200, 0.0, 20.0, 4800, 0.5, 4800.5);
    fHistHT2IDvsFiredEt = new TH2F("fHistHT2IDvsFiredEt", "HT2 tower ID vs fired transverse energy", 200, 0.0, 20.0, 4800, 0.5, 4800.5);
    fHistHT3IDvsFiredEt = new TH2F("fHistHT3IDvsFiredEt", "HT3 tower ID vs fired transverse energy", 200, 0.0, 20.0, 4800, 0.5, 4800.5);

    fHistNFiredHT0vsFlag = new TH1F("fHistNFiredHT0vsFlag", "NTowers fired HT0 vs Flag", 125, -0.5, 124.5);
    fHistNFiredHT1vsFlag = new TH1F("fHistNFiredHT1vsFlag", "NTowers fired HT1 vs Flag", 125, -0.5, 124.5);
    fHistNFiredHT2vsFlag = new TH1F("fHistNFiredHT2vsFlag", "NTowers fired HT2 vs Flag", 125, -0.5, 124.5);
    fHistNFiredHT3vsFlag = new TH1F("fHistNFiredHT3vsFlag", "NTowers fired HT3 vs Flag", 125, -0.5, 124.5);
    fHistNFiredJP0vsFlag = new TH1F("fHistNFiredJP0vsFlag", "NTowers fired JP0 vs Flag", 125, -0.5, 124.5);
    fHistNFiredJP1vsFlag = new TH1F("fHistNFiredJP1vsFlag", "NTowers fired JP1 vs Flag", 125, -0.5, 124.5);
    fHistNFiredJP2vsFlag = new TH1F("fHistNFiredJP2vsFlag", "NTowers fired JP2 vs Flag", 125, -0.5, 124.5);

    fHistNFiredHT0vsADC = new TH1F("fHistNFiredHT0vsADC", "NTowers fired HT0 vs ADc", 100, 0., 100.);
    fHistNFiredHT1vsADC = new TH1F("fHistNFiredHT1vsADC", "NTowers fired HT1 vs ADC", 100, 0., 100.);
    fHistNFiredHT2vsADC = new TH1F("fHistNFiredHT2vsADC", "NTowers fired HT2 vs ADC", 100, 0., 100.);
    fHistNFiredHT3vsADC = new TH1F("fHistNFiredHT3vsADC", "NTowers fired HT3 vs ADC", 100, 0., 100.);
    fHistNFiredJP0vsADC = new TH1F("fHistNFiredJP0vsADC", "NTowers fired JP0 vs ADC", 100, 0., 100.);
    fHistNFiredJP1vsADC = new TH1F("fHistNFiredJP1vsADC", "NTowers fired JP1 vs ADC", 100, 0., 100.);
    fHistNFiredJP2vsADC = new TH1F("fHistNFiredJP2vsADC", "NTowers fired JP2 vs ADC", 100, 0., 100.);

    // set up track and tower sparse
    UInt_t bitcodeTrack = 0; // bit coded, see GetDimParams() below
    bitcodeTrack = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4; 
    fhnTrackQA = NewTHnSparseFTracks("fhnTrackQA", bitcodeTrack);

    // set up track and tower sparse
    UInt_t bitcodeTower = 0; // bit coded, see GetDimParams() below
    bitcodeTower = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4;                      
    fhnTowerQA = NewTHnSparseFTowers("fhnTowerQA", bitcodeTower);

    // Switch on Sumw2 for all histos - (except profiles)
    SetSumw2();
}
//
// write histograms
//________________________________________________________________________
void StPicoTrackClusterQA::WriteHistograms() {
  // track and tower histograms
  fHistNTrackvsPt->Write();
  fHistNTrackvsPhi->Write();
  fHistNTrackvsEta->Write();
  fHistNTrackvsPhivsEta->Write();
  fHistNHadCorrTowervsE->Write();
  fHistNHadCorrTowervsEt->Write();
  fHistNHadCorrTowervsPhi->Write();
  fHistNHadCorrTowervsEta->Write();
  fHistNHadCorrTowervsPhivsEta->Write();
  fHistNHadCorrTowerHOTvsTowID->Write();
  fHistNTowervsE->Write();
  fHistNTowervsEt->Write();
  fHistNTowervsPhi->Write();
  fHistNTowervsEta->Write();
  fHistNTowervsPhivsEta->Write();
  fHistNTowerHOTvsTowID->Write();

  // QA histograms
  fHistEventSelectionQA->Write();
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
  fHistTriggerIDs->Write();
  fProfEventRefMult->Write();
  fProfEventRanking->Write();
  fProfEventZvtx->Write();
  fProfEventYvtx->Write();
  fProfEventXvtx->Write();
  fProfEventVzVPD->Write();
  fProfEventBBCx->Write();
  fProfEventZDCx->Write();

  // trigger QA histograms
  fHistNZeroEHT1vsID->Write();
  fHistNZeroEHT2vsID->Write();
  fHistNZeroEHT3vsID->Write();
  fHistNNegEHT1vsID->Write();
  fHistNNegEHT2vsID->Write();
  fHistNNegEHT3vsID->Write();
  fHistNFiredHT0vsID->Write();
  fHistNFiredHT1vsID->Write();
  fHistNFiredHT2vsID->Write();
  fHistNFiredHT3vsID->Write();
  fHistHT0FiredEtvsID->Write();
  fHistHT1FiredEtvsID->Write();
  fHistHT2FiredEtvsID->Write();
  fHistHT3FiredEtvsID->Write();
  fHistHT0IDvsFiredEt->Write();
  fHistHT1IDvsFiredEt->Write();
  fHistHT2IDvsFiredEt->Write();
  fHistHT3IDvsFiredEt->Write();

  fHistNFiredHT0vsFlag->Write();
  fHistNFiredHT1vsFlag->Write();
  fHistNFiredHT2vsFlag->Write();
  fHistNFiredHT3vsFlag->Write();
  fHistNFiredJP0vsFlag->Write();
  fHistNFiredJP1vsFlag->Write();
  fHistNFiredJP2vsFlag->Write();

  fHistNFiredHT0vsADC->Write();
  fHistNFiredHT1vsADC->Write();
  fHistNFiredHT2vsADC->Write();
  fHistNFiredHT3vsADC->Write();
  fHistNFiredJP0vsADC->Write();
  fHistNFiredJP1vsADC->Write();
  fHistNFiredJP2vsADC->Write();

  // sparses
  //fhnTrackQA->Write();
  //fhnTowerQA->Write();
}
//
//
//________________________________________________________________________
void StPicoTrackClusterQA::Clear(Option_t *opt) {
}
//
// Main loop, called for each event.
//________________________________________________________________________
int StPicoTrackClusterQA::Make()
{
  fGoodTrackCounter = 0;

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

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  // get vertex 3-vector and declare variables
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  // Z-vertex cut - per the Aj analysis (-40, 40)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // ============================ CENTRALITY ============================== //
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  int RunId = mPicoEvent->runId();
  float fBBCCoincidenceRate = mPicoEvent->BBCx();
  float fZDCCoincidenceRate = mPicoEvent->ZDCx();
  int grefMult = mPicoEvent->grefMult();
  Int_t centbin, cent16;

  // for non-pp analyses - calculate centrality
  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate); }
    else{ grefmultCorr->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ grefmultCorr->getRefMultCorr(grefMult, zVtx, fZDCCoincidenceRate, 2); }

    cent16 = grefmultCorr->getCentralityBin16();
    centbin = GetCentBin(cent16, 16);
  } else {
    centbin = 0, cent16 = 0;
  }

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStWarn; // maybe kStOk; - this is for lowest multiplicity events 80%+ centrality, cut on them
  fCentralityScaled = centbin*5.0;

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }
  // ===================================================

  // ========================= Trigger Info =============================== //
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  FillEmcTriggersHist(hEmcTriggers);

  // fill trigger IDs
  FillTriggerIDs(fHistTriggerIDs);

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  //if(fDebugLevel == kDebugEmcTrigger) 
  //cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    //if(fDebugLevel == kDebugEmcTrigger) 
    //cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  }
  //if(fDebugLevel == kDebugEmcTrigger) 
  //cout<<endl;
  // ======================== end of Triggers ============================= //

  // Event QA print-out
  // printing available information from the PicoDst objects
/*
  StPicoTrack* t = mPicoDst->track(1);
  if(t) t->Print();
  StPicoEmcTrigger *emcTrig = mPicoDst->emcTrigger(0);
  if(emcTrig) emcTrig->Print();
  StPicoMtdTrigger *mtdTrig = mPicoDst->mtdTrigger(0);
  if(mtdTrig) mtdTrig->Print();
  StPicoBTowHit *btowhit = mPicoDst->btowHit(0); 
  if(btowhit) btowhit->Print();
  StPicoBTofHit *btofhit = mPicoDst->btofHit(0);
  if(btofhit) btofhit->Print();
  //StPicoMtdHit *mtdhit = mPicoDst->mtdHit(0);
  //mtdhit->Print();
  StPicoBEmcPidTraits *emcpid = mPicoDst->bemcPidTraits(0); // OLD NAME (StPicoEmcPidTraits, emcPidTraits) now its StPicoBEmcPidTraits
  if(emcpid) emcpid->Print();
  StPicoBTofPidTraits *tofpid = mPicoDst->btofPidTraits(0);
  if(tofpid) tofpid->Print();
  StPicoMtdPidTraits *mtdpid = mPicoDst->mtdPidTraits(0);
  if(mtdpid) mtdpid->Print();
*/

  // switches for QA analysis
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);
  bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu
  if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE;

  // run tower QA for specific conditions
  if(fDoTowerQAforHT && fHaveEmcTrigger)  {
    RunEventQA();
    RunFiredTriggerQA();  //cout<<"HT.."<<endl; } // HT trigger
    RunTrackQA(); 
    RunTowerQA();
    RunHadCorrTowerQA();
  }

  // look for firing trigger specifically requested
  //if(!fDoTowerQAforHT && fHaveMBevent) {
  if(!fDoTowerQAforHT && fRunForMB) { // updated MB type Dec4, 2018
    RunEventQA();
    RunFiredTriggerQA();
    RunTrackQA();
    RunTowerQA(); 
    RunHadCorrTowerQA();
  }

/*
  int nTracks = mPicoDst->numberOfTracks();
  int nTrigs = mPicoDst->numberOfEmcTriggers();
  int nBTowHits = mPicoDst->numberOfBTowHits();
  int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();
  cout<<"nTracks = "<<nTracks<<"  nTrigs = "<<nTrigs<<"  nBTowHits = "<<nBTowHits<<"  nBEmcPidTraits = "<<nBEmcPidTraits<<endl;
  cout<<"highTowerThreshold 0 = "<<mPicoEvent->highTowerThreshold(0)<<endl;
  cout<<"highTowerThreshold 1 = "<<mPicoEvent->highTowerThreshold(1)<<endl;
  cout<<"highTowerThreshold 2 = "<<mPicoEvent->highTowerThreshold(2)<<endl;
  cout<<"highTowerThreshold 3 = "<<mPicoEvent->highTowerThreshold(3)<<endl;
  cout<<"jetPatchThreshold 0 = "<<mPicoEvent->jetPatchThreshold(0)<<endl;
  cout<<"jetPatchThreshold 1 = "<<mPicoEvent->jetPatchThreshold(1)<<endl;
  cout<<"jetPatchThreshold 2 = "<<mPicoEvent->jetPatchThreshold(2)<<endl;
*/

  // Event / object PRINT INFO!!
  //mPicoDst->printTracks();
  //mPicoDst->printTriggers();
  ////mPicoDst->printBTowHits();
  ////mPicoDst->printBEmcPidTraits();

  // fill Event QA after cuts
  //FillEventTriggerQA(fHistEventSelectionQAafterCuts);

  return kStOK;
}
//
// 
//________________________________________________________________________
void StPicoTrackClusterQA::RunTrackQA()
{
  // assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // track variables
  unsigned int ntracks = mPicoDst->numberOfTracks();

  // loop over ALL tracks in PicoDst 
  for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
    StPicoTrack* trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
    if(!trk){ continue; }

    // TODO - double check this, was commenting out for test before (Feb21, 2018)
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
    double p = mTrkMom.Mag();
    //double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
    short charge = trk->charge();         
    int bemcIndex = trk->bemcPidTraitsIndex();

    // shift track phi (0, 2*pi)
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    if(fDebugLevel == 8) cout<<"iTracks = "<<iTracks<<"  p = "<<p<<"  charge = "<<charge<<"  eta = "<<eta<<"  phi = "<<phi;
    if(fDebugLevel == 8) cout<<"  nHitsFit = "<<trk->nHitsFit()<<"  BEmc Index = "<<bemcIndex<<endl;

    // fill some QA histograms
    fHistNTrackvsPt->Fill(pt);
    fHistNTrackvsPhi->Fill(phi);
    fHistNTrackvsEta->Fill(eta);
    fHistNTrackvsPhivsEta->Fill(phi, eta);

    // fill track sparse
    Double_t trackEntries[5] = {fCentralityScaled, pt, eta, phi, zVtx};
    Double_t trefficiency = 1.0; // TODO update
    fhnTrackQA->Fill(trackEntries, 1.0/trefficiency);

    fGoodTrackCounter++;
  } // track loop

  // looping over clusters - STAR: matching already done
  // get # of clusters and set position variables
  unsigned int nclus = mPicoDst->numberOfBEmcPidTraits();
  TVector3  towPosition, clusPosition;

  // print EMCal cluster info
  if(fDebugLevel == 7) mPicoDst->printBEmcPidTraits();

  // loop over ALL clusters in PicoDst and add to jet //TODO
  for(unsigned short iClus = 0; iClus < nclus; iClus++){
    StPicoBEmcPidTraits* cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
    if(!cluster){ continue; }

    // print index of associated track in the event (debug = 2)
    if(fDebugLevel == 8) cout<<"iClus = "<<iClus<<"  trackIndex = "<<cluster->trackIndex()<<"  nclus = "<<nclus<<endl;

    // cluster and tower ID: ID's are calculated as such:
    // mBtowId       = (ntow[0] <= 0 || ntow[0] > 4800) ? -1 : (Short_t)ntow[0];
    // mBtowId23 = (ntow[1] < 0 || ntow[1] >= 9 || ntow[2] < 0 || ntow[2] >= 9) ? -1 : (Char_t)(ntow[1] * 10 + ntow[2]);
    int clusID = cluster->bemcId();  // index in bemc point array
    int towID = cluster->btowId();   // projected tower Id: 1 - 4800
    int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

    // cluster and tower position - from vertex and ID
    towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
    clusPosition = mEmcPosition->getPosFromVertex(mVertex, clusID); // INNACCURATE for BEmcPidTraits objects

    // index of associated track in the event
    int trackIndex = cluster->trackIndex();

    // get track corresponding to EMC pidTraits in question from its index
    StPicoTrack* trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
    //StMuTrack* trkMu = (StMuTrack*)mPicoDst->track(trackIndex);
    if(!trk) continue;
    //if(doUsePrimTracks) { if(!(trk->isPrimary())) continue; } // check if primary 

    if(fDebugLevel == 8) cout<<"cluster bemcId = "<<cluster->bemcId()<<"  cluster btowId = "<<cluster->btowId();
    if(fDebugLevel == 8) cout<<"  cluster trackIndex = "<<cluster->trackIndex()<<"  trkId = "<<trk->id()<<endl;

    // April 9, 2018: removed ton of useless code that was here just to originally figure out tower - track matching and indexing
    // add some useful QA code back here soon!

    // fill tower sparse
    //Double_t towerEntries[5] = {fCentralityScaled, cluster->bemcE0(), towEta, towPhi, zVtx};
    //fhnTowerQA->Fill(towerEntries);

  } // cluster loop

}  // track/cluster QA
//
// Set sum weights for histogram uncertainties
//________________________________________________________________________________
void StPicoTrackClusterQA::SetSumw2() {
  fHistNTrackvsPt->Sumw2();
  fHistNTrackvsPhi->Sumw2();
  fHistNTrackvsEta->Sumw2();
  fHistNTrackvsPhivsEta->Sumw2();
  fHistNHadCorrTowervsE->Sumw2();
  fHistNHadCorrTowervsEt->Sumw2();
  fHistNHadCorrTowervsPhi->Sumw2();
  fHistNHadCorrTowervsEta->Sumw2();
  fHistNHadCorrTowervsPhivsEta->Sumw2();
  fHistNHadCorrTowerHOTvsTowID->Sumw2();
  fHistNTowervsE->Sumw2();
  fHistNTowervsEt->Sumw2();
  fHistNTowervsPhi->Sumw2();
  fHistNTowervsEta->Sumw2();
  fHistNTowervsPhivsEta->Sumw2();
  fHistNTowerHOTvsTowID->Sumw2();

  // QA histograms
  fHistEventSelectionQA->Sumw2();
  fHistEventSelectionQAafterCuts->Sumw2();
  hTriggerIds->Sumw2();
  hEmcTriggers->Sumw2();
  fHistTriggerIDs->Sumw2();
  //fProfEventRefMult-Sumw2();
  //fProfEventRanking->Sumw2();
  //fProfEventZvtx->Sumw2();
  //fProfEventYvtx->Sumw2();
  //fProfEventXvtx->Sumw2();
  //fProfEventVzVPD->Sumw2();
  //fProfEventBBCx->Sumw2();
  //fProfEventZDCx->Sumw2();

  // trigger QA histograms
  fHistNZeroEHT1vsID->Sumw2();
  fHistNZeroEHT2vsID->Sumw2();
  fHistNZeroEHT3vsID->Sumw2();
  fHistNNegEHT1vsID->Sumw2();
  fHistNNegEHT2vsID->Sumw2();
  fHistNNegEHT3vsID->Sumw2();
  fHistNFiredHT0vsID->Sumw2();
  fHistNFiredHT1vsID->Sumw2();
  fHistNFiredHT2vsID->Sumw2();
  fHistNFiredHT3vsID->Sumw2();
  fHistHT0FiredEtvsID->Sumw2();
  fHistHT1FiredEtvsID->Sumw2();
  fHistHT2FiredEtvsID->Sumw2();
  fHistHT3FiredEtvsID->Sumw2();
  fHistHT0IDvsFiredEt->Sumw2();
  fHistHT1IDvsFiredEt->Sumw2();
  fHistHT2IDvsFiredEt->Sumw2();
  fHistHT3IDvsFiredEt->Sumw2();

  fHistNFiredHT0vsFlag->Sumw2();
  fHistNFiredHT1vsFlag->Sumw2();
  fHistNFiredHT2vsFlag->Sumw2();
  fHistNFiredHT3vsFlag->Sumw2();
  fHistNFiredJP0vsFlag->Sumw2();
  fHistNFiredJP1vsFlag->Sumw2();
  fHistNFiredJP2vsFlag->Sumw2();

  fHistNFiredHT0vsADC->Sumw2();
  fHistNFiredHT1vsADC->Sumw2();
  fHistNFiredHT2vsADC->Sumw2();
  fHistNFiredHT3vsADC->Sumw2();
  fHistNFiredJP0vsADC->Sumw2();
  fHistNFiredJP1vsADC->Sumw2();
  fHistNFiredJP2vsADC->Sumw2();

  fhnTrackQA->Sumw2();
  fhnTowerQA->Sumw2();
}
//
//
//________________________________________________________________________
Int_t StPicoTrackClusterQA::GetCentBin(Int_t cent, Int_t nBin) const
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
//
//
//__________________________________________________________________________________________
Bool_t StPicoTrackClusterQA::SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut) {
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

  bool doAnalysis;
  doAnalysis = kFALSE; // set false by default, to make sure user chooses an available bin

  // switch on bin selection
  switch(fCentralitySelectionCut) {
    case StJetFrameworkPicoBase::kCent010 :  // 0-10%
      if((centbin>-1) && (centbin<2)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent020 :  // 0-20%
      if((centbin>-1) && (centbin<4)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent1020 : // 10-20%
      if((centbin>1) && (centbin<4)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent1030 : // 10-30%
      if((centbin>1) && (centbin<6)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent1040 : // 10-40%
      if((centbin>1) && (centbin<8)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent2030 : // 20-30%
      if((centbin>3) && (centbin<6)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent2040 : // 20-40%
      if((centbin>3) && (centbin<8)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent2050 : // 20-50%
      if((centbin>3) && (centbin<10)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent2060 : // 20-60%
      if((centbin>3) && (centbin<12)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent3050 : // 30-50%
      if((centbin>5) && (centbin<10)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent3060 : // 30-60%
      if((centbin>5) && (centbin<12)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent4060 : // 40-60%
      if((centbin>7) && (centbin<12)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent4070 : // 40-70%
      if((centbin>7) && (centbin<14)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent4080 : // 40-80%
      if((centbin>7) && (centbin<16)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent5080 : // 50-80%
      if((centbin>9) && (centbin<16)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    case StJetFrameworkPicoBase::kCent6080 : // 60-80%
      if((centbin>11) && (centbin<16)) { doAnalysis = kTRUE; }
      else { doAnalysis = kFALSE; }
      break;

    default : // wrong entry
      doAnalysis = kFALSE;

  }

  return doAnalysis;
}
// 
// Track Quality Cuts
//________________________________________________________________________
Bool_t StPicoTrackClusterQA::AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert) {
  // constants: assume neutral pion mass
  //double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0*TMath::Pi();

  // get momentum vector of track - global or primary track
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
  //double dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // track pt, eta, phi cut
  if(pt < fTrackPtMinCut) return kFALSE;
  if(pt > fTrackPtMaxCut) return kFALSE; // 20.0 STAR -> 30.0 GeV, 100.0 ALICE
  if((eta < fTrackEtaMinCut) || (eta > fTrackEtaMaxCut)) return kFALSE;
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  if((phi < fTrackPhiMinCut) || (phi > fTrackPhiMaxCut)) return kFALSE;

  // additional quality cuts for tracks
  if(dca > fTrackDCAcut)            return kFALSE;
  if(nHitsFit < fTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
// 
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StPicoTrackClusterQA::AcceptTower(StPicoBTowHit *tower, Int_t towerID) {
  // constants:
  double pi = 1.0*TMath::Pi();

  // tower ID - passed into function
  //int towerID = tower->id();

  // make sure some of these aren't still in event array
  if(towerID < 0) return kFALSE;

  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  double phi = towerPosition.Phi();
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  double eta = towerPosition.PseudoRapidity();
  //int towerADC = tower->adc();
  //double towerEunCorr = tower->energy();  // uncorrected energy

  // check for bad (and dead) towers
  bool TowerOK = IsTowerOK(towerID);      // kTRUE means GOOD
  bool TowerDead = IsTowerDead(towerID);  // kTRUE means BAD
  if(!TowerOK)  { return kFALSE; }
  if(TowerDead) { return kFALSE; }

  // jet track acceptance cuts now - after getting 3vector - hardcoded
  if((eta < fTowerEtaMinCut) || (eta > fTowerEtaMaxCut)) return kFALSE;
  if((phi < fTowerPhiMinCut) || (phi > fTowerPhiMaxCut)) return kFALSE;

  // passed all above cuts - keep tower and fill input vector to fastjet
  return kTRUE;
}
//
// fill histogram keeping track of fired EMC triggers: HT, JP
//_________________________________________________________________________
TH1* StPicoTrackClusterQA::FillEmcTriggersHist(TH1* h) {
  // zero out trigger array and get number of Emcal Triggers
  for(int i = 0; i < 8; i++) { fEmcTriggerArr[i] = 0; }
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  //if(fDebugLevel == kDebugEmcTrigger) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  // tower - HT trigger types array
  // zero these out - so they are refreshed for each event
  for(int i = 0; i < 4801; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // emc trigger parameters
    int emcTrigID = emcTrig->id();
    //fTowerToTriggerType[i] = -1;

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT0 = emcTrig->isHT0();
    bool isHT1 = emcTrig->isHT1();
    bool isHT2 = emcTrig->isHT2();
    bool isHT3 = emcTrig->isHT3();
    if(isHT1) fTowerToTriggerTypeHT1[emcTrigID] = kTRUE;
    if(isHT2) fTowerToTriggerTypeHT2[emcTrigID] = kTRUE;
    if(isHT3) fTowerToTriggerTypeHT3[emcTrigID] = kTRUE;

    // print some EMCal Trigger info
    if(fDebugLevel == kDebugEmcTrigger) {
      cout<<"i = "<<i<<"  id = "<<emcTrigID<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
      cout<<"  isHT0: "<<isHT0<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3;
      cout<<"  isJP0: "<<emcTrig->isJP0()<<"  isJP1: "<<emcTrig->isJP1()<<"  isJP2: "<<emcTrig->isJP2()<<endl;
    }

    // fill for valid triggers
    if(isHT0) { h->Fill(1); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(isHT1) { h->Fill(2); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(isHT2) { h->Fill(3); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(isHT3) { h->Fill(4); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { h->Fill(5); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { h->Fill(6); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { h->Fill(7); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
  // kAny trigger - filled once per event
  h->Fill(10);

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
// fill histogram for various triggered events: based on RunFlag
//_____________________________________________________________________________
// Trigger QA histogram, label bins 
TH1* StPicoTrackClusterQA::FillEventTriggerQA(TH1* h) {
  // check and fill a Event Selection QA histogram for different trigger selections after cuts

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
  // Run14 triggers:
    int arrHT1[] = {450201, 450211, 460201};
    int arrHT2[] = {450202, 450212, 460202, 460212};
    int arrHT3[] = {450203, 450213, 460203};
    int arrMB[] = {450014};
    int arrMB30[] = {450010, 450020};
    int arrCentral5[] = {450010, 450020};
    int arrCentral[] = {460101, 460111};
    int arrMB5[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { bin = 4; h->Fill(bin); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); } // MB 
    if(DoComparison(arrCentral5, sizeof(arrCentral5)/sizeof(*arrCentral5))) { bin = 7; h->Fill(bin); }// Central-5
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 8; h->Fill(bin); } // Central & Central-mon
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { bin = 10; h->Fill(bin); }// VPDMB-5 
    if(DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { bin = 11; h->Fill(bin); } // VPDMB-30

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
    h->GetXaxis()->SetBinLabel(2, "BHT1*VPDMB-30");
    h->GetXaxis()->SetBinLabel(3, "BHT2*VPDMB-30");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, "Central-5");
    h->GetXaxis()->SetBinLabel(8, "Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-30");
  }

  // Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    // hard-coded trigger Ids for run16
    //int arrHT0[] = {520606, 520616, 520626, 520636, 520646, 520656};
    int arrHT1[] = {520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
    int arrHT2[] = {530202, 540203};
    int arrHT3[] = {520203, 530213};
    int arrMB[] = {520021};
    int arrMB5[] = {520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
    int arrMB10[] = {520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};
    int arrCentral[] = {520101, 520111, 520121, 520131, 520141, 520103, 520113, 520123};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { bin = 4; h->Fill(bin); } // HT3
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); }  // MB
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 7; h->Fill(bin); }// Central-5 & Central-novtx
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { bin = 10; h->Fill(bin); } // VPDMB-5 
    if(DoComparison(arrMB10, sizeof(arrMB10)/sizeof(*arrMB10))) { bin = 11; h->Fill(bin); }// VPDMB-10

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
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

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}
//
// elems: sizeof(myarr)/sizeof(*myarr) prior to passing to function
// upon passing the array collapses to a pointer and can not get size anymore
//________________________________________________________________________
Bool_t StPicoTrackClusterQA::DoComparison(int myarr[], int elems) {
  //std::cout << "Length of array = " << (sizeof(myarr)/sizeof(*myarr)) << std::endl;
  bool match = kFALSE;

  // loop over specific physics selection array and compare to specific event trigger
  for(int i=0; i<elems; i++) {
    if(mPicoEvent->isTrigger(myarr[i])) match = kTRUE;
    if(match) break;
  }
  //cout<<"elems: "<<elems<<"  match: "<<match<<endl;

  return match;
}
//
//______________________________________________________________________
THnSparse* StPicoTrackClusterQA::NewTHnSparseFTracks(const char* name, UInt_t entries) {
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
      GetDimParamsTracks(i, label, nbins[c], xmin[c], xmax[c]);
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
void StPicoTrackClusterQA::GetDimParamsTracks(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   //stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

    case 0:
      label = "centrality 5% bin";
      nbins = 20; //16;
      xmin = 0.;
      xmax = 100.; //16.;
      break;

    case 1:
      label = "track p_{T}";
      nbins = 300;
      xmin = 0.;
      xmax = 30.;
      break;

    case 2:
      label = "track #eta";
      nbins = 24;
      xmin = -1.2;
      xmax =  1.2;
      break;

    case 3:
      label = "track #phi";
      nbins = 144;
      xmin = 0.0;
      xmax = 2.0*pi;
      break;

    case 4:
      label = "Z-vertex";
      nbins = 20;
      xmin = -40.;
      xmax = 40.;
      break;

   }// end of switch
} // end of track sparse
//
//______________________________________________________________________
THnSparse* StPicoTrackClusterQA::NewTHnSparseFTowers(const char* name, UInt_t entries) {
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
      GetDimParamsTowers(i, label, nbins[c], xmin[c], xmax[c]);
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
void StPicoTrackClusterQA::GetDimParamsTowers(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   //stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

    case 0:
      label = "centrality 5% bin";
      nbins = 20; //16;
      xmin = 0.;
      xmax = 100.; //16.;
      break;

    case 1:
      label = "tower E";
      nbins = 200;
      xmin = 0.;
      xmax = 20.;
      break;

    case 2:
      label = "tower #eta";
      nbins = 24;
      xmin = -1.2;
      xmax =  1.2;
      break;

    case 3:
      label = "tower #phi";
      nbins = 144;
      xmin = 0.0;
      xmax = 2.0*pi;
      break;

    case 4:
      label = "Z-vertex";
      nbins = 20;
      xmin = -40.;
      xmax = 40.;
      break;

   }// end of switch
} // end of tower sparse
//
//============================================================================
//
// - not used, taken from Kolja's framework 
// this is for a test... for now.. Feb12, 2018
//____________________________________________________________________________
Bool_t StPicoTrackClusterQA::MuProcessBEMC() {
  /* this is reworked, but a simple reimplementation
     of the old Maker code. It should be updated for
     future use - currently a track is matched to the
     first tower within dR of its BEMC projection
     where it should be matched to the tower
     with the minimum dR. This might require
     a rework of MuProcessPrimaryTracks as well.
   */
  /* define the number of modules in the BEMC */
  UInt_t   mBEMCModules         = 120;
  UInt_t   mBEMCTowPerModule    = 20;
  Int_t    mBEMCSubSections     = 2;
  
  /* and we will count the # of towers/tracks matched, as this
     is saved in the event header. Not sure why 
   */
  Int_t nMatchedTowers = 0;
  Int_t nMatchedTracks = 0;
 
  // are collections used in AuAu??? FIXME
  // not according to: http://www.star.bnl.gov/public/comp/meet/RM200311/MuDstTutorial.pdf 
  StEmcCollection* mEmcCollection = static_cast<StEmcCollection*>(mMuDst->emcCollection());
  StMuEmcCollection* mMuEmcCollection = static_cast<StMuEmcCollection*>(mMuDst->muEmcCollection());
  mBemcTables->loadTables(static_cast<StMaker*>(this));
  
  /* if no collections are found, exit assuming error */
  if (!mEmcCollection || !mMuEmcCollection) return kFALSE;

  StEmcDetector* mBEMC = static_cast<StEmcDetector*>(mEmcCollection->detector(kBarrelEmcTowerId));
  StSPtrVecEmcPoint& mEmcContainer =  mEmcCollection->barrelPoints();
  
  /* if we can't get the detector exit assuming error */
  if (!mBEMC) return kFALSE;

  /* if there are no hits, continue assuming everything
     is working */
  if(mEmcContainer.size() == 0) return kTRUE;

  ////TStarJetPicoTower jetTower;
  
  /* loop over all modules */
  for (UInt_t i = 1; i <= mBEMCModules; ++i) {
    StSPtrVecEmcRawHit& mEmcTowerHits = mBEMC->module(i)->hits();
    
    /* loop over all hits in the module */
    for (UInt_t j = 0; j < mEmcTowerHits.size(); ++j) {
      StEmcRawHit* tow = mEmcTowerHits[j];
      
      if (abs(tow->module()) > mBEMCModules  || abs(tow->eta()) > mBEMCTowPerModule || tow->sub() >  mBEMCSubSections) continue;
      
      Int_t towerID, towerStatus;
      mGeom->getId((int)tow->module(), (int)tow->eta(), (int)tow->sub(), towerID);
      mBemcTables->getStatus(BTOW, towerID, towerStatus);
      
      if (mTowerStatusMode == RejectBadTowerStatus && towerStatus != 1) continue;
      
      Float_t towerEnergy = tow->energy();
      Float_t towerADC = tow->adc();
     
      // mTowerEnergyMin = 0.2
      if(towerEnergy < 0.2) continue;
      
      Float_t towerEta, towerPhi;
      mGeom->getEtaPhi(towerID, towerEta, towerPhi);
      
      /* check for SMD hits in the tower */
      Int_t ehits = MuFindSMDClusterHits(mEmcCollection, towerEta, towerPhi, 2);
      Int_t phits = MuFindSMDClusterHits(mEmcCollection, towerEta, towerPhi, 3);
      
      /* correct eta for Vz position */
      Float_t theta;
      theta = 2 * atan(exp(-towerEta)); /* getting theta from eta */
      Double_t z = 0;
      if(towerEta != 0) z = 231.0 / tan(theta);  /* 231 cm = radius of SMD */
      //Double_t zNominal = z - mMuDst->event()->primaryVertexPosition().z(); /* shifted z */
      Double_t zNominal = z - mVertex.z();       /* shifted z*/    // should be fixed
      Double_t thetaCorr = atan2(231.0, zNominal); /* theta with respect to primary vertex */
      Float_t etaCorr =-log(tan(thetaCorr / 2.0)); /* eta with respect to primary vertex */
      
      /* now match tracks to towers */
      for (unsigned k = 0; k < mBemcMatchedTracks.size(); ++k) {
        BemcMatch match = mBemcMatchedTracks[k];
        if (match.globalId == -1) continue;
        
        Double_t halfTowerWidth = 0.025;
        Double_t dEta = match.matchEta - towerEta;
        Double_t dPhi = match.matchPhi - towerPhi;
        while (dPhi > TMath::Pi())  dPhi -= TMath::Pi();
        while (dPhi < -TMath::Pi()) dPhi += TMath::Pi();
        if (fabs(dEta) > halfTowerWidth || fabs(dPhi) > halfTowerWidth) continue;
        nMatchedTracks++;
        mBemcMatchedTracks[k].globalId = -1;
        ////jetTower.AddMatchedTrack(match.trackId);
        
        /* set dEta & dPhi for the matched track */
        ////TStarJetPicoPrimaryTrack* matchedTrack = (TStarJetPicoPrimaryTrack*) mEvent->GetPrimaryTracks()->At(match.trackId);
        ////matchedTrack->SetEtaDiffHitProjected(dEta);
        ////matchedTrack->SetPhiDiffHitProjected(dPhi);
      }
      
      ////if (jetTower.GetNAssocTracks() > 0) nMatchedTowers++;
      ////mEvent->AddTower(&jetTower);
    }
    
  }
  ////mEvent->GetHeader()->SetNOfMatchedTracks(nMatchedTracks);
  ////mEvent->GetHeader()->SetNOfMatchedTowers(nMatchedTowers);
  
  return kTRUE;
}
//
// - Not used, from Kolja's framework
//___________________________________________________________________________________________________
Int_t StPicoTrackClusterQA::MuFindSMDClusterHits(StEmcCollection* coll, Double_t eta, Double_t phi, Int_t detectorID) {
  StEmcCluster* smdCluster = nullptr;
  Float_t dRmin = 9999;
  Float_t dRmin_cut = 0.05;
  StDetectorId id = static_cast<StDetectorId>(detectorID + kBarrelEmcTowerId);
  
  StEmcDetector* detector = coll->detector(id);
  if (!detector) return 0;
  StEmcClusterCollection* clusters = detector->cluster();
  if (!clusters) return 0;
  StSPtrVecEmcCluster& cl=clusters->clusters();
  
  for (unsigned i = 0; i < cl.size(); ++i) {
    Float_t clEta = cl[i]->eta();
    Float_t clPhi = cl[i]->phi();
    Float_t dR = sqrt((eta - clEta) * (eta - clEta) + (phi - clPhi) * (phi - clPhi));
    
    if (dR < dRmin && dR < dRmin_cut) {
      dRmin = dR;
      smdCluster = cl[i];
    }
    
  }
  
  if (smdCluster) return smdCluster->nHits();
  else return 0;
}
//
// function to get tower QA for hadronically corrected towers
//________________________________________________________________________
void StPicoTrackClusterQA::RunHadCorrTowerQA()
{
  // set / initialize some variables
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // towerStatus array
  float mTowerMatchTrkIndex[4801] = { 0 };
  bool mTowerStatusArr[4801] = { 0 };
  int matchedTowerTrackCounter = 0;

  // print
  //int nTracks = mPicoDst->numberOfTracks();
  //int nTrigs = mPicoDst->numberOfEmcTriggers();
  //int nBTowHits = mPicoDst->numberOfBTOWHits();
  int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();
  //cout<<"nTracks = "<<nTracks<<"  nTrigs = "<<nTrigs<<"  nBTowHits = "<<nBTowHits<<"  nBEmcPidTraits = "<<nBEmcPidTraits<<endl;
  
  // loop over ALL clusters in PicoDst
  for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
    StPicoBEmcPidTraits* cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
    if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

    // cluster and tower ID
    int clusID = cluster->bemcId();  // index in bemc point array
    int towID = cluster->btowId();   // projected tower Id: 1 - 4800
    int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

    // cluster and tower position - from vertex and ID
    TVector3  towPosition = mEmcPosition->getPosFromVertex(mVertex, towID);
    double towPhi = towPosition.Phi();
    if(towPhi < 0.0)    towPhi += 2.0*pi;
    if(towPhi > 2.0*pi) towPhi -= 2.0*pi;
    //double towEta = towPosition.PseudoRapidity();
    //double detectorRadius = mGeom->Radius();

    // matched track index
    int trackIndex = cluster->trackIndex();
    StPicoTrack* trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
    if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; } 

    // tower status set - towerID is matched to track passing quality cuts
    mTowerMatchTrkIndex[towID] = trackIndex;
    mTowerStatusArr[towID] = kTRUE;
    matchedTowerTrackCounter++;

    // print tower and track info
    //cout<<"towers: towPhi = "<<towPhi<<"  towEta = "<<towEta<<"  etaCorr = "<<etaCorr;  //<<endl;
    //cout<<"tracks:  pt = "<<pt<<"  p = "<<p<<"  phi = "<<phi<<"  eta = "<<eta<<endl;
  } // loop over BEmcPidTraits

  // print statment on matches
  //cout<<"Matched Tracks passing cuts (with tower): "<<matchedTowerTrackCounter<<"  nBTowHits = ";
  //cout<<mPicoDst->numberOfBTowHits()<<"  unFiltered Tracks = "<<mPicoDst->numberOfTracks()<<"  Filtered Tracks = "<<fGoodTrackCounter<<endl;

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // quality/acceptance cuts - TODO may not want to run this for QA
    if(!AcceptTower(tower, itow+1)) { continue; }

    // tower ID - get from itow shift, +1 above array index (needed since fall 2018)
    //int towerID = tower->id();
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    double towerPhi = towerPosition.Phi();
    if(towerPhi < 0.0)    towerPhi += 2.0*pi;
    if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;
    double towerEta = towerPosition.PseudoRapidity();
    //int towerADC = tower->adc();
    //double towerEunCorr = tower->energy();
    double towerE = tower->energy();

    // perform tower cuts - if tower was not matched to an accepted track, use it for jet by itself if > 0.2 GeV (constituent cut)
    if(mTowerStatusArr[towerID]) {
      //if(mTowerMatchTrkIndex[towerID] > 0) 
      StPicoTrack* trk = static_cast<StPicoTrack*>(mPicoDst->track( mTowerMatchTrkIndex[towerID] ));
      if(!trk) { cout<<"No trk pointer...."<<endl; continue; } 
      if(!AcceptTrack(trk, Bfield, mVertex)) { cout<<"track matched back doesn't pass cuts"<<endl; continue; }  // TODO double check this

      // get track variables to matched tower
      TVector3 mTrkMom;
      if(doUsePrimTracks) { 
        // get primary track vector
        mTrkMom = trk->pMom(); 
      } else { 
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield); 
      }

      //double pt = mTrkMom.Perp();
      //double phi = mTrkMom.Phi();
      //double eta = mTrkMom.PseudoRapidity();
      double p = mTrkMom.Mag();
      double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);

      // apply hadronic correction
      towerE = towerE - (mHadronicCorrFrac * E);
    } 
    // else - no match so treat towers on their own

    // cut on transvere tower energy - corrected or not
    double towerEt = towerE / (1.0*TMath::CosH(towerEta));
    if(towerEt < 0) towerEt = 0.0;
    if(towerEt < mTowerEnergyMin) continue;

    // print
    //cout<<"itow: "<<itow<<"  towerID = "<<towerID<<"  towerPhi = "<<towerPhi<<"  towerEta = "<<towerEta<<"  towerADC = "<<towerADC<<"  towerE = "<<towerE<<"  towerEunCorr = "<<towerEunCorr<<"  mIndex = "<<mTowerMatchTrkIndex[towerID]<<endl;

    // fill QA histos for towers
    fHistNHadCorrTowervsE->Fill(towerE);
    fHistNHadCorrTowervsEt->Fill(towerEt);
    fHistNHadCorrTowervsPhi->Fill(towerPhi);
    fHistNHadCorrTowervsEta->Fill(towerEta);
    fHistNHadCorrTowervsPhivsEta->Fill(towerPhi, towerEta);
    if(towerEt > 2.0) fHistNHadCorrTowerHOTvsTowID->Fill(towerID);  // HOT tower histogram

  } // tower loop

  //} // cluster loop
} // cluster / tower QA
//
// function to get tower QA
//________________________________________________________________________
void StPicoTrackClusterQA::RunTowerQA()
{
  // set / initialize some variables
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // loop over towers
  int nTowers = mPicoDst->numberOfBTowHits();
  for(int itow = 0; itow < nTowers; itow++) {
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
    if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

    // quality/acceptance cuts - TODO may not want to run this for QA
    //if(!AcceptTower(tower, itow+1)) { continue; }  // TURN this off for RAW QA of towers

    // tower ID - get from itow shift, +1 above array index (needed since fall 2018)
    int towerID = itow + 1;
    if(towerID < 0) continue; // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    double towerPhi = towerPosition.Phi();
    if(towerPhi < 0.0)    towerPhi += 2.0*pi;
    if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;
    double towerEta = towerPosition.PseudoRapidity();
    //int towerADC = tower->adc();
    double towerE = tower->energy();

    // cut on transvere tower energy - corrected or not
    double towerEt = towerE / (1.0*TMath::CosH(towerEta));
    if(towerEt < 0) towerEt = 0.0;
    if(towerEt < mTowerEnergyMin) continue;

    // fill QA histos for towers
    fHistNTowervsE->Fill(towerE);
    fHistNTowervsEt->Fill(towerEt);
    fHistNTowervsPhi->Fill(towerPhi);
    fHistNTowervsEta->Fill(towerEta);
    fHistNTowervsPhivsEta->Fill(towerPhi, towerEta);
    if(towerEt > 2.0) fHistNTowerHOTvsTowID->Fill(towerID);  // HOT tower histogram

  } // tower loop

}   //  tower QA
//
// this function is used for filling histograms used to determine bad and dead towers
//________________________________________________________________________
void StPicoTrackClusterQA::RunFiredTriggerQA()
{
  // assume neutral pion mass
  //double pi0mass = Pico::mMass[0]; // GeV
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // set flags for triggers - HAVE TO DO IT RUN-BY-RUN by of problems with STAR and how it saves them
  UInt_t HT0flag, HT1flag, HT2flag, HT3flag, JP0flag, JP1flag, JP2flag;
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    HT1flag = 10;
    HT2flag = 14;
    HT3flag = 8;  // HT3 is BAD for this dataset
    JP0flag = 16;
    JP1flag = 48;
    JP2flag = 112;
  } 
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) { 
    HT1flag = 2;
    HT2flag = 6;
    HT3flag = 14;  
    JP0flag = 112;
    JP1flag = 112;
    JP2flag = 112;
  }
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    HT1flag = 2;
    HT2flag = 6;
    HT3flag = 14;
    JP0flag = 112;
    JP1flag = 112;
    JP2flag = 112;
  }

  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // emc trigger parameters
    int emcTrigID = emcTrig->id();
    int adc = emcTrig->adc();
    UInt_t flag = emcTrig->flag();
    // flags: HT1 = 2, HT2 = 6, HT3 = 14, JP = 112
    // Run12: HT1 = 10, HT2 = 14, HT3 = 8, JP0 = 16,48,112 JP2 = 48,112, JP3 = 112

    // check if i'th trigger fired HT triggers by meeting threshold
    bool isHT0 = emcTrig->isHT0(); // ADC > 11, pp: ADC > ?
    bool isHT1 = emcTrig->isHT1(); // ADC > 15, pp: ADC > 15
    bool isHT2 = emcTrig->isHT2(); // ADC > 18, pp: ADC > 18
    bool isHT3 = emcTrig->isHT3(); // ADC > 25, pp: ADC >  8
    bool isJP0 = emcTrig->isJP0(); //           pp: ADC low 21
    bool isJP1 = emcTrig->isJP1(); //           pp: ADC low 29
    bool isJP2 = emcTrig->isJP2(); //           pp: ADC low 37

    // tower flag QA
    if(isHT0) fHistNFiredHT0vsFlag->Fill(flag);
    if(isHT1) fHistNFiredHT1vsFlag->Fill(flag);
    if(isHT2) fHistNFiredHT2vsFlag->Fill(flag);
    if(isHT3) fHistNFiredHT3vsFlag->Fill(flag);
    if(isJP0) fHistNFiredJP0vsFlag->Fill(flag);
    if(isJP1) fHistNFiredJP1vsFlag->Fill(flag);
    if(isJP2) fHistNFiredJP2vsFlag->Fill(flag);

    // tower ADC QA
    if(isHT0) fHistNFiredHT0vsADC->Fill(adc);
    if(isHT1) fHistNFiredHT1vsADC->Fill(adc);
    if(isHT2) fHistNFiredHT2vsADC->Fill(adc);
    if(isHT3) fHistNFiredHT3vsADC->Fill(adc);
    if(isJP0) fHistNFiredJP0vsADC->Fill(adc);
    if(isJP1) fHistNFiredJP1vsADC->Fill(adc);
    if(isJP2) fHistNFiredJP2vsADC->Fill(adc);

    // print some EMCal Trigger info
    if(fDebugLevel == kDebugEmcTrigger) {
      cout<<"i = "<<i<<"  id = "<<emcTrigID<<"  flag = "<<flag<<"  adc = "<<emcTrig->adc();
      cout<<"  isHT0: "<<isHT0<<"  isHT1: "<<isHT1<<"  isHT2: "<<isHT2<<"  isHT3: "<<isHT3;
      cout<<"  isJP0: "<<isJP0<<"  isJP1: "<<isJP1<<"  isJP2: "<<isJP2<<endl;
    }

    // continue for no HT trigger or if JetPatch
    if(flag == 112) continue;                          // JetPatch triggers
    if(!isHT0 && !isHT1 && !isHT2 && !isHT3) continue; // have a HT trigger

    // get associated tower - minus 1 to get array element
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(emcTrigID-1));
    if(!tower) { cout<<"No tower pointer... iTow = "<<emcTrigID<<endl; continue; }

    // tower ID (1-4800) - get from emcTrigID directly
    //int towerID = tower->id();
    int towerID = -1;
    if( gROOT->GetClass("StPicoBTowHit")->GetClassVersion() < 3) {
      //towerID = tower->id();
    } else {
      //towerID = tower->numericIndex2SoftId(emcTrigID - 1);
    }

    // tower ID is directly related to referenced emcTrigID
    towerID = emcTrigID;
    if(towerID < 0) { cout<<"tower ID < 0, tower ID = "<<towerID<<endl; continue; } // double check these aren't still in the event list

    // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    //double towerPhi = towerPosition.Phi();
    double towerEta = towerPosition.PseudoRapidity();
    double towerE = tower->energy();
    double towerEt = towerE / (1.0*TMath::CosH(towerEta));

    // if(towerEt < 0) cout<<"emcTrigID = "<<emcTrigID<<"  towerID = "<<towerID<<"  towerEta = "<<towerEta<<"  towerE = "<<towerE<<"  towerEt = "<<towerEt<<"  ADC = "<<tower->adc()<<endl;

    // fill some histograms for QA when have a zero energy entry
    if(towerE == 0) {  // for ZERO energy 
      if(isHT1 && (flag == HT1flag)) fHistNZeroEHT1vsID->Fill(towerID);
      if(isHT2 && (flag == HT2flag)) fHistNZeroEHT2vsID->Fill(towerID);
      if(isHT3 && (flag == HT3flag)) fHistNZeroEHT3vsID->Fill(towerID);
    }

    // Fill histogram with towerID when we come across a negative energy entry
    if(towerE < 0) { // for negative energy
      if(isHT1 && (flag == HT1flag)) fHistNNegEHT1vsID->Fill(towerID);
      if(isHT2 && (flag == HT2flag)) fHistNNegEHT2vsID->Fill(towerID);
      if(isHT3 && (flag == HT3flag)) fHistNNegEHT3vsID->Fill(towerID);
      continue;
    }

    // fill for fired triggers
    if(isHT0) fHistNFiredHT0vsID->Fill(emcTrigID);
    if(isHT1 && (flag == HT1flag)) fHistNFiredHT1vsID->Fill(emcTrigID);
    if(isHT2 && (flag == HT2flag)) fHistNFiredHT2vsID->Fill(emcTrigID);
    if(isHT3 && (flag == HT3flag)) fHistNFiredHT3vsID->Fill(emcTrigID);

    if(isHT0) fHistHT0FiredEtvsID->Fill(emcTrigID, towerEt);
    if(isHT1 && (flag == HT1flag)) fHistHT1FiredEtvsID->Fill(emcTrigID, towerEt);
    if(isHT2 && (flag == HT2flag)) fHistHT2FiredEtvsID->Fill(emcTrigID, towerEt);
    if(isHT3 && (flag == HT3flag)) fHistHT3FiredEtvsID->Fill(emcTrigID, towerEt);

    if(isHT0) fHistHT0IDvsFiredEt->Fill(towerEt, emcTrigID);
    if(isHT1 && (flag == HT1flag)) fHistHT1IDvsFiredEt->Fill(towerEt, emcTrigID);
    if(isHT2 && (flag == HT2flag)) fHistHT2IDvsFiredEt->Fill(towerEt, emcTrigID);
    if(isHT3 && (flag == HT3flag)) fHistHT3IDvsFiredEt->Fill(towerEt, emcTrigID);

  } // trigger loop

}
//
// check if event is MB
//_________________________________________________________________________
Bool_t StPicoTrackClusterQA::CheckForMB(int RunFlag, int type) {
  // Run11 triggers:
  int arrMB_Run11[] = {13, 320000, 320001, 320011, 320021, 330021};

  // Run12 (200 GeV pp) triggers:
  // 1) VPDMB
  int arrMB_Run12[] = {370001, 370011, 370983};

  // Run13 triggers:
  int arrMB_Run13[] = {39, 430001, 430011, 430021, 430031};

  // Run14 triggers:
  int arrMB_Run14[] = {450014};
  int arrMB30_Run14[] = {450010, 450020};
  int arrMB5_Run14[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};
  // additional 30: 4, 5, 450201, 450202, 450211, 450212

  // Run16 triggers:
  int arrMB_Run16[] = {520021};
  int arrMB5_Run16[] = {520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
  int arrMB10_Run16[] = {520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};

  // Run17 triggers:
  int arrMB30_Run17[] = {570001, 590001};
  int arrMB100_Run17[] = {590002};
  int arrMBnovtx_Run17[] = {55, 570004};

  // run flag selection to check for MB firing
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        switch(type) {
          case StJetFrameworkPicoBase::kRun14main :
              if((DoComparison(arrMB_Run14, sizeof(arrMB_Run14)/sizeof(*arrMB_Run14)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMB5 :
              if((DoComparison(arrMB5_Run14, sizeof(arrMB5_Run14)/sizeof(*arrMB5_Run14)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMB30 :
              if((DoComparison(arrMB30_Run14, sizeof(arrMB30_Run14)/sizeof(*arrMB30_Run14)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run14, sizeof(arrMB_Run14)/sizeof(*arrMB_Run14)))) { return kTRUE; }
        }
        break; // added May20

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        switch(type) {
          case StJetFrameworkPicoBase::kRun16main :
              if((DoComparison(arrMB_Run16, sizeof(arrMB_Run16)/sizeof(*arrMB_Run16)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMB5 :
              if((DoComparison(arrMB5_Run16, sizeof(arrMB5_Run16)/sizeof(*arrMB5_Run16)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kVPDMB10 :
              if((DoComparison(arrMB10_Run16, sizeof(arrMB10_Run16)/sizeof(*arrMB10_Run16)))) { return kTRUE; }
              break;
          default :
              if((DoComparison(arrMB_Run16, sizeof(arrMB_Run16)/sizeof(*arrMB_Run16)))) { return kTRUE; }
        }
        break; // added May20

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
              if((DoComparison(arrMB_Run12, sizeof(arrMB_Run12)/sizeof(*arrMB_Run12)))) { return kTRUE; }
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

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17 pp
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
Bool_t StPicoTrackClusterQA::CheckForHT(int RunFlag, int type) {
  // Run12 (200 GeV pp) triggers:
  int arrHT1_Run12[] = {370511, 370546};
  int arrHT2_Run12[] = {370521, 370522, 370531, 370980};
  int arrHT3_Run12[] = {380206, 380216}; // NO HT3 triggers in this dataset

  // Run14 triggers:
  int arrHT1_Run14[] = {450201, 450211, 460201};
  int arrHT2_Run14[] = {450202, 450212, 460202, 460212};
  int arrHT3_Run14[] = {450203, 450213, 460203};

  // Run16 triggers:
  int arrHT1_Run16[] = {520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
  int arrHT2_Run16[] = {530202, 540203};
  int arrHT3_Run16[] = {520203, 530213};

  // Run17 triggers: (HT1 and HT2 not exclusive)
  int arrHT1_Run17[] = {29, 570204, 570214};
  int arrHT2_Run17[] = {30, 31, 570205, 570215};
  int arrHT3_Run17[] = {16, 570201, 590201};

  // run flag selection to check for MB firing
  switch(RunFlag) {
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14 AuAu
        switch(type) {
          case StJetFrameworkPicoBase::kIsHT1 :
              if((DoComparison(arrHT1_Run14, sizeof(arrHT1_Run14)/sizeof(*arrHT1_Run14)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT2 :
              if((DoComparison(arrHT2_Run14, sizeof(arrHT2_Run14)/sizeof(*arrHT2_Run14)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT3 :
              if((DoComparison(arrHT3_Run14, sizeof(arrHT3_Run14)/sizeof(*arrHT3_Run14)))) { return kTRUE; }
              break;
          default :  // default to HT2
              if((DoComparison(arrHT2_Run14, sizeof(arrHT2_Run14)/sizeof(*arrHT2_Run14)))) { return kTRUE; }
        }
        break; // added May20

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16 AuAu
        switch(type) {
          case StJetFrameworkPicoBase::kIsHT1 :
              if((DoComparison(arrHT1_Run16, sizeof(arrHT1_Run16)/sizeof(*arrHT1_Run16)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT2 :
              if((DoComparison(arrHT2_Run16, sizeof(arrHT2_Run16)/sizeof(*arrHT2_Run16)))) { return kTRUE; }
              break;
          case StJetFrameworkPicoBase::kIsHT3 :
              if((DoComparison(arrHT3_Run16, sizeof(arrHT3_Run16)/sizeof(*arrHT3_Run16)))) { return kTRUE; }
              break;
          default :  // Run16 only has HT1's
              if((DoComparison(arrHT1_Run16, sizeof(arrHT1_Run16)/sizeof(*arrHT1_Run16)))) { return kTRUE; }
        }
        break; // added May20

    case StJetFrameworkPicoBase::Run12_pp200 : // Run12 pp
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

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17 pp
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
//
//____________________________________________________________________________________________
Bool_t StPicoTrackClusterQA::IsTowerOK( Int_t mTowId ){
  //if( badTowers.size()==0 ){
  if( badTowers.empty() ){
    __ERROR("StPicoTrackClusterQA::IsTowerOK: WARNING: You're trying to run without a bad tower list. If you know what you're doing, deactivate this throw and recompile.");
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
//
//____________________________________________________________________________________________
Bool_t StPicoTrackClusterQA::IsTowerDead( Int_t mTowId ){
  //if( deadTowers.size()==0 ){
  if( deadTowers.empty() ){
    __ERROR("StPicoTrackClusterQA::IsTowerDead: WARNING: You're trying to run without a dead tower list. If you know what you're doing, deactivate this throw and recompile.");
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
//
//____________________________________________________________________________
void StPicoTrackClusterQA::ResetBadTowerList( ){
  badTowers.clear();
}

// Add bad towers from comma separated values file
// Can be split into arbitrary many lines
// Lines starting with # will be ignored
Bool_t StPicoTrackClusterQA::AddBadTowers(TString csvfile){
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
Bool_t StPicoTrackClusterQA::AddDeadTowers(TString csvfile){
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
//
//____________________________________________________________________________
void StPicoTrackClusterQA::ResetDeadTowerList( ){
  deadTowers.clear();
}
//
// Returns pt of hardest track in the event
//
Double_t StPicoTrackClusterQA::GetMaxTrackPt()
{
  // get # of tracks
  int nTrack = mPicoDst->numberOfTracks();
  double fMaxTrackPt = -99;

  // loop over all tracks
  for(int i=0; i<nTrack; i++) {
    // get tracks
    StPicoTrack* track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // primary track switch
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
// fill trigger ids of the dataset into a histogram
// only currently set up for run 12 pp (200 GeV)
//________________________________________________________________________
void StPicoTrackClusterQA::FillTriggerIDs(TH1 *h) {
  // All non-test triggers for Run12
  // 27 entries, array goes from 0 - 26
  unsigned int triggers[] = {370001, 370011, 370021, 370022, 370031, 370032, 370301, 370341, 370361, 370501, 370511, 370521, 370522, 370531, 370541, 370542, 370546, 370601, 370611, 370621, 370641, 370701, 370801, 370980, 370981, 370982, 370983};

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    // check for valid, non-test trigger ID
    if(mytriggers[i] > 1000) {
      for(int j=0; j<27; j++) {
        if(mytriggers[i] == triggers[j]) h->Fill(j + 1);

      } // loops over ID's
    }   // non-test trigger
  }     // loop over triggers

  // label bins of the analysis trigger selection summary
  for(int i = 0; i < 27; i++) {
    h->GetXaxis()->SetBinLabel(i+1, Form("%i", triggers[i]));
  }

  // set x-axis labels vertically
  h->LabelsOption("v");

}
//
// function to do some event QA
//________________________________________________________________________________________
void StPicoTrackClusterQA::RunEventQA() {
  // get run ID, transform to run order for filling histogram of corrections
  int fRunId = mPicoEvent->runId();
  int RunId_Order = GetRunNo(fRunId);// + 1;
  //if(RunId_Order < -1) return kStOK;

  // centrality, refmult
  int refmult = mPicoEvent->refMult();

  // vertex
  TVector3 fPrimaryVertex = mPicoEvent->primaryVertex();
  float fZVtx = fPrimaryVertex.z();
  float fYVtx = fPrimaryVertex.y();
  float fXVtx = fPrimaryVertex.x();
  float fVzVPD = mPicoEvent->vzVpd();
  float ranking = mPicoEvent->ranking();

  // coincidence rates
  float fZDCx = mPicoEvent->ZDCx();
  float fBBCx = mPicoEvent->BBCx();

  // print
  //if(TMath::Abs(fVzVPD) > 200)  cout<<"RefMult: "<<refmult<<"  Ranking: "<<ranking<<"  ZVtx: "<<fZVtx<<"  fVzVPD: "<<fVzVPD<<"  BBCx: "<<fBBCx<<"  ZDCx: "<<fZDCx<<endl;

  // fill histograms
  fProfEventRefMult->Fill(RunId_Order + 1., refmult);
  fProfEventRanking->Fill(RunId_Order + 1., ranking);
  fProfEventZvtx->Fill(RunId_Order + 1., fZVtx);
  fProfEventYvtx->Fill(RunId_Order + 1., fYVtx);
  fProfEventXvtx->Fill(RunId_Order + 1., fXVtx);
  if(fVzVPD > -900.) fProfEventVzVPD->Fill(RunId_Order + 1., fVzVPD); // sometimes this is not set, don't include in average then
  fProfEventBBCx->Fill(RunId_Order + 1., fBBCx);
  fProfEventZDCx->Fill(RunId_Order + 1., fZDCx);

}
//
// this function checks for the bin number of the run from a runlist header 
// in order to apply various corrections and fill run-dependent histograms
// _________________________________________________________________________________
Int_t StPicoTrackClusterQA::GetRunNo(int runid){
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
