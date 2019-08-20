// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// track and tower input to cluster together using fastjet wrapper and create jets
//      - leading jet tag
//      - access to jet constituents
//      - general QA
//      - constitutent subtractor
//
// ################################################################
// $Id$

#include "StJetMakerTask.h"

// ROOT includes
#include "TROOT.h"
#include <TChain.h>
#include <TClonesArray.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "TFile.h"
#include "TVector3.h"
#include <sstream>
#include <fstream>

// StRoot includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoArrays.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoConstants.h"

// for clusters
#include "StEmcUtil/geometry/StEmcGeom.h"
//#include "StEmcUtil/projection/StEmcPosition.h" // old - uses StThreeVectorF
//class StEmcPosition; // old - uses StThreeVectorF
#include "StEmcPosition2.h"
class StEmcPosition2;

#include "StCentMaker.h"

// jet class and fastjet wrapper and dataset (Run#'s) 
#include "StJet.h"
#include "StFJWrapper.h"
#include "StJetFrameworkPicoBase.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

// extra includes
#include "StJetPicoDefinitions.h"

// StRoot classes
class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;

// constants - definitions
const Int_t StJetMakerTask::fgkConstIndexShift = 100000;

ClassImp(StJetMakerTask)

//________________________________________________________________________
StJetMakerTask::StJetMakerTask() : 
  StMaker(),
  doWriteHistos(kFALSE),
  doUsePrimTracks(kFALSE), 
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE), 
  fRequireCentSelection(kFALSE),
  doConstituentSubtr(kFALSE),
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fMaxEventTrackPt(30.0),
  fMaxEventTowerEt(1000.0), // 30.0
  doRejectBadRuns(kFALSE),
  Bfield(0.0),
//  mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fEmcTriggerEventType(0), // see StJetFrameworkPicoBase::fEmcTriggerFlagEnum
  fMBEventType(2),         // kVPDMB5, see StJetFrameworkPicoBase::fMBFlagEnum
  fTriggerToUse(0),        // kTriggerAny, see StJetFrameworkPicoBase::fTriggerEventTypeEnum
  fCentralityScaled(0.0),
  ref16(-99),
  ref9(-99), 
  mOutName(""),
  fTracksName(""),
  fCaloName(""),
  fJetsName(""),
  fJetAlgo(1), 
  fJetType(0), 
  fRecombScheme(fastjet::BIpt2_scheme), // was BIpt_scheme
  fjw("StJetMakerTask", "StJetMakerTask"),
  fRadius(0.4),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10.),
  fJetPhiMax(+10.),
  fJetEtaMin(-0.6),
  fJetEtaMax(0.6),
  fGhostArea(0.005), 
  fMinJetTrackPt(0.2),
  fMaxJetTrackPt(30.0),
  fMinJetClusPt(0.15),
  fMinJetClusE(0.2),
  fMinJetTowerE(0.2),
  fTrackEtaMin(-1.0),
  fTrackEtaMax(1.0),
  fTrackPhiMin(0.0),
  fTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackEtaMin(-1.0),
  fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0),
  fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fJetTowerEMin(0.2),
  fJetTowerEMax(100.0),
  fJetTowerEtaMin(-1.0),
  fJetTowerEtaMax(1.0),
  fJetTowerPhiMin(0.0),
  fJetTowerPhiMax(2.0*TMath::Pi()),
  mTowerEnergyMin(0.2),
  mHadronicCorrFrac(1.),
  fJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks), // default is using all matched Tracks, Aug2019, per Hanseul
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0x0),
  fJetsBGsub(0x0),
  fFull_Event(0),
  fConstituents(0),
  fJetsConstit(0x0),
  fJetsConstitBGsub(0x0),
  mGeom(StEmcGeom::instance("bemc")),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0)
{
  // Default constructor.
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = kFALSE; }

  for(int i=0; i<4801; i++) { 
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE; 

    mTowerMatchTrkIndexLast[i] = -1; // 0
    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

}

//________________________________________________________________________
StJetMakerTask::StJetMakerTask(const char *name, double mintrackPt = 0.20, bool doHistos = kFALSE, const char* outName = "") : 
  StMaker(name),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  doppAnalysis(kFALSE),
  fRequireCentSelection(kFALSE),
  doConstituentSubtr(kFALSE),
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  fCentralitySelectionCut(-99),
  fMaxEventTrackPt(30.0),
  fMaxEventTowerEt(1000.0), // 30.0
  doRejectBadRuns(kFALSE),
  Bfield(0.0),
//  mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fEmcTriggerEventType(0), // see StJetFrameworkPicoBase::fEMCTriggerFlagEnum
  fMBEventType(2),         // kVPDMB5, see StJetFrameworkPicoBase::fMBFlagEnum
  fTriggerToUse(0),        // kTriggerAny, see StJetFrameworkPicoBase::fTriggerEventTypeEnum
  fCentralityScaled(0.0),
  ref16(-99),
  ref9(-99),
  mOutName(outName),
  fTracksName("Tracks"),
  fCaloName("Clusters"),
  fJetsName("Jets"),
  fJetAlgo(1), 
  fJetType(0),
  fRecombScheme(fastjet::BIpt2_scheme), // was BIpt2_scheme
  fjw(name, name),
  fRadius(0.4),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10), 
  fJetPhiMax(+10),
  fJetEtaMin(-0.6), 
  fJetEtaMax(0.6),
  fGhostArea(0.005),
  fMinJetTrackPt(mintrackPt),
  fMaxJetTrackPt(30.0), 
  fMinJetClusPt(0.15),
  fMinJetClusE(0.2),
  fMinJetTowerE(0.2),
  fTrackEtaMin(-1.0),
  fTrackEtaMax(1.0),
  fTrackPhiMin(0.0),
  fTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackEtaMin(-1.0), 
  fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0),
  fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fJetTowerEMin(0.2),
  fJetTowerEMax(100.0),
  fJetTowerEtaMin(-1.0),
  fJetTowerEtaMax(1.0),
  fJetTowerPhiMin(0.0),
  fJetTowerPhiMax(2.0*TMath::Pi()),
  mTowerEnergyMin(0.2),
  mHadronicCorrFrac(1.),
  fJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks), // default is using all matched Tracks, Aug2019, per Hanseul
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0x0),
  fJetsBGsub(0x0),
  fFull_Event(0),
  fConstituents(0),
  fJetsConstit(0x0),
  fJetsConstitBGsub(0x0),
  mGeom(StEmcGeom::instance("bemc")),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  mCentMaker(0x0),
  mBaseMaker(0x0),
  mEmcPosition(0x0),
  grefmultCorr(0x0)
{
  // Standard constructor.
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = kFALSE; }

  for(int i=0; i<4801; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;

    mTowerMatchTrkIndexLast[i] = -1; // 0
    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StJetMakerTask::~StJetMakerTask()
{
  // Destructor
  //fJets->Clear(); delete fJets;
  if(fHistMultiplicity)        delete fHistMultiplicity;
  if(fHistCentrality)          delete fHistCentrality;
  if(fHistCentralityPostCut)   delete fHistCentralityPostCut;
  if(fHistFJRho)               delete fHistFJRho;
  if(fProfEventBBCx)           delete fProfEventBBCx;
  if(fProfEventZDCx)           delete fProfEventZDCx;

  if(fHistNTrackvsPt)          delete fHistNTrackvsPt;
  if(fHistNTrackvsPhi)         delete fHistNTrackvsPhi;
  if(fHistNTrackvsEta)         delete fHistNTrackvsEta;
  if(fHistNTrackvsPhivsEta)    delete fHistNTrackvsPhivsEta;
  if(fHistNTowervsID)          delete fHistNTowervsID;
  if(fHistNTowervsADC)         delete fHistNTowervsADC;
  if(fHistNTowervsE)           delete fHistNTowervsE;
  if(fHistNTowervsEt)          delete fHistNTowervsEt;
  if(fHistNTowervsPhi)         delete fHistNTowervsPhi;
  if(fHistNTowervsEta)         delete fHistNTowervsEta;
  if(fHistNTowervsPhivsEta)    delete fHistNTowervsPhivsEta;

  for(int i = 0; i < 5; i++) {
    if(fHistNMatchTrack[i])       delete fHistNMatchTrack[i];
    if(fHistHadCorrComparison[i]) delete fHistHadCorrComparison[i];
    if(fHistHadCorrComparisonLast[i]) delete fHistHadCorrComparisonLast[i];
  }
  if(fHistJetNTrackvsPt)       delete fHistJetNTrackvsPt;
  if(fHistJetNTrackvsPhi)      delete fHistJetNTrackvsPhi;
  if(fHistJetNTrackvsEta)      delete fHistJetNTrackvsEta;
  if(fHistJetNTrackvsPhivsEta) delete fHistJetNTrackvsPhivsEta;
  if(fHistJetNTowervsID)       delete fHistJetNTowervsID;
  if(fHistJetNTowervsADC)      delete fHistJetNTowervsADC;
  if(fHistJetNTowervsE)        delete fHistJetNTowervsE;
  if(fHistJetNTowervsEt)       delete fHistJetNTowervsEt;
  if(fHistJetNTowervsPhi)      delete fHistJetNTowervsPhi;
  if(fHistJetNTowervsEta)      delete fHistJetNTowervsEta;
  if(fHistJetNTowervsPhivsEta) delete fHistJetNTowervsPhivsEta;

  if(fHistNJetsvsPt)           delete fHistNJetsvsPt;
  if(fHistNJetsvsPhi)          delete fHistNJetsvsPhi;
  if(fHistNJetsvsEta)          delete fHistNJetsvsEta;
  if(fHistNJetsvsPhivsEta)     delete fHistNJetsvsPhivsEta;
  if(fHistNJetsvsArea)         delete fHistNJetsvsArea;
  if(fHistNJetsvsNConstituents)delete fHistNJetsvsNConstituents;
  if(fHistNJetsvsNTracks)      delete fHistNJetsvsNTracks;
  if(fHistNJetsvsNTowers)      delete fHistNJetsvsNTowers;

  if(fHistQATowIDvsEta)        delete fHistQATowIDvsEta;
  if(fHistQATowIDvsPhi)        delete fHistQATowIDvsPhi;

  if(mEmcPosition)             delete mEmcPosition;
}
//
//
//________________________________________________________________________
Int_t StJetMakerTask::Init() {
  DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  fJets->SetName(fJetsName);

  fJetsBGsub = new TClonesArray("StJet");
  fJetsBGsub->SetName(fJetsName+"BGsub");

  // may need array (name hard-coded, Feb20, 2018)
  fJetsConstit = new TClonesArray("StPicoTrack");
  fJetsConstit->SetName("JetConstituents");

  // TODO maybe add BGsub?

  // position object for Emc
  mEmcPosition = new StEmcPosition2();

  // ============================ set jet parameters for fastjet wrapper  =======================
  // recombination schemes:
  // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme
  fastjet::RecombinationScheme    recombScheme;
  if (fRecombScheme == 0)     recombScheme = fastjet::E_scheme;
  if (fRecombScheme == 1)     recombScheme = fastjet::pt_scheme; 
  if (fRecombScheme == 2)     recombScheme = fastjet::pt2_scheme;
  if (fRecombScheme == 3)     recombScheme = fastjet::Et_scheme;
  if (fRecombScheme == 4)     recombScheme = fastjet::Et2_scheme;
  if (fRecombScheme == 5)     recombScheme = fastjet::BIpt_scheme;
  if (fRecombScheme == 6)     recombScheme = fastjet::BIpt2_scheme;
  if (fRecombScheme == 7)     recombScheme = fastjet::WTA_pt_scheme;
  if (fRecombScheme == 8)     recombScheme = fastjet::WTA_modp_scheme;
  if (fRecombScheme == 99)    recombScheme = fastjet::external_scheme;

  // jet algorithm
  fastjet::JetAlgorithm          algorithm;
  if (fJetAlgo == 1)      algorithm = fastjet::antikt_algorithm;
  if (fJetAlgo == 0)      algorithm = fastjet::kt_algorithm;
  // extra algorithms
  if (fJetAlgo == 2)      algorithm = fastjet::cambridge_algorithm;
  if (fJetAlgo == 3)      algorithm = fastjet::genkt_algorithm;
  if (fJetAlgo == 11)     algorithm = fastjet::cambridge_for_passive_algorithm;
  if (fJetAlgo == 13)     algorithm = fastjet::genkt_for_passive_algorithm;
  if (fJetAlgo == 99)     algorithm = fastjet::plugin_algorithm;
  if (fJetAlgo == 999)    algorithm = fastjet::undefined_jet_algorithm;
  fastjet::Strategy               strategy = fastjet::Best;

  // setup fj wrapper
  fjw.SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw.SetStrategy(strategy);
  fjw.SetGhostArea(fGhostArea);
  fjw.SetR(fRadius);
  fjw.SetAlgorithm(algorithm);        //fJetAlgo);
  fjw.SetRecombScheme(recombScheme);  //fRecombScheme);
  fjw.SetMaxRap(1.2);

  // ghost-area specifications
  double ghost_maxrap = 1.2;
  fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
  fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

  // setting legacy mode
  //if(fLegacyMode) { fjw.SetLegacyMode(kTRUE); }

  return kStOK;
}
//
//
//_______________________________________________________________________________________
Int_t StJetMakerTask::Finish() {
/*
  //  Summarize the run.
  cout<<"StJetMakerTask::Finish()\n";
  cout<<"##### Jet parameters overview #####"<<endl;
  cout<<"type = "<<"   algorithm = "<<"  recombination scheme = "<<endl;
  cout<<"R = "<<fRadius<<"   ghostArea = "<<fGhostArea<<"  minJetArea = "<<fMinJetArea<<endl;
  cout<<"minTrackPt = "<<fMinJetTrackPt<<"  minClusterPt = "<<fMinJetClusPt<<"  maxTrackPt = "<<fMaxJetTrackPt<<endl;
  cout<<"minJetPhi = "<<fJetPhiMin<<"  maxJetPhi = "<<fJetPhiMax<<"  minJetEta = "<<fJetEtaMin<<"  maxJetEta = "<<fJetEtaMax<<endl;
  cout<<"End of StJetMakerTask::Finish"<<endl;
*/

  if(doWriteHistos && mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  return kStOK;
}
//
// Function: declare histograms
//________________________________________________________________________
void StJetMakerTask::DeclareHistograms() {
  // constants
  double pi = 1.0*TMath::Pi();

  // histogram settings:     pp : AuAu
  double kHistMultMax = (doppAnalysis) ? 100. : 800.;
  int kHistMultBins = (doppAnalysis) ? 100. : 400.;

  // set binning for run based corrections - run dependent
  Int_t nRunBins = 1; // - just a default
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   nRunBins = 857 + 43;
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) nRunBins = 830; //1654;
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) nRunBins = 1359;
  Double_t nRunBinsMax = (Double_t)nRunBins + 0.5;

  //fHistEventCounter = new TH1F("fHistEventCounter", "Event counter", 10, 0.5, 10.5);
  fHistMultiplicity = new TH1F("fHistMultiplicity", "No. events vs multiplicity", kHistMultBins, 0, kHistMultMax);
  fHistCentrality = new TH1F("fHistCentrality", "No. events vs centrality", 20, 0, 100);    
  fHistCentralityPostCut = new TH1F("fHistCentralityPostCut", "No. events vs centrality, after cut", 20, 0, 100);
  fHistFJRho = new TH1F("fHistFJRho", "Underlying event energy density via FastJet", 200, 0, 50);
 
  // event QA
  fProfEventBBCx = new TProfile("fProfEventBBCx", "Event averaged BBC coincidence rate", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);
  fProfEventZDCx = new TProfile("fProfEventZDCx", "Event averaged ZDC coincidence rate", nRunBins, 0.5, nRunBinsMax);//, -100., 100.);

  // track and tower plots
  fHistNTrackvsPt = new TH1F("fHistNTrackvsPt", "# track vs p_{T}", 200, 0., 40.);
  fHistNTrackvsPhi = new TH1F("fHistNTrackvsPhi", "# track vs #phi", 72, 0., 2.0*pi);
  fHistNTrackvsEta = new TH1F("fHistNTrackvsEta", "# track vs #eta", 40, -1.0, 1.0);
  fHistNTrackvsPhivsEta = new TH2F("fHistNTrackvsPhivsEta", "# track vs #phi vs #eta", 144, 0, 2.0*pi, 20, -1.0, 1.0);
  fHistNTowervsID = new TH1F("fHistNTowervsID", "# tower vs tower id", 4800, 0.5, 4800.5);
  fHistNTowervsADC = new TH1F("fHistNTowervsADC", "# tower vs ADC", 250, 0., 250.);
  fHistNTowervsE = new TH1F("fHistNTowervsE", "# tower vs energy", 200, 0., 40.0);
  fHistNTowervsEt = new TH1F("fHistNTowervsEt", "# tower vs transverse energy", 200, 0., 40.0);
  fHistNTowervsPhi = new TH1F("fHistNTowervsPhi", "# tower vs #phi", 144, 0., 2.0*pi);
  fHistNTowervsEta = new TH1F("fHistNTowervsEta", "# tower vs #eta", 40, -1.0, 1.0);
  fHistNTowervsPhivsEta = new TH2F("fHistNTowervsPhivsEta", "# vs #phi vs #eta", 144, 0, 2.0*pi, 20, -1.0, 1.0);

  // jet track and jet tower plots
  for(int i = 0; i < 5; i++) {
    fHistNMatchTrack[i] = new TH1F(Form("fHistNMatchTrack_%i", i), "Number of matched tracks to tower, cent bin: %i", 11, -0.5, 10.5);
    fHistHadCorrComparison[i] = new TH2F(Form("fHistHadCorrComparison_%i", i), "max E_{T} vs sum E_{T} of matched tracks for had corr, cent bin: %i", 120, 0., 30., 160, -10., 30.); 
    fHistHadCorrComparisonLast[i] = new TH2F(Form("fHistHadCorrComparisonLast_%i", i), "last E_{T} vs max E_{T} of matched tracks for had corr, cent bin: %i", 120, 0., 30., 160, -10., 30.);
  }
  fHistJetNTrackvsPt = new TH1F("fHistJetNTrackvsPt", "Jet track constituents vs p_{T}", 150, 0., 30.);
  fHistJetNTrackvsPhi = new TH1F("fHistJetNTrackvsPhi", "Jet track constituents vs #phi", 72, 0., 2.0*pi);
  fHistJetNTrackvsEta = new TH1F("fHistJetNTrackvsEta", "Jet track constituents vs #eta", 40, -1.0, 1.0);
  fHistJetNTrackvsPhivsEta = new TH2F("fHistJetNTrackvsPhivsEta", "Jet track constituents vs #phi vs #eta", 144, 0, 2.0*pi, 20, -1.0, 1.0);
  fHistJetNTowervsID = new TH1F("fHistJetNTowervsID", "Jet tower vs tower id", 4800, 0.5, 4800.5);
  fHistJetNTowervsADC = new TH1F("fHistJetNTowervsADC", "Jet tower constituents vs ADC", 250, 0., 250.);
  fHistJetNTowervsE = new TH1F("fHistJetNTowervsE", "Jet tower constituents vs energy", 150, 0., 30.0);
  fHistJetNTowervsEt = new TH1F("fHistJetNTowervsEt", "Jet tower constituents vs transverse energy", 150, 0., 30.0);
  fHistJetNTowervsPhi = new TH1F("fHistJetNTowervsPhi", "Jet tower constituents vs #phi", 144, 0., 2.0*pi);
  fHistJetNTowervsEta = new TH1F("fHistJetNTowervsEta", "Jet tower constituents vs #eta", 40, -1.0, 1.0);
  fHistJetNTowervsPhivsEta = new TH2F("fHistJetNTowervsPhivsEta", "Jet tower constituents vs #phi vs #eta", 144, 0, 2.0*pi, 20, -1.0, 1.0);

  // jet plots
  fHistNJetsvsPt = new TH1F("fHistNJetsvsPt", "NJets vs p_{T}", 100, 0., 100.0);
  fHistNJetsvsPhi = new TH1F("fHistNJetsvsPhi", "NJets vs #phi", 144, 0., 2.0*pi);
  fHistNJetsvsEta = new TH1F("fHistNJetsvsEta", "NJets vs #eta", 40, -1.0, 1.0);
  fHistNJetsvsPhivsEta = new TH2F("fHistNJetsvsPhivsEta", "NJets vs #phi vs #eta", 144, 0, 2.0*pi, 20, -1.0, 1.0);
  fHistNJetsvsArea = new TH1F("fHistNJetsvsArea", "NJets vs jet area", 100, 0.0, 1.0);
  fHistNJetsvsNConstituents = new TH1F("fHistNJetsvsNConstituents", "NJets vs NConstit", 51, -0.5, 50.5);
  fHistNJetsvsNTracks = new TH1F("fHistNJetsvsNTracks", "NJets vs NTracks", 51, -0.5, 50.5);
  fHistNJetsvsNTowers = new TH1F("fHistNJetsvsNTowers", "NJets vs NTowers", 51, -0.5, 50.5);

  fHistQATowIDvsEta = new TH2F("fHistQATowIDvsEta", "Tower ID vs #eta", 4800, 0.5, 4800.5, 40, -1.0, 1.0);
  fHistQATowIDvsPhi = new TH2F("fHistQATowIDvsPhi", "Tower ID vs #phi", 4800, 0.5, 4800.5, 144, 0.0, 2.0*pi);

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}
//
// Write histograms
//________________________________________________________________________
void StJetMakerTask::WriteHistograms() {
  fHistMultiplicity->Write();
  fHistCentrality->Write();
  fHistCentralityPostCut->Write();
  fHistFJRho->Write();
  fProfEventBBCx->Write();
  fProfEventZDCx->Write();

  fHistNTrackvsPt->Write();
  fHistNTrackvsPhi->Write();
  fHistNTrackvsEta->Write();
  fHistNTrackvsPhivsEta->Write();
  fHistNTowervsID->Write();
  fHistNTowervsADC->Write();
  fHistNTowervsE->Write();
  fHistNTowervsEt->Write();
  fHistNTowervsPhi->Write();
  fHistNTowervsEta->Write();
  fHistNTowervsPhivsEta->Write();

  for(int i = 0; i < 5; i++) {
    fHistNMatchTrack[i]->Write();
    fHistHadCorrComparison[i]->Write();
    fHistHadCorrComparisonLast[i]->Write();
  }
  fHistJetNTrackvsPt->Write();
  fHistJetNTrackvsPhi->Write();
  fHistJetNTrackvsEta->Write();
  fHistJetNTrackvsPhivsEta->Write();
  fHistJetNTowervsID->Write();
  fHistJetNTowervsADC->Write();
  fHistJetNTowervsE->Write();
  fHistJetNTowervsEt->Write();
  fHistJetNTowervsPhi->Write();
  fHistJetNTowervsEta->Write();
  fHistJetNTowervsPhivsEta->Write();

  fHistNJetsvsPt->Write();
  fHistNJetsvsPhi->Write();
  fHistNJetsvsEta->Write();
  fHistNJetsvsPhivsEta->Write();
  fHistNJetsvsArea->Write();
  fHistNJetsvsNConstituents->Write();
  fHistNJetsvsNTracks->Write();
  fHistNJetsvsNTowers->Write();

  fHistQATowIDvsEta->Write();
  fHistQATowIDvsPhi->Write();
}
//
// Function: clear or delete objects after running
//________________________________________________________________________________
void StJetMakerTask::Clear(Option_t *opt) {
  fJets->Clear();
  fJetsBGsub->Clear();
}
//
// Function: main loop, called for each event
//________________________________________________________________________________
int StJetMakerTask::Make()
{
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // ZERO's out the jet array
  fJets->Delete();
  fJetsBGsub->Delete();

  // ZERO these out for double checking they aren't set
  for(int i = 0; i < 4801; i++) {
    mTowerMatchTrkIndexLast[i] = -1; // 0
    for(int j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

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

  // get base class pointer - this class does not inherit from base class: StJetFrameworkPicoBase, but we want to reduce redundancy
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

  // cut event on max tower E > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  // Z-vertex cut - - per the Aj analysis (-40, 40)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

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

  // multiplicity histogram - event activity - compensate for pp or AuAu
  double kEventActivity = (doppAnalysis) ? (double)grefMult : refCorr2;
  fHistMultiplicity->Fill(kEventActivity);

  // scaled centrality - based on 5% bins
  fHistCentrality->Fill(fCentralityScaled);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!mBaseMaker->SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOK; }
  fHistCentralityPostCut->Fill(fCentralityScaled);

  // ============================ Triggers ============================== //
  // check for MB and HT triggers - Type Flag corresponds to selected type of MB or EMC
  // NEED to ADD new triggers and runs to StJetFrameworkPicoBase class !!
  // different access method because this class doesn't inherit from base
  bool fHaveMBevent    = mBaseMaker->CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event   = mBaseMaker->CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event  = mBaseMaker->CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = mBaseMaker->CheckForHT(fRunFlag, fEmcTriggerEventType);
  bool fHaveAnyEvent   = kTRUE;

  // fill trigger array
  FillEmcTriggersArr();

  // run event QA function
  RunEventQA();

  // no need for switch for few checks - FIXME
  if((fTriggerToUse == StJetFrameworkPicoBase::kTriggerMB) && (!fHaveMB5event) && (!fHaveMB30event)) return kStOK;  // MB triggered event
  if((fTriggerToUse == StJetFrameworkPicoBase::kTriggerHT) && (!fHaveEmcTrigger))                    return kStOK;  // HT triggered event
  // else fTriggerToUse is ANY and we still want to run analysis

  // switches for QA analysis
  bool fRunForMB = kFALSE;  // used to differentiate pp and AuAu
  if(doppAnalysis)  fRunForMB = (fHaveMBevent) ? kTRUE : kFALSE;                    // pp analysis
  if(!doppAnalysis) fRunForMB = (fHaveMB5event || fHaveMB30event) ? kTRUE : kFALSE; // NON-pp analysis
  bool fHaveHT2   = mBaseMaker->CheckForHT(fRunFlag, StJetFrameworkPicoBase::kIsHT2);  // HT2
  bool fHaveHT3   = mBaseMaker->CheckForHT(fRunFlag, StJetFrameworkPicoBase::kIsHT3);  // HT3
  bool fHaveMB30HT2HT3 = (fHaveMB30event || fHaveHT2 || fHaveHT3) ? kTRUE : kFALSE;
  bool fHaveMBHT2HT3 = (doppAnalysis && (fRunForMB || fHaveHT2 || fHaveHT3)) ? kTRUE : kFALSE;

  // used for test comparison (July 2019)
  //bool continueON = (doppAnalysis) ? fHaveMBHT2HT3 : fHaveMB30HT2HT3;
  //if(!continueON) return kStOK;

  // Find jets:  deprecated version -> FindJets(tracks, clus, fJetAlgo, fRadius);
  FindJets();

  // Fill jet branch
  if(!doConstituentSubtr) FillJetBranch();
  if( doConstituentSubtr) FillJetBGBranch();

  return kStOK;
}
//
// old class to FindJets - it is deprecated, but kept for backwards compatibility
// the parameters are global so they don't do anything here
//________________________________________________________________________
void StJetMakerTask::FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius)
{
  // call main (new) FindJets function
  FindJets();
}
//
// Main class to FindJets from tracks + towers
//________________________________________________________________________
void StJetMakerTask::FindJets()
{
  // clear out existing wrapper object
  fjw.Clear();
  fFull_Event.clear();

  // get cent bin for some histograms: TODO event activity - compensate for pp or AuAu
  Int_t cbin = -1;
  if (fCentralityScaled >= 0 && fCentralityScaled < 10)       cbin = 1; //  0-10%
  else if (fCentralityScaled >= 10 && fCentralityScaled < 20) cbin = 2; // 10-20%
  else if (fCentralityScaled >= 20 && fCentralityScaled < 30) cbin = 3; // 20-30%
  else if (fCentralityScaled >= 30 && fCentralityScaled < 50) cbin = 4; // 30-50%
  else if (fCentralityScaled >= 50 && fCentralityScaled < 80) cbin = 5; // 50-80%
  
  // assume neutral pion mass
  // additional parameters constructed
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();

  // loop over ALL tracks in PicoDst and add to jet, after acceptance and quality cuts 
  if((fJetType == kFullJet) || (fJetType == kChargedJet)) {
    for(unsigned short iTracks = 0; iTracks < ntracks; iTracks++){
      // get track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(iTracks));
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

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
      if(phi < 0.0)    phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;
      double eta = mTrkMom.PseudoRapidity();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      double p = mTrkMom.Mag();
      double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
      short charge = trk->charge();         
      //cout<<"Charge: "<<charge<<"  nHitsFit: "<<trk->nHitsFit()<<endl;

      // fill some QA histograms
      fHistNTrackvsPt->Fill(pt);
      fHistNTrackvsPhi->Fill(phi);
      fHistNTrackvsEta->Fill(eta);
      fHistNTrackvsPhivsEta->Fill(phi, eta);

      // send track info to FJ wrapper
      //fjw.AddInputVector(px, py, pz, p, iTracks);    // p -> E
      fjw.AddInputVector(px, py, pz, energy, iTracks); // includes E

      // if running constituent subtractor, do this in parallel
      if(doConstituentSubtr) {
        fastjet::PseudoJet particle_Track(px, py, pz, energy);
        particle_Track.set_user_index(iTracks);
        fFull_Event.push_back(particle_Track);
      }
    } // track loop

  }   // if full/charged jets

  // full or neutral jets - get towers and apply hadronic correction
  if((fJetType == kFullJet) || (fJetType == kNeutralJet)) {
    // towerStatus array - globally declared
    //double mTowerMatchTrkIndex[4801][7] = { -1 };
    //int mTowerStatusArr[4801] = { 0 };
    //double mTowerMatchTrkIndexLast[4801] = { -1 };
    int matchedTowerTrackCounter = 0;

    // print
    //int nTracks = mPicoDst->numberOfTracks();               // number of tracks
    //int nTrigs = mPicoDst->numberOfEmcTriggers();           // number of Emc triggers
    //int nBTowHits = mPicoDst->numberOfBTowHits();           // barrel tower hits, always 4800
    int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits(); // number of BEMC matched tracks
    //cout<<"nTracks = "<<nTracks<<"  nTrigs = "<<nTrigs<<"  nBTowHits = "<<nBTowHits<<"  nBEmcPidTraits = "<<nBEmcPidTraits<<endl;

    // loop over ALL clusters in PicoDst to get track<->tower matches saved to arrays for hadronic correction
    for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
      StPicoBEmcPidTraits *cluster = static_cast<StPicoBEmcPidTraits*>(mPicoDst->bemcPidTraits(iClus));
      if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

      // cluster and tower ID - ID's are calculated as such:
      // mBtowId       = (ntow[0] <= 0 || ntow[0] > 4800) ? -1 : (Short_t)ntow[0];
      // mBtowId23 = (ntow[1] < 0 || ntow[1] >= 9 || ntow[2] < 0 || ntow[2] >= 9) ? -1 : (Char_t)(ntow[1] * 10 + ntow[2]);
      // commented out clusID, towID2, towID3 as they are not used
      //int clusID = cluster->bemcId();  // index in bemc point array - associated BEMC cluster ID
      int towID = cluster->btowId();   // projected track matched tower Id: 1 - 4800 - track matched tower ID
      //int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
      //int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
      if(towID < 0) continue; // THIS IS NEEDED!

      // tower check
      StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towID-1));
      if(!tower) { cout<<"No tower pointer... iTow = "<<towID-1<<endl; continue; }
      if(!AcceptJetTower(tower, towID)) continue; // also eliminates bad towers

      // matched track index
      int trackIndex = cluster->trackIndex();
      if(trackIndex < 0) { continue; }
    
      // get pointer to StPicoTrack: with trackIndex
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(trackIndex));
      if(!trk) { cout<<"No trk pointer...."<<endl; continue; }

      // apply quality cut to matched tracks (as low as 0.2 GeV - no pt cut)
      if(!AcceptTrack(trk, Bfield, mVertex)) continue;

      // tower status set - towerID is matched to track passing quality cuts
      mTowerMatchTrkIndexLast[towID] = trackIndex;
      mTowerMatchTrkIndex[towID][ mTowerStatusArr[towID] ] = trackIndex;
      mTowerStatusArr[towID] = mTowerStatusArr[towID] + 1;
      matchedTowerTrackCounter++;
    } // PIDTraits loop

    // loop over towers and add input vectors to fastjet
    int nTowers = mPicoDst->numberOfBTowHits();
    for(int itow = 0; itow < nTowers; itow++) {
      // get tower pointer
      StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(itow));
      if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

      // tower ID - get from itow shift: which is 1 more than array index
      int towerID = itow + 1;
      if(towerID < 0) continue; // double check these aren't still in the event list

      // tower acceptance cuts - cuts on *bad* towers also
      if(!AcceptJetTower(tower, towerID)) continue;

      // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
      TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
      double towerPhi = towerPosition.Phi();
      if(towerPhi < 0.0)    towerPhi += 2.0*pi;
      if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;
      double towerEta = towerPosition.PseudoRapidity();
      int towerADC = tower->adc();
      double towerEunCorr = tower->energy();  // uncorrected energy
      double LastIndexTrktowerE = tower->energy();  // for use with using single (last indexed - OLD method) match for hadronic correction
      double towerE = tower->energy();        // corrected energy (hadronically - done below)
      double towEtunCorr = towerEunCorr / (1.0*TMath::CosH(towerEta));

      // fill QA histos for jet towers  
      fHistNTowervsID->Fill(towerID);
      fHistNTowervsADC->Fill(towerADC);
      fHistNTowervsE->Fill(towerEunCorr);
      fHistNTowervsEt->Fill(towEtunCorr);
      fHistNTowervsPhi->Fill(towerPhi);
      fHistNTowervsEta->Fill(towerEta);
      fHistNTowervsPhivsEta->Fill(towerPhi, towerEta);

      // fill a NmatchTrack histogram here (centrality dependent)
      if(cbin > 0 && cbin < 6) {
        int arrayBin = cbin - 1; // value -> array
        fHistNMatchTrack[arrayBin]->Fill( (double)mTowerStatusArr[towerID] );
      }

      // cut on min tower energy after filling histos
      if(towerEunCorr < mTowerEnergyMin) continue; // if we don't have enough E to start with, why mess around

      // =======================================================================
      // HADRONIC CORRECTION
      double maxEt = 0.;
      double sumEt = 0.;

      // if tower was is matched to a track or multiple, add up the matched track energies 
      //     (mult opt.) to then subtract from the corresponding tower
      // August 15: if *have* 1+ matched trk-tow AND uncorrected energy of tower is at least your tower constituent cut, then CONTINUE 
      if(mTowerStatusArr[towerID] > 0.5 && towerEunCorr > mTowerEnergyMin) {
        double maxE = 0.0;
        double sumE = 0.0;
// =======================================================================================================================
        // --- last indexed track matched to tower
        // get track pointer
        StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track( mTowerMatchTrkIndexLast[towerID] ));
        if(!trk) { cout<<"No trk pointer...."<<endl; continue; }

        // quality check on track should already be done above
        if(AcceptTrack(trk, Bfield, mVertex)) {
          // momentum vector of track
          TVector3 mTrkMom;
          if(doUsePrimTracks) { 
            // get primary track vector
            mTrkMom = trk->pMom(); 
          } else { 
            // get global track vector
            mTrkMom = trk->gMom(mVertex, Bfield); 
          }

          // track variables
          double p = mTrkMom.Mag();
          double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);

          // apply hadronic correction to tower
          LastIndexTrktowerE = towerEunCorr - (mHadronicCorrFrac * E);
        }  
// =======================================================================================================================
        // --- finds max E track matched to tower *AND* the sum of all matched track E and subtract from said tower
        //     USER provides readMacro.C which method to use for their analysis via SetJetHadCorrType(type);
        // loop over ALL matched tracks
        for(int itrk = 0; itrk < mTowerStatusArr[towerID]; itrk++) {
          StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track( mTowerMatchTrkIndex[towerID][itrk] ));
          if(!trk) { cout<<"No track pointer..."<<endl; continue; }
          if(!AcceptTrack(trk, Bfield, mVertex)) { cout<<"track matched back doesn't pass cuts"<<endl; continue; }  

          // get track variables to matched tower from 3-vector
          TVector3 mTrkMom;
          if(doUsePrimTracks) { 
            // get primary track vector
            mTrkMom = trk->pMom();
          } else { 
            // get global track vector
            mTrkMom = trk->gMom(mVertex, Bfield); 
          }

          // track variables
          double p = mTrkMom.Mag();
          double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
          if(E > maxE) maxE = E;
          sumE = sumE + E;

          // test statement
          //cout<<itrk<<"  itrkE: "<<E<<"  sumtrkE: "<<sumE<<endl;
        } // track loop

        // apply hadronic correction to tower
        maxEt  = (towerEunCorr - (mHadronicCorrFrac * maxE))/(1.0*TMath::CosH(towerEta));
        sumEt  = (towerEunCorr - (mHadronicCorrFrac * sumE))/(1.0*TMath::CosH(towerEta));        

        // fill histogram to compare 2 hadronic correction methods (centrality dependent)
        if(cbin > 0 && cbin < 6) {
          int arrayBin = cbin -1; // value -> array
          fHistHadCorrComparison[arrayBin]->Fill(maxEt, sumEt);
          fHistHadCorrComparisonLast[arrayBin]->Fill(LastIndexTrktowerE / (1.0*TMath::CosH(towerEta)), maxEt);
        }   
//=================================================================================================================
      }  // have a track-tower match
      // else - no match so treat towers on their own. Must meet constituent cut

      // Et - correction comparison
      // test area =============================
      double fMaxEt = (maxEt == 0) ? towerEunCorr / (1.0*TMath::CosH(towerEta)) : maxEt;
      double fSumEt = (sumEt == 0) ? towerEunCorr / (1.0*TMath::CosH(towerEta)) : sumEt;
      double fLastEt = LastIndexTrktowerE / (1.0*TMath::CosH(towerEta));
      //if(mTowerStatusArr[towerID] > 0.5) cout<<"towerEunCorr = "<<towerEunCorr<<"  CosH: "<<1.0*TMath::CosH(towerEta)<<"  LastEt: "<<fLastEt<<"   fMaxEt: "<<fMaxEt<<"   fSumEt: "<<fSumEt<<endl;
      // test area =============================

      // cut on transverse tower energy (more uniform)
      double towerEt = 0.0;
      if(mTowerStatusArr[towerID] < 1) { // no matches, use towers uncorrected energy
        towerEt = towerEunCorr / (1.0*TMath::CosH(towerEta)); 
      } else { 
          if(fJetHadCorrType == StJetFrameworkPicoBase::kLastMatchedTrack)     {  towerEt = fLastEt; towerE = fLastEt * 1.0*TMath::CosH(towerEta); }
          if(fJetHadCorrType == StJetFrameworkPicoBase::kHighestEMatchedTrack) {  towerEt = fMaxEt;  towerE = fMaxEt  * 1.0*TMath::CosH(towerEta); }
          if(fJetHadCorrType == StJetFrameworkPicoBase::kAllMatchedTracks)     {  towerEt = fSumEt;  towerE = fSumEt  * 1.0*TMath::CosH(towerEta); } 
      }
      if(towerEt == 0) { cout<<"fJetHadCorrType - or - towerE actually 0"<<endl; }  // it was unset, because you provided wrong fJetHadCorrType
      if(towerEt < 0) towerEt = 0.0;
      if(towerEt < mTowerEnergyMin) continue;

      // get components from Energy (p - momentum) - the below lines 'converts' the tower energy to momentum:
      TVector3 mom;
      GetMomentum(mom, tower, pi0mass, towerID, towerE);
      double towerPx = mom.x();
      double towerPy = mom.y();
      double towerPz = mom.z();

      // add towers to fastjet - shift tower index (tracks 0+, ghosts = -1, towers < -1)
      int uidTow = -(itow + 2);  
      fjw.AddInputVector(towerPx, towerPy, towerPz, towerE, uidTow); // includes E

      // if running constituent subtractor, do this in parallel
      if(doConstituentSubtr) {
        fastjet::PseudoJet particle_Tower(towerPx, towerPy, towerPz, towerE);
        particle_Tower.set_user_index(uidTow);
        fFull_Event.push_back(particle_Tower);
      }

    } // tower loop

  }   // neutral/full jets

  // run jet finder
  fjw.Run();

}
//
/**
 * This method fills the jet output branch (TClonesArray) with the jet found by the FastJet
 * wrapper. Before filling the jet branch, the utilities are prepared. Then the utilities are
 * called for each jet and finally after jet finding the terminate method of all utilities is called.
 */
void StJetMakerTask::FillJetBranch()
{
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  // loop over FastJet jets
  __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("%d jets found", (Int_t)jets_incl.size()));
  //for(UInt_t ij = 0, jetCount = 0; ij < jets_incl.size(); ++ij) {
  for(UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets,Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

    // PERFORM CUTS ON inclusive JETS before saving
    // cut on min jet pt
    if(jets_incl[ij].perp() < fMinJetPt) continue;
    // cut on min jet area
    if(fjw.GetJetArea(ij) < fMinJetArea*TMath::Pi()*fRadius*fRadius) continue;
    // cut on eta acceptance
    if((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax)) continue;
    // cut on phi acceptance 
    if((jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax)) continue;

    // need to figure out how to get m or E from STAR tracks
    StJet *jet = new ((*fJets)[jetCount])
      StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());

    jet->SetLabel(ij);

    // area vector and components
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());  // same as fjw.GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    // March 8, 2018 - probably don't need this anymore after coming up with method to pass and extract constituents via index!
    // get constituents of jets
    fConstituents = fjw.GetJetConstituents(ij);
    jet->SetJetConstituents(fConstituents);

    // fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    FillJetConstituents(jet, constituents, constituents);

    __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents()));

    jetCount++;
  } // jet loop 

}
//
/**
 * This method is called for each jet. It loops over the jet constituents and
 * adds them to the jet object.
 * @param jet Pointer to the AliEmcalJet object where the jet constituents will be added
 * @param constituents List of the jet constituents returned by the FastJet wrapper
 * @param constituents_unsub List of jet constituents before background subtraction
 * @param flag If kTRUE it means that the argument "constituents" is a list of subtracted constituents
 * @param particles_sub Array containing subtracted constituents - not used
 */
void StJetMakerTask::FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet>& constituents,
    std::vector<fastjet::PseudoJet>& constituents_unsub, Int_t flag, TString particlesSubName)
{
  // initialize some variables/counters
  Double_t neutralE = 0, maxTrack = 0, maxTower = 0;
  Int_t nt = 0; // track counter
  Int_t nc = 0; // cluster counter
  Int_t ng = 0; // ghost counter  
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  // initially set track and cluster constituent sizes
  jet->SetNumberOfTracks(constituents.size());
  jet->SetNumberOfClusters(constituents.size());

  // loop over constituents for ij'th jet
  for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
    // get user defined index
    Int_t uid = constituents[ic].user_index();

    // CHARGED COMPONENT (tracks)
    if(uid >= 0) {
      jet->AddTrackAt(uid, nt);

      // get jet track pointer
      StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track(uid));
      if(!trk) continue;

      // acceptance and kinematic quality cuts - probably not needed (done before passing to fastjet)
      if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

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
      double pt = mTrkMom.Perp();
      double phi = mTrkMom.Phi();
      double eta = mTrkMom.PseudoRapidity();

      // adjust phi value:  0 < phi < 2pi
      if(phi < 0.0)    phi += 2.0*pi;
      if(phi > 2.0*pi) phi -= 2.0*pi;

      // find max track pt
      if(pt > maxTrack) maxTrack = pt;

      // fill some QA histograms
      fHistJetNTrackvsPt->Fill(pt);
      fHistJetNTrackvsPhi->Fill(phi);
      fHistJetNTrackvsEta->Fill(eta);
      fHistJetNTrackvsPhivsEta->Fill(phi, eta);

      // increase track counter
      nt++;
    } else { // uid < 0

      // TODO - finish this if we ever care to look at ghosts.. (not a priority)
      // ghosts
      if(uid == -1) {
        //if(fFillGhost) jet->AddGhost(constituents[ic].px(), constituents[ic].py(), constituents[ic].pz(), constituents[ic].e());
        ng++;
      }

      // neutral towers: start at (index = -2, ghosts are -1, and tracks 0+)
      if(uid < -1) {
        // inialize Et variables
        double maxEt = 0.0;
        double sumEt = 0.0;

        // convert uid to tower index (index of tower - in BTowHit array)
        Int_t towIndex = -(uid + 2);   // 1 less than towerID
        jet->AddClusterAt(towIndex, nc);

        // get jet tower pointer
        StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(towIndex));
        if(!tower) continue;

        // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
        // tower ID: get from tower index shift: +1 above tower index in the BTowHit array
        int towerID = towIndex + 1;
        if(towerID < 0) continue;

        // tower acceptance cuts, also cuts bad towers - already done before input to FJ
        if(!AcceptJetTower(tower, towerID)) continue;

        // get vector for tower position
        TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
        double towerPhi = towerPosition.Phi();
        if(towerPhi < 0.0)    towerPhi += 2.0*pi;
        if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;
        double towerEta = towerPosition.PseudoRapidity();
        double towEuncorr = tower->energy();
        double LastIndexTrktowE = tower->energy();
        double towE = tower->energy();
        int towADC = tower->adc();

        // cut on min tower energy - shouldn't matter as this tower was *saved* to the jet 
        if(towEuncorr < mTowerEnergyMin) continue; // if we don't have enough E to start with, why mess around

        // need to perform hadronic correction again since StBTowHit object is not updated
        // if tower was not matched to an accepted track, use it for jet by itself if > constit cut
        if(mTowerStatusArr[towerID] > 0.5) {
          double maxE = 0.0;
          double sumE = 0.0;

// =======================================================================================================================
          // --- last indexed track matched to tower
          // get last matched track pointer
          StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track( mTowerMatchTrkIndexLast[towerID] ));
          if(!track) { cout<<"No track pointer...."<<endl; continue; }

          // want to check if pass cuts, but if not then don't correct tower, but use it
          if(AcceptTrack(track, Bfield, mVertex)) {
            // get track momentum vector
            TVector3 mTrkMom;
            if(doUsePrimTracks) {
              // get primary track vector
              mTrkMom = track->pMom();
            } else {
              // get global track vector
              mTrkMom = track->gMom(mVertex, Bfield);
            }

            // track variables
            double p = mTrkMom.Mag();
            double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);

            // apply hadronic correction to tower
            LastIndexTrktowE = towEuncorr - (mHadronicCorrFrac * E);
          }

// ============================================================================================================================
          // loop over ALL matched tracks
          for(int itrk = 0; itrk < mTowerStatusArr[towerID]; itrk++) {
            StPicoTrack *trk = static_cast<StPicoTrack*>(mPicoDst->track( mTowerMatchTrkIndex[towerID][itrk] ));
            if(!trk) { cout<<"No track pointer..."<<endl; continue; }
            if(!AcceptTrack(trk, Bfield, mVertex)) { 
              cout<<"track matched back doesn't pass cuts"<<endl; continue; 
            }

            // get track variables to matched tower from 3-vector
            TVector3 mTrkMom;
            if(doUsePrimTracks) {
              // get primary track vector
              mTrkMom = trk->pMom();
            } else {
              // get global track vector
              mTrkMom = trk->gMom(mVertex, Bfield);
            }

            // track variables
            double p = mTrkMom.Mag();
            double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);

            if(E > maxE) maxE = E;
            sumE = sumE + E;
          } // track loop

          // apply hadronic correction to tower
          maxEt  = (towEuncorr - (mHadronicCorrFrac * maxE))/(1.0*TMath::CosH(towerEta));
          sumEt  = (towEuncorr - (mHadronicCorrFrac * sumE))/(1.0*TMath::CosH(towerEta));
// ============================================================================================================================
        } // hadronic correction - have matched track
        // else - no match so treat towers on their own

        // set tower transverse energy, various options available
        double fMaxEt = (maxEt == 0) ? towEuncorr / (1.0*TMath::CosH(towerEta)) : maxEt;
        double fSumEt = (sumEt == 0) ? towEuncorr / (1.0*TMath::CosH(towerEta)) : sumEt;
        double fLastEt = LastIndexTrktowE / (1.0*TMath::CosH(towerEta));

        // cut on transverse tower energy (more uniform)
        double towEt = 0.0;
        if(mTowerStatusArr[towerID] < 1) { // no matches, use towers uncorrected energy
          towEt = towEuncorr / (1.0*TMath::CosH(towerEta));
        } else {
          if(fJetHadCorrType == StJetFrameworkPicoBase::kLastMatchedTrack)     { towEt = fLastEt; towE = fLastEt * 1.0*TMath::CosH(towerEta); }
          if(fJetHadCorrType == StJetFrameworkPicoBase::kHighestEMatchedTrack) { towEt = fMaxEt;  towE = fMaxEt  * 1.0*TMath::CosH(towerEta); }
          if(fJetHadCorrType == StJetFrameworkPicoBase::kAllMatchedTracks)     { towEt = fSumEt;  towE = fSumEt  * 1.0*TMath::CosH(towerEta); }
        }
        if(towEt == 0) { cout<<"fJetHadCorrType - or - towerE actually 0"<<endl; }  // it was unset, because you provided wrong fJetHadCorrType
        if(towEt < 0.0) towEt = 0.0;
        if(towEt < mTowerEnergyMin) continue;

        // find max tower Et and neutral E total
        if(towEt > maxTower) maxTower = towEt;
        neutralE += towE;

        // fill QA histos for jet towers
        fHistJetNTowervsID->Fill(towerID);
        fHistJetNTowervsADC->Fill(towADC);
        fHistJetNTowervsE->Fill(towE);
        fHistJetNTowervsEt->Fill(towEt);
        fHistJetNTowervsPhi->Fill(towerPhi);
        fHistJetNTowervsEta->Fill(towerEta);
        fHistJetNTowervsPhivsEta->Fill(towerPhi, towerEta);
        fHistQATowIDvsEta->Fill(towerID, towerEta);
        fHistQATowIDvsPhi->Fill(towerID, towerPhi);

        // increase tower counter
        nc++;
      } // towers

    } // uid < 0

  }  // end of constituent loop

  // set some jet properties
  jet->SetNumberOfTracks(nt);
  jet->SetNumberOfClusters(nc);
  jet->SetMaxTrackPt(maxTrack);
  jet->SetMaxTowerEt(maxTower);  
  jet->SetNEF(neutralE/jet->E());  // should this be Et? FIXME
  //jet->SortConstituents(); // TODO see how this works - sorts ClusterIds() and TrackIds() by index (increasing)

  // fill jets histograms
  fHistNJetsvsPt->Fill(jet->Pt()); 
  fHistNJetsvsPhi->Fill(jet->Phi());
  fHistNJetsvsEta->Fill(jet->Eta());
  fHistNJetsvsPhivsEta->Fill(jet->Phi(), jet->Eta());
  fHistNJetsvsArea->Fill(jet->Area());  // same as fjw.GetJetArea(ij)
  fHistNJetsvsNConstituents->Fill(nt + nc);
  fHistNJetsvsNTracks->Fill(nt);
  fHistNJetsvsNTowers->Fill(nc);

}
/**
 * Sorts jets by pT (decreasing)
 * @param[out] indexes This array is used to return the indexes of the jets ordered by pT
 * @param[in] array Vector containing the list of jets obtained by the FastJet wrapper
 * @return kTRUE if at least one jet was found in array; kFALSE otherwise
 */
Bool_t StJetMakerTask::GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const
{
  static Float_t pt[9999] = {0};
  const Int_t n = (Int_t)array.size();
  if(n < 1) return kFALSE;

  for(Int_t i = 0; i < n; i++)
    pt[i] = array[i].perp();

  TMath::Sort(n, pt, indexes);

  return kTRUE;
}
/**
* An instance of this class can be "locked". Once locked, it cannot be unlocked.
 * If the instance is locked, attempting to change the configuration will throw a
 * fatal and stop the execution of the program. This method checks whether the instance
 * is locked and throw a fatal if it is locked.
 */
Bool_t StJetMakerTask::IsLocked() const
{
  if (fLocked) {
    Form("Jet finder task is locked! Changing properties is not allowed."); 
    return kStFatal;
  } 
  else {
    return kStOK;
  }
}
/**
 * Converts the internal enum values representing jet algorithms in
 * the corresponding values accepted by the FastJet wrapper.
 * @param algo Algorithm represented in the EJetAlgo_t enum
 * @return Algortithm represented in the fastjet::JetAlgorithm enum
 */
/*
fastjet::JetAlgorithm StJetMakerTask::ConvertToFJAlgo(EJetAlgo_t algo)
{
  switch(algo) {
  case StJetFrameworkPicoBase::kt_algorithm:
    return fastjet::kt_algorithm;
  case StJetFrameworkPicoBase::antikt_algorithm:
    return fastjet::antikt_algorithm;
  case StJetFrameworkPicoBase::cambridge_algorithm:
    return fastjet::cambridge_algorithm;
  case StJetFrameworkPicoBase::genkt_algorithm:
    return fastjet::genkt_algorithm;
  case StJetFrameworkPicoBase::cambridge_for_passive_algorithm:
    return fastjet::cambridge_for_passive_algorithm;
  case StJetFrameworkPicoBase::genkt_for_passive_algorithm:
    return fastjet::genkt_for_passive_algorithm;
  case StJetFrameworkPicoBase::plugin_algorithm:
    return fastjet::plugin_algorithm;
  case StJetFrameworkPicoBase::undefined_jet_algorithm:
    return fastjet::undefined_jet_algorithm;

  default:
    ::Error("StJetMakerTask::ConvertToFJAlgo", "Jet algorithm %d not recognized!!!", algo);
    return fastjet::undefined_jet_algorithm;
  }
}
*/
/**
 * Converts the internal enum values representing jet recombination schemes in
 * the corresponding values accepted by the FastJet wrapper.
 * @param reco Recombination scheme represented in the EJetAlgo_t enum
 * @return Recombination scheme represented in the fastjet::JetAlgorithm enum
 */
/*
fastjet::RecombinationScheme StJetMakerTask::ConvertToFJRecoScheme(ERecoScheme_t reco)
{
  switch(reco) {
  case StJetFrameworkPicoBase::E_scheme:
    return fastjet::E_scheme;
  case StJetFrameworkPicoBase::pt_scheme:
    return fastjet::pt_scheme;
  case StJetFrameworkPicoBase::pt2_scheme:
    return fastjet::pt2_scheme;
  case StJetFrameworkPicoBase::Et_scheme:
    return fastjet::Et_scheme;
  case StJetFrameworkPicoBase::Et2_scheme:
    return fastjet::Et2_scheme;
  case StJetFrameworkPicoBase::BIpt_scheme:
    return fastjet::BIpt_scheme;
  case StJetFrameworkPicoBase::BIpt2_scheme:
    return fastjet::BIpt2_scheme;
  case StJetFrameworkPicoBase::external_scheme:
    return fastjet::external_scheme;

  default:
    ::Error("StJetMakerTask::ConvertToFJRecoScheme", "Recombination scheme %d not recognized!!!", reco);
    return fastjet::external_scheme;
  }
}
*/
//
// Function: track quality cuts
//________________________________________________________________________
Bool_t StJetMakerTask::AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert) {
  // constants: assume neutral pion mass
  ///double pi0mass = Pico::mMass[0]; // GeV
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
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // jet track acceptance cuts now - note difference from AcceptJetTrack()
  if((eta < fTrackEtaMin) || (eta > fTrackEtaMax)) return kFALSE;
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  if((phi < fTrackPhiMin) || (phi > fTrackPhiMax)) return kFALSE;

  // additional quality cuts for tracks
  if(dca > fJetTrackDCAcut)            return kFALSE;
  if(nHitsFit < fJetTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fJetTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
//
// Function: jet track quality cuts
// - this function should only be used for jet constituent tracks
// 	NOT when considering track-tower matches  (Aug 19')
//________________________________________________________________________
Bool_t StJetMakerTask::AcceptJetTrack(StPicoTrack *trk, Float_t B, TVector3 Vert) {
  // constants: assume neutral pion mass
  ///double pi0mass = Pico::mMass[0]; // GeV
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
  //double p = mTrkMom.Mag();
  //double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
  //short charge = trk->charge();
  //double dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // jet track acceptance cuts now
  if(pt < fMinJetTrackPt) return kFALSE;
  if(pt > fMaxJetTrackPt) return kFALSE; // 20.0 STAR, 100.0 ALICE
  if((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax)) return kFALSE;
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  if((phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) return kFALSE;
      
  // additional quality cuts for tracks
  if(dca > fJetTrackDCAcut)            return kFALSE;
  if(nHitsFit < fJetTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fJetTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
//
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StJetMakerTask::AcceptJetTower(StPicoBTowHit *tower, Int_t towerID) {
  // constants:
  double pi = 1.0*TMath::Pi();

  // tower ID - passed into function: make sure some of these aren't still in event array
  if(towerID < 0) return kFALSE; 

  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  double phi = towerPosition.Phi();
  if(phi < 0.0)    phi += 2.0*pi;
  if(phi > 2.0*pi) phi -= 2.0*pi;
  double eta = towerPosition.PseudoRapidity();

  // check for bad (and dead) towers
  bool TowerOK = mBaseMaker->IsTowerOK(towerID);      // kTRUE means GOOD
  bool TowerDead = mBaseMaker->IsTowerDead(towerID);  // kTRUE means BAD
  if(!TowerOK)  { return kFALSE; }
  if(TowerDead) { return kFALSE; }

  // jet track acceptance cuts njow
  if((eta < fJetTowerEtaMin) || (eta > fJetTowerEtaMax)) return kFALSE;
  if((phi < fJetTowerPhiMin) || (phi > fJetTowerPhiMax)) return kFALSE;

  // passed all above cuts - keep tower and fill input vector to fastjet
  return kTRUE;
}
//
// Function: get cent bin on 0-8 or 0-15 scale
//________________________________________________________________________
Int_t StJetMakerTask::GetCentBin(Int_t cent, Int_t nBin) const
{
  Int_t centbin = -1;

  if(nBin == 16) { centbin = nBin - 1 - cent; }
  if(nBin == 9)  { centbin = nBin - 1 - cent; }

  return centbin;
}
//
// Function: calculate momentum of a tower
//________________________________________________________________________________________________________
Bool_t StJetMakerTask::GetMomentum(TVector3 &mom, const StPicoBTowHit* tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const {
  // vertex components - only need if below method is used
  // mGeom3->getEtaPhi(towerID,tEta,tPhi);
  //double xVtx = mVertex.x();
  //double yVtx = mVertex.y();
  //double zVtx = mVertex.z();

  // get mass, E, P, ID
  if(mass < 0) mass = 0.0;
  //Double_t energy = tower->energy();
  Double_t energy = CorrectedEnergy;   // USE THIS, to use the corrected energy to get the momentum components
  Double_t p = TMath::Sqrt(energy*energy - mass*mass);

  // tower ID - passed into function
  // get tower position
  TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
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
  } else { return kFALSE; }

  return kTRUE;
}
//
//
//_________________________________________________________________________
void StJetMakerTask::FillEmcTriggersArr() {
  // zero out trigger array and get number of Emcal Triggers
  for(int i = 0; i < 8; i++) { fEmcTriggerArr[i] = kFALSE; }
  int nEmcTrigger = mPicoDst->numberOfEmcTriggers();

  // tower - HT trigger types arrays: zero these out - so they are refreshed for each event
  for(int i = 0; i < 4801; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }

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

    // fill for valid triggers
    if(isHT1) fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = kTRUE;
    if(isHT2) fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = kTRUE;
    if(isHT3) fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = kTRUE;
  }

}
//========
/*
      // get tower location - from ID: need to correct eta (when using this method)
      float tEta, tPhi;    // need to be floats
      StEmcGeom *mGeom3 = (StEmcGeom::instance("bemc"));
      double radius = mGeom3->Radius();
      mGeom3->getEtaPhi(towerID,tEta,tPhi);
      if(tPhi < 0.0)    tPhi += 2.0*pi;
      if(tPhi > 2.0*pi) tPhi -= 2.0*pi;

     // correct eta for Vz position 
     double theta;
     theta = 2 * atan(exp(-tEta)); // getting theta from eta 
     double z = 0;
     if(tEta != 0) z = radius / tan(theta);  // 231 cm = radius of SMD 
     double zNominal = z - mVertex.z();
     double thetaCorr = atan2(radius, zNominal); // theta with respect to primary vertex
     double etaCorr =-log(tan(thetaCorr / 2.0)); // eta with respect to primary vertex 

     cout<<"tEta = "<<tEta<<"  etaCorr = "<<etaCorr<<"  towerEta = "<<towerEta<<"  tPhi = "<<tPhi<<"  towerPhi = "<<towerPhi<<endl;
*/
//==========
//
// Function: Returns pt of hardest track in the event
//_________________________________________________________________________________________________
Double_t StJetMakerTask::GetMaxTrackPt()
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

    // primary track switch: get momentum vector of track - global or primary track
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
Double_t StJetMakerTask::GetMaxTowerEt()
{
  // get # of towers
  int nTowers = mPicoDst->numberOfBTowHits();
  double fMaxTowerEt = -99;

  // loop over all towers
  for(int i = 0; i < nTowers; i++) {
    // get tower pointer
    StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(i));
    if(!tower) continue;

    int towerID = i+1;
    if(towerID < 0) continue;

    // tower position - from vertex and ID: shouldn't need additional eta correction
    TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
    //double towerPhi = towerPosition.Phi();
    //if(towerPhi < 0.0)    towerPhi += 2.0*pi;
    //if(towerPhi > 2.0*pi) towerPhi -= 2.0*pi;
    double towerEta = towerPosition.PseudoRapidity();
    double towerEuncorr = tower->energy();               // uncorrected energy
    double towEt = towerEuncorr / (1.0*TMath::CosH(towerEta)); // should this be used instead?

    // get max tower
    if(towEt > fMaxTowerEt) { fMaxTowerEt = towEt; }
  }

  return fMaxTowerEt;
}
//
/**
 * This method fills the jet output branch (TClonesArray) with the jet found by the FastJet wrapper.
 * This is for constituent subtractor performed to jets
 */
void StJetMakerTask::FillJetBGBranch()
{
   // E_scheme, pt_scheme, pt2_scheme, Et_scheme, Et2_scheme, BIpt_scheme, BIpt2_scheme, WTA_pt_scheme, WTA_modp_scheme
   fastjet::RecombinationScheme    recombScheme;
   if (fRecombScheme == 0)     recombScheme = fastjet::E_scheme;
   if (fRecombScheme == 1)     recombScheme = fastjet::pt_scheme;
   if (fRecombScheme == 2)     recombScheme = fastjet::pt2_scheme;
   if (fRecombScheme == 3)     recombScheme = fastjet::Et_scheme;
   if (fRecombScheme == 4)     recombScheme = fastjet::Et2_scheme;
   if (fRecombScheme == 5)     recombScheme = fastjet::BIpt_scheme;
   if (fRecombScheme == 6)     recombScheme = fastjet::BIpt2_scheme;

   // jet algorithm
   fastjet::JetAlgorithm          algorithm;
   if (fJetAlgo == 1)      algorithm = fastjet::antikt_algorithm;
   if (fJetAlgo == 0)      algorithm = fastjet::kt_algorithm;
   if (fJetAlgo == 2)      algorithm = fastjet::cambridge_algorithm;
   if (fJetAlgo == 11)     algorithm = fastjet::cambridge_for_passive_algorithm;
   fastjet::Strategy              strategy = fastjet::Best;

   double jetAbsRapMax = 1.0;

   // create a jet definition for the clustering: We use the anti-kt algorithm with a radius of 0.5
   //----------------------------------------------------------
   fastjet::JetDefinition jet_def(algorithm, fRadius, recombScheme, strategy);

   // create an area definition for the clustering
   // ----------------------------------------------------------
   // ghosts should go up to the acceptance of the detector or (with infinite acceptance) at least 2R beyond the region where you plan to investigate jets.
   double ghost_maxrap = 1.2;
   fastjet::GhostedAreaSpec area_spec(ghost_maxrap, 1, fGhostArea);
   fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

   // run the jet clustering with the above jet and area definitions for both the hard and full event
   //
   // We retrieve the jets above 7 GeV in both case (note that the 7-GeV cut we be applied again later on after we subtract the jets from the full event)
   // ----------------------------------------------------------
   fastjet::ClusterSequenceArea clust_seq_full(fFull_Event, jet_def, area_def);

   // minimum jet pt for inclusive jets 
   double ptmin = fMinJetPt;
   //Selector sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(3.0);
   vector<fastjet::PseudoJet> full_jets = sorted_by_pt(clust_seq_full.inclusive_jets(ptmin));
   //vector<fastjet::PseudoJet> full_jets = fjw.GetInclusiveJets();

   // Now turn to the estimation of the background (for the full event)
   //
   // There are different ways to do that. In general, this also requires clustering the particles that will be handled internally in FastJet. 
   //
   // The suggested way to proceed is to use a BackgroundEstimator constructed from the following 3 arguments:
   //  - a jet definition used to cluster the particles.
   //    . We strongly recommend using the kt or Cambridge/Aachen algorithm (a warning will be issued otherwise)
   //    . The choice of the radius is a bit more subtle. R=0.4 has been chosen to limit the impact of hard jets; in samples of
   //      dominantly sparse events it may cause the UE/pileup to be underestimated a little, a slightly larger value (0.5 or 0.6) may be better.
   //  - An area definition for which we recommend the use of explicit ghosts (i.e. active_area_explicit_ghosts) As mentionned in the area example (06-area.cc), ghosts should
   //    extend sufficiently far in rapidity to cover the jets used in the computation of the background (see also the comment below)
   //  - A Selector specifying the range over which we will keep the jets entering the estimation of the background (you should
   //    thus make sure the ghosts extend far enough in rapidity to cover the range, a warning will be issued otherwise).
   //    In this particular example, the two hardest jets in the event are removed from the background estimation
   // ----------------------------------------------------------

   // create what we need for the background estimation
   //----------------------------------------------------------
   fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, fRadius, recombScheme, strategy);
   //fastjet::Selector selector = fastjet::SelectorAbsRapMax(jetAbsRapMax) * (!fastjet::SelectorNHardest(2));
   fastjet::Selector rho_range =  fastjet::SelectorAbsRapMax(3.0); // 3.0
   fastjet::ClusterSequenceArea clust_seq_rho(fFull_Event, jet_def, area_def); // not used FIXME 

   fastjet::JetMedianBackgroundEstimator bge_rho(rho_range, jet_def_for_rho, area_def);
   // TODO next 2 lines commented out to suppress warnings, doesn't affect results - Sept26, 2018
   //fastjet::BackgroundJetScalarPtDensity *scalarPtDensity = new fastjet::BackgroundJetScalarPtDensity();
   //bge_rho.set_jet_density_class(scalarPtDensity); // this changes computation of pt of patches from vector sum to scalar sum. Theor., the scalar sum seems more reasonable.
   bge_rho.set_particles(fFull_Event);

   // subtractor:
   //----------------------------------------------------------
   fastjet::contrib::ConstituentSubtractor subtractor(&bge_rho);
   subtractor.set_common_bge_for_rho_and_rhom(true); // TODO - need this to omit warning, same results Sept26, 2018

   // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
   //subtractor.use_common_bge_for_rho_and_rhom(true); // for massless input particles it does not make any difference (rho_m is always zero)
   ////cout << subtractor.description() << endl;
   ////cout << "  Giving, for the full event" << "    rho     = " << bge_rho.rho() << "    sigma   = " << bge_rho.sigma() << endl;

   // fill histogram with FastJet calculated rho
   fHistFJRho->Fill(bge_rho.rho());

   std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
   // sort jets according to jet pt
   static Int_t indexes[9999] = {-1};
   GetSortedArray(indexes, jets_incl);

   // ===================================
   // clear constituent array
   fConstituents.clear();

   // loop over FastJet jets
   __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("%d jets found", (Int_t)jets_incl.size()));
   for(UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet) {
     Int_t ij = indexes[ijet];
     __DEBUG(StJetFrameworkPicoBase::kDebugFillJets,Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)));

     // PERFORM CUTS ON inclusive JETS before saving
     // cut on min jet pt
     if(jets_incl[ij].perp() < fMinJetPt) continue;
     // cut on min jet area
     if(fjw.GetJetArea(ij) < fMinJetArea*TMath::Pi()*fRadius*fRadius) continue;
     // cut on eta acceptance
     if((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax)) continue;
     // cut on phi acceptance 
     if((jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax)) continue;

     // apply subtractor here
     fastjet::PseudoJet subtracted_jet = subtractor(jets_incl[ij]);

     // need to figure out how to get m or E from STAR tracks - FIXME - may want to just make this 'fJets'
     StJet *jet = new ((*fJetsBGsub)[jetCount])
       StJet(subtracted_jet.perp(), subtracted_jet.eta(), subtracted_jet.phi(), subtracted_jet.m());

     // set label
     jet->SetLabel(ij);

     // area vector and components
     fastjet::PseudoJet area(fjw.GetJetAreaVector(ij)); // FIXME - this MAY not be correct after subtraction, double check
     jet->SetArea(area.perp());  // same as fjw.GetJetArea(ij)
     jet->SetAreaEta(area.eta());
     jet->SetAreaPhi(area.phi());
     jet->SetAreaE(area.E());

     // get constituents of jets - these should have the background (fraction) subtracted
     //fConstituents = fjw.GetJetConstituents(ij);
     //jet->SetJetConstituents(fConstituents);

     // fill jet constituents - these are identified by their index
     vector<fastjet::PseudoJet> constituents = subtracted_jet.constituents();
     jet->SetJetConstituents(constituents);  // constituents are the corrected pseudojet objects
     FillJetConstituents(jet, constituents, constituents);

     __DEBUG(StJetFrameworkPicoBase::kDebugFillJets, Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents()));

     jetCount++;
   } // jet loop 

}
//
// sets errors on histograms up before filling
// set sum weights
//________________________________________________________________________
void StJetMakerTask::SetSumw2() {
  fHistMultiplicity->Sumw2();
  fHistCentrality->Sumw2();
  fHistCentralityPostCut->Sumw2();
  fHistFJRho->Sumw2();
  //fProfEventBBCx->Sumw2();
  //fProfEventZDCx->Sumw2();

  fHistNTrackvsPt->Sumw2();
  fHistNTrackvsPhi->Sumw2();
  fHistNTrackvsEta->Sumw2();
  fHistNTrackvsPhivsEta->Sumw2();
  fHistNTowervsID->Sumw2();
  fHistNTowervsADC->Sumw2();
  fHistNTowervsE->Sumw2();
  fHistNTowervsEt->Sumw2();
  fHistNTowervsPhi->Sumw2();
  fHistNTowervsEta->Sumw2();
  fHistNTowervsPhivsEta->Sumw2();

  for(int i = 0; i < 5; i++) { 
    fHistNMatchTrack[i]->Sumw2();
    fHistHadCorrComparison[i]->Sumw2();
    fHistHadCorrComparisonLast[i]->Sumw2();
  }
  fHistJetNTrackvsPt->Sumw2();
  fHistJetNTrackvsPhi->Sumw2();
  fHistJetNTrackvsEta->Sumw2();
  fHistJetNTrackvsPhivsEta->Sumw2();
  fHistJetNTowervsID->Sumw2();
  fHistJetNTowervsADC->Sumw2();
  fHistJetNTowervsE->Sumw2();
  fHistJetNTowervsEt->Sumw2();
  fHistJetNTowervsPhi->Sumw2();
  fHistJetNTowervsEta->Sumw2();
  fHistJetNTowervsPhivsEta->Sumw2();

  fHistNJetsvsPt->Sumw2();
  fHistNJetsvsPhi->Sumw2();
  fHistNJetsvsEta->Sumw2();
  fHistNJetsvsPhivsEta->Sumw2();
  fHistNJetsvsArea->Sumw2();
  fHistNJetsvsNConstituents->Sumw2();
  fHistNJetsvsNTracks->Sumw2();
  fHistNJetsvsNTowers->Sumw2();

  fHistQATowIDvsEta->Sumw2();
  fHistQATowIDvsPhi->Sumw2();
}
//
// function to do some event QA
//________________________________________________________________________________________
void StJetMakerTask::RunEventQA() {
  // get run ID, transform to run order for filling histogram of corrections
  int fRunId = mPicoEvent->runId();
  int RunId_Order = GetRunNo(fRunId);// + 1;
  //if(RunId_Order < -1) return kStOK;

  // vertex
  //TVector3 fPrimaryVertex = mPicoEvent->primaryVertex();
  //float fZVtx = fPrimaryVertex.z();
  //float fVzVPD = mPicoEvent->vzVpd();

  // coincidence rates
  float fZDCx = mPicoEvent->ZDCx();
  float fBBCx = mPicoEvent->BBCx();

  // fill event QA histograms
  fProfEventBBCx->Fill(RunId_Order + 1., fBBCx);
  fProfEventZDCx->Fill(RunId_Order + 1., fZDCx);

}
//
// this function checks for the bin number of the run from a runlist header 
// in order to apply various corrections and fill run-dependent histograms
// 1287 - Liang
// _________________________________________________________________________________
Int_t StJetMakerTask::GetRunNo(int runid){
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
  // 1654 for Run14 AuAu (OLD), new picoDst production is 830
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
