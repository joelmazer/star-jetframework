// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################
// $Id$

#include "StPicoTrackClusterQA.h"

// root classes
#include <TChain.h>
#include <TClonesArray.h>
#include "TH1F.h"
#include "TH2F.h"
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include "TFile.h"

// general StRoot classes
#include "StThreeVectorF.hh"

// StRoot classes
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoDstMaker/StPicoArrays.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StPicoConstants.h"

// test for clusters:
#include "StMuDSTMaker/COMMON/StMuTrack.h"
// StEmc:
#include "StEmcClusterCollection.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcPoint.h"
//#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcUtil/projection/StEmcPosition.h"

// extra includes
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StJetFrameworkPicoBase.h"

// centrality
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StPicoTrackClusterQA)

//________________________________________________________________________
StPicoTrackClusterQA::StPicoTrackClusterQA() : 
//  StJetFrameworkPicoBase(),
  StMaker(),
  doWriteHistos(kFALSE),
  doUsePrimTracks(kFALSE), 
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum
  fDoEffCorr(kFALSE),
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  mOutName(""),
  fTracksName(""),
  fCaloName(""),
  fTrackPtMinCut(0.2),
  fTrackPtMaxCut(20.0),
  fClusterPtMinCut(0.2),
  fClusterPtMaxCut(1000.0),
  fTrackPhiMinCut(0.0),
  fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0),
  fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15),
  fTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fTriggerEventType(0),
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
  mu(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  grefmultCorr(0x0)
{
  // Default constructor.
  for(int i=0; i<7; i++) { fEmcTriggerArr[i] = 0; }

}

//________________________________________________________________________
StPicoTrackClusterQA::StPicoTrackClusterQA(const char *name, bool doHistos = kFALSE, const char* outName = "") : 
//  StJetFrameworkPicoBase(name),
  StMaker(name),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum
  fDoEffCorr(kFALSE),
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  mOutName(outName),
  fTracksName("Tracks"),
  fCaloName("Clusters"),
  fTrackPtMinCut(0.2), //0.20
  fTrackPtMaxCut(20.0), 
  fClusterPtMinCut(0.2),
  fClusterPtMaxCut(1000.0),
  fTrackPhiMinCut(0.0),
  fTrackPhiMaxCut(2.0*TMath::Pi()),
  fTrackEtaMinCut(-1.0), 
  fTrackEtaMaxCut(1.0),
  fTrackDCAcut(3.0),
  fTracknHitsFit(15),
  fTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fCentralityScaled(0.),
  ref16(-99), ref9(-99),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  fRunNumber(0),
  fTriggerEventType(0),
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
  mu(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  grefmultCorr(0x0)
{
  // Standard constructor.
  if (!name) return;
  SetName(name);

  for(int i=0; i<7; i++) { fEmcTriggerArr[i] = 0; }
}

//________________________________________________________________________
StPicoTrackClusterQA::~StPicoTrackClusterQA()
{
  // Destructor
  if(fHistNTrackvsPt)   delete fHistNTrackvsPt;
  if(fHistNTrackvsPhi)  delete fHistNTrackvsPhi;
  if(fHistNTrackvsEta)  delete fHistNTrackvsEta;
  if(fHistNTrackvsPhivsEta) delete fHistNTrackvsPhivsEta;
  if(fHistNTowervsE)    delete fHistNTowervsE;
  if(fHistNTowervsPhi)  delete fHistNTowervsPhi;
  if(fHistNTowervsEta)  delete fHistNTowervsEta;
  if(fHistNTowervsPhivsEta) delete fHistNTowervsPhivsEta;

  delete fHistEventSelectionQA;
  delete fHistEventSelectionQAafterCuts;
  delete hTriggerIds;
  delete hEmcTriggers;
}

//-----------------------------------------------------------------------------
Int_t StPicoTrackClusterQA::Init() {
  DeclareHistograms();

  // may not need, used for old RUNS
  // StRefMultCorr* getgRefMultCorr() ; // For grefmult //Run14 AuAu200GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr(); }
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    if(fCentralityDef == StJetFrameworkPicoBase::kgrefmult) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr(); }
    if(fCentralityDef == StJetFrameworkPicoBase::kgrefmult_P16id) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id(); }
    if(fCentralityDef == StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx(); }
    if(fCentralityDef == StJetFrameworkPicoBase::kgrefmult_VpdMB30) { grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_VpdMB30(); }
  }

  return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StPicoTrackClusterQA::Finish() {
  //  Summarize the run.
/*
  cout<<"StPicoTrackClusterQA::Finish()\n";
  //cout << "\tProcessed " << mEventCounter << " events." << endl;
  //cout << "\tWithout PV cuts: "<< mAllPVEventCounter << " events"<<endl;
  //cout << "\tInput events: "<<mInputEventCounter<<endl;
*/

  if(doWriteHistos && mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "RECREATE");
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

//________________________________________________________________________
void StPicoTrackClusterQA::DeclareHistograms() {
    // declare histograms
    double pi = 1.0*TMath::Pi();
    fHistNTrackvsPt = new TH1F("fHistNTrackvsPt", "Ntracks vs p_{T}", 100, 0., 20.);
    fHistNTrackvsPhi = new TH1F("fHistNTrackvsPhi", "Ntracks vs #phi", 72, 0., 2*pi);
    fHistNTrackvsEta = new TH1F("fHistNTrackvsEta", "Ntracks vs #eta", 40, -1.0, 1.0);
    fHistNTrackvsPhivsEta = new TH2F("fHistNTrackvsPhivsEta", "Ntrackvs #phi vs #eta", 144, 0, 2*pi, 20, -1.0, 1.0);
    fHistNTowervsE = new TH1F("fHistNTowervsE", "Ntowers vs energy", 100, 0., 20.0);
    fHistNTowervsPhi = new TH1F("fHistNTowervsPhi", "Ntowers vs #phi", 144, 0., 2*pi);
    fHistNTowervsEta = new TH1F("fHistNTowervsEta", "Ntowers vs #eta", 40, -1.0, 1.0);
    fHistNTowervsPhivsEta = new TH2F("fHistNTowervsPhivsEta", "Ntowers vs #phi vs #eta", 144, 0, 2*pi, 20, -1.0, 1.0);

    // Event Selection QA histo
    fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
    fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
    hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
    hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);

}

//________________________________________________________________________

void StPicoTrackClusterQA::WriteHistograms() {
  // write histograms
  fHistNTrackvsPt->Write();
  fHistNTrackvsPhi->Write();
  fHistNTrackvsEta->Write();
  fHistNTrackvsPhivsEta->Write();
  fHistNTowervsE->Write();
  fHistNTowervsPhi->Write();
  fHistNTowervsEta->Write();
  fHistNTowervsPhivsEta->Write();

  // QA histos
  fHistEventSelectionQA->Write();
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
}

//-----------------------------------------------------------------------------
void StPicoTrackClusterQA::Clear(Option_t *opt) {
}

//________________________________________________________________________
int StPicoTrackClusterQA::Make()
{  // Main loop, called for each event.
  bool fHaveEmcTrigger = kFALSE;
  bool fHaveMBevent = kFALSE;

/*
  // get muDst 
  mu = (StMuDst*) GetInputDS("MuDst");         
  if(!mu) {
    LOG_WARN << " No muDst! Skip! " << endm;
    return kStWarn;
  }

  // get EmcCollection
  mEmcCol = (StEmcCollection*)mu->emcCollection();
  if(!mEmcCol) return kStWarn;
*/

  // Get PicoDstMaker
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

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();

  // get vertex 3-vector and declare variables
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  // Z-vertex cut - per the Aj analysis (-40, 40)
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // ============================ CENTRALITY ============================== //
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  Int_t RunId = mPicoEvent->runId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  int grefMult = mPicoEvent->grefMult();
  grefmultCorr->init(RunId);
  grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate);
  grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
  Int_t cent16 = grefmultCorr->getCentralityBin16();
  if(cent16 == -1) return kStOk; // maybe kStOk; - this is for lowest multiplicity events 80%+ centrality, cut on them
  Int_t centbin = GetCentBin(cent16, 16);
  fCentralityScaled = centbin*5.0;

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  // trigger information:  // cout<<"istrigger = "<<mPicoEvent->isTrigger(450021)<<endl; // NEW
  FillEmcTriggersHist(hEmcTriggers);

  // Run16 triggers:
  int arrBHT1[] = {7, 15, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
  int arrBHT2[] = {4, 16, 17, 530202, 540203};
  int arrBHT3[] = {17, 520203, 530213};
  int arrMB[] = {520021};
  int arrMB5[] = {1, 43, 45, 520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
  int arrMB10[] = {7, 8, 56, 520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  if(fDebugLevel == kDebugEmcTrigger) cout<<"EventTriggers: ";
  for(unsigned int i=0; i<mytriggers.size(); i++) {
    if(fDebugLevel == kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
    if((fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) && (mytriggers[i] == 450014)) { fHaveMBevent = kTRUE; }
    // FIXME Hard-coded for now
    if((fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) && (DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1)))) { fHaveEmcTrigger = kTRUE; }
    if((fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) && (DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5)))) { fHaveMBevent = kTRUE; }
  }

  if(fDebugLevel == kDebugEmcTrigger) cout<<endl;
  // ======================== end of Triggers ============================= //

//  mPicoDst->printTriggers();
//  mPicoDst->printBEmcPidTraits();
//  mPicoDst->printBTOWHits();

  // Run - Trigger Selection to process jets from
  //if((fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) && (!fEmcTriggerArr[fTriggerEventType])) continue;
  //if((fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) && (fHaveEmcTrigger)) continue; //FIXME//

  // TClonesArrays - not needed anymore
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;

  // note that parameter are globally set, may not need to have params
  RunQA(tracks, clus);


  // fill Event QA after cuts
  //FillEventTriggerQA(fHistEventSelectionQAafterCuts);

  return kStOK;
}

//________________________________________________________________________
void StPicoTrackClusterQA::RunQA(TObjArray *tracks, TObjArray *clus)
{
  // assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  unsigned int ntracks = mPicoDst->numberOfTracks();
  double phi, eta, p, pt, energy;

  // loop over ALL tracks in PicoDst 
  for(unsigned short iTracks=0;iTracks<ntracks;iTracks++){
    // get tracks
    StPicoTrack* trk = (StPicoTrack*)mPicoDst->track(iTracks);
    if(!trk){ continue; }

    // acceptance and kinematic quality cuts
    if(!AcceptTrack(trk, Bfield, mVertex)) { continue; }

    // get momentum
    if(doUsePrimTracks) {
      StThreeVectorF mPMomentum = trk->pMom();
      pt = mPMomentum.perp();
      phi = mPMomentum.phi();
      eta = mPMomentum.pseudoRapidity();
      p = mPMomentum.mag();
      energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
    } else {
      StThreeVectorF mGMomentum = trk->gMom(mVertex, Bfield);
      pt = mGMomentum.perp();
      phi = mGMomentum.phi();
      eta = mGMomentum.pseudoRapidity();
      p = mGMomentum.mag();
      energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
    }

    // test - this is only here to occassionally test track variables
    Short_t charge = trk->charge();         
    cout<<"iTracks = "<<iTracks<<"  P = "<<pt<<"  charge = "<<charge<<"  eta = "<<eta<<"  phi = "<<phi;
    cout<<"  nHitsFit = "<<trk->nHitsFit()<<"  BEmc Index = "<<trk->bemcPidTraitsIndex()<<endl;

    // fill some QA histograms
    fHistNTrackvsPt->Fill(pt);
    fHistNTrackvsPhi->Fill(phi);
    fHistNTrackvsEta->Fill(eta);
    fHistNTrackvsPhivsEta->Fill(phi, eta);

  } // track loop

  // looping over clusters - STAR: matching already done
  // get # of clusters and set variables
  unsigned int nclus = mPicoDst->numberOfBEmcPidTraits();
  int towID, towID2, towID3, clusID;
  StThreeVectorF  towPosition, clusPosition;
  StEmcPosition *mPosition = new StEmcPosition();
  StEmcPosition *mPosition2 = new StEmcPosition();

  // print EMCal cluster info
  if(fDebugLevel == 7) mPicoDst->printBEmcPidTraits();
  if(fDebugLevel == 2) cout<<"nClus = "<<nclus<<endl;  

  // loop over ALL clusters in PicoDst and add to jet //TODO
  for(unsigned short iClus=0;iClus<nclus;iClus++){
    StPicoBEmcPidTraits* cluster = mPicoDst->bemcPidTraits(iClus);
    if(!cluster){ continue; }

    // print index of associated track in the event (debug = 2)
    if(fDebugLevel == 8) cout<<"iClus = "<<iClus<<"  trackIndex = "<<cluster->trackIndex()<<"  nclus = "<<nclus<<endl;

    // use StEmcDetector to get position information
    //StEmcDetector* detector;
    //detector=mEmcCol->detector(kBarrelEmcTowerId);
    //if(!detector) cout<<"don't have detector object"<<endl;

    // cluster and tower ID
    // ID's are calculated as such:
    // mBtowId       = (ntow[0] <= 0 || ntow[0] > 4800) ? -1 : (Short_t)ntow[0];
    // mBtowId23 = (ntow[1] < 0 || ntow[1] >= 9 || ntow[2] < 0 || ntow[2] >= 9) ? -1 : (Char_t)(ntow[1] * 10 + ntow[2]);
    clusID = cluster->bemcId();  // index in bemc point array
    towID = cluster->btowId();   // projected tower Id: 1 - 4800
    towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    if(towID < 0) continue;

/*
    // get tower location - from ID
    float towerEta, towerPhi;    // need to be floats
    mGeom->getEtaPhi(towID,towerEta,towerPhi);
    if(towerPhi < 0)    towerPhi += 2*pi;
    if(towerPhi > 2*pi) towerPhi -= 2*pi;
    if(fDebugLevel == 8) {
      cout<<"towID1 = "<<towID<<"  towID2 = "<<towID2<<"  towID3 = "<<towID3;
      cout<<"   towerEta = "<<towerEta<<"  towerPhi = "<<towerPhi<<endl;
    }
*/

    // cluster and tower position - from vertex and ID
    towPosition = mPosition->getPosFromVertex(mVertex, towID);
    clusPosition = mPosition2->getPosFromVertex(mVertex, clusID);

    // matched track index
    int trackIndex = cluster->trackIndex();
    StPicoTrack* trk = (StPicoTrack*)mPicoDst->track(trackIndex);
    StMuTrack* trkMu = (StMuTrack*)mPicoDst->track(trackIndex);
    if(!trk) continue;
    //if(doUsePrimTracks) { if(!(trk->isPrimary())) continue; } // check if primary 

    StThreeVectorD position, momentum;
    double kilogauss = 1000;
    double tesla = 0.00001;
    double magneticField = Bfield * kilogauss  / tesla; // in Tesla
    bool ok = kFALSE;
    if(mPosition) { 
      ok = mPosition->projTrack(&position,&momentum,trkMu,(Double_t) Bfield); 
    }

    double z, eta, phi, theta;
    eta = position.pseudoRapidity(); 
    phi = position.phi();
    z   = position.z();
    theta = 2*TMath::ATan(exp(-1.0*eta));

    // TEST comparing track position with cluster and tower
    double towPhi = towPosition.phi();
    double towEta = towPosition.pseudoRapidity();
    double pmatchPhi = trk->pMom().phi();
    double gmatchPhi = trk->gMom().phi();
    double pmatchEta = trk->pMom().pseudoRapidity();
    double gmatchEta = trk->gMom().pseudoRapidity();
    double clusPhi = clusPosition.phi();
    double clusEta = clusPosition.pseudoRapidity();
    if(towPhi < 0) towPhi += 2*pi;
    if(towPhi > 2*pi) towPhi -= 2*pi;
    if(gmatchPhi < 0) gmatchPhi += 2*pi;
    if(gmatchPhi > 2*pi) gmatchPhi -= 2*pi;
    if(clusPhi < 0) clusPhi += 2*pi;
    if(clusPhi > 2*pi) clusPhi -= 2*pi;

    // fill QA histos for towers
    fHistNTowervsE->Fill(cluster->bemcE0());
    fHistNTowervsPhi->Fill(towPhi);
    fHistNTowervsEta->Fill(towEta);
    fHistNTowervsPhivsEta->Fill(towPhi, towEta);

    // print some tower / cluster / track debug info
    //if(fDebugLevel == 8) {
    if(fDebugLevel == 8 && cluster->bemcE0() > 1.0) {
      cout<<"tID = "<<towID<<" cID = "<<clusID<<" iTrk = "<<trackIndex<<" TrkID = "<<trk->id();
      cout<<" tEta = "<<towEta<<" tPhi = "<<towPhi<<"  towE = "<<cluster->bemcE0();
      cout<<" mTowE = "<<cluster->btowE()<<"  mTowE2 = "<<cluster->btowE2()<<"  mTowE3 = "<<cluster->btowE3()<<endl;
      cout<<"Track: -> trkPhi = "<<gmatchPhi<<" trkEta = "<<gmatchEta<<"  trkP = "<<trk->gMom().mag()<<endl;
      cout<<"Project trk -> eta = "<<eta<<"  phi = "<<phi<<"  z = "<<z;
      cout<<"  etaDist = "<<cluster->btowEtaDist()<<"  phiDist = "<<cluster->btowPhiDist()<<endl;
      cout<<"Cluster: cID = "<<clusID<<"  iClus = "<<iClus<<"  cEta = "<<clusEta<<"  cPhi = "<<clusPhi<<"  clusE = "<<cluster->bemcE()<<endl<<endl;
    } // debug statement
  } // cluster loop

}  // track/cluster QA

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

  Bool_t doAnalysis;
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

//________________________________________________________________________
Bool_t StPicoTrackClusterQA::AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert) {
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

  // do pt cut here to accommadate either type of track
  if(doUsePrimTracks) { // primary  track
    if(pt < fTrackPtMinCut) return kFALSE;
  } else { // global track
    if(pt < fTrackPtMinCut) return kFALSE;
  }

  // track pt, eta, phi cuts
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

//_________________________________________________________________________
TH1* StPicoTrackClusterQA::FillEmcTriggersHist(TH1* h) {
  // number of Emcal Triggers
  for(int i=0; i<7; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  //if(fDebugLevel == kDebugEmcTrigger) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  //static StPicoEmcTrigger* emcTrigger(int i) { return (StPicoEmcTrigger*)picoArrays[picoEmcTrigger]->UncheckedAt(i); }
  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    StPicoEmcTrigger *emcTrig = mPicoDst->emcTrigger(i);
    if(!emcTrig) continue;

    // print some EMCal Trigger info
    if(fDebugLevel == kDebugEmcTrigger) {
      cout<<"i = "<<i<<"  id = "<<emcTrig->id()<<"  flag = "<<emcTrig->flag()<<"  adc = "<<emcTrig->adc();
      cout<<"  isHT0: "<<emcTrig->isHT0()<<"  isHT1: "<<emcTrig->isHT1();
      cout<<"  isHT2: "<<emcTrig->isHT2()<<"  isHT3: "<<emcTrig->isHT3();
      cout<<"  isJP0: "<<emcTrig->isJP0()<<"  isJP1: "<<emcTrig->isJP1()<<"  isJP2: "<<emcTrig->isJP2()<<endl;
    }

    // fill for valid triggers
    if(emcTrig->isHT0()) { h->Fill(1); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { h->Fill(2); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { h->Fill(3); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { h->Fill(4); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
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

//_____________________________________________________________________________
// Trigger QA histogram, label bins 
TH1* StPicoTrackClusterQA::FillEventTriggerQA(TH1* h) {
  // check and fill a Event Selection QA histogram for different trigger selections after cuts

  // Run14 AuAu 200 GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    int arrBHT1[] = {450201, 450211, 460201};
    int arrBHT2[] = {450202, 450212};
    int arrBHT3[] = {460203, 6, 10, 14, 31, 450213};
    int arrMB[] = {450014};
    int arrMB30[] = {20, 450010, 450020};
    int arrCentral5[] = {20, 450010, 450020};
    int arrCentral[] = {15, 460101, 460111};
    int arrMB5[] = {1, 4, 16, 32, 450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};

    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
    int bin = 0;
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrBHT3, sizeof(arrBHT3)/sizeof(*arrBHT3))) { bin = 4; h->Fill(bin); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); } // MB 
    //if() { bin = 6; h->Fill(bin); } 
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
    // get vector of triggers
    vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
    int bin = 0;

    // hard-coded trigger Ids for run16
    int arrBHT0[] = {520606, 520616, 520626, 520636, 520646, 520656};
    int arrBHT1[] = {7, 15, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
    int arrBHT2[] = {4, 16, 17, 530202, 540203};
    int arrBHT3[] = {17, 520203, 530213};
    int arrMB[] = {520021};
    int arrMB5[] = {1, 43, 45, 520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
    int arrMB10[] = {7, 8, 56, 520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};
    int arrCentral[] = {6, 520101, 520111, 520121, 520131, 520141, 520103, 520113, 520123};

    // fill for kAny
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrBHT1, sizeof(arrBHT1)/sizeof(*arrBHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrBHT2, sizeof(arrBHT2)/sizeof(*arrBHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrBHT3, sizeof(arrBHT3)/sizeof(*arrBHT3))) { bin = 4; h->Fill(bin); } // HT3
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); }  // MB
    //if(mytriggers[i] == 999999) { bin = 6; h->Fill(bin); }
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 7; h->Fill(bin); }// Central-5 & Central-novtx
    //if(mytriggers[i] == 999999) { bin = 8; h->Fill(bin); } 
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

