// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// track and tower input to cluster over and create jets
//      - leading jet tag
//      - access to jet constituents
//      - general QA
//
// ################################################################
// $Id$

#include "StJetMakerTask.h"

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
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"

#include "StPicoConstants.h"

// test for clusters: TODO
#include "StMuDSTMaker/COMMON/StMuTrack.h"
//StEmc
#include "StEmcClusterCollection.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
////#include "StEmcUtil/others/emcDetectorName.h" ?
#include "StEmcUtil/projection/StEmcPosition.h"
class StEmcPosition;
class StEEmcCluster;

// jet class and fastjet wrapper
#include "StJet.h"
#include "StFJWrapper.h"
#include "StJetFrameworkPicoBase.h"

// centrality
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;

const Int_t StJetMakerTask::fgkConstIndexShift = 100000;

ClassImp(StJetMakerTask)

//________________________________________________________________________
StJetMakerTask::StJetMakerTask() : 
  StMaker(),
  doWriteHistos(kFALSE),
  doUsePrimTracks(kFALSE), 
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  mOutName(""),
  fTracksName(""),
  fCaloName(""),
  fJetsName(""),
  fJetAlgo(1), 
  fJetType(0), 
  fRecombScheme(fastjet::BIpt_scheme),
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
  fMaxJetTrackPt(20.0),
  fMinJetClusPt(0.15),
  fMinJetClusE(0.2),
  fMinJetTowerE(0.2),
  fJetTrackEtaMin(-1.0),
  fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0),
  fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  mTowerEnergyMin(0.2),
  mHadronicCorrFrac(1.),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0x0),
  fConstituents(0),
  fJetsConstit(0x0),
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
//  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
  mu(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  grefmultCorr(0x0)
//  fJetMakerName("")
{
  // Default constructor.
  for(int i=0; i<4801; i++) { 
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE; 
  }

}

//________________________________________________________________________
StJetMakerTask::StJetMakerTask(const char *name, double mintrackPt = 0.20, bool doHistos = kFALSE, const char* outName = "") : 
  StMaker(name),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  fRunFlag(0),       // see StJetFrameworkPicoBase::fRunFlagEnum
  fCentralityDef(4), // see StJetFrameworkPicoBase::fCentralityDefEnum
  fEventZVtxMinCut(-40.0), 
  fEventZVtxMaxCut(40.0),
  Bfield(0.0),
  mVertex(0x0),
  zVtx(0.0),
  mOutName(outName),
  fTracksName("Tracks"),
  fCaloName("Clusters"),
  fJetsName("Jets"),
  fJetAlgo(1), 
  fJetType(0),
  fRecombScheme(fastjet::BIpt_scheme),
  fjw(name, name),
  fRadius(0.4),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10), 
  fJetPhiMax(+10),
  fJetEtaMin(-0.6), 
  fJetEtaMax(0.6),
  fGhostArea(0.005),
  fMinJetTrackPt(mintrackPt), //0.20
  fMaxJetTrackPt(20.0), 
  fMinJetClusPt(0.15),
  fMinJetClusE(0.2),
  fMinJetTowerE(0.2),
  fJetTrackEtaMin(-1.0), 
  fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0),
  fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  mTowerEnergyMin(0.2),
  mHadronicCorrFrac(1.),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0x0),
  fConstituents(0),
  fJetsConstit(0x0),
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
//  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
  mu(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0),
  grefmultCorr(0x0)
//  fJetMakerName("")
{
  // Standard constructor.
  for(int i=0; i<4801; i++) {
    fTowerToTriggerTypeHT1[i] = kFALSE;
    fTowerToTriggerTypeHT2[i] = kFALSE;
    fTowerToTriggerTypeHT3[i] = kFALSE;
  }

  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StJetMakerTask::~StJetMakerTask()
{
  // Destructor
  //fJets->Clear(); delete fJets;
  if(fHistJetNTrackvsPt)       delete fHistJetNTrackvsPt;
  if(fHistJetNTrackvsPhi)      delete fHistJetNTrackvsPhi;
  if(fHistJetNTrackvsEta)      delete fHistJetNTrackvsEta;
  if(fHistJetNTrackvsPhivsEta) delete fHistJetNTrackvsPhivsEta;
  if(fHistJetNTowervsE)        delete fHistJetNTowervsE;
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
}

//-----------------------------------------------------------------------------
Int_t StJetMakerTask::Init() {
  DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  fJets->SetName(fJetsName);

  // may need array (name hard-coded, Feb20, 2018)
  fJetsConstit = new TClonesArray("StPicoTrack");
  fJetsConstit->SetName("JetConstituents");

  // ============================ Do some jet stuff =======================
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
  if (fRecombScheme == 7)     recombScheme = fastjet::WTA_modp_scheme;
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
  fjw.SetAlgorithm(algorithm); //fJetAlgo);
  fjw.SetRecombScheme(recombScheme);  //fRecombScheme);
  fjw.SetMaxRap(1);

  // setting legacy mode
  //if(fLegacyMode) {
  //  fjw.SetLegacyMode(kTRUE);
  //}

  // Setup container utils. Must be called after Init() so that the
  // containers' arrays are setup.
  //fClusterContainerIndexMap.CopyMappingFrom(StClusterContainer::GetEmcalContainerIndexMap(), fClusterCollArray);
  //fParticleContainerIndexMap.CopyMappingFrom(StParticleContainer::GetEmcalContainerIndexMap(), fParticleCollArray);

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

//
//_______________________________________________________________________________________
Int_t StJetMakerTask::Finish() {
  //  Summarize the run.
/*
  cout<<"StJetMakerTask::Finish()\n";
  //cout << "\tProcessed " << mEventCounter << " events." << endl;
  //cout << "\tWithout PV cuts: "<< mAllPVEventCounter << " events"<<endl;
  //cout << "\tInput events: "<<mInputEventCounter<<endl;

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

/*   ===== test ======
  if(mOutName!="") {
    cout<<"checking output file in StJetMakerTask::Finish().."<<endl;
    TFile *fout = new TFile(mOutName.Data(),"UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Close();
  }
*/

  return kStOK;
}

//________________________________________________________________________
void StJetMakerTask::DeclareHistograms() {
    // declare histograms
    double pi = 1.0*TMath::Pi();
    fHistJetNTrackvsPt = new TH1F("fHistJetNTrackvsPt", "Jet track constituents vs p_{T}", 100, 0., 20.);
    fHistJetNTrackvsPhi = new TH1F("fHistJetNTrackvsPhi", "Jet track constituents vs #phi", 72, 0., 2*pi);
    fHistJetNTrackvsEta = new TH1F("fHistJetNTrackvsEta", "Jet track constituents vs #eta", 40, -1.0, 1.0);
    fHistJetNTrackvsPhivsEta = new TH2F("fHistJetNTrackvsPhivsEta", "Jet track constituents vs #phi vs #eta", 144, 0, 2*pi, 20, -1.0, 1.0);
    fHistJetNTowervsE = new TH1F("fHistJetNTowervsE", "Jet tower constituents vs energy", 100, 0., 20.0);
    fHistJetNTowervsPhi = new TH1F("fHistJetNTowervsPhi", "Jet tower constituents vs #phi", 144, 0., 2*pi);
    fHistJetNTowervsEta = new TH1F("fHistJetNTowervsEta", "Jet tower constituents vs #eta", 40, -1.0, 1.0);
    fHistJetNTowervsPhivsEta = new TH2F("fHistJetNTowervsPhivsEta", "Jet tower constituents vs #phi vs #eta", 144, 0, 2*pi, 20, -1.0, 1.0);

    fHistNJetsvsPt = new TH1F("fHistNJetsvsPt", "NJets vs p_{T}", 100, 0., 100.0);
    fHistNJetsvsPhi = new TH1F("fHistNJetsvsPhi", "NJets vs #phi", 144, 0., 2*pi);
    fHistNJetsvsEta = new TH1F("fHistNJetsvsEta", "NJets vs #eta", 40, -1.0, 1.0);
    fHistNJetsvsPhivsEta = new TH2F("fHistNJetsvsPhivsEta", "NJets vs #phi vs #eta", 144, 0, 2*pi, 20, -1.0, 1.0);
    fHistNJetsvsArea = new TH1F("fHistNJetsvsArea", "NJets vs jet area", 100, 0.0, 1.0);
    fHistNJetsvsNConstituents = new TH1F("fHistNJetsvsNConstituents", "NJets vs NConstit", 51, -0.5, 50.5);
    fHistNJetsvsNTracks = new TH1F("fHistNJetsvsNTracks", "NJets vs NTracks", 50, -0.5, 50.5);
    fHistNJetsvsNTowers = new TH1F("fHistNJetsvsNTowers", "NJets vs NTowers", 50, -0.5, 50.5);
}

//________________________________________________________________________

void StJetMakerTask::WriteHistograms() {
  // write histograms
  fHistJetNTrackvsPt->Write();
  fHistJetNTrackvsPhi->Write();
  fHistJetNTrackvsEta->Write();
  fHistJetNTrackvsPhivsEta->Write();
  fHistJetNTowervsE->Write();
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
}

//-----------------------------------------------------------------------------
void StJetMakerTask::Clear(Option_t *opt) {
  // clear or delete objects after running
  fJets->Clear();
}

//________________________________________________________________________
int StJetMakerTask::Make()
{
  // Main loop, called for each event.
  // ZERO's out the jet array
  fJets->Delete();

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

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  // Z-vertex cut
  // per the Aj analysis (-40, 40)
  // TEST: kStOk -> kStErr // Error, drop this and go to the next event
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk; //Pico::kSkipThisEvent; //kStOk; //kStWarn;

  // ============================ CENTRALITY ============================== //
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  int RunId = mPicoEvent->runId();
  float fBBCCoincidenceRate = mPicoEvent->BBCx();
  int grefMult = mPicoEvent->grefMult();
  grefmultCorr->init(RunId);
  grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate);
  grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
  int cent16 = grefmultCorr->getCentralityBin16();
  if(cent16 == -1) return kStOk; // - this is for lowest multiplicity events 80%+ centrality, cut on them
  int centbin = GetCentBin(cent16, 16);

  // cut on centrality for analysis before doing anything
  // TEST: kStOk -> kStErr // Error, drop this and go to the next event
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; } // Pico::kSkipThisEvent; }

  // TClonesArrays - not needed anymore
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;

/*
//  tracks = dynamic_cast<TClonesArray*>
//  tracks = mPicoDst->picoArray(StPicoArrays::picoArrayNames("Tracks"));
//  tracks = StPicoDst::picoArrays[StPicoArrays::Track];
  //tracks = mPicoDst->picoArray(StPicoArrays::Tracks);  // Track -> Tracks Aug17
  if(!tracks) { return kStWarn; }
*/

/*
  // TRACKS: for FULL or Charged jets
  if ((fJetType==0)||(fJetType==1)) { }

  // NEUTRAL: for FULL or Neutral jets
  if ((fJetType==0)||(fJetType==2)) { }
*/      

  // Find jets
  // note that parameter are globally set, may not need to have params
  FindJets(tracks, clus, fJetAlgo, fRadius);

  // I think this is already working
  //FillJetBranch();

  return kStOK;
}

//________________________________________________________________________
void StJetMakerTask::FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius)
{
  // Find jets.
  // clear out existing wrapper object
  fjw.Clear();

  // initialize Emc position objects
  StEmcPosition *mPosition = new StEmcPosition();
  StEmcPosition *mPosition2 = new StEmcPosition();

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3-vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // assume neutral pion mass
  // additional parameters constructed
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV
  unsigned int ntracks = mPicoDst->numberOfTracks();

  // loop over ALL tracks in PicoDst and add to jet, after acceptance and quality cuts 
  if((fJetType == kFullJet) || (fJetType == kChargedJet)) {
    for(unsigned short iTracks=0; iTracks < ntracks; iTracks++){
      StPicoTrack* trk = (StPicoTrack*)mPicoDst->track(iTracks);
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

      // get momentum vector of track - global or primary track
      StThreeVectorF mTrkMom;
      if(doUsePrimTracks) { 
        // get primary track vector
        mTrkMom = trk->pMom(); 
      } else { 
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield); 
      }

      // track variables
      double pt = mTrkMom.perp();
      double phi = mTrkMom.phi();
      double eta = mTrkMom.pseudoRapidity();
      double px = mTrkMom.x();
      double py = mTrkMom.y();
      double pz = mTrkMom.z();
      double p = mTrkMom.mag();
      double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);

      // test - this is only here to occassionally test track variables
      short charge = trk->charge();         
      //cout<<"iTracks = "<<iTracks<<"  P = "<<pt<<"  charge = "<<charge<<"  eta = "<<eta<<"  phi = "<<phi<<"  nHitsFit = "<<trk->nHitsFit()<<"BEmc Index = "<<trk->bemcPidTraitsIndex()<<endl;

      // send track info to FJ wrapper
      //fjw.AddInputVector(px, py, pz, p, iTracks);    // p -> E
      fjw.AddInputVector(px, py, pz, energy, iTracks); // includes E

    } // track loop
  } // if full/charged jets

/*
  // looping over cluster to add to jet - Example
  if (clus) {
    Double_t vertex[3] = {0, 0, 0};
    const Int_t Nclus = clus->GetEntries();
    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(iClus)); //FIXME
      if (!c->IsEMCAL()) continue;
      TLorentzVector nPart;
      c->GetMomentum(nPart, vertex);
      Double_t energy = nPart.P();
      if (energy<fMinJetClusPt) continue;
      fjw.AddInputVector(nPart.Px(), nPart.Py(), nPart.Pz(), energy, -iClus-1);
    }
  }
*/

  // full or neutral jets - get towers and apply hadronic correction
  if((fJetType == kFullJet) || (fJetType == kNeutralJet)) {
    // looping over clusters to add to jet - STAR: matching already done
    // get # of clusters and set variables
    unsigned int nclus = mPicoDst->numberOfBEmcPidTraits();

/*
    // print EMCal cluster info
    if(fDebugLevel == 7) mPicoDst->printBEmcPidTraits();
    if(fDebugLevel == 2) cout<<"nClus = "<<nclus<<endl;  

    // loop over ALL clusters in PicoDst and add to jet //TODO
    for(unsigned short iClus = 0; iClus < nclus; iClus++){
      StPicoBEmcPidTraits* cluster = mPicoDst->bemcPidTraits(iClus); // NEW usage
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
      int clusID = cluster->bemcId();  // index in bemc point array
      int towID = cluster->btowId();   // projected tower Id: 1 - 4800
      int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
      int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
      if(towID < 0) continue;

      // get tower location - from ID
      float towerEta, towerPhi;    // need to be floats
      mGeom->getEtaPhi(towID,towerEta,towerPhi);
      if(towerPhi < 0)    towerPhi += 2*pi;
      if(towerPhi > 2*pi) towerPhi -= 2*pi;
      if(fDebugLevel == 8) {
        cout<<"towID1 = "<<towID<<"  towID2 = "<<towID2<<"  towID3 = "<<towID3<<"   towerEta = "<<towerEta<<"  towerPhi = "<<towerPhi<<endl;
      }

      // cluster and tower position - from vertex and ID
      StThreeVectorF  towPosition = mPosition->getPosFromVertex(mVertex, towID);
      StThreeVectorF  clusPosition = mPosition2->getPosFromVertex(mVertex, clusID);

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

      double eta = position.pseudoRapidity(); 
      double phi = position.phi();
      double z   = position.z();
      double theta = 2*TMath::ATan(exp(-1.0*eta));

      // TEST comparing track position with cluster and tower
      double towPhi = towPosition.phi();
      double towEta = towPosition.pseudoRapidity();

      // get track variables to matched tower
      StThreeVectorF mTrkMom;
      if(doUsePrimTracks) {
        // get primary track vector
        mTrkMom = trk->pMom();
      } else {
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield);
      }
 
      double matchPhi = mTrkMom.phi();
      double matchEta = mTrkMom.pseudoRapidity();
      double matchP = mTrkMom.mag();
      double clusPhi = clusPosition.phi();
      double clusEta = clusPosition.pseudoRapidity();
      if(towPhi < 0) towPhi += 2*pi;
      if(towPhi > 2*pi) towPhi -= 2*pi;
      if(matchPhi < 0)    matchPhi += 2*pi;
      if(matchPhi > 2*pi) matchPhi -= 2*pi;
      if(clusPhi < 0) clusPhi += 2*pi;
      if(clusPhi > 2*pi) clusPhi -= 2*pi;

      // print some tower / cluster / track debug info
      //if(fDebugLevel == 8) {
      if(fDebugLevel == 8 && cluster->bemcE0() > 1.0) {
        cout<<"tID = "<<towID<<" cID = "<<clusID<<" iTrk = "<<trackIndex<<" TrkID = "<<trk->id();
        cout<<" tEta = "<<towEta<<" tPhi = "<<towPhi<<"  towE = "<<cluster->bemcE0();
        cout<<" mTowE = "<<cluster->btowE()<<"  mTowE2 = "<<cluster->btowE2()<<"  mTowE3 = "<<cluster->btowE3()<<endl;
        cout<<"Track: -> trkPhi = "<<matchPhi<<" trkEta = "<<matchEta<<"  trkP = "<<matchP<<endl;
        cout<<"Project trk -> eta = "<<eta<<"  phi = "<<phi<<"  z = "<<z;
        cout<<"  etaDist = "<<cluster->btowEtaDist()<<"  phiDist = "<<cluster->btowPhiDist()<<endl;
        cout<<"Cluster: cID = "<<clusID<<"  iClus = "<<iClus<<"  cEta = "<<clusEta<<"  cPhi = "<<clusPhi<<"  clusE = "<<cluster->bemcE()<<endl<<endl;
      }

    } // 'cluster' loop
*/

// ==================== March 6th, 2018
    // set / initialize some variables
    double pi = 1.0*TMath::Pi();
    double pi0mass = Pico::mMass[0]; // GeV

    // towerStatus array
    float mTowerMatchTrkIndex[4801] = { 0 };
    bool mTowerStatusArr[4801] = { 0 };
    int matchedTowerTrackCounter = 0;

    // print
    int nTracks = mPicoDst->numberOfTracks();
    int nTrigs = mPicoDst->numberOfEmcTriggers();
    int nBTowHits = mPicoDst->numberOfBTOWHits();
    int nBEmcPidTraits = mPicoDst->numberOfBEmcPidTraits();
//    cout<<"nTracks = "<<nTracks<<"  nTrigs = "<<nTrigs<<"  nBTowHits = "<<nBTowHits<<"  nBEmcPidTraits = "<<nBEmcPidTraits<<endl;
  
    // loop over ALL clusters in PicoDst and add to jet //TODO
    for(unsigned short iClus = 0; iClus < nBEmcPidTraits; iClus++){
      StPicoBEmcPidTraits* cluster = mPicoDst->bemcPidTraits(iClus);
      if(!cluster){ cout<<"Cluster pointer does not exist.. iClus = "<<iClus<<endl; continue; }

      // cluster and tower ID
      int clusID = cluster->bemcId();  // index in bemc point array
      int towID = cluster->btowId();   // projected tower Id: 1 - 4800
      int towID2 = cluster->btowId2(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
      int towID3 = cluster->btowId3(); // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
      if(towID < 0) continue;

      ////////
      // get tower location - from ID: via this method, need to correct eta (when using this method)
      float tEta, tPhi;    // need to be floats
      StEmcGeom *mGeom2 = (StEmcGeom::instance("bemc"));
      //cout<<"radius = "<<mGeom2->Radius()<<endl;
      mGeom2->getEtaPhi(towID,tEta,tPhi);
      if(tPhi < 0)    tPhi += 2*pi;
      if(tPhi > 2*pi) tPhi -= 2*pi;

      // cluster and tower position - from vertex and ID - shouldn't need to correct via this method
      StThreeVectorF towPosition = mPosition->getPosFromVertex(mVertex, towID);
      double towPhi = towPosition.phi();
      double towEta = towPosition.pseudoRapidity();

      // matched track index
      int trackIndex = cluster->trackIndex();
      StPicoTrack* trk = (StPicoTrack*)mPicoDst->track(trackIndex);
      //StMuTrack* mutrk = (StMuTrack*)mPicoDst->track(trackIndex);
      if(!trk) { cout<<"No trk pointer...."<<endl; continue; }
      //if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; } // FIXME - do I want to apply quality cuts to matched track?

      // tower status set - towerID is matched to track passing quality cuts
      mTowerMatchTrkIndex[towID] = trackIndex;
      mTowerStatusArr[towID] = kTRUE;
      matchedTowerTrackCounter++;

      // get track variables to matched tower
      StThreeVectorF mTrkMom;
      if(doUsePrimTracks) { 
        // get primary track vector
        mTrkMom = trk->pMom();
      } else { 
        // get global track vector
        mTrkMom = trk->gMom(mVertex, Bfield); 
      }

      // track properties
      double pt = mTrkMom.perp();
      double phi = mTrkMom.phi();
      double eta = mTrkMom.pseudoRapidity();
      double p = mTrkMom.mag();

    }

    // loop over towers
    int nTowers = mPicoDst->numberOfBTOWHits();
    for(int itow = 0; itow < nTowers; itow++) {
      StPicoBTowHit *tower = mPicoDst->btowHit(itow);
      if(!tower) { cout<<"No tower pointer... iTow = "<<itow<<endl; continue; }

      // tower ID
      int towerID = tower->id();
      if(towerID < 0) continue; // double check these aren't still in the event list

      // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
      StThreeVectorF towerPosition = mPosition->getPosFromVertex(mVertex, towerID);
      double towerPhi = towerPosition.phi();
      if(towerPhi < 0)    towerPhi += 2.0*pi;
      if(towerPhi > 2*pi) towerPhi -= 2.0*pi;
      double towerEta = towerPosition.pseudoRapidity();
      double towX = towerPosition.x();
      double towY = towerPosition.y();
      double towZ = towerPosition.z();
      double towMag = towerPosition.mag();
      int towerADC = tower->adc();
      double towerEunCorr = tower->energy();
      double towerE = tower->energy();

      // tower matched to firing trigger - TODO
      //if(fTowerToTriggerTypeHT1[emcTrigID])
      //if(fTowerToTriggerTypeHT2[emcTrigID])
      //if(fTowerToTriggerTypeHT3[emcTrigID])

      // perform tower cuts
      // if tower was not matched to an accepted track, use it for jet by itself if > 0.2 GeV
      if(mTowerStatusArr[towerID]) {
        //if(mTowerMatchTrkIndex[towerID] > 0) 
        StPicoTrack* trk = (StPicoTrack*)mPicoDst->track( mTowerMatchTrkIndex[towerID] );
        if(!trk) { cout<<"No trk pointer...."<<endl; continue; } // March5, 2018 commented back in TODO
        //if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

        // get track variables to matched tower
        StThreeVectorF mTrkMom;
        if(doUsePrimTracks) { 
          // get primary track vector
          mTrkMom = trk->pMom(); 
        } else { 
          // get global track vector
          mTrkMom = trk->gMom(mVertex, Bfield); 
        }

        // track variables
        double pt = mTrkMom.perp();
        double phi = mTrkMom.phi();
        double eta = mTrkMom.pseudoRapidity();
        double p = mTrkMom.mag();
        double pi0mass = Pico::mMass[0]; // GeV
        double E = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);

        // apply hadronic correction to tower
        towerE = towerEunCorr - (mHadronicCorrFrac * E);
      } 
      // else - no match so treat towers on their own

/*
    // Feb26, 2018: don't think I need this if using getPosFromVertex(vert, id)
    // correct eta for Vz position 
    float theta;
    theta = 2 * atan(exp(-towEta)); // getting theta from eta 
    double z = 0;
    if(towEta != 0) z = 231.0 / tan(theta);  // 231 cm = radius of SMD 
    double zNominal = z - mVertex.z();
    double thetaCorr = atan2(231.0, zNominal); // theta with respect to primary vertex
    float etaCorr =-log(tan(thetaCorr / 2.0)); // eta with respect to primary vertex 
*/

      // cut on tower energy
      if(towerE < 0) towerE = 0.0;
      if(towerE < mTowerEnergyMin) continue;

      // get components from Energy (p - momentum)
      double towerPx = towerE * (towX / towMag);
      double towerPy = towerE * (towY / towMag);
      double towerPz = towerE * (towZ / towMag);
      //cout<<"x = "<<towX<<" y = "<<towY<<" z = "<<towZ<<endl;
      //cout<<"px = "<<towerPx<<"  py = "<<towerPy<<"  pz = "<<towerPz<<"  p = "<<endl;
      //double angle = 1.0*TMath::ATan2(towPy, towPx);
      //cout<<"angle = "<<angle<<"  towerPhi = "<<towerPhi<<endl;

      // add towers to fastjet
      // tower index - TODO check this!
      int uidTow = -(itow + 2);  
      //fjw.AddInputVector(towerPx, towerPy, towerPz, towerE, itow); // includes E
      fjw.AddInputVector(towerPx, towerPy, towerPz, towerE, uidTow); // includes E

      // fill QA histos for towers
      fHistJetNTowervsE->Fill(towerE);
      fHistJetNTowervsPhi->Fill(towerPhi);
      fHistJetNTowervsEta->Fill(towerEta);
      fHistJetNTowervsPhivsEta->Fill(towerPhi, towerEta);

    } // tower loop

// ====================
  } // neutral/full jets

  // run jet finder
  fjw.Run();

  // ======= The below can be a filljetbranch function =======
  // loop over FastJet jets  
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  for(UInt_t ij=0, jetCount=0; ij<jets_incl.size(); ++ij) {
    // PERFORM CUTS ON JETS before saving
    // cut on min jet pt
    if(jets_incl[ij].perp() < fMinJetPt) continue;
    // cut on min jet area
    if(fjw.GetJetArea(ij) < fMinJetArea*TMath::Pi()*fRadius*fRadius) continue;
    // cut on eta acceptance
    if((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax)) continue;
    // cut on phi acceptance 
    if((jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax)) continue;

    // need to figure out how to get m or E from STAR tracks
    //StJet *jet = new ((*fJets)[jetCount]) 
    //  StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());
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
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    fConstituents = fjw.GetJetConstituents(ij); 
    jet->SetJetConstituents(fConstituents);
///////////////////////////////////////////////



///////////////////////////////////////////////
    Double_t neutralE = 0, maxTrack = 0, maxTower = 0;
    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfClusters(constituents.size());
    Int_t nt = 0; // track counter
    Int_t nc = 0; // cluster counter
    Int_t ng = 0; // ghost counter  

    // loop over constituents for ij'th jet
    for(UInt_t ic = 0; ic < constituents.size(); ++ic) {
      // get user defined index
      Int_t uid = constituents[ic].user_index();

      // CHARGED COMPONENT (tracks)
      if(uid >= 0) {
	jet->AddTrackAt(uid, nt);
        StPicoTrack* trk = mPicoDst->track(uid);
        if(!trk) continue;

        // acceptance and kinematic quality cuts
        if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

        // primary track switch
        // get momentum vector of track - global or primary track
        StThreeVectorF mTrkMom;
        if(doUsePrimTracks) {
          // get primary track vector
          mTrkMom = trk->pMom();
        } else {
          // get global track vector
          mTrkMom = trk->gMom(mVertex, Bfield);
        }

        // track variables
        double pt = mTrkMom.perp();
        double phi = mTrkMom.phi();
        double eta = mTrkMom.pseudoRapidity();

        // adjust phi value:  0 < phi < 2pi
        if(phi < 0)    phi+= 2*pi;
        if(phi > 2*pi) phi-= 2*pi;

        // find max track pt
        if(pt > maxTrack) maxTrack = pt;

        // fill some QA histograms
        fHistJetNTrackvsPt->Fill(pt);
        fHistJetNTrackvsPhi->Fill(phi);
        fHistJetNTrackvsEta->Fill(eta);
        fHistJetNTrackvsPhivsEta->Fill(phi, eta);

	nt++;

//============================
/*
      // example usage for container index mapping
      else if (uid >= fgkConstIndexShift) { // track constituent
      Int_t iColl = uid / fgkConstIndexShift;
      Int_t tid = uid - iColl * fgkConstIndexShift;
      iColl--;
      StParticleContainer* partCont = GetParticleContainer(iColl); // FIXME
      if (!partCont) { continue; }
      StVParticle *t = partCont->GetParticle(tid); //FIXME
      if (!t) { continue; }
      jet->AddTrackAt(fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, tid), nt);
*/
// ===========================

      } else { // uid < 0
      // NEUTRAL componenet - TODO

/*
	// example usage
        jet->AddClusterAt(-(uid+1),nc);
        AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(-(uid+1)));
        TLorentzVector nP;
        c->GetMomentum(nP, vertex);
        neutralE += nP.P();
        if (nP.P()>maxCluster)
          maxCluster=nP.P();
	nc++;
*/

        // TODO - finish this if we ever care to look at ghosts.. (not a priority)
        // ghosts
        if(uid == -1) {
          //if(fFillGhost) jet->AddGhost(constituents[ic].px(), constituents[ic].py(), constituents[ic].pz(), constituents[ic].e());
          ng++;
        }

        // neutral towers: start at (index = -2, ghosts are -1, and tracks 0+)
        if(uid < -1) {
          // convert uid to tower id (index of tower)
          Int_t towid = -(uid + 2);
          jet->AddClusterAt(towid, nc);
          StPicoBTowHit *tower = mPicoDst->btowHit(towid);
          if(!tower) continue;

          // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
          StThreeVectorF towerPosition = mPosition->getPosFromVertex(mVertex, towid);
          double towerPhi = towerPosition.phi();
          double towerEta = towerPosition.pseudoRapidity();
          //int towerADC = tower->adc();
          double towE = tower->energy();
          neutralE += towE;
        
          // find max tower E
          if(towE > maxTower) maxTower = towE;

          nc++;
        } // towers

      } // uid < 0

    }  // end of constituent loop

    // set some jet properties
    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SetMaxTrackPt(maxTrack);
    jet->SetMaxClusterPt(maxTower);
    jet->SetNEF(neutralE/jet->E());
    jetCount++;

    // fill jets histograms
    fHistNJetsvsPt->Fill(jet->Pt()); 
    fHistNJetsvsPhi->Fill(jet->Phi());
    fHistNJetsvsEta->Fill(jet->Eta());
    fHistNJetsvsPhivsEta->Fill(jet->Phi(), jet->Eta());
    fHistNJetsvsArea->Fill(jet->Area());  // same as fjw.GetJetArea(ij)
    fHistNJetsvsNConstituents->Fill(nt + nc);
    fHistNJetsvsNTracks->Fill(nt);
    fHistNJetsvsNTowers->Fill(nc);

  }  // end of jet loop

}

// this function is not used anymore! - maybe more stuff 
/**
 * This method fills the jet output branch (TClonesArray) with the jet found by the FastJet
 * wrapper. Before filling the jet branch, the utilities are prepared. Then the utilities are
 * called for each jet and finally after jet finding the terminate method of all utilities is called.
 */
void StJetMakerTask::FillJetBranch()
{
  // loop over fastjet jets
  std::vector<fastjet::PseudoJet> jets_incl = fjw.GetInclusiveJets();
  // sort jets according to jet pt
  static Int_t indexes[9999] = {-1};
  GetSortedArray(indexes, jets_incl);

  Form("%d jets found", (Int_t)jets_incl.size()); //FIXME
  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_incl.size(); ++ijet) {
    Int_t ij = indexes[ijet];
    Form("Jet pt = %f, area = %f", jets_incl[ij].perp(), fjw.GetJetArea(ij)); //FIXME

    if (jets_incl[ij].perp() < fMinJetPt) continue;
    if (fjw.GetJetArea(ij) < fMinJetArea*TMath::Pi()*fRadius*fRadius) continue;
    if ((jets_incl[ij].eta() < fJetEtaMin) || (jets_incl[ij].eta() > fJetEtaMax) ||
        (jets_incl[ij].phi() < fJetPhiMin) || (jets_incl[ij].phi() > fJetPhiMax))
      continue;

    StJet *jet = new ((*fJets)[jetCount])
    		          StJet(jets_incl[ij].perp(), jets_incl[ij].eta(), jets_incl[ij].phi(), jets_incl[ij].m());
    //jet->SetLabel(ij);

    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    // Fill constituent info
    std::vector<fastjet::PseudoJet> constituents(fjw.GetJetConstituents(ij));
    ////FillJetConstituents(jet, constituents, constituents); //FIXME

    Form("Added jet n. %d, pt = %f, area = %f, constituents = %d", jetCount, jet->Pt(), jet->Area(), jet->GetNumberOfConstituents()); //FIXME
    jetCount++;
  }

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

// The below is only useful if I eventually figure out the container mapping indexes to work between the initial
// TClonesArray / Container(Class) to the index created by FastJet and then back later on after jet-finding

/**
 * This method is called for each jet. It loops over the jet constituents and
 * adds them to the jet object.
 * @param jet Pointer to the AliEmcalJet object where the jet constituents will be added
 * @param constituents List of the jet constituents returned by the FastJet wrapper
 * @param constituents_unsub List of jet constituents before background subtraction
 * @param flag If kTRUE it means that the argument "constituents" is a list of subtracted constituents
 * @param particles_sub Array containing subtracted constituents
 */
/*
// this might not be even worth implementing
void StJetMakerTask::FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet>& constituents,
    std::vector<fastjet::PseudoJet>& constituents_unsub, Int_t flag, TString particlesSubName)
{
  Int_t nt            = 0;
  Int_t nc            = 0;
  Double_t neutralE   = 0.;
  Double_t maxCh      = 0.;
  Double_t maxNe      = 0.;
  Int_t ncharged      = 0;
  Int_t nneutral      = 0;
  Double_t mcpt       = 0.;
  TClonesArray * particles_sub = 0;

  Int_t uid   = -1;

  jet->SetNumberOfTracks(constituents.size());
  jet->SetNumberOfClusters(constituents.size());

  for (UInt_t ic = 0; ic < constituents.size(); ++ic) {
    if (flag == 0) {
      uid = constituents[ic].user_index();
    } else {
      if (constituents[ic].perp()<1.e-10) {
        uid=-1;
      } else {
        uid = constituents[ic].user_index();
      }
      if (uid==0) {
        //Form("correspondence between un/subtracted constituent not found");
        continue;
      }
    }

    Form("Processing constituent %d", uid); //FIXME
    if (uid == -1) { //ghost particle

      if (fFillGhost) jet->AddGhost(constituents[ic].px(),
          constituents[ic].py(),
          constituents[ic].pz(),
          constituents[ic].e());
    }	
    else if (uid >= fgkConstIndexShift) { // track constituent
      Int_t iColl = uid / fgkConstIndexShift;
      Int_t tid = uid - iColl * fgkConstIndexShift;
      iColl--;
      //Form("Constituent %d is a track from collection %d and with ID %d", uid, iColl, tid);
      AliParticleContainer* partCont = GetParticleContainer(iColl); // FIXME
      if (!partCont) {
        //Form("Could not find particle container %d",iColl);
        continue;
      }
      StVParticle *t = partCont->GetParticle(tid); //FIXME
      if (!t) {
        //Form("Could not find track %d",tid);
        continue;
      }

      Double_t cEta = t->Eta();
      Double_t cPhi = t->Phi();
      Double_t cPt  = t->Pt();
      Double_t cP   = t->P();
      if (t->Charge() == 0) {
        neutralE += cP;
        ++nneutral;
        if (cPt > maxNe) maxNe = cPt;
      } else {
        ++ncharged;
        if (cPt > maxCh) maxCh = cPt;
      }

      // check if MC particle
      if (TMath::Abs(t->GetLabel()) > fMinMCLabel) mcpt += cPt;

      if (flag == 0 || particlesSubName == "") {
        jet->AddTrackAt(fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, tid), nt);
      }
      else {
        // Get the particle container and array corresponding to the subtracted particles
        partCont = GetParticleContainer(particlesSubName);
        particles_sub = partCont->GetArray();
        // Create the new particle in the particles_sub array and add it to the jet
        Int_t part_sub_id = particles_sub->GetEntriesFast();
        AliEmcalParticle* part_sub = new ((*particles_sub)[part_sub_id]) AliEmcalParticle(dynamic_cast<AliVTrack*>(t));   // SA: probably need to be fixed!!
        part_sub->SetPtEtaPhiM(constituents[ic].perp(),constituents[ic].eta(),constituents[ic].phi(),constituents[ic].m());
        jet->AddTrackAt(fParticleContainerIndexMap.GlobalIndexFromLocalIndex(partCont, part_sub_id), nt);
      }

      ++nt;
    } 
    else if (uid <= -fgkConstIndexShift) { // cluster constituent
      Int_t iColl = -uid / fgkConstIndexShift;
      Int_t cid = -uid - iColl * fgkConstIndexShift;
      iColl--;
      //sprintf(3,Form("Constituent %d is a cluster from collection %d and with ID %d", uid, iColl, cid));
      StClusterContainer* clusCont = GetClusterContainer(iColl);
      AliVCluster *c = clusCont->GetCluster(cid);
      if (!c) continue;

      AliTLorentzVector nP;
      clusCont->GetMomentum(nP, cid);

      Double_t cEta = nP.Eta();
      Double_t cPhi = nP.Phi_0_2pi();
      Double_t cPt  = nP.Pt();
      Double_t cP   = nP.P();

      neutralE += cP;
      if (cPt > maxNe) maxNe = cPt;

      // MC particle
      if (TMath::Abs(c->GetLabel()) > fMinMCLabel) mcpt += c->GetMCEnergyFraction() > 1e-6 ? cPt * c->GetMCEnergyFraction() : cPt;

      if (flag == 0 || particlesSubName == "") {
        jet->AddClusterAt(fClusterContainerIndexMap.GlobalIndexFromLocalIndex(clusCont, cid), nc);
      }
      else {
        // Get the cluster container and array corresponding to the subtracted particles
        clusCont = GetClusterContainer(particlesSubName);
        particles_sub = clusCont->GetArray();
        // Create the new particle in the particles_sub array and add it to the jet
        Int_t part_sub_id = particles_sub->GetEntriesFast();
        AliEmcalParticle* part_sub = new ((*particles_sub)[part_sub_id]) AliEmcalParticle(c);
        part_sub->SetPtEtaPhiM(constituents[ic].perp(),constituents[ic].eta(),constituents[ic].phi(),constituents[ic].m());
        jet->AddClusterAt(fClusterContainerIndexMap.GlobalIndexFromLocalIndex(clusCont, part_sub_id), nc);
      }

      ++nc;
      ++nneutral;
    } 
    else {
      //Form("%s: No logical way to end up here (uid = %d).", GetName(), uid);
      continue;
    }
  }

  jet->SetNumberOfTracks(nt);
  jet->SetNumberOfClusters(nc);
  jet->SetNEF(neutralE / jet->E());
  jet->SetMaxChargedPt(maxCh);
  jet->SetMaxNeutralPt(maxNe);
  jet->SetNumberOfCharged(ncharged);
  jet->SetNumberOfNeutrals(nneutral);
  jet->SetMCPt(mcpt);
  jet->SortConstituents();
}
*/

//________________________________________________________________________
Bool_t StJetMakerTask::AcceptJetTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert) {
  // constants: assume neutral pion mass
  double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0*TMath::Pi();

  // get momentum vector of track - global or primary track
  StThreeVectorF mTrkMom;
  if(doUsePrimTracks) { 
    if(!(trk->isPrimary())) return kFALSE; // check if primary

    // get primary track vector
    mTrkMom = trk->pMom(); 
  } else { 
    // get global track vector
    mTrkMom = trk->gMom(Vert, B); 
  }

  // track variables
  double pt = mTrkMom.perp();
  double phi = mTrkMom.phi();
  double eta = mTrkMom.pseudoRapidity();
  //double px = mTrkMom.x();
  //double py = mTrkMom.y();
  //double pz = mTrkMom.z();
  //double p = mTrkMom.mag();
  //double energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
  //short charge = trk->charge();
  double dca = (trk->dcaPoint() - mPicoEvent->primaryVertex()).mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // do pt cut here to accommadate either type
  if(doUsePrimTracks) { // primary  track
    if(pt < fMinJetTrackPt) return kFALSE;
  } else { // global track
    if(pt < fMinJetTrackPt) return kFALSE;
  }

  // jet track acceptance cuts now - after getting 3vector - hardcoded
  if(pt > fMaxJetTrackPt) return kFALSE; // 20.0 STAR, 100.0 ALICE
  if((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax)) return kFALSE;
  if(phi < 0) phi+= 2*pi;
  if(phi > 2*pi) phi-= 2*pi;
  if((phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) return kFALSE;
      
  // additional quality cuts for tracks
  if(dca > fJetTrackDCAcut) return kFALSE;
  if(nHitsFit < fJetTracknHitsFit) return kFALSE;
  if(nHitsRatio < fJetTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}

//________________________________________________________________________
Int_t StJetMakerTask::GetCentBin(Int_t cent, Int_t nBin) const
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
Bool_t StJetMakerTask::SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut) {
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

  Bool_t doAnalysis = kFALSE; // set false by default, to make sure user chooses an available bin

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

