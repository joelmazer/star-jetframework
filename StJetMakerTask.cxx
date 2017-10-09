// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This JetMakerTask gives access to the FastJetWrapper to give 
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

#include "StPicoConstants.h"

// test for clusters: TODO
#include "StMuDSTMaker/COMMON/StMuTrack.h"
//StEmc
#include "StEmcClusterCollection.h"
#include "StEmcCollection.h"
#include "StEmcCluster.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcUtil/projection/StEmcPosition.h"
class StEmcPosition;
class StEEmcCluster;

// jet class and fastjet wrapper
#include "StJet.h"
#include "StFJWrapper.h"

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
  mOutName(""),
  fTracksName(""),
  fCaloName(""),
  fJetsName(""),
  fRadius(0.4),
  fJetAlgo(1), 
  fJetType(0), 
  fRecombScheme(fastjet::BIpt_scheme),
  fjw("StJetMakerTask", "StJetMakerTask"),
  fMinJetArea(0.001),
  fMinJetPt(1.0),
  fJetPhiMin(-10.0),
  fJetPhiMax(+10.0),
  fJetEtaMin(-0.6),
  fJetEtaMax(0.6),
  fGhostArea(0.005), 
  fMinJetTrackPt(0.2),
  fMaxJetTrackPt(20.0),
  fMinJetClusPt(0.15),
  fJetTrackEtaMin(-1.0),
  fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0),
  fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fConstituents(0),
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
//  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
  mu(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0)
//  fJetMakerName("")
{
  // Default constructor.

}

//________________________________________________________________________
StJetMakerTask::StJetMakerTask(const char *name, double mintrackPt = 0.20, bool doHistos = kFALSE, const char* outName = "") : 
  StMaker(name),
  doWriteHistos(doHistos),
  doUsePrimTracks(kFALSE),
  fDebugLevel(0),
  mOutName(outName),
  fTracksName("Tracks"),
  fCaloName("Clusters"),
  fJetsName("Jets"),
  fRadius(0.4),
  fJetAlgo(1), 
  fJetType(0),
  fRecombScheme(fastjet::BIpt_scheme),
  fjw(name, name),
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
  fJetTrackEtaMin(-1.0), 
  fJetTrackEtaMax(1.0),
  fJetTrackPhiMin(0.0),
  fJetTrackPhiMax(2.0*TMath::Pi()),
  fJetTrackDCAcut(3.0),
  fJetTracknHitsFit(15),
  fJetTracknHitsRatio(0.52),
  fTrackEfficiency(1.),
  fLegacyMode(kFALSE),
  fFillGhost(kFALSE),
  fJets(0),
  fConstituents(0),
  mGeom(StEmcGeom::instance("bemc")),
  mEmcCol(0),
//  fClusterContainerIndexMap(),
  fParticleContainerIndexMap(),
  mu(0x0),
  mPicoDstMaker(0x0),
  mPicoDst(0x0),
  mPicoEvent(0x0)
//  fJetMakerName("")
{
  // Standard constructor.
  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StJetMakerTask::~StJetMakerTask()
{
  // Destructor
  //fJets->Clear(); delete fJets;
  if(fHistJetNTrackvsPt)   delete fHistJetNTrackvsPt;
  if(fHistJetNTrackvsPhi)  delete fHistJetNTrackvsPhi;
  if(fHistJetNTrackvsEta)  delete fHistJetNTrackvsEta;
  if(fHistJetNTowervsE)    delete fHistJetNTowervsE;
  if(fHistJetNTowervsPhi)  delete fHistJetNTowervsPhi;
  if(fHistJetNTowervsEta)  delete fHistJetNTowervsEta;
  if(fHistJetNTowervsPhivsEta) delete fHistJetNTowervsPhivsEta;
}

//-----------------------------------------------------------------------------
Int_t StJetMakerTask::Init() {
  DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  fJets->SetName(fJetsName);

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

  return kStOK;
}

//-----------------------------------------------------------------------------
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

  return kStOK;
}

//________________________________________________________________________
void StJetMakerTask::DeclareHistograms() {
    // declare histograms
    double pi = 1.0*TMath::Pi();
    fHistJetNTrackvsPt = new TH1F("fHistJetNTrackvsPt", "Jet track constituents vs p_{T}", 100, 0., 20.);
    fHistJetNTrackvsPhi = new TH1F("fHistJetNTrackvsPhi", "Jet track constituents vs #phi", 72, 0., 2*pi);
    fHistJetNTrackvsEta = new TH1F("fHistJetNTrackvsEta", "Jet track constituents vs #eta", 40, -1.0, 1.0);
    fHistJetNTowervsE = new TH1F("fHistJetNTowervsE", "Jet tower constituents vs energy", 100, 0., 20.0);
    fHistJetNTowervsPhi = new TH1F("fHistJetNTowervsPhi", "Jet tower constituents vs #phi", 144, 0., 2*pi);
    fHistJetNTowervsEta = new TH1F("fHistJetNTowervsEta", "Jet tower constituents vs #eta", 40, -1.0, 1.0);
    fHistJetNTowervsPhivsEta = new TH2F("fHistJetNTowervsPhivsEta", "Jet track constituents vs #phi vs #eta", 144, 0, 2*pi, 20, -1.0, 1.0);
}

//________________________________________________________________________

void StJetMakerTask::WriteHistograms() {
  // write histograms
  fHistJetNTrackvsPt->Write();
  fHistJetNTrackvsPhi->Write();
  fHistJetNTrackvsEta->Write();
  fHistJetNTowervsE->Write();
  fHistJetNTowervsPhi->Write();
  fHistJetNTowervsEta->Write();
  fHistJetNTowervsPhivsEta->Write();
}

//-----------------------------------------------------------------------------
void StJetMakerTask::Clear(Option_t *opt) {
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

  // get vertex 3 vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // zvertex cut - per the Aj analysis - User can make finer cut in Analysis Code
  if((1.0*TMath::Abs(zVtx)) > 40) return kStOk; // kStWarn;

  // TClonesArrays
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
//  tracks = dynamic_cast<TClonesArray*>
//  tracks = mPicoDst->picoArray(StPicoArrays::picoArrayNames("Tracks"));
//  tracks = StPicoDst::picoArrays[StPicoArrays::Track];
  //tracks = mPicoDst->picoArray(picoTrack); // COMMENTED OUT AUG 6th

/* test out
  tracks = mPicoDst->picoArray(StPicoArrays::Tracks);  // Track -> Tracks Aug17
  if(!tracks) {
    return kStWarn;
  }
*/

/*
  // TRACKS: for FULL or Charged jets
  if ((fJetType==0)||(fJetType==1)) { 
  }

  // NEUTRAL: for FULL or Neutral jets
  if ((fJetType==0)||(fJetType==2)) {
  }
*/      

  // Find jets
  // not that parameter are globally set, may not need to have params
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

  // get event B (magnetic) field
  Float_t Bfield = mPicoEvent->bField();

  // get vertex 3-vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // assume neutral pion mass
  double pi = 1.0*TMath::Pi();
  double pi0mass = Pico::mMass[0]; // GeV

  unsigned int ntracks = mPicoDst->numberOfTracks();
  double phi, eta, px, py, pz, p, pt, energy;

  // loop over ALL tracks in PicoDst and add to jet 
  if((fJetType == kFullJet) || (fJetType == kChargedJet)) {
    for(unsigned short iTracks=0;iTracks<ntracks;iTracks++){
      // get tracks
      StPicoTrack* trk = (StPicoTrack*)mPicoDst->track(iTracks);
      if(!trk){ continue; }

      // acceptance and kinematic quality cuts
      if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

      // get momentum
      if(doUsePrimTracks) {
        StThreeVectorF mPMomentum = trk->pMom();
        pt = mPMomentum.perp();
        phi = mPMomentum.phi();
        eta = mPMomentum.pseudoRapidity();
        px = mPMomentum.x();
        py = mPMomentum.y();
        pz = mPMomentum.z();
        p = mPMomentum.mag();
        energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
      } else {
        StThreeVectorF mGMomentum = trk->gMom(mVertex, Bfield);
        pt = mGMomentum.perp();
        phi = mGMomentum.phi();
        eta = mGMomentum.pseudoRapidity();
        px = mGMomentum.x();
        py = mGMomentum.y();
        pz = mGMomentum.z();
        p = mGMomentum.mag();
        energy = 1.0*TMath::Sqrt(p*p + pi0mass*pi0mass);
      }

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

  if((fJetType == kFullJet) || (fJetType == kNeutralJet)) {
    // looping over clusters to add to jet - STAR: matching already done
    // get # of clusters and set variables
    unsigned int nclus = mPicoDst->numberOfBEmcPidTraits();
    int towID, clusID;
    StThreeVectorF  towPosition, clusPosition;
    StEmcPosition *mPosition = new StEmcPosition();
    StEmcPosition *mPosition2 = new StEmcPosition();

    // print EMCal cluster info
    if(fDebugLevel == 7) mPicoDst->printBEmcPidTraits();
    if(fDebugLevel == 2) cout<<"nClus = "<<nclus<<endl; 
    // loop over ALL clusters in PicoDst and add to jet //TODO
    for(unsigned short iClus=0;iClus<nclus;iClus++){
      StPicoBEmcPidTraits* cluster = mPicoDst->bemcPidTraits(iClus); // NEW usage
      if(!cluster){ continue; }

     if(fDebugLevel==2) cout<<"iClus = "<<iClus<<"trackIndex = "<<cluster->trackIndex()<<endl;

      // use StEmcDetector to get position information
      //StEmcDetector* detector;
      //detector=mEmcCol->detector(kBarrelEmcTowerId);
      //if(!detector) cout<<"don't detector object"<<endl;

      // cluster and tower ID
      clusID = cluster->bemcId();  // index in bemc point array
      towID = cluster->btowId();   // projected tower Id: 1 - 4800
      if(towID < 0) continue;
      float towerEta, towerPhi;    // need to be floats

      // get tower location - from ID
      mGeom->getEtaPhi(towID,towerEta,towerPhi);
      if(towerPhi < 0) towerPhi += 2*pi;
      if(towerPhi > 2*pi) towerPhi -= 2*pi;
      if(fDebugLevel > 8) { 
        cout<<"towerID = "<<towID<<"  towerEta = "<<towerEta<<"  towerPhi = "<<towerPhi<<endl;
      }

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

      double z,eta,phi;
      eta = position.pseudoRapidity(); 
      phi = position.phi();
      z   = position.z();
      double theta = 2*TMath::ATan(exp(-1.0*eta));
      //cout<<"theta = "<<theta<<endl;      

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
      fHistJetNTowervsE->Fill(cluster->bemcE0());
      fHistJetNTowervsPhi->Fill(towPhi);
      fHistJetNTowervsEta->Fill(towEta);
      fHistJetNTowervsPhivsEta->Fill(towPhi, towEta);

      // print some tower / cluster / track debug info
      if(fDebugLevel > 8) {
        cout<<"tID = "<<towID<<" cID = "<<clusID<<" iTrk = "<<trackIndex<<" TrkID = "<<trk->id();
        cout<<" tEta = "<<towEta<<" tPhi = "<<towPhi<<"  towE = "<<cluster->bemcE0();
        cout<<" MatTowE = "<<cluster->btowE()<<"  MatTowE2 = "<<cluster->btowE2()<<"  MatTowE3 = "<<cluster->btowE3()<<endl;
        cout<<"Track: -> trkPhi = "<<gmatchPhi<<" trkEta = "<<gmatchEta<<"  trkP = "<<trk->gMom().mag()<<endl;
        cout<<"Project trk -> eta = "<<eta<<"  phi = "<<phi<<"  z = "<<z;
        cout<<"  etaDist = "<<cluster->btowEtaDist()<<"  phiDist = "<<cluster->btowPhiDist()<<endl;
        cout<<"Cluster: cID = "<<clusID<<"  iClus = "<<iClus<<"  cEta = "<<clusEta<<"  cPhi = "<<clusPhi<<"  clusE = "<<cluster->bemcE()<<endl<<endl;
      }

    }
  }

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

    // area vector added July26
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());  // same as fjw.GetJetArea(ij)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());

    // get constituents of jets
    vector<fastjet::PseudoJet> constituents = fjw.GetJetConstituents(ij);
    fConstituents = fjw.GetJetConstituents(ij); 
    jet->SetJetConstituents(fConstituents);

    Double_t neutralE = 0;Double_t maxTrack = 0;Double_t maxCluster=0;
    jet->SetNumberOfTracks(constituents.size());
    jet->SetNumberOfClusters(constituents.size());
    Int_t nt = 0; // track counter
    Int_t nc = 0; // cluster counter
    Double_t pt, phi, eta;
  
    // loop over constituents for ij'th jet
    for(UInt_t ic=0; ic<constituents.size(); ++ic) {
      Int_t uid = constituents[ic].user_index();

      if (uid>=0){
	jet->AddTrackAt(uid, nt);
        StPicoTrack* trk = mPicoDst->track(uid);
        if(!trk) continue;

        // acceptance and kinematic quality cuts
        if(!AcceptJetTrack(trk, Bfield, mVertex)) { continue; }

        // get momentum
        if(doUsePrimTracks) {
          StThreeVectorF mPMomentum = trk->pMom();
          pt = mPMomentum.perp();
          phi = mPMomentum.phi();
          eta = mPMomentum.pseudoRapidity();
        } else {
          StThreeVectorF mGMomentum = trk->gMom(mVertex, Bfield);
          pt = mGMomentum.perp();
          phi = mGMomentum.phi();
          eta = mGMomentum.pseudoRapidity();
        }

        // adjust phi value:  0 < phi < 2pi
        if(phi < 0) phi+= 2*pi;
        if(phi > 2*pi) phi-= 2*pi;

        // find max track pt
        if (pt>maxTrack) maxTrack=pt;

        // fill some QA histograms
        fHistJetNTrackvsPt->Fill(pt);
        fHistJetNTrackvsPhi->Fill(phi);
        fHistJetNTrackvsEta->Fill(eta);

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
      }
    }

    jet->SetNumberOfTracks(nt);
    jet->SetNumberOfClusters(nc);
    jet->SetMaxTrackPt(maxTrack);
    jet->SetMaxClusterPt(maxCluster); // clusters not added yet
    jet->SetNEF(neutralE/jet->E()); // FIXME : currently not set properly - neutralE = 0
    jetCount++;
  }
}

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
