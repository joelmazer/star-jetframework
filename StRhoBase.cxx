// $Id$
//
// Base class for rho calculation.
// Calculates parameterized rho for given centrality independent of input.
//
// adapated form the AliROOT class AliAnalysisTaskRhoBase.cxx for STAR
// original Author: S.Aiola

#include "StRhoBase.h"

// ROOT includes
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TClonesArray.h>
#include <TGrid.h>
#include <THnSparse.h>
#include "TVector3.h"

// jet-framework includes
#include "StRhoParameter.h"
#include "StJet.h"
#include "StJetMakerTask.h"
#include "StCentMaker.h"

// STAR includes
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// STAR centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StRhoBase)

//________________________________________________________________________
StRhoBase::StRhoBase() : 
  StJetFrameworkPicoBase(),
  fCentralityScaled(0.0),
  ref16(-99),
  ref9(-99),
  fMaxEventTrackPt(30.0),
  fMaxEventTowerEt(1000.0), // 30.0
  fOutRhoName(),
  fOutRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fIsAuAu(kTRUE),
  fOutRho(0),
  fOutRhoScaled(0),
  fCompareRho(0),
  fCompareRhoScaled(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistJetRhovsCent(0),
  fHistNjetvsCent(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fHistRhovsCent(0),
  fHistRhoScaledvsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistRhovsNtrackvsMult(0),
  fHistRhoScaledvsNtrackvsMult(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0),
  mCentMaker(0x0),
  mBaseMaker(0x0)
{
  ;
  mOutName = "";
  fJetMakerName = "";
  fRhoMakerName = "";

  // Constructor.
  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  for (Int_t i = 0; i < 12; i++) {
    fHistNjUEoverNjVsNj[i] = 0;
  }
}

//________________________________________________________________________
StRhoBase::StRhoBase(const char *name, Bool_t histo, const char *outName, const char *jetMakerName) :
  StJetFrameworkPicoBase(name),
  fCentralityScaled(0.0),
  ref16(-99),
  ref9(-99),
  fMaxEventTrackPt(30.0),
  fMaxEventTowerEt(1000.0), // 30.0
  fOutRhoName(),
  fOutRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fIsAuAu(kTRUE),
  fOutRho(0),
  fOutRhoScaled(0),
  fCompareRho(0),
  fCompareRhoScaled(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistJetRhovsCent(0),
  fHistNjetvsCent(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fHistRhovsCent(0),
  fHistRhoScaledvsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistRhovsNtrackvsMult(0),
  fHistRhoScaledvsNtrackvsMult(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0),
  mCentMaker(0x0),
  mBaseMaker(0x0)
{
  ;
  // Constructor.
  mOutName = outName;
  fJetMakerName =jetMakerName;
  fRhoMakerName = name;

  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  for (Int_t i = 0; i < 12; i++) {
    fHistNjUEoverNjVsNj[i] = 0;
  }
  //SetMakeGeneralHistograms(histo);
}

//________________________________________________________________________
StRhoBase::~StRhoBase()
{ /*  */
  // destructor
  delete fHistJetAreavsCent;
  delete fHistJetRhovsCent;
  delete fHistJetPtvsCent;
  delete fHistNjetvsCent;
  delete fHistJetPtvsNtrack;
  delete fHistJetAreavsNtrack;
  delete fHistNjetvsNtrack;
  delete fHistRhovsCent;
  delete fHistRhoScaledvsCent;
  delete fHistDeltaRhovsCent;
  delete fHistDeltaRhoScalevsCent;
  delete fHistRhovsNtrackvsMult;
  delete fHistRhoScaledvsNtrackvsMult;
  delete fHistDeltaRhovsNtrack;
  delete fHistDeltaRhoScalevsNtrack;
  delete fHistRhovsNcluster;
  delete fHistRhoScaledvsNcluster;
  delete fJets;
  delete fBGJets;

  for (Int_t i = 0; i < 4; i++) {
    delete fHistJetNconstVsPt[i];
    delete fHistJetRhovsEta[i];
  }
  for (Int_t i = 0; i < 12; i++) {
    delete fHistNjUEoverNjVsNj[i];
  }
}

//________________________________________________________________________
Int_t StRhoBase::Init()
{
  StJetFrameworkPicoBase::Init();

  // declare histograms
  DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);
 
  fBGJets = new TClonesArray("StJet");
  //fBGJets->SetName("fBGJetsName); // need too add name

  // Init the analysis.
  if(!fOutRho) { fOutRho = new StRhoParameter(fOutRhoName, 0);  }
  if(fScaleFunction && !fOutRhoScaled) { fOutRhoScaled = new StRhoParameter(fOutRhoScaledName, 0); }

/*
  if(!fCompareRhoName.IsNull() && !fCompareRho) {
    fCompareRho = dynamic_cast<StRhoParameter*>(InputEvent()->FindListObject(fCompareRhoName));
  }

  if(!fCompareRhoScaledName.IsNull() && !fCompareRhoScaled) {
    fCompareRhoScaled = dynamic_cast<StRhoParameter*>(InputEvent()->FindListObject(fCompareRhoScaledName));
  }
*/

  return kStOk;
}

//________________________________________________________________________
Int_t StRhoBase::Finish() {
  cout << "StRhoBase::Finish()\n";

/*
  //  Write histos to file and close it.
  if(mOutName!="") {
    // use "RECREATE" to create new file, if it exists it will be overwritten
    // use "UPDATE" to open existing file for writing, if no file exists, it is created
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->cd(fRhoMakerName);
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }
*/

  return kStOK;
}
//
// Function: write histograms and objects to file
//_________________________________________________________________________
void StRhoBase::WriteHistograms() 
{
  fHistRhovsCent->Write();
  fHistRhovsNtrackvsMult->Write();
  fHistRhovsNcluster->Write();
  fHistJetPtvsCent->Write();
  fHistJetAreavsCent->Write();
  fHistJetRhovsCent->Write();
  fHistNjetvsCent->Write();

  fHistJetPtvsNtrack->Write();
  fHistJetAreavsNtrack->Write();
  fHistNjetvsNtrack->Write();

  for(Int_t i=0; i<4; i++) {
    fHistJetNconstVsPt[i]->Write();
    fHistJetRhovsEta[i]->Write();

    for(Int_t j=0; j<3; j++) {
      fHistNjUEoverNjVsNj[i*3+j]->Write();
    }
  }

  fHistDeltaRhovsCent->Write();
  fHistDeltaRhovsNtrack->Write();
  fHistRhoScaledvsCent->Write();
  fHistRhoScaledvsNtrackvsMult->Write();
  fHistRhoScaledvsNcluster->Write();
  fHistDeltaRhoScalevsCent->Write();
  fHistDeltaRhoScalevsNtrack->Write();
} 
//
// Function: user create output objects, called at the beginning of the analysis
//________________________________________________________________________
void StRhoBase::DeclareHistograms()
{
  //if (!fCreateHisto) return; //FIXME

  // ranges for AuAu - TODO update range
  Float_t Ntrackrange[2] = {0, 6000};
  Float_t Mult[2] = {0.,25000.};
  if(!fIsAuAu){
    // set multiplicity related axes to a smaller max value
    Ntrackrange[1] = 200.;
    Mult[1] = 2000.;
  }
 
  int fNbins = 1;
  double fMinBinPt = 0.0;
  double fMaxBinPt = 10.0;
 
  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");

  fHistRhovsNtrackvsMult = new TH3F("fHistRhovsNtrackvsMult", "fHistRhovsNtrackvsMult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, Mult[0], Mult[1]);
  fHistRhovsNtrackvsMult->GetXaxis()->SetTitle("No. of tracks");
  fHistRhovsNtrackvsMult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
  fHistRhovsNtrackvsMult->GetZaxis()->SetTitle("mult");

  fHistRhovsNcluster = new TH2F("fHistRhovsNcluster", "fHistRhovsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
  fHistRhovsNcluster->GetXaxis()->SetTitle("No. of clusters");
  fHistRhovsNcluster->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");

  //if (fJetCollArray.GetEntriesFast()>0) {
    fHistJetPtvsCent = new TH2F("fHistJetPtvsCent", "fHistJetPtvsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetPtvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetPtvsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");

    fHistJetAreavsCent = new TH2F("fHistJetAreavsCent", "fHistJetAreavsCent", 101, -1, 100, 100, 0, 1);
    fHistJetAreavsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetAreavsCent->GetYaxis()->SetTitle("Jet area");

    fHistJetRhovsCent = new TH2F("fHistJetRhovsCent", "fHistJetRhovsCent", 101, -1, 100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistJetRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetRhovsCent->GetYaxis()->SetTitle("Jet #rho (GeV/c)");

    fHistNjetvsCent = new TH2F("fHistNjetvsCent",  "fHistNjetvsCent", 101, -1, 100, 150, -0.5, 149.5);
    fHistNjetvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNjetvsCent->GetYaxis()->SetTitle("No. of jets");

    fHistJetPtvsNtrack = new TH2F("fHistJetPtvsNtrack", "fHistJetPtvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
    fHistJetPtvsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistJetPtvsNtrack->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");

    fHistJetAreavsNtrack = new TH2F("fHistJetAreavsNtrack", "fHistJetAreavsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 100, 0, 1);
    fHistJetAreavsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistJetAreavsNtrack->GetYaxis()->SetTitle("Jet area");

    fHistNjetvsNtrack = new TH2F("fHistNjetvsNtrack", "fHistNjetvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 150, -0.5, 149.5);
    fHistNjetvsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistNjetvsNtrack->GetYaxis()->SetTitle("No. of jets");

    TString name;
    for (Int_t i = 0; i < 4; i++) {
      name = Form("fHistJetNconstVsPt_%d",i);
      fHistJetNconstVsPt[i] = new TH2F(name, name, 150, -0.5, 149.5, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetNconstVsPt[i]->GetXaxis()->SetTitle("No. of constituents");
      fHistJetNconstVsPt[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");

      name = Form("fHistJetRhovsEta_%d",i);
      fHistJetRhovsEta[i] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt*2, 16, -0.8, 0.8);
      fHistJetRhovsEta[i]->GetXaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fHistJetRhovsEta[i]->GetYaxis()->SetTitle("#eta");

      for (Int_t j = 0; j < 3; j++) {
	name = Form("NjUEoverNjVsNj_%d_%d",i,j+1);
	fHistNjUEoverNjVsNj[i*3+j] = new TH2F(name, name, 150, -0.5, 149.5, 120, 0.01, 1.21);
	fHistNjUEoverNjVsNj[i*3+j]->GetXaxis()->SetTitle("N_{jet}");
	fHistNjUEoverNjVsNj[i*3+j]->GetYaxis()->SetTitle("N_{jet_{UE}} / N_{jet}");
      }
    }
  //}
  
  //if (!fCompareRhoName.IsNull()) {
    fHistDeltaRhovsCent = new TH2F("fHistDeltaRhovsCent", "fHistDeltaRhovsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistDeltaRhovsCent->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");

    fHistDeltaRhovsNtrack = new TH2F("fHistDeltaRhovsNtrack", "fHistDeltaRhovsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaRhovsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistDeltaRhovsNtrack->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
  //}

  //if(fScaleFunction) {
    fHistRhoScaledvsCent = new TH2F("fHistRhoScaledvsCent", "fHistRhoScaledvsCent", 101, -1, 100, fNbins, fMinBinPt , fMaxBinPt*2);
    fHistRhoScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledvsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");

    fHistRhoScaledvsNtrackvsMult = new TH3F("fHistRhoScaledvsNtrackvsMult", "fHistRhoScaledvsNtrackvsMult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, Mult[0], Mult[1]);
    fHistRhoScaledvsNtrackvsMult->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoScaledvsNtrackvsMult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fHistRhoScaledvsNtrackvsMult->GetZaxis()->SetTitle("mult");

    fHistRhoScaledvsNcluster = new TH2F("fHistRhoScaledvsNcluster", "fHistRhoScaledvsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistRhoScaledvsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistRhoScaledvsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");

    //if (!fCompareRhoScaledName.IsNull()) {
      fHistDeltaRhoScalevsCent = new TH2F("fHistDeltaRhoScalevsCent", "fHistDeltaRhoScalevsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
      fHistDeltaRhoScalevsCent->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
      
      fHistDeltaRhoScalevsNtrack = new TH2F("fHistDeltaRhoScalevsNtrack", "fHistDeltaRhoScalevsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistDeltaRhoScalevsNtrack->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
    //}
  //}

}

//_______________________________________________________________________
void StRhoBase::Clear(Option_t *opt) 
{

}
//
// Runs the analysis for each event
//________________________________________________________________________
Int_t StRhoBase::Make() 
{
  // zero out these global variables
  fCentralityScaled = 0.0, ref9 = 0, ref16 = 0;

  // get the PicoDstMaker
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

  // get run number, check bad runs list if desired (kFALSE if bad)
  int fRunNumber = mPicoEvent->runId();
  if(doRejectBadRuns) {
    if( !mBaseMaker->IsRunOK(fRunNumber) ) return kStOK;
  }

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

  // get vertex 3 vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // zvertex cut - per the Aj analysis
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOK;

  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  if(!JetMaker) {
    LOG_WARN << " No JetMakerTask! Skip! " << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(JetMaker) {
    fJets =  static_cast<TClonesArray*>(JetMaker->GetJets());
    //fJets->SetName("BGJetsRho"); // name will be that set by specific Maker task
  }
  if(!fJets) return kStWarn;

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
  int fCentBin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  fCentralityScaled = mCentMaker->GetCentScaled();
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 

  // test for now FIXME
  double fCent = 0.0;

  // set Rho value
  Double_t rho = GetRhoFactor(fCent);
  fOutRho->SetVal(rho);

  // if we have a scale function - TODO scale factor for neutral if needed?
  if(fScaleFunction) {
    Double_t rhoScaled = rho * GetScaleFactor(fCent);
    fOutRhoScaled->SetVal(rhoScaled);
  }

  // test
  FillHistograms();

  return kStOk;
}

//________________________________________________________________________
Bool_t StRhoBase::FillHistograms() 
{
  // Fill histograms.

  // test for now FIXME
  double fCent = 0.0;
  int fCentBin = 0; 

  Int_t Ntracks   = 0;
  Int_t Nclusters = 0;
  //if (GetParticleContainer(0)) Ntracks = GetParticleContainer(0)->GetNAcceptedParticles(); // FIXME
  //if (GetClusterContainer(0)) Nclusters = GetClusterContainer(0)->GetNAcceptedClusters(); // FIXME

  if(fJets) {
    Int_t    Njets         = fJets->GetEntries();
    Int_t    NjetAcc       = 0;
    Int_t    NjetUE1Sigma  = 0;
    Int_t    NjetUE2Sigma  = 0;
    Int_t    NjetUE3Sigma  = 0;
    Double_t rhoPlus1Sigma = fOutRho->GetVal() + fInEventSigmaRho;
    Double_t rhoPlus2Sigma = fOutRho->GetVal() + 2*fInEventSigmaRho;
    Double_t rhoPlus3Sigma = fOutRho->GetVal() + 3*fInEventSigmaRho;

    // loop over jets
    for (Int_t i = 0; i < Njets; ++i) { 
      // get pointer to jet
      StJet *jet = static_cast<StJet*>(fJets->At(i));
      if(!jet) { continue; } 
      //if (!AcceptJet(jet)) continue; //FIXME
      
      fHistJetPtvsCent->Fill(fCent, jet->Pt());
      fHistJetAreavsCent->Fill(fCent, jet->Area());
      fHistJetRhovsCent->Fill(fCent, jet->Pt() / jet->Area());
      fHistJetRhovsEta[fCentBin]->Fill(jet->Pt() / jet->Area(), jet->Eta());

      fHistJetPtvsNtrack->Fill(Ntracks, jet->Pt());
      fHistJetAreavsNtrack->Fill(Ntracks, jet->Area());

      fHistJetNconstVsPt[fCentBin]->Fill(jet->GetNumberOfConstituents(), jet->Pt());

      if (jet->Pt() < rhoPlus1Sigma * jet->Area()) NjetUE1Sigma++;
      if (jet->Pt() < rhoPlus2Sigma * jet->Area()) NjetUE2Sigma++;
      if (jet->Pt() < rhoPlus3Sigma * jet->Area()) NjetUE3Sigma++;
      NjetAcc++;
    }
    
    if(NjetAcc > 0) {
      fHistNjUEoverNjVsNj[fCentBin*3  ]->Fill(NjetAcc,1.*NjetUE1Sigma/NjetAcc);
      fHistNjUEoverNjVsNj[fCentBin*3+1]->Fill(NjetAcc,1.*NjetUE2Sigma/NjetAcc);
      fHistNjUEoverNjVsNj[fCentBin*3+2]->Fill(NjetAcc,1.*NjetUE3Sigma/NjetAcc);
    }

    fHistNjetvsCent->Fill(fCent, NjetAcc);
    fHistNjetvsNtrack->Fill(Ntracks, NjetAcc);
  }
  
  fHistRhovsCent->Fill(fCent, fOutRho->GetVal());

  //fHistRhovsNtrackvsMult->Fill(Ntracks, fOutRho->GetVal(), multA+multC);
  fHistRhovsNcluster->Fill(Nclusters, fOutRho->GetVal());
  if(fCompareRho) {
    fHistDeltaRhovsCent->Fill(fCent, fOutRho->GetVal() - fCompareRho->GetVal());
    fHistDeltaRhovsNtrack->Fill(Ntracks, fOutRho->GetVal() - fCompareRho->GetVal());
  }

  // scaled Rho
  if(fOutRhoScaled) {
    fHistRhoScaledvsCent->Fill(fCent, fOutRhoScaled->GetVal());
    //fHistRhoScaledvsNtrackvsMult->Fill(Ntracks, fOutRhoScaled->GetVal(), multA+multC);
    //fHistRhoScaledvsNcluster->Fill(Nclusters,  fOutRhoScaled->GetVal());
    if (fCompareRhoScaled) {
      fHistDeltaRhoScalevsCent->Fill(fCent, fOutRhoScaled->GetVal() - fCompareRhoScaled->GetVal());
      fHistDeltaRhoScalevsNtrack->Fill(Ntracks, fOutRhoScaled->GetVal() - fCompareRhoScaled->GetVal());
    }
  }

  return kTRUE;
}      

//________________________________________________________________________
Double_t StRhoBase::GetRhoFactor(Double_t cent)
{
  // Return rho per centrality.
  Double_t rho = 0;
  if(fRhoFunction) rho = fRhoFunction->Eval(cent);

  return rho;
}

//________________________________________________________________________
Double_t StRhoBase::GetScaleFactor(Double_t cent)
{
  // Get scale factor.
  Double_t scale = 1;
  if(fScaleFunction) scale = fScaleFunction->Eval(cent);

  return scale;
}
//
// Load the scale function from a file.
//________________________________________________________________________
TF1* StRhoBase::LoadRhoFunction(const char* path, const char* name)
{
  // "STARfileLocation" needs to be updated if loading rho function from file - dummy now
  TString fname(path);
  if(fname.BeginsWith("STARfileLocation://")) {
    TGrid::Connect("STARfileLocation://");
  }

  // open file
  TFile *file = TFile::Open(path);

  // does file exist
  if(!file || file->IsZombie()) {
    ::Error("StRhoBase", "Could not open scale function file");
    return 0;
  }

  // get scale function from file
  TF1 *sfunc = dynamic_cast<TF1*>(file->Get(name));
  if(sfunc) {
    ::Info("StRhoBase::LoadRhoFunction", "Scale function %s loaded from file %s.", name, path);
  }
  else {
    ::Error("StRhoBase::LoadRhoFunction", "Scale function %s not found in file %s.", name, path);
    return 0;
  }

  // clone scale function
  fScaleFunction = static_cast<TF1*>(sfunc->Clone());

  // close and delete input file
  file->Close();
  delete file;

  return fScaleFunction;
}
