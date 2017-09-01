// $Id$
//
// Base class for rho calculation.
// Calculates parameterized rho for given centrality independent of input.
//
// adapated form the AliROOT class AliAnalysisTaskRhoBase.cxx for STAR
// original
// Author: S.Aiola

// ROOT includes
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TClonesArray.h>
#include <TGrid.h>
#include "TFile.h"
#include <THnSparse.h>

// STAR includes
#include "StRhoParameter.h"
#include "StJet.h"
#include "StRhoBase.h"
#include "StJetMakerTask.h"
#include "StMyAnalysisMaker.h"
#include "StMaker.h"

#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
//#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// STAR centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StRhoBase)

//________________________________________________________________________
StRhoBase::StRhoBase() : 
  //StJet("StRhoBase", kFALSE),
  fOutRhoName(),
  fOutRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fIsPbPb(kTRUE),
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
  fHistRhovsNtrackvsV0Mult(0),
  fHistRhoScaledvsNtrackvsV0Mult(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0),
  fJets(0),
  mPicoDstMaker(0),
  mPicoDst(0),
  mPicoEvent(0),
  grefmultCorr(0),
  JetMaker(0), 
  mOutName(""),
  fJetMakerName(""),
  fRhoMakerName("")
{
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
  //StJet(name, histo),
  fOutRhoName(),
  fOutRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fIsPbPb(kTRUE),
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
  fHistRhovsNtrackvsV0Mult(0),
  fHistRhoScaledvsNtrackvsV0Mult(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0),
  fJets(0),
  mPicoDstMaker(0),
  mPicoDst(0),
  mPicoEvent(0),
  grefmultCorr(0),
  JetMaker(0),
  mOutName(outName),
  fJetMakerName(jetMakerName),
  fRhoMakerName(name)
{
  // Constructor.

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
  delete fHistRhovsNtrackvsV0Mult;
  delete fHistRhoScaledvsNtrackvsV0Mult;
  delete fHistDeltaRhovsNtrack;
  delete fHistDeltaRhoScalevsNtrack;
  delete fHistRhovsNcluster;
  delete fHistRhoScaledvsNcluster;
  delete fJets;
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
  DeclareHistograms();

  // Create user objects.
  fJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);
 
  // Init the analysis.
  if(!fOutRho) { fOutRho = new StRhoParameter(fOutRhoName, 0);  }
  if(fScaleFunction && !fOutRhoScaled) { fOutRhoScaled = new StRhoParameter(fOutRhoScaledName, 0); }

/*
  if (!fCompareRhoName.IsNull() && !fCompareRho) {
    fCompareRho = dynamic_cast<StRhoParameter*>(InputEvent()->FindListObject(fCompareRhoName));
    if (!fCompareRho) {
      Form("%s: Could not retrieve rho %s!", GetName(), fCompareRhoName.Data());
    }
  }

  if (!fCompareRhoScaledName.IsNull() && !fCompareRhoScaled) {
    fCompareRhoScaled = dynamic_cast<StRhoParameter*>(InputEvent()->FindListObject(fCompareRhoScaledName));
    if (!fCompareRhoScaled) {
      Form("%s: Could not retrieve rho %s!", GetName(), fCompareRhoScaledName.Data());
    }
  }
*/

  // initialize centrality correction
  grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();

  return kStOk;
}

//________________________________________________________________________
Int_t StRhoBase::Finish() {
  //  Summarize the run.
  cout << "StRhoBase::Finish()\n";

/*
  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Close();
  }
*/

/*
  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->cd(fRhoMakerName);
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }
*/

  //cout<<"End of StRhoBase::Finish"<<endl;

  return kStOK;
}

//_________________________________________________________________________
void StRhoBase::WriteHistograms() 
{
  // added Jul17, 2017
  fHistRhovsCent->Write();
  fHistRhovsNtrackvsV0Mult->Write();
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
  fHistRhoScaledvsNtrackvsV0Mult->Write();
  fHistRhoScaledvsNcluster->Write();
  fHistDeltaRhoScalevsCent->Write();
  fHistDeltaRhoScalevsNtrack->Write();
} 

//________________________________________________________________________
void StRhoBase::DeclareHistograms()
{
  // User create output objects, called at the beginning of the analysis.
  //if (!fCreateHisto) return; //FIXME

  //ranges for PbPb, change for AuAug -FIXME
  Float_t Ntrackrange[2] = {0, 6000};
  Float_t V0Mult[2] = {0.,25000.};
  if(!fIsPbPb){
     //set multiplicity related axes to a smaller max value
     Ntrackrange[1] = 200.;
     V0Mult[1] = 2000.;
  }
 
  int fNbins = 1;
  double fMinBinPt = 0.0;
  double fMaxBinPt = 10.0;
 
  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");

  //if (fParticleCollArray.GetEntriesFast()>0) {
    fHistRhovsNtrackvsV0Mult = new TH3F("fHistRhovsNtrackvsV0Mult", "fHistRhovsNtrackvsV0Mult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, V0Mult[0], V0Mult[1]);
    fHistRhovsNtrackvsV0Mult->GetXaxis()->SetTitle("No. of tracks");
    fHistRhovsNtrackvsV0Mult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fHistRhovsNtrackvsV0Mult->GetZaxis()->SetTitle("V0 mult");
  //}

  //if (fClusterCollArray.GetEntriesFast()>0) {
    fHistRhovsNcluster = new TH2F("fHistRhovsNcluster", "fHistRhovsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistRhovsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistRhovsNcluster->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
  //}

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

    //if (fParticleCollArray.GetEntriesFast()>0) {
      fHistJetPtvsNtrack = new TH2F("fHistJetPtvsNtrack", "fHistJetPtvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
      fHistJetPtvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetPtvsNtrack->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");

      fHistJetAreavsNtrack = new TH2F("fHistJetAreavsNtrack", "fHistJetAreavsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 100, 0, 1);
      fHistJetAreavsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetAreavsNtrack->GetYaxis()->SetTitle("Jet area");

      fHistNjetvsNtrack = new TH2F("fHistNjetvsNtrack", "fHistNjetvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 150, -0.5, 149.5);
      fHistNjetvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistNjetvsNtrack->GetYaxis()->SetTitle("No. of jets");
    //}

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

    //if (fParticleCollArray.GetEntriesFast()>0) {
      fHistDeltaRhovsNtrack = new TH2F("fHistDeltaRhovsNtrack", "fHistDeltaRhovsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhovsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistDeltaRhovsNtrack->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
    //}
  //}

  //if (fScaleFunction) {
    fHistRhoScaledvsCent = new TH2F("fHistRhoScaledvsCent", "fHistRhoScaledvsCent", 101, -1, 100, fNbins, fMinBinPt , fMaxBinPt*2);
    fHistRhoScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledvsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");

    //if (fParticleCollArray.GetEntriesFast()>0) {
      fHistRhoScaledvsNtrackvsV0Mult = new TH3F("fHistRhoScaledvsNtrackvsV0Mult", "fHistRhoScaledvsNtrackvsV0Mult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, V0Mult[0], V0Mult[1]);
      fHistRhoScaledvsNtrackvsV0Mult->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledvsNtrackvsV0Mult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fHistRhoScaledvsNtrackvsV0Mult->GetZaxis()->SetTitle("V0 mult");
    //}

    //if (fClusterCollArray.GetEntriesFast()>0) {
      fHistRhoScaledvsNcluster = new TH2F("fHistRhoScaledvsNcluster", "fHistRhoScaledvsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
      fHistRhoScaledvsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoScaledvsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");
    //}

    //if (!fCompareRhoScaledName.IsNull()) {
      fHistDeltaRhoScalevsCent = new TH2F("fHistDeltaRhoScalevsCent", "fHistDeltaRhoScalevsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
      fHistDeltaRhoScalevsCent->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
      
      //if (fParticleCollArray.GetEntriesFast()>0) {
	fHistDeltaRhoScalevsNtrack = new TH2F("fHistDeltaRhoScalevsNtrack", "fHistDeltaRhoScalevsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
	fHistDeltaRhoScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
	fHistDeltaRhoScalevsNtrack->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
      //}
    //}
  //}

}

//_______________________________________________________________________
void StRhoBase::Clear(Option_t *opt) 
{

}

//________________________________________________________________________
Int_t StRhoBase::Make() 
{
  // Run the analysis.
  // test for now //FIXME
  double fCent = 0.0;  

  // get the PicoDstMaker
  mPicoDstMaker = (StPicoDstMaker*)GetMaker("picoDst");
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
    //return kStFatal;
  }

  // construct PicoDst object from maker
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
    //return kStFatal;
  }

  // create pointer to PicoEvent
  mPicoEvent = mPicoDst->event();
  if(!mPicoEvent) {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
    //return kStFatal;
  }

  // get vertex 3 vector and declare variables
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // zvertex cut - per the Aj analysis
  if((1.0*TMath::Abs(zVtx)) > 40) return kStWarn; //kStFatal;

  // get JetMaker
  JetMaker = (StJetMakerTask*)GetMaker(fJetMakerName);
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn; //kStFatal;
  }

  // if we have JetMaker, get jet collection associated with it
  if(JetMaker) {
    fJets =  JetMaker->GetJets();
    //fJets->SetName("BGJetsRho"); // name will be that set by specific Maker task
  }
  if(!fJets) return kStWarn; //kStFatal;

  // get run # for centrality correction
  Int_t RunId = mPicoEvent->runId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Float_t fZDCCoincidenceRate = mPicoEvent->ZDCx();

  // Centrality correction calculation
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  int grefMult = mPicoEvent->grefMult();
  int refMult = mPicoEvent->refMult();
  grefmultCorr->init(RunId);
  grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate);
  Int_t cent16 = grefmultCorr->getCentralityBin16();
  Int_t cent9 = grefmultCorr->getCentralityBin9();
  Int_t centbin = GetCentBin(cent16, 16);

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

  // test for now // FIXME
  double fCent = 0.0;
  int fCentBin = 0; 

  Int_t Ntracks   = 0;
  Int_t Nclusters = 0;
  //if (GetParticleContainer(0)) Ntracks = GetParticleContainer(0)->GetNAcceptedParticles(); //FIXME
  //if (GetClusterContainer(0)) Nclusters = GetClusterCon<F12>tainer(0)->GetNAcceptedClusters(); //FIXME

  if(fJets) {
    Int_t    Njets         = fJets->GetEntries();
    Int_t    NjetAcc       = 0;
    Int_t    NjetUE1Sigma  = 0;
    Int_t    NjetUE2Sigma  = 0;
    Int_t    NjetUE3Sigma  = 0;
    Double_t rhoPlus1Sigma = fOutRho->GetVal() + fInEventSigmaRho;
    Double_t rhoPlus2Sigma = fOutRho->GetVal() + 2*fInEventSigmaRho;
    Double_t rhoPlus3Sigma = fOutRho->GetVal() + 3*fInEventSigmaRho;

    for (Int_t i = 0; i < Njets; ++i) { 
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
    
    if(NjetAcc>0) {
      fHistNjUEoverNjVsNj[fCentBin*3  ]->Fill(NjetAcc,1.*NjetUE1Sigma/NjetAcc);
      fHistNjUEoverNjVsNj[fCentBin*3+1]->Fill(NjetAcc,1.*NjetUE2Sigma/NjetAcc);
      fHistNjUEoverNjVsNj[fCentBin*3+2]->Fill(NjetAcc,1.*NjetUE3Sigma/NjetAcc);
    }

    fHistNjetvsCent->Fill(fCent, NjetAcc);
    fHistNjetvsNtrack->Fill(Ntracks, NjetAcc);
  }
  
  fHistRhovsCent->Fill(fCent, fOutRho->GetVal());

  //fHistRhovsNtrackvsV0Mult->Fill(Ntracks, fOutRho->GetVal(),multV0A+multV0C);
  fHistRhovsNcluster->Fill(Nclusters, fOutRho->GetVal());
  if(fCompareRho) {
    fHistDeltaRhovsCent->Fill(fCent, fOutRho->GetVal() - fCompareRho->GetVal());
    fHistDeltaRhovsNtrack->Fill(Ntracks, fOutRho->GetVal() - fCompareRho->GetVal());
  }

  // scaled Rho
  if(fOutRhoScaled) {
    fHistRhoScaledvsCent->Fill(fCent, fOutRhoScaled->GetVal());
    //fHistRhoScaledvsNtrackvsV0Mult->Fill(Ntracks, fOutRhoScaled->GetVal(),multV0A+multV0C);
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

//________________________________________________________________________
TF1* StRhoBase::LoadRhoFunction(const char* path, const char* name)
{
  // Load the scale function from a file.
  TString fname(path);
  if(fname.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }

  TFile* file = TFile::Open(path);

  if(!file || file->IsZombie()) {
    ::Error("AddTaskRho", "Could not open scale function file");
    return 0;
  }

  TF1* sfunc = dynamic_cast<TF1*>(file->Get(name));
  if(sfunc) {
    ::Info("StRhoBase::LoadRhoFunction", "Scale function %s loaded from file %s.", name, path);
  }
  else {
    ::Error("StRhoBase::LoadRhoFunction", "Scale function %s not found in file %s.", name, path);
    return 0;
  }

  fScaleFunction = static_cast<TF1*>(sfunc->Clone());

  file->Close();
  delete file;

  return fScaleFunction;
}

//________________________________________________________________________
Int_t StRhoBase::GetCentBin(Int_t cent, Int_t nBin) const
{  // Get centrality bin.
  Int_t centbin = -1;

  if(nBin == 16) {
    if(cent == 0)  centbin = 15;
    if(cent == 1)  centbin = 14;
    if(cent == 2)  centbin = 13;
    if(cent == 3)  centbin = 12;
    if(cent == 4)  centbin = 11;
    if(cent == 5)  centbin = 10;
    if(cent == 6)  centbin = 9;
    if(cent == 7)  centbin = 8;
    if(cent == 8)  centbin = 7;
    if(cent == 9)  centbin = 6;
    if(cent == 10) centbin = 5;
    if(cent == 11) centbin = 4;
    if(cent == 12) centbin = 3;
    if(cent == 13) centbin = 2;
    if(cent == 14) centbin = 1;
    if(cent == 15) centbin = 0;
  }
  if(nBin == 9) {
    if(cent == 0)  centbin = 8;
    if(cent == 1)  centbin = 7;
    if(cent == 2)  centbin = 6;
    if(cent == 3)  centbin = 5;
    if(cent == 4)  centbin = 4;
    if(cent == 5)  centbin = 3;
    if(cent == 6)  centbin = 2;
    if(cent == 7)  centbin = 1;
    if(cent == 8)  centbin = 0;
  }

  return centbin;
}
