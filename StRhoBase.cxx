// $Id$
//
// Base class for rho calculation.
// Calculates parameterized rho for given centrality independent of input.
//
// adapated form the AliROOT class AliAnalysisTaskRhoBase.cxx for STAR
// 
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
#include "TFile.h"
#include <THnSparse.h>
#include "TVector3.h"

// STAR includes
#include "StRhoParameter.h"
#include "StJet.h"
#include "StJetMakerTask.h"
#include "StMyAnalysisMaker.h"
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
  doUseBBCCoincidenceRate(kFALSE),
  fMaxEventTrackPt(30.0),
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
  fHistRhoScaledvsNcluster(0)
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
  doUseBBCCoincidenceRate(kFALSE),
  fMaxEventTrackPt(30.0),
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
  fHistRhoScaledvsNcluster(0)
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
        break;

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

  return kStOk;
}

//________________________________________________________________________
Int_t StRhoBase::Finish() {
  //  Summarize the run.
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

  //cout<<"End of StRhoBase::Finish"<<endl;

  return kStOK;
}

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

//________________________________________________________________________
void StRhoBase::DeclareHistograms()
{
  // User create output objects, called at the beginning of the analysis.
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

  //if (fParticleCollArray.GetEntriesFast()>0) {
    fHistRhovsNtrackvsMult = new TH3F("fHistRhovsNtrackvsMult", "fHistRhovsNtrackvsMult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, Mult[0], Mult[1]);
    fHistRhovsNtrackvsMult->GetXaxis()->SetTitle("No. of tracks");
    fHistRhovsNtrackvsMult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fHistRhovsNtrackvsMult->GetZaxis()->SetTitle("mult");
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

  //if(fScaleFunction) {
    fHistRhoScaledvsCent = new TH2F("fHistRhoScaledvsCent", "fHistRhoScaledvsCent", 101, -1, 100, fNbins, fMinBinPt , fMaxBinPt*2);
    fHistRhoScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledvsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");

    //if (fParticleCollArray.GetEntriesFast()>0) {
      fHistRhoScaledvsNtrackvsMult = new TH3F("fHistRhoScaledvsNtrackvsMult", "fHistRhoScaledvsNtrackvsMult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, Mult[0], Mult[1]);
      fHistRhoScaledvsNtrackvsMult->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledvsNtrackvsMult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fHistRhoScaledvsNtrackvsMult->GetZaxis()->SetTitle("mult");
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
//
// Runs the analysis for each event
//________________________________________________________________________
Int_t StRhoBase::Make() 
{
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

  // cut event on max track pt > 30.0 GeV
  if(GetMaxTrackPt() > fMaxEventTrackPt) return kStOK;

  // get vertex 3 vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // zvertex cut - per the Aj analysis
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStWarn;

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

  // get run # for centrality correction
  Int_t RunId = mPicoEvent->runId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Float_t fZDCCoincidenceRate = mPicoEvent->ZDCx();

  // Centrality correction calculation
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  int grefMult = mPicoEvent->grefMult();
  //int refMult = mPicoEvent->refMult();
  Int_t fCentBin, cent16; // cent9

  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorr->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }

    cent16 = grefmultCorr->getCentralityBin16();
    //cent9 = grefmultCorr->getCentralityBin9();
    fCentBin = GetCentBin(cent16, 16); // centbin
  } else { // for pp
    fCentBin = 0, cent16 = 0; 
  }

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStWarn; // maybe kStOk; - this is for lowest multiplicity events 80%+ centrality, cut on them

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

  //fHistRhovsNtrackvsMult->Fill(Ntracks, fOutRho->GetVal(),multA+multC);
  fHistRhovsNcluster->Fill(Nclusters, fOutRho->GetVal());
  if(fCompareRho) {
    fHistDeltaRhovsCent->Fill(fCent, fOutRho->GetVal() - fCompareRho->GetVal());
    fHistDeltaRhovsNtrack->Fill(Ntracks, fOutRho->GetVal() - fCompareRho->GetVal());
  }

  // scaled Rho
  if(fOutRhoScaled) {
    fHistRhoScaledvsCent->Fill(fCent, fOutRhoScaled->GetVal());
    //fHistRhoScaledvsNtrackvsMult->Fill(Ntracks, fOutRhoScaled->GetVal(),multA+multC);
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
  TFile* file = TFile::Open(path);

  // does file exist
  if(!file || file->IsZombie()) {
    ::Error("StRhoBase", "Could not open scale function file");
    return 0;
  }

  // get scale function from file
  TF1* sfunc = dynamic_cast<TF1*>(file->Get(name));
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
