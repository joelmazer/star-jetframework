// $Id$
// Adapation from the AliROOT class AliAnalysisTalkRho.cxx
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".Apppend("_Scaled").
//
// original
// Authors: R.Reed, S.Aiola

#include "StRho.h"

// ROOT includes
#include <TClonesArray.h>
#include <TMath.h>
#include "TH2.h"
#include "TH2F.h"

class TH2;
class TH2F;

// STAR includes
#include "StJet.h"
#include "StRhoParameter.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
///#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

#include "StJetMakerTask.h"

//class TH2;
class StJetMakerTask;

ClassImp(StRho)

//________________________________________________________________________
StRho::StRho() : 
  StRhoBase("StRho"),
  fNExclLeadJets(0),
  fJets(0),
  mPicoDstMaker(0),
  mPicoDst(0),
  mPicoEvent(0),
  JetMaker(0),
  fHistMultvsRho(0),
  mOutName(""), 
  fJetMakerName(""),
  fRhoMakerName("")
{
  // Standard constructor.
}

//________________________________________________________________________
StRho::StRho(const char *name, Bool_t histo, const char *outName, const char *jetMakerName) :
  StRhoBase(name, histo, jetMakerName),
  fNExclLeadJets(0),
  fJets(0),
  mPicoDstMaker(0),
  mPicoDst(0),
  mPicoEvent(0),
  JetMaker(0),
  fHistMultvsRho(0),
  mOutName(outName),
  fJetMakerName(jetMakerName),
  fRhoMakerName(name)
{
  // Constructor.
  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StRho::~StRho()
{ /*  */
  // destructor
  delete fHistMultvsRho;

  //fJets->Clear(); delete fJets;
}

//________________________________________________________________________
Int_t StRho::Init()
{
  // nothing done - base class should take care of that
  StRhoBase::Init();

  DeclareHistograms();

  //Create user objects.
  fJets = new TClonesArray("StJet");
  //fJets->SetName(fJetsName);      

  return kStOk;
}

//________________________________________________________________________
Int_t StRho::Finish() {
  //  Summarize the run.
  //cout << "StRho::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(fRhoMakerName);
    fout->cd(fRhoMakerName);
    StRhoBase::WriteHistograms();
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  //cout<<"End of StRho::Finish"<<endl;
  return kStOK;
}

//________________________________________________________________________
void StRho::DeclareHistograms() {
    // declare histograms
    fHistMultvsRho = new TH2F("fHistMultvsRho", "fHistMultvsRho", 150, 0., 1500., 100, 0., 100.);
    fHistMultvsRho->GetXaxis()->SetTitle("Charged track multiplicity");
    fHistMultvsRho->GetYaxis()->SetTitle("#rho (GeV/c)/A");

}

//________________________________________________________________________
void StRho::WriteHistograms() {
  // write histograms
  fHistMultvsRho->Write();
}

//________________________________________________________________________
void StRho::Clear(Option_t *opt) {
  StRhoBase::Clear();

  fJets->Clear();
}

//________________________________________________________________________
Int_t StRho::Make() 
{
  // Run the analysis.

  // get PicoDstMaker
  mPicoDstMaker = (StPicoDstMaker*)GetMaker("picoDst");
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    //return kStWarn;
    return kStFatal;
  }

  // construct PicoDst object from maker
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    //return kStWarn;
    return kStFatal;
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
  if((1.0*TMath::Abs(zVtx)) > 40) return kStOk;   //kStWarn; //kStFatal;

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
    fJets->SetName("BGJetsRho");
  }
  if(!fJets) return kStWarn; //kStFatal;

  // get event multiplicity -TODO should define this differently
  const int multiplicity = mPicoDst->numberOfTracks();

  // initialize Rho and scaled Rho
  fOutRho->SetVal(0);
  if(fOutRhoScaled) fOutRhoScaled->SetVal(0);

  // get number of jets
  const Int_t Njets   = fJets->GetEntries();

  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  // exclude leading jets
  if(fNExclLeadJets > 0) {
    for (Int_t ij = 0; ij < Njets; ++ij) {
      StJet *jet = static_cast<StJet*>(fJets->At(ij));
      if(!jet) {
	//Form("%s: Could not receive jet %d", GetName(), ij); //FIXME
	continue;
      } 

      // NEED TO CHECK DEFAULTS // FIXME
      //if (!AcceptJet(jet)) continue;
      // get some jet parameters
      double jetpt = jet->Pt();
      double jetEta = jet->Eta();
      double jetPhi = jet->Phi();
      // some threshold cuts for tests
      if(jetpt < 0) continue;

      // get ID and pt of the leading and sub-leading jet   
      if(jet->Pt() > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jet->Pt();
	maxJetIds[0] = ij;
      } else if (jet->Pt() > maxJetPts[1]) {
	maxJetPts[1] = jet->Pt();
	maxJetIds[1] = ij;
      }
    }
    // only set to remove leading jet 
    if(fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = 0;
    }
  }

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;

  // push all jets within selected acceptance into stack
  for(Int_t iJets = 0; iJets < Njets; ++iJets) {

    // excluding lead jets
    if(iJets == maxJetIds[0] || iJets == maxJetIds[1])
      continue;

    // pointer to jets
    StJet *jet = static_cast<StJet*>(fJets->At(iJets));
    if(!jet) {
      continue;
    } 

    // NEED TO CHECK FOR DEFAULTS - cuts are done at the jet finder level
    //if(!AcceptJet(jet)) continue; //FIXME
    // get some get parameters
    double jetpt = jet->Pt();
    double jetEta = jet->Eta();
    double jetPhi = jet->Phi();
    // some threshold cuts for tests
    if(jetpt < 0) continue;
 
    rhovec[NjetAcc] = jet->Pt() / jet->Area();
    ++NjetAcc;
  }

  // when we have accepted Jets
  if(NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);
    fOutRho->SetVal(rho);

    // print statement
    //cout<<"Rho = "<<rho<<endl;
    fHistMultvsRho->Fill(multiplicity, rho);

    // if we want scaled Rho
    if(fOutRhoScaled) {
      //Double_t rhoScaled = rho * GetScaleFactor(fCent); //FIXME
      Double_t rhoScaled = rho * 1.0;
      fOutRhoScaled->SetVal(rhoScaled);
    }
  }

  StRhoBase::FillHistograms();
  //fHistMultvsRho->Fill(multiplicity, rho);

  return kStOk;
} 
