// $Id$
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".apppend("_Scaled").
//
#include "StRhoSparse.h"

#include <TClonesArray.h>
#include <TMath.h>
#include "TH2.h"
#include "TH2F.h"

// STAR includes
#include "StThreeVectorF.hh"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"

#include "StJet.h"
#include "StRhoParameter.h"
#include "StJetMakerTask.h"

ClassImp(StRhoSparse)

//________________________________________________________________________
StRhoSparse::StRhoSparse() : 
  StRhoBase("StRhoSparse"),
  fNExclLeadJets(0),
  fCreateHisto(0),
  fRhoCMS(0),
  fHistOccCorrvsCent(0),
  mOutName("")
{
  // Constructor.
}

//________________________________________________________________________
StRhoSparse::StRhoSparse(const char *name, Bool_t histo, const char *outName) :
//  StRhoBase(name, histo),
  StRhoBase(name, histo, "JetsMakerBG"),
  fNExclLeadJets(0),
  fCreateHisto(histo),
  fRhoCMS(0),
  fHistOccCorrvsCent(0),
  mOutName(outName)
{
  // Constructor.
  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StRhoSparse::~StRhoSparse()
{ /*  */
  // destructor
  delete fHistOccCorrvsCent;
}

//________________________________________________________________________
Int_t StRhoSparse::Init()
{
  // nothing done - base class should take care of that
  StRhoBase::Init();

  DeclareHistograms();

  //Create user objects.
  fBGJets = new TClonesArray("StJet");

  return kStOk;
}

//________________________________________________________________________
Int_t StRhoSparse::Finish() {

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile("test2.root", "UPDATE");
    fout->cd();
    fout->mkdir("RhoSparse");
    fout->cd("RhoSparse");
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

  return kStOK;
}

//________________________________________________________________________
void StRhoSparse::DeclareHistograms() {
    // declare histograms
  //fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "OccCorrvsCent", 101, -1, 100, 2000, 0 , 2);
  fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "OccCorrvsCent", 20, 0.0, 100, 2000, 0 , 2);
}

//________________________________________________________________________
void StRhoSparse::WriteHistograms() {
  // write histograms
  fHistOccCorrvsCent->Write();
}

//________________________________________________________________________
void StRhoSparse::Clear(Option_t *opt) {

}

//________________________________________________________________________
Bool_t StRhoSparse::IsJetOverlapping(StJet* jet1, StJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
  {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
    {
      Int_t jet2Track = jet2->TrackAt(j);
      if (jet1Track == jet2Track)
        return kTRUE;
    }
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t StRhoSparse::IsJetSignal(StJet* jet)
{
  if(jet->Pt()>5){
      return kTRUE;
  }else{
    return kFALSE;
  }
}


//________________________________________________________________________
Int_t StRhoSparse::Make() 
{
  // Run the analysis.
  fOutRho->SetVal(0);
  if (fOutRhoScaled)  fOutRhoScaled->SetVal(0);

  // get PicoDstMaker
  mPicoDstMaker = (StPicoDstMaker*)GetMaker("picoDst");
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStFatal;
  }

  // construct PicoDst object from maker
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStFatal;
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

  // zvertex cut - per the Aj analysis
  if((1.0*TMath::Abs(zVtx)) > 40) return kStOk;   //kStWarn; //kStFatal;

  // get JetMaker - (fJetMakerName) - TODO set name globally
  StJetMakerTask *JetMakerBG = (StJetMakerTask*)GetMaker("JetMakerBG");
  if(!JetMakerBG) { return kStWarn; }

  // if we have background JetMaker, get jet collection associated with it
  if(JetMakerBG) {
    fBGJets =  JetMakerBG->GetJets();
  }
  if(!fBGJets) { return kStWarn; }

  // # of bg jets
  const Int_t Njets = fBGJets->GetEntries();

  // get Signal JetMaker    - (fJetMakerName)
  JetMaker = (StJetMakerTask*)GetMaker("JetMaker");
  if(!JetMaker) { return kStWarn; }

  // if we have signal JetMaker, get jet collection associated with it
  fJets = JetMaker->GetJets();
  if(!fJets) { return kStWarn; }

  Int_t NjetsSig = 0;
  if (fJets) NjetsSig = fJets->GetEntries();

  // initialize leading jet arrays
  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  // leading jet exclusion
  if (fNExclLeadJets > 0) {
    for (Int_t ij = 0; ij < Njets; ++ij) {
      StJet *jet = static_cast<StJet*>(fBGJets->At(ij));
      if (!jet) { continue; } 
      //if (!AcceptJet(jet)) continue;

      if (jet->Pt() > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jet->Pt();
	maxJetIds[0] = ij;
      } else if (jet->Pt() > maxJetPts[1]) {
	maxJetPts[1] = jet->Pt();
	maxJetIds[1] = ij;
      }
    }
    if (fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = 0;
    }
  }

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;
  Double_t TotaljetArea=0;
  Double_t TotaljetAreaPhys=0;

  // push all jets within selected acceptance into stack
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
    // exlcuding lead jets
    if (iJets == maxJetIds[0] || iJets == maxJetIds[1]) continue;

    StJet *jet = static_cast<StJet*>(fBGJets->At(iJets));
    if (!jet) {  continue; } 

    TotaljetArea+=jet->Area();
    if(jet->Pt()>0.1){
      TotaljetAreaPhys+=jet->Area();
    }

    //if (!AcceptJet(jet)) continue;

   // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    if (fJets) {
      for(Int_t j=0;j<NjetsSig;j++) {
        StJet* signalJet = static_cast<StJet*>(fJets->At(j)); // GetAcceptJet(j)
        if(!signalJet) continue;
        if(!IsJetSignal(signalJet)) continue;
	  
        if(IsJetOverlapping(signalJet, jet)) {
          isOverlapping = kTRUE;
          break;
        }
      }
    }

    if(isOverlapping) continue;

    if(jet->Pt()>0.1){
      rhovec[NjetAcc] = jet->Pt() / jet->Area();
      ++NjetAcc;
    }
  }

  // fCent just initialized for test
  double fCent = 1.0;

  Double_t OccCorr=0.0;
  if(TotaljetArea>0) OccCorr=TotaljetAreaPhys/TotaljetArea;

  cout<<"NjetsSig = "<<NjetsSig<<"  NBGjets = "<<Njets<<"  TotaljetArea = "<<TotaljetArea<<"  OccCorr = "<<OccCorr<<endl;
 
  if (fCreateHisto) fHistOccCorrvsCent->Fill(fCent, OccCorr);

  // when we have accepted jets
  if (NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);

    if(fRhoCMS){
      rho = rho * OccCorr;
    }

    fOutRho->SetVal(rho);

    // if we need to scale rho from charged -> ch+ne
    if (fOutRhoScaled) {
      Double_t rhoScaled = rho * GetScaleFactor(fCent);
      fOutRhoScaled->SetVal(rhoScaled);
    }
  }

  return kTRUE;
} 
