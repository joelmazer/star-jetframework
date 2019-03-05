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
#include "TVector3.h"

// STAR includes
///#include "StThreeVectorF.hh"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"

// jet-framework STAR includes
#include "StJet.h"
#include "StRhoParameter.h"
#include "StJetMakerTask.h"
#include "StJetFrameworkPicoBase.h"

// STAR centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StRhoSparse)

//________________________________________________________________________
StRhoSparse::StRhoSparse() : 
  StRhoBase("StRhoSparse"),
  fNExclLeadJets(0),
  fCreateHisto(0),
  fRhoCMS(0),
  fHistOccCorrvsCent(0),
  fHistOccCorrvsMult(0),
  fHistMultvsUnCorrRho(0),
  fHistMultvsCorrRho(0),
  mOutName(""),
  fRhoSparseMakerName("")
{
  // Constructor.
}

//________________________________________________________________________
StRhoSparse::StRhoSparse(const char *name, Bool_t histo, const char *outName, const char *jetMakerName) :
//  StRhoBase(name, histo),
  StRhoBase(name, histo, jetMakerName), // FIXME?
  fNExclLeadJets(0),
  fCreateHisto(histo),
  fRhoCMS(0),
  fHistOccCorrvsCent(0),
  fHistOccCorrvsMult(0),
  fHistMultvsUnCorrRho(0),
  fHistMultvsCorrRho(0),
  mOutName(outName),
  fRhoSparseMakerName(name)
{
  // Constructor.
  if(!name) return;
  SetName(name);
}

//________________________________________________________________________
StRhoSparse::~StRhoSparse()
{ /*  */
  // destructor
  if(fHistOccCorrvsCent)   delete fHistOccCorrvsCent;
  if(fHistOccCorrvsMult)   delete fHistOccCorrvsMult;
  if(fHistMultvsUnCorrRho) delete fHistMultvsUnCorrRho;
  if(fHistMultvsCorrRho)   delete fHistMultvsCorrRho;
}

//________________________________________________________________________
Int_t StRhoSparse::Init()
{
  // nothing done - base class should take care of that
  StRhoBase::Init();

  DeclareHistograms();

  // Create user objects.
  fBGJets = new TClonesArray("StJet");

  return kStOk;
}

//________________________________________________________________________
Int_t StRhoSparse::Finish() {

  //  Write histos to file and close it.
  // "test.root"  
  if(mOutName!="") {
    //TFile *fout = new TFile("test.root", "UPDATE");
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
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
  // weird that this next line is needed to remove cppcheck warning
  delete fHistOccCorrvsCent;
  //fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "Occupancy correction vs centrality", 101, -1, 100, 2000, 0 , 2);
  fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "Occupancy correction vs centrality", 20, 0.0, 100, 200, 0, 2); //2000
  fHistOccCorrvsMult = new TH2F("OccCorrvsMult", "Occupancy correction vs multiplicity", 160, 0.0, 800, 200, 0, 2); //2000

  fHistMultvsUnCorrRho = new TH2F("fHistMultvsUncorrRho", "Multiplicity vs #rho", 160, 0., 800., 100, 0., 100.);
  fHistMultvsUnCorrRho->GetXaxis()->SetTitle("Charged track multiplicity");
  fHistMultvsUnCorrRho->GetYaxis()->SetTitle("#rho (GeV/c)/A");

  fHistMultvsCorrRho = new TH2F("fHistMultvsCorrRho", "Multiplicity vs corrected #rho", 160, 0., 800., 100, 0., 100.);
  fHistMultvsCorrRho->GetXaxis()->SetTitle("Charged track multiplicity");
  fHistMultvsCorrRho->GetYaxis()->SetTitle("corrected #rho (GeV/c)/A");
}

//________________________________________________________________________
void StRhoSparse::WriteHistograms() {
  // write histograms
  fHistOccCorrvsCent->Write();
  fHistOccCorrvsMult->Write();
  fHistMultvsUnCorrRho->Write();
  fHistMultvsCorrRho->Write();
}

//________________________________________________________________________
void StRhoSparse::Clear(Option_t *opt) {

}

//________________________________________________________________________
Bool_t StRhoSparse::IsJetOverlapping(StJet* jet1, StJet* jet2) {
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i) {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j) {
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
  if(fOutRhoScaled)  fOutRhoScaled->SetVal(0);

  if(fDebugLevel == 1) cout<<"fJetMakerName = "<<fJetMakerName<<"  fJetBGMakerName = "<<fJetBGMakerName<<endl;

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

  // get vertex 3 vector and declare variables
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // z-vertex cut - per the Aj analysis (-40,40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // get JetMaker of background jets ("JetMakerBG")
  JetMakerBG = static_cast<StJetMakerTask*>(GetMaker(fJetBGMakerName));
  if(!JetMakerBG) { return kStWarn; }

  // if we have background JetMaker, get jet collection associated with it
  fBGJets =  static_cast<TClonesArray*>(JetMakerBG->GetJets());
  if(!fBGJets) { return kStWarn; }

  // get Signal JetMaker  - ("JetMaker")
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  if(!JetMaker) { return kStWarn; }

  // if we have signal JetMaker, get jet collection associated with it
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) { return kStWarn; }

  // # of jets
  const Int_t Njets = fBGJets->GetEntries();
  const Int_t NjetsSig = fJets->GetEntries();

  // get run # for centrality correction
  Int_t RunId = mPicoEvent->runId();
  Float_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Float_t fZDCCoincidenceRate = mPicoEvent->ZDCx();

  // Centrality correction calculation
  // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning
  int grefMult = mPicoEvent->grefMult();
  Int_t centbin, cent16;
  Double_t refCorr2;

  // check for AuAu analyses
  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorr->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }

    // get centrality bin: either 0-7 or 0-15
    cent16 = grefmultCorr->getCentralityBin16();

    // re-order binning to be from central -> peripheral
    centbin = GetCentBin(cent16, 16);  // 0-16

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fZDCCoincidenceRate, 2); }

  } else { // for pp
    centbin = 0, cent16 = 0, refCorr2 = 0.0;
  } 

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStWarn; // maybe kStOk; - this is for lowest multiplicity events 80%+ centrality, cut on them
  Double_t fCent = centbin * 5.0;

  // cut on centrality for analysis before doing anything
  if(fDebugLevel == 3) { cout<<"fCentralitySelectionCut = "<<fCentralitySelectionCut<<endl; }
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // to limit filling unused entries in sparse, only fill for certain centrality ranges
  // ranges can be different than functional cent bin setter
  Int_t cbin = -1;
  if (centbin>-1 && centbin < 2)    cbin = 1; // 0-10%
  else if (centbin>1 && centbin<4)  cbin = 2; // 10-20%
  else if (centbin>3 && centbin<6)  cbin = 3; // 20-30%
  else if (centbin>5 && centbin<10) cbin = 4; // 30-50%
  else if (centbin>9 && centbin<16) cbin = 5; // 50-80%
  else cbin = -99;
  // ============================ end of CENTRALITY ============================== //

  // initialize leading jet arrays
  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  // leading jet exclusion
  if(fNExclLeadJets > 0) {
    for(Int_t ij = 0; ij < Njets; ++ij) {
      StJet *jet = static_cast<StJet*>(fBGJets->At(ij));
      if(!jet) { continue; } 
      //if (!AcceptJet(jet)) continue;

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
    if(fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = 0;
    }
  }

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;
  Double_t TotaljetArea=0;
  Double_t TotaljetAreaPhys=0;

  // push all jets within selected acceptance into stack
  for(Int_t iJets = 0; iJets < Njets; ++iJets) {
    // excluding lead jets
    if(iJets == maxJetIds[0] || iJets == maxJetIds[1]) continue;

    // get background jets
    StJet *jet = static_cast<StJet*>(fBGJets->At(iJets));
    if(!jet) { continue; } 

    // add total jet area
    TotaljetArea+=jet->Area();
    if(jet->Pt() > 0.1) {
      TotaljetAreaPhys+=jet->Area();
    }

    //if (!AcceptJet(jet)) continue;

    // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    if(fJets) {
      for(Int_t j = 0; j < NjetsSig; j++) {
        // get signal jets
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

  // calculate occupancy correction
  Double_t OccCorr = 0.0;
  if(TotaljetArea > 0) OccCorr = TotaljetAreaPhys / TotaljetArea;

  if(fDebugLevel == 2) cout<<"NjetsSig = "<<NjetsSig<<"  NBGjets = "<<Njets<<"  TotaljetArea = "<<TotaljetArea<<"  OccCorr = "<<OccCorr<<endl;
 
  // fill histos
  if(fCreateHisto) fHistOccCorrvsCent->Fill(fCent, OccCorr);
  if(fCreateHisto) fHistOccCorrvsMult->Fill(refCorr2, OccCorr);

  // when we have accepted jets
  if(NjetAcc > 0) {
    // find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);
    Double_t uncorrho = TMath::Median(NjetAcc, rhovec);
    if(fCreateHisto) fHistMultvsUnCorrRho->Fill(refCorr2, uncorrho);

    // correct rho
    if(fRhoCMS){
      rho = rho * OccCorr;
    }
    fOutRho->SetVal(rho);

    // fill histo
    if(fCreateHisto) fHistMultvsCorrRho->Fill(refCorr2, rho);

    // if we need to scale rho from charged -> ch+ne
    if(fOutRhoScaled) {
      //Double_t rhoScaled = rho * GetScaleFactor(fCent); // fCent in 5% bins
      Double_t rhoScaled = rho * 1.0;
      fOutRhoScaled->SetVal(rhoScaled);
    }
  }

  return kStOk;
} 
