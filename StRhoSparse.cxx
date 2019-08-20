// $Id$
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".apppend("_Scaled").
//
#include "StRhoSparse.h"

// ROOT includes
#include <TClonesArray.h>
#include <TMath.h>
#include "TH2.h"
#include "TH2F.h"
#include "TVector3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"

// jet-framework STAR includes
#include "StJet.h"
#include "StRhoParameter.h"
#include "StJetMakerTask.h"
#include "StJetFrameworkPicoBase.h"
#include "StCentMaker.h"

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

  // declare histograms
  DeclareHistograms();

  // Create user objects.
  fBGJets = new TClonesArray("StJet");

  return kStOk;
}
//
// write histos to file and close it
//________________________________________________________________________
Int_t StRhoSparse::Finish() {
  if(mOutName!="") {
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
//
// Function: declare histograms and user output objects
//________________________________________________________________________
void StRhoSparse::DeclareHistograms() {
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
//
// Function to write histograms
//________________________________________________________________________
void StRhoSparse::WriteHistograms() {
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
  // jet 1 tracks:
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i) {
    Int_t jet1Track = jet1->TrackAt(i);

    // jet 2 tracks:
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
  if(jet->Pt() > 5) {
      return kTRUE;
  } else {
    return kFALSE;
  }
}

//
// Functions that runs over the analysis for each event
//________________________________________________________________________
Int_t StRhoSparse::Make() 
{
  // re-initialize the Rho objects
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
  int centbin = mCentMaker->GetRef16();
  double refCorr2 = mCentMaker->GetRefCorr2();
  double fCent = mCentMaker->GetCentScaled(); // integer bins scaled up by 5% per
  //double refCorr = mCentMaker->GetCorrectedMultiplicity(refMult, zVtx, zdcCoincidenceRate, 0); // example usage
  // for pp analyses:    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;

/*
  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStOk; // this is for lowest multiplicity events 80%+ centrality, cut on them 
*/

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
      // get jet pointer
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

    // get background jet pointer
    StJet *jet = static_cast<StJet*>(fBGJets->At(iJets));
    if(!jet) { continue; } 

    // add total jet area
    TotaljetArea += jet->Area();
    if(jet->Pt() > 0.1) {
      TotaljetAreaPhys += jet->Area();
    }

    //if (!AcceptJet(jet)) continue;

    // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    if(fJets) {
      for(Int_t j = 0; j < NjetsSig; j++) {
        // get signal jets
        StJet *signalJet = static_cast<StJet*>(fJets->At(j)); // GetAcceptJet(j)
        if(!signalJet) continue;
        if(!IsJetSignal(signalJet)) continue;
	  
        if(IsJetOverlapping(signalJet, jet)) {
          isOverlapping = kTRUE;
          break;
        }
      }
    }

    if(isOverlapping) continue;

    if(jet->Pt() > 0.1){
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
