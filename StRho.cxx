// $Id$
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fOutRhoName".Apppend("_Scaled").
//

#include "StRho.h"

// ROOT includes
#include <TClonesArray.h>
#include <TMath.h>
#include "TH2.h"
#include "TH2F.h"

class TH2;
class TH2F;

// JetFramework includes
#include "StJet.h"
#include "StRhoParameter.h"
#include "StJetMakerTask.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"

// STAR centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

//class StJetMakerTask;

ClassImp(StRho)

//________________________________________________________________________
StRho::StRho() : StRhoBase("")
{
  fNExclLeadJets = 0;
  fJets = 0x0;
  fHistMultvsRho = 0x0;
  mOutName = ""; 
  fJetMakerName = "";
  fRhoMakerName = "";
}

//________________________________________________________________________
StRho::StRho(const char *name, Bool_t histo, const char *outName, const char *jetMakerName) :
  StRhoBase(name, histo, jetMakerName)
{
  fNExclLeadJets = 0;
  fJets = 0x0;
  fHistMultvsRho = 0x0;
  mOutName = outName;
  fJetMakerName = jetMakerName;
  fRhoMakerName = name;

  // Constructor.
  if (!name) return;
  SetName(name);
}

//________________________________________________________________________
StRho::~StRho()
{ /*  */
  // destructor
  if(fHistMultvsRho) delete fHistMultvsRho;

}

//________________________________________________________________________
Int_t StRho::Init()
{
  // nothing done - base class should take care of that
  // this in effect inherits from StJetFrameworkPicoBase - check it out!
  StRhoBase::Init();

  DeclareHistograms();

  // Create user objects.
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
    ///StRhoBase::WriteHistograms();
    WriteHistograms();
    fout->cd();
    fout->Write();
    fout->Close();
  }

/*   ===== test ======
  if(mOutName!="") {
    cout<<"checking output file in StRhoMaker::Finish().."<<endl;
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
    fout->cd();
    fout->Close();
  }
*/

  //cout<<"End of StRho::Finish"<<endl;
  return kStOK;
}

//________________________________________________________________________
void StRho::DeclareHistograms() {
    // declare histograms
    // weird that this next line is needed to remove cppcheck warning
    delete fHistMultvsRho;
    // mult was 150, 0, 1500
    fHistMultvsRho = new TH2F("fHistMultvsRho", "fHistMultvsRho", 160, 0., 800., 100, 0., 100.);
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
//  StRhoBase::Clear();

//  fJets->Clear();
}

//________________________________________________________________________
Int_t StRho::Make() 
{
  // Run the analysis - for each event

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
  StThreeVectorF mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();

  // z-vertex cut - per the Aj analysis (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // if we have JetMaker, get jet collection associated with it
  if(JetMaker) {
    fJets =  JetMaker->GetJets();
    //fJets->SetName("BGJetsRho");  // name is set by Maker who created it
  }
  if(!fJets) return kStWarn;

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

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }

  // get event multiplicity - TODO is this correct? same as that used for centrality
  //const int multiplicity = mPicoDst->numberOfTracks(); // this is total tracks not multiplicity
  const int multiplicity = refCorr2;

  // ============================ end of CENTRALITY ============================== //

  // initialize Rho and scaled Rho
  fOutRho->SetVal(0);
  if(fOutRhoScaled) fOutRhoScaled->SetVal(0);

  // get number of jets, initialize arrays
  const Int_t Njets   = fJets->GetEntries();
  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  // exclude leading jets
  if(fNExclLeadJets > 0) {
    // loop over jets
    for(Int_t ij = 0; ij < Njets; ++ij) {
      StJet *jet = static_cast<StJet*>(fJets->At(ij));
      if(!jet) { continue; } 

      // NEED TO CHECK DEFAULTS // FIXME
      //if (!AcceptJet(jet)) continue;
      // get some jet parameters
      double jetPt = jet->Pt();
      // some threshold cuts for tests
      if(jetPt < 0) continue;

      // get ID and pt of the leading and sub-leading jet   
      if(jetPt > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jetPt;
	maxJetIds[0] = ij;
      } else if (jetPt > maxJetPts[1]) {
	maxJetPts[1] = jetPt;
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

    // pointer to jet
    StJet *jet = static_cast<StJet*>(fJets->At(iJets));
    if(!jet) { continue; } 

    // NEED TO CHECK FOR DEFAULTS - cuts are done at the jet finder level
    //if(!AcceptJet(jet)) continue; //FIXME
    // get some get parameters
    double jetPt = jet->Pt();
    double jetArea = jet->Area();
    // some threshold cuts for tests
    if(jetPt < 0) continue;
 
    rhovec[NjetAcc] = jetPt / jetArea;
    ++NjetAcc;
  }

  // when we have accepted Jets - calculate and set rho
  if(NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);
    fOutRho->SetVal(rho);

    // fill histo
    fHistMultvsRho->Fill(multiplicity, rho);

    // if we want scaled Rho from charged -> ch+ne
    if(fOutRhoScaled) {
      //Double_t rhoScaled = rho * GetScaleFactor(fCent); //Don't need yet, fCent is in 5% bins
      Double_t rhoScaled = rho * 1.0;
      fOutRhoScaled->SetVal(rhoScaled);
    }
  }

  StRhoBase::FillHistograms();

  return kStOk;
} 
