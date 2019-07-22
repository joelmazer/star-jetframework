// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// Centrality QA
//
// ################################################################

#include "StCentMaker.h"
#include "StRoot/StarRoot/StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TFile.h"
#include "TParameter.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"

// jet-framework includes
#include "StJetFrameworkPicoBase.h"

// old file kept
#include "StPicoConstants.h"

// extra includes
#include "StJetPicoDefinitions.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StCentMaker)

//______________________________________________________________________________
StCentMaker::StCentMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName = "", bool mDoComments = kFALSE)
  : StJetFrameworkPicoBase(name)
{
  fDebugLevel = 0;
  fRunFlag = 0;
  doppAnalysis = kFALSE;
  fCentralityDef = 4; // see StJetFrameworkPicoBase::fCentralityDefEnum
  doUseBBCCoincidenceRate = kFALSE; // kFALSE = use ZDC
  fMaxEventTrackPt = 30.0;
  fMaxEventTowerE = 1000.0; // 30.0
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  mOutName = outName;
  doRejectBadRuns = kFALSE;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  //////////////////////////////
  kgrefMult = -99;
  krefMult = -99;
  krefCorr2 = 0.0;
  kref16 = -99; kref9 = -99;
  kcent9 = -99; kcent16 = -99;
  kCentralityScaled = -99;

  fEmcTriggerEventType = 0; 
  fMBEventType = 2;
  doComments = mDoComments;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;
}

//_____________________________________________________________________________
StCentMaker::~StCentMaker()
{ /*  */
  // destructor
  if(hEventZVertex)         delete hEventZVertex;
  if(hCentrality)           delete hCentrality;
  if(hMultiplicity)         delete hMultiplicity;
  if(fHistEventSelectionQA) delete fHistEventSelectionQA;
}

//_____________________________________________________________________________
Int_t StCentMaker::Init() {
  // initialize the histograms
  DeclareHistograms();

  // switch on Run Flag specifically requested for given run period for centrality definition setup
  switch(fRunFlag) {
    case StJetFrameworkPicoBase::Run11_pp500 : // Run11: 500 GeV pp
        break;

    case StJetFrameworkPicoBase::Run12_pp200 : // Run12: 200 GeV pp
        break;

    case StJetFrameworkPicoBase::Run12_pp500 : // Run12: 500 GeV pp
        break;

    case StJetFrameworkPicoBase::Run13_pp510 : // Run13: 510 (500) GeV pp
        break;
  
    case StJetFrameworkPicoBase::Run14_AuAu200 : // Run14: 200 GeV AuAu
        switch(fCentralityDef) {
          case StJetFrameworkPicoBase::kgrefmult :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P17id_VpdMB30();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30 :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P16id :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
              break;
          default:
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
        }
        break;

    case StJetFrameworkPicoBase::Run15_pp200 : // Run15: 200 GeV pp
        break;

    case StJetFrameworkPicoBase::Run16_AuAu200 : // Run16: 200 GeV AuAu
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

    case StJetFrameworkPicoBase::Run17_pp510 : // Run17: 510 (500) GeV pp
        // this is the default for Run17 pp - don't set anything for pp
        break;

    default :
        // for any non specified Run, or any pp Run, this default is set and as such not used, as parameters are set to 0
        cout<<"DEFAULT centrality setup!"<<endl;
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  }

  return kStOK;
}

//____________________________________________________________________________
Int_t StCentMaker::Finish() { 
  cout << "StCentMaker::Finish()\n";

  //  Write histos to file and close it.
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    fout->mkdir(GetName());
    fout->cd(GetName());
    WriteHistograms();
   
    fout->cd();
    fout->Write();
    fout->Close();
  }

  cout<<"End of StCentMaker::Finish"<<endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

//_____________________________________________________________________________
void StCentMaker::DeclareHistograms() {
  // setup for centrality binning
  int nHistCentBins = 20;

  // QA histos
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 200, -100., 100.);
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", 160, 0, 800);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);

  // Switch on Sumw2 for all histos - (except profiles)
  SetSumw2();
}

//
// write histograms
//_____________________________________________________________________________
void StCentMaker::WriteHistograms() {
  // default histos
  hEventZVertex->Write();
  hCentrality->Write();
  hMultiplicity->Write();

  // QA histos
  fHistEventSelectionQA->Write(); 
}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StCentMaker::Clear(Option_t *opt) {

}
// 
//  This method is called every event.
//_____________________________________________________________________________
Int_t StCentMaker::Make() {
  // reset global parameters
  kgrefMult = -99, krefMult = -99, kref9 = -99, kref16 = -99, kcent9 = -99, kcent16 = -99;
  kCentralityScaled = -99., krefCorr2 = 0.0;

  // get PicoDstMaker 
  mPicoDstMaker = static_cast<StPicoDstMaker*>(GetMaker("picoDst"));
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  // get PicoDst object from maker
  mPicoDst = static_cast<StPicoDst*>(mPicoDstMaker->picoDst());
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // get pointer to PicoEvent 
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

  // cut event on max tower E > 30.0 GeV
  //if(GetMaxTowerE() > fMaxEventTowerE) return kStOK;

  // get vertex 3-vector and z-vertex component
  TVector3 mVertex = mPicoEvent->primaryVertex();
  double zVtx = mVertex.z();
 
  // commented out for now, but centrality was configured with 30 cm z-vertex data 
  // Z-vertex cut - the Aj analysis cut on (-40, 40) for reference
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  hEventZVertex->Fill(zVtx);

  // coincidence rates used for corrected multiplicity calculation
  double fBBCCoincidenceRate = mPicoEvent->BBCx();
  double fZDCCoincidenceRate = mPicoEvent->ZDCx();

  // ============================ CENTRALITY ============================== //
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/
  Int_t centbin;
  kgrefMult = mPicoEvent->grefMult(); // grefMult
  krefMult = mPicoEvent->refMult();   // refMult

  // for non-pp analyses
  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(fRunNumber);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(kgrefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorr->initEvent(kgrefMult, zVtx, fZDCCoincidenceRate); }
//    if(grefmultCorr->isBadRun(fRunNumber)) cout << "Run is bad" << endl; 
//    if(grefmultCorr->isIndexOk()) cout << "Index Ok" << endl;
//    if(grefmultCorr->isZvertexOk()) cout << "Zvertex Ok" << endl;
//    if(grefmultCorr->isRefMultOk()) cout << "RefMult Ok" << endl;

    // get centrality bin: either 0-7 or 0-15 - but decreasing % with increasing value (poor practice, use below)
    kcent16 = grefmultCorr->getCentralityBin16();
    kcent9 = grefmultCorr->getCentralityBin9();

    // re-order binning to be from central -> peripheral
    kref9 = GetCentBin(kcent9, 9);
    kref16 = GetCentBin(kcent16, 16);
    centbin = GetCentBin(kcent16, 16);  // 0-16

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { krefCorr2 = grefmultCorr->getRefMultCorr(kgrefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ krefCorr2 = grefmultCorr->getRefMultCorr(kgrefMult, zVtx, fZDCCoincidenceRate, 2); }

    // calculate corrected multiplicity: 
    // Double_t getRefMultCorr(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag=2) const ;
    // flag=0:  Luminosity only
    // flag=1:  z-vertex only
    // flag=2:  full correction (default)
    //
    //grefmultCorr->isCentralityOk(cent16)
  } else { // for pp
    centbin = 0, kcent9 = 0, kcent16 = 0;
    krefCorr2 = 0.0, kref9 = 0, kref16 = 0;
  }

  // cut on unset centrality, > 80%
  if(kcent16 == -1) return kStOk; // this is for the lowest multiplicity events 80%+ centrality, cut on them here

  // multiplicity histogram
  hMultiplicity->Fill(krefCorr2);

  // centrality histogram
  kCentralityScaled = centbin*5.0;
  hCentrality->Fill(kCentralityScaled);

  // debug
  if(fDebugLevel == kDebugCentrality) { if(centbin > 15) cout<<"centbin = "<<centbin<<"  mult = "<<krefCorr2<<"   kCentralityScaled = "<<kCentralityScaled<<"  cent16 = "<<kcent16<<endl; }

  // ========================= Trigger Info =============================== //
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  return kStOK;
}

//
// Set the bin errors on histograms, set sum weights
// __________________________________________________________________________________
void StCentMaker::SetSumw2() {
  hEventZVertex->Sumw2();
  hCentrality->Sumw2();
  hMultiplicity->Sumw2();
  fHistEventSelectionQA->Sumw2();
}
//
// Get corrected multiplicity for different correction flags
// - luminosity only, z-vtx only, full correction (default)
//_______________________________________________________________________________________________
Double_t StCentMaker::GetCorrectedMultiplicity(const UShort_t RefMult, const Double_t z,
const Double_t zdcCoincidenceRate, const UInt_t flag) {
    // Event-by-event initialization. Call this function event-by-event
    //   * Default ZDC coincidence rate = 0 to make the function backward compatible 
    //   --> i.e. no correction will be applied unless users set the values for 200 GeV  

    //   * Comment for luminosity (zdc coincidence rate) correction
    //     - Luminosity correction is only valid for 200 GeV
    //     - The default argument is 0 for zdc coincidence rate in initEvent() function, see header StRefMultCorr.h,
    //       so that you can still use previous initEvent() function like
    //         void StRefMultCorr::initEvent(refmult, vz);   
    //       without specifying zdc coincidence rate for lower beam energies
    //     - (Very important) You should use BBC coincidence rate for refmult2

    // Corrected multiplity
    // flag=0:  Luminosity only
    // flag=1:  z-vertex only
    // flag=2:  full correction (default)

    double refCorr = -99.;
    if(flag == 0) refCorr = grefmultCorr->getRefMultCorr(RefMult, z, zdcCoincidenceRate, flag);
    if(flag == 1) refCorr = grefmultCorr->getRefMultCorr(RefMult, z, zdcCoincidenceRate, flag);
    if(flag == 2) refCorr = grefmultCorr->getRefMultCorr(RefMult, z, zdcCoincidenceRate, flag);

    return refCorr;
}
