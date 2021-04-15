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
#include <THnSparse.h>
#include "TFile.h"
#include "TParameter.h"

#include <sstream>
#include <fstream>

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
  fMaxEventTowerEt = 1000.0; // 30.0
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  grefmultCorr = 0x0;
  refmultCorr = 0x0;
  refmult2Corr = 0x0;
  grefmultCorrMB5 = 0x0;
  grefmultCorrUtil = 0x0;
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
  kMB5toMB30ReWeight = 1.0;
  kReWeight = 1.0;

  fEmcTriggerEventType = 0; 
  fMBEventType = 2;
  doComments = mDoComments;
  mBaseMaker = 0x0;
  fAnalysisMakerName = name;

  fhnCentQA = 0x0; 
}

//_____________________________________________________________________________
StCentMaker::~StCentMaker()
{ /*  */
  // destructor
  if(hEventZVertex)         delete hEventZVertex;
  if(hCentrality)           delete hCentrality;
  if(hMultiplicity)         delete hMultiplicity;
  if(fHistEventSelectionQA) delete fHistEventSelectionQA;
  if(hMB5refCorr2ReWeight)  delete hMB5refCorr2ReWeight;
  if(hMB5refCorr2Raw)       delete hMB5refCorr2Raw;
  if(hMB5CentScaledReWeight)delete hMB5CentScaledReWeight;
  if(hMB5CentScaledRaw)     delete hMB5CentScaledRaw;
  if(hMB30refCorr2ReWeight) delete hMB30refCorr2ReWeight;
  if(hMB30refCorr2Raw)      delete hMB30refCorr2Raw;
  if(hMB30CentScaledReWeight) delete hMB30CentScaledReWeight;
  if(hMB30CentScaledRaw)    delete hMB30CentScaledRaw;

  for(int i=0; i<30; i++) {
    if(hMB5refCorr2[i])     delete hMB5refCorr2[i];
    if(hMB5refCorr2ReWt[i])    delete hMB5refCorr2ReWt[i];
    if(hMB5refCorr2ReWtInv[i]) delete hMB5refCorr2ReWtInv[i];
    if(hMB5grefMult[i])     delete hMB5grefMult[i];
    if(hMB5CentScaled[i])   delete hMB5CentScaled[i];
    if(hMB5refCorr2Weight[i])    delete hMB5refCorr2Weight[i];
    if(hMB5refCorr2WeightInv[i]) delete hMB5refCorr2WeightInv[i];

    if(hMB30refCorr2[i])    delete hMB30refCorr2[i];
    if(hMB30refCorr2ReWt[i])    delete hMB30refCorr2ReWt[i];
    if(hMB30refCorr2ReWtInv[i]) delete hMB30refCorr2ReWtInv[i];
    if(hMB30grefMult[i])    delete hMB30grefMult[i];
    if(hMB30CentScaled[i])  delete hMB30CentScaled[i];

    if(hHT2refCorr2[i])     delete hHT2refCorr2[i];
    if(hHT2refCorr2ReWt[i])    delete hHT2refCorr2ReWt[i];
    if(hHT2refCorr2ReWtInv[i]) delete hHT2refCorr2ReWtInv[i];
    if(hHT2grefMult[i])     delete hHT2grefMult[i];
    if(hHT2CentScaled[i])   delete hHT2CentScaled[i];

    if(hMB5onlyrefCorr2Raw[i])  delete hMB5onlyrefCorr2Raw[i];
    if(hMB30onlyrefCorr2Raw[i]) delete hMB30onlyrefCorr2Raw[i];
    if(hMB5MB30refCorr2Raw[i])  delete hMB5MB30refCorr2Raw[i];
  }
  if(hMB5onlyrefCorr2RawTotal)  delete hMB5onlyrefCorr2RawTotal;
  if(hMB30onlyrefCorr2RawTotal) delete hMB30onlyrefCorr2RawTotal;
  if(hMB5MB30refCorr2RawTotal)  delete hMB5MB30refCorr2RawTotal;

  if(fhnCentQA)             delete fhnCentQA;
}

//_____________________________________________________________________________
Int_t StCentMaker::Init() {
  // initialize the histograms
  DeclareHistograms();

  // ====================================================================================================================
  // centrality definitions loaded here
  //    - these are generated for mid luminosity runs
  //  
  //  Run 14: 
  //  	getgRefMultCorr_P18ih_VpdMB30_AllLumi() - used for all luminosities
  //    getgRefMultCorr_P18ih_VpdMB30()         - used for mid-lumi (older set of definitions)
  // ====================================================================================================================

  // switch on Run Flag specifically requested for given run period for centrality definition setup
  switch(fRunFlag) {
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
          case StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30_AllLumi();
              break;
          case StJetFrameworkPicoBase::kgrefmult_P16id :
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
              break;
          default:
              grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
        }
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

    default :
        // for any non specified Run, or any pp Run, this default is set and as such not used, as parameters are set to 0
        cout<<"DEFAULT centrality setup!"<<endl;
        grefmultCorr = CentralityMaker::instance()->getgRefMultCorr();
  }

  grefmultCorrMB5  = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30_AllLumi();

  //grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30_AllLumi();
  grefmultCorrUtil = new StRefMultCorr("grefmult_P18ih_VpdMB30_AllLumi_MB5sc");
  //grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0); // OLD range
  //grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");  // OLD file
  grefmultCorrUtil->setVzForWeight(16, -16.0, 16.0); // NEW range
  //grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P18ih.txt"); // NEW file
  grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P18ih_set1.txt"); // NEW file

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
  int nMultBins = 800;

  // QA histos
  hEventZVertex = new TH1F("hEventZVertex", "z-vertex distribution", 200, -100., 100.);
  hCentrality   = new TH1F("hCentrality", "No. events vs centrality", nHistCentBins, 0, 100); 
  hMultiplicity = new TH1F("hMultiplicity", "No. events vs multiplicity", 160, 0, 800);

  hMB5refCorr2ReWeight = new TH1F("hMB5refCorr2ReWeight", "MB5 events reweighted to MB30", nMultBins, 0., 800.);
  hMB5refCorr2Raw      = new TH1F("hMB5refCorr2Raw", "MB5 events to MB5 - raw", nMultBins, 0., 800.);
  hMB5CentScaledReWeight = new TH1F("hMB5CentScaledReWeight", "MB5: Centrality Scaled to percent, reweighted", 20, 0., 100.);
  hMB5CentScaledRaw    = new TH1F("hMB5CentScaledRaw", "MB5: Centrality Scaled to percent, raw", 20, 0., 100.);  
  hMB30refCorr2ReWeight = new TH1F("hMB30refCorr2ReWeight", "MB30 events to MB30 - reweighted", nMultBins, 0., 800.);
  hMB30refCorr2Raw     = new TH1F("hMB30refCorr2Raw", "MB30 events to MB30 - raw", nMultBins, 0., 800.);
  hMB30CentScaledReWeight = new TH1F("hMB30CentScaledReWeight", "MB30: Centrality Scaled to percent, reweighted", 20, 0., 100.);
  hMB30CentScaledRaw   = new TH1F("hMB30CentScaledRaw", "MB30: Centrality Scaled to percent, raw", 20, 0., 100.);

  for(int i=0; i<30; i++) { // z-vtx bins (-30, 30)
    hMB5refCorr2[i]    = new TH1F(Form("hMB5refCorr2_%i", i), Form("MB5: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB5refCorr2ReWt[i] = new TH1F(Form("hMB5refCorr2ReWt_%i", i), Form("MB5: refCorr2 re-weighted, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB5refCorr2ReWtInv[i] = new TH1F(Form("hMB5refCorr2ReWtInv_%i", i), Form("MB5: refCorr2 re-weighted inverse, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB5grefMult[i]    = new TH1F(Form("hMB5grefMult%i", i), Form("MB5: grefMult, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB5CentScaled[i]  = new TH1F(Form("hMB5CentScaled%i", i), Form("MB5: Centrality Scaled to percent, z-vtx bin=%i", i), 20, 0., 100.);
    hMB5refCorr2Weight[i]    = new TH1F(Form("hMB5refCorr2Weight_%i", i), Form("MB5 weighted: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB5refCorr2WeightInv[i] = new TH1F(Form("hMB5refCorr2WeightInv_%i", i), Form("MB5 weighted inv: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);

    hMB30refCorr2[i]   = new TH1F(Form("hMB30refCorr2_%i", i), Form("MB30: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB30refCorr2ReWt[i] = new TH1F(Form("hMB30refCorr2ReWt_%i", i), Form("MB30: refCorr2 re-weighted, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB30refCorr2ReWtInv[i] = new TH1F(Form("hMB30refCorr2ReWtInv_%i", i), Form("MB30: refCorr2 re-weighted inverse, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB30grefMult[i]   = new TH1F(Form("hMB30grefMult%i", i), Form("MB30: grefMult, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB30CentScaled[i] = new TH1F(Form("hMB30CentScaled%i", i), Form("MB30: Centrality Scaled to percent, z-vtx bin=%i", i), 20, 0., 100.);

    hHT2refCorr2[i]    = new TH1F(Form("hHT2refCorr2_%i", i), Form("HT2: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hHT2refCorr2ReWt[i] = new TH1F(Form("hHT2refCorr2ReWt_%i", i), Form("HT2: refCorr2 re-weighted, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hHT2refCorr2ReWtInv[i] = new TH1F(Form("hHT2refCorr2ReWtInv_%i", i), Form("HT2: refCorr2 re-weighted inverse, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hHT2grefMult[i]    = new TH1F(Form("hHT2grefMult%i", i), Form("HT2: grefMult, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hHT2CentScaled[i]  = new TH1F(Form("hHT2CentScaled%i", i), Form("HT2: Centrality Scaled to percent, z-vtx bin=%i", i), 20, 0., 100.);

    hMB5onlyrefCorr2Raw[i]  = new TH1F(Form("hMB5onlyrefCorr2Raw_%i", i), Form("MB5 only: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB30onlyrefCorr2Raw[i] = new TH1F(Form("hMB30onlyrefCorr2Raw_%i", i), Form("MB30 only: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
    hMB5MB30refCorr2Raw[i]  = new TH1F(Form("hMB5MB30refCorr2Raw_%i", i), Form("MB5 && MB30: refCorr2, z-vtx bin=%i", i), nMultBins, 0., 800.);
  }
  hMB5onlyrefCorr2RawTotal  = new TH1F("hMB5onlyrefCorr2RawTotal", "MB5 only: refCorr2, all z-vtx", nMultBins, 0., 800.);
  hMB30onlyrefCorr2RawTotal = new TH1F("hMB30onlyrefCorr2RawTotal", "MB30 only: refCorr2, all z-vtx", nMultBins, 0., 800.);
  hMB5MB30refCorr2RawTotal  = new TH1F("hMB5MB30refCorr2RawTotal", "MB5 && MB30: refCorr2, all z-vtx", nMultBins, 0., 800.);

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);

  // set up jet-hadron sparse
  UInt_t bitcode = 0; // bit coded, see GetDimParams() below
  bitcode = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8; // | 1<<9 | 1<<10;
  fhnCentQA = NewTHnSparseF("fhnCentQA", bitcode);

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
  hMB5refCorr2ReWeight->Write();
  hMB5refCorr2Raw->Write();
  hMB5CentScaledReWeight->Write();
  hMB5CentScaledRaw->Write();
  hMB30refCorr2ReWeight->Write();
  hMB30refCorr2Raw->Write();
  hMB30CentScaledReWeight->Write();
  hMB30CentScaledRaw->Write();

  for(int i=0; i<30; i++) {
    hMB5refCorr2[i]->Write();
    hMB5refCorr2ReWt[i]->Write();
    hMB5refCorr2ReWtInv[i]->Write();
    hMB5grefMult[i]->Write();
    hMB5CentScaled[i]->Write();
    hMB5refCorr2Weight[i]->Write();
    hMB5refCorr2WeightInv[i]->Write();

    hMB30refCorr2[i]->Write();
    hMB30refCorr2ReWt[i]->Write();
    hMB30refCorr2ReWtInv[i]->Write();
    hMB30grefMult[i]->Write();
    hMB30CentScaled[i]->Write();

    hHT2refCorr2[i]->Write();
    hHT2refCorr2ReWt[i]->Write();
    hHT2refCorr2ReWtInv[i]->Write();
    hHT2grefMult[i]->Write();
    hHT2CentScaled[i]->Write();

    hMB5onlyrefCorr2Raw[i]->Write();
    hMB30onlyrefCorr2Raw[i]->Write();
    hMB5MB30refCorr2Raw[i]->Write();
  }
  hMB5onlyrefCorr2RawTotal->Write();
  hMB30onlyrefCorr2RawTotal->Write();
  hMB5MB30refCorr2RawTotal->Write();

  // sparse - set up before writing..
  // fhnCentQA->Write();
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
  kCentralityScaled = -99., krefCorr2 = 0.0, kMB5toMB30ReWeight = 1.0, kReWeight = 1.0;

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

  // cut event on max tower Et > 30.0 GeV
  //if(GetMaxTowerEt() > fMaxEventTowerEt) return kStOK;

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
  krefMult  = mPicoEvent->refMult();  // refMult

  // for non-pp analyses
  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(fRunNumber);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(kgrefMult, zVtx, fBBCCoincidenceRate); } 
    else{ grefmultCorr->initEvent(kgrefMult, zVtx, fZDCCoincidenceRate); } // default
//    if(grefmultCorr->isBadRun(fRunNumber)) cout << "Run is bad" << endl; 
//    if(grefmultCorr->isIndexOk()) cout << "Index Ok" << endl;
//    if(grefmultCorr->isZvertexOk()) cout << "Zvertex Ok" << endl;
//    if(grefmultCorr->isRefMultOk()) cout << "RefMult Ok" << endl;

    // get centrality bin: either 0-7 or 0-15 - but decreasing % with increasing value (poor practice, use below)
    kcent16 = grefmultCorr->getCentralityBin16();
    kcent9  = grefmultCorr->getCentralityBin9();

    // re-order binning to be from central -> peripheral
    kref9   = GetCentBin(kcent9, 9);
    kref16  = GetCentBin(kcent16, 16);
    centbin = GetCentBin(kcent16, 16);  // 0-16

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { krefCorr2 = grefmultCorr->getRefMultCorr(kgrefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ krefCorr2 = grefmultCorr->getRefMultCorr(kgrefMult, zVtx, fZDCCoincidenceRate, 2); }

    // get ReWeight for peripheral bins
    kReWeight = grefmultCorr->getWeight();

    // calculate corrected multiplicity: 
    // Double_t getRefMultCorr(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag=2) const ;
    // flag=0:  Luminosity only
    // flag=1:  z-vertex only
    // flag=2:  full correction (default)
    //
    //grefmultCorr->isCentralityOk(cent16)
  } else { // for pp
    centbin = 0, kcent9 = 0, kcent16 = 0;
    krefCorr2 = 0.0, kref9 = 0, kref16 = 0, kReWeight = 0;
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

  //================================================================================================
  //================================================================================================
  // initialize by RunID
  grefmultCorrMB5->init(fRunNumber);

  // initialize for the event
  grefmultCorrMB5->initEvent(kgrefMult, zVtx, fZDCCoincidenceRate); // default
  //const Double_t reweight = grefmultCorrMB5->getWeight();

  // get centrality bin: either 0-7 or 0-15 - but decreasing % with increasing value (poor practice, use below)
  int kcent16MB5 = grefmultCorrMB5->getCentralityBin16();
  int kcent9MB5  = grefmultCorrMB5->getCentralityBin9();

  // re-order binning to be from central -> peripheral
  int kref9MB5   = GetCentBin(kcent9MB5, 9);
  int kref16MB5  = GetCentBin(kcent16MB5, 16);
  int centbinMB5 = GetCentBin(kcent16MB5, 16);  // 0-16

  // calculate corrected multiplicity
  double krefCorr2MB5 = grefmultCorrMB5->getRefMultCorr(kgrefMult, zVtx, fZDCCoincidenceRate, 2);
  double kCentralityScaledMB5 = centbinMB5*5.0;

  // cut on unset centrality, > 80%
  if(kcent16MB5 == -1) {
    cout<<"MB5 definiton crashed as cent > 80%"<<endl;
    return kStOk; // this is for the lowest multiplicity events 80%+ centrality, cut on them here
  }  
  
  //================================================================================================
  // ========================= Trigger Info =============================== //
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);

  // get trigger IDs from PicoEvent class and loop over them
  vector<unsigned int> mytriggers = mPicoEvent->triggerIds();
  for(unsigned int i = 0; i < mytriggers.size(); i++) {
    if(fDebugLevel == kDebugEmcTrigger) cout<<"i = "<<i<<": "<<mytriggers[i] << ", ";
  }

  // check for MB/HT event
  bool fHaveMB5event   = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event  = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30); // MAIN
  bool fHaveHT1Trigger = CheckForHT(fRunFlag, StJetFrameworkPicoBase::kIsHT1);
  bool fHaveHT2Trigger = CheckForHT(fRunFlag, StJetFrameworkPicoBase::kIsHT2);   // MAIN
  bool fHaveHT3Trigger = CheckForHT(fRunFlag, StJetFrameworkPicoBase::kIsHT3);

  //===========================================================================================================================
//  StRefMultCorr *grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P18ih_VpdMB30_AllLumi();
//  StRefMultCorr* grefmultCorrUtil = new StRefMultCorr("grefmult_P18ih_VpdMB30_AllLumi_MB5sc");
  grefmultCorrUtil->init(fRunNumber);  // You need to specify the run number you are going to process

  // scale factor: This factor is needed for Run14 VPDMB5
  // Call initEvent(const UShort_t RefMult, const Double_t z) function event-by-event at the beginning before using any other functions
  grefmultCorrUtil->initEvent(kgrefMult, zVtx, fZDCCoincidenceRate);

  if(fHaveMB5event && !fHaveMB30event) {
    // grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
    // grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");

    //for(Int_t i = 0; i < 6; i++) { 
    for(Int_t i = 0; i < 16; i++) {
      //cout<<i<<" "<<grefmultCorrUtil->get(i, (double)grefmultCorrUtil->getRefMultCorr(kgrefMult, zVtx, fZDCCoincidenceRate, 2) )<<endl;  
    }
  } // MB5-event

  // Corrected refmult (with z-vertex dependent correction and luminositiy correction)
  // NOTE: type should be double or float, not integer
  const Double_t grefmultCorTMP2 = grefmultCorrUtil->getRefMultCorr(kgrefMult, zVtx, fZDCCoincidenceRate, 2);

  // Re-weighting corrections for peripheral bins
  const Double_t tmpWEIGHT = grefmultCorrUtil->getWeight();

  // corrected MB5 to MB30 level - with scaling to weight - FIXME THIS PART BELOW
  //	-  (-99) for unset - meaning it is not a MB5 only event
  //	- IN YOUR ANALYSIS MAKER: use this weight for a MB5 event and use kReWeight for MB30
  //	- if(fHaveMB5event && weight < 0.) return kStOk // this will happen WHEN we have a MB5 AND a MB30 event
  kMB5toMB30ReWeight = (fHaveMB5event && !fHaveMB30event) ? grefmultCorrUtil->getWeight() : 1.0;

  // 2 cm z-vtx bins: starting from -30 to +30cm = 30 bins 
  int ZvtxBin = GetZvtxBin(zVtx);
  if(ZvtxBin < 0) return kStOk; // outside the range (-30, 30)

  // MB5!!!
  if(fHaveMB5event && !fHaveMB30event) {
    hMB5refCorr2ReWeight->Fill(grefmultCorTMP2, (double)kMB5toMB30ReWeight);
    hMB5refCorr2Raw->Fill(grefmultCorTMP2);
    hMB5CentScaledReWeight->Fill(kCentralityScaled, (double)kMB5toMB30ReWeight);
    hMB5CentScaledRaw->Fill(kCentralityScaled);

    hMB5refCorr2[ZvtxBin]->Fill(grefmultCorTMP2);
    hMB5grefMult[ZvtxBin]->Fill(kgrefMult);
    hMB5CentScaled[ZvtxBin]->Fill(kCentralityScaled);

    hMB5refCorr2Weight[ZvtxBin]->Fill(grefmultCorTMP2, kReWeight);
    hMB5refCorr2WeightInv[ZvtxBin]->Fill(grefmultCorTMP2, 1.0/kReWeight);

    hMB5onlyrefCorr2Raw[ZvtxBin]->Fill(grefmultCorTMP2);
    hMB5onlyrefCorr2RawTotal->Fill(grefmultCorTMP2);

    hMB5refCorr2ReWt[ZvtxBin]->Fill(grefmultCorTMP2, (double)kMB5toMB30ReWeight);
    hMB5refCorr2ReWtInv[ZvtxBin]->Fill(grefmultCorTMP2, 1.0/kMB5toMB30ReWeight);
  }

  // MB30!!!
  if(fHaveMB30event && !fHaveMB5event) {
    hMB30refCorr2ReWeight->Fill(krefCorr2, kReWeight);
    hMB30refCorr2Raw->Fill(krefCorr2);
    hMB30CentScaledReWeight->Fill(kCentralityScaled, kReWeight);
    hMB30CentScaledRaw->Fill(kCentralityScaled);

    hMB30refCorr2[ZvtxBin]->Fill(krefCorr2);
    hMB30grefMult[ZvtxBin]->Fill(kgrefMult);
    hMB30CentScaled[ZvtxBin]->Fill(kCentralityScaled);

    hMB30onlyrefCorr2Raw[ZvtxBin]->Fill(krefCorr2);
    hMB30onlyrefCorr2RawTotal->Fill(krefCorr2);

    hMB30refCorr2ReWt[ZvtxBin]->Fill(krefCorr2, kReWeight);
    hMB30refCorr2ReWtInv[ZvtxBin]->Fill(krefCorr2, 1.0/kReWeight);
  }

  // MB5 && MB30
  if(fHaveMB5event && fHaveMB30event) {
    hMB5MB30refCorr2Raw[ZvtxBin]->Fill(krefCorr2);
    hMB5MB30refCorr2RawTotal->Fill(krefCorr2);
  }

  // HT2 trigger!
  if(fHaveHT2Trigger) {
    hHT2refCorr2[ZvtxBin]->Fill(krefCorr2);
    hHT2refCorr2ReWt[ZvtxBin]->Fill(krefCorr2, kReWeight);
    hHT2refCorr2ReWtInv[ZvtxBin]->Fill(krefCorr2, 1.0/kReWeight);
    hHT2grefMult[ZvtxBin]->Fill(kgrefMult);
    hHT2CentScaled[ZvtxBin]->Fill(kCentralityScaled);
  }

  // ========================= fill QA histograms ============================ //
  // double triggerEntries[9] = {(double)kgrefMult, (double)grefMult, refCorr2, fZDCCoincidenceRate, zVtx, 
  //double triggerEntries[9] = {centBinToUse, jetPtselected, pt, deta, dphijh, dEP, zVtx, (double)charge, (double)assocPtBin};
  //fhnCentQA->Fill(triggerEntries, 1.0/trkEfficiency);  // fill Sparse Histo with trigger entries

  // ========================= Print out event info ============================ //
  if(doComments) {
    cout<<"Event Overview: "<<endl;
    cout<<"MB5: "<<fHaveMB5event<<"   MB30: "<<fHaveMB30event<<"   fHaveHT2Trigger: "<<fHaveHT2Trigger<<endl;
    cout<<"  krefMult: "<<krefMult<<"   kgrefMult: "<<kgrefMult<<"   Z-vtx: "<<zVtx<<"   ZDCx: "<<fZDCCoincidenceRate<<"  kref9: "<<kref9<<"   kref16: "<<kref16<<endl;
    cout<<"  refCorr2: "<<krefCorr2<<"  refCorr2MB5: "<<krefCorr2MB5<<"   grefmultCorTMP2: "<<grefmultCorTMP2<<"   REWEIGHT: "<<kReWeight<<"  tmpWEIGHT-MB5: "<<tmpWEIGHT<<endl;
    cout<<"  Cent-scaled: "<<kCentralityScaled<<"  Cent-scaled MB5: "<<kCentralityScaledMB5<<endl<<endl;
  }

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
  hMB5refCorr2ReWeight->Sumw2();
  hMB5refCorr2Raw->Sumw2();
  hMB5CentScaledReWeight->Sumw2();
  hMB5CentScaledRaw->Sumw2();
  hMB30refCorr2ReWeight->Sumw2();
  hMB30refCorr2Raw->Sumw2();
  hMB30CentScaledReWeight->Sumw2();
  hMB30CentScaledRaw->Sumw2();  

  for(int i=0; i<30; i++) {
    hMB5refCorr2[i]->Sumw2();
    hMB5refCorr2ReWt[i]->Sumw2();
    hMB5refCorr2ReWtInv[i]->Sumw2();
    hMB5grefMult[i]->Sumw2();
    hMB5CentScaled[i]->Sumw2();
    hMB5refCorr2Weight[i]->Sumw2();
    hMB5refCorr2WeightInv[i]->Sumw2();

    hMB30refCorr2[i]->Sumw2();
    hMB30refCorr2ReWt[i]->Sumw2();
    hMB30refCorr2ReWtInv[i]->Sumw2();
    hMB30grefMult[i]->Sumw2();
    hMB30CentScaled[i]->Sumw2();

    hHT2refCorr2[i]->Sumw2();
    hHT2refCorr2ReWt[i]->Sumw2();
    hHT2refCorr2ReWtInv[i]->Sumw2();
    hHT2grefMult[i]->Sumw2();
    hHT2CentScaled[i]->Sumw2();

    hMB5onlyrefCorr2Raw[i]->Sumw2();
    hMB30onlyrefCorr2Raw[i]->Sumw2();
    hMB5MB30refCorr2Raw[i]->Sumw2();
  }
  hMB5onlyrefCorr2RawTotal->Sumw2();
  hMB30onlyrefCorr2RawTotal->Sumw2();
  hMB5MB30refCorr2RawTotal->Sumw2();

  fhnCentQA->Sumw2();
}
//
// Get corrected multiplicity for different correction flags
// - luminosity only, z-vtx only, full correction (default)
//_______________________________________________________________________________________________
Double_t StCentMaker::GetCorrectedMultiplicity(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag) {
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
//______________________________________________________________________
THnSparse* StCentMaker::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()
   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){
         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }

      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF

//
//________________________________________________________________________
void StCentMaker::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "centrality";
      nbins = 10;
      xmin = 0.;
      xmax = 100.;     
      break;

   case 1:
        label = "Jet p_{T}";
        nbins = 20;
        xmin = 0.;
        xmax = 100.;
      break;

   case 2:
      label = "Track p_{T}";
      nbins = 80; 
      xmin = 0.;
      xmax = 20.;
      break;

   case 3:
      label = "Relative Eta";
      nbins = 72; // 48
      xmin = -1.8;
      xmax = 1.8;
      break;

   case 4: 
      label = "Relative Phi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;

   case 5:
      label = "Relative angle of Jet and Reaction Plane";
      nbins = 3; // (12) 72
      xmin = 0;
      xmax = 0.5*pi;
      break;

   case 6:
      label = "z-vertex";
      nbins = 20; // 10
      xmin = -40; //-10
      xmax =  40; //+10
      break;

   case 7:
      label = "track charge";
      nbins = 3;
      xmin = -1.5;
      xmax = 1.5;
      break;

   case 8:
      label = "TPC associated pt EPbin";
      nbins = 5;
      xmin = -0.5;
      xmax =  4.5;
      break;

   case 9:
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;

   } // end of switch
} // end of getting dim-params

//
// function to converte "zvertex" to a bin or index from 0-19 
// - assuming -40 to 40 cm setup in 4 cm bins
//__________________________________________________________________________________
Int_t StCentMaker::GetZvtxBin(Double_t zvertex) const
{
  // cut on +/- 30cm
  if(TMath::Abs(zvertex) >= 30.0) return -99;

  // initialize z-vtx bin
  Int_t zBin = -99;

  // get zBin
  zBin = int((zvertex + 30.) / 2.);  // bin width is equal to 2 centimeters 

  return zBin;
}
