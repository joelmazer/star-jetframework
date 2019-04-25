// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// This code is set as an AnalysisMaker task, where it can perform:
// 1) Event Plane reconstruction and corrections
//	- BBC recentering and shifting
//	- ZDC recentering and shifting
// 	- TPC recentering and shifting
//      - Rho (underlying event) subtraction to jets
//      - leading jet tag
//      - general QA
//      
// can get a pointer to:
// 1) collection of jets  	
// 2) event wise rho parameter
// ################################################################

#include "StEventPlaneMaker.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"

// my STAR includes
#include "StJetFrameworkPicoBase.h"
#include "StRhoParameter.h"
#include "StRho.h"
#include "StJetMakerTask.h"
#include "StFemtoTrack.h"
#include "runlistP12id.h" // Run12 pp
#include "runlistP16ij.h"
#include "runlistP17id.h" // SL17i - Run14, now SL18b (March20)
#include "runlistRun14AuAu_P18ih.h" // new Run14 AuAu
#include "StPicoEPCorrectionsIncludes.h"

// old file kept
#include "StPicoConstants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

ClassImp(StEventPlaneMaker)

//_____________________________________________________________________________
StEventPlaneMaker::StEventPlaneMaker(const char* name, StPicoDstMaker *picoMaker, const char* jetMakerName = "", const char* rhoMakerName = "")
  : StJetFrameworkPicoBase(name)
{
  fLeadingJet = 0x0; fSubLeadingJet = 0x0; fExcludeLeadingJetsFromFit = 1.0;
  fTrackWeight = 1; //StJetFrameworkPicoBase::kPtLinear2Const5Weight; // see StJetFrameworkPicoBase::EPtrackWeightType 
  fEventPlaneMaxTrackPtCut = 5.0;
  fTPCEPmethod = 1;
  phi_shift_switch = kFALSE;         // keep off! 
  tpc_recenter_read_switch = kFALSE; // TPC recenter switch
  tpc_shift_read_switch = kFALSE;    // TPC shift switch
  tpc_apply_corr_switch = kFALSE;    // TPC finall corrections
  zdc_recenter_read_switch = kFALSE; // ZDC recenter switch
  zdc_shift_read_switch = kFALSE;    // ZDC shift switch
  zdc_apply_corr_switch = kFALSE;    // ZDC final corrections
  bbc_recenter_read_switch = kFALSE; // BBC recenter switch
  bbc_shift_read_switch = kFALSE;    // BBC shift switch 
  bbc_apply_corr_switch = kFALSE;    // BBC final corrections
  fHistCentBinMin = 0;
  fHistCentBinMax = 9;               // 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80
  fHistZvertBinMin = 0;
  fHistZvertBinMax = 20;             // (-40, 40) 4cm bins
  Q2x_raw = 0.; Q2y_raw = 0.;
  Q2x_p = 0.; Q2x_m = 0.;
  Q2y_p = 0.; Q2y_m = 0.;
  Q2x = 0.; Q2y = 0.;
  TPC_PSI2 = -999;
  TPCA_PSI2 = -999; // subevent A
  TPCB_PSI2 = -999; // subevent B
  BBC_PSI2 = -999; ZDC_PSI2 = -999;
  BBC_PSI1 = -999; ZDC_PSI1 = -999;
  PSI2 = -999;
  RES = -999;
  TPC_raw_comb = 0.; TPC_raw_neg = 0.; TPC_raw_pos = 0.;
  BBC_raw_comb = 0.; BBC_raw_east = 0.; BBC_raw_west = 0.;
  ZDC_raw_comb = 0.; ZDC_raw_east = 0.; ZDC_raw_west = 0.;
  fJets = 0x0;
  fRunNumber = 0;
  fEPcalibFileName = "$STROOT_CALIB/eventplaneFlat.root";
  fEPTPCn = 0.; fEPTPCp = 0.; fEPTPC = 0.; fEPBBC = 0.; fEPZDC = 0.;
  fCalibFile = 0x0; fCalibFile2 = 0x0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  JetMaker = 0;
  RhoMaker = 0;
  grefmultCorr = 0x0;
  mOutNameEP = "";
  doUsePrimTracks = kFALSE;
  fDebugLevel = 0;
  fRunFlag = 0;       // see StMyAnalysisMaker::fRunFlagEnum
  doppAnalysis = kFALSE;
  fCentralityDef = 4; // see StJetFrameworkPicoBase::fCentralityDefEnum
  fRequireCentSelection = kFALSE;
  fCentralitySelectionCut = -99;
  doUseBBCCoincidenceRate = kFALSE; // kFALSE = use ZDC
  fMaxEventTrackPt = 30.0;
  fDoEffCorr = kFALSE;
  fCorrJetPt = kFALSE;
  doEventPlaneRes = kFALSE;
  doTPCptassocBin = kFALSE;
  fTPCptAssocBin = -99;
  doReadCalibFile = kFALSE;
  fEmcTriggerEventType = 0; // see StJetFrameworkPicoBase::fEmcTriggerFlagEnum
  fMBEventType = 2;         // kVPDMB5
  fTriggerToUse = 0;        // kTriggerAny, see StJetFrameworkPicoBase::fTriggerEventTypeEnum
  fJetType = 0;             // full is default
  fMinPtJet = 10.0;         // perhaps this should be lower?
  fJetConstituentCut = 2.0;
  fTrackBias = 0.2;
  fTowerBias = 0.2;
  fJetRad = 0.4;
  fEventZVtxMinCut = -40.0; fEventZVtxMaxCut = 40.0;
  fTrackPtMinCut = 0.2; fTrackPtMaxCut = 30.0;
  fTrackPhiMinCut = 0.0; fTrackPhiMaxCut = 2.0*TMath::Pi();
  fTrackEtaMinCut = -1.0; fTrackEtaMaxCut = 1.0;
  fTrackDCAcut = 3.0;
  fTracknHitsFit = 15; fTracknHitsRatio = 0.52;
  fCentralityScaled = 0.;
  ref16 = -99; ref9 = -99;
  Bfield = 0.0;
//  mVertex = 0x0;
  zVtx = 0.0;
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  fAnalysisMakerName = name;
  fJetMakerName = jetMakerName;
  fRhoMakerName = rhoMakerName;

  InitParameters();
}
//
//
//_________________________________________________________________________________
StEventPlaneMaker::~StEventPlaneMaker()
{ /*  */
  // destructor
  if(hCentrality)    delete hCentrality;
  if(hCentralityEP)   delete hCentralityEP;
  if(hEventPlane)    delete hEventPlane;
  if(fHistEPTPCn)    delete fHistEPTPCn;
  if(fHistEPTPCp)    delete fHistEPTPCp;
  if(fHistEPBBC)     delete fHistEPBBC;
  if(fHistEPZDC)     delete fHistEPZDC;
  for(int i=0; i<9; i++){ // centrality
    if(hTrackPhi[i]) delete hTrackPhi[i];
    if(hTrackPt[i])  delete hTrackPt[i];
  }

  if(fHistEventSelectionQA)          delete fHistEventSelectionQA;
  if(fHistEventSelectionQAafterCuts) delete fHistEventSelectionQAafterCuts;
  if(hTriggerIds)    delete hTriggerIds;
  if(hEmcTriggers)   delete hEmcTriggers;
  if(hTPCepDebug)    delete hTPCepDebug;
  if(hBBCepDebug)    delete hBBCepDebug;
  if(hZDCepDebug)    delete hZDCepDebug;

  if(hZDCDis_W)      delete hZDCDis_W;
  if(hZDCDis_E)      delete hZDCDis_E;
  if(hBBCDis_W)      delete hBBCDis_W;
  if(hBBCDis_E)      delete hBBCDis_E;

  for(int i=0; i<9; i++){ // centrality
    for(int j=0; j<20; j++){ // vz (15)
      if(tpc_recenter_read_switch) {
        if(Q2_p[i][j]) delete Q2_p[i][j];
        if(Q2_m[i][j]) delete Q2_m[i][j];
      }

      if(tpc_shift_read_switch) { 
        if(hTPC_shift_N[i][j]) delete hTPC_shift_N[i][j];
        if(hTPC_shift_P[i][j]) delete hTPC_shift_P[i][j];
      }

      if(bbc_shift_read_switch) { 
        if(hBBC_shift_A[i][j]) delete hBBC_shift_A[i][j];
        if(hBBC_shift_B[i][j]) delete hBBC_shift_B[i][j];
      }

      if(zdc_shift_read_switch) {
        if(hZDC_shift_A[i][j]) delete hZDC_shift_A[i][j];
        if(hZDC_shift_B[i][j]) delete hZDC_shift_B[i][j];
      }

    } // vz
  } // cen

  //if(res_cen) delete res_cen;
  if(hZDC_center_ex) delete hZDC_center_ex;
  if(hZDC_center_ey) delete hZDC_center_ey;
  if(hZDC_center_wx) delete hZDC_center_wx;
  if(hZDC_center_wy) delete hZDC_center_wy;
  if(hBBC_center_ex) delete hBBC_center_ex;
  if(hBBC_center_ey) delete hBBC_center_ey;
  if(hBBC_center_wx) delete hBBC_center_wx;
  if(hBBC_center_wy) delete hBBC_center_wy;

  if(bbc_res)        delete bbc_res;
  if(zdc_psi)        delete zdc_psi;
  if(checkbbc)       delete checkbbc;
  if(psi2_tpc_bbc)   delete psi2_tpc_bbc;
  if(bbc_psi_e)      delete bbc_psi_e;
  if(bbc_psi_w)      delete bbc_psi_w;
  if(bbc_psi_evw)    delete bbc_psi_evw;
  if(bbc_psi1_raw)   delete bbc_psi1_raw;
  if(bbc_psi_raw)    delete bbc_psi_raw;
  if(bbc_psi_rcd)    delete bbc_psi_rcd;
  if(bbc_psi_sft)    delete bbc_psi_sft;
  if(bbc_psi_fnl)    delete bbc_psi_fnl;

  if(zdc_res)        delete zdc_res;
  if(zdc_psi_e)      delete zdc_psi_e;
  if(zdc_psi_w)      delete zdc_psi_w;
  if(zdc_psi_evw)    delete zdc_psi_evw;
  if(zdc_psi_raw)    delete zdc_psi_raw;
  if(zdc_psi_rcd)    delete zdc_psi_rcd;
  if(zdc_psi_sft)    delete zdc_psi_sft;
  if(zdc_psi_fnl)    delete zdc_psi_fnl;

  if(tpc_res)        delete tpc_res;
  if(tpc_psi_N)      delete tpc_psi_N;
  if(tpc_psi_P)      delete tpc_psi_P;
  if(tpc_psi_NvP)    delete tpc_psi_NvP;
  if(tpc_psi_raw)    delete tpc_psi_raw;
  if(tpc_psi_rcd)    delete tpc_psi_rcd;
  if(tpc_psi_sft)    delete tpc_psi_sft;
  if(tpc_psi_fnl)    delete tpc_psi_fnl;

  if(Psi2)           delete Psi2;
  if(Psi2m)          delete Psi2m;
  if(Psi2p)          delete Psi2p;
  if(Delta_Psi2)     delete Delta_Psi2;
  if(Delta_Psi2cyc)  delete Delta_Psi2cyc;
  if(Delta_Psi2old)  delete Delta_Psi2old;
  if(Shift_delta_psi2) delete Shift_delta_psi2;
  if(Psi2_rcd)       delete Psi2_rcd;
  if(Psi2_final)     delete Psi2_final;
  if(Psi2_final_folded) delete Psi2_final_folded;
  if(Psi2_final_raw) delete Psi2_final_raw;

  if(hTPCvsBBCep)    delete hTPCvsBBCep;
  if(hTPCvsZDCep)    delete hTPCvsZDCep;
  if(hBBCvsZDCep)    delete hBBCvsZDCep;

  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      if(fProfV2Resolution[i]) delete fProfV2Resolution[i];
      if(fProfV3Resolution[i]) delete fProfV3Resolution[i];
      if(fProfV4Resolution[i]) delete fProfV4Resolution[i];
      if(fProfV5Resolution[i]) delete fProfV5Resolution[i];
    }
  }

}
//
//
//_______________________________________________________________________________
Int_t StEventPlaneMaker::Init() {
  //StJetFrameworkPicoBase::Init();

  // initialize the histograms
  DeclareHistograms();

  // initialize calibration file for event plane
  //TFile *fCalibFile = new TFile("recenter_calib_file.root", "READ");
  ////fCalibFile = new TFile(fEPcalibFileName.Data(), "READ");

//  fCalibFile = new TFile("StRoot/StMyAnalysisMaker/corrections/recenter_calib_file.root", "READ");
  fCalibFile = new TFile("./recenter_calib_file.root", "READ");
  if(!fCalibFile) cout<<"recenter_calib_file.root does not exist..."<<endl;

//  fCalibFile2 = new TFile("StRoot/StMyAnalysisMaker/corrections/shift_calib_file.root", "READ");
  fCalibFile2 = new TFile("./shift_calib_file.root", "READ");
  if(!fCalibFile2) cout<<"shift_calib_file.root does not exist.."<<endl;

  // Jet TClonesArray
  fJets = new TClonesArray("StJet"); // will have name correspond to the Maker which made it

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

  return kStOK;
}
//
//
//_______________________________________________________________________________________
Int_t StEventPlaneMaker::Finish() { 
  // close event plane calibration files (if open)
  if(fCalibFile->IsOpen())  fCalibFile->Close();
  if(fCalibFile2->IsOpen()) fCalibFile2->Close();

/*
  //  Write event plane histos to file and close it.
  if(mOutNameEP!="") {
    TFile *foutEP = new TFile(mOutNameEP.Data(), "RECREATE");
    foutEP->cd();
    //foutEP->mkdir(GetName());
    //foutEP->cd(GetName());
    WriteEventPlaneHistograms();

    foutEP->cd();
    foutEP->Write();
    foutEP->Close();
  }
*/

  // ============================
  //  Write histos to file and close it.
  if(mOutNameEP!="") {
    TFile *foutEP = new TFile(mOutNameEP.Data(), "UPDATE");
    //fout->ls();
    foutEP->cd();
    foutEP->mkdir(GetName());
    foutEP->cd(GetName());
    WriteEventPlaneHistograms();

    foutEP->cd();
    foutEP->Write();
    foutEP->Close();
  }
  // ============================

  return kStOK;
}
//
//
//______________________________________________________________________________________
void StEventPlaneMaker::DeclareHistograms() {
  // constants
  double pi = 1.0*TMath::Pi();

  // QA histos
  hCentrality = new TH1F("hCentrality", "No. events vs centrality", 20, 0, 100);
  hCentralityEP = new TH1F("hCentralityEP", "No. events vs centrality for EP res", 20, 0, 100);

  hEventPlane = new TH1F("hEventPlane", "Event plane distribution", 72, 0.0, 1.0*pi);
  fHistEPTPCn = new TH2F("fHistEPTPCn", "", 20, 0., 100., 72, -pi, pi);
  fHistEPTPCp = new TH2F("fHistEPTPCp", "", 20, 0., 100., 72, -pi, pi);
  fHistEPBBC = new TH2F("fHistEPBBC", "", 20, 0., 100., 72, -pi, pi);
  fHistEPZDC = new TH2F("fHistEPZDC", "", 20, 0., 100., 72, -pi, pi);

  // track phi distribution for centrality
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i] = new TH1F(Form("hTrackPhi%d", i), Form("track distribution vs #phi, centr%d", i), 144, 0, 2*pi);
    hTrackPt[i] = new TH1F(Form("hTrackPt%d", i), Form("track distribution vs p_{T}, centr%d", i), 120, 0., 30.0);
  }

  // Event Selection QA histo
  fHistEventSelectionQA = new TH1F("fHistEventSelectionQA", "Trigger Selection Counter", 20, 0.5, 20.5);
  fHistEventSelectionQAafterCuts = new TH1F("fHistEventSelectionQAafterCuts", "Trigger Selection Counter after Cuts", 20, 0.5, 20.5);
  hTriggerIds = new TH1F("hTriggerIds", "Trigger Id distribution", 100, 0.5, 100.5);
  hEmcTriggers = new TH1F("hEmcTriggers", "Emcal Trigger counter", 10, 0.5, 10.5);

  hTPCepDebug = new TH1F("hTPCepDebug", "TPC event plane debugging", 10, 0.5, 10.5);
  hTPCepDebug->GetXaxis()->SetBinLabel(1, "Q2x_m <= 0");
  hTPCepDebug->GetXaxis()->SetBinLabel(2, "Q2y_m <= 0");
  hTPCepDebug->GetXaxis()->SetBinLabel(3, "Q2x_p <= 0");
  hTPCepDebug->GetXaxis()->SetBinLabel(4, "Q2y_p <= 0");
  hTPCepDebug->GetXaxis()->SetBinLabel(5, "Q2x <= 0");
  hTPCepDebug->GetXaxis()->SetBinLabel(6, "Q2y <= 0");

  hBBCepDebug = new TH1F("hBBCepDebug", "BBC event plane debugging", 10, 0.5, 10.5);
  hBBCepDebug->GetXaxis()->SetBinLabel(1, "sumcos_E <= 0");
  hBBCepDebug->GetXaxis()->SetBinLabel(2, "sumsin_E <= 0");
  hBBCepDebug->GetXaxis()->SetBinLabel(3, "sumcos_W <= 0");
  hBBCepDebug->GetXaxis()->SetBinLabel(4, "sumsin_W <= 0");
  hBBCepDebug->GetXaxis()->SetBinLabel(6, "sum_E <= 0");
  hBBCepDebug->GetXaxis()->SetBinLabel(7, "sum_W <= 0");

  hZDCepDebug = new TH1F("hZDCepDebug", "ZDC event plane debugging", 10, 0.5, 10.5);
  hZDCepDebug->GetXaxis()->SetBinLabel(1, "mQey <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(2, "mQex <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(3, "mQwy <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(4, "mQwx <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(6, "w_ev <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(7, "w_wv <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(8, "w_eh <= 0");
  hZDCepDebug->GetXaxis()->SetBinLabel(9, "w_wh <= 0");

  float range = 12.;
  float range_b = 2.;
  hZDCDis_W = new TH2F("hZDCDis_W","", 400, -range, range, 400, -range, range);
  hZDCDis_E = new TH2F("hZDCDis_E","", 400, -range, range, 400, -range, range);
  hBBCDis_W = new TH2F("hBBCDis_W","", 400, -range_b, range_b, 400, -range_b, range_b);
  hBBCDis_E = new TH2F("hBBCDis_E","", 400, -range_b, range_b, 400, -range_b, range_b);

  // TPC recentering Q-vector
  if(tpc_recenter_read_switch){
    for(int i=0; i<9; i++){ // centrality
      for(int j=0; j<20; j++){ // vz (15)
        Q2_p[i][j] = new TProfile(Form("Q2_p%d_%d",i,j), Form("Q2 for postive eta,centr%d,vz%d",i,j), 2, 0, 2);
        Q2_m[i][j] = new TProfile(Form("Q2_m%d_%d",i,j), Form("Q2 for minus eta,centr%d,vz%d",i,j), 2, 0, 2);
      }
    }
  }

  // TPC shift of event plane
  if(tpc_shift_read_switch){
    for(int i=0; i<9; i++){ // centrality
      for(int j=0; j<20; j++){ // vz (15)
        hTPC_shift_N[i][j] = new TProfile(Form("hTPC_shift_N%d_%d",i,j), "Nn", 20, 0, 20, -100, 100);
        hTPC_shift_P[i][j] = new TProfile(Form("hTPC_shift_P%d_%d",i,j), "Pn", 20, 0, 20, -100, 100);
      }
    }
  }

  // BBC shift
  if(bbc_shift_read_switch){
    for(int i=0; i<9; i++){ // centrality
      for(int j=0; j<20; j++){ //vz (15)
        hBBC_shift_A[i][j] = new TProfile(Form("hBBC_shift_A%d_%d",i,j), "An", 20, 0, 20, -100, 100);
        hBBC_shift_B[i][j] = new TProfile(Form("hBBC_shift_B%d_%d",i,j), "Bn", 20, 0, 20, -100, 100);
      }
    }
  }

  // ZDC shift
  if(zdc_shift_read_switch){
    for(int i=0; i<9; i++){ // centrality
      for(int j=0; j<20; j++){ //vz (15)
        hZDC_shift_A[i][j] = new TProfile(Form("hZDC_shift_A%d_%d",i,j), "An", 20, 0, 20, -100, 100);
        hZDC_shift_B[i][j] = new TProfile(Form("hZDC_shift_B%d_%d",i,j), "Bn", 20, 0, 20, -100, 100);
      }
    }
  }

  // set binning for run based corrections - run dependent
  Int_t nRunBins = 1; // - just a default
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) nRunBins = 830; //1654;
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) nRunBins = 1359;
  Double_t nRunBinsMax = (Double_t)nRunBins;

  // ZDC centering
  hZDC_center_ex = new TProfile("hZDC_center_ex", "", nRunBins, 0, nRunBinsMax, -100, 100);
  hZDC_center_ey = new TProfile("hZDC_center_ey", "", nRunBins, 0, nRunBinsMax, -100, 100);
  hZDC_center_wx = new TProfile("hZDC_center_wx", "", nRunBins, 0, nRunBinsMax, -100, 100);
  hZDC_center_wy = new TProfile("hZDC_center_wy", "", nRunBins, 0, nRunBinsMax, -100, 100);
  // BBC centering
  hBBC_center_ex = new TProfile("hBBC_center_ex", "", nRunBins, 0, nRunBinsMax, -100, 100);
  hBBC_center_ey = new TProfile("hBBC_center_ey", "", nRunBins, 0, nRunBinsMax, -100, 100);
  hBBC_center_wx = new TProfile("hBBC_center_wx", "", nRunBins, 0, nRunBinsMax, -100, 100);
  hBBC_center_wy = new TProfile("hBBC_center_wy", "", nRunBins, 0, nRunBinsMax, -100, 100);

  //// res_cen=new TProfile("res_cen","res vs. cen",10,0,10,-2,2);
  bbc_res = new TProfile("bbc_res", "", 10, 0, 10, -100, 100);
  zdc_psi = new TH1F("zdc_psi", "zdc psi1", 288, -0.5*pi, 2.5*pi);
  checkbbc = new TH1F("checkbbc", "difference between psi2 and bbc ps2", 288, -0.5*pi, 1.5*pi);
  psi2_tpc_bbc = new TH2F("psi2_tpc_bbc", "tpc psi2 vs. bbc psi2", 288, -0.5*pi, 1.5*pi, 288, -0.5*pi, 1.5*pi);
  bbc_psi_e = new TH1F("bbc_psi_e", "bbc psi2", 288, -0.5*pi, 1.5*pi); // TODO - check order
  bbc_psi_w = new TH1F("bbc_psi_w", "bbc psi2", 288, -0.5*pi, 1.5*pi); // TODO - check order
  bbc_psi_evw = new TH2F("bbc_psi_evw", "bbc psi2 east vs. bbc psi2 west", 288, -0.5*pi, 1.5*pi, 288, -0.5*pi, 1.5*pi);
  bbc_psi1_raw = new TH1F("bbc_psi1_raw", "bbc psi1 raw", 288, -0.5*pi, 2.5*pi);
  bbc_psi_raw = new TH1F("bbc_psi_raw", "bbc psi2 raw", 288, -0.5*pi, 1.5*pi); 
  bbc_psi_rcd = new TH1F("bbc_psi_rcd", "bbc psi2 centered", 288, -0.5*pi, 1.5*pi);  // TODO - check order
  bbc_psi_sft = new TH1F("bbc_psi_sft", "bbc psi2 shifted", 288, -0.5*pi, 1.5*pi);   // TODO - check order
  bbc_psi_fnl = new TH1F("bbc_psi_fnl", "bbc psi2 corrected", 288, -0.5*pi, 1.5*pi); // TODO - check order

  zdc_res = new TProfile("zdc_res", "", 10, 0, 10, -100, 100);
  //zdc_psi = new TH1F("zdc_psi", "zdc psi1", 288, -0.5*pi, 1.5*pi);
  zdc_psi_e = new TH1F("zdc_psi_e", "zdc psi1", 288, -0.5*pi, 1.5*pi);
  zdc_psi_w = new TH1F("zdc_psi_w", "zdc psi1", 288, -0.5*pi, 1.5*pi);
  zdc_psi_evw = new TH2F("zdc_psi_evw", "zdc psi2 east vs. zdc psi2 west", 288, -0.5*pi, 1.5*pi, 288, -0.5*pi, 1.5*pi);
  zdc_psi_raw = new TH1F("zdc_psi_raw", "zdc psi1 raw", 288, -0.5*pi, 1.5*pi);
  zdc_psi_rcd = new TH1F("zdc_psi_rcd", "zdc psi1 centered", 288, -0.5*pi, 1.5*pi);
  zdc_psi_sft = new TH1F("zdc_psi_sft", "zdc psi1 shifted", 288, -0.5*pi, 1.5*pi);
  zdc_psi_fnl = new TH1F("zdc_psi_fnl", "zdc psi1 corrected", 288, -0.5*pi, 1.5*pi);

  tpc_res = new TProfile("tpc_res", "", 10, 0, 10, -100, 100);
  tpc_psi = new TH1F("tpc_psi", "tpc psi2", 288, -0.5*pi, 1.5*pi);
  tpc_psi_N = new TH1F("tpc_psi_N", "tpc psi2 negative", 288, -0.5*pi, 1.5*pi);
  tpc_psi_P = new TH1F("tpc_psi_P", "tpc psi2 positive", 288, -0.5*pi, 1.5*pi);
  tpc_psi_NvP = new TH2F("tpc_psi_NvP", "tpc psi2 negative vs. tpc psi2 positive", 288, -0.5*pi, 1.5*pi, 288, -0.5*pi, 1.5*pi);
  tpc_psi_raw = new TH1F("tpc_psi_raw", "tpc psi2 raw", 288, -0.5*pi, 1.5*pi);
  tpc_psi_rcd = new TH1F("tpc_psi_rcd", "tpc psi2 centered", 288, -0.5*pi, 1.5*pi);
  tpc_psi_sft = new TH1F("tpc_psi_sft", "tpc psi2 shifted", 288, -0.5*pi, 1.5*pi);
  tpc_psi_fnl = new TH1F("tpc_psi_fnl", "tpc psi2 corrected", 288, -0.5*pi, 1.5*pi);

  // these were changed from [-2pi, 2pi] to [0, pi]
  Psi2 = new TH1F("Psi2", "raw #Psi_{2} distribution", 144, 0., 1.*pi);
  Psi2m = new TH1F("Psi2m", "minus eta raw #Psi_{2} distribution", 144, 0., 1.*pi);
  Psi2p = new TH1F("Psi2p", "positive eta raw #Psi_{2} distribution", 144, 0., 1.*pi);
  Delta_Psi2 = new TH1F("Delta_Psi2", "#Delta #Psi_{2} distribution", 144, -1.*pi, 1.0*pi);
  Delta_Psi2cyc = new TH1F("Delta_Psi2cyc", "#Delta #Psi_{2} distribution - cyc", 144, -0.5*pi, 0.5*pi);
  Delta_Psi2old = new TH1F("Delta_Psi2old", "#Delta #Psi_{2} distribution - old", 144, 0., 1.0*pi);
  Shift_delta_psi2 = new TH1F("Shift_delta_psi2", "shift_delta_psi2 distribution", 4000, -8.*pi, 8.*pi);
  Psi2_rcd = new TH1F("Psi2_rcd", "recentered #Psi_{2} distribution", 144, 0., 1.*pi);
  Psi2_final = new TH1F("Psi2_final", "final(shifted and recentered) #Psi_{2} distribution",2000,-4.*pi, 4.*pi);
  Psi2_final_folded = new TH1F("Psi2_final_folded", "final(shifted and recentered) #Psi_{2} distribution", 2000, -4.*pi, 4.*pi);
  Psi2_final_raw = new TH1F("Psi2_final_raw", "final(shifted and recentered) #Psi_{2} distribution", 2000, -4.*pi, 4.*pi);

  // 2-D event plane differences
  hTPCvsBBCep = new TH2F("hTPCvsBBCep", "TPC vs BBC 2nd order event plane", 144, -1*pi, 1*pi, 144, -1*pi, 1*pi);
  hTPCvsZDCep = new TH2F("hTPCvsZDCep", "TPC vs ZDC 2nd order event plane", 144, -1*pi, 1*pi, 144, -1*pi, 1*pi);
  hBBCvsZDCep = new TH2F("hBBCvsZDCep", "BBC vs ZDC 2nd order event plane", 144, -1*pi, 1*pi, 144, -1*pi, 1*pi);

  if(doEventPlaneRes){
    // Reaction Plane resolution as function of centrality - corrected for 2nd order event plane
    for (Int_t i=0; i<9; i++){
      // 2nd order correction to second order event plane
      fProfV2Resolution[i] = new TProfile(Form("fProfV2Resolution_%i", i), Form("fProfV2Resolution_%i", i), 25, 0.5, 25.5);
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(2(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(2(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(2(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(2(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(2(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(2(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(2(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(2(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(2(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(2(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(2(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(2(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(2(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(2(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(2(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(2(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(2(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV2Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(2(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 3rd order correction to second order event plane
      fProfV3Resolution[i] = new TProfile(Form("fProfV3Resolution_%i", i), Form("fProfV3Resolution_%i", i), 25, 0.5, 25.5);
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(3(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(3(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(3(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(3(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(3(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(3(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(3(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(3(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(3(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(3(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(3(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(3(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(3(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(3(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(3(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(3(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(3(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV3Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(3(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 4th order correction to second order event plane
      fProfV4Resolution[i] = new TProfile(Form("fProfV4Resolution_%i", i), Form("fProfV4Resolution_%i", i), 25, 0.5, 25.5);
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(4(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(4(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(4(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(4(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(4(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(4(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(4(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(4(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(4(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(4(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(4(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(4(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(4(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(4(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(4(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(4(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(4(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV4Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(4(#Psi_{ZDC1} - #Psi_{TPCp}))>");

      // 5th order correction to second order event plane
      fProfV5Resolution[i] = new TProfile(Form("fProfV5Resolution_%i", i), Form("fProfV5Resolution_%i", i), 25, 0.5, 25.5);
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(2, "<cos(5(#Psi_{BBC} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(3, "<cos(5(#Psi_{BBC} - #Psi_{ZDC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(4, "<cos(5(#Psi_{TPC} - #Psi_{BBC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(5, "<cos(5(#Psi_{TPC} - #Psi_{ZDC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(6, "<cos(5(#Psi_{ZDC} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(7, "<cos(5(#Psi_{ZDC} - #Psi_{BBC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(8, "<cos(5(#Psi_{BBC} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(9, "<cos(5(#Psi_{BBC} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(10, "<cos(5(#Psi_{ZDC} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(11, "<cos(5(#Psi_{ZDC} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(12, "<cos(5(#Psi_{TPCp} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(17, "<cos(5(#Psi_{BBC1} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(18, "<cos(5(#Psi_{BBC1} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(19, "<cos(5(#Psi_{BBC1} - #Psi_{TPCp}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(20, "<cos(5(#Psi_{BBC1} - #Psi_{ZDC1}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(21, "<cos(5(#Psi_{ZDC1} - #Psi_{TPC}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(22, "<cos(5(#Psi_{ZDC1} - #Psi_{TPCn}))>");
      fProfV5Resolution[i]->GetXaxis()->SetBinLabel(23, "<cos(5(#Psi_{ZDC1} - #Psi_{TPCp}))>");
    }
  }

  // Switch on Sumw2 for all histos
  SetEPSumw2();
}
//
// write event plane histograms
//_____________________________________________________________________________
void StEventPlaneMaker::WriteEventPlaneHistograms() {
  hZDCDis_W->Write();
  hZDCDis_E->Write();
  hBBCDis_W->Write();
  hBBCDis_E->Write();

  // default histos
  hCentrality->Write();
  hCentralityEP->Write();

  hEventPlane->Write();
  fHistEPTPCn->Write();
  fHistEPTPCp->Write();
  fHistEPBBC->Write();
  fHistEPZDC->Write();

  // track phi distribution for centrality
  for(int i=0; i<9; i++){ // centrality
    //hTrackPhi[i]->Write();  // comment out May 17 TODO
    //hTrackPt[i]->Write();   // comment out May 17 TODO
  }

  // QA histos
  fHistEventSelectionQA->Write();
  fHistEventSelectionQAafterCuts->Write();
  hTriggerIds->Write();
  hEmcTriggers->Write();
  hTPCepDebug->Write();
  hBBCepDebug->Write();
  hZDCepDebug->Write();

  for(int i=0; i<9; i++){ // centrality
    for(int j=0; j<20; j++){ // vz (15)
      if(tpc_recenter_read_switch) Q2_p[i][j]->Write();
      if(tpc_recenter_read_switch) Q2_m[i][j]->Write();

      if(tpc_shift_read_switch) hTPC_shift_N[i][j]->Write();
      if(tpc_shift_read_switch) hTPC_shift_P[i][j]->Write();
      
      if(bbc_shift_read_switch) hBBC_shift_A[i][j]->Write();
      if(bbc_shift_read_switch) hBBC_shift_B[i][j]->Write();

      if(zdc_shift_read_switch) hZDC_shift_A[i][j]->Write();
      if(zdc_shift_read_switch) hZDC_shift_B[i][j]->Write();
    }
  }

  if(zdc_recenter_read_switch){
    hZDC_center_ex->Write();
    hZDC_center_ey->Write();
    hZDC_center_wx->Write();
    hZDC_center_wy->Write();
  }
 
  if(bbc_recenter_read_switch){
    hBBC_center_ex->Write();
    hBBC_center_ey->Write();
    hBBC_center_wx->Write();
    hBBC_center_wy->Write();
  }

  bbc_res->Write();
  zdc_psi->Write();
  checkbbc->Write();
  psi2_tpc_bbc->Write();
  bbc_psi_e->Write();
  bbc_psi_w->Write();
  bbc_psi_evw->Write();
  bbc_psi1_raw->Write();
  bbc_psi_raw->Write();
  bbc_psi_rcd->Write();
  bbc_psi_sft->Write();
  bbc_psi_fnl->Write();

  zdc_res->Write();
  zdc_psi_e->Write();
  zdc_psi_w->Write();
  zdc_psi_evw->Write();
  zdc_psi_raw->Write();
  zdc_psi_rcd->Write();
  zdc_psi_sft->Write();
  zdc_psi_fnl->Write();

  tpc_res->Write();
  tpc_psi_N->Write();
  tpc_psi_P->Write();
  tpc_psi_NvP->Write();
  tpc_psi_raw->Write();
  tpc_psi_rcd->Write();
  tpc_psi_sft->Write();
  tpc_psi_fnl->Write();

  Psi2->Write();
  Psi2m->Write();
  Psi2p->Write();
  Delta_Psi2->Write();
  Delta_Psi2cyc->Write();
  Delta_Psi2old->Write();
  Shift_delta_psi2->Write();
  Psi2_rcd->Write();
  Psi2_final->Write();
  Psi2_final_folded->Write();
  Psi2_final_raw->Write();

  hTPCvsBBCep->Write();
  hTPCvsZDCep->Write();
  hBBCvsZDCep->Write();

  // (perhaps temp - save resolution hists to main output file instead of event plane calibration file)
  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      fProfV2Resolution[i]->Write();
      fProfV3Resolution[i]->Write();
      fProfV4Resolution[i]->Write();
      fProfV5Resolution[i]->Write();
    }
  }

}
//
// OLD user code says: //  Called every event after Make(). 
//_____________________________________________________________________________
void StEventPlaneMaker::Clear(Option_t *opt) {
  //fJets->Clear();
}
//
//  This method is called every event.
//_____________________________________________________________________________
Int_t StEventPlaneMaker::Make() {
  // just get out of Dodge if we are trying to run on pp collisions
  if(doppAnalysis) return kStOK; // use kStOK to not create tons of warning printouts

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

  // get event B (magnetic) field
  Bfield = mPicoEvent->bField(); 

  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();
  
  // Z-vertex cut 
  if((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;

  // let me know the Run #, fill, and event ID
  Int_t RunId = mPicoEvent->runId();
  fRunNumber = mPicoEvent->runId();
  Int_t fillId = mPicoEvent->fillId();
  Int_t eventId = mPicoEvent->eventId();
  Double_t fBBCCoincidenceRate = mPicoEvent->BBCx();
  Double_t fZDCCoincidenceRate = mPicoEvent->ZDCx();
  if(fDebugLevel == kDebugGeneralEvt) cout<<"RunID = "<<RunId<<"  fillID = "<<fillId<<"  eventID = "<<eventId<<endl; // what is eventID?

  // ============================ CENTRALITY ============================== //
  // for only 14.5 GeV collisions from 2014 and earlier runs: refMult, for AuAu run14 200 GeV: grefMult 
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_refmult.txt
  // https://github.com/star-bnl/star-phys/blob/master/StRefMultCorr/Centrality_def_grefmult.txt
  int grefMult = mPicoEvent->grefMult();
  //int refMult = mPicoEvent->refMult();
  Int_t centbin, cent9, cent16;
  Double_t refCorr2;

  if(!doppAnalysis) {
    // initialize event-by-event by RunID
    grefmultCorr->init(RunId);
    if(doUseBBCCoincidenceRate) { grefmultCorr->initEvent(grefMult, zVtx, fBBCCoincidenceRate); } // default
    else{ grefmultCorr->initEvent(grefMult, zVtx, fZDCCoincidenceRate); }
    // 10 14 21 29 40 54 71 92 116 145 179 218 263 315 373 441  // RUN 14 AuAu binning

    // get centrality bin: either 0-7 or 0-15
    cent16 = grefmultCorr->getCentralityBin16();
    cent9 = grefmultCorr->getCentralityBin9();

    // re-order binning to be from central -> peripheral
    ref9 = GetCentBin(cent9, 9);
    ref16 = GetCentBin(cent16, 16);
    centbin = GetCentBin(cent16, 16);  // 0-16

    // calculate corrected multiplicity
    if(doUseBBCCoincidenceRate) { refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fBBCCoincidenceRate, 2);
    } else{ refCorr2 = grefmultCorr->getRefMultCorr(grefMult, zVtx, fZDCCoincidenceRate, 2); }

  } else { // for pp
    centbin = 0, cent9 = 0, cent16 = 0, refCorr2 = 0.0, ref9 = 0, ref16 = 0;
  }

  // cut on unset centrality, > 80%
  if(cent16 == -1) return kStWarn; // maybe kStOk; - this is for lowest multiplicity events 80%+ centrality, cut on them

  // centrality / multiplicity histograms
  if(fDebugLevel == kDebugCentrality) { if(centbin > 15) cout<<"centbin = "<<centbin<<"  mult = "<<refCorr2<<"  Centbin*5.0 = "<<centbin*5.0<<"  cent16 = "<<cent16<<endl; }
  fCentralityScaled = centbin*5.0;
  hCentrality->Fill(fCentralityScaled);

  // cut on centrality for analysis before doing anything
  if(fRequireCentSelection) { if(!SelectAnalysisCentralityBin(centbin, fCentralitySelectionCut)) return kStOk; }
  // ============================ end of CENTRALITY ============================== //

  // ========================= Trigger Info =============================== //
  // looking at the EMCal triggers - used for QA and deciding on HT triggers
  // fill Event Trigger QA
  FillEventTriggerQA(fHistEventSelectionQA);
  FillEmcTriggersHist(hEmcTriggers);
 
  // check for MB and HT triggers - Type Flag corresponds to selected type of MB or EMC
  bool fHaveMBevent = CheckForMB(fRunFlag, fMBEventType);
  bool fHaveMB5event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB5);
  bool fHaveMB30event = CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
  bool fHaveEmcTrigger = CheckForHT(fRunFlag, fEmcTriggerEventType);

  // switches for Event Plane analysis
  Bool_t doEPAnalysis = kFALSE;  // set false by default

  // if trigger fired, turn on doEPAnalysis
  if(fHaveEmcTrigger) doEPAnalysis = kTRUE;

  // no need for switch for few checks - TODO will need to update depending on dataset used
  if((fTriggerToUse == StJetFrameworkPicoBase::kTriggerMB) && (!fHaveMB5event) && (!fHaveMB30event)) return kStOK;  // MB triggered event
  if((fTriggerToUse == StJetFrameworkPicoBase::kTriggerHT) && (!fHaveEmcTrigger))                    return kStOK;  // HT triggered event
  // else fTriggerToUse is ANY and we still want to run analysis
  // ======================== end of Triggers ============================= //

  // ================= JetMaker ================ //
  // get JetMaker
  JetMaker = static_cast<StJetMakerTask*>(GetMaker(fJetMakerName));
  const char *fJetMakerNameCh = fJetMakerName;
  if(!JetMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fJetMakerNameCh) << endm;
    return kStWarn;
  }

  // get jet collection associated with JetMaker
  fJets = static_cast<TClonesArray*>(JetMaker->GetJets());
  if(!fJets) {
    LOG_WARN << Form(" No fJets object! Skip! ") << endm;
    return kStWarn;
  }

  // ============== RhoMaker =============== //
  // get RhoMaker from event: old names "StRho_JetsBG", "OutRho", "StMaker#0"
  RhoMaker = static_cast<StRho*>(GetMaker(fRhoMakerName));
  const char *fRhoMakerNameCh = fRhoMakerName;
  if(!RhoMaker) {
    LOG_WARN << Form(" No %s! Skip! ", fRhoMakerNameCh) << endm;
    return kStWarn;
  }

  // set rho object, alt fRho = GetRhoFromEvent(fRhoName);
  fRho = static_cast<StRhoParameter*>(RhoMaker->GetRho());
  if(!fRho) {
    LOG_WARN << Form("Couldn't get fRho object! ") << endm;
    return kStWarn;
  }

  // =========== Event Plane Angle ============= //
  // cache the leading jet within acceptance
  if(fCorrJetPt) {
    fLeadingJet = GetLeadingJet(fJetMakerName, fRho);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName, fRho);
  } else {
    fLeadingJet = GetLeadingJet(fJetMakerName);
    fSubLeadingJet = GetSubLeadingJet(fJetMakerName);
  }

  //cout<<"print1:  TPC_PSI2: "<<TPC_PSI2<<"  TPCA_PSI2: "<<TPCA_PSI2<<"  TPCB_PSI2: "<<TPCB_PSI2<<endl;
  // ==================================================================
  // 1) get z-vertex binning
  // 2) call BBC_EP_Cal to calculate corrections for BBC event plane
  // 3) call ZDC_EP_Cal to calculate corrections for ZDC event plane
  // 4) call EventPlaneCal to calculate corrections for TPC event plane

  // get z-vertex region and cut on it
  int region_vz = GetVzRegion(zVtx);
  if(region_vz > 900) return kStOK;

  // get BBC, ZDC, TPC event planes
  BBC_EP_Cal(ref9, region_vz, 2);
  ZDC_EP_Cal(ref9, region_vz, 2);  // will probably want n=1 for ZDC
  EventPlaneCal(ref9, region_vz, 2, fTPCptAssocBin);
  hEventPlane->Fill(TPC_PSI2);
  //cout<<"print2:  TPC_PSI2: "<<TPC_PSI2<<"  TPCA_PSI2: "<<TPCA_PSI2<<"  TPCB_PSI2: "<<TPCB_PSI2<<endl;

  // compare BBC, ZDC, TPC event planes
  hTPCvsBBCep->Fill(BBC_PSI2, TPC_PSI2);
  hTPCvsZDCep->Fill(ZDC_PSI2, TPC_PSI2);
  hBBCvsZDCep->Fill(ZDC_PSI2, BBC_PSI2);
 
  // test statements (debug)
  if(fDebugLevel == kDebugEventPlaneCalc) {
    cout<<"BBC = "<<BBC_PSI2<<"  ZDC = "<<ZDC_PSI2<<"  TPC_PSI2 = "<<TPC_PSI2<<"  TPCneg = "<<fEPTPCn<<"  TPCpos = "<<fEPTPCp<<"  Multiplicity = "<<refCorr2<<endl;
    cout<<"BBCrawcomb = "<<BBC_raw_comb<<"  BBCrawE = "<<BBC_raw_east<<"  BBCrawW = "<<BBC_raw_west<<endl;
    cout<<"TPCrawcomb = "<<TPC_raw_comb<<"  TPCrawN = "<<TPC_raw_neg<<"  TPCrawP = "<<TPC_raw_pos<<endl;
    cout<<"ZDCrawcomb = "<<ZDC_raw_comb<<"  ZDCrawE = "<<ZDC_raw_east<<"  ZDCrawW = "<<ZDC_raw_west<<endl;
  }

/*
    // my original method for TPC (gets neg + pos)
    //cout<<"Method: kRemoveEtaPhiCone, ptbin = "<<fTPCptAssocBin<<endl;
    //GetEventPlane(kFALSE, 2, kRemoveEtaPhiCone, 2.0, fTPCptAssocBin);  // last param not used (ptcut)
    //cout<<"  bin = "<<fTPCptAssocBin<<"  TPCa = "<<fEPTPCn<<"  TPCb = "<<fEPTPCp<<"  RES = "<<TMath::Cos(2.*(fEPTPCn - fEPTPCp))<<endl;
    ////CalculateEventPlaneResolution(BBC_PSI2, ZDC_PSI2, TPC_PSI2, fEPTPCn, fEPTPCp, BBC_PSI1, ZDC_PSI1);
*/

  // event plane analysis for resolution calculation
  if(doEPAnalysis) {
    hCentralityEP->Fill(fCentralityScaled);

    // calculate / fill event plane resolution histograms
    if(doEventPlaneRes){
      // this version uses my method for TPCn and TPCp with no corrections
      ////CalculateEventPlaneResolution(BBC_PSI2, ZDC_PSI2, TPC_PSI2, fEPTPCn, fEPTPCp, BBC_PSI1, ZDC_PSI1);

      // this version uses the current iteration with random subevent A and B - TODO check this..
      //CalculateEventPlaneResolution(BBC_PSI2, ZDC_PSI2, TPC_PSI2, TPC_raw_neg, TPC_raw_pos, BBC_PSI1, ZDC_PSI1);

      // corrected event planes used for Resolution
      // MAKE sure for meaningful results to be in mode: STEP3
      CalculateEventPlaneResolution(BBC_PSI2, ZDC_PSI2, TPC_PSI2, TPCA_PSI2, TPCB_PSI2, BBC_PSI1, ZDC_PSI1);
    }

    // debug statement
    if(fDebugLevel == kDebugEventPlaneCalc) {
      cout<<"kRemoveEtaPhiCone: "<<"TPC = "<<fEPTPC<<"  TPCn = "<<fEPTPCn<<"  TPCp = "<<fEPTPCp<<" RES = "<<TMath::Cos(2.*(fEPTPCn - fEPTPCp))<<endl;
      cout<<"BBCtpcN Res = "<<1.0*TMath::Cos(2.*(BBC_PSI2 - fEPTPCn))<<"   BBCtpcP Res = "<<1.0*TMath::Cos(2.*(BBC_PSI2 - fEPTPCp))<<endl;
      cout<<endl;
    }

  } // have Emc HT trigger - process event plane resolutions

  return kStOK;
}
//
/*
//________________________________________________________________________
Bool_t StEventPlaneMaker::AcceptJet(StJet *jet) { // for jets
  // applies all jet cuts except pt
  if ((jet->Phi() < fPhimin) || (jet->Phi() > fPhimax)) return kFALSE;
  if ((jet->Eta() < fEtamin) || (jet->Eta() > fEtamax)) return kFALSE;
  if (jet->Area() < fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area() == 0) return kFALSE;
  // exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt() > 20) return kFALSE;
  // jet passed all above cuts
  return kTRUE;
}
*/
//
//_________________________________________________________________________
TH1* StEventPlaneMaker::FillEmcTriggersHist(TH1* h) {
  // number of Emcal Triggers
  for(int i=0; i<8; i++) { fEmcTriggerArr[i] = 0; }
  Int_t nEmcTrigger = mPicoDst->numberOfEmcTriggers();
  //if(fDebugLevel == kDebugEmcTrigger) { cout<<"nEmcTrigger = "<<nEmcTrigger<<endl; }

  // set kAny true to use of 'all' triggers
  fEmcTriggerArr[StJetFrameworkPicoBase::kAny] = 1;  // always TRUE, so can select on all event (when needed/wanted) 

  //static StPicoEmcTrigger* emcTrigger(int i) { return (StPicoEmcTrigger*)picoArrays[picoEmcTrigger]->UncheckedAt(i); }
  // loop over valid EmcalTriggers
  for(int i = 0; i < nEmcTrigger; i++) {
    // get trigger pointer
    StPicoEmcTrigger *emcTrig = static_cast<StPicoEmcTrigger*>(mPicoDst->emcTrigger(i));
    if(!emcTrig) continue;

    // fill for valid triggers
    if(emcTrig->isHT0()) { h->Fill(1); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT0] = 1; }
    if(emcTrig->isHT1()) { h->Fill(2); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT1] = 1; }
    if(emcTrig->isHT2()) { h->Fill(3); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT2] = 1; }
    if(emcTrig->isHT3()) { h->Fill(4); fEmcTriggerArr[StJetFrameworkPicoBase::kIsHT3] = 1; }
    if(emcTrig->isJP0()) { h->Fill(5); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP0] = 1; }
    if(emcTrig->isJP1()) { h->Fill(6); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP1] = 1; }
    if(emcTrig->isJP2()) { h->Fill(7); fEmcTriggerArr[StJetFrameworkPicoBase::kIsJP2] = 1; }
  }
  // kAny trigger - filled once per event
  h->Fill(10); 

  h->GetXaxis()->SetBinLabel(1, "HT0");
  h->GetXaxis()->SetBinLabel(2, "HT1");
  h->GetXaxis()->SetBinLabel(3, "HT2");
  h->GetXaxis()->SetBinLabel(4, "HT3");
  h->GetXaxis()->SetBinLabel(5, "JP0");
  h->GetXaxis()->SetBinLabel(6, "JP1");
  h->GetXaxis()->SetBinLabel(7, "JP2");
  h->GetXaxis()->SetBinLabel(10, "Any");

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");

  return h;
}
//
// Trigger QA histogram, label bins
// check and fill a Event Selection QA histogram for different trigger selections after cuts
//_____________________________________________________________________________
TH1* StEventPlaneMaker::FillEventTriggerQA(TH1* h) {
  // Run12 pp 200 GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    // Run12 (200 GeV pp) triggers:
    int arrHT1[] = {370511, 370546};
    int arrHT2[] = {370521, 370522, 370531, 370980};
    //int arrHT3[] = {380206, 380216}; // NO HT3 triggered events
    int arrMB[] = {370001, 370011, 370983};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { bin = 3; h->Fill(bin); } // HT2
    //if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { bin = 4; h->Fill(bin); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 10; h->Fill(bin); } // VPDMB

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
    h->GetXaxis()->SetBinLabel(2, "BHT1");
    h->GetXaxis()->SetBinLabel(3, "BHT2");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, ""); //"VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, ""); //"Central-5");
    h->GetXaxis()->SetBinLabel(8, ""); //"Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB");
  }

  // Run14 AuAu 200 GeV
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    int arrHT1[] = {450201, 450211, 460201};
    int arrHT2[] = {450202, 450212, 460202, 460212};
    int arrHT3[] = {460203, 450213, 460203};
    int arrMB[] = {450014};
    int arrMB30[] = {450010, 450020};
    int arrCentral5[] = {450010, 450020};
    int arrCentral[] = {460101, 460111};
    int arrMB5[] = {450005, 450008, 450009, 450014, 450015, 450018, 450024, 450025, 450050, 450060};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos
    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { bin = 4; h->Fill(bin); } // HT3 
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); } // MB 
    if(DoComparison(arrCentral5, sizeof(arrCentral5)/sizeof(*arrCentral5))) { bin = 7; h->Fill(bin); }// Central-5
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 8; h->Fill(bin); } // Central & Central-mon
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { bin = 10; h->Fill(bin); }// VPDMB-5 
    if(DoComparison(arrMB30, sizeof(arrMB30)/sizeof(*arrMB30))) { bin = 11; h->Fill(bin); } // VPDMB-30
 
    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
    h->GetXaxis()->SetBinLabel(2, "BHT1*VPDMB-30");
    h->GetXaxis()->SetBinLabel(3, "BHT2*VPDMB-30");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-nobsmd");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, "Central-5");
    h->GetXaxis()->SetBinLabel(8, "Central or Central-mon");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-30");
  }

  // Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    // hard-coded trigger Ids for run16
    //int arrHT0[] = {520606, 520616, 520626, 520636, 520646, 520656};
    int arrHT1[] = {520201, 520211, 520221, 520231, 520241, 520251, 520261, 520605, 520615, 520625, 520635, 520645, 520655, 550201, 560201, 560202, 530201, 540201};
    int arrHT2[] = {530202, 540203};
    int arrHT3[] = {520203, 530213};
    int arrMB[] = {520021};
    int arrMB5[] = {520001, 520002, 520003, 520011, 520012, 520013, 520021, 520022, 520023, 520031, 520033, 520041, 520042, 520043, 520051, 520822, 520832, 520842, 570702};
    int arrMB10[] = {520007, 520017, 520027, 520037, 520201, 520211, 520221, 520231, 520241, 520251, 520261, 520601, 520611, 520621, 520631, 520641};
    int arrCentral[] = {520101, 520111, 520121, 520131, 520141, 520103, 520113, 520123};

    // fill for kAny
    int bin = 0;
    bin = 1; h->Fill(bin);

    // check if event triggers meet certain criteria and fill histos    
    if(DoComparison(arrHT1, sizeof(arrHT1)/sizeof(*arrHT1))) { bin = 2; h->Fill(bin); } // HT1
    if(DoComparison(arrHT2, sizeof(arrHT2)/sizeof(*arrHT2))) { bin = 3; h->Fill(bin); } // HT2
    if(DoComparison(arrHT3, sizeof(arrHT3)/sizeof(*arrHT3))) { bin = 4; h->Fill(bin); } // HT3
    if(DoComparison(arrMB, sizeof(arrMB)/sizeof(*arrMB))) { bin = 5; h->Fill(bin); }  // MB
    if(DoComparison(arrCentral, sizeof(arrCentral)/sizeof(*arrCentral))) { bin = 7; h->Fill(bin); }// Central-5 & Central-novtx
    if(DoComparison(arrMB5, sizeof(arrMB5)/sizeof(*arrMB5))) { bin = 10; h->Fill(bin); } // VPDMB-5 
    if(DoComparison(arrMB10, sizeof(arrMB10)/sizeof(*arrMB10))) { bin = 11; h->Fill(bin); }// VPDMB-10

    // label bins of the analysis trigger selection summary
    h->GetXaxis()->SetBinLabel(1, "un-identified trigger");
    h->GetXaxis()->SetBinLabel(2, "BHT1");
    h->GetXaxis()->SetBinLabel(3, "BHT2");
    h->GetXaxis()->SetBinLabel(4, "BHT3");
    h->GetXaxis()->SetBinLabel(5, "VPDMB-5-p-sst");
    h->GetXaxis()->SetBinLabel(6, "");
    h->GetXaxis()->SetBinLabel(7, "Central");
    h->GetXaxis()->SetBinLabel(8, "");
    h->GetXaxis()->SetBinLabel(10, "VPDMB-5");
    h->GetXaxis()->SetBinLabel(11, "VPDMB-10");
  }

  // set x-axis labels vertically
  h->LabelsOption("v");
  //h->LabelsDeflate("X");
  
  return h;
}
//
// set sumw2 for event plane histograms
// __________________________________________________________________________________
void StEventPlaneMaker::SetEPSumw2() {
  // set sum weights
  //hCentrality->Sumw2();
  //hCentralityEP->Sumw2();

  hEventPlane->Sumw2();
  fHistEPTPCn->Sumw2();
  fHistEPTPCp->Sumw2();
  fHistEPBBC->Sumw2();
  fHistEPZDC->Sumw2();
  for(int i=0; i<9; i++){ // centrality
    hTrackPhi[i]->Sumw2();
    hTrackPt[i]->Sumw2();
  }

  fHistEventSelectionQA->Sumw2();
  fHistEventSelectionQAafterCuts->Sumw2();
  hTriggerIds->Sumw2();
  hEmcTriggers->Sumw2();

  hTPCepDebug->Sumw2();
  hBBCepDebug->Sumw2();
  hZDCepDebug->Sumw2();

  hZDCDis_W->Sumw2();
  hZDCDis_E->Sumw2();
  hBBCDis_W->Sumw2();
  hBBCDis_E->Sumw2();

  // don't turn on Sumw2 for profile plots
/*
  for(int i=0; i<9; i++){ // centrality
    for(int j=0; j<20; j++){ // vz (15)
      if(tpc_recenter_read_switch) Q2_p[i][j]->Sumw2();
      if(tpc_recenter_read_switch) Q2_m[i][j]->Sumw2();
      if(tpc_shift_read_switch) hTPC_shift_N[i][j]->Sumw2();
      if(tpc_shift_read_switch) hTPC_shift_P[i][j]->Sumw2();
      
      if(bbc_shift_read_switch) hBBC_shift_A[i][j]->Sumw2();
      if(bbc_shift_read_switch) hBBC_shift_B[i][j]->Sumw2();
      if(zdc_shift_read_switch) hZDC_shift_A[i][j]->Sumw2();
      if(zdc_shift_read_switch) hZDC_shift_B[i][j]->Sumw2();
    }
  }
  ////res_cen->Sumw2();
  if(zdc_recenter_read_switch){
    hZDC_center_ex->Sumw2();
    hZDC_center_ey->Sumw2();
    hZDC_center_wx->Sumw2();
    hZDC_center_wy->Sumw2();
  }
  if(bbc_recenter_read_switch){
    hBBC_center_ex->Sumw2();
    hBBC_center_ey->Sumw2();
    hBBC_center_wx->Sumw2();
    hBBC_center_wy->Sumw2();
  }
*/

  bbc_res->Sumw2();
  zdc_psi->Sumw2();
  checkbbc->Sumw2();
  psi2_tpc_bbc->Sumw2();
  bbc_psi_e->Sumw2();
  bbc_psi_w->Sumw2();
  bbc_psi_evw->Sumw2();
  bbc_psi1_raw->Sumw2();
  bbc_psi_raw->Sumw2();
  bbc_psi_rcd->Sumw2();
  bbc_psi_sft->Sumw2();
  bbc_psi_fnl->Sumw2();

  zdc_res->Sumw2();
  zdc_psi_e->Sumw2();
  zdc_psi_w->Sumw2();
  zdc_psi_evw->Sumw2();
  zdc_psi_raw->Sumw2();
  zdc_psi_rcd->Sumw2();
  zdc_psi_sft->Sumw2();
  zdc_psi_fnl->Sumw2();

  tpc_res->Sumw2();
  tpc_psi_N->Sumw2();
  tpc_psi_P->Sumw2();
  tpc_psi_NvP->Sumw2();
  tpc_psi_raw->Sumw2();
  tpc_psi_rcd->Sumw2();
  tpc_psi_sft->Sumw2();
  tpc_psi_fnl->Sumw2();

  Psi2->Sumw2();
  Psi2m->Sumw2();
  Psi2p->Sumw2();
  Delta_Psi2->Sumw2();
  Delta_Psi2cyc->Sumw2();
  Delta_Psi2old->Sumw2();
  Shift_delta_psi2->Sumw2();
  Psi2_rcd->Sumw2();
  Psi2_final->Sumw2();
  Psi2_final_folded->Sumw2();
  Psi2_final_raw->Sumw2();

  hTPCvsBBCep->Sumw2();
  hTPCvsZDCep->Sumw2();
  hBBCvsZDCep->Sumw2();

/*
  if(doEventPlaneRes){
    for(Int_t i=0; i<9; i++){
      fProfV2Resolution[i]->Sumw2();
      fProfV3Resolution[i]->Sumw2();
      fProfV4Resolution[i]->Sumw2();
      fProfV5Resolution[i]->Sumw2();
    }
  }
*/

}
//
// Initialize parameters
// __________________________________________________________________________________
void StEventPlaneMaker::InitParameters()
{
}  
//
//
// ________________________________________________________________________________________
void StEventPlaneMaker::GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin)
{ 
  // local variables
  TVector2 mQtpcn, mQtpcp, mQtpc;
  double mQtpcnx = 0., mQtpcny = 0., mQtpcpx = 0., mQtpcpy = 0., mQtpcX = 0., mQtpcY = 0.;
  int order = n;
  double pi = 1.0*TMath::Pi();
  int ntracksNEG = 0, ntracksPOS = 0;

  // leading jet check and removal
  double excludeInEta = -999., excludeInPhi = -999.;
  double excludeInEtaSub = -999., excludeInPhiSub = -999.;
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from EP estimate
    if(fLeadingJet) {
      excludeInEta = fLeadingJet->Eta();
      excludeInPhi = fLeadingJet->Phi();
      //cout<<"leading: pt = "<<fLeadingJet->Pt()<<"  eta = "<<fLeadingJet->Eta()<<"  phi = "<<fLeadingJet->Phi()<<endl;
    }

    // new Jan26
    if(fSubLeadingJet) {
      excludeInEtaSub = fSubLeadingJet->Eta();
      excludeInPhiSub = fSubLeadingJet->Phi();
      //cout<<"subleading: pt = "<<fSubLeadingJet->Pt()<<"  eta = "<<fSubLeadingJet->Eta()<<"  phi = "<<fSubLeadingJet->Phi()<<endl;
    }
  } // leading jets

  // loop over tracks
  TRandom3 *rand = new TRandom3();
  //TRandom *rand = new TRandom();
  int nTOT = 0, nA = 0, nB = 0;
  int nTrack = mPicoDst->numberOfTracks();
  for(int i = 0; i < nTrack; i++) {
    // get track pointer
    StPicoTrack* track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // primary track switch: get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    // more acceptance cuts now - after getting 3-vector
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    ////if(pt > ptcut) continue; // == TEST == //

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    if(doTPCptassocBin) {
      if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
      if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(ptbin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(ptbin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(ptbin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    // remove strip only when we have a leading jet
    // Method1: kRemoveEtaStrip
    //if(fTPCEPmethod == 1){
    if(method == 1){
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) && 
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    //} else if(fTPCEPmethod == 2){
    } else if(method == 2){
      // remove cone (in eta and phi) around leading jet
      // Method2: kRemoveEtaPhiCone
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    //} else if(fTPCEPmethod == 3){
    } else if(method == 3){ 
      // remove tracks above 2 GeV in cone around leading jet
      // Method3: kRemoveLeadingJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
    //} else if(fTPCEPmethod == 4){
    } else if(method == 4){
      // remove strip only when we have a leading + subleading jet
      // Method4: kRemoveEtaStripLeadSubLead
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEtaSub) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    //} else if(fTPCEPmethod == 5){
    } else if(method == 5){
      // remove cone (in eta and phi) around leading + subleading jet
      // Method5: kRemoveEtaPhiConeLeadSubLead
      double deltaR    = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi-excludeInPhiSub)*(phi-excludeInPhiSub));
      if((fLeadingJet)    && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR    < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaRSub < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    //} else if(fTPCEPmethod == 6){
    } else if(method == 6){ 
      // remove tracks above 2 GeV in cone around leading + subleading jet
      // Method6: kRemoveLeadingSubJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi - excludeInPhiSub)*(phi - excludeInPhiSub));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaRSub < fJetRad)) continue;
    } else {
      // DO NOTHING! nothing is removed...
    }

    // configure track weight when performing Q-vector summation
    double trackweight;
    if(fTrackWeight == kNoWeight) {
      trackweight = 1.0;
    } else if(fTrackWeight == kPtLinearWeight) {
      trackweight = pt;
    } else if(fTrackWeight == kPtLinear2Const5Weight) {
      if(pt <= 2.0) trackweight = pt;
      if(pt >  2.0) trackweight = 2.0;
    } else {
      // nothing choosen, so don't use weight
      trackweight = 1.0;
    }

    // test - Jan15 for random subevents
    // generate random distribution from 0 -> 1: and split subevents for [0,0.5] and [0.5, 1]
    double randomNum = rand->Rndm();
    ////double randomNum = gRandom->Rndm();  // > 0.5?

    // split up Q-vectors into 2 random TPC sub-events (A and B)
    if(doTPCptassocBin) {
      if(randomNum >= 0.5) {
        // sum up q-vectors for random event A
        mQtpcpx += trackweight * cos(phi * order);
        mQtpcpy += trackweight * sin(phi * order);
        nA++;
      }
      if(randomNum < 0.5) {
        // sum up q-vectors for random event B
        mQtpcnx += trackweight * cos(phi * order);
        mQtpcny += trackweight * sin(phi * order);
        nB++;
      }
    }  // pt-dependent mode

    // non-pt dependent mode: split up Q-vectors to +/- eta regions
    if(!doTPCptassocBin) {
      if(eta < 0) {
        // sum up q-vectors on negative eta side
        mQtpcnx += trackweight * cos(phi * order);
        mQtpcny += trackweight * sin(phi * order);
        ntracksNEG++;
      }
      if(eta >= 0) {
        // sum up q-vectors on positive eta side
        mQtpcpx += trackweight * cos(phi * order);
        mQtpcpy += trackweight * sin(phi * order);
        ntracksPOS++;
      }
    }  

    // combined TPC event plane q-vectors
    mQtpcX += trackweight * cos(phi * order);
    mQtpcY += trackweight * sin(phi * order);
    nTOT++;

  } // end of track loop

  // test statements
  //cout<<"ntracksNEG = "<<ntracksNEG<<"  ntracksPOS = "<<ntracksPOS<<endl;
  //cout<<"mQtpcpx = "<<mQtpcpx<<"  mQtpcpy = "<<mQtpcpy<<"  mQtpcnx = "<<mQtpcnx<<"  mQtpcny = "<<mQtpcny<<endl;
  //cout<<"nTOT = "<<nTOT<<"  nA = "<<nA<<"  nB = "<<nB<<endl;

  // set q-vector components 
  mQtpcn.Set(mQtpcnx, mQtpcny);
  mQtpcp.Set(mQtpcpx, mQtpcpy);
  mQtpc.Set(mQtpcX, mQtpcY);

  // Calculate the Event Plane
  fEPTPCn = mQtpcn.Phi() / order;
  fEPTPCp = mQtpcp.Phi() / order;
  fEPTPC = mQtpc.Phi() / order;

  // make sure event plane is contrained from [0, pi]
  if((fEPTPCn > pi) || (fEPTPCn < 0)) cout<<"fEPTPCn out of range..... "<<fEPTPCn<<endl;
  if((fEPTPCp > pi) || (fEPTPCp < 0)) cout<<"fEPTPCp out of range..... "<<fEPTPCp<<endl;
  if((fEPTPC  > pi) || (fEPTPC  < 0)) cout<<"fEPTPC out of range..... "<<fEPTPC<<endl;
  if(fEPTPCn > pi) fEPTPCn -= pi;
  if(fEPTPCn <  0) fEPTPCn += pi;
  if(fEPTPCp > pi) fEPTPCp -= pi;
  if(fEPTPCp <  0) fEPTPCp += pi;
  if(fEPTPC  > pi) fEPTPC  -= pi;
  if(fEPTPC  <  0) fEPTPC  += pi;

  // combine x-y vectors for neg and pos Eta ranges (cross-check)
  //double tpcn2 = (0.5*TMath::ATan2(mQtpcny, mQtpcnx));
  //double tpcp2 = (0.5*TMath::ATan2(mQtpcpy, mQtpcpx));
  //double tpc2 = (0.5*TMath::ATan2(mQtpcY, mQtpcX));

  // standard event plane distributions as function of centrality
  //if (fEPTPCResolution!=-1) fHistEPTPCResolution->Fill(fCentrality, fEPTPCResolution);
  fHistEPTPCn->Fill(fCentralityScaled, fEPTPCn);
  fHistEPTPCp->Fill(fCentralityScaled, fEPTPCp);
  fHistEPBBC->Fill(fCentralityScaled, fEPBBC);
  fHistEPZDC->Fill(fCentralityScaled, fEPZDC);

}
// 
// BBC event plane calculation
// ______________________________________________________________________
Int_t StEventPlaneMaker::BBC_EP_Cal(int ref9, int region_vz, int n) { //refmult, the region of vz, and order of EP
  // constant
  double pi = 1.0*TMath::Pi();
  
  // get run ID, transform to run order for filling histogram of corrections
  int RunId = mPicoEvent->runId();
  int RunId_Order = GetRunNo(RunId);// + 1;
  if(RunId_Order < -1) return kStOK;

  // initialize some BBC parameters
  const int N_B = 16;
  double bbc_E[N_B] = {0.};
  double bbc_W[N_B] = {0.};
  double sum_E = 0.;
  double sum_W = 0.;
  for(int i = 0; i < N_B; i++){
    bbc_E[i] = mPicoEvent->bbcAdcEast(i);
    bbc_W[i] = mPicoEvent->bbcAdcWest(i);
    sum_E += bbc_E[i];
    sum_W += bbc_W[i];
  }

  // TODO - double check
//  if(sum_E < 1e-6 || sum_W < 1e-6) { cout<<"BBC: sumE or sumW < 1e-6!!!!"<<endl; return kStOK; } // added test statement Jan5

  // loop over BBC tiles (are these weights?)
  for(int i = 0; i < N_B; i++){
    bbc_E[i] /= sum_E;
    bbc_W[i] /= sum_W;
  }

  // initialize BBC sum and angle parameters
  double sumsin_E = 0., sumcos_E = 0.;
  double sumsin_W = 0., sumcos_W = 0.;

  // loop over BBC tiles
  for(int i = 0; i < N_B; i++){
    double phi_pE = BBC_GetPhi(0, i);
    double phi_pW = BBC_GetPhi(1, i);
    sumsin_E += bbc_E[i]*sin(n*phi_pE); 
    sumcos_E += bbc_E[i]*cos(n*phi_pE);
    sumsin_W += bbc_W[i]*sin(n*phi_pW); 
    sumcos_W += bbc_W[i]*cos(n*phi_pW); 
  }

  // STEP1: for re-centering the BBC event plane angle
  if(bbc_recenter_read_switch){
    hBBC_center_ex->Fill(RunId_Order + 0.5, sumcos_E);
    hBBC_center_ey->Fill(RunId_Order + 0.5, sumsin_E);
    hBBC_center_wx->Fill(RunId_Order + 0.5, sumcos_W);
    hBBC_center_wy->Fill(RunId_Order + 0.5, sumsin_W); // FIXED bug Jan5, 2017 (typo)
  }

  // TEST - debug BBC
  if(fabs(sumcos_E) < 1e-6) { //cout<<"BBC sumcos_E < 1e-6, "<<sumcos_E<<endl; 
    hBBCepDebug->Fill(1.); }
  if(fabs(sumsin_E) < 1e-6) { //cout<<"BBC sumsin_E < 1e-6, "<<sumsin_E<<endl; 
    hBBCepDebug->Fill(2.); }
  if(fabs(sumcos_W) < 1e-6) { //cout<<"BBC sumcos_W < 1e-6, "<<sumcos_W<<endl; 
    hBBCepDebug->Fill(3.); }
  if(fabs(sumsin_W) < 1e-6) { //cout<<"BBC sumsin_W < 1e-6, "<<sumsin_W<<endl; 
    hBBCepDebug->Fill(4.); }
  if(fabs(sum_E) < 1e-6) {    //cout<<"BBC sum_E < 1e-6, "<<sum_E<<endl; 
    hBBCepDebug->Fill(6.); }
  if(fabs(sum_W) < 1e-6) {    //cout<<"BBC sum_W < 1e-6, "<<sum_W<<endl; 
    hBBCepDebug->Fill(7.); }

  // create Q-vectors
  TVector2 bQE, bQW, bQ, bQ_raw, bQ_raw1; 
  bQ_raw1.Set((sumcos_E - sumcos_W), (sumsin_E - sumsin_W)); // this is raw Q vector (1st order event plane)
  bQ_raw.Set((sumcos_E + sumcos_W), (sumsin_E + sumsin_W));  // raw Q vector (2nd order event plane)

  // test statement:
  if(fDebugLevel == kDebugEventPlaneCalc) {
    cout<<"BBC 1st order = "<<bQ_raw1.Phi()<<"  2nd order = "<<bQ_raw.Phi() / n<<endl;
  }
  BBC_PSI1 = bQ_raw1.Phi();

  // raw BBC angles
  TVector2 bbc_comb_raw_vec, bbc_east_raw_vec, bbc_west_raw_vec;
  bbc_comb_raw_vec.Set((sumcos_E + sumcos_W), (sumsin_E + sumsin_W));  // 2nd order
  bbc_east_raw_vec.Set(sumcos_E, sumsin_E);
  bbc_west_raw_vec.Set(sumcos_W, sumsin_W);
  BBC_raw_comb = bbc_comb_raw_vec.Phi(); // TODO - double check
  BBC_raw_east = bbc_east_raw_vec.Phi();
  BBC_raw_west = bbc_west_raw_vec.Phi();
  BBC_raw_comb /= n;
  BBC_raw_east /= n;
  BBC_raw_west /= n;

  // STEP2: correct angles: re-centering and then do shift below
  if(bbc_shift_read_switch) {
    // Method 1: reading values from a function in a *.h file
    // recentering procedure
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
      sumcos_E -= bbc_center_ex_Run14[RunId_Order];
      sumsin_E -= bbc_center_ey_Run14[RunId_Order];
      sumcos_W -= bbc_center_wx_Run14[RunId_Order];
      sumsin_W -= bbc_center_wy_Run14[RunId_Order];
    } else if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
      sumcos_E -= bbc_center_ex[RunId_Order];
      sumsin_E -= bbc_center_ey[RunId_Order];
      sumcos_W -= bbc_center_wx[RunId_Order];
      sumsin_W -= bbc_center_wy[RunId_Order];
    }

/*
    // read from file
    if(fCalibFile && doReadCalibFile){
      // Method 2: reading values from *.root files
      TProfile* htempBBC_center_ex = (TProfile*)fCalibFile->Get("hBBC_center_ex");
      htempBBC_center_ex->SetName("htempBBC_center_ex");
      sumcos_E -= htempBBC_center_ex->GetBinContent(RunId_Order + 1);

      TProfile* htempBBC_center_ey = (TProfile*)fCalibFile->Get("hBBC_center_ey");
      htempBBC_center_ey->SetName("htempBBC_center_ey");
      sumsin_E -= htempBBC_center_ey->GetBinContent(RunId_Order + 1);

      TProfile* htempBBC_center_wx = (TProfile*)fCalibFile->Get("hBBC_center_wx");
      htempBBC_center_wx->SetName("htempBBC_center_wx");
      sumcos_W -= htempBBC_center_wx->GetBinContent(RunId_Order + 1);

      TProfile* htempBBC_center_wy = (TProfile*)fCalibFile->Get("hBBC_center_wy");
      htempBBC_center_wy->SetName("htempBBC_center_wy");
      sumsin_W -= htempBBC_center_wy->GetBinContent(RunId_Order + 1);

      // test statement
      //cout<<"RunID_Order = "<<RunId_Order + 1<<"  hBBC_center_ex max = "<<htempBBC_center_ex->GetMaximum()<<"  bin# = "<<htempBBC_center_ex->GetMaximumBin()<<endl;
      // delete temp histos
      delete htempBBC_center_ex;
      delete htempBBC_center_ey;
      delete htempBBC_center_wx;
      delete htempBBC_center_wy;
    }
*/

  }

  // fill east/west BBC angles
  hBBCDis_W->Fill(sumcos_W, sumsin_W);
  hBBCDis_E->Fill(sumcos_E, sumsin_E);

  // set vector components
  bQE.Set(sumcos_E, sumsin_E);
  bQW.Set(sumcos_W, sumsin_W);
  bQ.Set((sumcos_E + sumcos_W), (sumsin_E + sumsin_W)); // 2nd order

  // set BBC angles
  double bPhi_East = bQE.Phi();
  double bPhi_West = bQW.Phi();
  double bPsi1_raw = bQ_raw1.Phi(); // 1st order
  double bPhi_raw = bQ_raw.Phi();
  double bPhi_rcd = bQ.Phi();
  bPhi_East /= n;
  bPhi_West /= n;
  bPhi_rcd /= n;
  bPhi_raw /= n;

  // STEP2: calculate shift correction and fill histos
  // shift correction to BBC event plane angle
  if(bbc_shift_read_switch){
    for(int s = 1; s < 21; s++) {
      double btimes = double(1./s); //2/(order*2)
      double bBn =   btimes*cos(2*s*bPhi_rcd);
      double bAn = -(btimes*sin(2*s*bPhi_rcd));
      hBBC_shift_A[ref9][region_vz]->Fill(s - 0.5, bAn);
      hBBC_shift_B[ref9][region_vz]->Fill(s - 0.5, bBn);

      // add east/west shifts here... TODO
    }	
  }

  // test statements (debug) east and west
  if(fDebugLevel == kDebugEventPlaneCalc) {
    cout<<"bPhi_East = "<<bPhi_East<<"  bPhi_West = "<<bPhi_West<<endl;
  }

  // STEP3: read shift correction for BBC event plane (read in from file)
  double bbc_delta_psi = 0.;	
  double bbc_shift_Aval = 0., bbc_shift_Bval = 0.; // comment in with code chunk below
  if(bbc_apply_corr_switch) { // need to have ran recentering + shift prior
    for(int nharm = 1; nharm < 21; nharm++){
      // Method 1: load from *.h file function
      // perform 'shift' to BBC event plane angle
      ////bbc_delta_psi += (bbc_shift_A[ref9][region_vz][z-1]*cos(2*z*bPhi_rcd) + bbc_shift_B[ref9][region_vz][z-1]*sin(2*z*bPhi_rcd)); 
      ///bbc_delta_psi += (bbc_shift_A[ref9][region_vz][nharm-1] * cos(2*nharm*bPhi_rcd) + 
      ///                  bbc_shift_B[ref9][region_vz][nharm-1] * sin(2*nharm*bPhi_rcd));

      // recentering procedure
      if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
        bbc_delta_psi += (bbc_shift_A_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*bPhi_rcd) +
                          bbc_shift_B_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*bPhi_rcd));
      } else if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
        bbc_delta_psi += (bbc_shift_A[ref9][region_vz][nharm-1] * cos(2*nharm*bPhi_rcd) +
                          bbc_shift_B[ref9][region_vz][nharm-1] * sin(2*nharm*bPhi_rcd));
      }

/*
      // reading in from file
      if(fCalibFile2 && doReadCalibFile){
        // Method 2: load from *.root calibration file
        TProfile* htempBBC_ShiftA = (TProfile*)fCalibFile2->Get(Form("hBBC_shift_A%i_%i", ref9, region_vz));
        htempBBC_ShiftA->SetName(Form("htempBBC_ShiftA%i_%i", ref9, region_vz));
        bbc_shift_Aval = htempBBC_ShiftA->GetBinContent(nharm);

        TProfile* htempBBC_ShiftB = (TProfile*)fCalibFile2->Get(Form("hBBC_shift_B%i_%i", ref9, region_vz));
        htempBBC_ShiftB->SetName(Form("htempBBC_ShiftB%i_%i", ref9, region_vz));
        bbc_shift_Bval = htempBBC_ShiftB->GetBinContent(nharm);

        // perform 'shift' to BBC event plane angle
        bbc_delta_psi += (bbc_shift_Aval * cos(2*nharm*bPhi_rcd) + bbc_shift_Bval * sin(2*nharm*bPhi_rcd));
        // delete temp histos
        delete htempBBC_ShiftA;
        delete htempBBC_ShiftB;
      }
*/

    }

  }

  int ns = 0;
  ns = int(fabs(bbc_delta_psi) / pi);
  if(bbc_delta_psi > 0) bbc_delta_psi -= ns*pi;
  if(bbc_delta_psi < 0) bbc_delta_psi += ns*pi;

  // shifted BBC event plane angle
  double bPhi_sft = bPhi_rcd + bbc_delta_psi; //(0, pi) + (-pi, pi)= (-pi, 2pi);
  //double b_res = cos(2*(bPhi_East - bPhi_West - pi));
  double b_res = cos(2*(bPhi_East - bPhi_West));

  // Corrected ANGLE: make shifted event plane from {0, pi}
  double bPhi_fnl = bPhi_rcd + bbc_delta_psi;
  if(bPhi_fnl < 0.0)    bPhi_fnl += pi;
  if(bPhi_fnl > 1.0*pi) bPhi_fnl -= pi;

  // fill a bunch of histograms
  ////checkbbc->Fill(PSI2-bPhi_sft);// FIXME
  //psi2_tpc_bbc->Fill(PSI2, bPhi_rcd);
  bbc_res->Fill((ref9) + 0.5, b_res);
  bbc_psi_e->Fill(bPhi_East);
  bbc_psi_w->Fill(bPhi_West);
  bbc_psi_evw->Fill(bPhi_East, bPhi_West);
  bbc_psi1_raw->Fill(bPsi1_raw);// raw 1st order event plane
  bbc_psi_raw->Fill(bPhi_raw);  // raw event plane (TODO order DOUBLE CHECK!)
  bbc_psi_rcd->Fill(bPhi_rcd);  // re-centered event plane
  bbc_psi_sft->Fill(bPhi_sft);  // shifted event plane
  bbc_psi_fnl->Fill(bPhi_fnl);  // final BBC event plane {0, pi}

  // set psi for n=2
  ////PSI2 = bPhi_rcd; // FIXME
  BBC_PSI2 = bPhi_fnl;

/*
  if(TMath::IsNaN(BBC_PSI2)) {
    cout<<endl;
    cout<<"fnl = "<<bPhi_fnl<<"  bPhiE = "<<bPhi_East<<"  bPhiW = "<<bPhi_West<<"  bbcdeltaPsi = "<<bbc_delta_psi<<endl;
    cout<<"sum_E = "<<sum_E<<"  sum_W = "<<sum_W<<"  sumcos_E = "<<sumcos_E<<"  sumcos_W = "<<sumcos_W<<"  sumsin_E = "<<sumsin_E<<"  sumsin_W = "<<sumsin_W<<endl;
    cout<<"BBC_raw_comb = "<<BBC_raw_comb<<"  BBC_raw_east = "<<BBC_raw_east<<"  BBC_raw_west = "<<BBC_raw_west<<endl;
    cout<<"bbc_E = ";
    for(int i = 0; i<16; i++) cout<<bbc_E[i]<<", ";
    cout<<endl<<"bbc_W = ";
    for(int i = 0; i<16; i++) cout<<bbc_W[i]<<", ";
    cout<<endl;
  }
*/

  return kStOk;
}
//
//
// ZDC event plane calculation
// ______________________________________________________________________
Int_t StEventPlaneMaker::ZDC_EP_Cal(int ref9, int region_vz, int n) {
  // constants
  double pi = 1.0*TMath::Pi();

  // get the runID
  int RunId = mPicoEvent->runId();
  int RunId_Order = GetRunNo(RunId);
  if(RunId_Order < -1) return kStOK;

  // initialize some east/west horizontal and vertical values
  double zdc_EH[8] = {0.}; // y
  double zdc_WH[8] = {0.}; // y
  double zdc_EV[7] = {0.}; // x
  double zdc_WV[7] = {0.}; // x

  // loop over horizontal - read in the value of the adc deposition
  for(int i = 0; i < 8; i++){
    zdc_EH[i] = mPicoEvent->ZdcSmdEastHorizontal(i); // y
    zdc_WH[i] = mPicoEvent->ZdcSmdWestHorizontal(i); // y
  }
	
  // loop over vertical - read in the value of the adc deposition
  for(int i = 0; i < 7; i++){
    zdc_EV[i] = mPicoEvent->ZdcSmdEastVertical(i); // x
    zdc_WV[i] = mPicoEvent->ZdcSmdWestVertical(i); // x
  }

  // initialize the east and west horizontal and vertical components
  double eh = 0., wh = 0., ev = 0., wv = 0.;
  double w_eh = 0., w_wh = 0., w_ev = 0., w_wv = 0.;

  // h: horizontal - Y, v: vertical - X     March 20, 2018: this is correct! 
  // https://pdfs.semanticscholar.org/9499/5cee9e50bc55027a8b187681ce8d74c735af.pdf
  // loop over horizontal tiles for ZDCSMD
  for(int i = 0; i < 8; i++){ // y
    eh += zdc_EH[i]*ZDCSMD_GetPosition(RunId_Order, 0, 1, i); // east=0, west=1
    wh += zdc_WH[i]*ZDCSMD_GetPosition(RunId_Order, 1, 1, i); // vertical=0, horizontal=1
    w_eh += zdc_EH[i];
    w_wh += zdc_WH[i];
  }

  // loop over vertical tiles for ZDCSMD
  for(int i = 0; i < 7; i++){ // x
    ev += zdc_EV[i]*ZDCSMD_GetPosition(RunId_Order, 0, 0, i); // east=0, west=1
    wv += zdc_WV[i]*ZDCSMD_GetPosition(RunId_Order, 1, 0, i); // vertical=0, horizontal=1
    w_ev += zdc_EV[i];
    w_wv += zdc_WV[i];
  }

  double mQex, mQey, mQwx, mQwy;
  // written as:
  // z = if(condition) then(?) <do this> else(:) <do this>  
  mQex = (w_ev>0. && w_wv>0. && w_eh>0. && w_wh>0.) ? (ev/w_ev):0.;
  mQey = (w_ev>0. && w_wv>0. && w_eh>0. && w_wh>0.) ? (eh/w_eh):0.;
  mQwx = (w_ev>0. && w_wv>0. && w_eh>0. && w_wh>0.) ? (wv/w_wv):0.;
  mQwy = (w_ev>0. && w_wv>0. && w_eh>0. && w_wh>0.) ? (wh/w_wh):0.;

  // STEP1: fill histograms with recentering corrections for ZDC
  // fill some histograms with Q vectors and after re-centering event plane
  if(w_ev>0. && w_wv>0. && w_eh>0. && w_wh>0.){
    hZDCDis_W->Fill(mQwx, mQwy);
    hZDCDis_E->Fill(mQex, mQey);

    // fill re-centering info for ZDC
    if(zdc_recenter_read_switch){
      hZDC_center_ex->Fill(RunId_Order + 0.5, mQex);
      hZDC_center_ey->Fill(RunId_Order + 0.5, mQey);
      hZDC_center_wx->Fill(RunId_Order + 0.5, mQwx);
      hZDC_center_wy->Fill(RunId_Order + 0.5, mQwy);
    }
  }

  //cout<<"w_ev = "<<w_ev<<"  w_wv = "<<w_wv<<"  w_eh = "<<w_eh<<"  w_wh = "<<w_wh<<endl;
  // TEST - debug ZDC
  if(fabs(mQey) < 1e-6) { //cout<<"ZDC mQey < 1e-6, "<<mQey<<endl; 
    hZDCepDebug->Fill(1.); }
  if(fabs(mQex) < 1e-6) { //cout<<"ZDC mQex < 1e-6, "<<mQex<<endl; 
    hZDCepDebug->Fill(2.); }
  if(fabs(mQwy) < 1e-6) { //cout<<"ZDC mQwy < 1e-6, "<<mQwy<<endl; 
    hZDCepDebug->Fill(3.); }
  if(fabs(mQwx) < 1e-6) { //cout<<"ZDC mQwx < 1e-6, "<<mQwx<<endl; 
    hZDCepDebug->Fill(4.); }
  if(fabs(w_ev) < 1e-6) { //cout<<"ZDC w_ev < 1e-6, "<<w_ev<<endl; 
    hZDCepDebug->Fill(6.); }
  if(fabs(w_wv) < 1e-6) { //cout<<"ZDC w_wv < 1e-6, "<<w_wv<<endl; 
    hZDCepDebug->Fill(7.); }
  if(fabs(w_eh) < 1e-6) { //cout<<"ZDC w_eh < 1e-6, "<<w_eh<<endl; 
    hZDCepDebug->Fill(8.); }
  if(fabs(w_wh) < 1e-6) { //cout<<"ZDC w_wh < 1e-6, "<<w_wh<<endl; 
    hZDCepDebug->Fill(9.); }

  // initialize vectors
  TVector2 mQ1, mQw, mQe, mQtot, zQ_raw;
  mQe.Set(mQex, mQey);
  mQw.Set(mQwx, mQwy);
  mQ1 = mQe - mQw;      // 1st order event plane:  Dec7, why is this written this way? for n=1
  ZDC_PSI1 = mQ1.Phi(); // set global 1st order ZDC event plane
  mQtot.Set((mQex + mQwx), (mQey + mQwy));   // 2nd order - Jan9 found bug!
  zQ_raw.Set((mQex + mQwx), (mQey + mQwy));  // 2nd order
  if(fDebugLevel == kDebugEventPlaneCalc) {
    cout<<"ZDC 1st order: "<<mQ1.Phi()<<"  ZDC 2nd order = "<<mQtot.Phi() / n<<endl;
  }

  // raw ZDC angles
  TVector2 zdc_comb_raw_vec, zdc_east_raw_vec, zdc_west_raw_vec;
  zdc_comb_raw_vec.Set((mQex + mQwx), (mQey + mQwy)); // 2nd order
  zdc_east_raw_vec.Set(mQex, mQey);
  zdc_west_raw_vec.Set(mQwx, mQwy);
  ZDC_raw_comb = zdc_comb_raw_vec.Phi();
  ZDC_raw_east = zdc_east_raw_vec.Phi();
  ZDC_raw_west = zdc_west_raw_vec.Phi();
  ZDC_raw_comb /= n;
  ZDC_raw_east /= n;
  ZDC_raw_west /= n;

  double zPhi_East = -999; // added
  double zPhi_West = -999; // added
  double zPhi_rcd = mQtot.Phi();  // 2nd order - TODO check that a '-999' doesn't get filled and affect average
  double zPhi_raw = zQ_raw.Phi(); // 2nd order

  // get east and west ZDC phi - if weights are valid
  if(w_ev > 1e-6 && w_eh > 1e-6) zPhi_East = mQe.Phi();
  if(w_wv > 1e-6 && w_wh > 1e-6) zPhi_West = mQw.Phi();
  
  // calculate full ZDC event plane - when east and west angles are set
  if(zPhi_East > -900. && zPhi_West > -900.) {
    double psi_full = mQ1.Phi(); // 1st order
    zdc_psi->Fill(psi_full);
  }

  // divide by reaction plane order (harmonic)
  zPhi_East /= n;
  zPhi_West /= n;
  zPhi_raw /= n;
  zPhi_rcd /= n;

  // ============================================
  // ================ NEW ========================
  // STEP2: calculate shift correction and fill histos
  // shift correction to ZDC event plane angle
  if(zdc_shift_read_switch){
    for(int s = 1; s < 21; s++) {
      double btimes = double(1./s); //2/(order*2) - TODO check with LIANG
      double zBn =   btimes*cos(2*s*zPhi_rcd);
      double zAn = -(btimes*sin(2*s*zPhi_rcd));
      hZDC_shift_A[ref9][region_vz]->Fill(s - 0.5, zAn);
      hZDC_shift_B[ref9][region_vz]->Fill(s - 0.5, zBn);
    }
  }

  // STEP3: read shift correction for BBC event plane (read in from file or header)
  double zdc_delta_psi = 0.;
  double zdc_shift_Aval = 0., zdc_shift_Bval = 0.; // comment in with below code chunk
  if(zdc_apply_corr_switch) { // need to have ran recentering + shift prior
    for(int nharm = 1; nharm < 21; nharm++){
      // Method 1: load from *.h file function
      // perform 'shift' to ZDC event plane angle
      ///zdc_delta_psi += (zdc_shift_A[ref9][region_vz][nharm-1] * cos(2*nharm*zPhi_rcd) +
      ///                  zdc_shift_B[ref9][region_vz][nharm-1] * sin(2*nharm*zPhi_rcd));

      if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
        zdc_delta_psi += (zdc_shift_A_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*zPhi_rcd) +
                          zdc_shift_B_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*zPhi_rcd));
      } else if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
        zdc_delta_psi += (zdc_shift_A[ref9][region_vz][nharm-1] * cos(2*nharm*zPhi_rcd) +
                          zdc_shift_B[ref9][region_vz][nharm-1] * sin(2*nharm*zPhi_rcd));
      }

/*
      // reading in from file
      if(fCalibFile2 && doReadCalibFile){
        // Method 2: load from *.root calibration file
        TProfile* htempZDC_ShiftA = (TProfile*)fCalibFile2->Get(Form("hZDC_shift_A%i_%i", ref9, region_vz));
        htempZDC_ShiftA->SetName(Form("htempZDC_ShiftA%i_%i", ref9, region_vz));
        zdc_shift_Aval = htempZDC_ShiftA->GetBinContent(nharm);

        TProfile* htempZDC_ShiftB = (TProfile*)fCalibFile2->Get(Form("hZDC_shift_B%i_%i", ref9, region_vz));
        htempZDC_ShiftB->SetName(Form("htempZDC_ShiftB%i_%i", ref9, region_vz));
        zdc_shift_Bval = htempZDC_ShiftB->GetBinContent(nharm);

        // perform 'shift' to ZDC event plane angle
        zdc_delta_psi += (zdc_shift_Aval * cos(2*nharm*zPhi_rcd) + zdc_shift_Bval * sin(2*nharm*zPhi_rcd));
        // delete temp histos
        delete htempZDC_ShiftA;
        delete htempZDC_ShiftB;
      }      
*/

    }
  }

  int ns = 0;
  ns = int(fabs(zdc_delta_psi) / pi);
  if(zdc_delta_psi > 0) zdc_delta_psi -= ns*pi;
  if(zdc_delta_psi < 0) zdc_delta_psi += ns*pi;

  // shifted ZDC event plane angle
  double zPhi_sft = zPhi_rcd + zdc_delta_psi; //(0, pi) + (-pi, pi)= (-pi, 2pi); // TODO - check
  double z_res = cos(2*(zPhi_East - zPhi_West - pi));

  // make shifted event plane from {0, pi}
  double zPhi_fnl = zPhi_rcd + zdc_delta_psi;
  if(zPhi_fnl < 0.0)    zPhi_fnl += pi;
  if(zPhi_fnl > 1.0*pi) zPhi_fnl -= pi;

  // fill a bunch of histograms  - Added
  zdc_res->Fill((ref9) + 0.5, z_res);
  zdc_psi_e->Fill(zPhi_East);
  zdc_psi_w->Fill(zPhi_West);
  zdc_psi_evw->Fill(zPhi_East, zPhi_West);
  zdc_psi_raw->Fill(zPhi_raw);  // raw event plane (2nd order)
  zdc_psi_rcd->Fill(zPhi_rcd);  // re-centered event plane
  zdc_psi_sft->Fill(zPhi_sft);  // shifted event plane
  zdc_psi_fnl->Fill(zPhi_fnl);  // final ZDC event plane {0, pi}

//  if(w_ev>0. && w_wv>0. && w_eh>0. && w_wh>0.){
    ZDC_PSI2 = zPhi_fnl;
//  } else { 
//    ZDC_PSI2 = -999;
//  }

  // ============================================	
  //  cout<<"zdc phi="<< psi_full<<endl;
  //  PSI2 = psi_full;

  return kStOk;
}
//
//
// get angle of BBC
// ______________________________________________________________________
Double_t StEventPlaneMaker::BBC_GetPhi(int e_w,int iTile){ //east == 0, (west == 1)
  double pi = 1.0*TMath::Pi();
  //TRandom* gRandom(1); // FIXME added for missing variable
  //TRandom rndm;
  double phi_div = pi/6.;
  double bbc_phi = phi_div;

  switch(iTile) {
    case 0: bbc_phi = 3*phi_div;
        break;
    case 1: bbc_phi = phi_div;
        break;
    case 2: bbc_phi = -1*phi_div;
        break;
    case 3: bbc_phi = -3*phi_div;
        break;
    case 4: bbc_phi = -5*phi_div;
        break;
    case 5: bbc_phi = 5*phi_div;
        break;
    case 6: bbc_phi = (gRandom->Rndm() > 0.5) ? 2*phi_div:4*phi_div;
        break;
    case 7: bbc_phi = 3*phi_div;
        break;
    case 8: bbc_phi = phi_div;
        break;
    case 9: bbc_phi = 0.;
        break;
    case 10: bbc_phi = -phi_div;
        break;
    case 11: bbc_phi = (gRandom->Rndm() > 0.5) ? -2*phi_div:-4*phi_div;
        break;
    case 12: bbc_phi = -3*phi_div;
        break;
    case 13: bbc_phi = -5*phi_div;
        break;
    case 14: bbc_phi = pi;
        break;
    case 15: bbc_phi = 5*phi_div;
        break;
    case 16: bbc_phi = 3*phi_div;
        break;
    case 17: bbc_phi = 0.;
        break;
    case 18: bbc_phi = -3*phi_div;
        break;
    case 19: bbc_phi = pi;
        break;
    case 20: bbc_phi = 3*phi_div;
        break;
    case 21: bbc_phi = 0.;
        break;
    case 22: bbc_phi = -3*phi_div;
        break;
    case 23: bbc_phi = pi;
        break;
  }

  if(e_w == 0) {
    if( bbc_phi > -0.001) {bbc_phi = pi - bbc_phi; }
    else{ bbc_phi = -pi - bbc_phi; }
  }
  if (bbc_phi < 0.) { bbc_phi += 2*pi;}

  return bbc_phi;
}
//
// 
// get position of ZDCSMD
// ______________________________________________________________________
Double_t StEventPlaneMaker::ZDCSMD_GetPosition(int id_order, int eastwest, int verthori, int strip) {
  Double_t zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25}; // pre-define the position of each slot.
  Double_t zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};

  // STEP2:
  // correct angles: re-centering ZDC event plane
  //TFile *fZDCcalibFile = new TFile("ZDC_recenter_calib_file.root", "READ");
  double mZDCSMDCenterex = 0.0, mZDCSMDCenterey = 0.0, mZDCSMDCenterwx = 0.0, mZDCSMDCenterwy = 0.0; 
  if(zdc_shift_read_switch || zdc_apply_corr_switch){
    // recentering procedure - read from a function in .h file 
    if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
      mZDCSMDCenterex = zdc_center_ex_Run14[id_order];
      mZDCSMDCenterey = zdc_center_ey_Run14[id_order];
      mZDCSMDCenterwx = zdc_center_wx_Run14[id_order];
      mZDCSMDCenterwy = zdc_center_wy_Run14[id_order];
    } else if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
      mZDCSMDCenterex = zdc_center_ex[id_order];
      mZDCSMDCenterey = zdc_center_ey[id_order];
      mZDCSMDCenterwx = zdc_center_wx[id_order];
      mZDCSMDCenterwy = zdc_center_wy[id_order];
    }

/*
    // reading in from file
    if(fCalibFile && doReadCalibFile){
      // get corrections from histograms of calibration .root file
      TProfile* htempZDC_center_ex = (TProfile*)fCalibFile->Get("hZDC_center_ex");
      htempZDC_center_ex->SetName("htempZDC_center_ex");
      mZDCSMDCenterex = htempZDC_center_ex->GetBinContent(id_order + 1);

      TProfile* htempZDC_center_ey = (TProfile*)fCalibFile->Get("hZDC_center_ey");
      htempZDC_center_ey->SetName("htempZDC_center_ey");
      mZDCSMDCenterey = htempZDC_center_ey->GetBinContent(id_order + 1);

      TProfile* htempZDC_center_wx = (TProfile*)fCalibFile->Get("hZDC_center_wx");
      htempZDC_center_wx->SetName("htempZDC_center_wx");
      mZDCSMDCenterwx = htempZDC_center_wx->GetBinContent(id_order + 1);

      TProfile* htempZDC_center_wy = (TProfile*)fCalibFile->Get("hZDC_center_wy");
      htempZDC_center_wy->SetName("htempZDC_center_wy");
      mZDCSMDCenterwy = htempZDC_center_wy->GetBinContent(id_order + 1);

      delete htempZDC_center_ex;
      delete htempZDC_center_ey;
      delete htempZDC_center_wx;
      delete htempZDC_center_wy;
    }
*/

  }

  // perform re-centering of ZDC event plane angle
  if(zdc_shift_read_switch || zdc_apply_corr_switch) { // TODO double check this - i think this is now FIXED Dec11, 2017
    if(strip < 7 && eastwest == 0 && verthori == 0) return zdcsmd_x[strip] - mZDCSMDCenterex;
    if(strip < 7 && eastwest == 1 && verthori == 0) return -mZDCSMDCenterwx - zdcsmd_x[strip];
    if(eastwest == 0 && verthori == 1) return (zdcsmd_y[strip])/(sqrt(2.)) - mZDCSMDCenterey;
    if(eastwest == 1 && verthori == 1) return (zdcsmd_y[strip])/(sqrt(2.)) - mZDCSMDCenterwy;
  }
// } else {
  if((zdc_recenter_read_switch) && (!zdc_shift_read_switch)){	
    if(strip < 7 && eastwest == 0 && verthori == 0) return zdcsmd_x[strip];
    if(strip < 7 && eastwest == 1 && verthori == 0) return -zdcsmd_x[strip];
    if(eastwest == 0 && verthori == 1) return zdcsmd_y[strip]/sqrt(2.);
    if(eastwest == 1 && verthori == 1) return zdcsmd_y[strip]/sqrt(2.);
  }

  return kStOk;
}
//
//
// this function checks for the bin number of the run from a runlist header 
// in order to apply various corrections and fill run-dependent histograms
// _________________________________________________________________________________
Int_t StEventPlaneMaker::GetRunNo(int runid){ 
  //1287 - Liang

  // Run12 pp (200 GeV)
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200) {
    for(int i = 0; i < 857; i++) {
      if(Run12pp_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run14 AuAu
  // Run14AuAu_IdNo: SL17id
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
    // 1654 for Run14 AuAu, new picoDst production is 830
    for(int i = 0; i < 830; i++) {
      if(Run14AuAu_P18ih_IdNo[i] == runid) {
        return i;
      }
    }
  }

  // Run16 AuAu
  if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
    // 1359 for Run16 AuAu
    for(int i = 0; i < 1359; i++){
      if(Run16AuAu_IdNo[i] == runid) {
        return i;
      }
    }
  }

  cout<<" *********** RunID not matched with list ************!!!! "<<endl;
  return -999;
}
//
//
// this is code from Liang to get the Vz region for event plane corrections
// __________________________________________________________________________________
Int_t StEventPlaneMaker::GetVzRegion(double Vz) // 0-14, 15          0-19, 20
{
  /////////if(fabs(Vz >= 30.)) return 999;  // THIS CUT IS ALREADY DONE
  //int regionvz=int((Vz+30)/6.);  // bin width is equal to 6 centi-meters
  //int regionvz = int((Vz+30.)/4.); // bin width is equal to 4 centi-meters
  //if(regionvz >= 15 || regionvz <= -1) return 999;   // FIXME! don't need this

  // 0-19, 20 bins (-40 to 40)
  int regionvz = int((Vz+40.)/4.); // bin width is equal to 4 centi-meters
  if(regionvz >= 20 || regionvz <= -1) return 999;

  return regionvz;
}
//
//
// Calculate TPC event plane angle with correction
// ___________________________________________________________________________________
Int_t StEventPlaneMaker::EventPlaneCal(int ref9, int region_vz, int n, int ptbin) {
  double res = 0.;
  double pi = 1.0*TMath::Pi();

  // initialize Q-vector components
  Q2x_raw = 0.;
  Q2y_raw = 0.;
  Q2x_p = 0.;
  Q2x_m = 0.;
  Q2y_p = 0.;
  Q2y_m = 0.;
  Q2x = 0.;
  Q2y = 0.;

  // function to calculate Q-vectors
  QvectorCal(ref9, region_vz, n, ptbin);

  // TEST - debug TPC
  if(fabs(Q2x_m) < 1e-6) { cout<<"TPC Q2x_m < 1e-6, "<<Q2x_m<<endl; 
    hTPCepDebug->Fill(1.); }
  if(fabs(Q2y_m) < 1e-6) { cout<<"TPC Q2y_m < 1e-6, "<<Q2y_m<<endl; 
    hTPCepDebug->Fill(2.); }
  if(fabs(Q2x_p) < 1e-6) { cout<<"TPC Q2x_p < 1e-6, "<<Q2x_p<<endl; 
    hTPCepDebug->Fill(3.); }
  if(fabs(Q2y_p) < 1e-6) { cout<<"TPC Q2y_p < 1e-6, "<<Q2y_p<<endl; 
    hTPCepDebug->Fill(4.); }
  if(fabs(Q2x) < 1e-6) {   cout<<"TPC Q2x < 1e-6, "<<Q2x<<endl; 
    hTPCepDebug->Fill(5.); }
  if(fabs(Q2y) < 1e-6) {   cout<<"TPC Q2y < 1e-6, "<<Q2y<<endl; 
    hTPCepDebug->Fill(6.); }

  if(fabs(Q2x_raw == 0.) && fabs(Q2y_raw == 0.)) { cout<<"Q2x_raw or Q2y_raw == 0"<<endl;  return kStOK; }

  // calculate event plane angles - CHECK RANGES TODO
  double psi2 = atan2(Q2y_raw, Q2x_raw); // (-pi, pi]
  double psi2p = atan2(Q2y_p, Q2x_p);    // (-pi, pi]
  double psi2m = atan2(Q2y_m, Q2x_m);    // (-pi, pi]
  double tPhi_rcd = atan2(Q2y, Q2x);     // (-pi, pi)

  // forces range {0, 2pi}
  if(psi2 < 0.)     psi2 += 2*pi;      // (0, 2*pi]
  if(psi2m < 0.)    psi2m += 2*pi;     // (0, 2*pi]
  if(psi2p < 0.)    psi2p += 2*pi;     // (0, 2*pi]
  if(tPhi_rcd < 0.) tPhi_rcd += 2*pi;  // (0, 2*pi]

  // divide by 2 (order) b/c  2*theta = tan-1(y/x)
  psi2 /= n;     // (0, pi)
  psi2m /= n;    // (0, pi)
  psi2p /= n;    // (0, pi)
  tPhi_rcd /= n; // (0, pi)

  // temp
  TVector2 tpc_comb_raw_vec, tpc_neg_raw_vec, tpc_pos_raw_vec;
  tpc_comb_raw_vec.Set(Q2x_raw, Q2y_raw);
  tpc_neg_raw_vec.Set(Q2x_m, Q2y_m); // TODO - might not be raw for STEP3
  tpc_pos_raw_vec.Set(Q2x_p, Q2y_p); // TODO - might not be raw for STEP3
  TPC_raw_comb = tpc_comb_raw_vec.Phi();
  TPC_raw_neg = tpc_neg_raw_vec.Phi(); // read above
  TPC_raw_pos = tpc_pos_raw_vec.Phi(); // read above
  TPC_raw_comb /= n;
  TPC_raw_neg /= n;
  TPC_raw_pos /= n;

  // fill TPC event plane histos
  Psi2->Fill(psi2);              // raw psi2
  Psi2_rcd->Fill(tPhi_rcd);      // recentered psi2
  Delta_Psi2->Fill(psi2m-psi2p); // raw delta psi2 - full range [-pi, pi]
  Delta_Psi2old->Fill(psi2m-psi2p); // raw delta psi2 - old
 
  // sanity check when doing unfolding - create angular difference [-pi/2, pi/2] - this will restrict the max range when used for smearing
  double deltaPsi2 = psi2m-psi2p; // angles from [-pi, pi]
  if(deltaPsi2 >  0.5*pi) deltaPsi2 -= pi;
  if(deltaPsi2 < -0.5*pi) deltaPsi2 += pi;
  Delta_Psi2cyc->Fill(deltaPsi2);

  // used for resolution calculation
  double t_res = cos(2*(psi2m - psi2p)); // added
  res = 2.*(cos(2*(psi2m - psi2p)));
  RES = res;

  // TODO - double check this! 
  if(fabs(Q2x_m) > 1e-6 || fabs(Q2y_m) > 1e-6) Psi2m->Fill(psi2m);
  if(fabs(Q2x_p) > 1e-6 || fabs(Q2y_p) > 1e-6) Psi2p->Fill(psi2p);

  //================================shift reading
  // STEP2: perform shift
  if(tpc_shift_read_switch){
    // loop over (s) harmonics
    for(int s = 1; s < 21; s++) {
      double times = double(1./s); //2/(order*2)
      double Bn =   times*cos(2*s*tPhi_rcd);
      double An = -(times*sin(2*s*tPhi_rcd));
      hTPC_shift_N[ref9][region_vz]->Fill(s - 0.5, An); // shift_A
      hTPC_shift_P[ref9][region_vz]->Fill(s - 0.5, Bn); // shift_B
    }  
  }

  //=================================shift correction
  // STEP3: read shift correction fo TPC event plane (read in from file)
  double tpc_delta_psi = 0.;
  double tpc_shift_Aval = 0., tpc_shift_Bval = 0.; // comment in with code chunk below
  if(tpc_apply_corr_switch) { // FIXME: file needs to exist and need to have ran recentering + shift prior
    // loop over harmonics
    for(int nharm = 1; nharm < 21; nharm++){
      // Method 1: load from *.h file function
      // perform 'shift' to TPC event plane angle
      // KEEP in mind, the naming convent here means nothing for functions: tpc_shift_N and tpc_shift_P,
      // they are corresponing to Bn and An components above!
      // so hTPC_shft_N and hTPC_shift_P also have misleading names
      //shift_delta_psi += (shift_A[ref9][region_vz][z-1]*cos(2*z*psi2_rcd) + shift_B[ref9][region_vz][z-1]*sin(2*z*psi2_rcd));
      //tpc_delta_psi += (tpc_shift_N[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
      //                  tpc_shift_P[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));

      if(doTPCptassocBin) {
        //tpc_delta_psi += GetTPCShiftingValue(tPhi_rcd, nharm, ref9, region_vz);
        tpc_delta_psi += GetTPCShiftingValueNEW(tPhi_rcd, nharm, ref9, region_vz);
      } else {
        // standard default method (all pt bins combined for reaction plane calculation)
        tpc_delta_psi += (tpc_shift_N[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
      }

/*
      // reading in from file
      // started correcting names..
      if(fCalibFile2 && doReadCalibFile){
        // Method 2: load from *.root calibration file
        TProfile* htempTPC_ShiftA = (TProfile*)fCalibFile2->Get(Form("hTPC_shift_N%i_%i", ref9, region_vz));
        htempTPC_ShiftA->SetName(Form("htempTPC_ShiftA%i_%i", ref9, region_vz));
        tpc_shift_Aval = htempTPC_ShiftA->GetBinContent(nharm);

        TProfile* htempTPC_ShiftB = (TProfile*)fCalibFile2->Get(Form("hTPC_shift_P%i_%i", ref9, region_vz));
        htempTPC_ShiftB->SetName(Form("htempTPC_ShiftB%i_%i", ref9, region_vz));
        tpc_shift_Bval = htempTPC_ShiftB->GetBinContent(nharm);

        // perform 'shift' to TPC event plane angle
        tpc_delta_psi += (tpc_shift_Aval * cos(2*nharm*tPhi_rcd) + tpc_shift_Bval * sin(2*nharm*tPhi_rcd));
        // delete temp histos
        delete htempTPC_ShiftA;
        delete htempTPC_ShiftB;
      }
*/

    } // nharm loop
  } // correction switch

  int ns = 0;
  //=====tpc_delta_psi (-pi, pi)
  ns = int(fabs(tpc_delta_psi) / pi);
  if(tpc_delta_psi > 0) tpc_delta_psi -= ns*pi;
  if(tpc_delta_psi < 0) tpc_delta_psi += ns*pi;

  // shifted TPC event plane angle
  double tPhi_sft = tPhi_rcd + tpc_delta_psi; //(0, pi) + (-pi, pi)= (-pi, 2pi);  TOOOOOOOOOOOOOOOOOOOOOOOO
  Psi2_final_raw->Fill(tPhi_sft);
  Shift_delta_psi2->Fill(tpc_delta_psi);

  // make shifted event plane from {0, pi}    -- ADDED
  double tPhi_fnl = tPhi_rcd + tpc_delta_psi;
  if(tPhi_fnl < 0.0)   tPhi_fnl += pi;
  if(tPhi_fnl > 1.0*pi) tPhi_fnl -= pi;

  double shifted_psi2_raw;
  if(tPhi_sft < 0.)  shifted_psi2_raw = tPhi_sft + pi;
  if(tPhi_sft >= pi) shifted_psi2_raw = tPhi_sft - pi;
  if(tPhi_sft > 0. && tPhi_sft < pi) shifted_psi2_raw = tPhi_sft;
  Psi2_final_folded->Fill(shifted_psi2_raw);
  Psi2_final->Fill(tPhi_sft);

  //================
  // fill a bunch of histograms  - Added
  tpc_res->Fill((ref9) + 0.5, t_res);
  tpc_psi_N->Fill(psi2m);          // TODO - update name (tPhi_N)
  tpc_psi_P->Fill(psi2p);          // TODO - update name (tPhi_P)
  tpc_psi_NvP->Fill(psi2m, psi2p); // TODO - update name (tPhi_N, tPhi_P)
  tpc_psi_raw->Fill(psi2);         // raw event plane TODO - update name (tPhi_raw)
  tpc_psi_rcd->Fill(tPhi_rcd);     // re-centered event plane
  tpc_psi_sft->Fill(tPhi_sft);     // shifted event plane
  tpc_psi_fnl->Fill(tPhi_fnl);     // final TPC event plane {0, pi}
  //================   

  //PSI2= tPhi_rcd;
  //PSI2= psi2;
  TPC_PSI2 = tPhi_fnl;
  TPCA_PSI2 = psi2p;
  TPCB_PSI2 = psi2m;

  return kStOk;
}
//
// this is a function for Qvector calculation for TPC event plane
// ______________________________________________________________________________________________
void StEventPlaneMaker::QvectorCal(int ref9, int region_vz, int n, int ptbin) {
  TVector2 mQtpcn, mQtpcp;
  //double mQtpcnx = 0., mQtpcny = 0., mQtpcpx = 0., mQtpcpy = 0.;
  int order = n; //2;
  double pi = 1.0*TMath::Pi();
  int ntracksNEG = 0, ntracksPOS = 0;
  int nTOT = 0, nA = 0, nB = 0; // counter for sub-event A & B

  // get random function to select sub-events
  TRandom3 *rand = new TRandom3();
  //TRandom *rand = new TRandom();

  // leading jet check and removal
  double excludeInEta = -999, excludeInPhi = -999;
  double excludeInEtaSub = -999., excludeInPhiSub = -999.;
  if(fExcludeLeadingJetsFromFit > 0 ) {    // remove the leading jet from EP estimate
    if(fLeadingJet) {
      excludeInEta = fLeadingJet->Eta();
      excludeInPhi = fLeadingJet->Phi();
    }

    if(fSubLeadingJet) {
      excludeInEtaSub = fSubLeadingJet->Eta();
      excludeInPhiSub = fSubLeadingJet->Phi();
    }
  } // get location of leading jets for removal

  // loop over tracks
  int Qtrack = mPicoDst->numberOfTracks();
  for(int i = 0; i < Qtrack; i++){
    // get track pointer
    StPicoTrack *track = static_cast<StPicoTrack*>(mPicoDst->track(i));
    if(!track) { continue; }

    // apply standard track cuts - (can apply more restrictive cuts below)
    if(!(AcceptTrack(track, Bfield, mVertex))) { continue; }

    // primary track switch: get momentum vector of track - global or primary track
    TVector3 mTrkMom;
    if(doUsePrimTracks) {
      // get primary track vector
      mTrkMom = track->pMom();
    } else {
      // get global track vector
      mTrkMom = track->gMom(mVertex, Bfield);
    }

    // track variables
    double pt = mTrkMom.Perp();
    double phi = mTrkMom.Phi();
    double eta = mTrkMom.PseudoRapidity();

    // should set a soft pt range (0.2 - 5.0?)
    if(pt > fEventPlaneMaxTrackPtCut) continue;   // 5.0 GeV
    if(phi < 0.0)    phi += 2.0*pi;
    if(phi > 2.0*pi) phi -= 2.0*pi;

    // 0.20-0.5, 0.5-1.0, 1.0-1.5, 1.5-2.0    - also added 2.0-3.0, 3.0-4.0, 4.0-5.0
    // when doing event plane calculation via pt assoc bin
    if(doTPCptassocBin) {
      if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations // FIXME - why is this 0.25??
      if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
      if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
      if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
      if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
      if(ptbin == 5) { if((pt > 2.00) && (pt <= 3.0)) continue; }  // 2.00 - 3.0 GeV assoc bin used for correlations
      if(ptbin == 6) { if((pt > 3.00) && (pt <= 4.0)) continue; }  // 3.00 - 4.0 GeV assoc bin used for correlations
      if(ptbin == 7) { if((pt > 4.00) && (pt <= 5.0)) continue; }  // 4.00 - 5.0 GeV assoc bin used for correlations
    }

    // FIXME FIXME bugs found here May25, 2018 - make sure they are fixed
    // remove strip only when we have a leading jet
    // Method1: kRemoveEtaStrip
    if(fTPCEPmethod == 1){
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    } else if(fTPCEPmethod == 2){
      // remove cone (in eta and phi) around leading jet
      // Method2: kRemoveEtaPhiCone - FIXME bug found May25th, 2018
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    } else if(fTPCEPmethod == 3){
      // remove tracks above 2 GeV in cone around leading jet
      // Method3: kRemoveLeadingJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
    } else if(fTPCEPmethod == 4){
      // remove strip only when we have a leading + subleading jet
      // Method4: kRemoveEtaStripLeadSubLead
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEta) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((TMath::Abs(eta - excludeInEtaSub) < fJetRad*fExcludeLeadingJetsFromFit ) ||
        ((TMath::Abs(eta) - fJetRad - 1.0 ) > 0) )) continue;
    } else if(fTPCEPmethod == 5){
      // remove cone (in eta and phi) around leading + subleading jet
      // Method5: kRemoveEtaPhiConeLeadSubLead
      double deltaR    = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi-excludeInPhiSub)*(phi-excludeInPhiSub));
      if((fLeadingJet)    && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaR    < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        ((deltaRSub < fJetRad) || (TMath::Abs(eta) - fJetRad - 1.0 > 0 ) )) continue;
    } else if(fTPCEPmethod == 6){
      // remove tracks above 2 GeV in cone around leading + subleading jet
      // Method6: kRemoveLeadingSubJetConstituents
      double deltaR = 1.0*TMath::Sqrt((eta - excludeInEta)*(eta - excludeInEta) + (phi - excludeInPhi)*(phi - excludeInPhi));
      double deltaRSub = 1.0*TMath::Sqrt((eta - excludeInEtaSub)*(eta - excludeInEtaSub) + (phi - excludeInPhiSub)*(phi - excludeInPhiSub));
      if((fLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaR < fJetRad)) continue;
      if((fSubLeadingJet) && (fExcludeLeadingJetsFromFit > 0) &&
        (pt > fJetConstituentCut) && (deltaRSub < fJetRad)) continue;
    } else {
      // DO NOTHING! nothing is removed...
    }

    // Liang cuts
    //if(pt >= 2.) continue;
    //if(fabs(eta)<=0.75) continue; //method 1
    //if(fabs(eta)>=0.25 || fabs(eta)<0.05) continue;   //method 2

    // configure track weight when performing Q-vector summation
    double trackweight;
    if(fTrackWeight == kNoWeight) {
      trackweight = 1.0;
    } else if(fTrackWeight == kPtLinearWeight) {
      trackweight = pt;
    } else if(fTrackWeight == kPtLinear2Const5Weight) {
      if(pt <= 2.0) trackweight = pt;
      if(pt >  2.0) trackweight = 2.0;
    } else {
      // nothing choosen, so don't use weight
      trackweight = 1.0;
    }

    // components (x and y)    (no segregation of minus and positive regions HERE - double check!)
    double x = trackweight * cos(order*phi);
    double y = trackweight * sin(order*phi);
    Q2x_raw += x;
    Q2y_raw += y;

    // test - Jan15 for random subevents
    // generate random distribution from 0 -> 1: and split subevents for [0,0.5] and [0.5, 1]
    double randomNum = rand->Rndm();
    //double randomNum = gRandom->Rndm();  // > 0.5?
    if(randomNum >= 0.5) nA++;
    if(randomNum <  0.5) nB++;
    nTOT++;

    // STEP1: calculate recentering for TPC event plane
    if(tpc_recenter_read_switch){
      if(doTPCptassocBin) {
        if(randomNum >= 0.5) { // subevent A
          Q2_p[ref9][region_vz]->Fill(0.5, x);
          Q2_p[ref9][region_vz]->Fill(1.5, y);
        }
        if(randomNum < 0.5) { // subevent B 
          Q2_m[ref9][region_vz]->Fill(0.5, x);
          Q2_m[ref9][region_vz]->Fill(1.5, y);
        }
      } else { // default: eta regions
        if(eta>0.){ // positive eta TPC region
          Q2_p[ref9][region_vz]->Fill(0.5, x);
          Q2_p[ref9][region_vz]->Fill(1.5, y);
        }
        if(eta<0.){ // negative eta TPC region
          Q2_m[ref9][region_vz]->Fill(0.5, x);
          Q2_m[ref9][region_vz]->Fill(1.5, y);
        }
      } // eta regions
    }

    //==================recentering procedure.
    // STEP2: read in recentering for TPC event plane
    //if(!tpc_recenter_read_switch){  // FIXME
    if(tpc_shift_read_switch){ 
      if(doTPCptassocBin) {
        //x -= GetTPCRecenterValue(randomNum, "x", ref9, region_vz);
        //y -= GetTPCRecenterValue(randomNum, "y", ref9, region_vz);

        x -= GetTPCRecenterValueNEW(randomNum, "x", ref9, region_vz);
        y -= GetTPCRecenterValueNEW(randomNum, "y", ref9, region_vz);
      } else {
        ////////////////////////////////////
        if(eta > 0){ // POSITIVE region
          // Method 1: from function
          x -= tpc_center_Qpx[ref9][region_vz];
          y -= tpc_center_Qpy[ref9][region_vz];

            // reading in from file
//          if(fCalibFile && doReadCalibFile){
//            // Method 2: from *.root calibration file
//            TProfile* hTPC_center_p = (TProfile*)fCalibFile->Get(Form("Q2_p%i_%i", ref9, region_vz));
//            hTPC_center_p->SetName("hTPC_center_p");
//            x -= hTPC_center_p->GetBinContent(1); // bin1 = x-vector
//            y -= hTPC_center_p->GetBinContent(2); // bin2 = y-vector
//            delete hTPC_center_p;
//          }

        } // positive eta
        if(eta < 0){ // NEGATIVE region
          // Method 1: from function
          x -= tpc_center_Qnx[ref9][region_vz];
          y -= tpc_center_Qny[ref9][region_vz];

            // reading in from file
//          if(fCalibFile && doReadCalibFile){
//            // Method 2: from *.root calibration file
//            TProfile* hTPC_center_m = (TProfile*)fCalibFile->Get(Form("Q2_m%i_%i", ref9, region_vz));
//            hTPC_center_m->SetName("hTPC_center_m");
//            x -= hTPC_center_m->GetBinContent(1); // bin1 = x-vector
//            y -= hTPC_center_m->GetBinContent(2); // bin2 = y-vector
//            delete hTPC_center_m;
//          }

        } // negative eta
      } // eta regions
    } // recenter tpc event plane

    // full TPC q-vectors
    Q2x += x;
    Q2y += y;

    if(doTPCptassocBin) {
      // construct TPC Q-vectors on a pt assoc bin basis
      if(randomNum >= 0.5) { // subevent A
        Q2x_p += x;
        Q2y_p += y;
        nA++;
      }
      if(randomNum < 0.5) { // subevent B
        Q2x_m += x;
        Q2y_m += y;
        nB++;
      }

      nTOT++;
    } else {
      // Positive / Negative Eta Regions
      if(eta>0.){ // positive eta region
        Q2x_p += x;
        Q2y_p += y;
        ntracksPOS++;
      }
      if(eta<0.){ // negative eta region
        Q2x_m += x;
        Q2y_m += y;
        ntracksNEG++;
      }
    } // eta regions

  } // track loop

  // test statements
  //cout<<"METHOD2... ntracksNEG = "<<ntracksNEG<<"  ntracksPOS = "<<ntracksPOS<<endl;
  //cout<<"Q2x_p = "<<Q2x_p<<"  Q2y_p = "<<Q2y_p<<"  Q2x_m = "<<Q2x_m<<"  Q2y_m = "<<Q2y_m<<endl;
  //cout<<"nA = "<<nA<<"  nB = "<<nB<<"  nTOT = "<<nTOT<<endl;
}

//
// Fill event plane resolution histograms
//_____________________________________________________________________________
void StEventPlaneMaker::CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1)
{
    //cout<<"bbc = "<<bbc<<"  zdc = "<<zdc<<"  tpc = "<<tpc<<"  tpcN = "<<tpcN<<"  tpcP = "<<tpcP<<"  bbc1 = "<<bbc1<<"  zdc1 = "<<zdc1<<endl;

    fProfV2Resolution[ref9]->Fill(12., TMath::Cos(2.*(tpcP - tpcN)));
    fProfV3Resolution[ref9]->Fill(12., TMath::Cos(3.*(tpcP - tpcN)));
    fProfV4Resolution[ref9]->Fill(12., TMath::Cos(4.*(tpcP - tpcN)));
    fProfV5Resolution[ref9]->Fill(12., TMath::Cos(5.*(tpcP - tpcN)));

    // fill the profiles for the resolution parameters
    // R2 resolution for 2nd order event plane
  if(!TMath::IsNaN(bbc)) {
    fProfV2Resolution[ref9]->Fill(2., TMath::Cos(2.*(bbc - tpc)));
    fProfV2Resolution[ref9]->Fill(3., TMath::Cos(2.*(bbc - zdc)));
    fProfV2Resolution[ref9]->Fill(4., TMath::Cos(2.*(tpc - bbc)));  // bin2
    fProfV2Resolution[ref9]->Fill(5., TMath::Cos(2.*(tpc - zdc)));  
    fProfV2Resolution[ref9]->Fill(6., TMath::Cos(2.*(zdc - tpc)));  // bin6
    fProfV2Resolution[ref9]->Fill(7., TMath::Cos(2.*(zdc - bbc)));  // bin3
    fProfV2Resolution[ref9]->Fill(8., TMath::Cos(2.*(bbc - tpcN)));
    fProfV2Resolution[ref9]->Fill(9., TMath::Cos(2.*(bbc - tpcP)));
    fProfV2Resolution[ref9]->Fill(10., TMath::Cos(2.*(zdc - tpcN)));
    fProfV2Resolution[ref9]->Fill(11., TMath::Cos(2.*(zdc - tpcP)));
    //fProfV2Resolution[ref9]->Fill(12., TMath::Cos(2.*(tpcP - tpcN)));
    fProfV2Resolution[ref9]->Fill(17., TMath::Cos(2.*(bbc1 - tpc)));
    fProfV2Resolution[ref9]->Fill(18., TMath::Cos(2.*(bbc1 - tpcN)));
    fProfV2Resolution[ref9]->Fill(19., TMath::Cos(2.*(bbc1 - tpcP)));
    fProfV2Resolution[ref9]->Fill(20., TMath::Cos(2.*(bbc1 - zdc1)));
    fProfV2Resolution[ref9]->Fill(21., TMath::Cos(2.*(zdc1 - tpc)));
    fProfV2Resolution[ref9]->Fill(22., TMath::Cos(2.*(zdc1 - tpcN)));
    fProfV2Resolution[ref9]->Fill(23., TMath::Cos(2.*(zdc1 - tpcP)));

    // R3 resolution for 2nd order event plane
    fProfV3Resolution[ref9]->Fill(2., TMath::Cos(3.*(bbc - tpc)));
    fProfV3Resolution[ref9]->Fill(3., TMath::Cos(3.*(bbc - zdc)));
    fProfV3Resolution[ref9]->Fill(4., TMath::Cos(3.*(tpc - bbc)));
    fProfV3Resolution[ref9]->Fill(5., TMath::Cos(3.*(tpc - zdc)));
    fProfV3Resolution[ref9]->Fill(6., TMath::Cos(3.*(zdc - tpc)));
    fProfV3Resolution[ref9]->Fill(7., TMath::Cos(3.*(zdc - bbc)));
    fProfV3Resolution[ref9]->Fill(8., TMath::Cos(3.*(bbc - tpcN)));
    fProfV3Resolution[ref9]->Fill(9., TMath::Cos(3.*(bbc - tpcP)));
    fProfV3Resolution[ref9]->Fill(10., TMath::Cos(3.*(zdc - tpcN)));
    fProfV3Resolution[ref9]->Fill(11., TMath::Cos(3.*(zdc - tpcP)));
    //fProfV3Resolution[ref9]->Fill(12., TMath::Cos(3.*(tpcP - tpcN)));
    fProfV3Resolution[ref9]->Fill(17., TMath::Cos(3.*(bbc1 - tpc)));
    fProfV3Resolution[ref9]->Fill(18., TMath::Cos(3.*(bbc1 - tpcN)));
    fProfV3Resolution[ref9]->Fill(19., TMath::Cos(3.*(bbc1 - tpcP)));
    fProfV3Resolution[ref9]->Fill(20., TMath::Cos(3.*(bbc1 - zdc1)));
    fProfV3Resolution[ref9]->Fill(21., TMath::Cos(3.*(zdc1 - tpc)));
    fProfV3Resolution[ref9]->Fill(22., TMath::Cos(3.*(zdc1 - tpcN)));
    fProfV3Resolution[ref9]->Fill(23., TMath::Cos(3.*(zdc1 - tpcP)));

    // R4 resolution for 2nd order event plane
    fProfV4Resolution[ref9]->Fill(2., TMath::Cos(4.*(bbc - tpc)));
    fProfV4Resolution[ref9]->Fill(3., TMath::Cos(4.*(bbc - zdc)));
    fProfV4Resolution[ref9]->Fill(4., TMath::Cos(4.*(tpc - bbc)));
    fProfV4Resolution[ref9]->Fill(5., TMath::Cos(4.*(tpc - zdc)));
    fProfV4Resolution[ref9]->Fill(6., TMath::Cos(4.*(zdc - tpc)));
    fProfV4Resolution[ref9]->Fill(7., TMath::Cos(4.*(zdc - bbc)));
    fProfV4Resolution[ref9]->Fill(8., TMath::Cos(4.*(bbc - tpcN)));
    fProfV4Resolution[ref9]->Fill(9., TMath::Cos(4.*(bbc - tpcP)));
    fProfV4Resolution[ref9]->Fill(10., TMath::Cos(4.*(zdc - tpcN)));
    fProfV4Resolution[ref9]->Fill(11., TMath::Cos(4.*(zdc - tpcP)));
    //fProfV4Resolution[ref9]->Fill(12., TMath::Cos(4.*(tpcP - tpcN)));
    fProfV4Resolution[ref9]->Fill(17., TMath::Cos(4.*(bbc1 - tpc)));
    fProfV4Resolution[ref9]->Fill(18., TMath::Cos(4.*(bbc1 - tpcN)));
    fProfV4Resolution[ref9]->Fill(19., TMath::Cos(4.*(bbc1 - tpcP)));
    fProfV4Resolution[ref9]->Fill(20., TMath::Cos(4.*(bbc1 - zdc1)));
    fProfV4Resolution[ref9]->Fill(21., TMath::Cos(4.*(zdc1 - tpc)));
    fProfV4Resolution[ref9]->Fill(22., TMath::Cos(4.*(zdc1 - tpcN)));
    fProfV4Resolution[ref9]->Fill(23., TMath::Cos(4.*(zdc1 - tpcP)));

    // R5 resolution for 2nd order event plane
    fProfV5Resolution[ref9]->Fill(2., TMath::Cos(5.*(bbc - tpc)));
    fProfV5Resolution[ref9]->Fill(3., TMath::Cos(5.*(bbc - zdc)));
    fProfV5Resolution[ref9]->Fill(4., TMath::Cos(5.*(tpc - bbc)));
    fProfV5Resolution[ref9]->Fill(5., TMath::Cos(5.*(tpc - zdc)));
    fProfV5Resolution[ref9]->Fill(6., TMath::Cos(5.*(zdc - tpc)));
    fProfV5Resolution[ref9]->Fill(7., TMath::Cos(5.*(zdc - bbc)));
    fProfV5Resolution[ref9]->Fill(8., TMath::Cos(5.*(bbc - tpcN)));
    fProfV5Resolution[ref9]->Fill(9., TMath::Cos(5.*(bbc - tpcP)));
    fProfV5Resolution[ref9]->Fill(10., TMath::Cos(5.*(zdc - tpcN)));
    fProfV5Resolution[ref9]->Fill(11., TMath::Cos(5.*(zdc - tpcP)));
    //fProfV5Resolution[ref9]->Fill(12., TMath::Cos(5.*(tpcP - tpcN)));
    fProfV5Resolution[ref9]->Fill(17., TMath::Cos(5.*(bbc1 - tpc)));
    fProfV5Resolution[ref9]->Fill(18., TMath::Cos(5.*(bbc1 - tpcN)));
    fProfV5Resolution[ref9]->Fill(19., TMath::Cos(5.*(bbc1 - tpcP)));
    fProfV5Resolution[ref9]->Fill(20., TMath::Cos(5.*(bbc1 - zdc1)));
    fProfV5Resolution[ref9]->Fill(21., TMath::Cos(5.*(zdc1 - tpc)));
    fProfV5Resolution[ref9]->Fill(22., TMath::Cos(5.*(zdc1 - tpcN)));
    fProfV5Resolution[ref9]->Fill(23., TMath::Cos(5.*(zdc1 - tpcP)));
  }

  // for the resolution of the combined vzero event plane, use two tpc halves as uncorrelated subdetectors
} 

//
//
//_____________________________________________________________________________
Double_t StEventPlaneMaker::CalculateEventPlaneChi(Double_t res) {
  // return chi for given resolution to combine event plane estimates from two subevents
  // see Phys. Rev. C no. CS6346 (http://arxiv.org/abs/nucl-ex/9805001)
  Double_t chi(2.), delta(1.), con((TMath::Sqrt(TMath::Pi()))/(2.*TMath::Sqrt(2)));
  for (Int_t i(0); i < 15; i++) {
    chi = ((con*chi*TMath::Exp(-chi*chi/4.)*(TMath::BesselI0(chi*chi/4.)+TMath::BesselI1(chi*chi/4.))) < res) ? chi + delta : chi - delta;
    delta = delta / 2.;
  }

  return chi;
}

//
//
//______________________________________________________________________
THnSparse* StEventPlaneMaker::NewTHnSparseEP(const char* name, UInt_t entries) {
  // generate new THnSparseD, axes are defined in GetDimParamsD()
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
      GetDimParamsEP(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }

    i++;
  }
  hnTitle += ";";

  return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseEP

//
// 
//______________________________________________________________________________________________
void StEventPlaneMaker::GetDimParamsEP(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

    case 0:
      label = "centrality 5% bin";
      nbins = 20;
      xmin = 0.;
      xmax = 100.;
      break;

    case 1:
      label = "track Qx-vector";
      nbins = 72;
      xmin = -90.;
      xmax = 90.;
      break;

    case 2:
      label = "track Qy-vector";
      nbins = 72;
      xmin = -90.;
      xmax = 90.;
      break;

    case 3:
      label = "event plane angle";
      nbins = 72;
      xmin = 0.;
      xmax = pi;
      break;

    case 4: // may delete this case
      label = "Jet p_{T} corrected with Rho";
      nbins = 50;
      xmin = -50.;
      xmax = 200.;
      break;

   }// end of switch
} // end of event plane sparse

//
//
//_____________________________________________________________________________
Double_t StEventPlaneMaker::GetEventPlaneAngle(TString det, Int_t order, Int_t correctin, TString subevt)
{
    // check for which type
    if(det.Contains("BBC")) {
      if(order == 1) return BBC_PSI1;
      if(order == 2) return BBC_PSI2;
    } else 
    if(det.Contains("ZDC")) { 
      if(order == 1) return ZDC_PSI1;
      if(order == 2) return ZDC_PSI2;
    } else 
    if(det.Contains("TPC")) {
      if(order == 2) {
        if(subevt.Contains("A")) return TPCA_PSI2;       
        if(subevt.Contains("B")) return TPCB_PSI2;
        return TPC_PSI2;
      }
    }

    // did not meet above criteria
    cout<<"Could not return event plane angle, check input !!"<<endl;
    return -999;
}
//
// return recentering value
//____________________________________________________________________________________________________________________
Double_t StEventPlaneMaker::GetTPCRecenterValue(Double_t randomNum, TString coordinate, Int_t ref9, Int_t region_vz) { 
    //==================recentering procedure.
    // STEP2: read in recentering for TPC event plane
        if(randomNum >= 0.5) { // subevent A
          switch(fTPCEPmethod) {
            // also have option of tpc_center_Qpx_bin0_Method1
            case kRemoveEtaStrip:
              if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
                if(fJetType == kFullJet) {
                  if(fJetRad == 0.4) { // Jet radius: R=0.4
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.4

                  if(fJetRad == 0.3) { // Jet radius: R=0.3
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1_R03_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.3

                  if(fJetRad == 0.2) { // Jet radius: R=0.2
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1_R02_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.2
                } // full jets

                if(fJetType == kChargedJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1ch_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // charged jets
              }
              break;

            case kRemoveEtaPhiCone:
              if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
                if(fJetType == kFullJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method2_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // full jets

                if(fJetType == kChargedJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method2ch_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // charged jets
              }

              if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qpx_bin0[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qpy_bin0[ref9][region_vz];
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qpx_bin1[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qpy_bin1[ref9][region_vz];
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qpx_bin2[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qpy_bin2[ref9][region_vz];
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qpx_bin3[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qpy_bin3[ref9][region_vz];
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qpx_bin4[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qpy_bin4[ref9][region_vz];
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              }
              break;

            case kRemoveLeadingJetConstituents:
              // dothis
              break;

            case kRemoveEtaPhiConeLeadSub:
              // dothis
              break;

            case kRemoveLeadingSubJetConstituents:
              // dothis
              break;

            default:
              // this is a default, but should never occur..
              break;

          } // METHOD switch

        }  // rand >= 0.5

        if(randomNum < 0.5) { // subevent B             
          switch(fTPCEPmethod) {
            case kRemoveEtaStrip:
              if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
                if(fJetType == kFullJet) {
                  if(fJetRad == 0.4) { // Jet radius: R=0.4
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  }

                  if(fJetRad == 0.3) { // Jet radius: R=0.3
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1_R03_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1_R03_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1_R03_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.3  

                  if(fJetRad == 0.2) { // Jet radius: R=0.2
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1_R02_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.2  
                } // full jets

                if(fJetType == kChargedJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1ch_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // charged jets       
              } // Run14
              break;

            case kRemoveEtaPhiCone:
              if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
                if(fJetType == kFullJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method2_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method2_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method2_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }

                } // full jets

                if(fJetType == kChargedJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method2ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method2ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method2ch_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // charged jets
              }

              if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qnx_bin0[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qny_bin0[ref9][region_vz];
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qnx_bin1[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qny_bin1[ref9][region_vz];
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qnx_bin2[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qny_bin2[ref9][region_vz];
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qnx_bin3[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qny_bin3[ref9][region_vz];
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  if(coordinate.Contains("x")) return tpc_center_Qnx_bin4[ref9][region_vz];
                  if(coordinate.Contains("y")) return tpc_center_Qny_bin4[ref9][region_vz];
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              }

              break;

            case kRemoveLeadingJetConstituents:
              // dothis
              break;

            case kRemoveEtaPhiConeLeadSub:
              // dothis
              break;

            case kRemoveLeadingSubJetConstituents:
              // dothis
              break;

            default:
              // this is a default, but should never occur..
              break;

          }  // METHOD switch
        } // pt < 0.5

    return 0.;
}
//
// return recentering value - new
//____________________________________________________________________________________________________________________
Double_t StEventPlaneMaker::GetTPCRecenterValueNEW(Double_t randomNum, TString coordinate, Int_t ref9, Int_t region_vz) {
    //==================recentering procedure.
    // STEP2: read in recentering for TPC event plane
        if(randomNum >= 0.5) { // subevent A
          switch(fTPCEPmethod) {
            // also have option of tpc_center_Qpx_bin0_Method1
            case kRemoveEtaStrip:
              if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
                if(fJetType == kFullJet) {
/*
                  if(fJetRad == 0.4) { // Jet radius: R=0.4
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.4
*/

                  if(fJetRad == 0.3) { // Jet radius: R=0.3
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.3

/*
                  if(fJetRad == 0.2) { // Jet radius: R=0.2
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1_R02_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.2
*/
                } // full jets

/*
                if(fJetType == kChargedJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin0_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin0_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin1_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin1_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin2_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin2_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin3_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin3_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qpx_bin4_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qpy_bin4_Method1ch_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // charged jets
*/

              }
              break;

            default:
              // this is a default, but should never occur..
              break;

          } // METHOD switch

        }  // rand >= 0.5

        if(randomNum < 0.5) { // subevent B             
          switch(fTPCEPmethod) {
            case kRemoveEtaStrip:
              if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
                if(fJetType == kFullJet) {

/*
                  if(fJetRad == 0.4) { // Jet radius: R=0.4
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  }
*/

                  if(fJetRad == 0.3) { // Jet radius: R=0.3
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1_R03_Run14_NEW[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1_R03_Run14_NEW[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.3  

/*
                  if(fJetRad == 0.2) { // Jet radius: R=0.2
                    if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1_R02_Run14[ref9][region_vz];
                    } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                      if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1_R02_Run14[ref9][region_vz];
                      if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1_R02_Run14[ref9][region_vz];
                    } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                  } // R=0.2  
*/

                } // full jets

/*
                if(fJetType == kChargedJet) {
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin0_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin0_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin1_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin1_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin2_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin2_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin3_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin3_Method1ch_Run14[ref9][region_vz];
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    if(coordinate.Contains("x")) return tpc_center_Qnx_bin4_Method1ch_Run14[ref9][region_vz];
                    if(coordinate.Contains("y")) return tpc_center_Qny_bin4_Method1ch_Run14[ref9][region_vz];
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // charged jets       
*/

              } // Run14
              break;

            default:
              // this is a default, but should never occur..
              break;

          }  // METHOD switch
        } // pt < 0.5

    return 0.;
}
//
// return shifting value
//____________________________________________________________________________________________________________________
Double_t StEventPlaneMaker::GetTPCShiftingValue(Double_t tPhi_rcd, Int_t nharm, Int_t ref9, Int_t region_vz) {
/*
  if(tpc_apply_corr_switch) { // FIXME: file needs to exist and need to have ran recentering + shift prior
    // loop over harmonics
    for(int nharm = 1; nharm < 21; nharm++){
      // Method 1: load from *.h file function
      // perform 'shift' to TPC event plane angle
      // KEEP in mind, the naming convent here means nothing for functions: tpc_shift_N and tpc_shift_P,
      // they are corresponing to Bn and An components above!
      // so hTPC_shft_N and hTPC_shift_P also have misleading names
      //shift_delta_psi += (shift_A[ref9][region_vz][z-1]*cos(2*z*psi2_rcd) + shift_B[ref9][region_vz][z-1]*sin(2*z*psi2_rcd));
      //tpc_delta_psi += (tpc_shift_N[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
      //                  tpc_shift_P[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
*/

        // TPC event plane calculation method
        switch(fTPCEPmethod) {
          case kRemoveEtaStrip: // added new usage to this switch on May31
            // run flag switch
            if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
              // jet type switch
              if(fJetType == kFullJet) {
                // jet size switch

                if(fJetRad == 0.4) { // Jet radius: R=0.4
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    return (tpc_shift_N_bin0_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin0_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    return (tpc_shift_N_bin1_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin1_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    return (tpc_shift_N_bin2_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin2_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    return (tpc_shift_N_bin3_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin3_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    return (tpc_shift_N_bin4_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin4_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // R=0.4

                if(fJetRad == 0.3) { // Jet radius: R=0.3
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    return (tpc_shift_N_bin0_Method1_R03_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin0_Method1_R03_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    return (tpc_shift_N_bin1_Method1_R03_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin1_Method1_R03_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    return (tpc_shift_N_bin2_Method1_R03_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin2_Method1_R03_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    return (tpc_shift_N_bin3_Method1_R03_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin3_Method1_R03_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    return (tpc_shift_N_bin4_Method1_R03_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin4_Method1_R03_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // R=0.3

                if(fJetRad == 0.2) { // Jet radius: R=0.2
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    return (tpc_shift_N_bin0_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin0_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    return (tpc_shift_N_bin1_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin1_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    return (tpc_shift_N_bin2_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin2_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    return (tpc_shift_N_bin3_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin3_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    return (tpc_shift_N_bin4_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin4_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // R=0.2
              } // full jets

              if(fJetType == kChargedJet) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  return (tpc_shift_N_bin0_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin0_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  return (tpc_shift_N_bin1_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin1_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  return (tpc_shift_N_bin2_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin2_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  return (tpc_shift_N_bin3_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin3_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  return (tpc_shift_N_bin4_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin4_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              } // charged jets
            } // RUN 14

            if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
              // FIXME - I believe I have these, or perhaps its standard usage and not pt associated bin
/*
              if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                return (tpc_shift_N_bin0_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin0_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                return (tpc_shift_N_bin1_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin1_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                return (tpc_shift_N_bin2_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin2_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                return (tpc_shift_N_bin3_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin3_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                return (tpc_shift_N_bin4_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin4_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
*/
            } // RUN 16

            break;

          case kRemoveEtaPhiCone:
            if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
              if(fJetType == kFullJet) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  return (tpc_shift_N_bin0_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin0_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  return (tpc_shift_N_bin1_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin1_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  return (tpc_shift_N_bin2_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin2_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  return (tpc_shift_N_bin3_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin3_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  return (tpc_shift_N_bin4_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin4_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              } // full jets

              if(fJetType == kChargedJet) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  return (tpc_shift_N_bin0_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin0_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  return (tpc_shift_N_bin1_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin1_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  return (tpc_shift_N_bin2_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin2_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  return (tpc_shift_N_bin3_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin3_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  return (tpc_shift_N_bin4_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin4_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              } // charged jets
            }

            if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
              if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                return (tpc_shift_N_bin0[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin0[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                return (tpc_shift_N_bin1[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin1[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                return (tpc_shift_N_bin2[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin2[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                return (tpc_shift_N_bin3[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin3[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                return (tpc_shift_N_bin4[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin4[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
            }
            break;

          case kRemoveLeadingJetConstituents:
            // dothis
            break;

          case kRemoveEtaPhiConeLeadSub:
            // dothis
            break;

          case kRemoveLeadingSubJetConstituents:
            // dothis
            break;

          default:
            // this is a default, but should never occur..
            break;

        } // METHOD switch

  return 0.0;
}
//
// return shifting value - NEW
//____________________________________________________________________________________________________________________
Double_t StEventPlaneMaker::GetTPCShiftingValueNEW(Double_t tPhi_rcd, Int_t nharm, Int_t ref9, Int_t region_vz) {
/*
  if(tpc_apply_corr_switch) { // FIXME: file needs to exist and need to have ran recentering + shift prior
    // loop over harmonics
    for(int nharm = 1; nharm < 21; nharm++){
      // Method 1: load from *.h file function
      // perform 'shift' to TPC event plane angle
      // KEEP in mind, the naming convent here means nothing for functions: tpc_shift_N and tpc_shift_P,
      // they are corresponing to Bn and An components above!
      // so hTPC_shft_N and hTPC_shift_P also have misleading names
      //shift_delta_psi += (shift_A[ref9][region_vz][z-1]*cos(2*z*psi2_rcd) + shift_B[ref9][region_vz][z-1]*sin(2*z*psi2_rcd));
      //tpc_delta_psi += (tpc_shift_N[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
      //                  tpc_shift_P[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
*/

        // TPC event plane calculation method
        switch(fTPCEPmethod) {
          case kRemoveEtaStrip: // added new usage to this switch on May31
            // run flag switch
            if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
              // jet type switch
              if(fJetType == kFullJet) {
                // jet size switch

/*
                if(fJetRad == 0.4) { // Jet radius: R=0.4
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    return (tpc_shift_N_bin0_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin0_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    return (tpc_shift_N_bin1_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin1_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    return (tpc_shift_N_bin2_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin2_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    return (tpc_shift_N_bin3_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin3_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    return (tpc_shift_N_bin4_Method1_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin4_Method1_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // R=0.4
*/

                if(fJetRad == 0.3) { // Jet radius: R=0.3
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    return (tpc_shift_N_bin0_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin0_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    return (tpc_shift_N_bin1_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin1_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    return (tpc_shift_N_bin2_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin2_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    return (tpc_shift_N_bin3_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin3_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    return (tpc_shift_N_bin4_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin4_Method1_R03_Run14_NEW[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // R=0.3

/*
                if(fJetRad == 0.2) { // Jet radius: R=0.2
                  if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                    return (tpc_shift_N_bin0_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin0_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                    return (tpc_shift_N_bin1_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin1_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                    return (tpc_shift_N_bin2_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin2_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                    return (tpc_shift_N_bin3_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin3_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                    return (tpc_shift_N_bin4_Method1_R02_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                            tpc_shift_P_bin4_Method1_R02_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                  } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
                } // R=0.2
*/

              } // full jets

/*
              if(fJetType == kChargedJet) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  return (tpc_shift_N_bin0_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin0_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  return (tpc_shift_N_bin1_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin1_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  return (tpc_shift_N_bin2_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin2_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  return (tpc_shift_N_bin3_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin3_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  return (tpc_shift_N_bin4_Method1ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin4_Method1ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              } // charged jets
*/
            } // RUN 14

            if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
/*
              if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                return (tpc_shift_N_bin0_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin0_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                return (tpc_shift_N_bin1_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin1_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                return (tpc_shift_N_bin2_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin2_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                return (tpc_shift_N_bin3_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin3_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                return (tpc_shift_N_bin4_Method1_Run16[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin4_Method1_Run16[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
*/
            } // RUN 16

            break;

/*
          case kRemoveEtaPhiCone:
            if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) {
              if(fJetType == kFullJet) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  return (tpc_shift_N_bin0_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin0_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  return (tpc_shift_N_bin1_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin1_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  return (tpc_shift_N_bin2_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin2_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  return (tpc_shift_N_bin3_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin3_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  return (tpc_shift_N_bin4_Method2_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin4_Method2_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              } // full jets

              if(fJetType == kChargedJet) {
                if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                  return (tpc_shift_N_bin0_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin0_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                  return (tpc_shift_N_bin1_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin1_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                  return (tpc_shift_N_bin2_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin2_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                  return (tpc_shift_N_bin3_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin3_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                  return (tpc_shift_N_bin4_Method2ch_Run14[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                          tpc_shift_P_bin4_Method2ch_Run14[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
                } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
              } // charged jets
            }

            if(fRunFlag == StJetFrameworkPicoBase::Run16_AuAu200) {
              if(fTPCptAssocBin == 0) {         // 0.20-0.50 GeV
                return (tpc_shift_N_bin0[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin0[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 1) {  // 0.50-1.00 GeV
                return (tpc_shift_N_bin1[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin1[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 2) {  // 1.00-1.50 GeV
                return (tpc_shift_N_bin2[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin2[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 3) {  // 1.50-2.00 GeV
                return (tpc_shift_N_bin3[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin3[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else if(fTPCptAssocBin == 4) {  // 2.00-20.0 GeV
                return (tpc_shift_N_bin4[ref9][region_vz][nharm-1] * cos(2*nharm*tPhi_rcd) +
                        tpc_shift_P_bin4[ref9][region_vz][nharm-1] * sin(2*nharm*tPhi_rcd));
              } else { cout<<"NOT CONFIGURED PROPERLY, please select pt assoc bin!"<<endl; }
            }
            break;

          case kRemoveLeadingJetConstituents:
            // dothis
            break;

          case kRemoveEtaPhiConeLeadSub:
            // dothis
            break;

          case kRemoveLeadingSubJetConstituents:
            // dothis
            break;
*/

          default:
            // this is a default, but should never occur..
            break;

        } // METHOD switch

  return 0.0;
}
