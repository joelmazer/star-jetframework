#ifndef StEventPlaneMaker_h
#define StEventPlaneMaker_h

#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class THnSparse;
class TProfile;
class TString;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// my STAR classes
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;

//class StEventPlaneMaker : public StMaker {
class StEventPlaneMaker : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugMixedEvents,
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality,
      kDebugEventPlaneCalc
    };

    // enumerator for TPC event plane method
    enum fTPCEPmethodEnum {
      kRemoveNothing,
      kRemoveEtaStrip,
      kRemoveEtaPhiCone,
      kRemoveLeadingJetConstituents, // greater than 2 GeV
      kRemoveEtaStripLeadSub,
      kRemoveEtaPhiConeLeadSub,
      kRemoveLeadingSubJetConstituents // greater than 2 GeV
    };

    // detector enum
    enum fDetectorType {
      kNoDetector = 0,
      kBBC = 1,
      kZDC = 2, 
      kTPC = 3
    };

    StEventPlaneMaker(const char *name, StPicoDstMaker *picoMaker, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StEventPlaneMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteEventPlaneHistograms();

    // ep stuff - Nov15
    void    InitParameters();
   
    // THnSparse Setup
    virtual THnSparse*      NewTHnSparseEP(const char* name, UInt_t entries);
    virtual void            GetDimParamsEP(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt) { fMaxEventTrackPt = mxpt; }

    // jet switches
    virtual void            SetJetType(Int_t jt)               { fJetType          = jt;}    // jet type (full, charged, neutral)
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    
    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }
    virtual void            SetTriggerToUse(UInt_t ttu)        { fTriggerToUse = ttu; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Int_t effcorr)        { fDoEffCorr = effcorr; }

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)        { fCorrJetPt = cpt; }

    // event plane
    virtual void            SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
    virtual void            SetEventPlaneTrackWeight(Int_t weight)          {fTrackWeight = weight; }
    virtual void            SetEventPlaneMaxTrackPtCut(Double_t m)          {fEventPlaneMaxTrackPtCut = m; }  
    virtual void            SetTPCEventPlaneMethod(Int_t tm)                {fTPCEPmethod = tm; }
    virtual void            SetPhiShift(Bool_t ps)                          {phi_shift_switch = ps; }
    virtual void            SetTPCRecenterRead(Bool_t trc)                  {tpc_recenter_read_switch = trc; }
    virtual void            SetTPCShiftRead(Bool_t ts)                      {tpc_shift_read_switch = ts; }
    virtual void            SetTPCApplyCorrections(Bool_t tac)              {tpc_apply_corr_switch = tac; }
    virtual void            SetZDCrecentering(Bool_t zrc)                   {zdc_recenter_read_switch = zrc; }
    virtual void            SetZDCShiftRead(Bool_t zs)                      {zdc_shift_read_switch = zs; }
    virtual void            SetZDCApplyCorrections(Bool_t zac)              {zdc_apply_corr_switch = zac; }
    virtual void            SetBBCrecentering(Bool_t brc)                   {bbc_recenter_read_switch = brc; }
    virtual void            SetBBCShiftRead(Bool_t bs)                      {bbc_shift_read_switch = bs; }
    virtual void            SetBBCApplyCorrections(Bool_t bac)              {bbc_apply_corr_switch = bac; }
    virtual void            SetHistBinLimitsCenZvert(Int_t cmin, Int_t cmax, Int_t zmin, Int_t zmax)   { fHistCentBinMin = cmin; fHistCentBinMax = cmax; fHistZvertBinMin = zmin; fHistZvertBinMax = zmax; }
    virtual void            SetdoEventPlaneRes(Bool_t depr)                 {doEventPlaneRes = depr; }
    virtual void            SetdoEPTPCptAssocMethod(Bool_t ptbin)           {doTPCptassocBin = ptbin; }
    virtual void            SetEPTPCptAssocBin(Int_t pb)                    {fTPCptAssocBin = pb; }

    // Where to read calib object with EP calibration if not default
    void                    SetEPcalibFileName(TString filename)            {fEPcalibFileName = filename; } 
    void                    SetOutFileNameEP(TString epout)                 {mOutNameEP = epout; }
    virtual void            SetdoReadCalibFile(Bool_t rc)                   {doReadCalibFile = rc; } 

    // get functions:
    Double_t                GetTPCEP()                { return TPC_PSI2; }

  protected:
    TH1                    *FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    TH1                    *FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
    void                    GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin);// get event plane / flatten and fill histos 
    Bool_t                  AcceptJet(StJet *jet);           // jets accept cuts function
    void                    SetEPSumw2(); // set errors weights for event plane histograms
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1);
    static Double_t         CalculateEventPlaneChi(Double_t res);
    Double_t                GetEventPlaneAngle(TString det, Int_t order, Int_t correctin, TString subevt);
    Double_t                GetTPCRecenterValue(Double_t randomNum, TString coordinate, Int_t ref9, Int_t region_vz);
    Double_t                GetTPCRecenterValueNEW(Double_t randomNum, TString coordinate, Int_t ref9, Int_t region_vz);
    Double_t                GetTPCShiftingValue(Double_t tPhi_rcd, Int_t nharm, Int_t ref9, Int_t region_vz);
    Double_t                GetTPCShiftingValueNEW(Double_t tPhi_rcd, Int_t nharm, Int_t ref9, Int_t region_vz);

    // Added from Liang
    void                    QvectorCal(int ref9, int region_vz, int n, int ptbin);
    Int_t                   EventPlaneCal(int ref9, int region_vz, int n, int ptbin);
    Int_t                   BBC_EP_Cal(int ref9, int region_vz, int n); //refmult, the region of vz, and order of EP
    Int_t                   ZDC_EP_Cal(int ref9, int region_vz, int n);
    Double_t                BBC_GetPhi(int e_w,int iTile); //east == 0
    Double_t                ZDCSMD_GetPosition(int id_order,int eastwest,int verthori,int strip);
    Int_t                   GetRunNo(int runid);
    Int_t                   GetVzRegion(double Vz);

    // switches
    Int_t                   fDoEffCorr;              // efficiency correction to tracks
    Bool_t                  doEventPlaneRes;         // event plane resolution switch
    Bool_t                  doTPCptassocBin;         // TPC event plane calculated on a pt assoc bin basis
    Int_t                   fTPCptAssocBin;          // pt associated bin to calculate event plane for
    Bool_t                  doReadCalibFile;         // read calibration file switch

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    UInt_t                  fTriggerToUse;               // trigger to use for analysis
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    // used for event plane calculation and resolution
    Double_t                fEventPlaneMaxTrackPtCut;    // max track pt cut for event plane calculation
    Int_t                   fTPCEPmethod;                // TPC event plane calculation method
    Bool_t                  phi_shift_switch;            // phi shift - for TPC: NOT USING!
    Bool_t                  tpc_recenter_read_switch;    // tpc recenter reader
    Bool_t                  tpc_shift_read_switch;       // tpc shift reader
    Bool_t                  tpc_apply_corr_switch;       // tpc apply final corrections
    Bool_t                  zdc_recenter_read_switch;    // zdc recentering switch
    Bool_t                  zdc_shift_read_switch;       // zdc shift reader
    Bool_t                  zdc_apply_corr_switch;       // zdc apply final corrections
    Bool_t                  bbc_recenter_read_switch;    // bbc recentering switch
    Bool_t                  bbc_shift_read_switch;       // bbc shift reader
    Bool_t                  bbc_apply_corr_switch;       // bbc apply final corrections
    Int_t                   fHistCentBinMin;             // min centrality bin for histogram loop
    Int_t                   fHistCentBinMax;             // max centrality bin for histogram loop
    Int_t                   fHistZvertBinMin;            // min z-vertex bin for histogram loop
    Int_t                   fHistZvertBinMax;            // min z-vertex bin for histogram loop

    // global variables used with TPC event plane corrections
    Double_t                Q2x_raw;
    Double_t                Q2y_raw;
    Double_t                Q2x_p;
    Double_t                Q2x_m;
    Double_t                Q2y_p;
    Double_t                Q2y_m;
    Double_t                Q2x;
    Double_t                Q2y;
    Double_t                TPC_PSI2;
    Double_t                TPCA_PSI2;
    Double_t                TPCB_PSI2;
    Double_t                BBC_PSI2;
    Double_t                ZDC_PSI2;
    Double_t                BBC_PSI1;
    Double_t                ZDC_PSI1;
    Double_t                PSI2;
    Double_t                RES;
    // temp (possibly)
    Double_t                TPC_raw_comb;
    Double_t                TPC_raw_neg;
    Double_t                TPC_raw_pos;
    Double_t                BBC_raw_comb;
    Double_t                BBC_raw_east;
    Double_t                BBC_raw_west;
    Double_t                ZDC_raw_comb;
    Double_t                ZDC_raw_east;
    Double_t                ZDC_raw_west;

  private:
    Int_t                   fRunNumber;
    TString                 fEPcalibFileName; 
    Double_t                fEPTPCResolution;
    Double_t                fEPTPCn;
    Double_t                fEPTPCp;
    Double_t                fEPTPC;
    Double_t                fEPBBC;
    Double_t                fEPZDC;

    TFile                  *fCalibFile;
    TFile                  *fCalibFile2;
    TFile                  *fBBCcalibFile;
    TFile                  *fZDCcalibFile;

    // histograms
    TH1F *hCentrality;//!
    TH1F *hCentralityEP;//!
    TH1F *hEventPlane;//!   
    TH2F *fHistEPTPCn;//!
    TH2F *fHistEPTPCp;//!
    TH2F *fHistEPBBC;//!
    TH2F *fHistEPZDC;//!
    TH1F *hTrackPhi[9];//!
    TH1F *hTrackPt[9];//!

    // QA histos
    TH1 *fHistEventSelectionQA;//! 
    TH1 *fHistEventSelectionQAafterCuts;//!
    TH1 *hTriggerIds;//!
    TH1 *hEmcTriggers;//!
    TH1 *hTPCepDebug;//!
    TH1 *hBBCepDebug;//!
    TH1 *hZDCepDebug;//!

    // event plane histograms for corrections and calculations
    TH2F *hZDCDis_W;//!
    TH2F *hZDCDis_E;//!
    TH2F *hBBCDis_W;//!
    TH2F *hBBCDis_E;//!
    TProfile *Q2_p[9][20];//!  // 15
    TProfile *Q2_m[9][20];//!  // 15
    TProfile *hTPC_shift_N[9][20];//! // 15
    TProfile *hTPC_shift_P[9][20];//! // 15
    TProfile *hBBC_shift_A[9][20];//! // 15
    TProfile *hBBC_shift_B[9][20];//! // 15
    TProfile *hZDC_shift_A[9][20];//! // new
    TProfile *hZDC_shift_B[9][20];//! // new
    TProfile *res_cen;//!
    TProfile *hZDC_center_ex;//!
    TProfile *hZDC_center_ey;//!
    TProfile *hZDC_center_wx;//!
    TProfile *hZDC_center_wy;//!
    TProfile *hBBC_center_ex;//!
    TProfile *hBBC_center_ey;//!
    TProfile *hBBC_center_wx;//!
    TProfile *hBBC_center_wy;//!
    TProfile *bbc_res;
    TH1F *zdc_psi;//!
    TH1F *checkbbc;//!
    TH2F *psi2_tpc_bbc;//!
    TH1F *bbc_psi_e;//!
    TH1F *bbc_psi_w;//!
    TH2F *bbc_psi_evw;//!
    TH1F *bbc_psi1_raw;//! // 1st order
    TH1F *bbc_psi_raw;//!  // 2nd order
    TH1F *bbc_psi_rcd;//!
    TH1F *bbc_psi_sft;//!
    TH1F *bbc_psi_fnl;//!

    //TH1F *zdc_psi;//! // declared already
    TProfile *zdc_res;//!
    TH1F *zdc_psi_e;//!
    TH1F *zdc_psi_w;//!
    TH2F *zdc_psi_evw;//!
    TH1F *zdc_psi_raw;//!
    TH1F *zdc_psi_rcd;//!
    TH1F *zdc_psi_sft;//!
    TH1F *zdc_psi_fnl;//!

    TProfile *tpc_res;//!
    TH1F *tpc_psi;//!
    TH1F *tpc_psi_N;//!
    TH1F *tpc_psi_P;//!
    TH2F *tpc_psi_NvP;//!
    TH1F *tpc_psi_raw;//!
    TH1F *tpc_psi_rcd;//!
    TH1F *tpc_psi_sft;//!
    TH1F *tpc_psi_fnl;//!

    TH1F *Psi2;//!
    TH1F *Psi2m;//!
    TH1F *Psi2p;//!
    TH1F *Delta_Psi2;//!
    TH1F *Delta_Psi2cyc;//!
    TH1F *Delta_Psi2old;//!
    TH1F *Shift_delta_psi2;//!
    TH1F *Psi2_rcd;//!
    TH1F *Psi2_final;//!
    TH1F *Psi2_final_raw;//!
    TH1F *Psi2_final_folded;//!

    TH2F *hTPCvsBBCep;//!
    TH2F *hTPCvsZDCep;//!
    TH2F *hBBCvsZDCep;//!

    // EP resoltuion profiles
    TProfile              *fProfV2Resolution[9];//! resolution parameters for v2
    TProfile              *fProfV3Resolution[9];//! resolution parameters for v3
    TProfile              *fProfV4Resolution[9];//! resolution parameters for v4
    TProfile              *fProfV5Resolution[9];//! resolution parameters for v5
//    TH1F                  *fDiffV2Resolution[9];//! difference of event plane angles for n=2
//    TH1F                  *fDiffV3Resolution[9];//! difference of event plane angles for n=3
//    TH1F                  *fDiffV4Resolution[9];//! difference of event plane angles for n=4
//    TH1F                  *fDiffV5Resolution[9];//! difference of event plane angles for n=5

    // maker names
    TString                fAnalysisMakerName;
//    TString                fJetMakerName;
//    TString                fRhoMakerName;
                
    ClassDef(StEventPlaneMaker, 2)
};
#endif
