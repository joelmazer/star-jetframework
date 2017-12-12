#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

// some includes
//#include <initializer_list>

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
class StEventPoolManager;
class StCalibContainer;
class StEPFlattener;

//class StMyAnalysisMaker : public StMaker {
class StMyAnalysisMaker : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugMixedEvents,
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality
    };

    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments, double minJetPtCut, double trkbias, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StMyAnalysisMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();
    void    WriteEventPlaneHistograms();

    // ep stuff - Nov15
    void    InitParameters();
   
    // THnSparse Setup
    virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
    virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
    virtual THnSparse*      NewTHnSparseFCorr(const char* name, UInt_t entries);
    virtual void            GetDimParamsCorr(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
    virtual THnSparse*      NewTHnSparseEP(const char* name, UInt_t entries);
    virtual void            GetDimParamsEP(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

    static TString GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray* partCont, TClonesArray* clusCont, TString tag);

    // TClonesArrays function returners of analysis objects
    TClonesArray* jets() const { return mJets; }
    TClonesArray* tracks() const { return mTracks; }
    TClonesArray* towers() const { return mTowers; }
    TClonesArray* particles() const { return mParticles; }

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }

    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }

    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

    // event mixing - setters
    virtual void            SetEventMixing(Int_t yesno)	       { fDoEventMixing=yesno; }
    virtual void            SetMixingTracks(Int_t tracks)      { fMixingTracks = tracks; }
    virtual void            SetNMixedTr(Int_t nmt)             { fNMIXtracks = nmt; }
    virtual void            SetNMixedEvt(Int_t nme)            { fNMIXevents = nme; }

    // mixed selection - setters
    virtual void            SetTriggerEventType(UInt_t te)       { fTriggerEventType = te; }
    virtual void            SetMixedEventType(UInt_t me)         { fMixingEventType = me; }
    virtual void            SetCentBinSize(Int_t centbins)       { fCentBinSize = centbins; }
    virtual void            SetReduceStatsCent(Int_t red)        { fReduceStatsCent = red; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Int_t effcorr)          { fDoEffCorr = effcorr; }

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    // event plane
    StJet*                  GetLeadingJet(StRhoParameter* eventRho = 0x0);
    StJet*                  GetSubLeadingJet(StRhoParameter* eventRho = 0x0);
    virtual void            SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
    virtual void            SetEventPlaneTrackWeight(Int_t weight)          {fTrackWeight = weight; }
    virtual void            SetEventPlaneMaxTrackPtCut(Double_t m)          {fEventPlaneMaxTrackPtCut = m; }  
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

    // Where to read calib object with EP calibration if not default
    void                   SetEPcalibFileName(const TString filename) {fEPcalibFileName = filename; } 
    void                   SetOutFileNameEP(const TString epout) {mOutNameEP = epout; }

  protected:
    Int_t                  GetCentBin(Int_t cent, Int_t nBin) const;             // centrality bin
    Double_t               RelativePhi(Double_t mphi,Double_t vphi) const;       // relative jet track angle
    Double_t               RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    TH1*                   FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    TH1*                   FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
    Double_t               GetReactionPlane();                                   // get reaction plane angle
    void                   GetEventPlane(Bool_t flattenEP);                      // get event plane / flatten and fill histos 
    Bool_t                 AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function
    Bool_t                 AcceptJet(StJet *jet);           // jets accept cuts function
    Bool_t                 DoComparison(int myarr[], int elems);
    void                   SetSumw2(); // set errors weights 
    void                   SetEPSumw2(); // set errors weights for event plane histograms
    //Double_t               EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                   CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP);

    // Added from Liang
    void                   QvectorCal(int ref9, int region_vz, int n);
    Int_t                  EventPlaneCal(int ref9, int region_vz, int n);
    Int_t                  BBC_EP_Cal(int ref9, int region_vz, int n); //refmult, the region of vz, and order of EP
    Int_t                  ZDC_EP_Cal(int ref9, int region_vz, int n);
    Double_t                BBC_GetPhi(int e_w,int iTile); //east == 0
    Double_t                ZDCSMD_GetPosition(int id_order,int eastwest,int verthori,int strip);
    Int_t                  GetRunNo(int runid);
    Int_t                  GetVzRegion(double Vz);

    // switches
    Bool_t                 doUsePrimTracks;         // primary track switch
    Int_t                  fDebugLevel;             // debug printout level
    Bool_t                 doPrintEventCounter;     // print event # switch
    Int_t                  fRunFlag;                // Run Flag numerator value
    Int_t                  fCentralityDef;          // Centrality Definition enumerator value
    Int_t                  fDoEffCorr;              // efficiency correction to tracks
    Bool_t                 fCorrJetPt;              // correct jet pt by rho
    Bool_t                 doEventPlaneRes;         // event plane resolution switch

    // cuts
    Double_t               fMinPtJet;               // min jet pt to keep jet in output
    Double_t               fTrackBias;              // high pt track in jet bias
    Double_t               fJetRad;                 // jet radius
    Double_t               fEventZVtxMinCut;        // min event z-vertex cut
    Double_t               fEventZVtxMaxCut;        // max event z-vertex cut
    Double_t               fTrackPtMinCut;          // min track pt cut
    Double_t               fTrackPtMaxCut;          // max track pt cut
    Double_t               fTrackPhiMinCut;         // min track phi cut
    Double_t               fTrackPhiMaxCut;         // max track phi cut
    Double_t               fTrackEtaMinCut;         // min track eta cut
    Double_t               fTrackEtaMaxCut;         // max track eta cut
    Double_t               fTrackDCAcut;            // max track dca cut
    Int_t                  fTracknHitsFit;          // requirement for track hits
    Double_t               fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax

    // event mixing
    Int_t          fDoEventMixing;              // switch ON/off event mixing
    Int_t          fMixingTracks;               // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t          fNMIXtracks;                 // MIN # of mixing track in pool before performing mixing
    Int_t          fNMIXevents;                 // MIN # of mixing events in pool before performing mixing
    Int_t          fCentBinSize;                // centrality bin size of mixed event pools
    Int_t          fReduceStatsCent;            // bins to use for reduced statistics of sparse

    // centrality    
    Double_t       fCentralityScaled;           // scaled by 5% centrality 
    Int_t          ref16;                       // multiplicity bin (16)
    Int_t          ref9;                        // multiplicity bin (9)

    // event
    Double_t        Bfield;                      // event Bfield
    StThreeVectorF  mVertex;                     // event vertex 3-vector
    Double_t        zVtx;                        // z-vertex component

    // event selection types
    UInt_t         fTriggerEventType;           // Physics selection of event used for signal
    UInt_t         fMixingEventType;            // Physics selection of event used for mixed event
    Int_t          fEmcTriggerArr[7];           // EMCal triggers array: used to select signal and do QA

    // used for event plane calculation and resolution
    StJet*         fLeadingJet;//! leading jet
    StJet*         fSubLeadingJet;//! subleading jet
    Float_t        fExcludeLeadingJetsFromFit;  // exclude n leading jets from fit
    Int_t          fTrackWeight;                // track weight for Q-vector summation
    Double_t       fEventPlaneMaxTrackPtCut;    // max track pt cut for event plane calculation
    Bool_t         phi_shift_switch;            // phi shift - for TPC: NOT USING!
    Bool_t         tpc_recenter_read_switch;    // tpc recenter reader
    Bool_t         tpc_shift_read_switch;       // tpc shift reader
    Bool_t         tpc_apply_corr_switch;       // tpc apply final corrections
    Bool_t         zdc_recenter_read_switch;    // zdc recentering switch
    Bool_t         zdc_shift_read_switch;       // zdc shift reader
    Bool_t         zdc_apply_corr_switch;       // zdc apply final corrections
    Bool_t         bbc_recenter_read_switch;    // bbc recentering switch
    Bool_t         bbc_shift_read_switch;       // bbc shift reader
    Bool_t         bbc_apply_corr_switch;       // bbc apply final corrections
    Int_t          fHistCentBinMin;             // min centrality bin for histogram loop
    Int_t          fHistCentBinMax;             // max centrality bin for histogram loop
    Int_t          fHistZvertBinMin;            // min z-vertex bin for histogram loop
    Int_t          fHistZvertBinMax;            // min z-vertex bin for histogram loop

    // global variables used with TPC event plane corrections
    Double_t       Q2x_raw;
    Double_t       Q2y_raw;
    Double_t       Q2x_p;
    Double_t       Q2x_m;
    Double_t       Q2y_p;
    Double_t       Q2y_m;
    Double_t       Q2x;
    Double_t       Q2y;
    Double_t       TPC_PSI2;
    Double_t       BBC_PSI2;
    Double_t       ZDC_PSI2;
    Double_t       PSI2;
    Double_t       RES;

    // event pool
    TClonesArray          *CloneAndReduceTrackList();
    StEventPoolManager    *fPoolMgr;//!  // event pool Manager object

    // clonesarray collections of tracks and jets
    TClonesArray          *fTracksME;//! track collection to slim down for mixed events
    TClonesArray          *fJets;//! jet collection

  private:
    //void                   GetVZEROEventPlane(Bool_t isFlatten);
    Double_t               ApplyFlatteningTPCn(Double_t phi, Double_t c);
    Double_t               ApplyFlatteningTPCp(Double_t phi, Double_t c);
    Double_t               ApplyFlatteningBBC(Double_t phi, Double_t c);
    Double_t               ApplyFlatteningZDC(Double_t phi, Double_t c);

    Int_t                  fRunNumber;
    TString                fEPcalibFileName; 
    StCalibContainer      *fFlatContainer;
    StEPFlattener         *fTPCnFlat;
    StEPFlattener         *fTPCpFlat;
    StEPFlattener         *fBBCFlat;
    StEPFlattener         *fZDCFlat;
    Double_t               fEPTPCResolution;
    Double_t               fEPTPCn;
    Double_t               fEPTPCp;
    Double_t               fEPBBC;
    Double_t               fEPZDC;

    // PicoDstMaker and PicoDst object pointer
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StJetMakerTask   *JetMaker;
    StRho          *RhoMaker;

    // centrality objects
    StRefMultCorr* grefmultCorr;

    // TCloneArray of analysis objects
    TClonesArray   *mJets;
    TClonesArray   *mTracks;
    TClonesArray   *mTowers;
    TClonesArray   *mParticles;
   
    // output file name string 
    TString      mOutName;
    TString      mOutNameEP;

    TFile        *fCalibFile;
    TFile        *fCalibFile2;
    TFile        *fBBCcalibFile;
    TFile        *fZDCcalibFile;

/* 
    Int_t        mEventCounter;//!
    Int_t        mAllPVEventCounter;//!
    Int_t        mInputEventCounter;//!
*/

    // switches
    bool         doComments;

    // histograms
    TH1F* hEventPlane;//!   
    TH1F* hEventPlane2pi;//!
    TH1F* hEventPlaneWeighted;//!
    TH2F* fHistEPTPCnAlt;//!
    TH2F* fHistEPTPCpAlt;//!
    TH2F* fHistEPTPCn;//!
    TH2F* fHistEPTPCp;//!
    TH2F* fHistEPBBC;//!
    TH2F* fHistEPZDC;//!
    TH2F* fHistEPTPCnFlatten;//!
    TH2F* fHistEPTPCpFlatten;//!
    TH2F* fHistEPBBCFlatten;//!
    TH2F* fHistEPZDCFlatten;//!
    TH1F* hEventZVertex;//!
    TH1F* hCentrality;//!
    TH1F* hMultiplicity;//!
    TH2F* hRhovsCent;//!
    
    // jet histos
    TH1F* hJetPt;//!
    TH1F* hJetCorrPt;//!
    TH1F* hJetPt2;//!
    TH1F* hJetE;//!
    TH1F* hJetEta;//!
    TH1F* hJetPhi;//!
    TH1F* hJetNEF;//!
    TH1F* hJetArea;//!
    TH1F* hJetTracksPt;//!
    TH1F* hJetTracksPhi;//!
    TH1F* hJetTracksEta;//!
    TH2F* hJetPtvsArea;//!

    // correlation histo
    TH2  *fHistJetHEtaPhi;//!

    // QA histos
    TH1  *fHistEventSelectionQA;//! 
    TH1  *fHistEventSelectionQAafterCuts;//!
    TH1  *hTriggerIds;//!
    TH1  *hEmcTriggers;//!
    TH1  *hMixEvtStatZVtx;//!
    TH1  *hMixEvtStatCent;//!
    TH2  *hMixEvtStatZvsCent;//!

    // event plane histograms for corrections and calculations
    TH2F *hZDCDis_W;//!
    TH2F *hZDCDis_E;//!
    TH2F *hBBCDis_W;//!
    TH2F *hBBCDis_E;//!
    TProfile *Q2_p[9][20];//!  // 15
    TProfile *Q2_m[9][20];//!  // 15
    TProfile *phishiftA[2][2][9];
    TProfile *phishiftB[2][2][9];
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
    TH1F *bbc_psi_raw;//!
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

    // THn Sparse's jet sparse
    THnSparse             *fhnJH;//!           // jet hadron events matrix
    THnSparse             *fhnMixedEvents;//!  // mixed events matrix
    THnSparse             *fhnCorr;//!         // sparse to get # jet triggers

    THnSparse             *fhnEP;//!           // event plane sparse

    // Rho objects
    StRhoParameter        *GetRhoFromEvent(const char *name);
    StRhoParameter        *fRho;//!<!          // event rho
    Double_t               fRhoVal;//!<!       // event rho value, same for local rho
    TString                fRhoName;///<       // rho name

    // maker names
    TString                fAnalysisMakerName;
    TString                fJetMakerName;
    TString                fRhoMakerName;
    TString                fEventMixerMakerName;

/*
    // counters
    Int_t GetEventCounter() {return mEventCounter;}
    Int_t GetAllPVEventCounter() {return mAllPVEventCounter;}
    Int_t GetInputEventCounter() {return mInputEventCounter;}
*/
                
    ClassDef(StMyAnalysisMaker, 1)
};
/*
template <typename T> bool is_in(const T& val, const std::initializer_list<T>& list)
{
    for (const auto& i : list) {
        if (val == i) {
            return true;
        }
    }
    return false;
}
*/
#endif
