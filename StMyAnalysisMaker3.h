#ifndef StMyAnalysisMaker3_h
#define StMyAnalysisMaker3_h

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
//class StEventPlaneMaker;

//class StMyAnalysisMaker3 : public StMaker {
class StMyAnalysisMaker3 : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugMixedEvents,
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality,
      kDebugEventPlaneCalc,
      kDebugJetvsEPtype
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

    StMyAnalysisMaker3(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments, double minJetPtCut, double trkbias, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StMyAnalysisMaker3();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();
    void    WriteTrackQAHistograms();
    void    WriteJetEPQAHistograms();

    // THnSparse Setup
    virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
    virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
    virtual THnSparse*      NewTHnSparseFCorr(const char* name, UInt_t entries);
    virtual void            GetDimParamsCorr(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetWriteTrackQAHistograms(Bool_t w){ doWriteTrackQAHist = w; }
    virtual void            SetWriteJetQAHistograms(Bool_t w)  { doWriteJetQAHist = w; }

    // jet setters
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetMaxTowerE(Double_t t)        { fTowerBias        = t; }    // tower bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }

    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

    // tower setters
    virtual void            SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
    virtual void            SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
    virtual void            SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

    // event mixing - setters
    virtual void            SetEventMixing(Int_t yesno)	       { fDoEventMixing=yesno; }
    virtual void            SetMixingTracks(Int_t tracks)      { fMixingTracks = tracks; }
    virtual void            SetNMixedTr(Int_t nmt)             { fNMIXtracks = nmt; }
    virtual void            SetNMixedEvt(Int_t nme)            { fNMIXevents = nme; }
    virtual void            SetCentBinSize(Int_t centbins)     { fCentBinSize = centbins; }
    virtual void            SetReduceStatsCent(Int_t red)      { fReduceStatsCent = red; }
    virtual void            SetDoFilterPtMixEvents(Int_t fil)  { fDoFilterPtMixEvents = fil; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }
    virtual void            SetMixedEventType(UInt_t me)       { fMixingEventType = me; }

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
    void                    SetOutFileNameQA(TString QAout)                 {mOutNameQA = QAout; }

    virtual void            SetEventPlaneMakerName(const char *epn)         {fEventPlaneMakerName = epn; }

  protected:
    TH1*                    FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    TH1*                    FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
    Double_t                GetReactionPlane();                                   // get reaction plane angle
    void                    GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin);// get event plane / flatten and fill histos 
    Bool_t                  AcceptJet(StJet *jet);           // jets accept cuts function
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1);
    static Double_t         CalculateEventPlaneChi(Double_t res);
    void                    TrackQA();

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
    Bool_t                  doPrintEventCounter;     // print event # switch
    Bool_t                  doWriteTrackQAHist;      // write track QA histograms
    Bool_t                  doWriteJetQAHist;        // write jet QA histograms
    Int_t                   fDoEffCorr;              // efficiency correction to tracks
    Bool_t                  doEventPlaneRes;         // event plane resolution switch
    Bool_t                  doTPCptassocBin;         // TPC event plane calculated on a pt assoc bin basis
    Int_t                   fTPCptAssocBin;          // pt associated bin to calculate event plane for

    // cuts
    //Double_t                fMinPtJet;               // min jet pt to keep jet in output
    //Double_t                fJetConstituentCut;      // min jet constituent

    // event mixing
    Int_t          fDoEventMixing;              // switch ON/off event mixing
    Int_t          fMixingTracks;               // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t          fNMIXtracks;                 // MIN # of mixing track in pool before performing mixing
    Int_t          fNMIXevents;                 // MIN # of mixing events in pool before performing mixing
    Int_t          fCentBinSize;                // centrality bin size of mixed event pools
    Int_t          fReduceStatsCent;            // bins to use for reduced statistics of sparse
    Bool_t         fDoFilterPtMixEvents;        // filter mixed event pool by pt (reduce memory) switch

    // event selection types
    UInt_t         fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t         fMBEventType;                // Physics selection of event used for MB
    UInt_t         fMixingEventType;            // Physics selection of event used for mixed event
    Int_t          fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    // used for event plane calculation and resolution
    //Float_t        fExcludeLeadingJetsFromFit;  // exclude n leading jets from fit
    //Int_t          fTrackWeight;                // track weight for Q-vector summation
    Double_t       fEventPlaneMaxTrackPtCut;    // max track pt cut for event plane calculation
    Int_t          fTPCEPmethod;                // TPC event plane calculation method
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
    Double_t       TPC_PSI2;
    Double_t       TPCA_PSI2;
    Double_t       TPCB_PSI2;
    Double_t       BBC_PSI2;
    Double_t       ZDC_PSI2;
    Double_t       BBC_PSI1;
    Double_t       ZDC_PSI1;
    Double_t       PSI2;
    Double_t       RES;
    // temp (possibly)
    Double_t       TPC_raw_comb;
    Double_t       TPC_raw_neg;
    Double_t       TPC_raw_pos;
    Double_t       BBC_raw_comb;
    Double_t       BBC_raw_east;
    Double_t       BBC_raw_west;
    Double_t       ZDC_raw_comb;
    Double_t       ZDC_raw_east;
    Double_t       ZDC_raw_west;

    // event pool
    TClonesArray          *CloneAndReduceTrackList();
    StEventPoolManager    *fPoolMgr;//!  // event pool Manager object

  private:
    Int_t                  fRunNumber;
    TString                fEPcalibFileName; 
    Double_t               fEPTPCResolution;
    Double_t               fEPTPCn;
    Double_t               fEPTPCp;
    Double_t               fEPTPC;
    Double_t               fEPBBC;
    Double_t               fEPZDC;

    // switches
    bool         doComments;

    // histograms
    TH1F* hEventPlane;//!   
    TH1F* hEventPlaneWeighted;//!
    TH2F* fHistEPTPCn;//!
    TH2F* fHistEPTPCp;//!
    TH2F* fHistEPBBC;//!
    TH2F* fHistEPZDC;//!
    TH1F* hEventZVertex;//!
    TH1F* hCentrality;//!
    TH1F* hMultiplicity;//!
    TH2F* hRhovsCent;//!
    TH1F* hTrackPhi[9];//!
    TH1F* hTrackEta[9];//!
    TH1F* hTrackPt[9];//!
    TH2F* hTrackEtavsPhi;//!

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
    TH1F* hJetTracksZ;//!
    TH2F* hJetPtvsArea;//!
    TH1F* hJetEventEP;//!
    TH2F* hJetPhivsEP;//!

    TH1F* hJetPtIn;//!
    TH1F* hJetPhiIn;//!
    TH1F* hJetEtaIn;//!
    TH1F* hJetEventEPIn;//!
    TH2F* hJetPhivsEPIn;//!
    TH1F* hJetPtMid;//!
    TH1F* hJetPhiMid;//!
    TH1F* hJetEtaMid;//!
    TH1F* hJetEventEPMid;//!
    TH2F* hJetPhivsEPMid;//!
    TH1F* hJetPtOut;//!
    TH1F* hJetPhiOut;//!
    TH1F* hJetEtaOut;//!
    TH1F* hJetEventEPOut;//!
    TH2F* hJetPhivsEPOut;//!

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

    // THn Sparse's jet sparse
    THnSparse             *fhnJH;//!           // jet hadron events matrix
    THnSparse             *fhnMixedEvents;//!  // mixed events matrix
    THnSparse             *fhnCorr;//!         // sparse to get # jet triggers

    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    ClassDef(StMyAnalysisMaker3, 1)
};
#endif
