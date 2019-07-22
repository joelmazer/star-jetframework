#ifndef StMyAnalysisMaker3_h
#define StMyAnalysisMaker3_h

//#include "StMaker.h"
//#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StJetFrameworkPicoBase.h"
#include <set>

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

// jet-framework classes
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;
class StEventPoolManager;
class StEventPool;
class StCentMaker;

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
      kDebugJetConstituents,
      kDebugJetvsEPtype,
      kDebugRhoEstimate,
      kDebugLeadSubLeadJets
    };

    // enumerator for jet analysis jet type
    enum fJetAnalysisJetTypeEnum {
      kInclusiveJets,
      kLeadingJets,
      kSubLeadingJets
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
    void    WriteJetShapeHistograms(Int_t option);

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
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp; }
    virtual void            SetdoRunAnalysis(Bool_t ra)        { doRunAnalysis     = ra; }
    virtual void            SetdoJetHadronCorrelationAnalysis(Bool_t jhc) { doJetHadronCorrelationAnalysis = jhc; }
    virtual void            SetdoJetShapeAnalysis(Bool_t js)   { doJetShapeAnalysis = js; }
    virtual void            SetJetAnalysisJetType(Int_t t)     { fJetAnalysisJetType  = t; }
    virtual void            SetdoRequireAjSelection(Bool_t d)  { doRequireAjSelection = d; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetWriteTrackQAHistograms(Bool_t w){ doWriteTrackQAHist = w; }
    virtual void            SetWriteJetQAHistograms(Bool_t w)  { doWriteJetQAHist = w; }
    virtual void            SetdoUseMainEPAngle(Bool_t m)      { doUseMainEPAngle = m; }

    // jet setters
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetMaxTowerE(Double_t t)        { fTowerBias        = t; }    // tower bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    virtual void            SetJetShapeTrackPtRange(Double_t min, Double_t max)  { fJetShapeTrackPtMin = min; fJetShapeTrackPtMax = max; }  // jet shape analysis pt range
    virtual void            SetJetShapePtAssocBin(Int_t p)     { fJetShapePtAssocBin = p; }  // pt associated bin used in jet shape analysis 
    virtual void            SetJetLJSubLJPtThresholds(Double_t lj, Double_t slj) { fLeadJetPtMin = lj; fSubLeadJetPtMin = slj; }

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerE(Double_t mxE)    { fMaxEventTowerE = mxE; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

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
    virtual void            SetCentBinSizeJS(Int_t centbins)   { fCentBinSizeJS = centbins; }
    virtual void            SetReduceStatsCent(Int_t red)      { fReduceStatsCent = red; }
    virtual void            SetDoFilterPtMixEvents(Bool_t fil) { fDoFilterPtMixEvents = fil; }
    virtual void            SetDoUseMultBins(Bool_t mult)      { fDoUseMultBins = mult; }
    virtual void            SetdoUseEPBins(Bool_t ep)          { doUseEPBins = ep; }
    virtual void            SetdoIgnoreExternalME(Bool_t ig)   { doIgnoreExternalME = ig; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Int_t effcorr)        { fDoEffCorr = effcorr; }

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)        { fCorrJetPt = cpt; }

    // event plane
    virtual void            SetExcludeLeadingJetsFromFit(Float_t n)         { fExcludeLeadingJetsFromFit = n; }
    virtual void            SetEventPlaneTrackWeight(Int_t weight)          { fTrackWeight = weight; }
    virtual void            SetEventPlaneMaxTrackPtCut(Double_t m)          { fEventPlaneMaxTrackPtCut = m; }  
    virtual void            SetTPCEventPlaneMethod(Int_t tm)                { fTPCEPmethod = tm; }
    virtual void            SetHistBinLimitsCenZvert(Int_t cmin, Int_t cmax, Int_t zmin, Int_t zmax)   { fHistCentBinMin = cmin; fHistCentBinMax = cmax; fHistZvertBinMin = zmin; fHistZvertBinMax = zmax; }
    virtual void            SetdoEventPlaneRes(Bool_t depr)                 { doEventPlaneRes = depr; }
    virtual void            SetdoEPTPCptAssocMethod(Bool_t ptbin)           { doTPCptassocBin = ptbin; }
    virtual void            SetEPTPCptAssocBin(Int_t pb)                    { fTPCptAssocBin = pb; }

    // Where to read calib object with EP calibration if not default
    void                    SetEPcalibFileName(TString filename)            { fEPcalibFileName = filename; } 
    void                    SetOutFileNameEP(TString epout)                 { mOutNameEP = epout; }
    void                    SetOutFileNameQA(TString QAout)                 { mOutNameQA = QAout; }
    void                    SetOutFileNameMixEvt(TString MEout)             { mOutNameME = MEout; }

    virtual void            SetEventPlaneMakerName(const char *epn)         { fEventPlaneMakerName = epn; }

    // ##### External event pool configuration
    void                    SetExternalEventPoolManager(StEventPoolManager* mgr) { fPoolMgr = mgr;}
    StEventPoolManager*     GetEventPoolManager()                                { return fPoolMgr;}
    void                    SetUsePtBinnedEventPool(Bool_t val)                  { fUsePtBinnedEventPool = val;}
    void                    SetCheckEventNumberInMixedEvent(Bool_t val)          { fCheckEventNumberInMixedEvent = val;}

    // Set which pools will be saved
    virtual void            AddEventPoolsToOutput(Double_t minCent, Double_t maxCent, Double_t minZvtx, Double_t maxZvtx, Double_t minPsi2, Double_t maxPsi2, Double_t minPt, Double_t maxPt);

  protected:
    TH1                    *FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    Double_t                GetReactionPlane();                                   // get reaction plane angle
    void                    GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin);// get event plane / flatten and fill histos 
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1);
    static Double_t         CalculateEventPlaneChi(Double_t res);
    void                    TrackQA();
    void                    FillTowerTriggersArr();
    Bool_t                  DidTowerConstituentFireTrigger(StJet *jet);
    Bool_t                  DidBadTowerFireTrigger();
    Int_t                   JetHadronCorrelationAnalysis(StJet *jet, StEventPool *pool, Int_t centbin, Int_t assocPtBin);
    Int_t                   JetShapeAnalysis(StJet *jet, StEventPool *pool, Double_t refCorr2, Int_t assocPtBin);
    void                    GetJetV2(StJet *jet, Double_t EPangle, Int_t ptAssocBin);

    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Int_t                   fJetAnalysisJetType;     // type of jets to use for jet analysis
    Bool_t                  doRequireAjSelection;    // requirement of Aj selection on jets for Jet Shape Analysis
    Bool_t                  doWriteTrackQAHist;      // write track QA histograms
    Bool_t                  doWriteJetQAHist;        // write jet QA histograms
    Int_t                   fDoEffCorr;              // efficiency correction to tracks
    Bool_t                  doEventPlaneRes;         // event plane resolution switch
    Bool_t                  doTPCptassocBin;         // TPC event plane calculated on a pt assoc bin basis
    Int_t                   fTPCptAssocBin;          // pt associated bin to calculate event plane for
    Bool_t                  doUseMainEPAngle;        // use 0.2-2.0 GeV charged tracks for event plane
    Bool_t                  doRejectBadRuns;         // switch to reject bad runs and thus skip from analysis
    Bool_t                  doIgnoreExternalME;      // does standared event mixing (without use of external approach)
    Bool_t                  doRunAnalysis;           // switch to run jet shape / jet-hadron correlation analyses
    Bool_t                  doJetHadronCorrelationAnalysis; // perform jet-hadron correlation analysis

    // cuts
    //Double_t                fMinPtJet;               // min jet pt to keep jet in output
    //Double_t                fJetConstituentCut;      // min jet constituent
    Double_t                fJetShapeTrackPtMin;     // jet shape analysis - min track pt
    Double_t                fJetShapeTrackPtMax;     // jet shape analysis - max track pt
    Int_t                   fJetShapePtAssocBin;     // jet shape analysis - pt associated bin
    Double_t                fLeadJetPtMin;           // leading jet pt min
    Double_t                fSubLeadJetPtMin;        // sub-leading jet pt min

    // event mixing
    Int_t                   fDoEventMixing;          // switch ON/off event mixing
    Int_t                   fMixingTracks;           // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t                   fNMIXtracks;             // MIN # of mixing track in pool before performing mixing
    Int_t                   fNMIXevents;             // MIN # of mixing events in pool before performing mixing
    Int_t                   fCentBinSize;            // centrality bin size of mixed event pools
    Int_t                   fCentBinSizeJS;          // centrality bin size of mixed event pools for jet shape analysis
    Int_t                   fReduceStatsCent;        // bins to use for reduced statistics of sparse
    Bool_t                  fDoFilterPtMixEvents;    // filter mixed event pool by pt (reduce memory) switch
    Bool_t                  fDoUseMultBins;          // use multiplicity bins instead of centrality bins - used for Jet Shape Analysis
    Bool_t                  doUseEPBins;             // use event plane bins: 0.2-2.0 GeV charged tracks

    // event selection types
    UInt_t                  fEmcTriggerEventType;    // Physics selection of event used for signal
    UInt_t                  fMBEventType;            // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];       // EMCal triggers array: used to select signal and do QA

    // tower to firing trigger type matched array
    Bool_t                  fTowerToTriggerTypeHT1[4801];// Tower with corresponding HT1 trigger type array
    Bool_t                  fTowerToTriggerTypeHT2[4801];// Tower with corresponding HT2 trigger type array
    Bool_t                  fTowerToTriggerTypeHT3[4801];// Tower with corresponding HT3 trigger type array

    // used for event plane calculation and resolution
    //Float_t                 fExcludeLeadingJetsFromFit;  // exclude n leading jets from fit
    //Int_t                   fTrackWeight;                // track weight for Q-vector summation
    Double_t                fEventPlaneMaxTrackPtCut;// max track pt cut for event plane calculation
    Int_t                   fTPCEPmethod;            // TPC event plane calculation method
    Int_t                   fHistCentBinMin;         // min centrality bin for histogram loop
    Int_t                   fHistCentBinMax;         // max centrality bin for histogram loop
    Int_t                   fHistZvertBinMin;        // min z-vertex bin for histogram loop
    Int_t                   fHistZvertBinMax;        // min z-vertex bin for histogram loop

    // global variables used with TPC event plane corrections
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

    // event pool
    TClonesArray           *CloneAndReduceTrackList();
    StEventPoolManager     *fPoolMgr;//!  // event pool Manager object

  private:
    Int_t                   fRunNumber;
    TString                 fEPcalibFileName; 
    TString                 mOutNameME;
    Double_t                fEPTPCResolution;
    Double_t                fEPTPCn;
    Double_t                fEPTPCp;
    Double_t                fEPTPC;
    Double_t                fEPTPCcomb;
    Double_t                fEPBBC;
    Double_t                fEPZDC;

    // switches
    bool                    doComments;

    // histograms
    TH1F *hdEPReactionPlaneFnc;//!
    TH1F *hdEPEventPlaneFncN2;//!
    TH1F *hdEPEventPlaneFncP2;//!
    TH1F *hdEPEventPlaneFnc2;//!
    TH1F *hdEPEventPlaneClass;//!
    TH1F *hReactionPlaneFnc;//!
    TH1F *hEventPlaneFncN2;//!
    TH1F *hEventPlaneFncP2;//!
    TH1F *hEventPlaneFnc2;//!
    TH1F *hEventPlaneClass;//!

    TH1F *hEventPlane;//!   
    TH2F *fHistEPTPCn;//!
    TH2F *fHistEPTPCp;//!
    TH2F *fHistEPBBC;//!
    TH2F *fHistEPZDC;//!
    TH1F *hEventZVertex;//!
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
    TH1F *hStats;//!
    TH2F *hRhovsCent;//!
    TH1F *hdEPtrk[5];//!
    TH1F *hTrackPhi[9];//!
    TH1F *hTrackEta[9];//!
    TH1F *hTrackPt[9];//!
    TH2F *hTrackEtavsPhi;//!

    // jet histos
    TH1F *hJetPt;//!
    TH1F *hJetCorrPt;//!
    TH1F *hJetLeadingPt;//!
    TH1F *hJetSubLeadingPt;//!
    TH1F *hJetLeadingPtAj;//!
    TH1F *hJetSubLeadingPtAj;//!
    TH1F *hJetDiJetAj;//!
    TH1F *hJetE;//!
    TH1F *hJetEta;//!
    TH1F *hJetPhi;//!
    TH1F *hJetNEF;//!
    TH1F *hJetArea;//!
    TH1F *hJetTracksPt;//!
    TH1F *hJetTracksPhi;//!
    TH1F *hJetTracksEta;//!
    TH1F *hJetTracksZ;//!
    TH2F *hJetPtvsArea;//!

    TH1F *hJetEventEP[5];//!
    TH2F *hJetPhivsEP[5];//!
    TH1F *hJetPtIn[5];//!
    TH1F *hJetPhiIn[5];//!
    TH1F *hJetEtaIn[5];//!
    TH1F *hJetEventEPIn[5];//!
    TH2F *hJetPhivsEPIn[5];//!
    TH1F *hJetPtMid[5];//!
    TH1F *hJetPhiMid[5];//!
    TH1F *hJetEtaMid[5];//!
    TH1F *hJetEventEPMid[5];//!
    TH2F *hJetPhivsEPMid[5];//!
    TH1F *hJetPtOut[5];//!
    TH1F *hJetPhiOut[5];//!
    TH1F *hJetEtaOut[5];//!
    TH1F *hJetEventEPOut[5];//!
    TH2F *hJetPhivsEPOut[5];//!

    // correlation histo
    TH2  *fHistJetHEtaPhi;//!

    // QA histos
    TH1  *fHistEventSelectionQA;//! 
    TH1  *fHistEventSelectionQAafterCuts;//!
    TH1  *hTriggerIds;//!
    TH1  *hEmcTriggers;//!
    TH1  *hBadTowerFiredTrigger;//!
    TH1  *hNGoodTowersFiringTrigger;//!
    TH1  *hMixEvtStatZVtx;//!
    TH1  *hMixEvtStatCent;//!
    TH2  *hMixEvtStatZvsCent;//!
    TH1  *hTriggerEvtStatZVtx;//!
    TH1  *hTriggerEvtStatCent;//!
    TH2  *hTriggerEvtStatZvsCent;//!
    TH1  *hMBvsMult;//!
    TH1  *hMB5vsMult;//!
    TH1  *hMB30vsMult;//!
    TH1  *hHTvsMult;//!
    TH1  *hNMixEvents;//!

    TH2F *hTPCvsBBCep;//!
    TH2F *hTPCvsZDCep;//!
    TH2F *hBBCvsZDCep;//!

    // jet shape histos - in jetpt and centrality arrays
    TH1F *hJetShape[4][4][4][9];//! jet shape histograms in annuli bins
    TH1F *hJetShapeCase1[4][4][4][9];//! jet shape case1 histograms in annuli bins
    TH1F *hJetShapeCase2[4][4][4][9];//! jet shape case2 histograms in annuli bins
    TH1F *hJetShapeBG[4][4][4][9];//! jet shape backround histograms in annuli bins
    TH1F *hJetShapeBGCase1[4][4][4][9];//! jet shape case1 histograms in annuli bins
    TH1F *hJetShapeBGCase2[4][4][4][9];//! jet shape case2 histograms in annuli bins
    TH1F *hJetShapeBGCase3[4][4][4][9];//! jet shape case3 histograms in annuli bins
    TH1F *hJetCounter[4][4][4];//! jet shape and pt profile histogram - jet counter
    TH1F *hJetCounterCase1[4][4][4];//! jet shape and pt profile case1 histogram - jet counter
    TH1F *hJetCounterCase2[4][4][4];//! jet shape and pt profile case2 histogram - jet counter
    TH1F *hJetCounterCase3BG[4][4][4];//! jet shape and pt profile case3 histograms - jet counter FOR BG only!

    // jet pt profile histos - in jetpt and centrality arrays
    TH1F *hJetPtProfile[4][4][4][9];//! jet pt profile histograms in annuli bins
    TH1F *hJetPtProfileCase1[4][4][4][9];//! jet pt profile case1 histograms in annuli bins
    TH1F *hJetPtProfileCase2[4][4][4][9];//! jet pt profile case2 histograms in annuli bins
    TH1F *hJetPtProfileBG[4][4][4][9];//! jet pt profile backround histograms in annuli bins
    TH1F *hJetPtProfileBGCase1[4][4][4][9];//! jet profile case1 histograms in annuli bins
    TH1F *hJetPtProfileBGCase2[4][4][4][9];//! jet profile case2 histograms in annuli bins
    TH1F *hJetPtProfileBGCase3[4][4][4][9];//! jet profile case3 histograms in annuli bins

    // EP resoltuion profiles
    TProfile              *fProfV2Resolution[9];//! resolution parameters for v2
    TProfile              *fProfV3Resolution[9];//! resolution parameters for v3
    TProfile              *fProfV4Resolution[9];//! resolution parameters for v4
    TProfile              *fProfV5Resolution[9];//! resolution parameters for v5

    // jet vn measurement
    TProfile              *fProfJetV2[4][4][4];//! jet v2 

    // THn Sparse's jet sparse
    THnSparse             *fhnJH;//!           // jet hadron events matrix
    THnSparse             *fhnMixedEvents;//!  // mixed events matrix
    THnSparse             *fhnCorr;//!         // sparse to get # jet triggers

    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    // centrality maker pointer/
    StCentMaker            *mCentMaker; 

    // base class pointer
    StJetFrameworkPicoBase *mBaseMaker;

    // bad and dead tower list functions and arrays
    //void                   ResetBadTowerList( );
    //void                   ResetDeadTowerList( );
    //Bool_t                 AddBadTowers(TString csvfile);
    //Bool_t                 AddDeadTowers(TString csvfile);
    //Bool_t                 IsTowerOK( Int_t mTowId );
    //Bool_t                 IsTowerDead( Int_t mTowId );
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list 
    //void                   ResetBadRunList( );
    //Bool_t                 AddBadRuns(TString csvfile);
    //Bool_t                 IsRunOK( Int_t mRunId );
    std::set<Int_t>        badRuns;

    // Event pool variables - TEST
    vector<vector<Double_t> >   fEventPoolOutputList; // vector representing a list of pools (given by value range) that will be saved
    Bool_t                      fUsePtBinnedEventPool; // uses event pool in pt bins
    Bool_t                      fCheckEventNumberInMixedEvent; // check event number before correlation in mixed event
    TList                      *fListOfPools; //  Output list of containers

    ClassDef(StMyAnalysisMaker3, 2)
};
#endif
