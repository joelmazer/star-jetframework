#ifndef StMyAnalysisMaker3_h
#define StMyAnalysisMaker3_h

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

// jet-framework classes
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;
class StEventPoolManager;
class StEventPool;
class StCentMaker;

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
      kDebugLeadSubLeadJets,
      kDebugTowersFiringTriggers
    };

    // enumerator for jet analysis jet type
    enum fJetAnalysisJetTypeEnum {
      kInclusiveJets,
      kLeadingJets,
      kSubLeadingJets
    };

    // enumerator for jet analysis jet type
    enum fSystematicUncTypeEnum {
      kDoNothing,
      kTrkEffMin,
      kTrkEffMax
    };

    StMyAnalysisMaker3(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments, double minJetPtCut, const char *jetMakerName, const char *rhoMakerName);
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
    virtual THnSparse      *NewTHnSparseF(const char *name, UInt_t entries);
    virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
    virtual THnSparse      *NewTHnSparseFCorr(const char *name, UInt_t entries);
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
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetWriteTrackQAHistograms(Bool_t w){ doWriteTrackQAHist = w; }
    virtual void            SetWriteJetQAHistograms(Bool_t w)  { doWriteJetQAHist = w; }
    virtual void            SetdoUseMainEPAngle(Bool_t m)      { doUseMainEPAngle = m; }
    virtual void            SetSystematicUncType(Int_t a)      { fSysUncType = a; }

    // jet setters
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetMaxTowerEt(Double_t t)       { fTowerBias        = t; }    // tower bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    virtual void            SetJetShapeTrackPtRange(Double_t min, Double_t max)  { fJetShapeTrackPtMin = min; fJetShapeTrackPtMax = max; }  // jet shape analysis pt range
    virtual void            SetJetLJSubLJPtThresholds(Double_t lj, Double_t slj) { fLeadJetPtMin = lj; fSubLeadJetPtMin = slj; }
    virtual void            SetdoSkip1ParticleJets(Bool_t sk)  { doSkip1ParticleJets = sk; } // skip 1 particle jets
    virtual void            SetdoBiasJetLeadConstituent(Bool_t j) { doBiasJetLeadConstituent = j; } // bias jet shape / jet-had jets
    virtual void            SetdoRequireJetTowFireTrig(Bool_t a) { doRequireJetTowFireTrig = a; }   // require jet tower to have fired HT trigger      

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt; } // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt; } // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d;         }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h;       }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r;     }

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
    virtual void            SetDoFilterPtMixEvents(Bool_t fil) { fDoFilterPtMixEvents = fil; }
    virtual void            SetDoUseMultBins(Bool_t mult)      { fDoUseMultBins = mult; }
    virtual void            SetdoUseEPBins(Bool_t ep)          { doUseEPBins = ep; }
    virtual void            SetnEPBins(Int_t nep)              { fnEPBins = nep; }
    virtual void            SetdoIgnoreExternalME(Bool_t ig)   { doIgnoreExternalME = ig; }
    virtual void            SetBGConeFractionCut(Double_t c)   { fBackgroundConeFractionCut = c; }
    virtual void            SetdoGenerateBadMixEventBGcone(Bool_t i) { doGenerateBadMixEventBGcone = i; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Bool_t effcorr)       { fDoEffCorr = effcorr; }
    virtual void            SetTrackEfficiencyType(Int_t t)    { fTrackEfficiencyType = t; }

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
    void                    SetExternalEventPoolManager(StEventPoolManager *mgr) { fPoolMgr = mgr;}
    StEventPoolManager     *GetEventPoolManager()                                { return fPoolMgr;}
    void                    SetUsePtBinnedEventPool(Bool_t val)                  { fUsePtBinnedEventPool = val;}
    void                    SetCheckEventNumberInMixedEvent(Bool_t val)          { fCheckEventNumberInMixedEvent = val;}

    // Set which pools will be saved
    virtual void            AddEventPoolsToOutput(Double_t minCent, Double_t maxCent, Double_t minZvtx, Double_t maxZvtx, Double_t minPsi2, Double_t maxPsi2, Double_t minPt, Double_t maxPt);

  protected:
    TH1                    *FillEmcTriggersHist(TH1 *h);                          // EmcTrigger counter histo
    Double_t                GetReactionPlane();                                   // get reaction plane angle
    void                    GetEventPlane(Bool_t flattenEP, Int_t n, Int_t method, Double_t ptcut, Int_t ptbin);// get event plane / flatten and fill histos 
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    CalculateEventPlaneResolution(Double_t bbc, Double_t zdc, Double_t tpc, Double_t tpcN, Double_t tpcP, Double_t bbc1, Double_t zdc1);
    static Double_t         CalculateEventPlaneChi(Double_t res);
    void                    TrackQA();
    void                    RunJetQA();
    void                    FillTowerTriggersArr();
    Bool_t                  DidTowerConstituentFireTrigger(StJet *jet);
    Bool_t                  DidBadTowerFireTrigger();
    Bool_t                  DidBadTowerFireHTTrigger(); // TEST - August 2019
    void                    JetHadronCorrelationAnalysis(StJet *jet, StEventPool *pool, Int_t centbin, Int_t assocPtBin);
    void                    JetShapeAnalysis(StJet *jet, StEventPool *pool, Double_t refCorr2, Int_t assocPtBin);
    void                    GetJetV2(StJet *jet, Double_t EPangle, Int_t ptAssocBin);
    void                    FillTriggerIDs(TH1 *h);
    void                    SetupMixEvtPool();
    Int_t                   GetZvtxBin(Double_t zvertex) const;

    // switches
    Bool_t                  doPrintEventCounter;     // print event # switch
    Int_t                   fJetAnalysisJetType;     // type of jets to use for jet analysis
    Bool_t                  doRequireAjSelection;    // requirement of Aj selection on jets for Jet Shape Analysis
    Bool_t                  doWriteTrackQAHist;      // write track QA histograms
    Bool_t                  doWriteJetQAHist;        // write jet QA histograms
    Bool_t                  fDoEffCorr;              // efficiency correction to tracks
    Int_t                   fTrackEfficiencyType;    // track efficiency type: pt-eta, pt, eta
    Bool_t                  doEventPlaneRes;         // event plane resolution switch
    Bool_t                  doTPCptassocBin;         // TPC event plane calculated on a pt assoc bin basis
    Int_t                   fTPCptAssocBin;          // pt associated bin to calculate event plane for
    Bool_t                  doUseMainEPAngle;        // use 0.2-2.0 GeV charged tracks for event plane
    Bool_t                  doIgnoreExternalME;      // does standared event mixing (without use of external approach)
    Bool_t                  doRunAnalysis;           // switch to run jet shape / jet-hadron correlation analyses
    Bool_t                  doJetHadronCorrelationAnalysis; // perform jet-hadron correlation analysis
    Bool_t                  doSkip1ParticleJets;     // switch to skip 1-particle jets
    Bool_t                  doBiasJetLeadConstituent;// switch to require leading bias to jet shape / jet-had jets
    Bool_t                  doRequireJetTowFireTrig; // switch to require jet constituent tower to have fired an event HT trigger  
    Int_t                   fSysUncType;             // systematic uncertainty type

    // cuts
    Double_t                fJetShapeTrackPtMin;     // jet shape analysis - min track pt
    Double_t                fJetShapeTrackPtMax;     // jet shape analysis - max track pt
    Double_t                fLeadJetPtMin;           // leading jet pt min
    Double_t                fSubLeadJetPtMin;        // sub-leading jet pt min

    // event mixing
    Int_t                   fDoEventMixing;          // switch ON/off event mixing
    Int_t                   fMixingTracks;           // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t                   fNMIXtracks;             // MIN # of mixing track in pool before performing mixing
    Int_t                   fNMIXevents;             // MIN # of mixing events in pool before performing mixing
    Int_t                   fCentBinSize;            // centrality bin size of mixed event pools
    Int_t                   fReduceStatsCent;        // bins to use for reduced statistics of sparse
    Bool_t                  fDoFilterPtMixEvents;    // filter mixed event pool by pt (reduce memory) switch
    Bool_t                  fDoUseMultBins;          // use multiplicity bins instead of centrality bins - used for Jet Shape Analysis
    Bool_t                  doUseEPBins;             // use event plane bins: 0.2-2.0 GeV charged tracks
    Int_t                   fnEPBins;                // number of event plane bins to use for event mixing (0, pi) range
    Double_t                fBackgroundConeFractionCut; // cut threshold to exclude mixed events, when pt sum of tracks in background cone reach 
    Bool_t                  doGenerateBadMixEventBGcone; // loops through mixed event prior to running analysis portion, to change normalization and eliminate bad events from BGcone

    // event selection types
    UInt_t                  fEmcTriggerEventType;    // Physics selection of event used for signal
    UInt_t                  fMBEventType;            // Physics selection of event used for MB
    Int_t                   fEmcTriggerArr[8];       // EMCal triggers array: used to select signal and do QA

    // tower to firing trigger type matched array
    Bool_t                  fTowerToTriggerTypeHT1[4800];// Tower with corresponding HT1 trigger type array
    Bool_t                  fTowerToTriggerTypeHT2[4800];// Tower with corresponding HT2 trigger type array
    Bool_t                  fTowerToTriggerTypeHT3[4800];// Tower with corresponding HT3 trigger type array

    // used for event plane calculation and resolution
    //Float_t                 fExcludeLeadingJetsFromFit;  // exclude n leading jets from fit
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
    
    // track efficiency file and function
    TFile                  *fEfficiencyInputFile;
    Double_t                ApplyTrackingEff(Bool_t applyEff, Double_t tpt, Double_t teta, Int_t cbin, Double_t ZDCx, Int_t effType, TFile *infile); // single-track reconstruction efficiency 
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
    TH1F *hCentralityPostCut;//!
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
    TH1F *hJetMass;//!
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
    TH1F *hJetHTrigMaxTowEt;//!
    TH1F *hJetHTrigMaxTrkPt;//!
    TH2 *fHistJetHEtaPhi;//!

    // QA histos
    TH1  *fHistEventSelectionQA;//! 
    TH1  *fHistEventSelectionQAafterCuts;//!
    TH1  *hEmcTriggers;//!
    TH1  *hEventTriggerIDs;//!
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
    TH1  *hNMixNormBefore[9];//!
    TH1  *hNMixNormAfter[9];//!
    TH1  *hBGconeFractionOfJetPt;//!
    TH2  *hJetPtvsBGconeFraction;//!
    TH2  *hJetPtvsBGconePt;//!

    TH1  *hMB5TrkPtRaw;//!
    TH1  *hMB5TrkPtReWeight;//!
    TH1  *hMB30TrkPtRaw;//!
    TH1  *hMB30TrkPtReWeight;//!

    TH1  *hNEventsvsZvtxMB5;//!
    TH1  *hNEventsvsCentMB5;//!
    TH1  *hNEventsvsMultMB5;//!
    TH2  *hNEventsvsZvsCentMB5;//!
    TH1  *hNEventsvsZvtxMB30;//!
    TH1  *hNEventsvsCentMB30;//!
    TH1  *hNEventsvsMultMB30;//!
    TH2  *hNEventsvsZvsCentMB30;//!
    TH1  *hNEventsvsZvtxMB5Wt;//!
    TH1  *hNEventsvsCentMB5Wt;//!
    TH1  *hNEventsvsMultMB5Wt;//!
    TH2  *hNEventsvsZvsCentMB5Wt;//!
    TH1  *hNEventsvsZvtxMB30Wt;//!
    TH1  *hNEventsvsCentMB30Wt;//!
    TH1  *hNEventsvsMultMB30Wt;//!
    TH2  *hNEventsvsZvsCentMB30Wt;//!
    TH1  *hNEventsvsZvtxHT2;//!
    TH1  *hNEventsvsCentHT2;//!
    TH1  *hNEventsvsMultHT2;//!
    TH2  *hNEventsvsZvsCentHT2;//!
    TH1  *hNPairsvsZvtxMB5;//!
    TH1  *hNPairsvsZvtxMB30;//!
    TH1  *hNPairsvsZvtxHT2;//!
    TH1  *hNPairsvsZvtxMB5Wt;//!
    TH1  *hNPairsvsZvtxMB30Wt;//!

    TH2F *hTPCvsBBCep;//!
    TH2F *hTPCvsZDCep;//!
    TH2F *hBBCvsZDCep;//!

    // jet shape histos - in jetpt and centrality arrays
    TH1F *hJetShape[4][4][4][9];//! jet shape histograms in annuli bins
    TH1F *hJetShapeCase1[4][4][4][9];//! jet shape case1 histograms in annuli bins
    TH1F *hJetShapeCase2[4][4][4][9];//! jet shape case2 histograms in annuli bins
    TH1F *hJetShapeCase3[4][4][4][9];//! jet shape case3 histograms in annuli bins
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
    TH1F *hJetPtProfileCase3[4][4][4][9];//! jet pt profile case3 histograms in annuli bins
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

    // base class pointer
    StJetFrameworkPicoBase *mBaseMaker;

    // bad and dead tower list functions and arrays
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list 
    std::set<Int_t>        badRuns;

    // Event pool variables - TEST
    vector<vector<Double_t> >   fEventPoolOutputList; // vector representing a list of pools (given by value range) that will be saved
    Bool_t                      fUsePtBinnedEventPool; // uses event pool in pt bins
    Bool_t                      fCheckEventNumberInMixedEvent; // check event number before correlation in mixed event
    TList                      *fListOfPools; //  Output list of containers

    ClassDef(StMyAnalysisMaker3, 3)
};
#endif
