#ifndef StJetShapeAnalysis_h
#define StJetShapeAnalysis_h

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

// jet-framework classes
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;
class StEventPoolManager;
class StEventPool;

//class StJetShapeAnalysis : public StMaker {
class StJetShapeAnalysis : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugMixedEvents,
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality,
      kDebugEventPlaneCalc,
      kDebugJetvsEPtype, 
      kDebugRhoEstimate
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

    // enumerator for TPC event plane method
    enum fJetShapeJetTypeEnum {
      kInclusiveJets,
      kLeadingJets,
      kSubLeadingJets
    };

    StJetShapeAnalysis(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments, double minJetPtCut, double trkbias, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StJetShapeAnalysis();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();
    void    WriteTrackQAHistograms();
    void    WriteJetShapeHistograms();

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp; }
    virtual void            SetdoJetShapeAnalysis(Bool_t js)   { doJetShapeAnalysis = js; }
    virtual void            SetJetShapeJetType(Int_t t)        { fJetShapeJetType  = t; }
    virtual void            SetdoRequireAjSelection(Bool_t d)  { doRequireAjSelection = d; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }
    virtual void            SetWriteTrackQAHistograms(Bool_t w){ doWriteTrackQAHist = w; }

    // jet setters
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetMaxTowerE(Double_t t)        { fTowerBias        = t; }    // tower bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    virtual void            SetJetShapeTrackPtRange(Double_t min, Double_t max)  { fJetShapeTrackPtMin = min; fJetShapeTrackPtMax = max; }  // jet shape analysis pt range
    virtual void            SetJetShapePtAssocBin(Int_t p)     { fJetShapePtAssocBin = p; }  // pt associated bin used in jet shape analysis 

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

    // track setters
    virtual void            SetMinarackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
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
    virtual void            SetCentBinSizeJS(Int_t centbins)   { fCentBinSizeJS = centbins; }
    virtual void            SetReduceStatsCent(Int_t red)      { fReduceStatsCent = red; }
    virtual void            SetDoFilterPtMixEvents(Int_t fil)  { fDoFilterPtMixEvents = fil; }
    virtual void            SetDoUseMultBins(Bool_t mult)      { fDoUseMultBins = mult; }

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
    virtual void            SetHistBinLimitsCenZvert(Int_t cmin, Int_t cmax, Int_t zmin, Int_t zmax)   { fHistCentBinMin = cmin; fHistCentBinMax = cmax; fHistZvertBinMin = zmin; fHistZvertBinMax = zmax; }
    virtual void            SetdoEPTPCptAssocMethod(Bool_t ptbin)           {doTPCptassocBin = ptbin; }
    virtual void            SetEPTPCptAssocBin(Int_t pb)                    {fTPCptAssocBin = pb; }

    // Where to read calib object with EP calibration if not default
    void                    SetEPcalibFileName(TString filename)            {fEPcalibFileName = filename; } 
    void                    SetOutFileNameEP(TString epout)                 {mOutNameEP = epout; }
    void                    SetOutFileNameQA(TString QAout)                 {mOutNameQA = QAout; }

    virtual void            SetEventPlaneMakerName(const char *epn)         {fEventPlaneMakerName = epn; }

  protected:
    TH1                    *FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    Double_t                GetReactionPlane();                                   // get reaction plane angle
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    TrackQA();
    void                    FillTowerTriggersArr();
    Bool_t                  DidTowerConstituentFireTrigger(StJet *jet);
    Int_t                   JetShapeAnalysis(StJet *jet, StEventPool *pool, Double_t refCorr2);

    // switches
    Int_t                   fJetShapeJetType;        // type of jets to use for jet shape analysis
    Bool_t                  doRequireAjSelection;    // requirement of Aj selection on jets for Jet Shape Analysis
    Bool_t                  doWriteTrackQAHist;      // write track QA histograms
    Int_t                   fDoEffCorr;              // efficiency correction to tracks
    Bool_t                  doTPCptassocBin;         // TPC event plane calculated on a pt assoc bin basis
    Int_t                   fTPCptAssocBin;          // pt associated bin to calculate event plane for
    Bool_t                  doRejectBadRuns;         // switch to reject bad runs and thus skip from analysis

    // cuts
    //Double_t                fMinPtJet;               // min jet pt to keep jet in output
    //Double_t                fJetConstituentCut;      // min jet constituent
    Double_t                fJetShapeTrackPtMin;     // jet shape analysis - min track pt
    Double_t                fJetShapeTrackPtMax;     // jet shape analysis - max track pt
    Int_t                   fJetShapePtAssocBin;     // jet shape analysis - pt associated bin

    // event mixing
    Int_t                   fDoEventMixing;          // switch ON/off event mixing
    Int_t                   fMixingTracks;           // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t                   fNMIXtracks;             // MIN # of mixing track in pool before performing mixing
    Int_t                   fNMIXevents;             // MIN # of mixing events in pool before performing mixing
    Int_t                   fCentBinSizeJS;          // centrality bin size of mixed event pools
    Int_t                   fReduceStatsCent;        // bins to use for reduced statistics of sparse
    Bool_t                  fDoFilterPtMixEvents;    // filter mixed event pool by pt (reduce memory) switch
    Bool_t                  fDoUseMultBins;          // use multiplicity bins instead of centrality bins - used for Jet Shape Analysis

    // event selection types
    UInt_t                  fEmcTriggerEventType;    // Physics selection of event used for signal
    UInt_t                  fMBEventType;            // Physics selection of event used for MB
    UInt_t                  fMixingEventType;        // Physics selection of event used for mixed event
    Int_t                   fEmcTriggerArr[8];       // EMCal triggers array: used to select signal and do QA

    // tower to firing trigger type matched array
    Bool_t                  fTowerToTriggerTypeHT1[4801];// Tower with corresponding HT1 trigger type array
    Bool_t                  fTowerToTriggerTypeHT2[4801];// Tower with corresponding HT2 trigger type array
    Bool_t                  fTowerToTriggerTypeHT3[4801];// Tower with corresponding HT3 trigger type array

    // used for event plane calculation and resolution
    Double_t                fEventPlaneMaxTrackPtCut;// max track pt cut for event plane calculation
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

    // event pool
    TClonesArray           *CloneAndReduceTrackList();
    StEventPoolManager     *fPoolMgr;//!  // event pool Manager object

  private:
    Int_t                   fRunNumber;
    TString                 fEPcalibFileName; 
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
    TH1F *hEventPlane;//!   
    TH1F *hEventZVertex;//!
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
    TH2F *hRhovsCent;//!
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
    TH1F *hJetEventEP;//!
    TH2F *hJetPhivsEP;//!

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
    TH1  *hTriggerEvtStatZVtx;//!
    TH1  *hTriggerEvtStatCent;//!
    TH2  *hTriggerEvtStatZvsCent;//!
    TH1  *hMBvsMult;//!
    TH1  *hMB5vsMult;//!
    TH1  *hMB30vsMult;//!
    TH1  *hHTvsMult;//!
    TH1  *hNMixEvents;//!

    // jet shape histos - in jetpt and centrality arrays
    TH1F *hJetShape[4][4][4];//! jet shape histograms in annuli bins
    TH1F *hJetShapeCase1[4][4][4];//! jet shape case1 histograms in annuli bins
    TH1F *hJetShapeCase2[4][4][4];//! jet shape case2 histograms in annuli bins
    TH1F *hJetShapeBG[4][4][4];//! jet shape backround histograms in annuli bins
    TH1F *hJetShapeBGCase1[4][4][4];//! jet shape case1 histograms in annuli bins
    TH1F *hJetShapeBGCase2[4][4][4];//! jet shape case2 histograms in annuli bins
    TH1F *hJetShapeBGCase3[4][4][4];//! jet shape case3 histograms in annuli bins
    TH1F *hJetCounter[4][4][4];//! jet shape and pt profile histogram - jet counter
    TH1F *hJetCounterCase1[4][4][4];//! jet shape and pt profile case1 histogram - jet counter
    TH1F *hJetCounterCase2[4][4][4];//! jet shape and pt profile case2 histogram - jet counter
    TH1F *hJetCounterCase3BG[4][4][4];//! jet shape and pt profile case3 histograms - jet counter FOR BG only!

    // jet pt profile histos - in jetpt and centrality arrays
    TH1F *hJetPtProfile[4][4][4];//! jet pt profile histograms in annuli bins
    TH1F *hJetPtProfileCase1[4][4][4];//! jet pt profile case1 histograms in annuli bins
    TH1F *hJetPtProfileCase2[4][4][4];//! jet pt profile case2 histograms in annuli bins
    TH1F *hJetPtProfileBG[4][4][4];//! jet pt profile backround histograms in annuli bins
    TH1F *hJetPtProfileBGCase1[4][4][4];//! jet profile case1 histograms in annuli bins
    TH1F *hJetPtProfileBGCase2[4][4][4];//! jet profile case2 histograms in annuli bins
    TH1F *hJetPtProfileBGCase3[4][4][4];//! jet profile case3 histograms in annuli bins

    // bad and dead tower list functions and arrays
    void                   ResetBadTowerList( );
    void                   ResetDeadTowerList( );
    Bool_t                 AddBadTowers(TString csvfile);
    Bool_t                 AddDeadTowers(TString csvfile);
    Bool_t                 IsTowerOK( Int_t mTowId );
    Bool_t                 IsTowerDead( Int_t mTowId );
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list 
    void                   ResetBadRunList( );
    Bool_t                 AddBadRuns(TString csvfile);
    Bool_t                 IsRunOK( Int_t mRunId );
    std::set<Int_t>        badRuns;

    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    ClassDef(StJetShapeAnalysis, 2)
};
#endif
