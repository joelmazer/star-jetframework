#ifndef StEventPoolMaker_h
#define StEventPoolMaker_h

#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH2;
class TString;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// jet-framework classes
class StEventPoolManager;
class StEventPool;

class StEventPoolMaker : public StJetFrameworkPicoBase {
  public:

    StEventPoolMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments);
    virtual ~StEventPoolMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt) { fMaxEventTrackPt = mxpt; }

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
    virtual void            SetDoUseMultBins(Bool_t mult)      { fDoUseMultBins = mult; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }
    virtual void            SetMixedEventType(UInt_t me)       { fMixingEventType = me; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Int_t effcorr)        { fDoEffCorr = effcorr; }

  protected:
    TH1                    *FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    TH1                    *FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    FillTowerTriggersArr();
    Bool_t                  DidTowerConstituentFireTrigger(StJet *jet);

    // switches
    Bool_t                  doPrintEventCounter;         // print event # switch
    Int_t                   fDoEffCorr;                  // efficiency correction to tracks

    // event mixing
    Int_t                   fDoEventMixing;              // switch ON/off event mixing
    Int_t                   fMixingTracks;               // MAX # of mixing tracks to keep in pool, before removing old to add new
    Int_t                   fNMIXtracks;                 // MIN # of mixing track in pool before performing mixing
    Int_t                   fNMIXevents;                 // MIN # of mixing events in pool before performing mixing
    Int_t                   fCentBinSize;                // centrality bin size of mixed event pools
    Int_t                   fReduceStatsCent;            // bins to use for reduced statistics of sparse
    Bool_t                  fDoUseMultBins;              // use multiplicity bins instead of centrality bins

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
    UInt_t                  fMixingEventType;            // Physics selection of event used for mixed event
    Int_t                   fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA

    // tower to firing trigger type matched array
    Bool_t                  fTowerToTriggerTypeHT1[4801];// Tower with corresponding HT1 trigger type array
    Bool_t                  fTowerToTriggerTypeHT2[4801];// Tower with corresponding HT2 trigger type array
    Bool_t                  fTowerToTriggerTypeHT3[4801];// Tower with corresponding HT3 trigger type array

    // used for event plane calculation and resolution
    Int_t                   fHistCentBinMin;             // min centrality bin for histogram loop
    Int_t                   fHistCentBinMax;             // max centrality bin for histogram loop
    Int_t                   fHistZvertBinMin;            // min z-vertex bin for histogram loop
    Int_t                   fHistZvertBinMax;            // min z-vertex bin for histogram loop

    // event pool
    TClonesArray           *CloneAndReduceTrackList();
    StEventPoolManager     *fPoolMgr;//!  // event pool Manager object

  private:
    Int_t                   fRunNumber;

    // switches
    bool                    doComments;

    // histograms
    TH1F *hEventZVertex;//!
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
    TH2F *hTrackEtavsPhi;//!

    // QA histos
    TH1  *fHistEventSelectionQA;//! 
    TH1  *fHistEventSelectionQAafterCuts;//!
    TH1  *hTriggerIds;//!
    TH1  *hEmcTriggers;//!
    TH1  *hMixEvtStatZVtx;//!
    TH1  *hMixEvtStatCent;//!
    TH2  *hMixEvtStatZvsCent;//!

    // maker names
    TString                fAnalysisMakerName;
    TString                fEventMixerMakerName;

    ClassDef(StEventPoolMaker, 2)
};
#endif
