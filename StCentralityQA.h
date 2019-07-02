#ifndef StCentralityQA_h
#define StCentralityQA_h

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

class StCentralityQA : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality,
    };

    StCentralityQA(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments);
    virtual ~StCentralityQA();
   
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
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerE(Double_t mxE)    { fMaxEventTowerE = mxE; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }
    virtual void            SetBadRunListVers(Int_t i)         { fBadRunListVers = i; }

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

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Int_t effcorr)        { fDoEffCorr = effcorr; }

  protected:
    TH1                    *FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function
    void                    FillTowerTriggersArr();

    // FIXME - so far for Run14
    StRefMultCorr *grefmultCorr;
    StRefMultCorr *grefmultCorrNEW;

    // switches
    Int_t                   fDoEffCorr;                  // efficiency correction to tracks

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB
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

  private:
    Int_t                   fRunNumber;

    // switches
    bool                    doComments;

    // histograms
    TH1F *hEventZVertex;//!
    TH1F *hCentrality;//!
    TH1F *hCentralityNEW;//!
    TH1F *hCentralityDiff;//!
    TH1F *hMultiplicity;//!
    TH1F *hMultiplicityNEW;//!
    TH1F *hMultiplicityDiff;//!

    // QA histos
    TH1  *fHistEventSelectionQA;//! 
    TH1  *fHistEventSelectionQAafterCuts;//!
    TH1  *hTriggerIds;//!
    TH1  *hEmcTriggers;//!

    // bad run list 
    void                   ResetBadRunList( );
    Bool_t                 AddBadRuns(TString csvfile);
    Bool_t                 IsRunOK( Int_t mRunId );
    std::set<Int_t>        badRuns;

    // maker names
    TString                fAnalysisMakerName;

    ClassDef(StCentralityQA, 2)
};
#endif
