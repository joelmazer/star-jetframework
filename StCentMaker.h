#ifndef StCentMaker_h
#define StCentMaker_h

#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TString;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

class StCentMaker : public StJetFrameworkPicoBase {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum {
      kDebugNothing, // don't want lowest elements to be used
      kDebugEmcTrigger,
      kDebugGeneralEvt,
      kDebugCentrality,
    };

    StCentMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName, bool mDoComments);
    virtual ~StCentMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    // switches
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerE(Double_t mxE)    { fMaxEventTowerE = mxE; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }
    virtual void            SetBadRunListVers(Int_t i)         { fBadRunListVers = i; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }

    Int_t                   GetgrefMult() const                { return kgrefMult; }
    Int_t                   GetrefMult() const                 { return krefMult; }
    Int_t                   GetRef9() const                    { return kref9; }
    Int_t                   GetRef16() const                   { return kref16; }    
    Double_t                GetRefCorr2() const                { return krefCorr2; }
    Double_t                GetCentScaled() const              { return kCentralityScaled; }

  protected:
    // functions
    void                    SetSumw2(); // set errors weights 
    Double_t                GetCorrectedMultiplicity(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag);

    // RefMultCorr object
    StRefMultCorr *grefmultCorr;

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB

    // member variables
    Int_t                   kgrefMult;                   // global ref mult
    Int_t                   krefMult;                    // ref mult
    Double_t                krefCorr2;                   // ref mult corrected
    Int_t                   kref9;                       // 9 bin reference multiplicity (0-5, 5-10, 10-20, 20-30, .. 70-80)  
    Int_t                   kref16;                      // 16 bin reference multiplicity (0-5, 5-10, ... 70-75, 75-80)
    Double_t                kCentralityScaled;           // scaled centrality in 5% ranges

  private:
    // switches
    bool                    doComments;

    // histograms
    TH1F *hEventZVertex;//!
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
    TH1F *fHistEventSelectionQA;//! 

    // bad run list 
    void                   ResetBadRunList( );
    Bool_t                 AddBadRuns(TString csvfile);
    Bool_t                 IsRunOK( Int_t mRunId );
    std::set<Int_t>        badRuns;

    // maker names
    TString                fAnalysisMakerName;

    ClassDef(StCentMaker, 1)
};
#endif
