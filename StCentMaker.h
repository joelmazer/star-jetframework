#ifndef StCentMaker_h
#define StCentMaker_h

#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class THnSparse;
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

    // THnSparse Setup
    virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
    virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

    // switches
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }

    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

    // event selection - setters
    virtual void            SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
    virtual void            SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }

    Int_t                   GetgrefMult() const                { return kgrefMult; }
    Int_t                   GetrefMult() const                 { return krefMult; }
    Int_t                   GetRef9() const                    { return kref9; }   // 9 bin increasing % with bin
    Int_t                   GetRef16() const                   { return kref16; }  // 16 bin increasing % with bin
    Int_t                   GetCent9() const                   { return kcent9; }  // 9 bin decreasing % with bin - not used by ME, but for reference
    Int_t                   GetCent16() const                  { return kcent16; } // 16 bin decreasing % with bin- not used by ME, but for reference
    Double_t                GetRefCorr2() const                { return krefCorr2; }
    Double_t                GetCentScaled() const              { return kCentralityScaled; }
    Double_t                GetMB5toMB30ReWeight() const       { return kMB5toMB30ReWeight; }
    Double_t                GetReWeight() const                { return kReWeight; }

  protected:
    // functions
    void                    SetSumw2(); // set errors weights 
    Double_t                GetCorrectedMultiplicity(const UShort_t RefMult, const Double_t z, const Double_t zdcCoincidenceRate, const UInt_t flag);
    Int_t                   GetZvtxBin(Double_t zvertex) const;

    // event selection types
    UInt_t                  fEmcTriggerEventType;        // Physics selection of event used for signal
    UInt_t                  fMBEventType;                // Physics selection of event used for MB

    // member variables
    Int_t                   kgrefMult;                   // global ref mult
    Int_t                   krefMult;                    // ref mult
    Double_t                krefCorr2;                   // ref mult corrected
    Int_t                   kref9;                       // 9 bin reference multiplicity (0-5, 5-10, 10-20, 20-30, .. 70-80)  
    Int_t                   kref16;                      // 16 bin reference multiplicity (0-5, 5-10, ... 70-75, 75-80)
    Int_t                   kcent9;                      // 9 bin reverse order (STAR has bad practice, hence the above added in increasing %)
    Int_t                   kcent16;                     // 16 bin reverse order (STAR has bad practice, hence the above added in increasing %)
    Double_t                kCentralityScaled;           // scaled centrality in 5% ranges
    Double_t                kMB5toMB30ReWeight;          // MB5 to MB30 conversion weight (includes re-weight)
    Double_t                kReWeight;                   // ReWeight for peripheral bins

    // centrality objects
    StRefMultCorr          *grefmultCorr;
    StRefMultCorr          *refmultCorr;
    StRefMultCorr          *refmult2Corr;
    StRefMultCorr          *grefmultCorrMB5;
    StRefMultCorr          *grefmultCorrUtil;

  private:
    // switches
    bool                    doComments;

    // histograms
    TH1F *hEventZVertex;//!
    TH1F *hCentrality;//!
    TH1F *hMultiplicity;//!
    TH1F *fHistEventSelectionQA;//! 

    TH1F *hMB5refCorr2ReWeight;//!
    TH1F *hMB5refCorr2Raw;//!
    TH1F *hMB5CentScaledReWeight;//!
    TH1F *hMB5CentScaledRaw;//!
    TH1F *hMB30refCorr2ReWeight;//!
    TH1F *hMB30refCorr2Raw;//!
    TH1F *hMB30CentScaledReWeight;//!
    TH1F *hMB30CentScaledRaw;//!
      
    TH1F *hMB5refCorr2[30];//!
    TH1F *hMB5refCorr2ReWt[30];//!
    TH1F *hMB5refCorr2ReWtInv[30];//!
    TH1F *hMB5grefMult[30];//!
    TH1F *hMB5CentScaled[30];//!
    TH1F *hMB5refCorr2Weight[30];//!
    TH1F *hMB5refCorr2WeightInv[30];//!

    TH1F *hMB30refCorr2[30];//!
    TH1F *hMB30refCorr2ReWt[30];//!
    TH1F *hMB30refCorr2ReWtInv[30];//!
    TH1F *hMB30grefMult[30];//!
    TH1F *hMB30CentScaled[30];//!

    TH1F *hHT2refCorr2[30];//!
    TH1F *hHT2refCorr2ReWt[30];//!
    TH1F *hHT2refCorr2ReWtInv[30];//!
    TH1F *hHT2grefMult[30];//!
    TH1F *hHT2CentScaled[30];//!

    TH1F *hMB5onlyrefCorr2Raw[30];//!
    TH1F *hMB30onlyrefCorr2Raw[30];//!
    TH1F *hMB5MB30refCorr2Raw[30];//!
    TH1F *hMB5onlyrefCorr2RawTotal;//!
    TH1F *hMB30onlyrefCorr2RawTotal;//!
    TH1F *hMB5MB30refCorr2RawTotal;//!

    // bad run list 
    std::set<Int_t>        badRuns;

    // base class pointer object
    StJetFrameworkPicoBase *mBaseMaker;

    // maker names
    TString                fAnalysisMakerName;

    // THn Sparse's jet sparse
    THnSparse             *fhnCentQA;//!           // centrality events matrix

    ClassDef(StCentMaker, 2)
};
#endif
