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
   
    // THnSparse Setup
    virtual THnSparse*      NewTHnSparseF(const char* name, UInt_t entries);
    virtual void            GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
    virtual THnSparse*      NewTHnSparseFCorr(const char* name, UInt_t entries);
    virtual void            GetDimParamsCorr(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

    static TString GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray* partCont, TClonesArray* clusCont, TString tag);

    // TClonesArrays function returners of analysis objects
    TClonesArray* jets() const { return mJets; }
    TClonesArray* tracks() const { return mTracks; }
    TClonesArray* towers() const { return mTowers; }
    TClonesArray* particles() const { return mParticles; }

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
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
    virtual void            SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
    virtual void            SetEventPlaneTrackWeight(Int_t weight)          {fTrackWeight = weight; }
    virtual void            SetEventPlaneMaxTrackPtCut(Double_t m)          {fEvenPlaneMaxTrackPtCut = m; }  

  protected:
    Int_t                  GetCentBin(Int_t cent, Int_t nBin) const;             // centrality bin
    Double_t               RelativePhi(Double_t mphi,Double_t vphi) const;       // relative jet track angle
    Double_t               RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    TH1*                   FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
    TH1*                   FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
    Double_t               GetReactionPlane();                                   // get reaction plane angle
    Bool_t                 AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function
    Bool_t                 AcceptJet(StJet *jet);           // jets accept cuts function
    Bool_t                 DoComparison(int myarr[], int elems);
    void                   SetSumw2(); // set errors 

    //Double_t               EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function

    // switches
    Bool_t                 doUsePrimTracks;         // primary track switch
    Int_t                  fDebugLevel;             // debug printout level
    Int_t                  fRunFlag;                // Run Flag numerator value
    Int_t                  fCentralityDef;          // Centrality Definition enumerator value
    Int_t                  fDoEffCorr;              // efficiency correction to tracks
    Bool_t                 fCorrJetPt;              // correct jet pt by rho

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

    // event selection types
    UInt_t         fTriggerEventType;           // Physics selection of event used for signal
    UInt_t         fMixingEventType;            // Physics selection of event used for mixed event
    Int_t          fEmcTriggerArr[7];           // EMCal triggers array: used to select signal and do QA

    // used for event plane calculation and resolution
    StJet*         fLeadingJet;//! leading jet
    Float_t        fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
    Int_t          fTrackWeight; // track weight for Q-vector summation
    Double_t       fEvenPlaneMaxTrackPtCut; // max track pt cut for event plane calculation

    // event pool
    TClonesArray         *CloneAndReduceTrackList();
    StEventPoolManager   *fPoolMgr;//!  // event pool Manager object

    // clonesarray collections of tracks and jets
    TClonesArray          *fTracksME;//! track collection to slim down for mixed events
    TClonesArray          *fJets;//! jet collection

  private:
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

/* 
    Int_t        mEventCounter;//!
    Int_t        mAllPVEventCounter;//!
    Int_t        mInputEventCounter;//!
*/

    // switches
    bool         doComments;

    // histograms
    TH1F* hEventPlane;//!   
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

    // THn Sparse's jet sparse
    THnSparse             *fhnJH;//!           // jet hadron events matrix
    THnSparse             *fhnMixedEvents;//!  // mixed events matrix
    THnSparse             *fhnCorr;//!         // sparse to get # jet triggers

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
