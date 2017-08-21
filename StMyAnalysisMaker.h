#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

// some includes
#include "StMaker.h"

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
class StRho;
class StRhoParameter;
class StEventPoolManager;

class StMyAnalysisMaker : public StMaker {
  public:

  // jet type enumerator
  enum EJetType_t {
    kFullJet,
    kChargedJet,
    kNeutralJet
  };

  // jet algorithm enumerator
  enum EJetAlgo_t {
    kt_algorithm                    = 0,
    antikt_algorithm                = 1,
    cambridge_algorithm             = 2,
    genkt_algorithm                 = 3,
    cambridge_for_passive_algorithm = 11,
    genkt_for_passive_algorithm     = 13,
    plugin_algorithm                = 99,
    undefined_jet_algorithm         = 999
  };

  // jet recombination schme enumerator
  enum ERecoScheme_t {
    E_scheme        = 0,
    pt_scheme       = 1,
    pt2_scheme      = 2,
    Et_scheme       = 3,
    Et2_scheme      = 4,
    BIpt_scheme     = 5,
    BIpt2_scheme    = 6,
    WTA_pt_scheme   = 7,
    WTA_modp_scheme = 8,
    external_scheme = 99
  };

  /**
   * @enum JetAcceptanceType
   * @brief Bit definition for jet geometry acceptance. Defined here for backwards compatibility. This will be 
   * removed. Please use StJet::JetAcceptanceType in your code.
   */
/*
  enum JetAcceptanceType {
    kTPC              = StJet::kTPC,          ///< TPC acceptance
    kTPCfid           = StJet::kTPCfid,       ///< TPC fiducial acceptance (each eta edge narrowed by jet R)
    kEMCAL            = StJet::kEMCAL,        ///< EMCal acceptance
    kEMCALfid         = StJet::kEMCALfid,     ///< EMCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kDCAL             = StJet::kDCAL,         ///< DCal acceptance -- spans entire rectangular region in eta-phi (including most of PHOS)
    kDCALfid          = StJet::kDCALfid,      ///< DCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kDCALonly         = StJet::kDCALonly,     ///< DCal acceptance -- spans ONLY DCal (no PHOS or gap)
    kDCALonlyfid      = StJet::kDCALonlyfid,  ///< DCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kPHOS             = StJet::kPHOS,         ///< PHOS acceptance
    kPHOSfid          = StJet::kPHOSfid,      ///< PHOS fiducial acceptance (each eta, phi edge narrowed by jet R)
    kUser             = StJet::kUser          ///< Full acceptance, i.e. no acceptance cut applied -- left to user
  };
*/

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
    void                    SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }

    // jet setters
    void                    SetMinJetPt(Double_t j)            { fMinPtJet         = j; }
    void                    SetMinJetTrackPt(Double_t t)       { fTrackBias        = t; }

    void                    SetMinTrackPt(Double_t tp)         { fTrackPtCut       = tp;}

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

    // use local rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    // don't use this: OLD from run11
    Int_t centrality(int);

  protected:
    Int_t                  GetCentBin(Int_t cent, Int_t nBin) const; // centrality bin
    Double_t               RelativePhi(Double_t mphi,Double_t vphi) const; // relative jet track angle
    Double_t               RelativeEPJET(Double_t jetAng, Double_t EPAng) const;  // relative jet event plane angle
    TH1*                   FillEventTriggerQA(TH1* h, UInt_t t); // filled event trigger QA plots

    //Double_t               EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function

    // switches
    Bool_t                 doUsePrimTracks;         // primary track switch
    Int_t                  fDoEffCorr; // efficiency correction to tracks
    Bool_t                 fCorrJetPt; // correct jet pt by rho

    // cuts
    Double_t               fMinPtJet;               // min jet pt to keep jet in output
    Double_t               fTrackBias;              // high pt track in jet bias

    Double_t               fTrackPtCut;             // min track pt cut

    Int_t      mCentrality;

    // event mixing
    Int_t          fDoEventMixing;
    Int_t          fMixingTracks;
    Int_t          fNMIXtracks;
    Int_t          fNMIXevents;
    Int_t          fCentBinSize; // centrality bin size of mixed event pools
    Int_t          fReduceStatsCent; // bins to use for reduced statistics of sparse

    // event selection types
    UInt_t         fTriggerEventType;
    UInt_t         fMixingEventType;

    // event pool
//    TObjArray      *CloneAndReduceTrackList(TObjArray* tracks);
    TClonesArray      *CloneAndReduceTrackList(TClonesArray* tracks);
    StEventPoolManager   *fPoolMgr;//!  // event pool Manager object

    // clonesarray collections of tracks and jets
    //TClonesArray          *tracksClone;//! mixed event track collection
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
    StRefMultCorr* refmultCorr;
    StRefMultCorr* refmult2Corr;

    // TCloneArray of analysis objects
    TClonesArray   *mJets;
    TClonesArray   *mTracks;
    TClonesArray   *mTowers;
    TClonesArray   *mParticles;
    
    TString      mOutName;
 
    // counters 
    Int_t        mEventCounter;//!
    Int_t        mAllPVEventCounter;//!
    Int_t        mInputEventCounter;//!
 
    // switches
    bool       doComments;

    // histograms
    TH1F*      mKsM;//!
    TH1F*      mLambdaM;//!
    TH1F*      mLbarM;//!  
    TH1F*      mKsRM;//!
    TH1F*      mLambdaRM;//!
    TH1F*      mLbarRM;//!
    
    TH2F*      mdedxvspt;//!
    TH2F*      mKsDecayL;//!
    TH2F*      mLambdaDecayL;//!
    TH2F*      mLbarDecayL;//!

    TH1F* hTriggerPt;//!
    
    // jet histos
    TH1F* hJetPt;//!
    TH1F* hJetCorrPt;//!
    TH1F* hJetPt2;//!
    TH1F* hJetE;//!
    TH1F* hJetEta;//!
    TH1F* hJetPhi;//!
    TH1F* hJetNEF;//!
    TH1F* hJetArea;//!

    TH2  *fHistJetHEtaPhi;//!

    // QA histos
    TH1  *fHistEventSelectionQA;//! 
    TH1  *fHistEventSelectionQAafterCuts;//!

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

    // counters
    Int_t GetEventCounter() {return mEventCounter;}
    Int_t GetAllPVEventCounter() {return mAllPVEventCounter;}
    Int_t GetInputEventCounter() {return mInputEventCounter;}
                
    ClassDef(StMyAnalysisMaker, 1)
};

#endif