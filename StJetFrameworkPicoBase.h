#ifndef StJetFrameworkPicoBase_h
#define StJetFrameworkPicoBase_h

// some includes
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"

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

class StJetFrameworkPicoBase : public StMaker {
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

    // event plane track weight type enumerator
    enum EPtrackWeightType {
      kNoWeight,
      kPtLinearWeight,
      kPtLinear2Const5Weight
    };

    // run flags for specifics - update this
    enum fRunFlagEnum {
      Run14_AuAu200,
      Run16_AuAu200
    };

    // run flags for specifics
    enum fTriggerFlagEnum {
      kAny,
      kIsHT0, kIsHT1, kIsHT2, kIsHT3,
      kIsJP0, kIsJP1, kIsJP2
    };

    // Centrality Interfaces:
    //getRefMult3Corr() ; // For refmult3
    //getTofTrayMultCorr() ; // For TOF tray multiplicity
    //getgRefMultCorr()  ; // For grefmult //Run14 AuAu200GeV
    //getgRefMultCorr_P16id()  ; // For grefmult //Run14 AuAu200GeV, P16id
    //getgRefMultCorr_VpdMB30()  ; // for VPDMB-30; |vz| < 30
    //getgRefMultCorr_VpdMBnoVtx()  ; //  for VPDMB-noVtx; |vz| < 100
    // centrality enum
    enum fCentralityDefEnum {
      krefmult, krefmult2, krefmult3,
      kgrefmult, kgrefmult_P16id, kgrefmult_VpdMBnoVtx, kgrefmult_VpdMB30
    };    

    StJetFrameworkPicoBase();
    StJetFrameworkPicoBase(const char *name);
    virtual ~StJetFrameworkPicoBase();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    static TString GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray* partCont, TClonesArray* clusCont, TString tag);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }

    // jet setters
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

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    // leading jet
    StJet*                  GetLeadingJet(StRhoParameter* eventRho = 0x0);
    virtual void            SetExcludeLeadingJetsFromFit(Float_t n)         {fExcludeLeadingJetsFromFit = n; }
    virtual void            SetEventPlaneTrackWeight(int weight)            {fTrackWeight = weight; }

  protected:
    Int_t                  GetCentBin(Int_t cent, Int_t nBin) const; // centrality bin
    Double_t               RelativePhi(Double_t mphi,Double_t vphi) const; // relative jet track angle
    Double_t               RelativeEPJET(Double_t jetAng, Double_t EPAng) const;  // relative jet event plane angle
    Bool_t                 AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function
    Double_t               GetReactionPlane(); // get reaction plane angle

    // switches
    Bool_t                 doUsePrimTracks;         // primary track switch
    Int_t                  fDebugLevel;             // debug printout level
    Int_t                  fRunFlag;                // Run Flag enumerator value
    Int_t                  fCentralityDef;          // Centrality Definition enumerator value
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

    // used for event plane calculation and resolution
    StJet*         fLeadingJet;//! leading jet
    Float_t        fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
    Int_t          fTrackWeight; // track weight for Q-vector summation

    TClonesArray      *CloneAndReduceTrackList(TClonesArray* tracks);

    // clonesarray collections of tracks and jets
    TClonesArray          *fTracksME;//! track collection to slim down for mixed events
    TClonesArray          *fJets;//! jet collection

    // PicoDstMaker and PicoDst object pointer
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent    *mPicoEvent;
    StJetMakerTask *JetMaker;
    StRho          *RhoMaker;

    // centrality objects
    StRefMultCorr* grefmultCorr;
    StRefMultCorr* refmultCorr;
    StRefMultCorr* refmult2Corr;

    // output file name string 
    TString      mOutName;

    // counters 
    Int_t        mEventCounter;//!
    Int_t        mAllPVEventCounter;//!
    Int_t        mInputEventCounter;//!

    // Rho objects
    StRhoParameter        *GetRhoFromEvent(const char *name);
    StRhoParameter        *fRho;//!<!          // event rho
    Double_t               fRhoVal;//!<!       // event rho value, same for local rho
    TString                fRhoName;///<       // rho name

    // maker names
    TString                fJetMakerName;
    TString                fRhoMakerName;

    // counters
    Int_t GetEventCounter() {return mEventCounter;}
    Int_t GetAllPVEventCounter() {return mAllPVEventCounter;}
    Int_t GetInputEventCounter() {return mInputEventCounter;}
                
  private:
    ClassDef(StJetFrameworkPicoBase, 1)
};

#endif
