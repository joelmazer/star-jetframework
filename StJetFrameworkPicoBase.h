#ifndef StJetFrameworkPicoBase_h
#define StJetFrameworkPicoBase_h

// some includes
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include <set>

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
class StPicoBTowHit;

// jet-framework classes
class StEmcPosition2;
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;
class StEventPlaneMaker;
class StCentMaker;

class StJetFrameworkPicoBase : public StMaker {
  public:

    // debug flags for specifics
    enum fDebugFlagEnum_t {
      kDebugNothing, // don't want lowest elements to be used
      kDebugFillJets = 11,
      kDebugMixedEvents = 12,
      kDebugEmcTrigger = 13
      //kDebugGeneralEvt,
      //kDebugCentrality,
      //kDebugEventPlaneCalc
    };

    // debug flags for tower lists
    enum fBadTowerListsEnum_t {
      kAltBadTow, 
      kBadTow1, kBadTow2, kBadTow3, kBadTow4, kBadTow5
    };

    // bad run list type enumerator
    enum fBadRunListsEnum_t {
      fBadRuns_w_missing_HT,
      fBadRuns_wo_missing_HT
    };

    // efficiency type enumerator
    enum fEfficiencyTypeEnum_t {
      kNormalPtEtaBased, // this is default
      kPtBased, kEtaBased, kHeaderArray
    };

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

    // jet recombination scheme enumerator
    enum ERecoScheme_t {
      E_scheme        = 0,
      pt_scheme       = 1, pt2_scheme      = 2,
      Et_scheme       = 3, Et2_scheme      = 4,
      BIpt_scheme     = 5, BIpt2_scheme    = 6,
      WTA_pt_scheme   = 7, WTA_modp_scheme = 8,
      external_scheme = 99
    };

    // jet hadronic correction type enumerator
    enum EHadCorrType_t {
      kLastMatchedTrack,
      kHighestEMatchedTrack,
      kAllMatchedTracks
    };

    // event plane track weight type enumerator
    enum EPtrackWeightTypeEnum {
      kNoWeight,
      kPtLinearWeight,
      kPtLinear2Const5Weight
    };

    // run flags for specifics - update this as needed: TODO
    enum fRunFlagEnum {
      Run11_pp500, // 500
      Run12_pp200, 
      Run12_pp500, // 500
      Run13_pp510, // 500
      Run14_AuAu200, Run14_AuAu200_MB,
      Run15_pp200,
      Run16_AuAu200,
      Run17_pp510  // 500
    };

    // trigger flags
    enum fEmcTriggerFlagEnum {
      kAny,
      kIsHT0, kIsHT1, kIsHT2, kIsHT3,
      kIsJP0, kIsJP1, kIsJP2
    };

    // MB flags for specifics
    enum fMBFlagEnum {
      kRun14main = 0,
      kRun16main = 1,
      kVPDMB     = 2,
      kVPDMB5    = 3,
      kVPDMB10   = 4,
      kVPDMB30   = 5,
      kVPDMB100  = 6,
      kVPDMBnovtx= 7,
      kRun12main = 8,
      kRun12alt  = 9
    };

    // trigger type used to run specific part of analysis
    enum fTriggerEventTypeEnum {
      kTriggerANY = 0,
      kTriggerMB  = 1, // corresponding to specific MB selection
      kTriggerHT  = 2, // corresponding to specific HT selection
      kTriggerJP  = 3, // corresponding to specific JP selection
      kTriggerMB30HT2HT3 = 4, // corresponding to selection used in some QA: MB30 || HT2 || HT3
      kTriggerMBHT2HT3 = 5    // corresponding to selection used in some QA: MB || HT2 || HT3  (for pp analysis where it lacks much MB30)
    };

    // Centrality Interfaces:
    //getRefMult3Corr();           // For refmult3
    //getTofTrayMultCorr();        // For TOF tray multiplicity
    //getgRefMultCorr();           // For grefmult //Run14 AuAu200GeV
    //getgRefMultCorr_P16id();     // For grefmult //Run14 AuAu200GeV, P16id
    //getgRefMultCorr_P17id_VpdMB30();        // For grefmult //Run14 AuAu200GeV, P17id for VPDMB-30; |vz| < 30
    //getgRefMultCorr_P18ih_VpdMB30();        // For grefmult //Run14 AuAu200GeV, P18ih for VPDMB-30; |vz| < 30 (new - June10, 2019)
    //getgRefMultCorr_P18ih_VpdMB30_AllLumi();// For grefmult //Run14 AuAu200GeV, P18ih for VPDMB-30; |vz| < 30 (new - Aug15, 2019), all luminosities def
    //getgRefMultCorr_VpdMB30();   // for VPDMB-30; |vz| < 30
    //getgRefMultCorr_VpdMBnoVtx();// for VPDMB-noVtx; |vz| < 100
    // centrality enum
    enum fCentralityDefEnum {
      krefmult, krefmult2, krefmult3,
      kgrefmult, kgrefmult_P16id, kgrefmult_VpdMBnoVtx, kgrefmult_VpdMB30,
      kgrefmult_P17id_VpdMB30, kgrefmult_P18ih_VpdMB30, kgrefmult_P18ih_VpdMB30_AllLumi    
    };    

    // centrality bin range
    enum fCentralityBinEnum {
      kCent010, kCent020,
      kCent1020, kCent1030, kCent1040,
      kCent2030, kCent2040, kCent2050, kCent2060,
      kCent3050, kCent3060,
      kCent4060, kCent4070, kCent4080,
      kCent5080,
      kCent6080
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

    // 'explicit' added below to constructor to remove cppcheck warning
    StJetFrameworkPicoBase();
    explicit StJetFrameworkPicoBase(const char *name);
    virtual ~StJetFrameworkPicoBase();
   
    // class required functions
    virtual Int_t Init();
    //virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    static TString GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, TClonesArray *partCont, TClonesArray *clusCont, TString tag);

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetdoppAnalysis(Bool_t pp)         { doppAnalysis      = pp;}
    virtual void            SetdoJetShapeAnalysis(Bool_t js)   { doJetShapeAnalysis = js; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }
    virtual void            SetTurnOnCentSelection(Bool_t o)   { fRequireCentSelection = o; }
    virtual void            SetCentralityBinCut(Int_t c)       { fCentralitySelectionCut = c; }

    // jet setters
    virtual void            SetJetType(Int_t jt)               { fJetType          = jt;}    // jet type (full, charged, neutral)
    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    virtual void            SetCorrectJetPt(Bool_t cpt)        { fCorrJetPt = cpt; }
    
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
    virtual void            SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
    virtual void            SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
    virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
    virtual void            SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }
    virtual void            SetBadRunListVers(Int_t i)         { fBadRunListVers = i; }
    virtual void            SetBadTowerListVers(UInt_t ibt)    { fBadTowerListVers = ibt; }

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

    // leading/subleading jets
    Double_t                GetDiJetAj(StJet *jet1, StJet *jet2, StRhoParameter *eventRho = 0x0, Bool_t doCorrJetPt = kFALSE);
    StJet                  *GetLeadingJet(TString fJetMakerNametemp, StRhoParameter *eventRho = 0x0);
    StJet                  *GetSubLeadingJet(TString fJetMakerNametemp, StRhoParameter *eventRho = 0x0);
 
    virtual void            SetExcludeLeadingJetsFromFit(Float_t n)   {fExcludeLeadingJetsFromFit = n; }
    virtual void            SetEventPlaneTrackWeight(int weight)      {fTrackWeight = weight; }

    // set names of makers for global use
    virtual void            SetOutputFileName(const char *on)         { mOutName = on; }
    void                    SetOutFileNameEP(TString epout)           { mOutNameEP = epout; }
    void                    SetOutFileNameQA(TString QAout)           { mOutNameQA = QAout; }
    // TODO add MIXED event name
    virtual void            SetJetMakerName(const char *jn)           { fJetMakerName = jn; }
    virtual void            SetJetBGMakerName(const char *bjn)        { fJetBGMakerName = bjn; }
    virtual void            SetRhoMakerName(const char *rn)           { fRhoMakerName = rn; }
    virtual void            SetRhoSparseMakerName(const char *rpn)    { fRhoSparseMakerName = rpn; }
    virtual void            SetEventPlaneMakerName(const char *epn)   { fEventPlaneMakerName = epn; }

    // add-to histogram name
    virtual void            AddToHistogramsName(TString add)           { fAddToHistogramsName = add  ; }
    virtual TString         GetAddedHistogramsStringToName() const     { return fAddToHistogramsName ; }

    // bad and dead tower list functions and arrays
    void                    ResetBadTowerList( );
    void                    ResetDeadTowerList( );
    Bool_t                  AddBadTowers(TString csvfile);
    Bool_t                  AddDeadTowers(TString csvfile);
    Bool_t                  IsTowerOK( Int_t mTowId );
    Bool_t                  IsTowerDead( Int_t mTowId );

    // bad run list 
    void                    ResetBadRunList( );
    Bool_t                  AddBadRuns(TString csvfile);
    Bool_t                  IsRunOK( Int_t mRunId );

    // return setup object for use in other classes
    std::set<Int_t>         GetBadTowers()                   { return badTowers          ; }
    std::set<Int_t>         GetDeadTowers()                  { return deadTowers         ; }
    std::set<Int_t>         GetBadRuns()                     { return badRuns            ; }

    Bool_t                  SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut); // centrality bin to cut on for analysis
    Bool_t                  DoComparison(int myarr[], int elems);
    Bool_t                  CheckForMB(Int_t RunFlag, Int_t type);
    Bool_t                  CheckForHT(Int_t RunFlag, Int_t type);

    // functions
    Double_t                ApplyTrackingEff(Bool_t applyEff, Double_t tpt, Double_t teta, Int_t cbin, Double_t ZDCx, Int_t effType, TFile *infile); // single-track reconstruction efficiency 
    Int_t                   GetRunNo(Int_t RunFlag, Int_t runid);
    Int_t                   GetNDataSetRuns(Int_t RunFlag);

  protected:
    TH1                    *FillEventTriggerQA(TH1 *h);               // filled event trigger QA plots
    Int_t                   GetCentBin(Int_t cent, Int_t nBin) const; // centrality bin
    Int_t                   GetCentBin10(Int_t cbin) const;           // centrality bin (10% size)
    Int_t                   Get4CentBin(Double_t scaledCent) const;
    Double_t                RelativePhi(Double_t mphi,Double_t vphi) const;               // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const;         // relative jet event plane angle
    Bool_t                  AcceptJet(StJet *jet);                                   // jets accept cuts function
    Bool_t                  AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert); // track accept cuts function
    //Bool_t                  AcceptTower(StPicoBTowHit *tower, TVector3 Vertex, Int_t towerID);     // tower accept cuts function
    Double_t                GetReactionPlane(); // get reaction plane angle
    Int_t                   EventCounter();     // when called, provides Event #
    Double_t                GetRhoValue(TString fRhoMakerNametemp);
    Bool_t                  GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, StPicoEvent *PicoEvent, Int_t towerID) const;
    Double_t                GetMaxTrackPt();               // find max track pt in event
    Double_t                GetMaxTowerEt();               // find max tower Et in event
    Int_t                   GetAnnuliBin(Double_t deltaR) const;
    Int_t                   GetJetPtBin(Double_t jetpt) const;
    Int_t                   GetJetEPBin(Double_t dEP) const;
    Int_t                   GetLuminosityBin(Double_t lumi) const;
    Double_t                GetDeltaR(StJet *jet, StPicoTrack *trk);
    Int_t                   GetVzRegion(double Vz);
    Int_t                   GetZVertex4cmBin(Double_t zvertex) const;


    static Double_t        *GenerateFixedBinArray(Int_t n, Double_t min, Double_t max);
    static void             GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t *array);

    // switches
    Bool_t                  doUsePrimTracks;         // primary track switch
    Int_t                   fDebugLevel;             // debug printout level
    Int_t                   fRunFlag;                // Run Flag enumerator value
    Bool_t                  doppAnalysis;            // use pp analysis data
    Bool_t                  doJetShapeAnalysis;      // perform jet shape analysis
    Bool_t                  fCorrJetPt;              // correct jet pt by rho
    Int_t                   fCentralityDef;          // Centrality Definition enumerator value
    Bool_t                  fRequireCentSelection;   // require particular centrality bin
    Bool_t                  doUseBBCCoincidenceRate; // use BBC or ZDC Coincidence Rate, kFALSE = ZDC

    // centrality    
    Double_t                fCentralityScaled;       // scaled by 5% centrality 
    Int_t                   ref16;                   // multiplicity bin (16)
    Int_t                   ref9;                    // multiplicity bin (9)

    // event
    Double_t                Bfield;                  // event Bfield
    TVector3                mVertex;                 // event vertex 3-vector
    Double_t                zVtx;                    // z-vertex component
    Double_t                fMaxEventTrackPt;        // max track pt in the event (to cut on)    
    Double_t                fMaxEventTowerEt;        // max tower Et in the event (to cut on)    
    Bool_t                  doRejectBadRuns;         // switch to reject bad runs and thus skip from analysis
    Int_t                   fBadRunListVers;         // version of bad runs file list to use
    UInt_t                  fBadTowerListVers;       // version of bad tower file list to use

    // cuts
    Int_t                   fJetType;                // jet type (full, charged, neutral)
    Double_t                fMinPtJet;               // min jet pt to keep jet in output
    Double_t                fJetConstituentCut;      // min jet constituent
    Double_t                fTrackBias;              // high pt track in jet bias
    Double_t                fTowerBias;              // high E tower in jet bias
    Double_t                fJetRad;                 // jet radius
    Double_t                fEventZVtxMinCut;        // min event z-vertex cut
    Double_t                fEventZVtxMaxCut;        // max event z-vertex cut
    Int_t                   fCentralitySelectionCut; // centrality selection cut
    Double_t                fTrackPtMinCut;          // min track pt cut
    Double_t                fTrackPtMaxCut;          // max track pt cut
    Double_t                fTrackPhiMinCut;         // min track phi cut
    Double_t                fTrackPhiMaxCut;         // max track phi cut
    Double_t                fTrackEtaMinCut;         // min track eta cut
    Double_t                fTrackEtaMaxCut;         // max track eta cut
    Double_t                fTrackDCAcut;            // max track dca cut
    Int_t                   fTracknHitsFit;          // requirement for track hits
    Double_t                fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax
    Double_t                fTowerEMinCut;           // min tower energy cut
    Double_t                fTowerEMaxCut;           // max tower energy cut
    Double_t                fTowerEtaMinCut;         // min tower eta cut
    Double_t                fTowerEtaMaxCut;         // max tower eta cut
    Double_t                fTowerPhiMinCut;         // min tower phi cut
    Double_t                fTowerPhiMaxCut;         // max tower phi cut

    // used for event plane calculation and resolution
    StJet                  *fLeadingJet;//! leading jet
    StJet                  *fSubLeadingJet;//! sub-leading jet
    Float_t                 fExcludeLeadingJetsFromFit;    // exclude n leading jets from fit
    Int_t                   fTrackWeight; // track weight for Q-vector summation

    TClonesArray           *CloneAndReduceTrackList(TClonesArray *tracks);

    // clonesarray collections of tracks and jets
    TClonesArray           *fTracksME;//! track collection to slim down for mixed events
    TClonesArray           *fJets;//!  jet array
    TClonesArray           *fBGJets;//! background jets array
    // added for Thomas - multiple jet collections
    TClonesArray           *fJets1;//! jet array
    TClonesArray           *fJets2;//! jet array

    // PicoDstMaker and PicoDst object pointer
    StPicoDstMaker         *mPicoDstMaker;
    StPicoDst              *mPicoDst;
    StPicoEvent            *mPicoEvent;
    StJetMakerTask         *JetMaker;
    StJetMakerTask         *JetMakerBG;
    StJetMakerTask         *JetMaker1; // for thomas, multiple jet collections
    StJetMakerTask         *JetMaker2; // for thomas, multiple jet collection
    StRho                  *RhoMaker;
    StRho                  *RhoMaker1; // for thomas, multiple jet collections
    StRho                  *RhoMaker2; // for thomas, multiple jet collections
    StEventPlaneMaker      *EventPlaneMaker;
    StCentMaker            *mCentMaker;

    // position object
    StEmcPosition2         *mEmcPosition;

    // centrality objects
    StRefMultCorr          *grefmultCorr;
    StRefMultCorr          *refmultCorr;
    StRefMultCorr          *refmult2Corr;

    // output file name string 
    TString                 mOutName;
    TString                 mOutNameEP;
    TString                 mOutNameQA;

    // maker names
    TString                 fJetMakerName;
    TString                 fJetBGMakerName;
    TString                 fRhoMakerName;
    TString                 fRhoSparseMakerName;
    TString                 fEventPlaneMakerName;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);
    StRhoParameter         *fRho;//!<!          // event rho
    StRhoParameter         *fRho1;//!<!         // temp for thomas
    StRhoParameter         *fRho2;//!<!         // temp for thomas
    Double_t                fRhoVal;//!<!       // event rho value, same for local rho
    TString                 fRhoName;///<       // rho name

    // add to histograms name
    TString                 fAddToHistogramsName; ///<  Add this string to histograms name.

    // counters 
    Int_t                   mEventCounter;//!
    Int_t                   mAllPVEventCounter;//!
    Int_t                   mInputEventCounter;//!

    // counters get functions
    Int_t GetEventCounter()      { return mEventCounter;}
    Int_t GetAllPVEventCounter() { return mAllPVEventCounter;}
    Int_t GetInputEventCounter() { return mInputEventCounter;}

  private:
    // bad and dead tower list functions and arrays
    std::set<Int_t>        badTowers;
    std::set<Int_t>        deadTowers;

    // bad run list 
    std::set<Int_t>        badRuns;

    ClassDef(StJetFrameworkPicoBase, 2)
};

/**
 * Generate array with fixed binning within min and max with n bins. The parameter array
 * will contain the bin edges set by this function. Attention, the array needs to be
 * provided from outside with a size of n+1
 * @param[in] n Number of bins
 * @param[in] min Minimum value for the binning
 * @param[in] max Maximum value for the binning
 * @param[out] array Array containing the bin edges
 */
inline void StJetFrameworkPicoBase::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t *array)
{
    Double_t binWidth = (max-min)/n;
    array[0] = min;
    for (Int_t i = 1; i <= n; i++) {
      array[i] = array[i-1]+binWidth;
    }
}

/**
 * Generate array with fixed binning within min and max with n bins. The array containing the bin
 * edges set will be created by this function. Attention, this function does not take care about
 * memory it allocates - the array needs to be deleted outside of this function
 * @param[in] n Number of bins
 * @param[in] min Minimum value for the binning
 * @param[in] max Maximum value for the binning
 * @return Array containing the bin edges created bu this function
 */
inline Double_t *StJetFrameworkPicoBase::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max)
{
    Double_t *array = new Double_t[n+1];
    GenerateFixedBinArray(n, min, max, array);
    return array;
}

#endif
