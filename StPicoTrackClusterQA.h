#ifndef STPICOTRACKCLUSTERQA_H
#define STPICOTRACKCLUSTERQA_H
// $Id$

// base class
//class StJetFrameworkPicoBase;
//class StJetMakerTask;

#include <set>

// includes
//#include "StJetFrameworkPicoBase.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH2;
class THnSparse;
class TProfile;

// STAR classes
class StMaker;
class StChain;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;

class StEmcGeom;
class StBemcTables; //v3.14
class StEmcCluster;
class StEmcCollection;
class StEmcPosition;

// centrality class
class StRefMultCorr;

#include "StMyAnalysisMaker.h"

/*  Used to store track & tower matching 
 *  information between computation steps     
 */
struct BemcMatch {
  Int_t globalId;
  Int_t trackId;
  Double_t trackEta;
  Double_t trackPhi;
  Double_t matchEta;
  Double_t matchPhi;
  
  BemcMatch() : globalId(-1), trackId(-1), trackEta(0.0), trackPhi(0.0), matchEta(0.0), matchPhi(0.0) {};
  BemcMatch(int id, int trkId, double trackEta, double trackPhi, double matchEta, double matchPhi) :
  globalId(id), trackId(trkId), trackEta(trackEta), trackPhi(trackPhi), matchEta(matchEta), matchPhi(matchPhi) {};
  
};


class StPicoTrackClusterQA : public StMaker {
//class StPicoTrackClusterQA : public StJetFrameworkPicoBase {
 public:

  // debug flags for specifics
  enum fDebugFlagEnum {
    kDebugNothing, // don't want lowest elements to be used
    kDebugEmcTrigger,
    kDebugGeneralEvt,
    kDebugCentrality,
  };

  enum towerMode{AcceptAllTowers=0, RejectBadTowerStatus=1};

  StPicoTrackClusterQA();
  StPicoTrackClusterQA(const char *name, bool dohistos, const char* outName);
  virtual ~StPicoTrackClusterQA();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  // Use one set to reject bad towers
  void ResetBadTowerList( );
  void ResetDeadTowerList( );
  Bool_t AddBadTowers(TString csvfile);
  Bool_t AddDeadTowers(TString csvfile);

  // THnSparse Setup
  virtual THnSparse*      NewTHnSparseFTracks(const char* name, UInt_t entries);
  virtual void GetDimParamsTracks(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual THnSparse*      NewTHnSparseFTowers(const char* name, UInt_t entries);
  virtual void GetDimParamsTowers(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  // switches
  virtual void         SetUsePrimaryTracks(Bool_t P)    { doUsePrimTracks       = P; } 
  virtual void         SetDebugLevel(Int_t l)           { fDebugLevel           = l; }
  virtual void         SetRunFlag(Int_t f)              { fRunFlag              = f; }
  virtual void         SetdoppAnalysis(Bool_t pp)       { doppAnalysis          = pp;}
  virtual void         SetCentralityDef(Int_t c)        { fCentralityDef        = c; }
  virtual void         SetTurnOnCentSelection(Bool_t o) { fRequireCentSelection = o; }
  virtual void         SetCentralityBinCut(Int_t c)     { fCentralitySelectionCut = c; }

  // event setters
  virtual void         SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
  virtual void         SetUseBBCCoincidenceRate(Bool_t b) { doUseBBCCoincidenceRate = b; }
  virtual void         SetMaxEventTrackPt(Double_t mxpt) { fMaxEventTrackPt = mxpt; }

  // track / cluster setters 
  virtual void         SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;}   // min track cut
  virtual void         SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;}   // max track cut
  virtual void         SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
  virtual void         SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
  virtual void         SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
  virtual void         SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
  virtual void         SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }
  virtual void         SetMinClusterPt(Double_t minpt)    { fClusterPtMinCut    = minpt;} // min cluster cut
  virtual void         SetMaxClusterPt(Double_t maxpt)    { fClusterPtMaxCut    = maxpt;} // max cluster cut

  // tower setters
  virtual void         SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
  virtual void         SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
  virtual void         SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

  // event selection
  virtual void         SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
  virtual void         SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }       
  virtual void         SetDoTowerQAforHT(Bool_t m)        { fDoTowerQAforHT = m; }

  // efficiency correction setter
  virtual void         SetDoEffCorr(Int_t effcorr)        { fDoEffCorr = effcorr; }

  // common setters
  void                 SetClusName(const char *n)       { fCaloName      = n;  }
  void                 SetTracksName(const char *n)     { fTracksName    = n;  }
  void                 SetTrackEfficiency(Double_t t)   { fTrackEfficiency  = t     ; }

  Double_t             GetTrackEfficiency()             { return fTrackEfficiency   ; }

  /* define if tower status should be used to reject towers, or if all
   * towers should be accepted - default is to accept all towers, then
   * generate a bad tower list for the entire data set.
  */
  void                 SetTowerAcceptMode(towerMode mode) { mTowerStatusMode = mode; }

  /* set the minimum tower energy to be reconstructed (default = 0.15) */
  void                 SetTowerEnergyMin(double mMin)     { mTowerEnergyMin = mMin; }

  // set hadronic correction fraction for matched tracks to towers
  void                 SetHadronicCorrFrac(float frac)    { mHadronicCorrFrac = frac; }

 protected:
  void                   RunQA();
  void                   RunTowerTest();
  void                   RunFiredTriggerQA();  
  Bool_t                 AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function
  Bool_t                 AcceptTower(StPicoBTowHit *tower);                              // tower accept cuts function
  Int_t                  GetCentBin(Int_t cent, Int_t nBin) const;                       // centrality bin
  Bool_t                 SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut); // centrality bin to cut on for analysis
  TH1*                   FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
  TH1*                   FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
  Bool_t                 DoComparison(int myarr[], int elems);
  Bool_t                 CheckForMB(int RunFlag, int type);
  Bool_t                 CheckForHT(int RunFlag, int type);
  Double_t               GetMaxTrackPt();

  void                   SetSumw2(); // set errors weights 

  // switches
  Bool_t                 doWriteHistos;           // write QA histos
  Bool_t                 doUsePrimTracks;         // primary track switch
  Int_t                  fDebugLevel;             // debug printout level
  Int_t                  fRunFlag;                // Run Flag numerator value
  Bool_t                 doppAnalysis;            // use pp analysis data
  Int_t                  fCentralityDef;          // Centrality Definition enumerator value
  Bool_t                 fDoEffCorr;              // efficiency correction to tracks
  Bool_t                 fDoTowerQAforHT;         // do tower QA for HT triggers (else do for MB) - temp

  // event cuts
  Double_t               fEventZVtxMinCut;        // min event z-vertex cut
  Double_t               fEventZVtxMaxCut;        // max event z-vertex cut
  Int_t                  fCentralitySelectionCut; // centrality selection cut
  Bool_t                 fRequireCentSelection;   // require particular centrality bin
  Bool_t                 doUseBBCCoincidenceRate; // use BBC or ZDC Coincidence Rate, kFALSE = ZDC
  Double_t               fMaxEventTrackPt;        // max track pt in the event (to cut on) 

  // names
  TString                mOutName;                // name of output file
  TString                fAnalysisMakerName;      // name of this analysis maker
  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection

  Double_t               fTrackPtMinCut;          // min track pt cut
  Double_t               fTrackPtMaxCut;          // max track pt cut
  Double_t               fClusterPtMinCut;        // min cluster pt cut
  Double_t               fClusterPtMaxCut;        // max cluster pt cut
  Double_t               fTrackPhiMinCut;         // min track phi cut
  Double_t               fTrackPhiMaxCut;         // max track phi cut
  Double_t               fTrackEtaMinCut;         // min track eta cut
  Double_t               fTrackEtaMaxCut;         // max track eta cut
  Double_t               fTrackDCAcut;            // max track dca cut
  Int_t                  fTracknHitsFit;          // requirement for track hits
  Double_t               fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)
  Int_t                  fGoodTrackCounter;       // good tracks - passed quality cuts
  Double_t               fTowerEMinCut;           // min tower energy cut
  Double_t               fTowerEMaxCut;           // max tower energy cut
  Double_t               fTowerEtaMinCut;         // min tower eta cut
  Double_t               fTowerEtaMaxCut;         // max tower eta cut
  Double_t               fTowerPhiMinCut;         // min tower phi cut
  Double_t               fTowerPhiMaxCut;         // max tower phi cut

  // centrality    
  Double_t        fCentralityScaled;           // scaled by 5% centrality 
  Int_t           ref16;                       // multiplicity bin (16)
  Int_t           ref9;                        // multiplicity bin (9)

  // event
  Float_t         Bfield;                      // event Bfield
  StThreeVectorF  mVertex;                     // event vertex 3-vector
  Double_t        zVtx;                        // z-vertex component
  Int_t           fRunNumber;                  // Run number

  // event selection types
  UInt_t          fEmcTriggerEventType;        // Physics selection of event used for signal
  UInt_t          fMBEventType;                // Physics selection of event used for MB
  Int_t           fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA
  Bool_t          fTowerToTriggerTypeHT1[4801];// Tower with corresponding HT1 trigger type array
  Bool_t          fTowerToTriggerTypeHT2[4801];// Tower with corresponding HT2 trigger type array
  Bool_t          fTowerToTriggerTypeHT3[4801];// Tower with corresponding HT3 trigger type array

  // Emc objects
  StEmcGeom             *mGeom;
  StEmcCollection       *mEmcCol; 
  StBemcTables          *mBemcTables; 
  std::vector<BemcMatch> mBemcMatchedTracks;

  towerMode              mTowerStatusMode;
  Double_t               mTowerEnergyMin;
  Float_t                mHadronicCorrFrac;

 private:
  Bool_t MuProcessBEMC();
  Bool_t PicoProcessBEMC();
  Int_t  MuFindSMDClusterHits(StEmcCollection* coll, Double_t eta, Double_t phi, Int_t detectorID);

  StMuDstMaker      *mMuDstMaker;   // MuDstMaker object
  StMuDst           *mMuDst;        // muDst object
  StMuEvent         *mMuInputEvent; // muDst event object
  StPicoDstMaker    *mPicoDstMaker; // PicoDstMaker object
  StPicoDst         *mPicoDst;      // PicoDst object
  StPicoEvent       *mPicoEvent;    // PicoEvent object

  //bool              *mTowerStatusArr; // tower status array

  // centrality objects
  StRefMultCorr   *grefmultCorr;

  // histograms
  TH1F           *fHistNTrackvsPt;//!
  TH1F           *fHistNTrackvsPhi;//!
  TH1F           *fHistNTrackvsEta;//!
  TH2F           *fHistNTrackvsPhivsEta;//!
  TH1F           *fHistNTowervsE;//!
  TH1F           *fHistNTowervsEt;//!
  TH1F           *fHistNTowervsPhi;//!
  TH1F           *fHistNTowervsEta;//!
  TH2F           *fHistNTowervsPhivsEta;//!

  // QA histos
  TH1            *fHistEventSelectionQA;//! 
  TH1            *fHistEventSelectionQAafterCuts;//!
  TH1            *hTriggerIds;//!
  TH1            *hEmcTriggers;//!

  // trigger histos for zero and negative energy
  TH1F           *fHistNZeroEHT1vsID;//!
  TH1F           *fHistNZeroEHT2vsID;//!
  TH1F           *fHistNZeroEHT3vsID;//!
  TH1F           *fHistNNegEHT1vsID;//!
  TH1F           *fHistNNegEHT2vsID;//!
  TH1F           *fHistNNegEHT3vsID;//!

  // trigger histos - firing towers QA
  TH1F           *fHistNFiredHT0vsID;//!
  TH1F           *fHistNFiredHT1vsID;//!
  TH1F           *fHistNFiredHT2vsID;//!
  TH1F           *fHistNFiredHT3vsID;//!
  TH1F           *fHistHT0FiredEtvsID;//!
  TH1F           *fHistHT1FiredEtvsID;//!
  TH1F           *fHistHT2FiredEtvsID;//!
  TH1F           *fHistHT3FiredEtvsID;//!
  TH2F           *fHistHT0IDvsFiredEt;//!
  TH2F           *fHistHT1IDvsFiredEt;//!
  TH2F           *fHistHT2IDvsFiredEt;//!
  TH2F           *fHistHT3IDvsFiredEt;//!

  // THn Sparse's
  THnSparse      *fhnTrackQA;//!      // sparse of track info
  THnSparse      *fhnTowerQA;//!      // sparse of tower info

  // bad and dead tower list functions and arrays
  Bool_t IsTowerOK( Int_t mTowId );
  Bool_t IsTowerDead( Int_t mTowId );
  std::set<Int_t> badTowers;
  std::set<Int_t> deadTowers;

  StPicoTrackClusterQA(const StPicoTrackClusterQA&);            // not implemented
  StPicoTrackClusterQA &operator=(const StPicoTrackClusterQA&); // not implemented

  ClassDef(StPicoTrackClusterQA, 1) // track/cluster QA task
};
#endif
