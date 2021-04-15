#ifndef STPICOTRACKCLUSTERQA_H
#define STPICOTRACKCLUSTERQA_H
// $Id$

// C++ includes
#include <set>

// STAR includes
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StEmcUtil/geometry/StEmcGeom.h"

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH1F;
class TH2;
class TH2F;
class THnSparse;
class TProfile;
class TVector3;

// STAR classes
class StMaker;
class StChain;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoBTowHit;

// EMC and tower related classes
class StEmcGeom;

// star jet-frameworks classes
class StJetFrameworkPicoBase;
class StEmcPosition2;
class StCentMaker;

// centrality class
class StRefMultCorr;

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
    kDebugTrackTowerMatching
  };

  StPicoTrackClusterQA();
  StPicoTrackClusterQA(const char *name, bool dohistos, const char *outName);
  virtual ~StPicoTrackClusterQA();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  // THnSparse Setup
  virtual THnSparse      *NewTHnSparseFTracks(const char *name, UInt_t entries);
  virtual void GetDimParamsTracks(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  virtual THnSparse      *NewTHnSparseFTowers(const char *name, UInt_t entries);
  virtual void GetDimParamsTowers(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);

  // switches
  virtual void         SetUsePrimaryTracks(Bool_t P)    { doUsePrimTracks       = P; } 
  virtual void         SetDebugLevel(Int_t l)           { fDebugLevel           = l; }
  virtual void         SetRunFlag(Int_t f)              { fRunFlag              = f; }
  virtual void         SetdoppAnalysis(Bool_t pp)       { doppAnalysis          = pp;}
  virtual void         SetTurnOnCentSelection(Bool_t o) { fRequireCentSelection = o; }
  virtual void         SetCentralityBinCut(Int_t c)     { fCentralitySelectionCut = c; }

  // event setters
  virtual void         SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
  virtual void         SetMaxEventTrackPt(Double_t mxpt)  { fMaxEventTrackPt = mxpt; }
  virtual void         SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
  virtual void         SetRejectBadRuns(Bool_t rj)        { doRejectBadRuns = rj; }

  // track / cluster setters 
  virtual void         SetTrackPtRange(Double_t ptmi, Double_t ptma) { fTrackPtMinCut = ptmi; fTrackPtMaxCut = ptma; }
  virtual void         SetTrackPhiRange(Double_t phimi, Double_t phima) { fTrackPhiMinCut = phimi; fTrackPhiMaxCut = phima; }
  virtual void         SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
  virtual void         SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
  virtual void         SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
  virtual void         SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

  // tower setters
  virtual void         SetTowerERange(Double_t enmi, Double_t enmx) { fTowerEMinCut = enmi; fTowerEMaxCut = enmx; }
  virtual void         SetTowerEtaRange(Double_t temi, Double_t temx) { fTowerEtaMinCut = temi; fTowerEtaMaxCut = temx; }
  virtual void         SetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fTowerPhiMinCut = tpmi; fTowerPhiMaxCut = tpmx; }

  // event selection
  virtual void         SetTriggerToUse(UInt_t ttu)        { fTriggerToUse = ttu; }
  virtual void         SetEmcTriggerEventType(UInt_t te)  { fEmcTriggerEventType = te; }
  virtual void         SetMBEventType(UInt_t mbe)         { fMBEventType = mbe; }       
  virtual void         SetDoTowerQAforHT(Bool_t m)        { fDoTowerQAforHT = m; }

  // efficiency correction setter
  virtual void         SetDoEffCorr(Bool_t effcorr)     { fDoEffCorr = effcorr; }

  // common setters
  void                 SetClusName(const char *n)       { fCaloName      = n;  }
  void                 SetTracksName(const char *n)     { fTracksName    = n;  }

  /* set the minimum tower energy to be reconstructed (default = 0.15) */
  void                 SetTowerEnergyMin(double mMin)     { mTowerEnergyMin = mMin; }

  // set hadronic correction fraction for matched tracks to towers
  void                 SetHadronicCorrFrac(float frac)    { mHadronicCorrFrac = frac; }
  void                 SetJetHadCorrType(Int_t hct)       { fJetHadCorrType = hct;}

 protected:
  // functions
  void                 RunEventQA();
  void                 RunTrackQA();
  void                 RunBEmcPidTraitsQA();
  void                 RunHadCorrTowerQA();
  void                 RunTowerQA();
  void                 RunFiredTriggerQA();  
  Bool_t               AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert);  // track accept cuts function
  Bool_t               AcceptTower(StPicoBTowHit *tower, Int_t towerID);         // tower accept cuts function
  Int_t                GetCentBin(Int_t cent, Int_t nBin) const;                 // centrality bin
  TH1                 *FillEmcTriggersHist(TH1 *h);                              // EmcTrigger counter histo
  TH1                 *FillEventTriggerQA(TH1 *h);                               // fill event trigger QA plots
  Bool_t               DoComparison(int myarr[], int elems);
  Double_t             GetMaxTrackPt();               // find max track pt in event
  Double_t             GetMaxTowerEt();               // find max tower Et in event
  void                 FillTriggerIDs(TH1 *h);
  void                 SetSumw2(); // set errors weights 

  // switches
  Bool_t               doWriteHistos;           // write QA histos
  Bool_t               doUsePrimTracks;         // primary track switch
  Int_t                fDebugLevel;             // debug printout level
  Int_t                fRunFlag;                // Run Flag numerator value
  Bool_t               doppAnalysis;            // use pp analysis data
  Bool_t               fDoEffCorr;              // efficiency correction to tracks
  Bool_t               fDoTowerQAforHT;         // do tower QA for HT triggers (else do for MB) - temp

  // event cuts
  Double_t             fMaxEventTrackPt;        // max track pt in the event (to cut on) 
  Double_t             fMaxEventTowerEt;        // max tower Et in the event (to cut on)    
  Bool_t               doRejectBadRuns;         // switch to reject bad runs and thus skip from analysis
  Double_t             fEventZVtxMinCut;        // min event z-vertex cut
  Double_t             fEventZVtxMaxCut;        // max event z-vertex cut
  Int_t                fCentralitySelectionCut; // centrality selection cut
  Bool_t               fRequireCentSelection;   // require particular centrality bin

  // names
  TString              mOutName;                // name of output file
  TString              fAnalysisMakerName;      // name of this analysis maker
  TString              fTracksName;             // name of track collection
  TString              fCaloName;               // name of calo cluster collection

  Double_t             fTrackPtMinCut;          // min track pt cut
  Double_t             fTrackPtMaxCut;          // max track pt cut
  Double_t             fTrackPhiMinCut;         // min track phi cut
  Double_t             fTrackPhiMaxCut;         // max track phi cut
  Double_t             fTrackEtaMinCut;         // min track eta cut
  Double_t             fTrackEtaMaxCut;         // max track eta cut
  Double_t             fTrackDCAcut;            // max track dca cut
  Int_t                fTracknHitsFit;          // requirement for track hits
  Double_t             fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax
  Int_t                fGoodTrackCounter;       // good tracks - passed quality cuts
  Double_t             fTowerEMinCut;           // min tower energy cut
  Double_t             fTowerEMaxCut;           // max tower energy cut
  Double_t             fTowerEtaMinCut;         // min tower eta cut
  Double_t             fTowerEtaMaxCut;         // max tower eta cut
  Double_t             fTowerPhiMinCut;         // min tower phi cut
  Double_t             fTowerPhiMaxCut;         // max tower phi cut

  // centrality    
  Double_t             fCentralityScaled;       // scaled by 5% centrality 
  Int_t                ref16;                   // multiplicity bin (16)
  Int_t                ref9;                    // multiplicity bin (9)

  // event
  Float_t              Bfield;                  // event Bfield
  TVector3             mVertex;                 // event vertex 3-vector
  Double_t             zVtx;                    // z-vertex component
  Int_t                fRunNumber;              // Run number

  // event selection types
  UInt_t               fTriggerToUse;               // trigger to use for analysis
  UInt_t               fEmcTriggerEventType;        // Physics selection of event used for signal
  UInt_t               fMBEventType;                // Physics selection of event used for MB
  Int_t                fEmcTriggerArr[8];           // EMCal triggers array: used to select signal and do QA
  Bool_t               fTowerToTriggerTypeHT1[4800];// Tower with corresponding HT1 trigger type array
  Bool_t               fTowerToTriggerTypeHT2[4800];// Tower with corresponding HT2 trigger type array
  Bool_t               fTowerToTriggerTypeHT3[4800];// Tower with corresponding HT3 trigger type array

  // Emc objects
  StEmcGeom           *mGeom;

  Double_t             mTowerEnergyMin;

  // hadronic correction
  Float_t              mHadronicCorrFrac;       // hadronic correction fraction to subtract
  Int_t                fJetHadCorrType;         // hadronic correction type to be used

 private:
  // OLD - part of tests
  Bool_t               PicoProcessBEMC();

  StPicoDstMaker      *mPicoDstMaker; // PicoDstMaker object
  StPicoDst           *mPicoDst;      // PicoDst object
  StPicoEvent         *mPicoEvent;    // PicoEvent object
  StCentMaker         *mCentMaker;    // Centrality maker object
  StJetFrameworkPicoBase *mBaseMaker; // Base maker object

  // position object
  StEmcPosition2      *mEmcPosition;

  // centrality objects
  StRefMultCorr       *grefmultCorr;

  // histograms
  TH1F           *fHistNTrackvsPt;//!
  TH1F           *fHistNTrackvsPhi;//!
  TH1F           *fHistNTrackvsEta;//!
  TH2F           *fHistNTrackvsPhivsEta;//!
  TH1F           *fHistTrackToTowerIndex;//!
  TH1F           *fHistNHadCorrTowervsE;//!
  TH1F           *fHistNHadCorrTowervsEt;//!
  TH1F           *fHistNHadCorrTowervsPhi;//!
  TH1F           *fHistNHadCorrTowervsEta;//!
  TH2F           *fHistNHadCorrTowervsPhivsEta;//!
  TH1F           *fHistNHadCorrTowerHOTvsTowID;//!
  TH1F           *fHistNTowervsADC;//!
  TH1F           *fHistNTowervsE;//!
  TH1F           *fHistNTowervsEt;//!
  TH1F           *fHistNTowervsPhi;//!
  TH1F           *fHistNTowervsEta;//!
  TH2F           *fHistNTowervsPhivsEta;//!
  TH1F           *fHistNTowerHOTvsTowID;//!

  // trigger / event selection QA histos
  TH1F           *fHistCentrality;//!
  TH1F           *fHistMultiplicity;//!
  TH1F           *fHistEventCounter;//!
  TH1F           *fHistEventSelectionQA;//! 
  TH1F           *fHistEventSelectionQAafterCuts;//!
  TH1F           *fHistEventSelectionTrg;//!
  TH1F           *hEmcTriggers;//!
  TH1F           *fHistTriggerIDs;//!

  // event QA
  TH1F           *fHistEventNTrig_MB30;//!
  TH1F           *fHistEventNTrig_HT;//!
  TH1F           *fHistRefMult_MB30;//!
  TH1F           *fHistVzVPDVz_MB30;//!
  TH2F           *fHistVyvsVx_MB30;//!
  TH1F           *fHistRvtx_MB30;//!
  TH1F           *fHistPerpvtx_MB30;//!
  TH1F           *fHistZvtx_MB30;//!
  TH1F           *fHistZDCx_MB30;//!
  TH1F           *fHistEventID_MB30;//!
  TH1F           *fHistRunID_MB30;//!
  TProfile       *fProfEventTrackPt_MB30;//!
  TProfile       *fProfEventRefMult_MB30;//!
  TProfile       *fProfEventXvtx_MB30;//!
  TProfile       *fProfEventYvtx_MB30;//!
  TProfile       *fProfEventZvtx_MB30;//!
  TProfile       *fProfEventRvtx_MB30;//!
  TProfile       *fProfEventPerpvtx_MB30;//!
  TProfile       *fProfEventBBCx_MB30;//!
  TProfile       *fProfEventZDCx_MB30;//!
  TProfile       *fProfEventTrackPt;//!
  TProfile       *fProfEventRefMult;//!
  TProfile       *fProfEventRanking;//!
  TProfile       *fProfEventZvtx;//!
  TProfile       *fProfEventYvtx;//!
  TProfile       *fProfEventXvtx;//!
  TProfile       *fProfEventVzVPD;//!
  TProfile       *fProfEventBBCx;//!
  TProfile       *fProfEventZDCx;//!
 
  // trigger histos for hot tower (threshold levels for varying bad tower lists)
  TProfile       *fProfTowerAvgEvsID;//!
  TProfile       *fProfTowerAvgEtvsID;//!
  TH1F           *fHistNFiredHT1vsIDEt200MeV;//!
  TH1F           *fHistNFiredHT2vsIDEt200MeV;//!
  TH1F           *fHistNFiredHT3vsIDEt200MeV;//!
  TH1F           *fHistNFiredHT1vsIDEt1000MeV;//!
  TH1F           *fHistNFiredHT2vsIDEt1000MeV;//!
  TH1F           *fHistNFiredHT3vsIDEt1000MeV;//!
  TH1F           *fHistNFiredHT1vsIDEt2000MeV;//!
  TH1F           *fHistNFiredHT2vsIDEt2000MeV;//!
  TH1F           *fHistNFiredHT3vsIDEt2000MeV;//!
  TH1F           *fHistNFiredvsIDEt200MeV;//!
  TH1F           *fHistNFiredvsIDEt1000MeV;//!
  TH1F           *fHistNFiredvsIDEt2000MeV;//!
 
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

  TH1F           *fHistNFiredHT0vsFlag;//!
  TH1F           *fHistNFiredHT1vsFlag;//!
  TH1F           *fHistNFiredHT2vsFlag;//!
  TH1F           *fHistNFiredHT3vsFlag;//!
  TH1F           *fHistNFiredJP0vsFlag;//!
  TH1F           *fHistNFiredJP1vsFlag;//!
  TH1F           *fHistNFiredJP2vsFlag;//!

  TH1F           *fHistNFiredHT0vsADC;//!
  TH1F           *fHistNFiredHT1vsADC;//!
  TH1F           *fHistNFiredHT2vsADC;//!
  TH1F           *fHistNFiredHT3vsADC;//!
  TH1F           *fHistNFiredJP0vsADC;//!
  TH1F           *fHistNFiredJP1vsADC;//!
  TH1F           *fHistNFiredJP2vsADC;//!

  // THn Sparse's
  THnSparse      *fhnTrackQA;//!      // sparse of track info
  THnSparse      *fhnTowerQA;//!      // sparse of tower info

  // bad and dead tower list
  std::set<Int_t>        badTowers;
  std::set<Int_t>        deadTowers;

  // bad run list 
  std::set<Int_t>        badRuns;

  StPicoTrackClusterQA(const StPicoTrackClusterQA&);            // not implemented
  StPicoTrackClusterQA &operator=(const StPicoTrackClusterQA&); // not implemented

  ClassDef(StPicoTrackClusterQA, 2) // track/cluster QA task
};
#endif
