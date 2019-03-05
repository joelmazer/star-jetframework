#ifndef STJETMAKERTASKBGSUB_H
#define STJETMAKERTASKBGSUB_H

// $Id$

#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"

#include <set>

// for clusters
#include "StEmcUtil/geometry/StEmcGeom.h"
///#include "StEmcUtil/projection/StEmcPosition.h"
///#include "StMuDSTMaker/COMMON/StMuDst.h"
class StEmcGeom;
class StEmcPosition2;

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH2;

// STAR classes
class StPicoDst;
class StPicoDstMaker;

// Jet classes
class StFJWrapper;
class StJetUtility;

// Centrality class
class StRefMultCorr;

// STAR includes
#include "StFJWrapper.h"
#include "FJ_includes.h"
#include "StJet.h"
#include "StMyAnalysisMaker.h"

namespace fastjet {
  class PseudoJet;
}

class StJetMakerTaskBGsub : public StMaker {
 public:

  // jet type enumerator
  enum EJetType_t {
    kFullJet,
    kChargedJet,
    kNeutralJet
  };

  StJetMakerTaskBGsub();
  StJetMakerTaskBGsub(const char *name, double mintrackPt, bool dohistos, const char* outName);
  virtual ~StJetMakerTaskBGsub();

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

  // switches
  virtual void         SetUsePrimaryTracks(Bool_t P)    { doUsePrimTracks       = P; } 
  virtual void         SetDebugLevel(Int_t l)           { fDebugLevel           = l; }
  virtual void         SetRunFlag(Int_t f)              { fRunFlag              = f; }
  virtual void         SetdoppAnalysis(Bool_t pp)       { doppAnalysis          = pp;}
  virtual void         SetCentralityDef(Int_t c)        { fCentralityDef        = c; }
  virtual void         SetTurnOnCentSelection(Bool_t o) { fRequireCentSelection = o; }
  virtual void         SetCentralityBinCut(Int_t c)     { fCentralitySelectionCut = c; }
  virtual void         SetdoConstituentSubtr(Bool_t c)  { doConstituentSubtr    = c; }

  // event setters
  virtual void         SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }
  virtual void         SetEmcTriggerEventType(UInt_t te)   { fEmcTriggerEventType = te; }
  virtual void         SetMBEventType(UInt_t mbe)       { fMBEventType = mbe; }   
  virtual void         SetTriggerToUse(UInt_t ttu)      { fTriggerToUse = ttu; }
  virtual void         SetBadTowerListVers(UInt_t ibt)  { fBadTowerListVers = ibt; }
  virtual void         SetMaxEventTrackPt(Double_t mxpt) { fMaxEventTrackPt = mxpt; }

  // common setters
  void         SetClusName(const char *n)                 { fCaloName      = n;  }
  void         SetTracksName(const char *n)               { fTracksName    = n;  }
  void         SetJetsName(const char *n)                 { fJetsName      = n;  }
  void         SetJetAlgo(Int_t a)                        { fJetAlgo          = a     ; }
  void         SetJetType(Int_t t)                        { fJetType          = t     ; }
  void         SetRecombScheme(Int_t scheme)              { fRecombScheme     = scheme; }
  void         SetMinJetArea(Double_t a)                  { fMinJetArea       = a     ; }
  void         SetMinJetPt(Double_t j)                    { fMinJetPt         = j     ; }
  void         SetRadius(Double_t r)                      { fRadius        = r;  }
  void         SetGhostArea(Double_t gharea)              { fGhostArea        = gharea; }
  void         SetJetEtaRange(Double_t emi, Double_t ema) { fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void         SetJetPhiRange(Double_t pmi, Double_t pma) { fJetPhiMin        = pmi   ; fJetPhiMax = pma; }

  void         SetMinJetTrackPt(Double_t min)             { fMinJetTrackPt = min;}
  void         SetMaxJetTrackPt(Double_t max)             { fMaxJetTrackPt = max;}
  void         SetJetTrackEtaRange(Double_t etmi, Double_t etma) { fJetTrackEtaMin = etmi; fJetTrackEtaMax = etma; }
  void         SetJetTrackPhiRange(Double_t ptmi, Double_t ptma) { fJetTrackPhiMax = ptmi; fJetTrackPhiMax = ptma; }
  void         SetJetTrackDCAcut(Double_t d)              { fJetTrackDCAcut   = d     ; }
  void         SetJetTracknHitsFit(Double_t h)            { fJetTracknHitsFit = h     ; }
  void         SetJetTracknHitsRatio(Double_t r)          { fJetTracknHitsRatio = r   ; }
  void         SetMinJetTowerE(Double_t min)              { mTowerEnergyMin = min;}
  void         SetJetTowerERange(Double_t enmi, Double_t enmx) { fJetTowerEMin = enmi; fJetTowerEMax = enmx; }
  void         SetJetTowerEtaRange(Double_t temi, Double_t temx) { fJetTowerEtaMin = temi; fJetTowerEtaMax = temx; }
  void         SetJetTowerPhiRange(Double_t tpmi, Double_t tpmx) { fJetTowerPhiMin = tpmi; fJetTowerPhiMax = tpmx; }
  void         SetMinJetClusPt(Double_t min)              { fMinJetClusPt  = min;}
  void         SetMinJetClusE(Double_t min)               { fMinJetClusE   = min;}

  void         SetLocked()                                { fLocked = kTRUE;}
  void         SetTrackEfficiency(Double_t t)             { fTrackEfficiency  = t     ; }
  void         SetLegacyMode(Bool_t mode)                 { fLegacyMode       = mode  ; }
  void         SetFillGhost(Bool_t b=kTRUE)               { fFillGhost        = b     ; }

  // for jet substructure routines
  StJetUtility*          AddUtility(StJetUtility* utility);
  TObjArray*             GetUtilities()                   { return fUtilities ; }

  // jets
  TClonesArray*          GetJets()                        { return fJets; }
  TClonesArray*          GetJetsBGsub()                   { return fJetsBGsub; }
  TClonesArray*          GetJetConstit()                  { return fJetsConstit; }
 
  // getters
  Double_t               GetGhostArea()                   { return fGhostArea         ; }
  const char*            GetJetsName()                    { return fJetsName.Data()   ; }
  const char*            GetJetsTag()                     { return fJetsTag.Data()    ; }
  Double_t               GetJetEtaMin()                   { return fJetEtaMin         ; }
  Double_t               GetJetEtaMax()                   { return fJetEtaMax         ; }
  Double_t               GetJetPhiMin()                   { return fJetPhiMin         ; }
  Double_t               GetJetPhiMax()                   { return fJetPhiMax         ; }
  UInt_t                 GetJetType()                     { return fJetType           ; }
  UInt_t                 GetJetAlgo()                     { return fJetAlgo           ; }
  Int_t                  GetRecombScheme()                { return fRecombScheme      ; }
  Bool_t                 GetLegacyMode()                  { return fLegacyMode        ; }
  Double_t               GetMinJetArea()                  { return fMinJetArea        ; }
  Double_t               GetMinJetPt()                    { return fMinJetPt          ; }
  Double_t               GetRadius()                      { return fRadius            ; }
  Double_t               GetJetTrackEtaMin()              { return fJetTrackEtaMin    ; }
  Double_t               GetJetTrackEtaMax()              { return fJetTrackEtaMax    ; }
  Double_t               GetJetTrackPhiMin()              { return fJetTrackPhiMin    ; }
  Double_t               GetJetTrackPhiMax()              { return fJetTrackPhiMax    ; }
  Double_t               GetTrackEfficiency()             { return fTrackEfficiency   ; }
  Double_t               GetJetTowerEtaMin()              { return fJetTowerEtaMin    ; }
  Double_t               GetJetTowerEtaMax()              { return fJetTowerEtaMax    ; }
  Double_t               GetJetTowerPhiMin()              { return fJetTowerPhiMin    ; }
  Double_t               GetJetTowerPhiMax()              { return fJetTowerPhiMax    ; }

  Bool_t                 IsLocked() const;

  // set hadronic correction fraction for matched tracks to towers
  void                   SetHadronicCorrFrac(float frac)    { mHadronicCorrFrac = frac; }

 protected:
  // this 1st version is deprecated as the parameters are global for the class and already set
  void                   FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius);
  void                   FindJets();
  void                   FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet>& constituents,
                            std::vector<fastjet::PseudoJet>& constituents_sub, Int_t flag = 0, TString particlesSubName = "");
  Bool_t                 AcceptJetTrack(StPicoTrack *trk, Float_t B, TVector3 Vert);// track accept cuts function
  Bool_t                 AcceptJetTower(StPicoBTowHit *tower, Int_t towerID);             // tower accept cuts function
  Int_t                  GetCentBin(Int_t cent, Int_t nBin) const;                        // centrality bin
  Bool_t                 SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut); // centrality bin to cut on for analysis
  Bool_t                 GetMomentum(TVector3 &mom, const StPicoBTowHit* tower, Double_t mass, Int_t towerID) const;
  Bool_t                 CheckForMB(int RunFlag, int type);
  Bool_t                 CheckForHT(int RunFlag, int type);
  Bool_t                 DoComparison(int myarr[], int elems);
  void                   FillEmcTriggersArr();
  Double_t               GetMaxTrackPt();
  Int_t                  FastJetBGsub();

  void                   SetSumw2(); // set errors weights 

  // may not need any of these except fill jet branch if I want 2 different functions
  void                   FillJetBranch();
  void                   FillJetBGBranch();
  void                   InitUtilities();
  void                   PrepareUtilities();
  void                   ExecuteUtilities(StJet* jet, Int_t ij);
  void                   TerminateUtilities();

  Bool_t                 GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  // switches
  Bool_t                 doWriteHistos;           // write QA histos
  Bool_t                 doUsePrimTracks;         // primary track switch
  Int_t                  fDebugLevel;             // debug printout level
  Int_t                  fRunFlag;                // Run Flag numerator value
  Bool_t                 doppAnalysis;            // use pp analysis data
  Int_t                  fCentralityDef;          // Centrality Definition enumerator value
  Bool_t                 fRequireCentSelection;   // require particular centrality bin
  Bool_t                 doConstituentSubtr;      // run constituent subtractor

  // event cuts
  Double_t               fEventZVtxMinCut;        // min event z-vertex cut
  Double_t               fEventZVtxMaxCut;        // max event z-vertex cut
  Int_t                  fCentralitySelectionCut; // centrality selection cut
  Bool_t                 doUseBBCCoincidenceRate; // use BBC or ZDC Coincidence Rate, kFALSE = ZDC
  Double_t               fMaxEventTrackPt;        // max track pt in the event (to cut on)    

  // event variables
  Double_t               Bfield;                  // event Bfield
  TVector3               mVertex;                 // event vertex 3-vector
  Double_t               zVtx;                    // z-vertex component

  // event selection types
  UInt_t                 fEmcTriggerEventType;    // Physics selection of event used for signal - HT or JP
  UInt_t                 fMBEventType;            // MB selection  
  UInt_t                 fTriggerToUse;           // trigger to use for analysis
  UInt_t                 fBadTowerListVers;       // version of bad tower file list to use
  Int_t                  fEmcTriggerArr[8];       // EMCal triggers array: used to select signal and do QA

  // tower to firing trigger type matched array
  Bool_t                 fTowerToTriggerTypeHT1[4801];// Tower with corresponding HT1 trigger type array
  Bool_t                 fTowerToTriggerTypeHT2[4801];// Tower with corresponding HT2 trigger type array
  Bool_t                 fTowerToTriggerTypeHT3[4801];// Tower with corresponding HT3 trigger type array

  // output file name string
  TString                mOutName;
  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection

  Int_t                  fJetAlgo;                // jet algorithm (kt, akt, etc)
  Int_t                  fJetType;                // jet type (full, charged, neutral)
  Int_t                  fRecombScheme;           // recombination scheme used by fastjet

  StFJWrapper            fjw;//!fastjet wrapper
  //StFJWrapper            fjwBG;//!fastjet wrapper for background

  // jet attributes
  Double_t               fRadius;                 // jet radius
  TString                fJetsTag;                // tag of jet collection (usually = "Jets")
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  Double_t               fJetPhiMin;              // minimum phi to keep jet in output
  Double_t               fJetPhiMax;              // maximum phi to keep jet in output
  Double_t               fJetEtaMin;              // minimum eta to keep jet in output
  Double_t               fJetEtaMax;              // maximum eta to keep jet in output
  Double_t               fGhostArea;              // ghost area

  // track attributes
  Double_t               fMinJetTrackPt;          // min jet track transverse momentum cut
  Double_t               fMaxJetTrackPt;          // max jet track transverse momentum cut
  Double_t               fMinJetClusPt;           // min jet cluster transverse momentum cut
  Double_t               fMinJetClusE;            // min jet cluster energy cut
  Double_t               fMinJetTowerE;           // min jet tower energy cut - not used (use mTowerEnergyMin)
  Double_t               fJetTrackEtaMin;         // min jet track eta cut
  Double_t               fJetTrackEtaMax;         // max jet track eta cut
  Double_t               fJetTrackPhiMin;         // min jet track phi cut
  Double_t               fJetTrackPhiMax;         // max jet track phi cut
  Double_t               fJetTrackDCAcut;         // max jet track dca cut
  Int_t                  fJetTracknHitsFit;       // requirement for track hits
  Double_t               fJetTracknHitsRatio;     // requirement for nHitsFit / nHitsMax
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)

  // tower attributes
  Double_t               fJetTowerEMin;           // min jet tower energy cut
  Double_t               fJetTowerEMax;           // max jet tower energy cut
  Double_t               fJetTowerEtaMin;         // min jet tower eta cut
  Double_t               fJetTowerEtaMax;         // max jet tower eta cut
  Double_t               fJetTowerPhiMin;         // min jet tower phi cut
  Double_t               fJetTowerPhiMax;         // max jet tower phi cut
  Double_t               mTowerEnergyMin;         // min jet tower energy cut
  Float_t                mHadronicCorrFrac;       // hadronic correction fraction from 0.0 to 1.0

  // may not need some of next bools
  TObjArray             *fUtilities;              // jet utilities (gen subtractor, constituent subtractor etc.)
  Bool_t                 fLocked;                 // true if lock is set
  Bool_t                 fIsInit;                 //!=true if already initialized
  Bool_t                 fLegacyMode;             //!=true to enable FJ 2.x behavior
  Bool_t                 fFillGhost;              //!=true ghost particles will be filled in StJet obj

  // jet and jet constituent objects
  TClonesArray          *fJets;                   //!jet collection
  TClonesArray          *fJetsBGsub;              //!jet background subtracted collection
  vector<fastjet::PseudoJet> fFull_Event;         //!jet input vectors
  vector<fastjet::PseudoJet> fConstituents;       //!jet constituents
  TClonesArray          *fJetsConstit;            //!jet constituents ClonesArray
  TClonesArray          *fJetsConstitBGsub;       //!jet constituents background subtracted ClonesArray  

  // fastjet definitions
  // might add..

  // TEST ---
  StEmcGeom       *mGeom;
  
  static const Int_t     fgkConstIndexShift;      //!contituent index shift

 private:
  StPicoDstMaker *mPicoDstMaker; // PicoDstMaker object
  StPicoDst      *mPicoDst;      // PicoDst object
  StPicoEvent    *mPicoEvent;    // PicoEvent object

  // position object
  StEmcPosition2 *mEmcPosition;

  // centrality objects
  StRefMultCorr* grefmultCorr;

  Float_t        mTowerMatchTrkIndex[4801];
  Bool_t         mTowerStatusArr[4801];

  // histograms
  TH1F           *fHistMultiplicity;//!
  TH1F           *fHistCentrality;//!
  TH1F           *fHistFJRho;//!

  TH1F           *fHistNTrackvsPt;//!
  TH1F           *fHistNTrackvsPhi;//!
  TH1F           *fHistNTrackvsEta;//!
  TH2F           *fHistNTrackvsPhivsEta;//!
  TH1F           *fHistNTowervsID;//!
  TH1F           *fHistNTowervsE;//!
  TH1F           *fHistNTowervsEt;//!
  TH1F           *fHistNTowervsPhi;//!
  TH1F           *fHistNTowervsEta;//!
  TH2F           *fHistNTowervsPhivsEta;//!

  TH1F           *fHistJetNTrackvsPt;//!
  TH1F           *fHistJetNTrackvsPhi;//!
  TH1F           *fHistJetNTrackvsEta;//!
  TH2F           *fHistJetNTrackvsPhivsEta;//!
  TH1F           *fHistJetNTowervsID;//!
  TH1F           *fHistJetNTowervsE;//!
  TH1F           *fHistJetNTowervsEt;//!
  TH1F           *fHistJetNTowervsPhi;//!
  TH1F           *fHistJetNTowervsEta;//!
  TH2F           *fHistJetNTowervsPhivsEta;//!

  TH1F           *fHistNJetsvsPt;//!
  TH1F           *fHistNJetsvsPhi;//!
  TH1F           *fHistNJetsvsEta;//!
  TH2F           *fHistNJetsvsPhivsEta;//!
  TH1F           *fHistNJetsvsArea;//!
  TH1F           *fHistNJetsvsNConstituents;//!
  TH1F           *fHistNJetsvsNTracks;//!
  TH1F           *fHistNJetsvsNTowers;//!

  TH2F           *fHistQATowIDvsEta;//!
  TH2F           *fHistQATowIDvsPhi;//!

  // bad and dead tower list functions and arrays
  Bool_t IsTowerOK( Int_t mTowId );
  Bool_t IsTowerDead( Int_t mTowId );
  std::set<Int_t> badTowers; 
  std::set<Int_t> deadTowers;

  // maker names
  //TString         fJetMakerName;

  StJetMakerTaskBGsub(const StJetMakerTaskBGsub&);            // not implemented
  StJetMakerTaskBGsub &operator=(const StJetMakerTaskBGsub&); // not implemented

  ClassDef(StJetMakerTaskBGsub, 1) // Jet producing task
};
#endif
