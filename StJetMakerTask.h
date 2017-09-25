#ifndef STJETTASKNEW_H
#define STJETTASKNEW_H

// $Id$

#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"

// TEST for clusters TODO
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
class StEmcGeom;
class StEmcCluster;
class StEmcCollection;
class StBemcTables; //v3.14

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH2;

// STAR classes
class StPicoDst;
class StPicoDstMaker;

class StFJWrapper;
class StJetUtility;

// STAR includes
#include "StFJWrapper.h"
#include "FJ_includes.h"
#include "StJet.h"
#include "StMyAnalysisMaker.h"

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "StIndexMap.h"
#endif

namespace fastjet {
  class PseudoJet;
}

/**
 * @brief General jet finder task implementing a wrapper for FastJet
 * 
 * after getting a bunch of the functionality working, have added some code 
 * directly from the ALICE version of AliEmcalJetTask written by:
 * 
 * @author Constantin Lozides <cloizides@lbl.gov>, Lawrence Berkeley National Laboratory
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * This class implements a wrapper for the FastJet jet finder. It allows
 * to set a jet definition (jet algorithm, recombination scheme) and the
 * list of jet constituents. The jet finding is delegated to
 * the class StFJWrapper which implements an interface to FastJet.
 *
 * The below is not functional yet:
 * The FastJet contrib utilities are available via the StJetUtility base class
 * and its derived classes. Utilities can be added via the AddUtility(StJetUtility*) method.
 * All the utilities added in the list will be executed. Users can implement new utilities
 * deriving a new class from StJetUtility to interface functionalities of the FastJet contribs.
 */


//class StJetMakerTask : public StMyAnalysisMaker {
class StJetMakerTask : public StMaker {
 public:

  // jet type enumerator
  enum EJetType_t {
    kFullJet,
    kChargedJet,
    kNeutralJet
  };

/*
  typedef StMyAnalysisMaker::EJetType_t EJetType_t;
  typedef StMyAnalysisMaker::EJetAlgo_t EJetAlgo_t;
  typedef StMyAnalysisMaker::ERecoScheme_t ERecoScheme_t;

#if !defined(__CINT__) && !defined(__MAKECINT__)
  typedef fastjet::JetAlgorithm FJJetAlgo;
  typedef fastjet::RecombinationScheme FJRecoScheme;
#endif
*/

  StJetMakerTask();
  StJetMakerTask(const char *name, double mintrackPt, bool dohistos, const char* outName);
  virtual ~StJetMakerTask();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  // switches
  void         SetUsePrimaryTracks(Bool_t P)    { doUsePrimTracks   = P; } 
  void         SetDebugLevel(Int_t l)           { fDebugLevel    = l;  }

  // common setters
  void         SetClusName(const char *n)       { fCaloName      = n;  }
  void         SetTracksName(const char *n)     { fTracksName    = n;  }
  void         SetJetsName(const char *n)       { fJetsName      = n;  }
  void         SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min;}
  void         SetMaxJetTrackPt(Double_t max)   { fMaxJetTrackPt = max;}
  void         SetMinJetClusPt(Double_t min)    { fMinJetClusPt  = min;}
  void         SetRadius(Double_t r)            { fRadius        = r;  }

  void         SetMinJetClusE(Double_t min);
  void         SetGhostArea(Double_t gharea)              { fGhostArea        = gharea; }
  void         SetJetEtaRange(Double_t emi, Double_t ema) { fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void         SetJetPhiRange(Double_t pmi, Double_t pma) { fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
  void         SetJetTrackEtaRange(Double_t etmi, Double_t etma) { fJetTrackEtaMin = etmi; fJetTrackEtaMax = etma; }
  void         SetJetTrackPhiRange(Double_t ptmi, Double_t ptma) { fJetTrackPhiMax = ptmi; fJetTrackPhiMax = ptma; }
  void         SetJetTracknHitsFit(Double_t h)            { fJetTracknHitsFit = h     ; }
  void         SetJetTracknHitsRatio(Double_t r)          { fJetTracknHitsRatio = r   ; }
  void         SetJetAlgo(Int_t a)                        { fJetAlgo          = a     ; }
  void         SetJetType(Int_t t)                        { fJetType          = t     ; }
  void         SetRecombScheme(Int_t scheme)              { fRecombScheme     = scheme; }
//  void                   SetJetAlgo(EJetAlgo_t a)                   { fJetAlgo          = a     ; }
//  void                   SetJetType(EJetType_t t)                   { fJetType          = t     ; }
//  void                   SetRecombScheme(ERecoScheme_t scheme)      { fRecombScheme     = scheme; }
  void         SetMinJetArea(Double_t a)                  { fMinJetArea       = a     ; }
  void         SetMinJetPt(Double_t j)                    { fMinJetPt         = j     ; }
  void         SetLocked()                                { fLocked = kTRUE;}
  void         SetTrackEfficiency(Double_t t)             { fTrackEfficiency  = t     ; }
  void         SetLegacyMode(Bool_t mode)                 { fLegacyMode       = mode  ; }
  void         SetFillGhost(Bool_t b=kTRUE)               { fFillGhost        = b     ; }

// ========
/*
  void                   SetGhostArea(Double_t gharea)              { if (IsLocked()) return; fGhostArea        = gharea; }
  void                   SetJetsName(const char *n)                 { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetsTag(const char *n)                  { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) { if (IsLocked()) return; fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) { if (IsLocked()) return; fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
  void                   SetJetTrackEtaRange(Double_t etmi, Double_t etma) { if (IsLocked()) return; fJetTrackEtaMin = etmi; fJetTrackEtaMax = etma; }
  void                   SetJetTrackPhiRange(Double_t ptmi, Double_t ptma) { if (IsLocked()) return; fJetTrackPhiMax = ptmi; fJetTrackPhiMax = ptma; }
  void                   SetJetAlgo(EJetAlgo_t a)                   { if (IsLocked()) return; fJetAlgo          = a     ; }
  void                   SetJetType(EJetType_t t)                   { if (IsLocked()) return; fJetType          = t     ; }
  void                   SetLocked()                                { fLocked = kTRUE;}
  void                   SetMinJetArea(Double_t a)                  { if (IsLocked()) return; fMinJetArea       = a     ; }
  void                   SetMinJetPt(Double_t j)                    { if (IsLocked()) return; fMinJetPt         = j     ; }
  void                   SetRecombScheme(ERecoScheme_t scheme)      { if (IsLocked()) return; fRecombScheme     = scheme; }
  void                   SetTrackEfficiency(Double_t t)             { if (IsLocked()) return; fTrackEfficiency  = t     ; }
  void                   SetLegacyMode(Bool_t mode)                 { if (IsLocked()) return; fLegacyMode       = mode  ; }
  void                   SetFillGhost(Bool_t b=kTRUE)               { if (IsLocked()) return; fFillGhost        = b     ; }
//  void                   SetRadius(Double_t r)                      { if (IsLocked()) return; fRadius           = r     ; }
*/
// =======

  // for jet substructure routines
  StJetUtility* AddUtility(StJetUtility* utility);
  TObjArray* GetUtilities() { return fUtilities ; }

  // jets
  TClonesArray* GetJets()                        { return fJets; }
 
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

  void                   FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet>& constituents,
                                             std::vector<fastjet::PseudoJet>& constituents_sub, Int_t flag = 0, TString particlesSubName = "");

  Bool_t                 IsLocked() const;

/*
#if !defined(__CINT__) && !defined(__MAKECINT__)
  static FJJetAlgo       ConvertToFJAlgo(EJetAlgo_t algo);
  static FJRecoScheme    ConvertToFJRecoScheme(ERecoScheme_t reco);
#endif
*/

 protected:
  void FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius);
  //Int_t FindJets();
  Bool_t                 AcceptJetTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function

  // switches
  Bool_t                 doWriteHistos;           // write QA histos
  Bool_t                 doUsePrimTracks;         // primary track switch
  Int_t                  fDebugLevel;             // debug printout level

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  Double_t               fRadius;                 // jet radius

  // need to tweak type of next 3
/*
  EJetType_t             fJetType;                // jet type (full, charged, neutral)
  EJetAlgo_t             fJetAlgo;                // jet algorithm (kt, akt, etc)
  ERecoScheme_t          fRecombScheme;           // recombination scheme used by fastjet
*/
  Int_t             fJetAlgo;                // jet algorithm (kt, akt, etc)
  Int_t             fJetType;                // jet type (full, charged, neutral)
  Int_t             fRecombScheme;           // recombination scheme used by fastjet

  StFJWrapper            fjw; //!fastjet wrapper

  // may not need any of these except fill jet branch if I want 2 different functions
  void                   FillJetBranch();
  void                   InitUtilities();
  void                   PrepareUtilities();
  void                   ExecuteUtilities(StJet* jet, Int_t ij); void                   TerminateUtilities();

  Bool_t                 GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  TString                fJetsTag;                // tag of jet collection (usually = "Jets")
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  Double_t               fJetPhiMin;              // minimum phi to keep jet in output
  Double_t               fJetPhiMax;              // maximum phi to keep jet in output
  Double_t               fJetEtaMin;              // minimum eta to keep jet in output
  Double_t               fJetEtaMax;              // maximum eta to keep jet in output
  Double_t               fGhostArea;              // ghost area
  Double_t               fMinJetTrackPt;          // min jet track transverse momentum
  Double_t               fMaxJetTrackPt;          // max jet track transverse momentum
  Double_t               fMinJetClusPt;           // min jet cluster transverse momentum
  Double_t               fJetTrackEtaMin;         // min jet track eta
  Double_t               fJetTrackEtaMax;         // max jet track eta
  Double_t               fJetTrackPhiMin;         // min jet track phi
  Double_t               fJetTrackPhiMax;         // max jet track phi
  Int_t                  fJetTracknHitsFit;       // requirement for track hits
  Double_t               fJetTracknHitsRatio;     // requirement for nHitsFit / nHitsMax
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)

  // may not need some of next bools
  TObjArray             *fUtilities;              // jet utilities (gen subtractor, constituent subtractor etc.)
  Bool_t                 fLocked;                 // true if lock is set
  Bool_t                 fIsInit;                 //!=true if already initialized
  Bool_t                 fLegacyMode;             //!=true to enable FJ 2.x behavior
  Bool_t                 fFillGhost;              //!=true ghost particles will be filled in StJet obj

  TClonesArray          *fJets;                   //!jet collection
  vector<fastjet::PseudoJet> fConstituents;       //!jet constituents

  // TEST ---
  StEmcGeom       *mGeom;
  StEmcCollection *mEmcCol;
  
  static const Int_t     fgkConstIndexShift;      //!contituent index shift

#if !(defined(__CINT__) || defined(__MAKECINT__))
// Handle mapping between index and containers
//StIndexMap <TClonesArray, StVCluster> fClusterContainerIndexMap;    //!<! Mapping between index and cluster containers  // FIXME
StIndexMap <TClonesArray, StVParticle> fParticleContainerIndexMap; //!<! Mapping between index and particle containers
#endif

 private:
  StMuDst        *mu; // muDst object
  StPicoDstMaker *mPicoDstMaker; // PicoDstMaker object
  StPicoDst      *mPicoDst; // PicoDst object
  StPicoEvent    *mPicoEvent; // PicoEvent object

  // histograms
  TH1F           *fHistJetNTrackvsPt;//!
  TH1F           *fHistJetNTrackvsPhi;//!
  TH1F           *fHistJetNTrackvsEta;//!
  TH1F           *fHistJetNTowervsE;//!
  TH1F           *fHistJetNTowervsPhi;//!
  TH1F           *fHistJetNTowervsEta;//!
  TH2F           *fHistJetNTowervsPhivsEta;//!

  // output file name string
  TString         mOutName;
     
  // maker names
  //TString         fJetMakerName;

  StJetMakerTask(const StJetMakerTask&);            // not implemented
  StJetMakerTask &operator=(const StJetMakerTask&); // not implemented

  ClassDef(StJetMakerTask, 1) // Jet producing task
};
#endif
