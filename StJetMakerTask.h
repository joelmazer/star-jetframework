#ifndef STJETTASKNEW_H
#define STJETTASKNEW_H

// $Id$

#include "StMaker.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"

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

namespace fastjet {
  class PseudoJet;
}

/**
 * @class StJetTask
 * @brief General jet finder task implementing a wrapper for FastJet
 * 
 * after getting a bunch of the functionality working, have added some code 
 * directly from the ALICE version of AliEmcalJetTask written by:
 * @author Constantin Lozides <cloizides@lbl.gov>, Lawrence Berkeley National Laboratory
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * This class implements a wrapper for the FastJet jet finder. It allows
 * to set a jet definition (jet algorithm, recombination scheme) and the
 * list of jet constituents. The jet finding is delegated to
 * the class StFJWrapper which implements an interface to FastJet.
 *
 * The below is not funcitonal yet:
 * The FastJet contrib utilities are available via the StJetUtility base class
 * and its derived classes. Utilities can be added via the AddUtility(StJetUtility*) method.
 * All the utilities added in the list will be executed. Users can implement new utilities
 * deriving a new class from StJetUtility to interface functionalities of the FastJet contribs.
 */


//class StJetMakerTask : public StMyAnalysisMaker {
class StJetMakerTask : public StMaker {
 public:
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
  StJetMakerTask(const char *name, double mintrackPt);
  virtual ~StJetMakerTask();

  // needed class functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // switches
  void                    SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; } 

  // common setters
  void         SetAlgo(Int_t a)                 { fAlgo          = a;  }
  void         SetClusName(const char *n)       { fCaloName      = n;  }
  void         SetJetsName(const char *n)       { fJetsName      = n;  }
  void         SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min;}
  void         SetMinJetClusPt(Double_t min)    { fMinJetClusPt  = min;}
  void         SetRadius(Double_t r)            { fRadius        = r;  }
  void         SetTracksName(const char *n)     { fTracksName    = n;  }
  void         SetType(Int_t t)                 { fType          = t;  }

  void                   SetEtaRange(Double_t emi, Double_t ema);
  void                   SetMinJetClusE(Double_t min);
  void                   SetPhiRange(Double_t pmi, Double_t pma);

  void                   SetGhostArea(Double_t gharea)              { fGhostArea        = gharea; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) { fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) { fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
  void                   SetJetAlgo(Int_t a)                   { fJetAlgo          = a     ; }
  void                   SetJetType(Int_t t)                   { fJetType          = t     ; }
  void                   SetRecombScheme(Int_t scheme)      { fRecombScheme     = scheme; }
//  void                   SetJetAlgo(EJetAlgo_t a)                   { fJetAlgo          = a     ; }
//  void                   SetJetType(EJetType_t t)                   { fJetType          = t     ; }
//  void                   SetRecombScheme(ERecoScheme_t scheme)      { fRecombScheme     = scheme; }
  void                   SetMinJetArea(Double_t a)                  { fMinJetArea       = a     ; }
  void                   SetMinJetPt(Double_t j)                    { fMinJetPt         = j     ; }
  void                   SetLocked()                                { fLocked = kTRUE;}
  void                   SetTrackEfficiency(Double_t t)             { fTrackEfficiency  = t     ; }
  void                   SetLegacyMode(Bool_t mode)                 { fLegacyMode       = mode  ; }
  void                   SetFillGhost(Bool_t b=kTRUE)               { fFillGhost        = b     ; }

// ========
/*
  void                   SetGhostArea(Double_t gharea)              { if (IsLocked()) return; fGhostArea        = gharea; }
  void                   SetJetsName(const char *n)                 { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetsTag(const char *n)                  { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) { if (IsLocked()) return; fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) { if (IsLocked()) return; fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
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
  Double_t               GetTrackEfficiency()             { return fTrackEfficiency   ; }


  void                   FillJetConstituents(StJet *jet, std::vector<fastjet::PseudoJet>& constituents,
                                             std::vector<fastjet::PseudoJet>& constituents_sub, Int_t flag = 0, TString particlesSubName = "");

  UInt_t FindJetAcceptanceType(Double_t eta, Double_t phi, Double_t r);

  Bool_t                 IsLocked() const;
  //void SetType(Int_t t);

/*
#if !defined(__CINT__) && !defined(__MAKECINT__)
  static FJJetAlgo       ConvertToFJAlgo(EJetAlgo_t algo);
  static FJRecoScheme    ConvertToFJRecoScheme(ERecoScheme_t reco);
#endif
*/

 protected:
  void FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius);
  //Int_t FindJets();

  // switches
  Bool_t                 doUsePrimTracks;         // primary track switch

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  Int_t                  fAlgo;                   // algo (0==kt, 1==antikt)
  Int_t                  fType;                   // jet type (0=all, 1=ch, 2=neutral)
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

//  StFJWrapper           fFastJetWrapper;         //!fastjet wrapper
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
  Double_t               fMinJetTrackPt;          // min jet track momentum
  Double_t               fMinJetClusPt;           // min jet cluster momentum
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)

  // may not need some of next bools
  TObjArray             *fUtilities;              // jet utilities (gen subtractor, constituent subtractor etc.)
  Bool_t                 fLocked;                 // true if lock is set
  Bool_t                 fIsInit;                 //!=true if already initialized
  Bool_t                 fLegacyMode;             //!=true to enable FJ 2.x behavior
  Bool_t                 fFillGhost;              //!=true ghost particles will be filled in StJet obj

  static const Int_t fgkConstIndexShift; //!contituent index shift

  TClonesArray          *fJets;                   //!jet collection

 private:
  StPicoDstMaker *mPicoDstMaker; // PicoDstMaker object
  StPicoDst      *mPicoDst; // PicoDst object
  StPicoEvent    *mPicoEvent; // PicoEvent object


  StJetMakerTask(const StJetMakerTask&);            // not implemented
  StJetMakerTask &operator=(const StJetMakerTask&); // not implemented

  ClassDef(StJetMakerTask, 1) // Jet producing task
};
#endif
