#ifndef STPICOTRACKCLUSTERQA_H
#define STPICOTRACKCLUSTERQA_H
// $Id$

// base class
//class StJetFrameworkPicoBase;
//class StJetMakerTask;

// includes
//#include "StJetFrameworkPicoBase.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"

// ROOT classes
class TClonesArray;
class TObjArray;
class TList;
class TH1;
class TH2;

// STAR classes
class StMaker;
class StChain;
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;

class StEmcGeom;
class StEmcCluster;
class StEmcCollection;
class StBemcTables; //v3.14
class StEmcPosition;

// centrality class
class StRefMultCorr;

// TEST for clusters TODO
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"

#include "StMyAnalysisMaker.h"

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

  // switches
  virtual void         SetUsePrimaryTracks(Bool_t P)    { doUsePrimTracks       = P; } 
  virtual void         SetDebugLevel(Int_t l)           { fDebugLevel           = l; }
  virtual void         SetRunFlag(Int_t f)              { fRunFlag              = f; }
  virtual void         SetCentralityDef(Int_t c)        { fCentralityDef        = c; }
  virtual void         SetTurnOnCentSelection(Bool_t o) { fRequireCentSelection = o; }
  virtual void         SetCentralityBinCut(Int_t c)     { fCentralitySelectionCut = c; }

  // event setters
  virtual void         SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }

  // track / cluster setters 
  virtual void         SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
  virtual void         SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
  virtual void         SetMinClusterPt(Double_t minpt)    { fClusterPtMinCut    = minpt;} // min cluster cut
  virtual void         SetMaxClusterPt(Double_t maxpt)    { fClusterPtMaxCut    = maxpt;} // max cluster cut
  virtual void         SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
  virtual void         SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
  virtual void         SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
  virtual void         SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
  virtual void         SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

  virtual void         SetTriggerEventType(UInt_t te)       { fTriggerEventType = te; }

  // efficiency correction setter
  virtual void         SetDoEffCorr(Int_t effcorr)          { fDoEffCorr = effcorr; }

  // common setters
  void                 SetClusName(const char *n)       { fCaloName      = n;  }
  void                 SetTracksName(const char *n)     { fTracksName    = n;  }
  void                 SetTrackEfficiency(Double_t t)   { fTrackEfficiency  = t     ; }

  Double_t             GetTrackEfficiency()             { return fTrackEfficiency   ; }

 protected:
  void                   RunQA(TObjArray *tracks, TObjArray *clus);
  Bool_t                 AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function
  Int_t                  GetCentBin(Int_t cent, Int_t nBin) const; // centrality bin
  Bool_t                 SelectAnalysisCentralityBin(Int_t centbin, Int_t fCentralitySelectionCut); // centrality bin to cut on for analysis
  TH1*                   FillEmcTriggersHist(TH1* h);                          // EmcTrigger counter histo
  TH1*                   FillEventTriggerQA(TH1* h);                           // filled event trigger QA plots
  Bool_t                 DoComparison(int myarr[], int elems);

  // switches
  Bool_t                 doWriteHistos;           // write QA histos
  Bool_t                 doUsePrimTracks;         // primary track switch
  Int_t                  fDebugLevel;             // debug printout level
  Int_t                  fRunFlag;                // Run Flag numerator value
  Int_t                  fCentralityDef;          // Centrality Definition enumerator value
  Bool_t                 fRequireCentSelection;   // require particular centrality bin
  Bool_t                 fDoEffCorr;              // efficiency correction to tracks

  // event cuts
  Double_t               fEventZVtxMinCut;        // min event z-vertex cut
  Double_t               fEventZVtxMaxCut;        // max event z-vertex cut
  Int_t                  fCentralitySelectionCut; // centrality selection cut

  // names
  TString                mOutName;                // name of output file
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
  UInt_t          fTriggerEventType;           // Physics selection of event used for signal
  Int_t           fEmcTriggerArr[7];           // EMCal triggers array: used to select signal and do QA

  StEmcGeom       *mGeom;
  StEmcCollection *mEmcCol; 
 
 private:
  StMuDst         *mu; // muDst object
  StPicoDstMaker  *mPicoDstMaker; // PicoDstMaker object
  StPicoDst       *mPicoDst; // PicoDst object
  StPicoEvent     *mPicoEvent; // PicoEvent object

  // centrality objects
  StRefMultCorr   *grefmultCorr;

  // histograms
  TH1F           *fHistNTrackvsPt;//!
  TH1F           *fHistNTrackvsPhi;//!
  TH1F           *fHistNTrackvsEta;//!
  TH2F           *fHistNTrackvsPhivsEta;//!
  TH1F           *fHistNTowervsE;//!
  TH1F           *fHistNTowervsPhi;//!
  TH1F           *fHistNTowervsEta;//!
  TH2F           *fHistNTowervsPhivsEta;//!

  // QA histos
  TH1  *fHistEventSelectionQA;//! 
  TH1  *fHistEventSelectionQAafterCuts;//!
  TH1  *hTriggerIds;//!
  TH1  *hEmcTriggers;//!


  StPicoTrackClusterQA(const StPicoTrackClusterQA&);            // not implemented
  StPicoTrackClusterQA &operator=(const StPicoTrackClusterQA&); // not implemented

  ClassDef(StPicoTrackClusterQA, 1) // track/cluster QA task
};
#endif
