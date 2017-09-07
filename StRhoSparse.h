#ifndef STRHOSPARSE_H
#define STRHOSPARSE_H

// $Id$
class TH2F;

#include "StRhoBase.h"

class StRhoSparse : public StRhoBase {

 public:
  StRhoSparse();
  StRhoSparse(const char *name, Bool_t histo=kFALSE, const char* outName="");
  virtual ~StRhoSparse();// {}

  // class required functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }
  void             SetCreateHistos(Bool_t cr)      { fCreateHisto = cr     ; }
  void             SetRhoCMS(Bool_t cms)           { fRhoCMS = cms ; }
  Bool_t           IsJetOverlapping(StJet* jet1, StJet* jet2);
  Bool_t           IsJetSignal(StJet* jet1);

 protected:
  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation
  Bool_t           fCreateHisto;                   // switch to create histograms
  Bool_t           fRhoCMS;                        // flag to run CMS method

  TH2F            *fHistOccCorrvsCent;             //!occupancy correction vs. centrality

  TString          mOutName;

  StRhoSparse(const StRhoSparse&);             // not implemented
  StRhoSparse& operator=(const StRhoSparse&);  // not implemented
  
  ClassDef(StRhoSparse, 1); // Rho task
};
#endif
