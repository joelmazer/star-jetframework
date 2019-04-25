#ifndef STRHOSPARSE_H
#define STRHOSPARSE_H

// $Id$
class TH2F;

#include "StRhoBase.h"

class StRhoSparse : public StRhoBase {

 public:
  StRhoSparse();
  StRhoSparse(const char *name, Bool_t histo=kFALSE, const char* outName="", const char* jetMakerName="");
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

  // set names of makers for global use
  virtual void     SetOutputFileName(const char *on)         { mOutName = on; }
  virtual void     SetJetMakerName(const char *jn)           { fJetMakerName = jn; }
  virtual void     SetJetBGMakerName(const char *bjn)        { fJetBGMakerName = bjn; }
  virtual void     SetRhoSparseMakerName(const char *rpn)    { fRhoSparseMakerName = rpn; }

 protected:
  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation
  Bool_t           fCreateHisto;                   // switch to create histograms
  Bool_t           fRhoCMS;                        // flag to run CMS method

 private:
  TH2F            *fHistOccCorrvsCent;//!      occupancy correction vs. centrality
  TH2F            *fHistOccCorrvsMult;//!      occupancy correction vs. multiplicity
  TH2F            *fHistMultvsUnCorrRho;//!    multiplicity vs. rho
  TH2F            *fHistMultvsCorrRho;//!      multiplicity vs. corrected rho 

  // maker names
  TString          mOutName;
  TString          fRhoSparseMakerName;

  StRhoSparse(const StRhoSparse&);             // not implemented
  StRhoSparse& operator=(const StRhoSparse&);  // not implemented
  
  ClassDef(StRhoSparse, 2); // Rho task
};
#endif
