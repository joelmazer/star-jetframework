#ifndef STRHO_H
#define STRHO_H

// $Id$
// adapted from the AliROOT class AliAnalysisTaskRho.h

#include "StRhoBase.h"

// additional includes
#include "StMaker.h"

// ROOT classes
class TH2;
class TH2F;
class TClonesArray;

// STAR classes
class StMaker;

class StRho : public StRhoBase {

 public:
  StRho();
  StRho(const char *name, Bool_t histo=kFALSE, const char* outName="", const char* jetMakerName="");
  //virtual ~StRho() {}
  virtual ~StRho();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt="");
  virtual Int_t Finish();

  // booking of histograms (optional)
  void    DeclareHistograms();
  void    WriteHistograms();

  void    SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

 protected:
  UInt_t            fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  TClonesArray     *fJets;//!jet collection

 private:
  TH2F             *fHistMultvsRho;//!

  StRho(const StRho&);             // not implemented
  StRho& operator=(const StRho&);  // not implemented
  
  ClassDef(StRho, 2); // Rho task
};
#endif
