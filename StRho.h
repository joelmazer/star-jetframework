#ifndef STRHO_H
#define STRHO_H

// $Id$
// adapted from the AliROOT class AliAnalysisTaskRho.h

class TH2;
class TH2F;
class TClonesArray;

#include "StRhoBase.h"

// some includes
#include "StMaker.h"

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StMaker;

//class StRho : virtual public StMaker, virtual public StRhoBase { //FIXME
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

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

 protected:
  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  TClonesArray          *fJets;//!jet collection

 private:
  // PicoDstMaker and PicoDst object pointer
  StPicoDstMaker *mPicoDstMaker;
  StPicoDst      *mPicoDst;
  StPicoEvent    *mPicoEvent;
  StJetMakerTask   *JetMaker;

  TH2F             *fHistMultvsRho;//!

  // maker names
  TString          mOutName;
  TString          fJetMakerName;
  TString          fRhoMakerName;

  StRho(const StRho&);             // not implemented
  StRho& operator=(const StRho&);  // not implemented
  
  ClassDef(StRho, 1); // Rho task
};
#endif
