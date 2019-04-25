#ifndef StRHOPARAMETER_H
#define StRHOPARAMETER_H

// $Id$
// Adapated from the AliROOT similar function: AliRhoParameter.h

class TString;
class TF1;

#include <TParameter.h>

class StRhoParameter : public TParameter<Double_t> {
 public: 
  StRhoParameter();
  StRhoParameter(const char *name, Double_t val=0);
  void        Clear(Option_t *option="");

private:
  StRhoParameter(const StRhoParameter&);             // not implemented
  StRhoParameter& operator=(const StRhoParameter&);  // not implemented
  
  ClassDef(StRhoParameter, 2); // Rho parameter
};
#endif
