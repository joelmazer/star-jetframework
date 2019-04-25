// $Id$
//
// Parameter class with a Clear.
//
// adapted from AliROOT class file AliRhoParameter.cxx for STAR
// original
// Author: C.Loizides

#include "StRhoParameter.h"

ClassImp(StRhoParameter)

//________________________________________________________________________
StRhoParameter::StRhoParameter() : 
  TParameter<Double_t>()
{
  // Dummy constructor.
}

//________________________________________________________________________
StRhoParameter::StRhoParameter(const char *name, Double_t val) :
  TParameter<Double_t>(name,val)
{
  // Constructor.
}

//________________________________________________________________________
void StRhoParameter::Clear(Option_t * /*option*/) 
{ 
  // Clear.
  SetVal(0);
}
