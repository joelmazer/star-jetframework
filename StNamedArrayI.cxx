// $Id: StNamedArrayI.cxx  $
//
// Named integer array.
//
// Author: S.Aiola

#include "StNamedArrayI.h"

ClassImp(StNamedArrayI)

//________________________________________________________________________
StNamedArrayI::StNamedArrayI() : 
  TNamed("StNamedArrayI","StNamedArrayI"),
  TArrayI()
{
  // Dummy constructor.

}

//________________________________________________________________________
StNamedArrayI::StNamedArrayI(const char *name, Int_t n) :
  TNamed(name,name),
  TArrayI(n)
{
  // Standard constructor.
  Clear();
}

//________________________________________________________________________
StNamedArrayI::StNamedArrayI(const char *name, Int_t n, const Int_t* array) :
  TNamed(name,name),
  TArrayI(n, array)
{
  // TArrayI copy c-style array constructor.

}

//________________________________________________________________________
StNamedArrayI::StNamedArrayI(const char *name, const TArrayI& array) :
  TNamed(name,name),
  TArrayI(array)
{
  // TArrayI copy constructor.
  
}

//________________________________________________________________________
void StNamedArrayI::Clear(Option_t * /*option*/) 
{ 
  // Clear.

  Reset(-1);
}
