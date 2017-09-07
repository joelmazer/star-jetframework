#ifndef STNAMEDARRAYI_H
#define STNAMEDARRAYI_H

// $Id: StNamedArrayI.h  $

#include <TArrayI.h>
#include <TNamed.h>

class StNamedArrayI : public TNamed, public TArrayI {
 public: 
  StNamedArrayI();
  StNamedArrayI(const char *name, Int_t n);
  StNamedArrayI(const char *name, Int_t n, const Int_t* array);
  StNamedArrayI(const char *name, const TArrayI& array);

  void Clear(Option_t *option="");

private:
  StNamedArrayI(const StNamedArrayI&);             // not implemented
  StNamedArrayI& operator=(const StNamedArrayI&);  // not implemented
  
  ClassDef(StNamedArrayI, 1); // Named integer array
};
#endif
