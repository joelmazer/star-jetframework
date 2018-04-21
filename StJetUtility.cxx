#include "StJetUtility.h"

ClassImp(StJetUtility)

//______________________________________________________________________________
StJetUtility::StJetUtility() :
  TNamed(),
  fJetTask(0),
  fInit(kFALSE)
{
  // Dummy constructor.
}

//______________________________________________________________________________
StJetUtility::StJetUtility(const char* name) :
  TNamed(name, name),
  fJetTask(0),
  fInit(kFALSE)
{
  // Default constructor.
}

//______________________________________________________________________________
StJetUtility::StJetUtility(const StJetUtility &other) :
  TNamed(other),
  fJetTask(other.fJetTask),
  fInit(other.fInit)
{
  // Copy constructor.
}

//______________________________________________________________________________
StJetUtility& StJetUtility::operator=(const StJetUtility &other)
{
  // Assignment
  if (&other == this) return *this;
  TNamed::operator=(other);
  fJetTask = other.fJetTask;
  fInit = other.fInit;
  return *this;
}
