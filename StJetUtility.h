#ifndef STJETUTILITY_H
#define STJETUTILITY_H

// ROOT includes
#include <TNamed.h>

// my includes
#include "StFJWrapper.h"

// classes
class StJetMakerTask;
class StJet;
class StFJWrapper;

class StJetUtility : public TNamed
{
 public:
  // explicit added in the below 2 lines to remove cppcheck warning for c++11 compilers

  StJetUtility();
  explicit StJetUtility(const char* name);
  explicit StJetUtility(const StJetUtility &jet);
  StJetUtility& operator=(const StJetUtility &jet);
  ~StJetUtility() {;}

  virtual void Init() = 0;                                                 // Executed only once at the end of AliEmcalJetTask::DoInit()
  virtual void InitEvent(StFJWrapper& fjw) = 0;                            //
  virtual void Prepare(StFJWrapper& fjw) = 0;                              // Executed for each event at the beginning of AliEmcalJetTask::FillJetBranch()
  virtual void ProcessJet(StJet* jet, Int_t ij, StFJWrapper& fjw) = 0;     // Executed for each jet in the loop in AliEmcalJetTask::FillJetBranch()
  virtual void Terminate(StFJWrapper& fjw) = 0;                            // Executed for each event at the end of AliEmcalJetTask::FillJetBranch()

  void SetJetTask(StJetMakerTask* jetTask) { fJetTask = jetTask; }

 protected:

  StJetMakerTask       *fJetTask     ;// pointer to the main jet task
  Bool_t                fInit        ;//! whether or not the utility has been initialized

  ClassDef(StJetUtility, 2) // Abstract Emcal jet utility class
};
#endif
