#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

#include "TRegexp.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjectSet.h"

#include "StChain/StChain.h"
#include "StChain/StChainOpt.h"
#include "St_base/StMessMgr.h"
#include "StarRoot/TAttr.h"

#include "StFileManagerMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"


//_____________________________________________________________________________
StFileManagerMaker::StFileManagerMaker(char const* name) : StMaker(name),
  mMuDst(nullptr), mPicoDst(new StPicoDst()),
  mInputFileName(), mOutputFileName(), mOutputFile(nullptr)
{
}
//_____________________________________________________________________________
StFileManagerMaker::StFileManagerMaker(PicoIoMode ioMode, char const* fileName, char const* name) : StFileManagerMaker(name)
{
  StMaker::m_Mode = ioMode;
  //mInputFileName = fileName;
  mOutputFileName = fileName;
}
//_____________________________________________________________________________
StFileManagerMaker::~StFileManagerMaker()
{
}
//_____________________________________________________________________________
Int_t StFileManagerMaker::Init()
{
  switch (StMaker::m_Mode)
  {
    case PicoIoMode::IoWrite:
      // may not need this chunk..
      if (mInputFileName.Length() == 0) {
        // No input file
        mOutputFileName = GetChainOpt()->GetFileOut();
        mOutputFileName.ReplaceAll(".root", ".picoDst.root");
      }
      else // have given an Input File Name
      {
        mInputFileName = mInputFileName(mInputFileName.Index("st_"), mInputFileName.Length());
        mOutputFileName = mInputFileName;
        mOutputFileName.ReplaceAll("MuDst.root", "picoDst.root");
      }
      openWrite();
      break;

    case PicoIoMode::IoRead:
      openRead();
      break;

    default:
      LOG_ERROR << "Pico IO mode is not set ... " << endm;
      return kStErr;
  }

  return kStOK;
}

//_____________________________________________________________________________
Int_t StFileManagerMaker::Finish()
{
  if (StMaker::m_Mode == PicoIoMode::IoRead)
  {
    closeRead();
  }
  else if (StMaker::m_Mode == PicoIoMode::IoWrite)
  {
    closeWrite();
  }
  return kStOK;
}
//_____________________________________________________________________________
void StFileManagerMaker::openWrite()
{
  mOutputFile = new TFile(mOutputFileName.Data(), "RECREATE");
  LOG_INFO << " Output file: " << mOutputFileName.Data() << " created." << endm;
}
//_____________________________________________________________________________
//_____________________________________________________________________________
void StFileManagerMaker::Clear(char const*)
{
  if (StMaker::m_Mode == PicoIoMode::IoRead)
    return;
}
//_____________________________________________________________________________
void StFileManagerMaker::closeWrite()
{
  if (StMaker::m_Mode == PicoIoMode::IoWrite)
  {
    if (mOutputFile)
    {
      mOutputFile->Write();
      mOutputFile->Close();
    }
  }
}
//_____________________________________________________________________________
int StFileManagerMaker::Make()
{
  int returnStarCode = kStOK;
  return returnStarCode;
}
