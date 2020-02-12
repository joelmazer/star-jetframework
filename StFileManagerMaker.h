#ifndef StFileManagerMaker_h
#define StFileManagerMaker_h

#include "StChain/StMaker.h"

class TClonesArray;
class TChain;
class TFile;
class StMuDst;
class StPicoDst;

class StFileManagerMaker : public StMaker
{
public:
  enum PicoIoMode {IoWrite=1, IoRead=2};

  // 'explicit' added below to constructor to remove cppcheck warning
  explicit StFileManagerMaker(char const* name = "PicoDst");
  StFileManagerMaker(PicoIoMode ioMode, char const* fileName = "", char const* name = "PicoDst");
  virtual ~StFileManagerMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void  Clear(Option_t* option = "");
  virtual Int_t Finish();

  /// Returns null pointer if no StPicoDst
  StPicoDst* picoDst();
  /// In read mode, returns pointer to the chain of .picoDst.root files
  TChain* chain();

private:

  void openWrite();
  void write();
  void closeWrite();
  Int_t openRead();
  void  read();
  void closeRead();

  /// A pointer to the main input source containing all muDst `TObjArray`s
  /// filled from corresponding muDst branches
  StMuDst*  mMuDst;

  /// A pointer to the main input/outpur picoDst structure containing all `TObjArray`s
  StPicoDst*  mPicoDst;

  TString   mInputFileName;        //! *.list - MuDst or picoDst
  TString   mOutputFileName;       //! FileName
  TFile*    mOutputFile;

  ClassDef(StFileManagerMaker, 0)
};
#endif
