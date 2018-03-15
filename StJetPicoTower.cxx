
#include "StJetPicoTower.h"
#include "TString.h"

#include "StJetPicoDefinitions.h"

ClassImp(StJetPicoTower)

#define __PICO_MAX_MATCHED_TRACKS 1

StJetPicoTower::StJetPicoTower()
  : TObject()
  , fSMDClusterP(0) 
  , fSMDClusterE(0) 
  , fTowerStatus(0)
  , fId(0)           
  , fEnergy(0)       
  , fEta(0)          
  , fPhi(0)          
  , fADC(0)          
  , fEtaCorrected(0) 
  , fPhiCorrected(0) 
  , fNAssocTracks(0)
  , fMatchedTracks(__PICO_MAX_MATCHED_TRACKS)
{
  //
  // Default constructor
  // All indexes of the matched tracks are set to -1
  //
  fMatchedTracks.Reset(-1);
}

StJetPicoTower::StJetPicoTower(StJetPicoTower &t)
  : TObject(t)
  , fSMDClusterP(t.fSMDClusterP) 
  , fSMDClusterE(t.fSMDClusterE) 
  , fTowerStatus(t.fTowerStatus)
  , fId(t.fId)           
  , fEnergy(t.fEnergy)       
  , fEta(t.fEta)          
  , fPhi(t.fPhi)          
  , fADC(t.fADC)          
  , fEtaCorrected(t.fEtaCorrected) 
  , fPhiCorrected(t.fPhiCorrected)   
  , fNAssocTracks(t.fNAssocTracks)
  , fMatchedTracks(t.fMatchedTracks)
{
  //
  // Copy constructor
  //
  ;
}

StJetPicoTower::~StJetPicoTower()
{
  //
  // Destructor
  // In fact it is safe not to clear at this point
  //
  Clear();
}

void StJetPicoTower::Clear(Option_t */*Option*/)
{
  //
  // Clear everything
  // 
  fSMDClusterP = 0; 
  fSMDClusterE = 0;
  fTowerStatus = 0;
  fId = 0;           
  fEnergy = 0;       
  fEta = 0;          
  fPhi = 0;          
  fADC = 0;          
  fEtaCorrected = 0; 
  fPhiCorrected = 0; 
  fNAssocTracks = 0;

  fMatchedTracks.Set(__PICO_MAX_MATCHED_TRACKS);
  fMatchedTracks.Reset(-1);
}

void StJetPicoTower::AddMatchedTrack(Int_t idx)
{
  //
  // Add new index
  // Resize the array if necessary
  //
  if (fNAssocTracks < fMatchedTracks.GetSize())
    {
      fMatchedTracks[fNAssocTracks++] = idx;      
    }
  else
    {
      //Int_t newsize = fMatchedTracks.GetSize() * 2;
      //since the N of matched tracks should drop steeply we can be modest
      //Int_t newsize = fMatchedTracks.GetSize() + 2;
      Int_t newsize = fMatchedTracks.GetSize() + 1;
      __DEBUG(3, Form("Resizing indexes array from %d to %d", fMatchedTracks.GetSize(), newsize));
      fMatchedTracks.Set(newsize);
      AddMatchedTrack(idx);
    }
}
