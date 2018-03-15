#ifndef __STJETPICOTOWER_H
#define __STJETPICOTOWER_H

#include <TObject.h>
#include <TMath.h>
#include <TArrayI.h>

class StJetPicoTower : public TObject
{
 public:
  StJetPicoTower();
  StJetPicoTower(StJetPicoTower &t);

  virtual ~StJetPicoTower();

  void    Clear(Option_t *Option = "");

  Int_t   GetId()     const {return fId;}
  Float_t GetEnergy() const {return fEnergy;}
  Float_t GetEta()    const {return fEta;}
  Float_t GetPhi()    const {return fPhi;}
  Int_t   GetADC()    const {return fADC;}
  Float_t GetEtaCorrected()    const {return fEtaCorrected;}
  Float_t GetEt()     const {return fEnergy / TMath::CosH(fEtaCorrected);}
  Float_t GetPhiCorrected()    const {return fPhiCorrected;}
  
  Int_t   GetSMDClusterP()   const {return fSMDClusterP;}
  Int_t   GetSMDClusterE()   const {return fSMDClusterE;}

  Int_t   GetTowerStatus()   const {return fTowerStatus;}
  Int_t   GetNAssocTracks() const {return fNAssocTracks;}

  Int_t   GetMatchedTrackIndex(Int_t idx) {return fMatchedTracks.At(idx);}

  const TArrayI *GetMatchedTracks() {return &fMatchedTracks;}
  const TArrayI *GetMatchedTrackIndexes() {return &fMatchedTracks;}
  
  // SETTERS
  void SetId(Int_t val)       {fId = val;}
  void SetEnergy(Float_t val) {fEnergy = val;}
  void SetEta(Float_t val)    {fEta = val;}
  void SetPhi(Float_t val)    {fPhi = val;}
  void SetADC(Int_t val)      {fADC = val;}
  void SetEtaCorrected(Float_t val)    {fEtaCorrected = val;}
  void SetPhiCorrected(Float_t val)    {fPhiCorrected = val;}
  
  void SetSMDClusterP(Int_t val)     {fSMDClusterP = val;}
  void SetSMDClusterE(Int_t val)     {fSMDClusterE = val;}
  
  void SetTowerStatus(Int_t val)     {fTowerStatus = val;}

  void AddMatchedTrack(Int_t idx);
  void AddMatchedTrackIndex(Int_t idx) {AddMatchedTrack(idx);}

  void SetNAssocTracks(Int_t val) {fNAssocTracks = val;}

 protected:

 private:
  Int_t           fSMDClusterP; // was SMDClustP[1111];   //[nCand]
  Int_t           fSMDClusterE; // was SMDClustE[1111];   //[nCand]
  Int_t           fTowerStatus; // wasTwrStatusMatched[1111];   //[nCand]

  Int_t           fId;           // was TwrAllId[4641];   //[nFlag]
  Float_t         fEnergy;       // was TwrAllEnergy[4641];   //[nFlag]
  Float_t         fEta;          // was TwrAllEta[4641];   //[nFlag]
  Float_t         fPhi;          // was TwrAllPhi[4641];   //[nFlag]
  Int_t           fADC;          // was TwrAllAdc[4641];   //[nFlag]
  Float_t         fEtaCorrected; // was TwrAllEtaCor[4641];   //[nFlag]
  Float_t         fPhiCorrected; // was TwrAllPhiCor[4641];   //[nFlag]
  Int_t           fNAssocTracks; //n of assoc tracks

  TArrayI         fMatchedTracks; //fMatchedTracks

  ClassDef(StJetPicoTower, 1)
};

#endif
