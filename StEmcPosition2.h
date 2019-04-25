//******************************************************************************
//                                                                            
// StEmcPosition2.h
//
// Authors: Joel Mazer
//
// Initial version: 12/15/18
// Added to on: 1/22/19
//******************************************************************************

#ifndef StEmcPosition2_H
#define StEmcPosition2_H

// ROOT includes 
#include "TObject.h"
#include "TVector3.h"

// StRoot classes
class StMuTrack;
class StEmcGeom;

class StEmcPosition2 : public TObject
{
   public:            
 
             StEmcPosition2();
    virtual  ~StEmcPosition2();

    Int_t             getTowerEtaPhi(Double_t eta, Double_t phi, Float_t* towerEta, Float_t* towerPhi) const; ///< Return tower eta/phi

    Int_t             getNextTowerId(Float_t eta, Float_t phi, Int_t nTowersdEta, Int_t nTowersdPhi) const; ///< Return neighbor tower id's
    Int_t             getNextTowerId(Int_t softId, Int_t nTowersdEta, Int_t nTowersdPhi) const; ///< Return neighbor tower id's
    Int_t             getNextTowerId(Int_t m, Int_t e, Int_t s, Int_t nTowersdEta, Int_t nTowersdPhi) const; ///< Return neighbor tower id's
    Int_t             getNextId(Int_t det, Int_t m, Int_t e, Int_t s, Int_t nEta, Int_t nPhi) const; ///< Return neighbor id (works for all detectors 1=bemc, 2=bprs, 3=bsmde, 4=bsmdp)
    Int_t             getNextId(Int_t det, Int_t softId, Int_t nEta, Int_t nPhi)const;///< Return neighbor id (works for all detectors 1=bemc, 2=bprs, 3=bsmde, 4=bsmdp)
    Float_t           getDistTowerToTrack(Double_t trackEta, Double_t trackPhi, Int_t nTowersdEta, Int_t nTowersdPhi) const; ///< Return distance from track to center of one tower

    TVector3          getPosFromVertex(const TVector3& position, Int_t TowerId) const; ///< Return Position from collision vertex
    Float_t           getThetaFromVertex(const TVector3& vertex, Int_t TowerId) const; ///< Return theta of the tower considering the collision vertex
    Float_t           getEtaFromVertex(const TVector3& vertex, Int_t TowerId) const; ///< Return eta of the tower considering the collision vertex
    Float_t           getPhiFromVertex(const TVector3& vertex, Int_t TowerId) const; ///< Return phi of the tower considering the collision vertex

   protected:     
 
     StEmcGeom *mGeom[4];   
 
   ClassDef(StEmcPosition2, 2)
};
#endif
