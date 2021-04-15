//******************************************************************************
//                                                                            
// StEmcPosition2.cxx
//
// Authors: Joel Mazer
//
// Initial version: 2018/12/20
// Added to on: 1/22/19
//******************************************************************************

#include "StEmcPosition2.h"

// C++ includes
#include <math.h>
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

// StRoot includes
#include "StEmcUtil/geometry/StEmcGeom.h"
 
// StMuDstMaker:
#include "StMuDSTMaker/COMMON/StMuTrack.h"
 
ClassImp(StEmcPosition2)
 
//_______________________________________________________________________________________________
StEmcPosition2::StEmcPosition2():TObject()
{  
  mGeom[0] = StEmcGeom::getEmcGeom("bemc");  
  mGeom[1] = StEmcGeom::getEmcGeom("bprs");  
  mGeom[2] = StEmcGeom::getEmcGeom("bsmde");  
  mGeom[3] = StEmcGeom::getEmcGeom("bsmdp");  
}
//_______________________________________________________________________________________________
StEmcPosition2::~StEmcPosition2()
{
}
//_______________________________________________________________________________________________
Int_t StEmcPosition2::getTowerEtaPhi(const Double_t eta, const Double_t phi,
                                  Float_t *towerEta, Float_t *towerPhi ) const
{
  *towerEta = 0; *towerPhi = 0;
  Float_t tempTowerEta = 0, tempTowerPhi = 0;
  Int_t m = 0, e = 0, s = 0, towerId = -1;

  mGeom[0]->getBin(phi, eta, m, e, s);
  if (m == 0) return -1;
  if (s < 0) s = 1;
  mGeom[0]->getId(m, e, s, towerId);
  mGeom[0]->getEtaPhi(towerId, tempTowerEta, tempTowerPhi);
  *towerEta = tempTowerEta;
  *towerPhi = tempTowerPhi;
  return 0;
}
//_______________________________________________________________________________________________
Int_t StEmcPosition2::getNextTowerId(const Float_t eta, const Float_t phi, const Int_t nTowersdEta, const Int_t nTowersdPhi) const
{
  Int_t m,e,s;
  mGeom[0]->getBin( phi, eta, m, e, s );
  if(m > 0 && m <= 120) {
    if(s < 0) s = 1;
    return getNextTowerId(m,e,s,nTowersdEta,nTowersdPhi);
  }
  return 0;
}
//_______________________________________________________________________________________________
Int_t StEmcPosition2::getNextTowerId(const Int_t softId, const Int_t nTowersdEta, const Int_t nTowersdPhi) const
{
  if(softId < 1 || softId > 4800) return 0;
  Int_t m,e,s;
  mGeom[0]->getBin(softId,m,e,s);
  return getNextTowerId(m,e,s,nTowersdEta,nTowersdPhi);
}
//_______________________________________________________________________________________________
Int_t StEmcPosition2::getNextTowerId(const Int_t m, const Int_t e, const Int_t s, const Int_t nTowersdEta, const Int_t nTowersdPhi) const
{
  if(m< 1 || m > 120) return 0;
  if(e< 1 || e > 20) return 0;
  if(s< 1 || s > 2) return 0;
  return getNextId(1,m,e,s,nTowersdEta,nTowersdPhi);
}
//_______________________________________________________________________________________________
Int_t StEmcPosition2::getNextId(const Int_t det, const Int_t m, const Int_t e, const Int_t s, const Int_t nEta, const Int_t nPhi) const
{
  if(det < 1 || det > 4) return 0;
  if(m < 1 || m > 120) return 0;
  if(s < 1 || s > mGeom[det-1]->NSub()) return 0;
  if(e < 1 || e > mGeom[det-1]->NEta()) return 0;

  Int_t ef = e+nEta;
  Int_t sf = s+nPhi;
  Int_t mf = m;

  Int_t NE = mGeom[det-1]->NEta();
  Int_t NS = mGeom[det-1]->NSub();

  if(abs(ef) > NE) return 0;

  do {
    if(sf <= 0) {
      sf += NS;
      mf--;
      if(mf == 60) mf  = 120;
      if(mf == 0)  mf  = 60;
    }
    if(sf > NS) {
      sf -= NS;
      mf++;
      if(mf == 61)  mf = 1;
      if(mf == 121) mf = 61;
    }
  } while(sf <= 0 || sf > NS);

  if(ef <= 0) {
    ef = 1-ef;
    sf = NS-sf+1;
    if(ef > NE) return 0;
    Int_t rid,etmp,stmp;
    Float_t eta,phi;
    mGeom[det-1]->getId(mf, ef, sf, rid);
    mGeom[det-1]->getEtaPhi(rid, eta, phi);
    mGeom[det-1]->getBin(phi,-eta,mf,etmp,stmp);
  }

  Int_t rid;
  if(mf < 1 || mf > 120) return 0;
  if(ef < 1 || ef > NE)  return 0;
  if(sf < 1 || sf > NS)  return 0;
  mGeom[det-1]->getId(mf, ef, sf, rid);
  return rid;

}
//_______________________________________________________________________________________________
Int_t StEmcPosition2::getNextId(const Int_t det, const Int_t softId, const Int_t nEta, const Int_t nPhi)const
{
  if(det < 1 || det > 4) return 0;
  Int_t m,e,s;
  if(softId < 1)return 0;
  if((det == 1 || det == 2) && softId > 4800)return 0;
  if((det == 3 || det == 4) && softId > 18000)return 0;
  mGeom[det-1]->getBin(softId,m,e,s);
  if(m > 0 && m <= 120)
    {
      if(s < 0) s = 1;
      return getNextId(det,m,e,s,nEta,nPhi);
    }
  return 0;
}
//_______________________________________________________________________________________________
Float_t StEmcPosition2::getDistTowerToTrack( Double_t trackEta, Double_t trackPhi,
                                         Int_t nTowersdEta, Int_t nTowersdPhi ) const

{
  Int_t towerId = 0;
  Float_t towerEta = 0, towerToTrackdEta = 0;
  Float_t towerPhi = 0, towerToTrackdPhi = 0;
  Float_t mdistTowerToTrack = 0;

  towerId = getNextTowerId( trackEta, trackPhi, nTowersdEta, nTowersdPhi );
  if (towerId != 0) {
    // Getting eta and phi of neighbour tower
    mGeom[0]->getEtaPhi(towerId, towerEta, towerPhi);
    towerToTrackdEta = towerEta-trackEta;
    towerToTrackdPhi = towerPhi-trackPhi;

    mdistTowerToTrack = ::sqrt( ::pow(towerToTrackdEta, 2) + ::pow(towerToTrackdPhi, 2) );

    return mdistTowerToTrack;
  } else return -1;
}
//_______________________________________________________________________________________________
TVector3 StEmcPosition2::getPosFromVertex( const TVector3& position, int TowerId ) const
{
  TVector3 Zero(0,0,0);
  if(TowerId < 1 || TowerId > 4800) return Zero;
         
  float xTower, yTower, zTower;
  //TVector3 position = vertex->position(); //modified to work with StMuDst instead of StEvent
  mGeom[0]->getXYZ(TowerId, xTower, yTower, zTower);
  TVector3 towerPosition(xTower, yTower, zTower);
  TVector3 PositionFromVertex = towerPosition - position;
         
  return PositionFromVertex;
}
//_______________________________________________________________________________________________
Float_t StEmcPosition2::getThetaFromVertex( const TVector3& vertex, int TowerId ) const
{
  TVector3 p = getPosFromVertex(vertex, TowerId );
  return p.Theta();
}
//_______________________________________________________________________________________________
Float_t StEmcPosition2::getEtaFromVertex( const TVector3& vertex, int TowerId ) const
{
  TVector3 p = getPosFromVertex(vertex, TowerId );
  return p.PseudoRapidity();
}
//_______________________________________________________________________________________________
Float_t StEmcPosition2::getPhiFromVertex( const TVector3& vertex, int TowerId ) const
{
  TVector3 p = getPosFromVertex(vertex, TowerId );
  return p.Phi();
}
