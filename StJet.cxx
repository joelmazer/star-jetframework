/* ***********************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// this class is adapted from the class AliEmcalJet
//
#include "StJet.h"
#include "StVParticle.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"

#include "Riostream.h"

//#include "AliParticleContainer.h"
//#include "AliClusterContainer.h"

/// \cond CLASSIMP
ClassImp(StJet);
/// \endcond

/**
 * Default constructor
 */
StJet::StJet() :
  StVParticle(),//FIXME
  fPt(0),
  fEta(0),
  fPhi(0),
  fM(0),
  fNEF(0),
  fArea(0),
  fAreaEta(0),
  fAreaPhi(0),
  fAreaE(0),
  fMaxCPt(0),
  fMaxNPt(0),
  fMaxTrackPt(0),
  fMCPt(0),
  fNn(0),
  fNch(0),
  fClusterIDs(),
  fTrackIDs(),
  fMatched(2),
  fMatchingType(0),
  fPtSub(0),
  fPtSubVect(0),
  fTriggers(0),
  fLabel(-1),
  fHasGhost(kFALSE),
  fGhosts(),
  fJetConstit()
{
  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;

}

/**
 * Constructor that uses the 3-momentum to define the jet axis.
 * It assumes zero mass for the jet.
 * @param px First transverse component of the jet momentum
 * @param px Second transverse component of the jet momentum
 * @param pz Longitudinal component of the jet momentum
 */
StJet::StJet(Double_t px, Double_t py, Double_t pz) :
  StVParticle(),//FIXME
  fPt(TMath::Sqrt(px * px + py* py)),
  fEta(TMath::ASinH(pz / fPt)),
  fPhi(0),
  fM(0),
  fNEF(0),
  fArea(0),
  fAreaEta(0),
  fAreaPhi(0),
  fAreaE(0),
  fMaxCPt(0),
  fMaxNPt(0),
  fMaxTrackPt(0),
  fMCPt(0),
  fNn(0),
  fNch(0),
  fClusterIDs(),
  fTrackIDs(),
  fMatched(2),
  fMatchingType(0),
  fPtSub(0),
  fPtSubVect(0),
  fTriggers(0),
  fLabel(-1),
  fHasGhost(kFALSE),
  fGhosts(),
  fJetConstit()
{
  if (fPt != 0) {
    fPhi = TVector2::Phi_0_2pi(TMath::ATan2(py, px));
  }

  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;

}

/**
 * Constructor that uses the 4-momentum to define the jet axis.
 * Coordinates are given in cylindrical system plus the mass.
 * @param pt Transverse component of the jet momentum
 * @param eta Pseudo-rapidity of the jet
 * @param phi Azimuthal angle of the jet axis
 * @param m Mass of the jet
 */
StJet::StJet(Double_t pt, Double_t eta, Double_t phi, Double_t m) :
  StVParticle(),//FIXME
  fPt(pt),
  fEta(eta),
  fPhi(phi),
  fM(m),
  fNEF(0),
  fArea(0),
  fAreaEta(0),
  fAreaPhi(0),
  fAreaE(0),
  fMaxCPt(0),
  fMaxNPt(0),
  fMaxTrackPt(0),
  fMCPt(0),
  fNn(0),
  fNch(0),
  fClusterIDs(),
  fTrackIDs(),
  fMatched(2),
  fMatchingType(0),
  fPtSub(0),
  fPtSubVect(0),
  fTriggers(0),
  fLabel(-1),
  fHasGhost(kFALSE),
  fGhosts(),
  fJetConstit()
{
  fPhi = TVector2::Phi_0_2pi(fPhi);

  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;

}

/**
 * Copy constructor.
 * @param jet Constant reference to copy from
 */
StJet::StJet(const StJet& jet) :
  StVParticle(jet),//FIXME
  fPt(jet.fPt),
  fEta(jet.fEta),
  fPhi(jet.fPhi),
  fM(jet.fM),
  fNEF(jet.fNEF),
  fArea(jet.fArea),
  fAreaEta(jet.fAreaEta),
  fAreaPhi(jet.fAreaPhi),
  fAreaE(jet.fAreaE),
  fMaxCPt(jet.fMaxCPt),
  fMaxNPt(jet.fMaxNPt),
  fMaxTrackPt(jet.fMaxTrackPt),
  fMCPt(jet.fMCPt),
  fNn(jet.fNn),
  fNch(jet.fNch),
  fClusterIDs(jet.fClusterIDs),
  fTrackIDs(jet.fTrackIDs),
  fMatched(jet.fMatched),
  fMatchingType(jet.fMatchingType),
  fPtSub(jet.fPtSub),
  fPtSubVect(jet.fPtSubVect),
  fTriggers(jet.fTriggers),
  fLabel(jet.fLabel),
  fHasGhost(jet.fHasGhost),
  fGhosts(jet.fGhosts),
  fJetConstit(jet.fJetConstit)
{
  // Copy constructor.
  fClosestJets[0]     = jet.fClosestJets[0];
  fClosestJets[1]     = jet.fClosestJets[1];
  fClosestJetsDist[0] = jet.fClosestJetsDist[0];
  fClosestJetsDist[1] = jet.fClosestJetsDist[1];

}

/**
 * Destructor.
 */
StJet::~StJet()
{
}

/**
 * Assignment operator
 * @param jet Constant reference to copy from
 * @return A reference to this
 */
StJet& StJet::operator=(const StJet& jet)
{
  if (this != &jet) {
    StVParticle::operator=(jet); //FIXME
    fPt                 = jet.fPt;
    fEta                = jet.fEta;
    fPhi                = jet.fPhi;
    fM                  = jet.fM;
    fNEF                = jet.fNEF;
    fArea               = jet.fArea;
    fAreaEta            = jet.fAreaEta;
    fAreaPhi            = jet.fAreaPhi;
    fAreaE              = jet.fAreaE;
    fMaxCPt             = jet.fMaxCPt;
    fMaxNPt             = jet.fMaxNPt;
    fMaxTrackPt         = jet.fMaxTrackPt;
    fMCPt               = jet.fMCPt;
    fNn                 = jet.fNn;
    fNch                = jet.fNch;
    fClusterIDs         = jet.fClusterIDs;
    fTrackIDs           = jet.fTrackIDs;
    fClosestJets[0]     = jet.fClosestJets[0];
    fClosestJets[1]     = jet.fClosestJets[1];
    fClosestJetsDist[0] = jet.fClosestJetsDist[0];
    fClosestJetsDist[1] = jet.fClosestJetsDist[1];
    fMatched            = jet.fMatched;
    fPtSub              = jet.fPtSub;
    fPtSubVect          = jet.fPtSubVect;
    fTriggers           = jet.fTriggers;
    fLabel              = jet.fLabel;
    fHasGhost = jet.fHasGhost;
    fGhosts   = jet.fGhosts;
    fJetConstit = jet.fJetConstit;
  }

  return *this;
}

/**
 * Compares two instances of StJet, ordering them based on their transverse momentum.
 * @param obj Pointer to another instance of StJet
 * @return -1 if this is smaller than obj, 1 if this is larger than obj, 0 if objects are equal or if obj is NULL or not an instance of StJet
 */
Int_t StJet::Compare(const TObject* obj) const
{
  //Return -1 if this is smaller than obj, 0 if objects are equal and 1 if this is larger than obj.

  if (obj == this) return 0;

  const StJet* jet = dynamic_cast<const StJet*>(obj);
  if (!jet) return 0;

  if (Pt() > jet->Pt()) return -1;
  else if (Pt() < jet->Pt()) return 1;
  else return 0;
}

/**
 * Builds a 4-momentum object using information contained in this instance of StJet
 * @param[out] vec Reference to TLorentzVector where the 4-momentum is returned
 */
void StJet::GetMomentum(TLorentzVector& vec) const
{
  vec.SetPtEtaPhiE(fPt, fEta, fPhi, E());
}

/**
 * Calculates transverse momentum after scalar average background subtraction.
 * The result can be negative. It saves the result if requested.
 * @param rho Scalar average background
 * @param save If kTRUE, stores the result in a class field
 * @return The subtracted transverse momentum
 */
Double_t StJet::PtSub(Double_t rho, Bool_t save)
{
  Double_t ptcorr = fPt - rho * fArea;
  if (save) fPtSub = ptcorr;
  return ptcorr;
}

/**
 * Calculates transverse momentum after vectorial average background subtraction.
 * The result can be negative. It saves the result if requested.
 * @param rho Vectorial average background
 * @param save If kTRUE, stores the result in a class field
 * @return The subtracted transverse momentum
 */
Double_t StJet::PtSubVect(Double_t rho, Bool_t save)
{
  Double_t dx = Px() - rho * fArea * TMath::Cos(fAreaPhi);
  Double_t dy = Py() - rho * fArea * TMath::Sin(fAreaPhi);
  //Double_t dz = Pz() - rho * fArea * TMath::SinH(fAreaEta);
  Double_t ptcorr = TMath::Sqrt(dx * dx + dy * dy);
  if (save) fPtSubVect = ptcorr;
  return ptcorr;
}

/**
 * Calculates 4-momentum after vectorial average background subtraction.
 * It saves the result if requested.
 * @param rho Vectorial average background
 * @param save If kTRUE, stores the result in a class field
 * @return The subtracted 4-momentum
 */
TLorentzVector StJet::SubtractRhoVect(Double_t rho, Bool_t save)
{
  TLorentzVector vecCorr;
  GetMomentum(vecCorr);
  TLorentzVector vecBg;
  vecBg.SetPtEtaPhiE(fArea, fAreaEta, fAreaPhi, fAreaE);
  vecBg *= rho;
  vecCorr -= vecBg;
  if (save) {
    Double_t dPhi = TMath::Abs(TVector2::Phi_mpi_pi(Phi() - vecCorr.Phi()));
    Int_t signum = dPhi <= TMath::PiOver2() ? 1 : -1;
    fPtSubVect = signum * vecCorr.Pt();
  }
  return vecCorr;
}

/**
 *  Sort constituent by index (increasing).
 *
 */
void StJet::SortConstituents()
{
  std::sort(fClusterIDs.GetArray(), fClusterIDs.GetArray() + fClusterIDs.GetSize());
  std::sort(fTrackIDs.GetArray(), fTrackIDs.GetArray() + fTrackIDs.GetSize());
}

/**
 * Helper function to calculate the distance between two jets or a jet and a particle
 * @param part Constant pointer to another particle
 * @return Distance in the eta-phi phase space
 */
Double_t StJet::DeltaR(const StVParticle* part) const // FIXME
{
  Double_t dPhi = Phi() - part->Phi();
  Double_t dEta = Eta() - part->Eta();
  dPhi = TVector2::Phi_mpi_pi(dPhi);
  return TMath::Sqrt(dPhi * dPhi + dEta * dEta);
}

/**
 * Sorting jet constituents by pT (decreasing)  
 * It returns a standard vector with the indexes of the constituents relative to fTrackIDs.
 * To retrieve the track do:
 * ~~~{.cxx}
 * TClonesArray* fTracksContArray = jetCont->GetParticleContainer()->GetArray();
 * std::vector< int > index_sorted_list = jet->GetPtSortedTrackConstituentIndexes(fTracksContArray);
 * for (std::size_t i = 0; i < jet->GetNumberOfTracks(); i++ ) {
 * track = jet->TrackAt ( index_sorted_list.at (i), fTracksContArray );
 * // use track;
 * }
 * ~~~
 * @param tracks Array containing pointers to the tracks from which jet constituents are drawn
 * @return Standard vector with the list of constituent indexes (relative to fTrackIDs)
 */

/*
std::vector<int> StJet::GetPtSortedTrackConstituentIndexes(TClonesArray* tracks) const
{
  typedef std::pair<Double_t, Int_t> ptidx_pair;

  // Create vector for Pt sorting
  std::vector<ptidx_pair> pair_list;

  for (Int_t i_entry = 0; i_entry < GetNumberOfTracks(); i_entry++) {
    StVParticle* track = TrackAt(i_entry, tracks); // FIXME
    if (!track) {
      Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)", i_entry, tracks->GetName(), TrackAt(i_entry), tracks->GetEntriesFast()); // FIXME
      continue;
    }
    pair_list.push_back(std::make_pair(track->Pt(), i_entry));
  }

  std::stable_sort(pair_list.begin() , pair_list.end() , sort_descend());

  // return a vector of indexes of constituents (sorted descending by pt)
  std::vector <int> index_sorted_list;

  // populating the return object with indexes of sorted tracks
  for (auto it : pair_list) index_sorted_list.push_back(it.second);

  return index_sorted_list;
}
*/

/**
 * Get the momentum fraction of a jet constituent
 * @param trkPx First transverse component of the momentum of the jet constituent
 * @param trkPy Second transverse component of the momentum of the jet constituent
 * @param trkPz Longitudinal component of the momentum of the jet constituent
 * @return Momentum fraction
 */
Double_t StJet::GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz) const
{
  Double_t pJetSq = P();
  pJetSq *= pJetSq;

  if (pJetSq > 1e-6) {
    return (trkPx * Px() + trkPy * Py() + trkPz * Pz()) / pJetSq;
  }
  else {
    Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f", GetName(), pJetSq); // FIXME
    return -1;
  }
}

/**
 * Get the momentum fraction of a jet constituent
 * @param trk Jet constituent
 * @return Momentum fraction
 */
Double_t StJet::GetZ(const StVParticle* trk) const // FIXME
{
  return GetZ(trk->Px(), trk->Py(), trk->Pz());
}

/**
 * Get Xi = Log(1 / z) of constituent track
 * @param trk Pointer to a constituent track
 * @return Xi of the constituent
 */
Double_t StJet::GetXi(const StVParticle* trk) const
{
  return TMath::Log(1 / GetZ(trk));
}

/**
 * Get Xi = Log(1 / z) of constituent track
 * @param trkPx First transverse component of the momentum of the jet constituent
 * @param trkPy Second transverse component of the momentum of the jet constituent
 * @param trkPz Longitudinal component of the momentum of the jet constituent
 * @return Xi of the constituent
 */
Double_t StJet::GetXi( const Double_t trkPx, const Double_t trkPy, const Double_t trkPz ) const
{
  return TMath::Log(1 / GetZ(trkPx, trkPy, trkPz));
}

/**
 * Find the leading track constituent of the jet.
 * @param tracks Array containing the pointers to the tracks from which jet constituents are drawn
 * @return Pointer to the leading track of the jet
 */

/*
StVParticle* StJet::GetLeadingTrack(TClonesArray* tracks) const //FIXME
{
  StVParticle* maxTrack = 0; // FIXME
  for (Int_t i = 0; i < GetNumberOfTracks(); i++) {
    StVParticle* track = TrackAt(i, tracks); // FIXME
    if (!track) {
      AliError(Form("Unable to find jet track %d in collection %s (pos in collection %d, max %d)",
          i, tracks->GetName(), TrackAt(i), tracks->GetEntriesFast()));
      continue;
    }
    if (!maxTrack || track->Pt() > maxTrack->Pt())
      maxTrack = track;
  }

  return maxTrack;
}
*/

/**
 * Reset jet matching information.
 */
void StJet::ResetMatching()
{
  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 999;
  fClosestJetsDist[1] = 999;
  fMatched = 2;
}

/**
 * Checks whether a certain track is among the jet constituents by looking for its index
 * @param Index of the track to search
 * @return The position of the track in the jet constituent array, if the track is found; -1 if the track is not a jet constituent
 */
Int_t StJet::ContainsTrack(Int_t it) const
{
  for (Int_t i = 0; i < fTrackIDs.GetSize(); i++) {
    if (it == fTrackIDs[i]) return i;
  }
  return -1;
}

/**
 * Create string representation of the Jet. The string respresentation contains
 * - \f$ p_{t} \f$ of the jet
 * - \f$ \eta \f$ of the jet
 * - \f$ \phi \f$ of the jet
 * - \f$ p_{t} \f$ of the maximum charged and neutral particle
 * - Number of constituent tracks
 * - Number of constituent clusters
 * - Jet area
 * - Neutral energy fraction
 * @return String representation of the jet
 */
TString StJet::toString() const {
  return TString::Format("Jet pT = %.2f, eta = %.2f, phi = %.2f, max charged pT = %.2f, max neutral pT = %.2f, N tracks = %d, N clusters = %d, Area = %.2f, NEF = %.2f",
         Pt(), Eta(), Phi(), MaxChargedPt(), MaxNeutralPt(), GetNumberOfTracks(), GetNumberOfClusters(), Area(), NEF());
}

/**
 * Print basic jet information using the string representation provided by
 * StJet::toString
 *
 * @param unused
 */
void StJet::Print(Option_t* /*opt*/) const
{
  Printf("%s\n", toString().Data());
}

/**
 * Print basic jet information on an output stream using the string representation provided by
 * StJet::toString. Used by operator<<
 * @param in output stream stream
 * @return reference to the output stream
 */
std::ostream &StJet::Print(std::ostream &in) const {
  in << toString().Data();
  return in;
}

/**
 * Prints the list of constituents in the standard output
 * @param tracks Array containing the pointers to tracks
 * @param clusters Array containing the pointers to the clusters
 */

/*
void StJet::PrintConstituents(TClonesArray* tracks, TClonesArray* clusters) const
{
  if (tracks) {
    for (Int_t i = 0; i < GetNumberOfTracks(); i++) {
      StVParticle* part = TrackAt(i, tracks); // FIXME
      if (part) {
        Printf("Track %d (index = %d) pT = %.2f, eta = %.2f, phi = %.2f, PDG code = %d", i, TrackAt(i), part->Pt(), part->Eta(), part->Phi(), part->PdgCode());
      }
    }
  }
}
*/


/**
 * Implementation of the output stream operator for StJet. Printing
 * basic jet information provided by function toString
 * @param in output stream
 * @param myjet Jet which will be printed
 * @return Reference to the output stream
 */
std::ostream &operator<<(std::ostream &in, const StJet &myjet) {
  std::ostream &result = myjet.Print(in);
  return result;
}

/**
 * Add a ghost particle to the ghost particle array.
 * This function should be called by the jet finder to fill the ghost particle array.
 * @param dPx First component of the transverse momentum of the particle
 * @param dPy First component of the transverse momentum of the particle
 * @param dPz Longitudinal component of the momentum of the particle
 * @param dE Energy of the particle
 */
void StJet::AddGhost(const Double_t dPx, const Double_t dPy, const Double_t dPz, const Double_t dE)
{
  TLorentzVector ghost(dPx, dPy, dPz, dE);
  fGhosts.push_back(ghost);
  if (!fHasGhost) fHasGhost = kTRUE;
  return;
}

/**
 * Clear this object: remove matching information, jet constituents, ghosts
 */
void StJet::Clear(Option_t */*option*/)
{
  fClusterIDs.Set(0);
  fTrackIDs.Set(0);
  fClosestJets[0] = 0;
  fClosestJets[1] = 0;
  fClosestJetsDist[0] = 0;
  fClosestJetsDist[1] = 0;
  fMatched = 0;
  fPtSub = 0;
  fGhosts.clear();
  fHasGhost = kFALSE;
  fJetConstit.clear();
}

/**
 * Retrieve the track constituent corresponding to the index found at a certain position.
 * Automatically retrieves the particle from the proper TClonesArray. This function is preferred to
 * TrackAt(Int_t, TClonesArray), which is only kept for backwards compatibility.
 *
 * @param idx Position of the track constituent
 * @return Pointer to the track constituent requested (if found)
 */

/* ==============
StVParticle* StJet::Track(Int_t idx) const // FIXME
{
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;

  return StParticleContainer::GetEmcalContainerIndexMap().GetObjectFromGlobalIndex(TrackAt(idx)); // FIXME
}
================ */

/**
 * Finds the track constituent corresponding to the index found at a certain position.
 * @param idx Position of the track constituent
 * @param ta Array with pointers to the tracks from which jet constituents are drawn
 * @return Pointer to the track constituent requested (if found)
 */

/* ======================
StVParticle* StJet::TrackAt(Int_t idx, TClonesArray *ta) const // FIXME
{
  if (!ta) return 0;
  auto res =  StParticleContainer::GetEmcalContainerIndexMap().LocalIndexFromGlobalIndex(TrackAt(idx)); // FIXME
  if (res.second != ta) {
    //AliWarning(Form("TClonesArray %s that was passed does not correspond to the passed index! The index belongs to a different TClonesArray named %s! Returning the object corresponding to the index (not the passed TClonesArray)! Consider fixing by updating to jet->Track(index).", ta->GetName(), res.second->GetName())); // FIXME
  }
  return dynamic_cast<StVParticle*>(res.second->At(res.first)); // FIXME
}
========================*/

/**
 * Checks whether a given track is among the jet constituents
 * @param track Pointer to the track to be searched
 * @param tracks Array with pointers to the tracks from which jet constituents are drawn
 * @return Position of the track among the jet constituents, if the track is found; -1 otherwise
 */
Int_t StJet::ContainsTrack(StVParticle* track, TClonesArray* tracks) const // FIXME
{
  if (!tracks || !track) return 0;
  return ContainsTrack(tracks->IndexOf(track));
}


// test functions for now = FIXME
/**
 * Checks whether a cluster is among the constituents of a jet.
 * @param jet Pointer to an StJet object
 * @param iclus Index of the cluster to look for
 * @param sorted If the constituents are sorted by index it will speed up computation
 * @return kTRUE if the cluster is among the constituents, kFALSE otherwise
 */
Bool_t StJet::IsJetCluster(StJet* jet, Int_t iclus, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfClusters(); ++i) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

/**
 * Checks whether a track is among the constituents of a jet.
 * @param jet Pointer to an StJet object
 * @param itrack Index of the track to look for
 * @param sorted If the constituents are sorted by index it will speed up computation
 * @return kTRUE if the track is among the constituents, kFALSE otherwise
 */
Bool_t StJet::IsJetTrack(StJet* jet, Int_t itrack, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfTracks(); ++i) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

// TODO TEST
void StJet::AddJetConstit(const Double_t dPx, const Double_t dPy, const Double_t dPz, const Double_t dE)
{
  TLorentzVector constit(dPx, dPy, dPz, dE);
  //fJetConstit.push_back(constit);
  return;
}
