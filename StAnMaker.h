#ifndef StAnMaker_h
#define StAnMaker_h

#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StJetFrameworkPicoBase.h"
class StJetFrameworkPicoBase;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class THnSparse;
class TProfile;
class TString;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;

// my STAR classes
class StJetMakerTask;
class StJet;
class StRho;
class StRhoParameter;

//class StAnMaker : public StMaker {
class StAnMaker : public StJetFrameworkPicoBase {
  public:

    StAnMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char *jetMakerName, const char *rhoMakerName);
    virtual ~StAnMaker();
   
    // class required functions
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    // booking of histograms (optional)
    void    DeclareHistograms();
    void    WriteHistograms();

    // switches
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetDebugLevel(Int_t l)             { fDebugLevel       = l; }
    virtual void            SetPrintEventCounter(Bool_t c)     { doPrintEventCounter = c; }
    virtual void            SetRunFlag(Int_t f)                { fRunFlag          = f; }
    virtual void            SetCentralityDef(Int_t c)          { fCentralityDef    = c; }

    virtual void            SetMinJetPt(Double_t j)            { fMinPtJet         = j; }    // min jet pt
    virtual void            SetJetConstituentCut(Double_t mc)  { fJetConstituentCut= mc;}    // min constituent pt cut
    virtual void            SetJetMaxTrackPt(Double_t t)       { fTrackBias        = t; }    // track bias
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    
    // event setters
    virtual void            SetEventZVtxRange(Double_t zmi, Double_t zma) { fEventZVtxMinCut = zmi; fEventZVtxMaxCut = zma; }

    // track setters
    virtual void            SetMinTrackPt(Double_t minpt)      { fTrackPtMinCut    = minpt;} // min track cut
    virtual void            SetMaxTrackPt(Double_t maxpt)      { fTrackPtMaxCut    = maxpt;} // max track cut
    virtual void            SetTrackPhiRange(Double_t ptmi, Double_t ptma) { fTrackPhiMinCut = ptmi; fTrackPhiMaxCut = ptma; }
    virtual void            SetTrackEtaRange(Double_t etmi, Double_t etma) { fTrackEtaMinCut = etmi; fTrackEtaMaxCut = etma; }
    virtual void            SetTrackDCAcut(Double_t d)         { fTrackDCAcut = d       ; }
    virtual void            SetTracknHitsFit(Double_t h)       { fTracknHitsFit = h     ; }
    virtual void            SetTracknHitsRatio(Double_t r)     { fTracknHitsRatio = r   ; }

    // mixed selection - setters
    virtual void            SetTriggerEventType(UInt_t te)       { fTriggerEventType = te; }

    // efficiency correction setter
    virtual void            SetDoEffCorr(Int_t effcorr)          { fDoEffCorr = effcorr; }

    // use rho to correct jet pt in correlation sparses
    virtual void            SetCorrectJetPt(Bool_t cpt)          { fCorrJetPt = cpt; }

    // leading jets
    StJet*                  GetLeadingJet(StRhoParameter* eventRho = 0x0);
    StJet*                  GetSubLeadingJet(StRhoParameter* eventRho = 0x0);

  protected:
    void                    RunTracks();
    void                    RunTowers();
    void                    RunJets();
    Int_t                   GetCentBin(Int_t cent, Int_t nBin) const;             // centrality bin
    Double_t                RelativePhi(Double_t mphi, Double_t vphi) const;      // relative jet track angle
    Double_t                RelativeEPJET(Double_t jetAng, Double_t EPAng) const; // relative jet event plane angle
    void                    FillEmcTriggers();                          // EmcTrigger counter histo
    void                    FillEventTriggerQA();                           // filled event trigger QA plots
    Bool_t                  AcceptTrack(StPicoTrack *trk, Float_t B, StThreeVectorF Vert);  // track accept cuts function
    Bool_t                  AcceptJet(StJet *jet);           // jets accept cuts function
    Bool_t                  DoComparison(int myarr[], int elems);
    void                    SetSumw2(); // set errors weights 
    //Double_t                EffCorrection(Double_t trkETA, Double_t trkPT, Int_t effswitch) const; // efficiency correction function

    // switches
    Bool_t                  doUsePrimTracks;         // primary track switch
    Int_t                   fDebugLevel;             // debug printout level
    Bool_t                  doPrintEventCounter;     // print event # switch
    Int_t                   fRunFlag;                // Run Flag numerator value
    Int_t                   fCentralityDef;          // Centrality Definition enumerator value
    Int_t                   fDoEffCorr;              // efficiency correction to tracks
    Bool_t                  fCorrJetPt;              // correct jet pt by rho

    // cuts
    Double_t                fMinPtJet;               // min jet pt to keep jet in output
    Double_t                fJetConstituentCut;      // min jet constituent
    Double_t                fTrackBias;              // high pt track in jet bias
    Double_t                fJetRad;                 // jet radius
    Double_t                fEventZVtxMinCut;        // min event z-vertex cut
    Double_t                fEventZVtxMaxCut;        // max event z-vertex cut
    Double_t                fTrackPtMinCut;          // min track pt cut
    Double_t                fTrackPtMaxCut;          // max track pt cut
    Double_t                fTrackPhiMinCut;         // min track phi cut
    Double_t                fTrackPhiMaxCut;         // max track phi cut
    Double_t                fTrackEtaMinCut;         // min track eta cut
    Double_t                fTrackEtaMaxCut;         // max track eta cut
    Double_t                fTrackDCAcut;            // max track dca cut
    Int_t                   fTracknHitsFit;          // requirement for track hits
    Double_t                fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax

    // centrality    
    Double_t                fCentralityScaled;           // scaled by 5% centrality 
    Int_t                   ref16;                       // multiplicity bin (16)
    Int_t                   ref9;                        // multiplicity bin (9)

    // event
    Double_t                Bfield;                      // event Bfield
    StThreeVectorF          mVertex;                     // event vertex 3-vector
    Double_t                zVtx;                        // z-vertex component

    // event selection types
    UInt_t                  fTriggerEventType;           // Physics selection of event used for signal
    Int_t                   fEmcTriggerArr[7];           // EMCal triggers array: used to select signal and do QA

    // used for event plane calculation and resolution
    StJet*                  fLeadingJet;//! leading jet
    StJet*                  fSubLeadingJet;//! subleading jet

    // clonesarray collections of tracks and jets
    TClonesArray           *fJets;//! jet collection

  private:
    Int_t                   fRunNumber;

    // PicoDstMaker and PicoDst object pointer
    StPicoDstMaker         *mPicoDstMaker;
    StPicoDst              *mPicoDst;
    StPicoEvent            *mPicoEvent;
    StJetMakerTask         *JetMaker;
    StRho                  *RhoMaker;

    // centrality objects
    StRefMultCorr          *grefmultCorr;
   
    // output file name string 
    TString                 mOutName;

    // Rho objects
    StRhoParameter         *GetRhoFromEvent(const char *name);
    StRhoParameter         *fRho;//!<!          // event rho
    Double_t                fRhoVal;//!<!       // event rho value, same for local rho
    TString                 fRhoName;///<       // rho name

    // maker names
    TString                 fAnalysisMakerName;
    TString                 fJetMakerName;
    TString                 fRhoMakerName;
    TString                 fEventMixerMakerName;

    ClassDef(StAnMaker, 1)
};
#endif
