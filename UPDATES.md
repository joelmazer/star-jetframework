# star-jetframework
changes to framework which may require a change to other code will be listed below


* May 22, 2019: moving calls of CheckForMB() and CheckForHT() from StJetMakerTask_ to frameworks BaseClass (StJetFrameworkPicoBase)
Add the following right after the class initialiation of StPicoDstMaker in your readPicoDst.C macro
```  
  // (updated August 20th, 2019)
  // create base class maker pointer
  StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
  baseMaker->SetRunFlag(RunFlag);                         // run flag (year)
  baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
  baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
  baseMaker->SetBadTowerListVers(TowerListToUse); 
```

This creates an instance of the base class to allow classes which don`t inherit from the base class to be able to access its member functions directly from the class object. This will also require adding the following at the beginning of a Make() function which does not inherit from the base class:
```
  // get base class pointer
  // this class does not inherit from base class: StJetFrameworkPicoBase, but we want to reduce redundancy
  StJetFrameworkPicoBase *baseMaker = static_cast<StJetFrameworkPicoBase*>(GetMaker("baseClassMaker"));
  if(!baseMaker) { cout<<"no baseMaker.. returning"<<endl;  return kStOK; }
```

# Other update notes...
* August 20, 2019
Changes needed from past update to be added to your readPicoDst.C macro after instance of base class:
```
        // create centrality class maker pointer
        StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, outputFile, doComments);
        CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
        CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
        CentMaker->SetRunFlag(RunFlag);                         // Run Flag
        CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
        CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
        CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used?
        CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        cout<<CentMaker->GetName()<<endl;  // print name of class instance
```
new changes or notes:
1) see tower flags section in readPicoDstDummyMaker.C, as there are now 3 bad tower lists, for different constituent cuts. And it is set from the 
```
  double fJetConstituentCut = 2.0;
```
The bad tower list used is decided based on your constituent cut.

2) remove instances of the SetCentralityDef() setter from your classes (except StCentMaker) if they exist and use the cent maker as used in StDummyMaker.cxx and readPicoDstDummyMaker.C
```
  jetTask->SetCentralityDef(CentralityDefinition); // run based centrality definition
```

3) new options have been added for how to perform that hadronic correction when using the neutral component for jets
  ..* using the last matched indexed track of the BEmcPidTraits cluster containing tower: StJetFrameworkPicoBase::kLastMatchedTrack
  ..* using the highest E matched indexed track of the BEmcPidTraits cluster containing tower:       StJetFrameworkPicoBase::kHighestEMatchedTrack 
  ..* using the sum of all matched tracks of the BEmcPidTraits cluster containing tower:
StJetFrameworkPicoBase::kAllMatchedTrack

The default is: StJetFrameworkPicoBase::kAllMatchedTrack
To change, add the corresponding setting in your readPicoDst.C macro:
```
jetTask->SetJetHadCorrType( type );
```
4) If one has in their readPicoDst.C macro or were using a class as their own from the framework, please update the Setters below.  Changes were made to the name from E to Et.
```
anaMaker[i]->SetJetMaxTowerE( val ); 
virtual void            SetMaxEventTowerE(Double_t mxE)  { fMaxEventTowerE = mxE; }

```
to
```
anaMaker[i]->SetJetMaxTowerEt( val ); 
virtual void            SetMaxEventTowerEt(Double_t mxEt)  { fMaxEventTowerEt = mxEt; }
```

5) There were changes to correct a problem with the GetMomentum() function in StJetMakerTask. It priorly used the tower energy instead of hadronically corrected tower energy when input to fastjet.  Nothing needs to be done on USERS end. 

6) many functions were ported to the base class to reduce redundancy, please see StJetMakerTask as an example on how to call things in classes that don't inherit from the base class. Or use/see StDummyMaker for a class that does inherit from the base class.  From StJetMakerTask:
```
  bool fHaveMB30event  = mBaseMaker->CheckForMB(fRunFlag, StJetFrameworkPicoBase::kVPDMB30);
```
7) please also update the centrality package to have updated defintions and files in your StRoot/StRefMultCorr directory.


*August 22, 2019
1) added new centrality definition to readMacros.C - can be used for all luminosities, and is set as default
2) added SetJetHadCorrType() to the readMacros.C so USERS can see it - set to default: kAllMatchedTracks



*February 12, 2020 - Part 1:
1) old(er) event plane correction (recentering and shifting) values were removed and NEW(er) values pushed to the repo
2) StAnMaker class removed along with its corresponding macro: readPicoDstTest.C - it has a replacement:  StDummyMaker and readPicoDstDummyMaker.C
3) trackingEfficiency_Run14.h was added, but not used - .root files are instead loaded - they may be added to the repo in upcoming updates, but are located on RCF


*February 12, 2020 - Part 2:
:::: Notable changes ::::
1) StCentMaker: method added for all luminosities for Run 14 
```
getgRefMultCorr_P18ih_VpdMB30_AllLumi()
```
2) StEventPlaneMaker: cleaned up, consolidation, added switches for additional jet resolutions (jet R),
3) StJetFrameworkPicoBase:
  - added another centrality bin function:  
  ```
  Int_t StJetFrameworkPicoBase::GetCentBin10(Int_t cbin) const
  ```
  - added luminosity bin function (used for Run 14 tracking efficiency)
  ```
  Int_t StJetFrameworkPicoBase::GetLuminosityBin(Double_t lumi) const  
  ```
  - updated and cleaned the   ApplyTrackingEff() function.  Tracking efficiency can be called in a class which inherits from the base class via:
 
4) StJetMakerTask:
  - added switches for tracking efficiency
  - loads file for tracking efficiency
  ```
    // input file - for tracking efficiency: Run14 AuAu and Run12 pp
  const char *input = "";
  if(fRunFlag == StJetFrameworkPicoBase::Run14_AuAu200) input=Form("./StRoot/StMyAnalysisMaker/Run14_efficiencySmaller2D.root");
  if(fRunFlag == StJetFrameworkPicoBase::Run12_pp200)   input=Form("./StRoot/StMyAnalysisMaker/Run12_efficiency_New.root"); // Oct 17, 2019 added
  if(fDoEffCorr) {
    fEfficiencyInputFile = new TFile(input, "READ");
    if(!fEfficiencyInputFile) cout<<Form("do not have input file: %s", input);
  }
  ```
  - switch added for Thomas to correct tracks for efficiency before passing to FastJet (done for charged jets currently) and thus giving back jets with efficiency corrected pt.  Need to set in your readPicoMacro.C
 ``` 
    virtual void         SetDoCorrectTracksforEffBeforeJetReco(Bool_t ec) {   doCorrectTracksforEffBeforeJetReco = ec;
 ```
where the following is done
  ```
      // - NEW TODO - for charged tracks only - CHARGED JETS only (January 22, 2020)
      // calcualte tracking efficency - calculate single particle track efficiency
      int effCent   = mCentMaker->GetRef16();
      double fZDCx  = mPicoEvent->ZDCx();
      double trkEff = mBaseMaker->ApplyTrackingEff(fDoEffCorr, pt, eta, effCent, fZDCx, fTrackEfficiencyType, fEfficiencyInputFile); // may not need to correct here, just analysis stage
      double pxCorr = px / trkEff;
      double pyCorr = py / trkEff;
      double pzCorr = pz / trkEff;
      double energyCorr = 1.0*TMath::Sqrt(pxCorr*pxCorr + pyCorr*pyCorr + pzCorr*pzCorr + pi0mass*pi0mass);

      // add track input vector to FastJet
      if(doCorrectTracksforEffBeforeJetReco) {
        fjw.AddInputVector(pxCorr, pyCorr, pzCorr, energyCorr, iTracks); // efficiency corrected components
      } else {
        // THIS IS DEFAULT !
        // send track info to FJ wrapper
        fjw.AddInputVector(px, py, pz, energy, iTracks); // includes E
      }
  ```
    
  and similarly within the constitutent subtractor switch
  - consolidation and clean-up
  - now call a couple methods from the base class
  - new histograms added for QA with track-tower matching and hadronic corrections
  - if YOU want to look into jet mass, you *WILL* need to set the recombination scheme in your readPicoMacro.C
  ```
  jetTask->SetRecombScheme(E_scheme); 
  ```
5) added capabilities to classes which inherit from the base class, to access and apply tracking efficiency for Run 12 pp and Run 14 AuAu.

6) StMyAnalysisMaker3: much has been added/tweaked in my combined class for jet-h and jet shape
  - added more histograms
  - updated added tracking efficiency switches / usage, loads efficiency file(s), 
  - added switches / members for:
    a) biasing jets
    b) filtering single particle jets
    c) removing mixed events with which have cone with specific fractions of energy of signal jet
    d) running with some systematic uncertainty variations
    e) jet towers firing trigger requirement function
  - cleaned
  - updated main class constructor
  - added more exclusion checks to hStats histogram
  - updated event plane option to mixed events
  - updated event plane type usage to follow scheme of enumerator
  

IF THERE IS ANYTHING ELSE - please me know or update this file yourself and push change.


# Task lists - *to be updated*
- [] port more functions from classes not inherting from base class to base class calls 
- [] find more straightforward way for other users to incorporate git framework into cons / StRoot setup instead of symbolic linking
- []
- []

## Author
**Joel Mazer**

## Notes
See README.md and class documentation for additional details
