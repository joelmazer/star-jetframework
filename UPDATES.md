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
1) see trigger flags section in readPicoDstDummyMaker.C, as there are now 3 bad tower lists, for different constituent cuts. And it is set from the 
```
  double fJetConstituentCut = 2.0;
```

2) remove instances of the SetCentralityDef() setter from your classes if they exist and use the cent maker as used in StDummyMaker.cxx and readPicoDstDummyMaker.C
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
