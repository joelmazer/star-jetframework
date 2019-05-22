# star-jetframework
changes to framework which may require a change to other code will be listed below


* May 22, 2019: moving calls of CheckForMB() and CheckForHT() from StJetMakerTask_ to frameworks BaseClass (StJetFrameworkPicoBase)
Add the following right after the class initialiation of StPicoDstMaker in your readPicoDst.C macro
```
  // create base class maker pointer
  StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
  baseMaker->SetRunFlag(RunFlag);           // run flag (year)
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

# Task lists - *to be updated*
- [] port more functions from classes not inherting from base class to base class calls 
- [] find more straightforward way for other users to incorporate git framework into cons / StRoot setup instead of symbolic linking
- []
- []

## Author
**Joel Mazer**

## Notes
See README.md and class documentation for additional details
