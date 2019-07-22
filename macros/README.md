# Descriptions for each macro:
*this will be updated as new macros are added...*

readPicoDstDummyMaker.C
* MAIN macro for NEW users - has the basic setup for jets, rho and runs a dummy class:  StDummyMaker, use this to start an analysis, should work right out of the box.  This is a NEW version of readPicoDstTest.C set up around the StDummyMaker class.

readPicoDstTest.C
* macro for NEW users - has the basic setup for jets, rho and runs a dummy class:  StAnMaker, use this to start an analysis

readPicoDstMakeTrees.C
* macro used to create Trees - may need updating 

readPicoDstMultPtBins.C
* used for the correlation and jet shape analysis - my main macro

readPicoDstQA.C
* used to run StPicoTrackClusterQA class which runs various QA, generates histograms used for bad/hot tower QA filtering

readPicoDst.C
* OLD version that needs to be updated

## Author
**Joel Mazer**

## Notes
May 6th, 2019: added switch for loading and omitting bad runs

June 10th, 2019: removed readPicoDstJoel.C, updated for new centrality definitions

July 22nd, 2019: added readPicoDstDummyMaker.C over past few days
