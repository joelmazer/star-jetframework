// ************************************** //
// Sept28, 2017
// some current notes:
// Run16 AuAu 200 GeV:  P16ij on RCAS has no EmcalTriggers and no clusters
//     ::P16ij has triggers that need to be accessed a different way
//
//

#include <TSystem>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
// my added STAR classes
class StJetMakerTask;
class StRho;
class StRhoBase;
class StRhoSparse;
class StMyAnalysisMaker;
class StJetFrameworkPicoBase;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();

// do more pt cuts of constituents
bool doAdditionalPtCuts = kFALSE;
bool usePrimaryTracks = kFALSE;

// additional switches
bool doEventPlaneCorrections = kTRUE;   // needs to be set to kTRUE to run corrections
bool phi_shift_switch = kFALSE;         // not used! (TPC)
bool tpc_recenter_read_switch = kFALSE; // STEP1: kTRUE = write output for recentering (shift)
bool tpc_shift_read_switch = kFALSE;    // STEP2: kTRUE = write output for shift (flattening)
bool tpc_apply_corr_switch = kFALSE;    // STEP3: kTRUE = read shift info for final TPC calculation
bool zdc_recenter_read_switch = kFALSE; // STEP1: kTRUE = write output for recentering of zdc (shift)
bool zdc_shift_read_switch = kFALSE;    // STEP2: kTRUE = write output for zdc shift (flattening)
bool zdc_apply_corr_switch = kFALSE;    // STEP3: kTRUE = read shift info for final ZDC calculation
bool bbc_recenter_read_switch = kFALSE; // STEP1: kTRUE = write output for recentering of bbc (shift)
bool bbc_shift_read_switch = kFALSE;    // STEP2: kTRUE = write output for bbc shift (flattening)
bool bbc_apply_corr_switch = kFALSE;    // STEP3: kTRUE = read shift info for final BBC calculation
bool doSTEP1 = kFALSE;
bool doSTEP2 = kFALSE;
bool doSTEP3 = kTRUE;

Double_t ZVtxMin = -20.0;
Double_t ZVtxMax = 20.0;

StChain *chain;

//void readPicoDst(const Char_t *inputFile="test2.list", const Char_t *outputFile="test2.root", Int_t nEv = 10)
//void readPicoDst(const Char_t *inputFile="newPicoDsts.list", const Char_t *outputFile="newtest.root", Int_t nEv = 10)
void readPicoDstJoel(const Char_t *inputFile="test_run17061011.list", const Char_t *outputFile="test.root", Int_t nEv = 10)
{
//        Int_t nEvents = 10000000;
//	Int_t nEvents = 25000;
//        Int_t nEvents = 10;
        Int_t nEvents = 10000;
//        Int_t nEvents = 1000;
        if(nEv > 100) nEvents = 10000000;

        if(doSTEP1) { bbc_recenter_read_switch = kTRUE; }
        if(doSTEP2) { bbc_shift_read_switch = kTRUE; }
        if(doSTEP3) { bbc_apply_corr_switch = kTRUE; bbc_shift_read_switch = kTRUE; bbc_recenter_read_switch = kTRUE; } 

        if(doSTEP1) { tpc_recenter_read_switch = kTRUE; }
        if(doSTEP2) { tpc_shift_read_switch = kTRUE; }
        if(doSTEP3) { tpc_apply_corr_switch = kTRUE; tpc_shift_read_switch = kTRUE; tpc_recenter_read_switch = kTRUE; }

        if(doSTEP1) { zdc_recenter_read_switch = kTRUE; }
        if(doSTEP2) { zdc_shift_read_switch = kTRUE; }
        if(doSTEP3) { zdc_apply_corr_switch = kTRUE; zdc_shift_read_switch = kTRUE; zdc_recenter_read_switch = kTRUE; }


        // set up Jet enumerators:
        // jet type
        enum EJetType_t {
          kFullJet,  // tracks + clusters
          kChargedJet,
          kNeutralJet
        };

        // jet algorithm
        enum EJetAlgo_t {
          kt_algorithm                    = 0, // background jets
          antikt_algorithm                = 1, // signal jets
          cambridge_algorithm             = 2,
          genkt_algorithm                 = 3,
          cambridge_for_passive_algorithm = 11,
          genkt_for_passive_algorithm     = 13,
          plugin_algorithm                = 99,
          undefined_jet_algorithm         = 999
        };

        // jet recombination scheme
        enum ERecoScheme_t {
          E_scheme        = 0,
          pt_scheme       = 1,
          pt2_scheme      = 2,
          Et_scheme       = 3,
          Et2_scheme      = 4,
          BIpt_scheme     = 5,
          BIpt2_scheme    = 6,
          WTA_pt_scheme   = 7,
          WTA_modp_scheme = 8,
          external_scheme = 99
        };

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        Bool_t doCentSelection = kTRUE; //kTRUE;
        Int_t centralitySelection = StJetFrameworkPicoBase::kCent2050;
	
        // open and close output .root file (so it exist and can be updated by Analysis Tasks)
        TFile *fout = new TFile(outputFile, "RECREATE");
        //fout->cd();
        fout->Close();

        // create chain
        StChain* chain = new StChain();

	// create the picoMaker maker
	//StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
        StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst"); // updated Aug6th
        picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

        // if(bFillGhost) jetTask->SetFillGhost();
        // create JetFinder first (JetMaker)
        // 0.2 GeV + tracks
        StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", 1.0, kTRUE, outputFile);
        jetTask->SetJetType(kFullJet);           //kChargedJet); //jetType
        jetTask->SetJetAlgo(antikt_algorithm); //jetAlgo
        jetTask->SetRecombScheme(BIpt2_scheme); //recomb
        jetTask->SetRadius(0.4);
        jetTask->SetJetsName("Jets");
        jetTask->SetMinJetPt(10.0); // 15.0 signal jets
        jetTask->SetMaxJetTrackPt(20.0);
        jetTask->SetGhostArea(0.005);
        jetTask->SetMinJetArea(0.0);
        jetTask->SetJetEtaRange(-0.6,0.6);
        jetTask->SetJetPhiRange(0,2*pi); 
        jetTask->SetUsePrimaryTracks(usePrimaryTracks);
        //jetTask->SetDebugLevel(2);
        jetTask->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        jetTask->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        jetTask->SetTurnOnCentSelection(doCentSelection);
        jetTask->SetCentralityBinCut(centralitySelection);

/*
        StJetMakerTask* jetTask = new StJetMakerTask(name);
        if (bFillGhosts) jetTask->SetFillGhost();
        if (lockTask) jetTask->SetLocked();
*/

/*
        // create JetFinder for background now (JetMakerBG)
//        StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
        StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 2.0, kTRUE, outputFile); // TEST
        jetTaskBG->SetJetType(kChargedJet); //jetType
        jetTaskBG->SetJetAlgo(kt_algorithm); //jetAlgo
        jetTaskBG->SetRecombScheme(BIpt2_scheme); //reco
        jetTaskBG->SetRadius(0.4); //0.4
        jetTaskBG->SetJetsName("JetsBG");
        jetTaskBG->SetMinJetPt(0.0);
        jetTaskBG->SetMaxJetTrackPt(20.0);
        jetTaskBG->SetGhostArea(0.005);
        jetTaskBG->SetMinJetArea(0.0);
        jetTaskBG->SetJetEtaRange(-0.6,0.6); //-0.5,0.5
        jetTaskBG->SetJetPhiRange(0,2*pi);  //0,pi
        jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
        jetTaskBG->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        jetTaskBG->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        jetTaskBG->SetTurnOnCentSelection(doCentSelection);
        jetTaskBG->SetCentralityBinCut(centralitySelection);
*/

        // this is the centrality dependent scaling for RHO as used in ALICE - don't use right now
        // s(Centrality) = 0.00015 ×Centrality^2 ? 0.016 ×Centrality + 1.91
        // scaled function for rho and setting parameters
        TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
        sfunc->SetParameter(2,1.80642);
        sfunc->SetParameter(1,-0.0112331);
        sfunc->SetParameter(0,0.0001139561);
  
        // name string for Rho
        //TString name(Form("StRho_%s_%s", nJets,cutType));
        bool dohisto = kFALSE;

        // Rho task, and scale it up to include neutral constituents
//        StRho *rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG");
        StRho *rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMaker");
        rhoTask->SetExcludeLeadJets(2);
        rhoTask->SetOutRhoName("OutRho");
        rhoTask->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        rhoTask->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTask->SetTurnOnCentSelection(doCentSelection);
        rhoTask->SetCentralityBinCut(centralitySelection);
        //rhoTask->SetScaleFunction(sfunc); // don't NEED
        //AddTaskRho(sBkgdJetsName, sTracksName, sClusName,"Rho", kJetRadius,"TPC",0.01,0,sfunc,2,kTRUE,"Rho");

        // Rho Base
        //StRhoBase *rhoTaskBase = dynamic_cast<StRhoBase*>rhoTask;

/*
        // Rho Sparse
        //StRhoSparse *rhoTaskSparse = dynamic_cast<StRhoSparse*>rhoTask;
        StRhoSparse *rhoTaskSparse = new StRhoSparse("StRhoSparse_JetsBG", kTRUE, outputFile, "JetMakerBG");
        rhoTaskSparse->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        rhoTaskSparse->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        rhoTaskSparse->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTaskSparse->SetExcludeLeadJets(2);
        rhoTaskSparse->SetJetMakerName("JetMaker");
        rhoTaskSparse->SetJetBGMakerName("JetMakerBG");
        rhoTaskSparse->SetDebugLevel(0); // 0 does nothing
        rhoTaskSparse->SetTurnOnCentSelection(doCentSelection);
        rhoTaskSparse->SetCentralityBinCut(centralitySelection);
        //rhoTaskSparse->SetScaleFunction(sfunc); // don't NEED
*/

	// create the analysis maker!
        bool doComments = kFALSE;
        //StMyAnalysisMaker *anaMaker = new StMyAnalysisMaker("AnalysisMaker", picoMaker, outputFile, doComments, 5.0, 0.15, "JetMaker", "StRho_JetsBG");
	StMyAnalysisMaker *anaMaker = new StMyAnalysisMaker("AnalysisMakerTrackBias1GeV", picoMaker, outputFile, doComments, 10.0, 4.0, "JetMaker", "StRho_JetsBG");
        anaMaker->SetEventMixing(kTRUE);
        anaMaker->SetCentBinSize(5);
        anaMaker->SetMixingTracks(25000); // 50,000 in ALICE
        anaMaker->SetNMixedTr(1500);      // was 2500
        anaMaker->SetNMixedEvt(2);        // was 5
        anaMaker->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
        anaMaker->SetCorrectJetPt(kFALSE); // subtract Rho BG
        anaMaker->SetMinTrackPt(0.2);
        anaMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        anaMaker->SetTrackPhiRange(0.0, 2*TMath::Pi());
        anaMaker->SetTrackEtaRange(-1.0, 1.0);
        //if(nEv<100) anaMaker->SetDebugLevel(7); // 0 turned off
        anaMaker->SetEventPlaneMaxTrackPtCut(5.0); // default
        anaMaker->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        anaMaker->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        anaMaker->SetTriggerEventType(StJetFrameworkPicoBase::kIsHT1);  // kIsHT1 or kIsHT2
        anaMaker->SetPrintEventCounter(kFALSE); 
        anaMaker->SetTurnOnCentSelection(doCentSelection);         // run analysis for specific centrality
        anaMaker->SetCentralityBinCut(centralitySelection);        // specific centrality range to run 
        //if(doSTEP1) anaMaker->SetEPcalibFileName("");                                                  // STEP1 - nothing
        //if(doSTEP2) anaMaker->SetEPcalibFileName("StRoot/StMyAnalysisMaker/recenter_calib_file.root"); // STEP2 read
        //if(doSTEP3) anaMaker->SetEPcalibFileName("StRoot/StMyAnalysisMaker/shift_calib_file.root");    // STEP3 read
        if(doSTEP1) anaMaker->SetOutFileNameEP("recenter_calib_file.root");  // set output file for BBC calibration of event plane
        if(doSTEP2) anaMaker->SetOutFileNameEP("shift_calib_file.root"); // set output file for BBC calibration of event plane
        if(doSTEP3) anaMaker->SetOutFileNameEP("final_eventplane.root"); // set output file for final BBC calculation

        // switches for event plane correction
        if(doEventPlaneCorrections) {
          //anaMaker->SetPhiShift(phi_shift_switch); // keep off, not used
          anaMaker->SetTPCRecenterRead(tpc_recenter_read_switch);   // TPC - STEP1
          anaMaker->SetTPCShiftRead(tpc_shift_read_switch);         // TPC - STEP2
          anaMaker->SetTPCApplyCorrections(tpc_apply_corr_switch);  // TPC - STEP3
          anaMaker->SetZDCrecentering(zdc_recenter_read_switch);    // ZDC - STEP1
          anaMaker->SetZDCShiftRead(zdc_shift_read_switch);         // ZDC - STEP2 
          anaMaker->SetZDCApplyCorrections(zdc_apply_corr_switch);  // ZDC - STEP3
          anaMaker->SetBBCrecentering(bbc_recenter_read_switch);    // BBC - STEP1
          anaMaker->SetBBCShiftRead(bbc_shift_read_switch);         // BBC - STEP2
          anaMaker->SetBBCApplyCorrections(bbc_apply_corr_switch);  // BBC - STEP3
        }


        // initialize chain
        chain->Init();
        cout<<"chain->Init();"<<endl;
        int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
        if(nEvents>total) nEvents = total;
  
        for (Int_t i=0; i<nEvents; i++){
          if(i%100==0) cout << "Working on eventNumber " << i << endl;

          chain->Clear();
          int iret = chain->Make(i);	
          if (iret) { cout << "Bad return code!" << iret << endl; break;}

          total++;		
	}
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;	

        // close output file if open
        if(fout->IsOpen()) fout->Close();

        //StMemStat::PrintMem("load StChain");
}

void LoadLibs()
{
  // load fastjet libraries 3.x
  //gSystem->Load("libCGAL"); - not installed 
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjettools");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");

  // add include path to use its functionality
  gSystem->AddIncludePath("-I/global/homes/j/jmazer/STAR/Y2017/mytestinstalldir/FastJet/fastjet-install/include");

  // load the system libraries - these were defaults
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // these are needed for new / additional classes
//  gSystem->Load("StPicoDstMaker");
  //gSystem->Load("StJetMakerTask");
  //gSystem->Load("StRho");
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  //gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StMyAnalysisMaker");

  gSystem->ListLibraries();
} 

void LoadMacros()
{
    //gROOT->LoadMacro("StRoot/StMyAnalysisMaker/ConfigJetFinders.C");
}
