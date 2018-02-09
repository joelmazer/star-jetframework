
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
class StPicoBase;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();

// do more pt cuts of constituents
bool doAdditionalPtCuts = kFALSE;
bool usePrimaryTracks = kTRUE;

StChain *chain;
void readPicoDst(const Char_t *inputFile="test2.list", const Char_t *outputFile="test2.root", Int_t nEv = 10)
{
        Int_t nEvents = 100;
        if(nEv > 100) nEvents = 10000000;

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
        // 0.15 GeV + tracks
        StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", 1.0);
        jetTask->SetJetType(kChargedJet); //jetType
        jetTask->SetJetAlgo(antikt_algorithm); //jetAlgo
        jetTask->SetRecombScheme(BIpt2_scheme); //recomb
        jetTask->SetRadius(0.4);
        jetTask->SetJetsName("Jets");
        jetTask->SetMinJetPt(5.0); // 15.0 signal jets
        jetTask->SetMaxJetTrackPt(20.0);
        jetTask->SetGhostArea(0.005);
        jetTask->SetMinJetArea(0.0);
        jetTask->SetJetEtaRange(-0.6,0.6);
        jetTask->SetJetPhiRange(0,2*pi); 
        jetTask->SetUsePrimaryTracks(usePrimaryTracks);

        // create JetFinder for background now (JetMakerBG)
        StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2);
        jetTaskBG->SetJetType(kChargedJet); //jetType
        jetTaskBG->SetJetAlgo(kt_algorithm); //jetAlgo
        jetTaskBG->SetRecombScheme(BIpt2_scheme); //reco
        jetTaskBG->SetRadius(0.4); //0.4
        jetTaskBG->SetJetsName("JetsBG");
        jetTaskBG->SetMinJetPt(1.0);
        jetTaskBG->SetMaxJetTrackPt(20.0);
        jetTaskBG->SetGhostArea(0.005);
        jetTaskBG->SetMinJetArea(0.0);
        jetTaskBG->SetJetEtaRange(-0.6,0.6); //-0.5,0.5
        jetTaskBG->SetJetPhiRange(0,2*pi);  //0,pi
        jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);

        // this is the centrality dependent scaling for RHO as used in ALICE - don't use right now
        // s(Centrality) = 0.00015 ×Centrality^2 ? 0.016 ×Centrality + 1.91
        // scaled function for rho and setting parameters
        TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
        sfunc->SetParameter(2,1.80642);
        sfunc->SetParameter(1,-0.0112331);
        sfunc->SetParameter(0,0.0001139561);
  
        bool dohisto = kFALSE;
        // Rho task, and scale it up to include neutral constituents
        StRho *rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG");
        rhoTask->SetExcludeLeadJets(2);
        rhoTask->SetOutRhoName("OutRho");
        //rhoTask->SetScaleFunction(sfunc); // don't NEED

        // Rho Sparse
        //StRhoSparse *rhoTaskSparse = new StRhoSparse("StRhoSparse", kTRUE, outputFile);
        //rhoTaskSparse->SetExcludeLeadJets(2);

	// create the analysis maker!
        bool doComments = kFALSE;
	StMyAnalysisMaker *anaMaker = new StMyAnalysisMaker("AnalysisMakerTrackBias1GeV", picoMaker, outputFile, doComments, 5.0, 1.0, "JetMaker", "StRho_JetsBG");
        anaMaker->SetEventMixing(kTRUE);
        anaMaker->SetCentBinSize(5);
        anaMaker->SetMixingTracks(25000); // 50,000 in ALICE
        anaMaker->SetNMixedTr(2500);
        anaMaker->SetNMixedEvt(5);
        anaMaker->SetUsePrimaryTracks(usePrimaryTracks); // kFALSE
        anaMaker->SetCorrectJetPt(kFALSE); // kTRUE
        anaMaker->SetMinTrackPt(0.2);

        // initialize the chain
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
/*
  gSystem->Load("libTable");
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities");        // new addition 22jul99
  gSystem->Load("StTreeMaker");
  gSystem->Load("StIOMaker");
  gSystem->Load("StarClassLibrary");
  gSystem->Load("StTriggerDataMaker"); // new starting from April 2003
  gSystem->Load("StBichsel");
  gSystem->Load("StEvent");
  gSystem->Load("StEventUtilities");
  gSystem->Load("StDbLib");
  gSystem->Load("StEmcUtil");
  gSystem->Load("StTofUtil");
  gSystem->Load("StPmdUtil");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("libStStrangeMuDstMaker");
  gSystem->Load("libStMuDSTMaker");
  gSystem->Load("libStarAgmlUtil");

  gSystem->Load("StTpcDb");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("libgen_Tables");
  gSystem->Load("libsim_Tables");
  gSystem->Load("libglobal_Tables");
  gSystem->Load("StEmcTriggerMaker");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");  // XXX: loaded twice !
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcSimulatorMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StTriggerUtilities");

  gSystem->Load("StMagF");
  gSystem->Load("StMtdUtil");
  gSystem->Load("StMtdMatchMaker");
  gSystem->Load("StMtdCalibMaker");
*/

  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  //gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StMyAnalysisMaker");
  gSystem->Load("StPicoBase");

  gSystem->ListLibraries();

} void LoadMacros()
{

}
