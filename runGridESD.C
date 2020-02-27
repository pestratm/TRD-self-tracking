class  AliAnalysisManager;
class  AliAnalysisAlien;


// modes: "test" to run over a small set of files (requires alien connection but everything stored locally),
//        "full" to run over everything on the grid,
//        "terminate" to merge results after "full"
void runGridESD(TString mode="terminate",Int_t sub=702, TString fname="Ali_AS_analysis_TRD_digits", Int_t alien=1)
{
    // Use root 5

    // aliroot runGridESD.C\(\"test\",702,\"Ali_AS_analysis_TRD_digits\",0\)
    // aliroot runGridESD.C\(\"test\",702,\"Ali_AS_analysis_TRD_digits\",1\)
    // aliroot runGridESD.C\(\"full\",702,\"Ali_AS_analysis_TRD_digits\",1\)
    // aliroot runGridESD.C\(\"terminate\",702,\"Ali_AS_analysis_TRD_digits\",1\)

    cout << "Start macro runGridESD, alien: " << alien << endl;

    AliLog::SetGlobalDebugLevel(0); // 5

    //__________________________________________________________________________
    // Use AliRoot includes to compile our task
    //gROOT->ProcessLine(".include $ALICE_PHYSICS/include $ALICE_ROOT/include");



    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    // Create and configure the alien handler plugin
    if(alien==1)
    {
	AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,sub,fname);
	if (!alienHandler) return;
    }

    // Create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
    if(alien==1) mgr->SetGridHandler(alienHandler);
    // Connect plug-in to the analysis manager

    //AliESDInputHandler* esdH = new AliESDInputHandler();
    //mgr->SetInputEventHandler(esdH);


#if 0
    AliVEventHandler* handler;
    // gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C"));
    // handler = AddMCHandler(kFALSE);
    // ((AliMCEventHandler*)handler)->SetReadTR(kFALSE);
    // mgr->SetMCtruthEventHandler(handler);
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));
    handler = AddESDHandler();
    ((AliESDInputHandler*)handler)->SetFriendFileName("AliESDfriends.root");
    ((AliESDInputHandler*)handler)->SetReadFriends(kTRUE);
    ((AliESDInputHandler*)handler)->SetNeedField();
    mgr->SetInputEventHandler(handler);
#endif

    AliESDInputHandler *esdH = new AliESDInputHandler();
    esdH->SetFriendFileName("AliESDfriends.root");
    esdH->SetReadFriends(kTRUE);
    mgr->SetInputEventHandler(esdH);


    //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); // No idea <-
    //  AddTaskPhysicsSelection(kTRUE);
    //AddTaskPhysicsSelection(kFALSE,kTRUE); // <-
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); // No idea <-
    AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kTRUE,"1"); // <-
    //AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse();


    //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    //AliCentralitySelectionTask *taskCentrality=AddTaskCentrality();
    //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    //AddTaskMultSelection(); // files need to have full correct path locally

    // to run a local task not in AliPhysics (there will be conflicts if a task in AliPhysics has the same name)

    cout << "" << endl;
    cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "----------------------- Load analysis macros -----------------------" << endl;
    cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    gROOT->LoadMacro("Ali_AS_analysis_TRD_digits.cxx+g");
    cout << "Loaded macro Ali_AS_analysis_TRD_digits.cxx" << endl;
    gROOT->LoadMacro("AddTask_aschmah.C");
    cout << "Loaded macro AddTask_aschmah.C" << endl;
    AddTask_aschmah();
    cout << "" << endl;
    cout << "----------------------- Added task -----------------------" << endl;
    cout << "" << endl;

    if( !mgr->InitAnalysis() ) return;

    mgr->PrintStatus();

    // Start analysis in grid.
    //  mgr->StartAnalysis("grid",5); // just some chunks for testing purposes

    cout << "" << endl;
    cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Start analysis on grid" << endl;
    cout << "-------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "" << endl;
    //mgr->StartAnalysis("grid",300); // all chunks, N events
    mgr->StartAnalysis("grid"); // all chunks, N events

    cout << "" << endl;
    cout << "Analysis done" << endl;

    // to run over files stored locally, uncomment this section,
    // and comment out the above lines related to alienHandler and StartAnalysis("grid")
    if(alien==0)
    {
        cout << "----------------------- Run local job -----------------------" << endl;
	TChain *chain = new TChain("esdTree");
	//chain->AddFile("/home/ceres/schmah/ALICE/TRD/Data/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1001/AliESDs.root");
        //chain->AddFile("/misc/alidata120/alice_u/schmah/alice/data/2015/LHC15o/000245450/pass1/15000245450036.1001/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1012/AliESDs.root");

	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1008/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1012/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1004/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1011/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000245450036.1008/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000246980039.9708/AliESDs.root");
	//chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000246980039.9709/AliESDs.root");
        //chain->AddFile("/misc/alidata120/alice_u/schmah/PbPb_2015/raw/15000246980039.9711/AliESDs.root");
        chain->AddFile("/misc/alidata120/alice_u/schmah/pPb_2016/16000265338039.6007/AliESDs.root");
	//chain->Print();
	mgr->StartAnalysis("local",chain);
    }


}



AliAnalysisGrid* CreateAlienHandler(TString mode="test",Int_t sub=0,TString fname="testName")
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    plugin->SetOverwriteMode();
    plugin->SetExecutableCommand("aliroot -q -b");
    plugin->SetRunMode(mode.Data());
    plugin->SetNtestFiles(1); // Number of files used in the testing case
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    //plugin->SetROOTVersion("v5-34-26");
    //plugin->SetAliROOTVersion("v5-08-18c-1 ");
    //plugin->SetAliPhysicsVersion("vAN-20170831-1"); // change to something up-to-date
    //plugin->SetAliPhysicsVersion("vAN-20171012-1"); // change to something up-to-date
    plugin->SetAliPhysicsVersion("vAN-20181112-1");
    // Declare input data to be processed.
    plugin->SetNrunsPerMaster(1);   // 1 folder per run; then 1 masterjob per run ; 1 merged file per run ; if 10 then 10 runs in 1 masterjob, merged all in one file; merging in stages if large files, might fail;
    plugin->SetMaxMergeStages(1);  // if >1 iterative merging; all subjobs of single masterjob in one file
    plugin->SetSplitMaxInputFileNumber(50); // 3 in the LEGO trains, 100, 1000 for more files combined together



    if(sub==702)
    {   // LHC16q p-Pb 5 TeV
	plugin->SetGridDataDir("/alice/data/2016/LHC16q/"); // OK
	plugin->SetGridWorkingDir(Form("%s/sub%d/",fname.Data(),sub)); // No idea
	plugin->SetDataPattern("*/pass1_CENT_wSDD/*/AliESDs.root"); // OK
	plugin->SetAnalysisMacro(Form("TaskTrackAna%d.C",sub));
	plugin->SetExecutable(Form("TaskTrackAna%d.sh",sub));
	plugin->SetJDLName(Form("TaskTrackAna%d.jdl",sub));
	plugin->SetRunPrefix("000");
	//Int_t runnumbers[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, // OK
	//265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377,
        //265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309}; // 32 runs in total

        Int_t runnumbers[] = {265338, 265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, // OK
	265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377,
	265344, 265343, 265342, 265339, 265336, 265334, 265332, 265309}; // 32 runs in total
	for(Int_t irun = 0; irun < 1; irun++)
	{
	    Printf("%d %d",irun,runnumbers[irun]);
	    plugin->AddRunNumber(runnumbers[irun]);
	}
    }

    if(sub==703)
    {   // LHC15o Pb-Pb 5 TeV
	plugin->SetGridDataDir("/alice/data/2015/LHC15o/"); // OK
	plugin->SetGridWorkingDir(Form("%s/sub%d/",fname.Data(),sub)); // No idea
	//plugin->SetDataPattern("/pass2_UD_CCUP/*/AliESDs.root"); // OK
        plugin->SetDataPattern("/pass1/*/AliESDs.root"); // OK
	plugin->SetAnalysisMacro(Form("TaskTrackAna%d.C",sub));
	plugin->SetExecutable(Form("TaskTrackAna%d.sh",sub));
	plugin->SetJDLName(Form("TaskTrackAna%d.jdl",sub));
	plugin->SetRunPrefix("000");
	Int_t runnumbers[] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};
	for(Int_t irun = 0; irun < 1; irun++)
	{
	    Printf("%d %d",irun,runnumbers[irun]);
	    plugin->AddRunNumber(runnumbers[irun]);
	}
    }



    TString extraLibs;
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

    //plugin->SetAdditionalLibs("Ali_AS_Event.cxx Ali_AS_Event.h Ali_AS_analysis_TRD_digits.cxx Ali_AS_analysis_TRD_digits.h");
    //plugin->SetAdditionalLibs("Ali_AS_Event.h Ali_AS_Event.cxx Ali_AS_analysis_TRD_digits.h Ali_AS_analysis_TRD_digits.cxx");
    plugin->SetAdditionalLibs("Ali_AS_Event.h Ali_AS_EventLinkDef.h Ali_AS_analysis_TRD_digits.h Ali_AS_analysis_TRD_digits.cxx");
    //plugin->SetAdditionalLibs("Ali_AS_Event_cxx.so Ali_AS_analysis_TRD_digits_cxx.so"); // doesn't work at all
    // needed if running over a local task, comment out if the task is in AliPhysics
    plugin->SetAnalysisSource("Ali_AS_analysis_TRD_digits.cxx");

    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    //plugin->SetDefaultOutputs(kFALSE); // bool, either 0 or 1 (<-)
    //plugin->SetDefaultOutputs(kTRUE); // bool, either 0 or 1 (<-)
    plugin->SetOutputToRunNo(1);
    //plugin->SetOutputFiles("AnalysisResults_blubb.root");
    //plugin->SetOutputFiles("TRDPIDTree_tree_T1.root TRDPIDTree_hists_T1.root");
    //plugin->SetOutputFiles("TRDPIDTree_hists_T0.root TRDPIDTree_hists_T1.root"); // Root output file
    //plugin->SetOutputFiles("TRDPIDTree_hists_T0.root"); // Root output file
    plugin->SetMergeViaJDL(kTRUE);
    plugin->SetOneStageMerging(kFALSE);
    //plugin->SetMaxMergeFiles(40);
    //plugin->SetMaxMergeStages(4);

    //plugin->SetTTL(86399); // Time to live
    plugin->SetTTL(39999); // Time to live in seconds
    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");
    plugin->SetKeepLogs(kTRUE); // keeps the logs
    // Optionally modify job price (default 1)
    plugin->SetPrice(1); // Grid price for the job;
    // Optionally modify split mode (default 'se')
    //plugin->SetSplitMaxInputFileNumber();
    plugin->SetSplitMode("se");
    return plugin;
}
