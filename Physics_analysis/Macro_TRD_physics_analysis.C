

// Main macro to control analysis

// Environment variables
R__LOAD_LIBRARY(TRD_physics_analysis_cxx.so);

// Environment variables
//#define ENV_PI
#define ENV_ALEX
//#define ENV_PI_SVEN

void Macro_TRD_physics_analysis(TString input_list = "List_physics.txt")
{
// load the shared libraries here gSystem ->Load("...
    gSystem ->Load("TRD_physics_analysis_cxx.so");

	TString out_file_name = input_list;
    out_file_name += "_physics_out.root";

    TString input_dir  = "../ST_out/";
    TString output_dir = "./";
    TString inlists_dir = "./";
#if defined(ENV_PI)
    inlists_dir = "/home/ceres/schmah/ALICE/TRD_self_tracking/Physics_analysis/";
    input_dir   = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/ST_out/";
    output_dir  = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/";
#endif

#if defined(ENV_PI_SVEN)
    inlists_dir = "/home/ceres/schmah/ALICE/TRD_self_tracking/Lists_tracklets/"; // need to create a folder in HD if want to use 
    input_dir   = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/ST_out/";
    output_dir  = "/misc/alidata120/alice_u/hoppner/TRD_self_tracking/ST_out/";
#endif

#if defined(ENV_ALEX)
    inlists_dir = "./";
    input_dir  = "../ST_out/";
    output_dir = "./";
#endif

    Ali_TRD_physics_analysis*  TRD_physics_analysis = new Ali_TRD_physics_analysis(output_dir,out_file_name);
    //Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze(output_dir,"test.root",graphics);

    TRD_physics_analysis ->set_input_lists(inlists_dir);
    TRD_physics_analysis ->set_input_dir(input_dir);
    TRD_physics_analysis ->Init_tree(input_list.Data());

// for(Long64_t event = start_event; event < stop_event; event++)
//    {
//Ana->Init();
//Ana->Do_something(); etc.
}

