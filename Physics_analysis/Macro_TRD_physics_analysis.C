

// Main macro to control analysis

// first need to compile physics class: .L TRD_physics_analysis.cxx++
// then run Macro:                      .x Macro_TRD_physics_analysis.C 

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
    
    Long64_t N_Events = TRD_physics_analysis ->get_N_Events();

    Int_t start_event = 0;
    Int_t stop_event  = (Int_t) N_Events;
    //Int_t stop_event  = 15;


    for(Long64_t event = start_event; event < stop_event; event++)
    {
        //-------> Tracklets = new Ali_TRD_ST_Tracklets*[NumTracklets]; what is this??
        //class _copy
        //all Floates?

        if (event != 0  &&  event % 50 == 0)
        cout << "." << flush;
        if (event != 0  &&  event % 500 == 0)
        {
            printf("event: %lld out of %lld, %4.2f%% total done \n",event,N_Events,((Double_t)event/(Double_t)N_Events)*100.0);
        }
        
        TRD_physics_analysis ->Loop_event(event);

        //printf("\n Event: %lld, Photons: %d, Nuclear interactions: %d \n",event,(Int_t)vec_PhotonVertex.size(),(Int_t)vec_NIVertex.size());

    }
}

