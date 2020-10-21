
#include "TRD_TFTracker.h"
#include "TRD_TFTracker.cxx"

R__LOAD_LIBRARY(TRD_Kalman_Tracking_cxx.so);
R__LOAD_LIBRARY(TRD_ST_Analyze_tracklets_cxx.so);

// Environment variables
//#define ENV_PI
#define ENV_ALEX
//#define ENV_PI_SVEN

void Macro_TRD_tracking(TString input_list = "List_data_ADC.txt")
{

    // First compile type root and then:
    // .L TRD_ST_Analyze_tracklets.cxx++
    // .L TRD_Kalman_Tracking.cxx++

    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(11);
    //gStyle->SetOptFit(1111);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFillColor(0);
    gStyle->SetPalette(27);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);


    gSystem ->Load("TRD_Kalman_Tracking_cxx.so");

    gSystem ->Load("TRD_ST_Analyze_tracklets_cxx.so");

    //------------------------------------
    // Define output file name and directory
    TString out_file_name = input_list;
    out_file_name += "_out.root";

    TString out_file_name_calib = input_list;
    out_file_name_calib += "_out_calib.root";

    TString input_dir  = "./Data/";
    TString output_dir = "./";
#if defined(ENV_PI)
    input_dir  = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/Calib_tracklets/";
    output_dir = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/ST_out/";
#endif

#if defined(ENV_PI_SVEN)
    input_dir  = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/Calib_tracklets/";
    output_dir = "/misc/alidata120/alice_u/hoppner/TRD_self_tracking/ST_out/";
#endif

#if defined(ENV_ALEX)
    input_dir  = "./Data/";
    output_dir = "./ST_out/";
#endif
    Int_t use_prim_vertex           = 1; // 0 = no primary vertex, 1 = primary vertex used
    Int_t KF_tracker                = 1; // Kalman filter tracker
    Int_t TF_tracker                = 0; // Tensorflow tracker

    Int_t graphics                  = 1; // 0 = no 3D graphics, 1 = 3D graphics (#define USEEVE in TRD_ST_Analyze_tracklets needs to be defined too)
    Int_t draw_tracklets_TPC_match  = 0; // Draw tracklets matched with TPC tracks
    Int_t draw_all_TPC_tracks       = 0; // Draw all TPC tracks
    Int_t draw_all_TRD_tracks       = 0; // Draw all TRD tracks
    Int_t draw_all_tracklets        = 0; // Draw all TRD tracklets
    Int_t draw_found_tracklets      = 1; // Draws tracklets found by tracker
    Int_t draw_matched_TPC_track    = 0; // Draw TPC to TRD matched TPC track
    Int_t draw_matched_TRD_track    = 0; // Draw TPC to TRD matched Kalman/TF track
    Int_t draw_secondary_vertices   = 0; // Draws tracks and secondary vertices

    //------------------------------------

    printf("TRD_ST_Analyze_tracklets started \n");
    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze(output_dir,out_file_name,graphics);
    TRD_ST_Analyze ->set_input_dir(input_dir);
    TRD_ST_Analyze ->Init_tree(input_list.Data());

    Long64_t N_Events = TRD_ST_Analyze ->get_N_Events();
    TH1D* h_layer_radii = TRD_ST_Analyze ->get_layer_radii_hist();
    //Long64_t event = 10;

    TRD_Kalman_Trackfinder kalid;
    kalid.set_layer_radii_hist(h_layer_radii);


    //------------------------------------------------------------
    //Tensorflow tracker
  //  TRD_TFTrackMaker tftracker;
    //------------------------------------------------------------



    // photon events: 88, 285, 29, 273
    // photon events thermal shield: 47,
    // pi0 event: 378
    // nuclear interaction event: 158, 168(!), 3741, 92, 328(!)
	//for(Long64_t event = 0; event < (Int_t) N_Events; event++) // 2,3
    for(Long64_t event = 0; event < 1; event++) // 2,3
    //for(Long64_t event = 0; event < 1; event++) // 2,3
    //for(Long64_t event = 0; event < N_Events; event++) // 2,3
    //for(Long64_t event = 0; event < 2; event++) // 2,3
    //Int_t event_plot = 555; // 168
    //for (Long64_t event = event_plot; event < (event_plot+1); event++) // 2,3   192
    {

        if (event != 0  &&  event % 50 == 0)
            cout << "." << flush;
        if (event != 0  &&  event % 500 == 0)
        {
            printf("event: %lld out of %lld, %4.2f%% total done \n",event,N_Events,((Double_t)event/(Double_t)N_Events)*100.0);
        }

        vector< vector<Ali_TRD_ST_Tracklets*> > tracker_found_tracklets_KF;
        vector<vector<Double_t>>                mHelices_tracker_KF;
        vector< vector<Ali_TRD_ST_Tracklets*> > tracker_found_tracklets_TF;
        vector<vector<Double_t>>                mHelices_tracker_TF;


        TRD_ST_Analyze ->Loop_event(event,graphics);
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;

        TRD_ST_Analyze ->Draw_event(event,graphics,draw_all_TPC_tracks,draw_all_tracklets);  // ->draws TPC tracks
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;


        TRD_ST_Analyze ->Do_TPC_TRD_matching(event,3.0,10.0,graphics*draw_tracklets_TPC_match); // last one is graphics  --> draws kalman TRD tracklets

        TRD_ST_Analyze ->set_TPC_helix_params(event);

        //------------------------------------------------------------
        //Kalman Filter tracker
        if(KF_tracker)
        {
            tracker_found_tracklets_KF = kalid.Kalman_Trackfind(TRD_ST_Analyze->Tracklets,TRD_ST_Analyze->Number_Tracklets,use_prim_vertex); // 0 = no primary vertex, 1 = primary vertex used

            vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks=TRD_ST_Analyze->matched_tracks; // TPC track matched tracklets
            vector< vector<Ali_TRD_ST_Tracklets*> > matched_beautiful_tracks;

            mHelices_tracker_KF = kalid.get_Kalman_helix_params();
            //printf("size of mHelices_tracker: %d \n",(Int_t)mHelices_tracker_KF.size());

            TRD_ST_Analyze ->set_Kalman_helix_params(mHelices_tracker_KF);
            TRD_ST_Analyze ->set_Kalman_TRD_tracklets(tracker_found_tracklets_KF);

            if(graphics && draw_found_tracklets) TRD_ST_Analyze ->Draw_Kalman_Tracklets(tracker_found_tracklets_KF); // Draws the Kalman matched TRD tracklets
            if(graphics && draw_all_TRD_tracks)  TRD_ST_Analyze ->Draw_Kalman_Helix_Tracks(-1,kGreen,3.0,500.0); // -1 -> all kalman tracks drawn
            TRD_ST_Analyze ->Match_kalman_tracks_to_TPC_tracks(graphics,draw_matched_TPC_track,draw_matched_TRD_track,kGreen+1);
        }
        //------------------------------------------------------------


        //------------------------------------------------------------
        //TensorFlow tracker
   /*     if(TF_tracker)
        {
            tftracker.Trackfind(TRD_ST_Analyze->Tracklets,TRD_ST_Analyze->Number_Tracklets);
            mHelices_tracker_TF = tftracker.get_Kalman_helix_params();
            tracker_found_tracklets_TF = tftracker.get_Tracklets();

            TRD_ST_Analyze ->set_Kalman_helix_params(mHelices_tracker_TF);
            TRD_ST_Analyze ->set_Kalman_TRD_tracklets(tracker_found_tracklets_TF);

            if(graphics && draw_found_tracklets) TRD_ST_Analyze ->Draw_Kalman_Tracklets(tracker_found_tracklets_TF); // Draws the Kalman matched TRD tracklets
            if(graphics && draw_all_TRD_tracks)  TRD_ST_Analyze ->Draw_Kalman_Helix_Tracks(-1,kRed,3.0,500.0); // -1 -> all kalman tracks drawn
            TRD_ST_Analyze ->Match_kalman_tracks_to_TPC_tracks(graphics,draw_matched_TPC_track,draw_matched_TRD_track,kRed+1);
        }
        //------------------------------------------------------------

 */      
        TRD_ST_Analyze ->Calibrate(graphics);
        TRD_ST_Analyze ->Calc_Kalman_efficiency();


        Int_t found_good_AP_vertex_TPC = TRD_ST_Analyze ->Calculate_secondary_vertices(graphics*draw_secondary_vertices,0); // (0 = no graphics), (0 = TRD, 1 = TPC)
        Int_t found_good_AP_vertex_TRD = TRD_ST_Analyze ->Calculate_secondary_vertices(graphics*draw_secondary_vertices,1); // (0 = no graphics), (0 = TRD, 1 = TPC)
        //if(found_good_AP_vertex) printf(" ----> Good AP vertex found in event: %lld \n",event);

        //vector< vector<Ali_TRD_ST_Tracklets*> > tracker_found_tracklets=kalid.found_tracks;
        //vector<Double_t> track_accuracy;

        //Ali_TRD_ST_Tracklets* last_tracklet;

    }
	
    TRD_ST_Analyze ->Draw_n_Save_Calibration(output_dir,out_file_name_calib);
    TRD_ST_Analyze ->Write();

}