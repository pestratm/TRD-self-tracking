
R__LOAD_LIBRARY(TRD_Kalman_Tracking_cxx.so);
R__LOAD_LIBRARY(TRD_ST_Analyze_tracklets_cxx.so);

// Environment variables
//#define ENV_PI
#define ENV_ALEX
//#define ENV_PI_SVEN

void drawhists(TString input_list = "List_data_ADC.txt")
{

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
    Int_t graphics = 1; // 0 = no 3D graphics, 1 = 3D graphics (#define USEEVE in TRD_ST_Analyze_tracklets needs to be defined too)
    Int_t use_prim_vertex = 0; // 0 = no primary vertex, 1 = primary vertex used
    //------------------------------------

    TH1F *histo = new TH1F("histogram","efficiency Kalman Trackfinder",20,0,1.2);
    TH1F *histo2 = new TH1F("histogram","purity Kalman Trackfinder",20,0,1.2);

    printf("TRD_ST_Analyze_tracklets started \n");
    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze(output_dir,out_file_name,graphics);
    TRD_ST_Analyze ->set_input_dir(input_dir);
    TRD_ST_Analyze ->Init_tree(input_list.Data());

    Long64_t N_Events = TRD_ST_Analyze ->get_N_Events();
    TH1D* h_layer_radii = TRD_ST_Analyze ->get_layer_radii_hist();
    //Long64_t event = 10;

    TRD_Kalman_Trackfinder kalid;
    kalid.set_layer_radii_hist(h_layer_radii);

    // photon events: 88
    // nuclear interaction event: 158

    //for (Long64_t event = 0; event < N_Events; event++) // 2,3
    for (Long64_t event = 225; event < 226; event++) // 2,3   192
    {

        if (event != 0  &&  event % 50 == 0)
            cout << "." << flush;
        if (event != 0  &&  event % 500 == 0)
        {
            printf("event: %lld out of %lld, %4.2f%% total done \n",event,N_Events,((Double_t)event/(Double_t)N_Events)*100.0);
        }


        TRD_ST_Analyze ->Loop_event(event,graphics);
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;


        //TRD_ST_Analyze ->Draw_event(event);  // ->draws TPC tracks
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;
        TRD_ST_Analyze ->Do_TPC_TRD_matching(event,3.0,10.0,graphics); // last one is graphics  --> draws kalman TRD tracklets
        //TRD_ST_Analyze ->Do_TPC_TRD_matching_allEvents(3.0,10.0);


#if 1
        vector< vector<Ali_TRD_ST_Tracklets*> > kalman_found_tracks = kalid.Kalman_Trackfind(TRD_ST_Analyze->Tracklets,TRD_ST_Analyze->Number_Tracklets,use_prim_vertex); // 0 = no primary vertex, 1 = primary vertex used
        if(graphics) TRD_ST_Analyze ->Draw_Kalman_Tracks(kalman_found_tracks);


        vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks=TRD_ST_Analyze->matched_tracks; // TPC track matched tracklets
        vector< vector<Ali_TRD_ST_Tracklets*> > matched_beautiful_tracks;

        vector<vector<Double_t>> mHelices_kalman = kalid.get_Kalman_helix_params();
        //printf("size of mHelices_kalman: %d \n",(Int_t)mHelices_kalman.size());

        TRD_ST_Analyze ->set_Kalman_helix_params(mHelices_kalman);
        TRD_ST_Analyze ->set_Kalman_TRD_tracklets(kalman_found_tracks);


        if(graphics) TRD_ST_Analyze ->Draw_Kalman_Helix_Tracks(-1,kRed); // -1 -> all kalman tracks drawn

        TRD_ST_Analyze ->Match_kalman_tracks_to_TPC_tracks(graphics);
        //TRD_ST_Analyze ->Match_kalman_tracks_to_TPC_tracks(1);

#endif
        Int_t found_good_AP_vertex = TRD_ST_Analyze ->Calculate_secondary_vertices(graphics); // 0 = no graphics
        //if(found_good_AP_vertex) printf(" ----> Good AP vertex found in event: %lld \n",event);

        //vector< vector<Ali_TRD_ST_Tracklets*> > kalman_found_tracks=kalid.found_tracks;
        //vector<Double_t> track_accuracy;

        //Ali_TRD_ST_Tracklets* last_tracklet;

    }
    //gRandom->Rndm();

    TRD_ST_Analyze ->Plot_AP();
    TRD_ST_Analyze ->Plot_pT_TPC_vs_Kalman();
    TRD_ST_Analyze ->Write();

    TCanvas * c1= new TCanvas("c1", "fitted data",5,5,800,600);

    TColor *c3 = new TColor(9003,0,1,0);
    histo->SetLineColor(9003);

    histo->SetTitle("Efficiency Of Tracklet Detection");
    histo->SetXTitle("efficiency");
    histo->SetYTitle("count");

    histo->Draw("histo");
    TCanvas * c2= new TCanvas("c2", "fitted data",5,5,800,600);

    //TColor *c3 = new TColor(9003,0,1,0);


    histo2->SetLineColor(9003);

    histo2->SetTitle("Purity Of Tracklet Detection");
    histo2->SetXTitle("purity");
    histo2->SetYTitle("count");

    histo2->Draw("histo");


    //TRD_ST_Analyze ->Draw_Kalman_Tracks(kalid.found_tracks,kalid.nbr_tracks);

}