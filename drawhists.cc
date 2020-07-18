//#include "TRD_ST_Analyze_tracklets.h"
//#include "TRD_Kalman_Tracking.h"
//#include "TRD_Kalman_Tracking.cxx"

R__LOAD_LIBRARY(TRD_Kalman_Tracking_cxx.so);
R__LOAD_LIBRARY(TRD_ST_Analyze_tracklets_cxx.so);

Bool_t fitting_track(Ali_TRD_ST_Tracklets* a,Ali_TRD_ST_Tracklets* b){
    //direction is not ok

    TVector3 a_angle = (TVector3){a->get_TV3_dir()[0], a->get_TV3_dir()[1], 0};
    TVector3 b_angle = (TVector3){b->get_TV3_dir()[0], b->get_TV3_dir()[1], 0};

    if (abs(a_angle.Angle(b_angle)*TMath::RadToDeg() )>7)
        return 0;
    if (abs(a->get_TV3_dir().Angle(b->get_TV3_dir())*TMath::RadToDeg() )>30)
        return 0;
    //position is not ok
    Double_t cylind_rad_dist=sqrt(b->get_TV3_offset()[0]*b->get_TV3_offset()[0] +b->get_TV3_offset()[1]*b->get_TV3_offset()[1]) -sqrt(a->get_TV3_offset()[0]*a->get_TV3_offset()[0] +a->get_TV3_offset()[1]*a->get_TV3_offset()[1]);
    Double_t cylind_rad_len=sqrt(a->get_TV3_dir()[0]*a->get_TV3_dir()[0] +a->get_TV3_dir()[1]*a->get_TV3_dir()[1]);
    TVector3 z=a->get_TV3_offset()+a->get_TV3_dir()*(cylind_rad_dist/cylind_rad_len );
    /*if (abs((a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] ) -b->get_TV3_offset())[1])>1)
     return 0;
     if (abs((a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] ) -b->get_TV3_offset())[2])>4)
     return 0;*/

    if ((z-b->get_TV3_offset())[0]>2.)
        return 0;

    if ((z-b->get_TV3_offset())[1]>2.)
        return 0;
    if ((z-b->get_TV3_offset())[2]>10.)
        return 0;

    return 1;
}
void drawhists(TString input_list = "List_data_ADC.txt")
{

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
    TString out_dir = "./";
    //------------------------------------

    TH1F *histo = new TH1F("histogram","efficiency Kalman Trackfinder",20,0,1.2);

    TH1F *histo2 = new TH1F("histogram","purity Kalman Trackfinder",20,0,1.2);


    printf("TRD_ST_Analyze_tracklets started \n");
    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze(out_dir,out_file_name);
    TRD_ST_Analyze ->Init_tree(input_list.Data());

    //TRD_ST_Analyze ->create_output_file(out_dir,out_file_name);

    Long64_t N_Events = TRD_ST_Analyze ->get_N_Events();
    TH1D* h_layer_radii = TRD_ST_Analyze ->get_layer_radii_hist();
    //Long64_t event = 10;

    Int_t graphics = 0;
    TRD_Kalman_Trackfinder kalid;
    kalid.set_layer_radii_hist(h_layer_radii);

    //for (Long64_t event = 0; event < N_Events; event++) // 2,3
    for (Long64_t event = 0; event < 1000; event++) // 2,3
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
        vector< vector<Ali_TRD_ST_Tracklets*> > kalman_found_tracks = kalid.Kalman_Trackfind(TRD_ST_Analyze->Tracklets,TRD_ST_Analyze->Number_Tracklets,1); // 0 = no primary vertex, 1 = primary vertex used
        if(graphics) TRD_ST_Analyze ->Draw_Kalman_Tracks(kalman_found_tracks);


        vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks=TRD_ST_Analyze->matched_tracks; // TPC track matched tracklets
        vector< vector<Ali_TRD_ST_Tracklets*> > matched_beautiful_tracks;

        vector<vector<Double_t>> mHelices_kalman = kalid.get_Kalman_helix_params();
        //printf("size of mHelices_kalman: %d \n",(Int_t)mHelices_kalman.size());

        TRD_ST_Analyze ->set_Kalman_helix_params(mHelices_kalman);
        TRD_ST_Analyze ->set_Kalman_TRD_tracklets(kalman_found_tracks);


        //if(graphics) TRD_ST_Analyze ->Draw_Kalman_Helix_Tracks(-1,kRed); // -1 -> all kalman tracks drawn

        TRD_ST_Analyze ->Match_kalman_tracks_to_TPC_tracks(graphics);

#endif
        //TRD_ST_Analyze ->Calculate_secondary_vertices(graphics); // 0 = no graphics

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