
#include "TRD_ST_Analyze_tracklets.h"

#if 1
// distance between circle and point 2D
static Double_t distance_circ_point_2D(Double_t x,Double_t y,Double_t *p)
{
    // D = |sqrt((xp-x0)^2 + (yp-y0)^2)-r|

    Double_t d2 = TMath::Abs(TMath::Sqrt(TMath::Power(x-p[0],2.0)+TMath::Power(y-p[1],2.0))-p[2]);
    return d2;
}


// function to be minimized - for circle fit
static void sum_distance_circ_point_2D(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t )
{

    sum = 0;

    for (Int_t i_layer = 0; i_layer < 6; ++i_layer)
    {
        //printf("i_layer: %d, vec_Dt_digit_pos_cluster[0]: %4.3f, vec_Dt_digit_pos_cluster[1]: %4.3f \n",i_layer,vec_Dt_digit_pos_cluster[i_layer][0],vec_Dt_digit_pos_cluster[i_layer][1]);
        
        if(vec_Dt_digit_pos_cluster[i_layer][0] == -999.9) continue;
        Double_t d = 100*distance_circ_point_2D(vec_Dt_digit_pos_cluster[i_layer][0],vec_Dt_digit_pos_cluster[i_layer][1],par);
        sum += d;
        //printf("i_layer: %d, vec_Dt_digit_pos_cluster[6][i_layer][0]: %4.3f, vec_Dt_digit_pos_cluster[6][i_layer][1]: %4.3f, sum: %4.8f \n",i_layer,vec_Dt_digit_pos_cluster[6][i_layer][0],vec_Dt_digit_pos_cluster[6][i_layer][1],sum);
        //printf("i_layer: %d, d: %4.8f \n",i_layer,d);

    }
}
#endif


//----------------------------------------------------------------------------------------
Ali_TRD_ST_Analyze::Ali_TRD_ST_Analyze(TString out_dir, TString out_file_name, Int_t graphics)
{
    //layer_radii_file = TFile::Open("./TRD_layer_radii.root");
    //h_layer_radii_det = (TH1D*)layer_radii_file ->Get("h_layer_radii_det");
	
	vec_tp_pvdca_vs_sector = new TProfile("vec_th1d_pvdca_vs_sector","vec_th1d_pvdca_vs_sector",90,0,89);
    vec_TH2D_pvdca_vs_sector = new TH2D("vec_th2d_pvdca_vs_sector","vec_th2d_pvdca_vs_sector",90,0,89,100,0,50);
 
    TPC_single_helix = new Ali_Helix();
    vec_h2D_pT_vs_TPC_TRD_residuals.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_h2D_pT_vs_TPC_TRD_residuals[i_det].resize(3); // x,y,z
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            HistName = "vec_h2D_pT_vs_TPC_TRD_residuals_det";
            HistName += i_det;
            HistName += "_xyz";
            HistName += i_xyz;
            vec_h2D_pT_vs_TPC_TRD_residuals[i_det][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),200,-10.0,10.0,10,-5,5);
        }
    }

    TFile* file_TRD_geometry = TFile::Open("./TRD_Geometry.root");
    TList* list_geom = (TList*)file_TRD_geometry->Get("TRD_Digits_output");
    h2D_TRD_det_coordinates = (TH2D*)list_geom ->FindObject("h2D_TRD_det_coordinates");
    h_layer_radii_det = new TH1D("h_layer_radii_det","h_layer_radii_det",540,0,540);

#if defined(USEEVE)
    TEveP_TRD_det_origin = new TEvePointSet();
#endif


    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        Double_t x_val = h2D_TRD_det_coordinates ->GetBinContent(i_det+1,1);
        Double_t y_val = h2D_TRD_det_coordinates ->GetBinContent(i_det+1,2);
        Double_t z_val = h2D_TRD_det_coordinates ->GetBinContent(i_det+1,3);

        Double_t radius_align = TMath::Sqrt(x_val*x_val + y_val*y_val);
        h_layer_radii_det ->SetBinContent(i_det+1,radius_align);

#if defined(USEEVE)
        //TEveP_TRD_det_origin ->SetPoint(i_det,x_val,y_val,z_val);
#endif
    }


    HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    outputfile = new TFile(HistName.Data(),"RECREATE");

    HistName = out_dir;
    HistName += "/";
    HistName += "gain_";
    HistName += out_file_name;
    
    out_gain = new TFile(HistName.Data(),"RECREATE");

    //------------------------------------------------
    outputfile ->cd();
    // TRD self tracking output data containers
    TRD_ST_Tracklet_out   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track_out  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event_out      = new Ali_TRD_ST_Event();

    //TRD self tracking output data containers for photons and nuclear interactions
    TRD_Self_Event = new Ali_TRD_Self_Event();
    TRD_Photon = new Ali_TRD_Photon();
    TRD_Nuclear_interaction = new Ali_TRD_Nuclear_interaction();
    TRD_Kalman_track = new Ali_Kalman_Track();
    TPC_track = new Ali_TPC_Track();

    //create new 


    //NT_secondary_vertices = new TNtuple("NT_secondary_vertices","NT_secondary_vertices Ntuple","x:y:z:ntracks:pT_AB:qpT_A:qpT_B:AP_pT:AP_alpha:dcaTPC:pathTPC:InvM:Eta:Phi:GlobEv:dotprod:TPCdEdx_A:dca_TPC_A:p_TPC_A:TPCdEdx_B:dca_TPC_B:p_TPC_B:InvMK0s:dcaAB:InvML:InvMaL");
    //NT_secondary_vertices ->SetAutoSave( 5000000 );

    //NT_secondary_vertex_cluster = new TNtuple("NT_secondary_vertex_cluster","NT_secondary_vertex_cluster Ntuple","x:y:z:nvertices:dcaTPC:tof:trklength:dEdx:dcaprim:pT:mom:layerbitmap");
    //NT_secondary_vertex_cluster ->SetAutoSave( 5000000 );

    Tree_TRD_ST_Event_out  = NULL;
    Tree_TRD_ST_Event_out  = new TTree("Tree_TRD_ST_Event_out" , "TRD_ST_Events_out" );
    Tree_TRD_ST_Event_out  ->Branch("Tree_TRD_ST_Event_branch_out"  , "TRD_ST_Event_out", TRD_ST_Event_out );
    Tree_TRD_ST_Event_out  ->SetAutoSave( 5000000 );

    Tree_TRD_Self_Event_out  = NULL;
    Tree_TRD_Self_Event_out  = new TTree("Tree_TRD_Self_Event" , "TRD_Self_Events" );
    Tree_TRD_Self_Event_out  ->Branch("Tree_TRD_Self_Event_branch"  , "TRD_Self_Event", TRD_Self_Event );
    Tree_TRD_Self_Event_out  ->SetAutoSave( 5000000 );
    //------------------------------------------------


    // constructor
#if defined(USEEVE)
    if(graphics)
    {
        TEveManager::Create();

        TEveP_primary_vertex = new TEvePointSet();

        TPL3D_helix = new TEveLine();
        TEveLine_beam_axis = new TEveLine();
        TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-650.0);
        TEveLine_beam_axis ->SetNextPoint(0.0,0.0,650.0);
        TEveLine_beam_axis ->SetName("beam axis");
        TEveLine_beam_axis ->SetLineStyle(1);
        TEveLine_beam_axis ->SetLineWidth(4);
        TEveLine_beam_axis ->SetMainAlpha(0.7);
        TEveLine_beam_axis ->SetMainColor(kBlue);
        gEve->AddElement(TEveLine_beam_axis);

        vec_TEveLine_tracklets.resize(6); // layers
        vec_TEveLine_tracklets_match.resize(6); // layers

        TEveP_TRD_det_origin ->SetMarkerSize(2);
        TEveP_TRD_det_origin ->SetMarkerStyle(20);
        TEveP_TRD_det_origin ->SetMainColor(kMagenta+1);
        //gEve->AddElement(TEveP_TRD_det_origin);

        Double_t tube_R = 1.0;
        Int_t TPC_color = kGray;

        auto TPC_A = new TEveArrow();
        TPC_A ->SetVector(0.0,0.0,250.0);
        TPC_A ->SetOrigin(0.0,0.0,-250.0);
        TPC_A ->SetConeL(0);
        //TPC_A ->SetConeR(50.0);
        TPC_A ->SetDrawQuality(100);
        TPC_A->SetMainColor(TPC_color);
        TPC_A ->SetMainAlpha(0.6);
        TPC_A->SetTubeR(tube_R);
        TPC_A->SetPickable(kTRUE);

        auto TPC_B = new TEveArrow();
        TPC_B ->SetVector(0.0,0.0,-250.0);
        TPC_B ->SetOrigin(0.0,0.0,250.0);
        TPC_B ->SetConeL(0);
        //TPC_B ->SetConeR(50.0);
        TPC_B ->SetDrawQuality(100);
        TPC_B->SetMainColor(TPC_color);
        TPC_B ->SetMainAlpha(0.6);
        TPC_B->SetTubeR(tube_R);
        TPC_B->SetPickable(kTRUE);

        //gEve->AddElement(TPC_A);
        //gEve->AddElement(TPC_B);


        TEveP_beamA = new TEvePointSet();
        TEveP_beamB = new TEvePointSet();
    }
#endif

    th1d_TRD_layer_radii = new TH1D("th1d_TRD_layer_radii","th1d_TRD_layer_radii",900,250,400.0);
    vec_th1d_TRD_layer_radii_det.resize(540);
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        HistName = "vec_th1d_TRD_layer_radii_det_";
        HistName += TRD_detector;
        vec_th1d_TRD_layer_radii_det[TRD_detector] = new TH1D(HistName.Data(),HistName.Data(),900,250,400.0);
    }


    tp_efficiency_matching_vs_pT = new TProfile("tp_efficiency_matching_vs_pT","tp_efficiency_matching_vs_pT",100,0,10);
    tp_efficiency_all_vs_pT      = new TProfile("tp_efficiency_all_vs_pT","tp_efficiency_all_vs_pT",100,0,10);
    vec_h2D_delta_pT_all_vs_pT.resize(6); // number of matched tracklets
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        HistName = "vec_h2D_delta_pT_all_vs_pT_";
        HistName += i_layer;
        vec_h2D_delta_pT_all_vs_pT[i_layer] = new TH2D(HistName.Data(),HistName.Data(),200,-10,10,400,-10,10);
    }



    //--------------------------
    // Open histogram which defines good and bad chambers
    //TFile* file_TRD_QA = TFile::Open("./Data/chamber_QC.root");
    //h_good_bad_TRD_chambers = (TH1D*)file_TRD_QA ->Get("all_defects_hist");

    h_good_bad_TRD_chambers = new TH1I("h_good_bad_TRD_chambers","h_good_bad_TRD_chambers",540,0,540);
    TFile* file_TRD_QA_flags = TFile::Open("./Data/chamber_QC_flags.root");
    vector<int> *t_flags;
    file_TRD_QA_flags ->GetObject("QC_flags", t_flags);

    // Its a 3 digit binary number. LSB is ADC good = 0 or bad = 1, next bit is anode HV good = 0, or bad = 1, and last bit is drift HV
    // so a 3 means that the ADC and the anode HV was bad, but the drift HV was okay

    // LSB = official QA, bit 1 = no fit, bit 2 = anode HV defect, bit 3 = drift HV defect, bit 4 = adc defect

    // number   adc defect   drift HV defect   anode HD defect    no fit   official QA
    //   0          0               0                0               0          0         --> all good
    //   1          0               0                0               0          1         --> official QA bad, rest good
    //  ...
    //   31         1               1                1               1          1         --> all bad

    Int_t i_chamber = 0;
    for(vector<int>::iterator it = t_flags->begin(); it != t_flags->end(); ++it)
    {
        //cout << "chamber: " << i_chamber << ", it: "  << *it << ", " << t_flags->at(i_chamber) << endl;
        h_good_bad_TRD_chambers ->SetBinContent(i_chamber+1,t_flags->at(i_chamber));
        i_chamber++;
    }
    //--------------------------

    //--------------------------
    // Load TRD geometry
    TFile* file_TRD_geom = TFile::Open("./Data/TRD_Geom.root");
    vec_TH1D_TRD_geometry.resize(3); // x,y,z
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TRD_geometry[i_xyz].resize(8); // 8 vertices
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            HistName = "vec_TH1D_TRD_geometry_xyz_";
            HistName += i_xyz;
            HistName += "_V";
            HistName += i_vertex;
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] = (TH1D*)file_TRD_geom->Get(HistName.Data());
        }

    }
#if defined(USEEVE)
    vec_eve_TRD_detector_box.resize(540);
#endif
    Int_t color_flag_QC[32];
    for(Int_t i_QC_flag = 0; i_QC_flag < 32; i_QC_flag++)
    {
        color_flag_QC[i_QC_flag] = kCyan;

        Int_t k_bit = 1; // fit
        Int_t bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // no fit
        {
            color_flag_QC[i_QC_flag] = kPink;
        }

        k_bit = 4; // ADC value
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // ADC low
        {
            color_flag_QC[i_QC_flag] = kMagenta;
        }

        k_bit = 2; // anode HV
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // anode HV low
        {
            color_flag_QC[i_QC_flag] = kYellow;
        }

        k_bit = 3; // drift HV bit
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // drift HV defect
        {
            color_flag_QC[i_QC_flag] = kOrange;
        }

        k_bit = 0; // official QA
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // official QA bad
        {
            color_flag_QC[i_QC_flag] = kRed;
        }
    }
    color_flag_QC[31] = kRed;
    //= {kCyan,kPink,kMagenta,kMagenta+2,kOrange,kOrange+2,kRed,kRed+2};

#if defined(USEEVE)
    if(graphics)
    {
        for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
        {
            vec_eve_TRD_detector_box[TRD_detector] = new TEveBox;

            HistName = "TRD_box_";
            HistName += TRD_detector;
            vec_eve_TRD_detector_box[TRD_detector] ->SetName(HistName.Data());
            Int_t flag_QC = h_good_bad_TRD_chambers ->GetBinContent(TRD_detector+1);
            if(!flag_QC) // chamber is OK flagged by QA
            {
                vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(kCyan);
                vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(95); // the higher the value the more transparent
            }
            else // bad chamber
            {
                vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(color_flag_QC[flag_QC]);
                vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(85); // the higher the value the more transparent
            }
            for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
            {
                Double_t arr_pos_glb[3] = {vec_TH1D_TRD_geometry[0][i_vertex]->GetBinContent(TRD_detector),vec_TH1D_TRD_geometry[1][i_vertex]->GetBinContent(TRD_detector),vec_TH1D_TRD_geometry[2][i_vertex]->GetBinContent(TRD_detector)};
                vec_eve_TRD_detector_box[TRD_detector]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);
            }

            gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
        }
    }
#endif
    //--------------------------


    TH2D_AP_plot          = new TH2D("TH2D_AP_plot","TH2D_AP_plot",200,-2.0,2.0,400,-0.1,4.0);
    TH2D_pT_TPC_vs_Kalman = new TH2D("TH2D_pT_TPC_vs_Kalman","TH2D_pT_TPC_vs_Kalman",2000,-10.0,10.0,2000,-10.0,10.0);
    vec_TH2D_AP_plot_radius.resize(N_AP_radii); // secondary vertex radius
    for(Int_t i_radius = 0; i_radius < N_AP_radii; i_radius++)
    {
        HistName = "vec_TH2D_AP_plot_radius_";
        HistName += i_radius;
        vec_TH2D_AP_plot_radius[i_radius] = new TH2D(HistName.Data(),HistName.Data(),200,-2.0,2.0,400,-0.1,4.0);
    }

    vec_TH2D_pT_TPC_vs_Kalman.resize(N_pT_resolution);
    vec_TH2D_one_over_pT_TPC_vs_Kalman.resize(N_pT_resolution);
    for(Int_t i_pt_res = 0; i_pt_res < N_pT_resolution; i_pt_res++)
    {
        HistName = "vec_TH2D_pT_TPC_vs_Kalman_";
        HistName += i_pt_res;
        vec_TH2D_pT_TPC_vs_Kalman[i_pt_res] = new TH2D(HistName.Data(),HistName.Data(),2000,-10.0,10.0,2000,-10.0,10.0);

        HistName = "vec_TH2D_one_over_pT_TPC_vs_Kalman_";
        HistName += i_pt_res;
        vec_TH2D_one_over_pT_TPC_vs_Kalman[i_pt_res] = new TH2D(HistName.Data(),HistName.Data(),2000,-10.0,10.0,2000,-10.0,10.0);
    }

    TFile* file_TRD_centers = TFile::Open("./TV3_TRD_center.root");
    vec_TV3_TRD_center.resize(540);
    Int_t i_vec;

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_TRD_center[i_det].resize(3);
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {

            i_vec = i_det*3 +i_xyz+1;

            vec_TV3_TRD_center[i_det][i_xyz] = (TVector3*)file_TRD_centers -> Get(Form("TVector3;%d",i_vec));
            //vec_TV3_TRD_center[i_det][i_xyz] = file_TRD_geometry ->GetObject(Form("TVector3;%d",i_vec), TVector3);
        }
    }

    vec_TV3_TRD_offset.clear();
    TFile* file_TRD_offsets = TFile::Open("./TRD_center_offset.root");
    vec_TV3_TRD_offset.resize(540);
    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_TRD_offset[i_det] = (TVector3*)file_TRD_offsets ->Get(Form("TVector3;%d",i_det+1));
    //printf("i_det: %d, vec_TV3_TRD_offset: {%4.3f, %4.3f, %4.3f} \n", i_det,(*vec_TV3_TRD_offset[i_det]).X(),(*vec_TV3_TRD_offset[i_det]).Y(),(*vec_TV3_TRD_offset[i_det]).Z());

    }

    vec_trd_TRD_pp_geom.clear();
    TFile* TRD_pp_geom = TFile::Open("./geom_out.root");
    vec_trd_TRD_pp_geom.resize(540);
    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        TGraph* geom = (TGraph*)TRD_pp_geom ->Get(Form("geom_param_%d;1",i_det));
        vec_trd_TRD_pp_geom[i_det].resize(6);
        for(Int_t i_param = 0; i_param < 6 ; i_param++)
        {
            vec_trd_TRD_pp_geom[i_det][i_param] = geom ->GetY()[i_param+1]; // n_col n_row col_start col_end row_start row_end
        }

        printf("i_det: %d, i_det from tg: %4.3f, n_cols: %4.3f, n_rows: %4.3f, col start: %4.3f, col end: %4.3frow start: %4.3f, row end: %4.3f \n", i_det,geom ->GetY()[0],vec_trd_TRD_pp_geom[i_det][0],vec_trd_TRD_pp_geom[i_det][1],vec_trd_TRD_pp_geom[i_det][2],vec_trd_TRD_pp_geom[i_det][3],vec_trd_TRD_pp_geom[i_det][4],vec_trd_TRD_pp_geom[i_det][5]);

    }

	       
    TFile* file_bethe_bloch = TFile::Open("Data/Bethe_Bloch.root");
    tg_bethe_bloch.resize(5);
    tg_bethe_bloch[0] = (TGraph*)file_bethe_bloch ->Get("TSPLINE3_DATA_ELECTRON_LHC10H_PASS2_PBPB_MEAN");
    tg_bethe_bloch[1] = (TGraph*)file_bethe_bloch ->Get("TSPLINE3_DATA_MUON_LHC10H_PASS2_PBPB_MEAN");
    tg_bethe_bloch[2] = (TGraph*)file_bethe_bloch ->Get("TSPLINE3_DATA_PION_LHC10H_PASS2_PBPB_MEAN");
    tg_bethe_bloch[3] = (TGraph*)file_bethe_bloch ->Get("TSPLINE3_DATA_KAON_LHC10H_PASS2_PBPB_MEAN");
    tg_bethe_bloch[4] = (TGraph*)file_bethe_bloch ->Get("TSPLINE3_DATA_PROTON_LHC10H_PASS2_PBPB_MEAN");

    tp_gain = new TProfile("tp_gain","tp_gain",540,0,539,0,10);
    tp_gain_uncor = new TProfile("tp_gain_uncor","tp_gain_uncor",540,0,539,0,10);
    th2d_gain = new TH2D("th2d_gain","th2d_gain",540,0,539,50,0,20);

    vec_tp2d_gain_vs_xz.clear();
    vec_tp2d_gain_vs_xz.resize(540);

    // n_col n_row col_start col_end row_start row_end

    for(Int_t i_det = 0; i_det < 540; i_det++) 
    {
        Double_t n_col = vec_trd_TRD_pp_geom[i_det][0]; // == x == 144
        Double_t n_row = vec_trd_TRD_pp_geom[i_det][1]; // == z

        if (n_row == 0 || n_col == 0) 
        {
            vec_tp2d_gain_vs_xz[i_det] = new TProfile2D(Form("vec_tp2d_gain_vs_xz%d",i_det),Form("vec_tp2d_gain_vs_xz%d",i_det),144,-60.0,60.0,16,-60.0,60.0);
        }

        else
        {
            Double_t half_x = fabs(vec_trd_TRD_pp_geom[i_det][3] - vec_trd_TRD_pp_geom[i_det][2])/2 + 0.5; //x == col
            //Double_t half_x = 60.0;
            Double_t half_z = fabs(vec_trd_TRD_pp_geom[i_det][5] - vec_trd_TRD_pp_geom[i_det][4])/2 + 5.0; //z == row

            vec_tp2d_gain_vs_xz[i_det] = new TProfile2D(Form("vec_tp2d_gain_vs_xz%d",i_det),Form("vec_tp2d_gain_vs_xz%d",i_det),n_col,-half_x,half_x,n_row,-half_z,half_z);
        }
    }

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::fHelixAtoPointdca(TVector3 space_vec, Ali_Helix* helixA, Float_t &pathA, Float_t &dcaAB)
{
    // V1.1
    Float_t pA[2] = {100.0,-100.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    Double_t helix_point[3];
    for(Int_t r = 0; r < 2; r++)
    {
        helixA ->Evaluate(pA[r],helix_point);  // 3D-vector of helixB point at path pB[r]
        testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.01 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
        //    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pA[1]-pA[0])*scale; // the next length interval
            pA[0]     = pA[1] + scale_length; // the new path

            helixA ->Evaluate(pA[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            distarray[0] = (testA-space_vec).Mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pA[0]-pA[1])*scale;
            pA[1]     = pA[0] + scale_length;
            helixA ->Evaluate(pA[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[1] = (testA-space_vec).Mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathA = pA[1];
        dcaAB = distarray[1];
    }
    //cout << "pathA = " << pathA << ", dcaAB = " << dcaAB << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,
                           Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
    // (x1,y1) -> center of circle 1, r1 -> radius of circle 1
    // (x2,y2) -> center of circle 2, r2 -> radius of circle 2
    // (x1_c,y1_c) -> crossing point 1
    // (x2_c,y2_c) -> crossing point 2
    // Solution see Wikipedia

    Double_t s = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));  // distance between the center of the circles

    x1_c = 0;
    y1_c = 0;
    x2_c = 0;
    y2_c = 0;

    //printf("s: %4.3f, r1+r2: %4.3f, pos1: {%4.3f, %4.3f}, pos2: {%4.3f, %4.3f}, r1: %4.3f, r2: %4.3f \n",s,r1+r2,x1,y1,x2,y2,r1,r2);

    if(x1 != x2 && y1 != y2 && s < (r1+r2))
    {
        Double_t m  = (x1-x2)/(y2-y1);
        Double_t n  = (-r2*r2 + r1*r1 + y2*y2 - y1*y1 + x2*x2 - x1*x1)/(2.0*(y2-y1));
        Double_t p  = (2.0*(-x1 + m*(n-y1)))/(1.0 + m*m);
        Double_t q  = (x1*x1 + (n-y1)*(n-y1) -r1*r1)/(1.0 + m*m);
        Double_t p2 = (p/2.0)*(p/2.0);


        if(p2 >= q)
        {
            x1_c = (-p/2.0) + sqrt(p2 - q);
            x2_c = (-p/2.0) - sqrt(p2 - q);
            y1_c = m*x1_c + n;
            y2_c = m*x2_c + n;
            return 1;
        }
        else return 0;
    }
    else return 0;

} 
//----------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::fCircle_Interception(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,
                                               Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
	Double_t dif_x = -(x1 - x2);
	Double_t dif_y = -(y1 - y2);
	Double_t dist=TMath::Sqrt( (dif_x*dif_x) + (dif_y*dif_y) );
        if(dist==0) // circles overlap -> most likely identical helices
        {
            x1_c = 0.0;
            y1_c = 0.0;
            x2_c = 0.0;
            y2_c = 0.0;
            return 0;
        };
        Double_t dist_inv =1/dist;
	
        if((dist<(r1+r2))&& (dist> TMath::Abs(r1-r2)))
        {
	//2 intersections
//		8<i<i<i<iii<i<iiii>ii


		Double_t x_loc 		=((dist*dist +r1*r1 -r2*r2)*(0.5*dist_inv));
		Double_t y_loc		=TMath::Sqrt(r1*r1 -x_loc*x_loc);
		
                x1_c	=(x1*dist+x_loc*dif_x -y_loc*dif_y)*dist_inv;
                y1_c	=(y1*dist+x_loc*dif_y +y_loc*dif_x)*dist_inv;
		
                x2_c	=(x1*dist+x_loc*dif_x +y_loc*dif_y)*dist_inv;
                y2_c	=(y1*dist+x_loc*dif_y -y_loc*dif_x)*dist_inv;

                return 1;
		}
		else if (dist<= TMath::Abs(r1-r2))
		{
			Double_t normrad	=(r1-dist+r2)*dist_inv*0.5;
            
            //Double_t normrad	= (r1 + 0.5*(dist - r1 - r2))*dist_inv;
			if(r2>r1)
			{	
                x1_c		=x1 - normrad*dif_x;
                y1_c		=y1 - normrad*dif_y;
			}
			else
			{
				x1_c		=x2 + normrad*dif_x;
                y1_c		=y2 + normrad*dif_y;
			
			}	
                x2_c		=x1_c;
                y2_c		=y1_c;
        
            	return 2;
		}	
        else
        {
	//no intersection (maybe 1)
            Double_t normrad	=(r1+dist-r2)*dist_inv*0.5;
            
            //Double_t normrad	= (r1 + 0.5*(dist - r1 - r2))*dist_inv;
           	x1_c		=x1 + normrad*dif_x;
       	    y1_c		=y1 + normrad*dif_y;
		
  	        x2_c		=x1_c;
            y2_c		=y1_c;
        
         	return 2;
		}	
	
	return 3;
}	
//------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
pair<Double_t,Double_t> Ali_TRD_ST_Analyze::fpathLength(Double_t r,Ali_Helix* helixA) const
{
	//taken from https://www.star.bnl.gov/webdata/dox/html/StHelix_8cc_source.html
	
	pair<Double_t,Double_t> value;
	pair<Double_t,Double_t> VALUE(999999999.,999999999.);
	Double_t curvature		= fabs(helixA->getHelix_param(4));
	Double_t radius			= fabs(1/curvature);
	Double_t x0				= helixA->getHelix_param(5) +radius*TMath::Sin(helixA->getHelix_param(2));
	Double_t y0				= helixA->getHelix_param(0) -radius*TMath::Cos(helixA->getHelix_param(2));
	Double_t z0				= helixA->getHelix_param(1);

	Double_t phase			= helixA->getHelix_param(2) -TMath::Pi()/2;
	//Double_t phase=atan2f(-(x0-helixA->getHelix_param(5)),(y0-helixA->getHelix_param(0)));
	
	while(phase<-TMath::Pi()) phase	+= 2.*TMath::Pi() ;
	while(phase>TMath::Pi()) phase	-= 2.*TMath::Pi() ;
	
	Double_t dipangle		= TMath::ATan(helixA->getHelix_param(3));
	Double_t cosdipangle	= TMath::Cos(dipangle);
	Double_t sindipangle	= TMath::Sin(dipangle);
	Double_t h				= TMath::Sign(1,helixA->getHelix_param(4));
	if(phase> TMath::Pi()) phase	-= 2.*TMath::Pi() ;
	Double_t cosphase		= TMath::Cos(phase);
	Double_t sinphase		= TMath::Sin(phase);

        //printf("h: %4.3f, phase: %4.3f, dipangle: %4.3f \n",h,phase,dipangle);

	Double_t t1 	= y0*curvature;
	Double_t t2 	= sinphase;
  	Double_t t3 	= curvature*curvature;
  	Double_t t4 	= y0*t2;
	Double_t t5 	= cosphase;
	Double_t t6 	= x0*t5;
	Double_t t8 	= x0*x0;
	Double_t t11	= y0*y0;
	Double_t t14 	= r*r;
	Double_t t15 	= t14*curvature;
	Double_t t17 	= t8*t8;
	Double_t t19 	= t11*t11;
	Double_t t21 	= t11*t3;
	Double_t t23 	= t5*t5;
	Double_t t32 	= t14*t14;
	Double_t t35 	= t14*t3;
	Double_t t38 	= 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*curvature*t6 +
				4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
				4.0*t8*x0*curvature*t5 - 4.0*t11*t23 -
				4.0*t11*y0*curvature*t2 + 4.0*t11 - 4.0*t14 +
				t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
	//cout<<t38<<endl;
	Double_t t40 	= (-t3*t38);
	if (t40<0.) return VALUE;
        t40 = ::sqrt(t40);

        if (h<0) phase += TMath::Pi();
	
	if(phase>TMath::Pi()) phase	-= 2.*TMath::Pi();

	Double_t t43 	= x0*curvature;
	Double_t t45 	= 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
	Double_t t46 	= h*cosdipangle*curvature;

	value.first 	= (-phase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
	value.second 	= -(phase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

	//
	//   Solution can be off by +/- one period, select smallest
	//
	Double_t p 		= fabs(2*TMath::Pi()/(h*curvature*cosdipangle));
	if (! std::isnan(value.first)) {
		if (fabs(value.first-p) < fabs(value.first)) value.first = value.first-p;
	   	else if (fabs(value.first+p) < fabs(value.first)) value.first = value.first+p;
	}
	if (! std::isnan(value.second)) {
	   	if (fabs(value.second-p) < fabs(value.second)) value.second = value.second-p;
	   	else if (fabs(value.second+p) < fabs(value.second)) value.second = value.second+p;
	}
	
	if (value.first > value.second)
	swap(value.first,value.second);
	return(value);

	
}
//----------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
TVector3 Ali_TRD_ST_Analyze::fDCA_Helix_Estimate(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{

    // Calculates the 2D crossing point, calculates the corresponding 3D point and returns pathA and pathB
    // TV3_photon is a result of a second calculation which assumes its a photon conversion in TRD

    TVector3 TV3_photon;
    TV3_photon.SetXYZ(-999.0,-999.0,-999.0);
    Double_t helix_point[3];

    Double_t x1 = helixA->getHelix_param(5);
    Double_t y1 = helixA->getHelix_param(0);
    Double_t x2 = helixB->getHelix_param(5);
    Double_t y2 = helixB->getHelix_param(0);
    Double_t c1 = helixA->getHelix_param(4);
    Double_t c2 = helixB->getHelix_param(4);
    Double_t r1 = 0.0;
    Double_t r2 = 0.0;
    if(c1 != 0 && c2 != 0)
    {
        r1 = fabs(1.0/c1);
        r2 = fabs(1.0/c2);
    } else return TV3_photon;

    Double_t x1_c = 0.0;
    Double_t y1_c = 0.0;
    Double_t x2_c = 0.0;
    Double_t y2_c = 0.0;

    //Int_t bCross_points = fCross_points_Circles(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);
    //printf("2D circle cross points (Alex): {%4.3f, %4.3f}, {%4.3f, %4.3f} \n",x1_c,y1_c,x2_c,y2_c);

    pathA = 0.0;
    pathB = 0.0;
    dcaAB = 0.0;


    Int_t cCross_points = fCircle_Interception(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);
    //printf("2D circle cross points, circle: {x1: %4.3f, y1: %4.3f, r1: %4.3f} {x2: %4.3f, y2: %4.3f, r2: %4.3f} (Sven): {%4.3f, %4.3f}, {%4.3f, %4.3f}, return: %d \n",x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c,cCross_points);

    Double_t radiusA = sqrt(x1_c*x1_c+y1_c*y1_c);
    Double_t radiusB = sqrt(x2_c*x2_c+y2_c*y2_c);
    //printf("radiusA: %4.3f, radiusB: %4.3f \n",radiusA,radiusB);

    //cout << "bCross_points = " << bCross_points << ", xyr(1) = {" << x1 << ", " << y1 << ", " << r1
    //    << "}, xyr(2) = {"  << x2 << ", " << y2 << ", " << r2 << "}, p1 = {" << x1_c << ", " << y1_c << "}, p2 = {" << x2_c << ", " << y2_c << "}" << endl;

    //if(bCross_points == 0) return 0;

    TVector3 pointA,pointB,pointA1,pointB1,pointA2,pointB2;

    Double_t path_lengthA_c1,path_lengthA_c2,path_lengthB_c1,path_lengthB_c2;

    // first crossing point for helix A
    pair< double, double > path_lengthA = fpathLength(radiusA,helixA);
    Double_t path_lengthA1 = path_lengthA.first;
    Double_t path_lengthA2 = path_lengthA.second;

    helixA ->Evaluate(path_lengthA1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x1_c-pointA1.x())*(x1_c-pointA1.x()) + (y1_c-pointA1.y())*(y1_c-pointA1.y())) <
       ((x1_c-pointA2.x())*(x1_c-pointA2.x()) + (y1_c-pointA2.y())*(y1_c-pointA2.y())))
    {
        path_lengthA_c1 = path_lengthA1;
    }
    else
    {
        path_lengthA_c1 = path_lengthA2;
    }

    // second crossing point for helix A
    path_lengthA = fpathLength(radiusB,helixA);
    path_lengthA1 = path_lengthA.first;
    path_lengthA2 = path_lengthA.second;

    helixA ->Evaluate(path_lengthA1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x2_c-pointA1.x())*(x2_c-pointA1.x()) + (y2_c-pointA1.y())*(y2_c-pointA1.y())) <
       ((x2_c-pointA2.x())*(x2_c-pointA2.x()) + (y2_c-pointA2.y())*(y2_c-pointA2.y())))
    {
        path_lengthA_c2 = path_lengthA1;
    }
    else
    {
        path_lengthA_c2 = path_lengthA2;
    }

    // first crossing point for helix B
    pair< double, double > path_lengthB = fpathLength(radiusA,helixB);
    Double_t path_lengthB1 = path_lengthB.first;
    Double_t path_lengthB2 = path_lengthB.second;

    helixB ->Evaluate(path_lengthB1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x1_c-pointB1.x())*(x1_c-pointB1.x()) + (y1_c-pointB1.y())*(y1_c-pointB1.y())) <
       ((x1_c-pointB2.x())*(x1_c-pointB2.x()) + (y1_c-pointB2.y())*(y1_c-pointB2.y())))
    {
        path_lengthB_c1 = path_lengthB1;
    }
    else
    {
        path_lengthB_c1 = path_lengthB2;
    }

    // second crossing point for helix B
    path_lengthB = fpathLength(radiusB,helixB);
    path_lengthB1 = path_lengthB.first;
    path_lengthB2 = path_lengthB.second;

    helixB ->Evaluate(path_lengthB1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x2_c-pointB1.x())*(x2_c-pointB1.x()) + (y2_c-pointB1.y())*(y2_c-pointB1.y())) <
       ((x2_c-pointB2.x())*(x2_c-pointB2.x()) + (y2_c-pointB2.y())*(y2_c-pointB2.y())))
    {
        path_lengthB_c2 = path_lengthB1;
    }
    else
    {
        path_lengthB_c2 = path_lengthB2;
    }

    helixA ->Evaluate(path_lengthA_c1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA_c2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    helixB ->Evaluate(path_lengthB_c1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB_c2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    TVector3 TV3_pointC1;
    TVector3 TV3_pointC2;

    if((pointA1 - pointB1).Mag() < (pointA2 - pointB1).Mag())
    {
        TV3_pointC1 = pointA1 + pointB1;
        TV3_pointC2 = pointA2 + pointB2;
    }
    else
    {
        TV3_pointC1 = pointA2 + pointB1;
        TV3_pointC2 = pointA1 + pointB2;
    }

    TV3_pointC1 *= 0.5;
    TV3_pointC2 *= 0.5;


    TV3_photon = pointA1 + pointA2 + pointB1 + pointB2;
    TV3_photon *= 0.25;

    TVector3 TV3_dist_AB = TV3_pointC1 - TV3_pointC2;
    if(TV3_dist_AB.Mag() > 50.0)
    {
        if(TV3_pointC1.Perp() < TV3_pointC2.Perp())
        {
            TV3_photon = TV3_pointC1;
        }
        else
        {
            TV3_photon = TV3_pointC2;
        }
    }


    Double_t dcaAB1 = (pointA1-pointB1).Mag();
    Double_t dcaAB2 = (pointA2-pointB2).Mag();
    Double_t dcaAB3 = (pointA1-pointB2).Mag();
    Double_t dcaAB4 = (pointA2-pointB1).Mag();

#if 0
    printf("pointA1: {%4.3f, %4.3f, %4.3f} \n",pointA1.X(),pointA1.Y(),pointA1.Z());
    printf("pointA2: {%4.3f, %4.3f, %4.3f} \n",pointA2.X(),pointA2.Y(),pointA2.Z());
    printf("pointB1: {%4.3f, %4.3f, %4.3f} \n",pointB1.X(),pointB1.Y(),pointB1.Z());
    printf("pointB2: {%4.3f, %4.3f, %4.3f} \n",pointB2.X(),pointB2.Y(),pointB2.Z());
    printf("dcaAB1: %4.3f, dcaAB2: %4.3f, dcaAB3: %4.3f, dcaAB4: %4.3f \n",dcaAB1,dcaAB2,dcaAB3,dcaAB4);
#endif

    if(dcaAB1 < dcaAB2)
    {
        pathA = path_lengthA_c1;
        pathB = path_lengthB_c1;
        dcaAB = dcaAB1;
    }
    else
    {
        pathA = path_lengthA_c2;
        pathB = path_lengthB_c2;
        dcaAB = dcaAB2;
    }


    return TV3_photon;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// V1.0
void Ali_TRD_ST_Analyze::fHelixABdca(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB,
                                    Float_t pathA_in, Float_t pathB_in)
{
    //cout << "Standard fHelixABdca called..." << endl;
    Float_t pA[2];
    Float_t pB[2]; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    if(pathA_in < -9990.0 && pathB_in < -9990.0)
    {
        pA[0] = 0.0;
        pA[1] = 0.0;
        pB[0] = 0.0;
        pB[1] = -70.0;
    }
    else
    {
        pA[0] = pathA_in+5.0;
        pA[1] = pathA_in-5.0;
        pB[0] = pathB_in+5.0;
        pB[1] = pathB_in-5.0;
    }
    Float_t distarray[2];
    TVector3 testA, testB;
    Double_t helix_point[3];
    for(Int_t r = 0; r < 2; r++)
    {
        helixB ->Evaluate(pB[r],helix_point);  // 3D-vector of helixB point at path pB[r]
        testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation

        helixA ->Evaluate(pathA_dca,helix_point);  // 3D-vector of helixB point at path pB[r]
        testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        distarray[r] = (testA-testB).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            helixB ->Evaluate(pB[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            helixA ->Evaluate(pA[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[0] = (testA-testB).Mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            helixB ->Evaluate(pB[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            helixA ->Evaluate(pA[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[1] = (testA-testB).Mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Float_t Ali_TRD_ST_Analyze::Calc_nuclev_bitmap(vector<Int_t> vec_idx_kalman_tracks_nuclev_in)
{
    // Calculated how many independent tracklets are used for the nuclear interaction vertex candidate.
    // Returns a bitmap (6 layers, each layer has 5 bits to store the number of tracklets)

#if 0
    for(Int_t idx_i_kalman_track = 0; idx_i_kalman_track < (Int_t)vec_idx_kalman_tracks_nuclev_in.size() - 1; idx_i_kalman_track++)
    {
        for(Int_t idx_i_kalman_trackB = (idx_i_kalman_track+1); idx_i_kalman_trackB < (Int_t)vec_idx_kalman_tracks_nuclev_in.size(); idx_i_kalman_trackB++)
        {
            Int_t i_kalman_track_B = vec_idx_kalman_tracks_nuclev_in[idx_i_kalman_trackB];
            if(i_kalman_track_A == i_kalman_track_B)
            {
                //printf("%s WARNING! Same kalman track used twice in nuclev, %s trackA: %d(%d), trackB: %d(%d) \n",KRED,KNRM,i_kalman_track_A,idx_i_kalman_track,i_kalman_track_B,idx_i_kalman_trackB);
            }
        }
    }
#endif

    vector< vector<Int_t> > vec_nuclev_kalman_tracklets_per_layer;
    vec_nuclev_kalman_tracklets_per_layer.resize(6);
    for(Int_t i_lay = 0; i_lay < 6; i_lay++)
    {
        for(Int_t idx_i_kalman_track = 0; idx_i_kalman_track < (Int_t)vec_idx_kalman_tracks_nuclev_in.size(); idx_i_kalman_track++)
        {
            Int_t i_kalman_track_A = vec_idx_kalman_tracks_nuclev_in[idx_i_kalman_track];
            if(vec_kalman_TRD_trackets[i_kalman_track_A][i_lay] != NULL)
            {
                Int_t idx_TRD_tracklet = vec_kalman_TRD_trackets[i_kalman_track_A][i_lay]->get_TRD_index();
                Int_t flag_used = 0;
                for(Int_t idx = 0; idx < (Int_t)vec_nuclev_kalman_tracklets_per_layer[i_lay].size(); idx++)
                {
                    Int_t idx_TRD_tracklet_stored = vec_nuclev_kalman_tracklets_per_layer[i_lay][idx];
                    if(idx_TRD_tracklet_stored == idx_TRD_tracklet)
                    {
                        flag_used = 1;
                        break;
                    }
                }
                if(!flag_used) vec_nuclev_kalman_tracklets_per_layer[i_lay].push_back(idx_TRD_tracklet);
            }
        }
    }

    Int_t bitmap_return = 0;
    for(Int_t i_lay = 0; i_lay < 6; i_lay++)
    {
        Int_t N_independent_tracklets = (Int_t)vec_nuclev_kalman_tracklets_per_layer[i_lay].size();

        for(Int_t i_bit = 0; i_bit < 5; i_bit++) // asume only 5 bits are filled -> 32 tracklets per layer at maximum
        {
            if(((Int_t)N_independent_tracklets >> i_bit) & 1)
            {
                bitmap_return |= 1 << (i_bit + 5*i_lay); // set bit i_layer to 1
            }
        }

#if 0
        printf("%s i_lay: %d %s, N_independent_tracklets: %d \n",KMAG,i_lay,KNRM,N_independent_tracklets);
        for(Int_t i_trk = 0; i_trk < N_independent_tracklets; i_trk++)
        {
            printf(" ----> i_trk: %d, idx: %d \n",i_trk,vec_nuclev_kalman_tracklets_per_layer[i_lay][i_trk]);
        }
#endif
    }

    return (Float_t)bitmap_return;
}
//----------------------------------------------------------------------------------------

void Ali_TRD_ST_Analyze::set_self_event_info()
{
// ???? will it work here  ? Fill event information
    TRD_Self_Event ->clearPhotonList();
    TRD_Self_Event ->clearNucInteractionsList();

    TRD_Self_Event ->setTriggerWord( TRD_ST_Event->getTriggerWord() );
    TRD_Self_Event ->setx( TRD_ST_Event->getx() );
    TRD_Self_Event ->sety( TRD_ST_Event->gety() );
    TRD_Self_Event ->setz( TRD_ST_Event->getz() );
    TRD_Self_Event ->setid( TRD_ST_Event->getid() );
    TRD_Self_Event ->setN_TPC_tracks( TRD_ST_Event->getN_tracks() );
    TRD_Self_Event ->setN_TRD_tracklets( TRD_ST_Event->getN_TRD_tracklets() );
    TRD_Self_Event ->setBeamIntAA( TRD_ST_Event->getBeamIntAA() );
    TRD_Self_Event ->setT0zVertex( TRD_ST_Event->getT0zVertex() );
    TRD_Self_Event ->setcent_class_ZNA( TRD_ST_Event->getcent_class_ZNA() );
    TRD_Self_Event ->setcent_class_ZNC( TRD_ST_Event->getcent_class_ZNC() );
    TRD_Self_Event ->setcent_class_V0A( TRD_ST_Event->getcent_class_V0A() );
    TRD_Self_Event ->setcent_class_V0C( TRD_ST_Event->getcent_class_V0C() );
    TRD_Self_Event ->setcent_class_V0M( TRD_ST_Event->getcent_class_V0M() );
    TRD_Self_Event ->setcent_class_CL0( TRD_ST_Event->getcent_class_CL0() );
    TRD_Self_Event ->setcent_class_CL1( TRD_ST_Event->getcent_class_CL1() );
    TRD_Self_Event ->setcent_class_SPD( TRD_ST_Event->getcent_class_SPD() );
    TRD_Self_Event ->setcent_class_V0MEq( TRD_ST_Event->getcent_class_V0MEq() );
    TRD_Self_Event ->setcent_class_V0AEq( TRD_ST_Event->getcent_class_V0AEq() );
    TRD_Self_Event ->setcent_class_V0CEq( TRD_ST_Event->getcent_class_V0CEq() );
    //printf("Event information filled \n");
}


//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::Calculate_secondary_vertices(Int_t graphics, Int_t flag_TRD_TPC_tracks)
{
    

    Double_t min_radius_cut = 250.0;
    Double_t max_radius_cut = 356.0;
    if(flag_TRD_TPC_tracks == 0) // TRD tracks
    {
        //printf("Secondary vertex finding: TRD tracks used \n");
        vec_helices = vec_helices_TRD;
    }
    else // TPC tracks
    {
        //printf("Secondary vertex finding: TPC tracks used \n");
        vec_helices = vec_helices_TPC;
        min_radius_cut = 10.0;
        max_radius_cut = 85.0;
    }

    // XALEX
    //min_radius_cut = 0.0;
    //max_radius_cut = 356.0;

    Int_t flag_found_good_AP_vertex = 0;

    Float_t Arr_seconary_params[26];

    Double_t helix_pointA[3];
    Double_t helix_pointB[3];
    Double_t helix_pointAs[3];
    Double_t helix_pointBs[3];
    Double_t vertex_point[3];
    Double_t helix_pointA_est[3];
    Double_t helix_pointB_est[3];
    Double_t vertex_point_est[3];
    TVector3 TV3_dirA;
    TVector3 TV3_dirB;
    TVector3 TV3_dirAB;
    TLorentzVector TLV_A;
    TLorentzVector TLV_B;
    TLorentzVector TLV_AB;
    Double_t EnergyA, EnergyB, Energy_AB;
    Double_t Inv_mass_AB, Inv_mass_AB_K0s, Inv_mass_AB_Lambda, Inv_mass_AB_antiLambda, Momentum_AB, pT_AB, Eta_AB, Phi_AB;

    TVector3 TV3_prim_vertex(EventVertexX,EventVertexY,EventVertexZ);
    TVector3 TV3_sec_vertex;
    TVector3 TV3_dir_prim_sec_vertex;

    vector< vector<Int_t> > vec_idx_kalman_sec_vert;
    vec_idx_kalman_sec_vert.clear();
    vec_idx_kalman_sec_vert.resize(2);
    vec_TV3_secondary_vertices.clear();

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_photon_vertices     = new TEvePointSet();
        TEveP_sec_vertices        = new TEvePointSet();
        TEveP_close_TPC_photon    = new TEvePointSet();
        TEveP_nucl_int_vertices   = new TEvePointSet();
        TEveP_first_point_helix   = new TEvePointSet();
        TEveP_second_point_helix  = new TEvePointSet();
    }
#endif

    Int_t i_vertex           = 0;
    Int_t i_comb             = 0;
    Int_t i_vertex_acc       = 0;
    Int_t i_close_TPC_photon = 0;

    //printf("N_helices: %d \n",(Int_t)vec_helices.size());
    for(Int_t i_track_A = 0; i_track_A < ((Int_t)vec_helices.size() - 1); i_track_A++) // TRD Kalman track A
    {        
        for(Int_t i_track_B = (i_track_A+1); i_track_B < (Int_t)vec_helices.size(); i_track_B++) // TRD Kalman track B
        {            
            Float_t pathA_est, pathB_est, dcaAB_est;
            //printf("i_track_A: %d, i_track_B: %d, i_comb: %d \n",i_track_A,i_track_B,i_comb);
            TVector3 TV3_photon = fDCA_Helix_Estimate(vec_helices[i_track_A],vec_helices[i_track_B],pathA_est,pathB_est,dcaAB_est);

            //printf("track A,B: {%d, %d}, dcaAB_est: %4.3f \n",i_track_A,i_track_B,dcaAB_est);
            //printf(" \n");

            if(dcaAB_est > 15.0) continue;

            vec_helices[i_track_A] ->Evaluate(pathA_est,helix_pointA_est);
            vec_helices[i_track_B] ->Evaluate(pathB_est,helix_pointB_est);
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                vertex_point_est[i_xyz] = (helix_pointA_est[i_xyz] + helix_pointB_est[i_xyz])/2.0;
            }
            Double_t radius_est = TMath::Sqrt(vertex_point_est[0]*vertex_point_est[0] + vertex_point_est[1]*vertex_point_est[1]);
            if(radius_est > 400.0) continue;


            Float_t pathA, pathB, dcaAB;
            if(TV3_photon.Mag() > 1000)
            {
                pathA_est = -9999.0;
                pathB_est = -9999.0;
            }
            fHelixABdca(vec_helices[i_track_A],vec_helices[i_track_B],pathA,pathB,dcaAB,pathA_est,pathB_est);
            //printf("  --> track A,B: {%d, %d}, dcaAB: %4.3f \n",i_track_A,i_track_B,dcaAB);

            if(dcaAB < 5.0)
            {
                //------------------------------------------------------------
                // This is for photons, replaces the fHelixABdca calculation. It takes the average of all intersections.
                Float_t pathA_photon, pathB_photon, dcaA_photon, dcaB_photon;
                fHelixAtoPointdca(TV3_photon,vec_helices[i_track_A],pathA_photon,dcaA_photon); // new helix to point dca calculation
                fHelixAtoPointdca(TV3_photon,vec_helices[i_track_B],pathB_photon,dcaB_photon); // new helix to point dca calculation
                pathA = pathA_photon;
                pathB = pathB_photon;
                dcaAB = (dcaA_photon + dcaB_photon)/2.0;
                vertex_point[0] = TV3_photon.X();
                vertex_point[1] = TV3_photon.Y();
                vertex_point[2] = TV3_photon.Z();
                //------------------------------------------------------------



                //------------------------------------------------------------
                // Calculate number of shared and independent tracklets
                Int_t N_independent_AB[2] = {0};
                Int_t N_shared_AB         = 0;

                // number |= 1 << bit_to_set; // set bit bit_to_set to 1
                // Int_t bit_status = (TRD_ADC_time_layer[i_layer] >> bitcheck) & 1; // check bit bitcheck
                // Int_t     HasITShit_on_layer(Int_t ilayer) { return ((NITScls >> ilayer) & 1);}  // ITShit -> LOL
                Int_t bit_TRD_layer_shared = 0;
                if(flag_TRD_TPC_tracks != 0) bit_TRD_layer_shared = -1; // TPC track



                //------------------------
                if(flag_TRD_TPC_tracks == 0) // TRD track
                {
                    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
                    {
                        TVector3 kalman_TV3_offset_A;
                        Int_t flag_use_A = 0;
                        if(vec_kalman_TRD_trackets[i_track_A][i_layer] != NULL)
                        {
                            kalman_TV3_offset_A = vec_kalman_TRD_trackets[i_track_A][i_layer] ->get_TV3_offset();
                            flag_use_A = 1;
                        }

                        TVector3 kalman_TV3_offset_B;
                        Int_t flag_use_B = 0;
                        if(vec_kalman_TRD_trackets[i_track_B][i_layer] != NULL)
                        {
                            kalman_TV3_offset_B = vec_kalman_TRD_trackets[i_track_B][i_layer] ->get_TV3_offset();
                            flag_use_B = 1;
                        }

                        Int_t flag_shared = 0;
                        if(flag_use_A && flag_use_B)
                        {
                            TVector3 TV3_diff_AB = kalman_TV3_offset_A - kalman_TV3_offset_B;
                            Double_t diff_AB = TV3_diff_AB.Mag();

                            if(diff_AB < 0.1)
                            {
                                N_shared_AB++;
                                flag_shared = 1;
                                bit_TRD_layer_shared |= 1 << (i_layer + 12); // set bit i_layer to 1
                                //printf("Shared layer: %d \n",i_layer);
                            }
                        }
                        if(!flag_shared) // only one tracklet available per track or tracklets are not shared
                        {
                            if(flag_use_A)
                            {
                                N_independent_AB[0]++;
                                bit_TRD_layer_shared |= 1 << i_layer; // set bit i_layer to 1
                                //printf("Independent layer track A: %d \n",i_layer);
                            }
                            if(flag_use_B)
                            {
                                N_independent_AB[1]++;
                                bit_TRD_layer_shared |= 1 << (i_layer + 6); // set bit i_layer to 1
                                //printf("Independent layer track B: %d \n",i_layer);
                            }
                        }
                    }

#if 0
                    // Test the bit map
                    Float_t test = (Float_t)bit_TRD_layer_shared;
                    for(Int_t i_bit = 0; i_bit < 18; i_bit++)
                    {
                        cout << "i_bit: " << i_bit << ", value: " << (((Int_t)test >> i_bit) & 1) << endl;
                    }
#endif
                    //printf("%s Number of shared tracklets: %s %d, independent A,B: {%d, %d} \n",KGRN,KNRM,N_shared_AB, N_independent_AB[0], N_independent_AB[1]);
                    //------------------------------------------------------------
                }
                //------------------------



                //if(i_comb == 38)
                {
                    vec_helices[i_track_A] ->Evaluate(pathA,helix_pointA);
                    vec_helices[i_track_B] ->Evaluate(pathB,helix_pointB);
                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        vertex_point[i_xyz] = (helix_pointA[i_xyz] + helix_pointB[i_xyz])/2.0;
                    }

                    TV3_sec_vertex.SetXYZ(vertex_point[0],vertex_point[1],vertex_point[2]);
                    vec_TV3_secondary_vertices.push_back(TV3_sec_vertex);
                    vec_idx_kalman_sec_vert[0].push_back(i_track_A);
                    vec_idx_kalman_sec_vert[1].push_back(i_track_B);

                    TV3_dir_prim_sec_vertex = TV3_sec_vertex - TV3_prim_vertex;
                    if(TV3_dir_prim_sec_vertex.Mag() <= 0.0) continue;
                    TV3_dir_prim_sec_vertex *= 1.0/TV3_dir_prim_sec_vertex.Mag();

                    Double_t radius_vertex = TMath::Sqrt(vertex_point[0]*vertex_point[0] + vertex_point[1]*vertex_point[1]);
                    Int_t i_radius = (Int_t)(radius_vertex/Delta_AP_radius);


                    Double_t pTA, pzA, pA, pTB, pzB, pB, CA, CB;
                    Double_t Chi_2A = -1.0;
                    Double_t Chi_2B = -1.0;
                    Double_t TPCdEdx_A = -1.0;
                    Double_t dca_TPC_A = -1.0;
                    Double_t p_TPC_A   = -1.0;
                    Double_t TPCdEdx_B = -1.0;
                    Double_t dca_TPC_B = -1.0;
                    Double_t p_TPC_B   = -1.0;

                    if(flag_TRD_TPC_tracks == 0) // TRD track
                    {
                        pTA = mHelices_kalman[i_track_A][6]; // pT
                        pzA = mHelices_kalman[i_track_A][7]; // pz
                        pA  = TMath::Sqrt(pTA*pTA + pzA*pzA); // p
                        pTB = mHelices_kalman[i_track_B][6]; // pT
                        pzB = mHelices_kalman[i_track_B][7]; // pz
                        pB  = TMath::Sqrt(pTB*pTB + pzB*pzB); // p
                        CA  = mHelices_kalman[i_track_A][4]; // curvature
                        CB  = mHelices_kalman[i_track_B][4]; // curvature

                        Chi_2A = mChi_2s_kalman[i_track_A];
                        Chi_2B = mChi_2s_kalman[i_track_B];

                        //n_(Int_t)vec_kalman_TRD_trackets[i_track_A].size();

                    }
                    else // TPC track
                    {
                        pTA = mHelices_TPC[i_track_A][6]; // pT
                        pzA = mHelices_TPC[i_track_A][7]; // pz
                        pA  = TMath::Sqrt(pTA*pTA + pzA*pzA); // p
                        pTB = mHelices_TPC[i_track_B][6]; // pT
                        pzB = mHelices_TPC[i_track_B][7]; // pz
                        pB  = TMath::Sqrt(pTB*pTB + pzB*pzB); // p
                        CA  = mHelices_TPC[i_track_A][4]; // curvature
                        CB  = mHelices_TPC[i_track_B][4]; // curvature


                        TPCdEdx_A = PID_params_TPC[i_track_A][0];
                        dca_TPC_A = PID_params_TPC[i_track_A][1];
                        p_TPC_A   = PID_params_TPC[i_track_A][2];

                        TPCdEdx_B = PID_params_TPC[i_track_B][0];
                        dca_TPC_B = PID_params_TPC[i_track_B][1];
                        p_TPC_B   = PID_params_TPC[i_track_B][2];

                        //printf("CA: {%4.5f, %4.5f}, p: {%4.5f, %4.5f} \n",CA,vec_helices[i_track_A]->getHelix_param(4),pA,p_TPC_A);
                        //printf("CA: %4.3f, CB: %4.5f, pTA: %4.5f, pTB: %4.3f \n",CA,CB,pTA,pTB);
                    }

                    //if(radius_vertex > 240.0 && radius_vertex < 360.0)
#if defined(USEEVE)
                    //if(graphics) TEveP_sec_vertices ->SetPoint(i_vertex,vertex_point[0],vertex_point[1],vertex_point[2]);
#endif
                    //-------------------------------------------------
                    // Armenteros-Podolanski
                    //printf("radius_vertex: %4.3f \n",radius_vertex);
                    if(radius_vertex > min_radius_cut && radius_vertex < max_radius_cut) // TRD acceptance in R-direction
                    {
                        vec_helices[i_track_A] ->Evaluate(pathA+0.1,helix_pointAs);
                        vec_helices[i_track_B] ->Evaluate(pathB+0.1,helix_pointBs);
                        TV3_dirA.SetXYZ(helix_pointAs[0] - helix_pointA[0],helix_pointAs[1] - helix_pointA[1],helix_pointAs[2] - helix_pointA[2]);
                        TV3_dirB.SetXYZ(helix_pointBs[0] - helix_pointB[0],helix_pointBs[1] - helix_pointB[1],helix_pointBs[2] - helix_pointB[2]);

                        TV3_dirA *= pA/TV3_dirA.Mag();
                        TV3_dirB *= pB/TV3_dirB.Mag();
                        TV3_dirAB = TV3_dirA + TV3_dirB;
                        TV3_dirAB *= 1.0/TV3_dirAB.Mag();
                        Double_t projA = TV3_dirA.Dot(TV3_dirAB);
                        Double_t projB = TV3_dirB.Dot(TV3_dirAB);
                        Double_t AP_pT = TMath::Sqrt(TMath::Power(TV3_dirA.Mag(),2) - projA*projA);
                        Double_t AP_pTB = TMath::Sqrt(TMath::Power(TV3_dirB.Mag(),2) - projB*projB);
                        Double_t AP_alpha = (projA - projB)/(projA + projB);
                        if(CA < CB) AP_alpha *= -1.0;

                        TLV_A.SetXYZM(TV3_dirA.X(),TV3_dirA.Y(),TV3_dirA.Z(),0.13957);
                        TLV_B.SetXYZM(TV3_dirB.X(),TV3_dirB.Y(),TV3_dirB.Z(),0.13957);
                        TLV_AB = TLV_A + TLV_B;
                        Inv_mass_AB_K0s = TLV_AB.M();

                        TLV_A.SetXYZM(TV3_dirA.X(),TV3_dirA.Y(),TV3_dirA.Z(),0.938272);
                        TLV_B.SetXYZM(TV3_dirB.X(),TV3_dirB.Y(),TV3_dirB.Z(),0.13957);
                        TLV_AB = TLV_A + TLV_B;
                        Inv_mass_AB_Lambda = TLV_AB.M();

                        TLV_A.SetXYZM(TV3_dirA.X(),TV3_dirA.Y(),TV3_dirA.Z(),0.13957);
                        TLV_B.SetXYZM(TV3_dirB.X(),TV3_dirB.Y(),TV3_dirB.Z(),0.938272);
                        TLV_AB = TLV_A + TLV_B;
                        Inv_mass_AB_antiLambda = TLV_AB.M();


                        TLV_A.SetXYZM(TV3_dirA.X(),TV3_dirA.Y(),TV3_dirA.Z(),0.000511);
                        TLV_B.SetXYZM(TV3_dirB.X(),TV3_dirB.Y(),TV3_dirB.Z(),0.000511);
                        TLV_AB = TLV_A + TLV_B;
                        Inv_mass_AB = TLV_AB.M();
                        Eta_AB      = TLV_AB.Eta();
                        Phi_AB      = TLV_AB.Phi();
                        Energy_AB   = TLV_AB.Energy();
                        Momentum_AB = TLV_AB.P();
                        pT_AB 	    = TLV_AB.Pt();

                        Double_t dot_product_dir_vertex = TV3_dirAB.Dot(TV3_dir_prim_sec_vertex);


                        TVector3 TV3_close_TPC_photon;
                        Int_t flag_close_TPC_photon = 0;

                        // AP cuts from https://www.physi.uni-heidelberg.de//Publications/PhDThesis_Leardini.pdf (eqn. 5.2)
                        Double_t AP_alpha_max = 0.95;
                        Double_t AP_qT_max    = 0.05;
                        Double_t AP_cut_value = 1.0;
                        Double_t AP_value = TMath::Power(AP_alpha/AP_alpha_max,2.0) + TMath::Power(AP_pT/AP_qT_max,2.0);

                        //if(fabs(AP_alpha) < 0.2 && AP_pT > 0.0 && AP_pT < 0.02 && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.5 && pTB < 0.5 && dot_product_dir_vertex > 0.9) // TRD photon conversion
                        //printf("AP_value: %4.3f, dot_product_dir_vertex: %4.3f, CA*CB: %4.3f, pTA: %4.3f, pTB: %4.3f \n",AP_value,dot_product_dir_vertex,CA*CB,pTA,pTB);

                        // XALEX
                        if(AP_value < AP_cut_value && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.8 && pTB < 0.8 && dot_product_dir_vertex > 0.9)
                        {                            
                            Double_t dca_min  = 999.0;
                            Double_t path_min = -999.0;
                            Int_t    i_track_min = -1;

                            //------------
                            if(flag_TRD_TPC_tracks == 0) // TRD track
                            {
                                UShort_t NumTracks = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
                                for(Int_t i_track = 0; i_track < NumTracks; i_track++)
                                {
                                    Float_t pathA_dca = -1.0;
                                    Float_t dcaAB_dca = -1.0;
                                    TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);
                                    TPC_single_helix ->setHelix(TRD_ST_TPC_Track->getHelix_param(0),TRD_ST_TPC_Track->getHelix_param(1),TRD_ST_TPC_Track->getHelix_param(2),TRD_ST_TPC_Track->getHelix_param(3),TRD_ST_TPC_Track->getHelix_param(4),TRD_ST_TPC_Track->getHelix_param(5));
                                    fHelixAtoPointdca(TV3_sec_vertex,TPC_single_helix,pathA_dca,dcaAB_dca); // new helix to point dca calculation

                                    Double_t helix_point_TPC_photon[3];
                                    TPC_single_helix ->Evaluate(pathA_dca,helix_point_TPC_photon);  // 3D-vector of helixB point at path pB[r]
                                    TV3_close_TPC_photon.SetXYZ(helix_point_TPC_photon[0],helix_point_TPC_photon[1],helix_point_TPC_photon[2]);
                                    flag_close_TPC_photon = 1;

                                    if(dcaAB_dca < dca_min)
                                    {
                                        dca_min     = dcaAB_dca;
                                        path_min    = pathA_dca;
                                        i_track_min = i_track;
                                    }
                                }
                            }
                            //------------

                            // XALEX
                            //flag_close_TPC_photon = 0;

                            //-----------------------------------
                            // Topology photon conversion cuts - for print out
                            Int_t independent_layer_A[6] = {0};
                            Int_t independent_layer_B[6] = {0};
                            Int_t shared_layer[6]        = {0};
                            for(Int_t i_bit = 0; i_bit < 18; i_bit++)
                            {
                                if(((Int_t)bit_TRD_layer_shared >> i_bit) & 1)
                                {
                                    if(i_bit < 6)                independent_layer_A[i_bit]   = 1;
                                    if(i_bit >= 6 && i_bit < 12) independent_layer_B[i_bit-6] = 1;
                                    if(i_bit >= 12)              shared_layer[i_bit-12]       = 1;
                                }
                            }


                            if(flag_TRD_TPC_tracks == 0) // TRD track
                            {
                                if(
                                   (
                                    (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 1 &&
                                    (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3]) > 1 &&
                                    (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3]) > 1
                                   )
                                   ||
                                   (
                                    (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 2 &&
                                    (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3]) > 0 &&
                                    (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3]) > 0
                                   )
                                   ||
                                   (
                                    (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 0 &&
                                    (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3] + independent_layer_A[2]) > 2 &&
                                    (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3] + independent_layer_B[2]) > 2
                                   )
                                  )
                                {
                                    if(dca_min > 10.0 || (dca_min <= 10.0 && path_min < 0.0)) // no close-by TPC track
                                    {
#if defined(USEEVE)
                                        if(graphics && flag_close_TPC_photon) TEveP_close_TPC_photon ->SetPoint(i_close_TPC_photon,TV3_close_TPC_photon[0],TV3_close_TPC_photon[1],TV3_close_TPC_photon[2]);
#endif
                                        i_close_TPC_photon++;
                                        printf("%s Number of shared tracklets: %s %d, independent A,B: {%d, %d}: dca_min: %4.3f, i_track_min: %d \n",KGRN,KNRM,N_shared_AB, N_independent_AB[0], N_independent_AB[1],dca_min,i_track_min);
                                        printf("      --> Found vertex for AP in %s event: %lld %s at radius: %4.3f, pos: {%4.3f, %4.3f, %4.3f}, AP_pT: %4.3f, AP_alpha: %4.3f, pTA: %4.3f, pTB: %4.3f, dot: %4.3f, Inv_mass_AB: %4.3f, Energy_AB: %4.3f, Momentum_AB: %4.3f \n",KBLU,Global_Event,KNRM,radius_vertex,TV3_sec_vertex.X(),TV3_sec_vertex.Y(),TV3_sec_vertex.Z(),AP_pT,AP_alpha,pTA,pTB,dot_product_dir_vertex,Inv_mass_AB,Energy_AB,Momentum_AB);
                                    }
                                }
                            }
                            //-----------------------------------

                            //-----------------------------------
                            // Ntuple for photon candidates
                            Arr_seconary_params[0]  = (Float_t)vertex_point[0];
                            Arr_seconary_params[1]  = (Float_t)vertex_point[1];
                            Arr_seconary_params[2]  = (Float_t)vertex_point[2];
                            Arr_seconary_params[3]  = (Float_t)bit_TRD_layer_shared;
                            Arr_seconary_params[4]  = (Float_t)pT_AB;
                            Arr_seconary_params[5]  = (Float_t)pTA*TMath::Sign(1,CA);
                            Arr_seconary_params[6]  = (Float_t)pTB*TMath::Sign(1,CB);
                            Arr_seconary_params[7]  = (Float_t)AP_pT;
                            Arr_seconary_params[8]  = (Float_t)AP_alpha;
                            Arr_seconary_params[9]  = (Float_t)dca_min;
                            Arr_seconary_params[10] = (Float_t)path_min;
                            Arr_seconary_params[11] = (Float_t)Inv_mass_AB;
                            Arr_seconary_params[12] = (Float_t)Eta_AB;
                            Arr_seconary_params[13] = (Float_t)Phi_AB;
                            Arr_seconary_params[14] = (Float_t)Global_Event;
                            Arr_seconary_params[15] = (Float_t)dot_product_dir_vertex;
                            Arr_seconary_params[16] = (Float_t)TPCdEdx_A;
                            Arr_seconary_params[17] = (Float_t)dca_TPC_A;
                            Arr_seconary_params[18] = (Float_t)p_TPC_A;
                            Arr_seconary_params[19] = (Float_t)TPCdEdx_B;
                            Arr_seconary_params[20] = (Float_t)dca_TPC_B;
                            Arr_seconary_params[21] = (Float_t)p_TPC_B;
                            Arr_seconary_params[22] = (Float_t)Inv_mass_AB_K0s;
                            Arr_seconary_params[23] = (Float_t)dcaAB;
                            Arr_seconary_params[24] = (Float_t)Inv_mass_AB_Lambda;
                            Arr_seconary_params[25] = (Float_t)Inv_mass_AB_antiLambda;

                            //printf("vertex pos: {%4.3f, %4.3f, %4.3f}, ntracks: %4.3f \n",Arr_seconary_params[0],Arr_seconary_params[1],Arr_seconary_params[2],Arr_seconary_params[3]);
                            //NT_secondary_vertices ->Fill(Arr_seconary_params);

                            //printf("bit_TRD_layer_shared: %d \n",bit_TRD_layer_shared);
                            //-----------------------------------

                            //-----FILL NEW CLASS PHOTON-----------------------------------

                            TRD_Photon = TRD_Self_Event ->createPhoton();

                            TRD_Photon ->set_vertex_point((Float_t)vertex_point[0],(Float_t)vertex_point[1],(Float_t)vertex_point[2]);
                            TRD_Photon ->set_bit_TRD_layer_shared((Float_t)bit_TRD_layer_shared);
                            TRD_Photon ->set_pT_AB((Float_t)pT_AB);
                            TRD_Photon ->set_AP_pT((Float_t)AP_pT);
                            TRD_Photon ->set_AP_alpha((Float_t)AP_alpha);
                            TRD_Photon ->set_dca_min((Float_t)dca_min);
                            TRD_Photon ->set_path_min((Float_t)path_min);
                            TRD_Photon ->set_Inv_mass_AB((Float_t)Inv_mass_AB);
                            TRD_Photon ->set_Eta_AB((Float_t)Eta_AB);
                            TRD_Photon ->set_Phi_AB((Float_t)Phi_AB);
                            TRD_Photon ->set_dot_product_dir_vertex((Float_t)dot_product_dir_vertex);
                            TRD_Photon ->set_Inv_mass_AB_K0s((Float_t)Inv_mass_AB_K0s);
                            TRD_Photon ->set_dcaAB((Float_t)dcaAB);
                            TRD_Photon ->set_Inv_mass_AB_Lambda((Float_t)Inv_mass_AB_Lambda);
                            TRD_Photon ->set_Inv_mass_AB_antiLambda((Float_t)Inv_mass_AB_antiLambda);
                            TRD_Photon ->set_TLV_part_A(TLV_A);
                            TRD_Photon ->set_TLV_part_B(TLV_B);

                            //--------- fill TRD/TPC track info --------------------------

                            if (flag_TRD_TPC_tracks == 0) //Kalman
                            {
                                for (Int_t i_trackAB = 0; i_trackAB < 2; i_trackAB++)
                                {
                                    //printf("test1 \n");
                                    TRD_Kalman_track = TRD_Photon -> createKalman_Track();
                                    //printf("test2 \n");


                                    if (i_trackAB == 0 ) { 
                                        TRD_Kalman_track ->set_Chi2(Chi_2A); 
                                        TRD_Kalman_track ->setKalmanHelix_param(vec_helices[i_track_A]->getHelix_param(0),vec_helices[i_track_A]->getHelix_param(1),
                                            vec_helices[i_track_A]->getHelix_param(2),vec_helices[i_track_A]->getHelix_param(3),
                                            vec_helices[i_track_A]->getHelix_param(4),vec_helices[i_track_A]->getHelix_param(5)); 

                                    }
                                    if (i_trackAB == 1 ) { 
                                        TRD_Kalman_track ->set_Chi2(Chi_2B);
                                        TRD_Kalman_track ->setKalmanHelix_param(vec_helices[i_track_B]->getHelix_param(0),vec_helices[i_track_B]->getHelix_param(1),
                                            vec_helices[i_track_B]->getHelix_param(2),vec_helices[i_track_B]->getHelix_param(3),
                                            vec_helices[i_track_B]->getHelix_param(4),vec_helices[i_track_B]->getHelix_param(5)); 
                                    }
                                }
                            }

                            if (flag_TRD_TPC_tracks == 1) //TPC
                            {
                                for (Int_t i_trackAB = 0; i_trackAB < 2; i_trackAB++)
                                {
                                    TPC_track = TRD_Photon -> createTPC_Track();

                                    if (i_trackAB == 0 ) { 
                                        TPC_track ->setHelix(vec_helices[i_track_A]->getHelix_param(0),vec_helices[i_track_A]->getHelix_param(1),
                                            vec_helices[i_track_A]->getHelix_param(2),vec_helices[i_track_A]->getHelix_param(3),
                                            vec_helices[i_track_A]->getHelix_param(4),vec_helices[i_track_A]->getHelix_param(5)); 

                                    }
                                    if (i_trackAB == 1 ) { 
                                        TPC_track ->setHelix(vec_helices[i_track_B]->getHelix_param(0),vec_helices[i_track_B]->getHelix_param(1),
                                            vec_helices[i_track_B]->getHelix_param(2),vec_helices[i_track_B]->getHelix_param(3),
                                            vec_helices[i_track_B]->getHelix_param(4),vec_helices[i_track_B]->getHelix_param(5)); 
                                    }
                                }
                            }

                            //--------- end fill Photon and TRD/TPC track info --------------------------



                            //-------------------------------------------------------------


                            if(flag_TRD_TPC_tracks == 0) // TRD tracks
                            {
#if defined(USEEVE)
                                //if(AP_pT > 0.05 && AP_pT < 0.25 && CA*CB < 0.0 && pTA > 0.2 && pTB > 0.2 && dot_product_dir_vertex > 0.95)
                                //if(fabs(AP_alpha) < 0.2 && AP_pT > 0.0 && AP_pT < 0.02 && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.5 && pTB < 0.5 && dot_product_dir_vertex > 0.9) // TRD photon conversion
                                //if(fabs(AP_alpha) < 0.2 && AP_pT > 0.0 && AP_pT < 0.02 && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.2 && pTB < 0.2 && dot_product_dir_vertex > 0.9) // TRD photon conversion
                                //if(fabs(AP_alpha) < 0.2 && AP_pT > 0.0 && AP_pT < 0.02 && CA*CB < 0.0 && pTA > 0.3 && pTB > 0.3 && pTA < 0.5 && pTB < 0.5 && dot_product_dir_vertex > 0.9) // TPC photon conversion
                                {
                                    //printf("-----> Found vertex for AP in event: %lld at radius: %4.3f, AP_pT: %4.3f, AP_alpha: %4.3f, pTA: %4.3f, pTB: %4.3f, dot: %4.3f, Inv_mass_AB: %4.3f, Energy_AB: %4.3f, Momentum_AB: %4.3f \n",Global_Event,radius_vertex,AP_pT,AP_alpha,pTA,pTB,dot_product_dir_vertex,Inv_mass_AB,Energy_AB,Momentum_AB);
                                    flag_found_good_AP_vertex = 1;
                                    if(graphics)
                                    {
                                        // Draw kalman tracks

                                        // XALEX
                                        Draw_Kalman_Helix_Tracks(i_track_A,kGreen,280.0,500.0); // 280.0, 500.0
                                        Draw_Kalman_Helix_Tracks(i_track_B,kRed,280.0,500.0);   // 280.0, 500.0

                                        // Draw kalman tracklets
                                        Draw_matched_Kalman_Tracklets(i_track_A);
                                        Draw_matched_Kalman_Tracklets(i_track_B);

                                        TEveP_sec_vertices    ->SetPoint(i_vertex_acc,vertex_point[0],vertex_point[1],vertex_point[2]);
                                        TEveP_photon_vertices ->SetPoint(i_vertex_acc,TV3_photon.X(),TV3_photon.Y(),TV3_photon.Z());
                                        //printf("i_comb: %d, A p,pT,pz: {%4.3f, %4.3f, %4.3f}, B p,pT,pz: {%4.3f, %4.3f, %4.3f}, AP pT, alpha: {%4.3f, %4.3f}, TV3_dirAB: {%4.3f, %4.3f, %4.3f} \n",i_comb,pA,pTA,pzA,pB,pTB,pzB,AP_pT,AP_alpha,TV3_dirAB[0],TV3_dirAB[1],TV3_dirAB[2]);


                                        TEveP_first_point_helix  ->SetPoint(i_vertex_acc,helix_pointA[0],helix_pointA[1],helix_pointA[2]);
                                        TEveP_second_point_helix ->SetPoint(i_vertex_acc,helix_pointAs[0],helix_pointAs[1],helix_pointAs[2]);
                                        TEveP_first_point_helix  ->SetPoint(i_vertex_acc+2,helix_pointB[0],helix_pointB[1],helix_pointB[2]);
                                        TEveP_second_point_helix ->SetPoint(i_vertex_acc+2,helix_pointBs[0],helix_pointBs[1],helix_pointBs[2]);

                                        TEveLine_mother.resize(i_vertex_acc+1);
                                        TEveLine_mother[i_vertex_acc] = new TEveLine();

                                        //TV3_sec_vertex, TV3_dirAB, TV3_prim_vertex // XLINE
                                        Double_t distance_to_prim_vertex = (TV3_prim_vertex - TV3_sec_vertex).Mag();

                                        TEveLine_mother[i_vertex_acc] ->SetNextPoint(vertex_point[0],vertex_point[1],vertex_point[2]);
                                        TEveLine_mother[i_vertex_acc] ->SetNextPoint(vertex_point[0] - distance_to_prim_vertex*TV3_dirAB[0],vertex_point[1] - distance_to_prim_vertex*TV3_dirAB[1],vertex_point[2] - distance_to_prim_vertex*TV3_dirAB[2]);

                                        i_vertex_acc++;
                                    }
                                }
#endif

                                if(pTA > 0.3 && pTB > 0.3 && CA*CB < 0.0 && (pzA + pzB) != 0.0)
                                {
                                    //printf("CA: %4.3f, CB: %4.3f \n",CA,CB);


                                    //printf("   ----> AP_alpha: %4.3f, AP_pT: %4.3f, AP_pTB: %4.3f, projA: %4.3f, projB: %4.3f \n",AP_alpha,AP_pT,AP_pTB,projA,projB);
                                    TH2D_AP_plot ->Fill(AP_alpha,AP_pT);
                                    if(i_radius < N_AP_radii) vec_TH2D_AP_plot_radius[i_radius] ->Fill(AP_alpha,AP_pT);
                                }
                            }

                        }

                    }
                    //-------------------------------------------------
                    i_vertex++;
                }
                i_comb++;
            }
        }
    }


    //printf("Number of secondary vertices found: %d \n",i_vertex);


    if(flag_TRD_TPC_tracks == 0) // TRD tracks
    {        
        //------------------------
        Int_t N_sec_vertices = (Int_t)vec_TV3_secondary_vertices.size();
        Double_t radius_sec_vertex;
        TVector3 TV3_diff_sec_vertices;
        TVector3 TV3_avg_sec_vertex;
        Float_t Arr_cluster_params[12];

        //Bool_t visited[N_sec_vertices];
        vector<Int_t> visited;
        visited.resize(N_sec_vertices);
        for(Int_t i_vis = 0; i_vis < N_sec_vertices; i_vis++) visited[i_vis] = 0;

        Int_t i_vertex_nucl_int = 0;
        for(Int_t i_sec_vtx_A = 0; i_sec_vtx_A < (N_sec_vertices - 1); i_sec_vtx_A++)
        {
            //printf("i_sec_vtx_A: %d \n",i_sec_vtx_A);
            radius_sec_vertex = vec_TV3_secondary_vertices[i_sec_vtx_A].Perp();
            if(radius_sec_vertex < 240.0 || radius_sec_vertex > 370.0) continue;
            if(visited[i_sec_vtx_A]) continue;
            visited[i_sec_vtx_A] = 1;
            Int_t N_close_vertex = 1;
            TV3_avg_sec_vertex = vec_TV3_secondary_vertices[i_sec_vtx_A]*1;
            vector<Int_t> idx_acc_vtx_B;
            idx_acc_vtx_B.clear();
            for(Int_t i_sec_vtx_B = (i_sec_vtx_A + 1); i_sec_vtx_B < N_sec_vertices; i_sec_vtx_B++)
            {
                if(visited[i_sec_vtx_B]) continue;
                TV3_diff_sec_vertices = vec_TV3_secondary_vertices[i_sec_vtx_A] - vec_TV3_secondary_vertices[i_sec_vtx_B];
                if(TV3_diff_sec_vertices.Mag() < 3.0)
                {
                    TV3_avg_sec_vertex += vec_TV3_secondary_vertices[i_sec_vtx_B];
                    N_close_vertex++;
                    visited[i_sec_vtx_B] = 1;
                    idx_acc_vtx_B.push_back(i_sec_vtx_B);
                }
            }

            if(N_close_vertex > 2)
            {
#if 0
                for(Int_t i_sec_vtx = 0; i_sec_vtx < N_sec_vertices; i_sec_vtx++)
                {
                    printf("i_sec_vtx: %d, visited: %d \n",i_sec_vtx,visited[i_sec_vtx]);
                }
#endif

                TV3_avg_sec_vertex *= 1.0/(Float_t)(N_close_vertex);


                // Determine TPC track(s) which are close to nuclear interaction vertex
                Float_t pathA_dca = -1.0;
                Float_t dcaAB_dca = -1.0;
                Float_t dcaAB_min        = 999.0;
                Float_t TOFsignal_min    = 999.0;
                Float_t Track_length_min = 999.0;
                Float_t TPCdEdx_min      = 999.0;
                Float_t dca_to_prim      = 999.0;
                Float_t pT_min           = 999.0;
                Float_t momentum_min     = 999.0;
                Float_t path_min         = 999.0;


                //-------------------------------
                // Get TPC track information
                UShort_t NumTracks = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
                Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
                Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
                Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
                Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
                Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
                Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
                Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
                Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
                TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
                UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
                UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
                UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
                Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
                Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
                Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
                Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();

                Float_t momentum        = TLV_part.P();
                Float_t eta_track       = TLV_part.Eta();
                Float_t pT_track        = TLV_part.Pt();
                Float_t theta_track     = TLV_part.Theta();
                Float_t phi_track       = TLV_part.Phi();
                //-------------------------------


                Int_t flag_close_TPC_track = 0;
                Int_t idx_close_TPC_track = -1;
                for(Int_t i_track = 0; i_track < NumTracks; i_track++)
                {
                    TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);
                    TPC_single_helix ->setHelix(TRD_ST_TPC_Track->getHelix_param(0),TRD_ST_TPC_Track->getHelix_param(1),TRD_ST_TPC_Track->getHelix_param(2),TRD_ST_TPC_Track->getHelix_param(3),TRD_ST_TPC_Track->getHelix_param(4),TRD_ST_TPC_Track->getHelix_param(5));
                    fHelixAtoPointdca(TV3_avg_sec_vertex,TPC_single_helix,pathA_dca,dcaAB_dca); // new helix to point dca calculation
                    if(dcaAB_dca < dcaAB_min)
                    {
                        dcaAB_min        = dcaAB_dca;
                        TOFsignal_min    = TOFsignal;
                        Track_length_min = Track_length;
                        TPCdEdx_min      = TPCdEdx;
                        dca_to_prim      = dca;
                        pT_min           = pT_track;
                        momentum_min     = momentum;
                        path_min         = pathA_dca;
                    }
                    if(dcaAB_dca < 3.0)
                    {
                        idx_close_TPC_track = i_track;
                        flag_close_TPC_track = 1;
                    }
                }

                //printf("dcaAB_min: %4.3f \n",dcaAB_min);


                //-----------------------------------
                // Store the kalman tracks which are used in nuclear event candidate
                vector<Int_t> vec_idx_kalman_tracks_nuclev;
                vec_idx_kalman_tracks_nuclev.push_back(vec_idx_kalman_sec_vert[0][i_sec_vtx_A]);
                vec_idx_kalman_tracks_nuclev.push_back(vec_idx_kalman_sec_vert[1][i_sec_vtx_A]);
                for(Int_t idx_B = 0; idx_B < (Int_t)idx_acc_vtx_B.size(); idx_B++)
                {
                    vec_idx_kalman_tracks_nuclev.push_back(vec_idx_kalman_sec_vert[0][idx_acc_vtx_B[idx_B]]);
                    vec_idx_kalman_tracks_nuclev.push_back(vec_idx_kalman_sec_vert[1][idx_acc_vtx_B[idx_B]]);
                }

                ///if(TV3_avg_sec_vertex.Perp() > 297.0 && flag_close_TPC_track) printf("%s ----> Nuclear interaction vertex at radius:  %s %4.3f cm, pos: {%4.3f, %4.3f, %4.3f} at event: %lld, i_vertex_nucl_int: %d, N_close_vertex: %d, path_min: %4.3f \n",KRED,KNRM,TV3_avg_sec_vertex.Perp(),TV3_avg_sec_vertex[0],TV3_avg_sec_vertex[1],TV3_avg_sec_vertex[2],Global_Event,i_vertex_nucl_int,N_close_vertex,path_min);
                Float_t nuclev_bitmap = Calc_nuclev_bitmap(vec_idx_kalman_tracks_nuclev); // Calculates the number of shared/independent tracklets per layer

                // Read bitmap
                Int_t Arr_tracklets_layer[6] = {0};
                for(Int_t i_lay = 0; i_lay < 6; i_lay++)
                {
                    Int_t N_tracklets_layer = 0;
                    for(Int_t i_bit = 0; i_bit < 5; i_bit++) // asume only 5 bits are filled -> 32 tracklets per layer at maximum
                    {
                        if(((Int_t)nuclev_bitmap >> (i_bit + i_lay*5)) & 1) // read the bits
                        {
                            N_tracklets_layer |= 1 << i_bit; // set bit i_bit to 1
                        }
                    }
                    Arr_tracklets_layer[i_lay] = N_tracklets_layer;
                    //printf(" --> bitmap read: %d, i_lay: %d, N_tracklets: %d \n",(Int_t)nuclev_bitmap,i_lay,N_tracklets_layer);
                }
                //-----------------------------------

                //-----------------------------------
                // Toplogy selection for nuclear events
                //if(Arr_tracklets_layer[5] > 3 || Arr_tracklets_layer[4] > 3 || Arr_tracklets_layer[3] > 3)
#if 0
                if(
                   (Arr_tracklets_layer[5] > 2 && Arr_tracklets_layer[4] > 1 && Arr_tracklets_layer[3] > 0)
                   || (Arr_tracklets_layer[4] > 2 && Arr_tracklets_layer[3] > 1 && Arr_tracklets_layer[2] > 0)
                   || (Arr_tracklets_layer[3] > 2 && Arr_tracklets_layer[2] > 1 && Arr_tracklets_layer[1] > 0)
                   || (Arr_tracklets_layer[2] > 2 && Arr_tracklets_layer[1] > 1 && Arr_tracklets_layer[0] > 0)
                   //
                   || (Arr_tracklets_layer[5] > 2 && Arr_tracklets_layer[3] > 1 && Arr_tracklets_layer[2] > 0)
                   || (Arr_tracklets_layer[4] > 2 && Arr_tracklets_layer[2] > 1 && Arr_tracklets_layer[1] > 0)
                   || (Arr_tracklets_layer[3] > 2 && Arr_tracklets_layer[1] > 1 && Arr_tracklets_layer[0] > 0)
                   //
                   || (Arr_tracklets_layer[4] > 2 && Arr_tracklets_layer[2] > 1 && Arr_tracklets_layer[1] > 0)
                   || (Arr_tracklets_layer[4] > 2 && Arr_tracklets_layer[1] > 1 && Arr_tracklets_layer[0] > 0)
                   || (Arr_tracklets_layer[3] > 2 && Arr_tracklets_layer[2] > 1 && Arr_tracklets_layer[0] > 0)
                  )
#endif
                    if(
                       (Arr_tracklets_layer[5] > 2 && Arr_tracklets_layer[4] > 2 && Arr_tracklets_layer[3] > 2)
                      )
                    {
                        if(TV3_avg_sec_vertex.Perp() > 297.0 && flag_close_TPC_track)
                        {
                            printf("%s ----> Nuclear interaction vertex at radius:  %s %4.3f cm, pos: {%4.3f, %4.3f, %4.3f} at event: %lld, i_vertex_nucl_int: %d, N_close_vertex: %d, path_min: %4.3f, trkl: {%d, %d, %d, %d, %d, %d} \n",KRED,KNRM,TV3_avg_sec_vertex.Perp(),TV3_avg_sec_vertex[0],TV3_avg_sec_vertex[1],TV3_avg_sec_vertex[2],Global_Event,i_vertex_nucl_int,N_close_vertex,path_min,Arr_tracklets_layer[0],Arr_tracklets_layer[1],Arr_tracklets_layer[2],Arr_tracklets_layer[3],Arr_tracklets_layer[4],Arr_tracklets_layer[5]);
                        }
                    }

                // Ntuple for nculear interaction event candidates
                Arr_cluster_params[0]	= (Float_t)TV3_avg_sec_vertex[0];
                Arr_cluster_params[1]	= (Float_t)TV3_avg_sec_vertex[1];
                Arr_cluster_params[2]	= (Float_t)TV3_avg_sec_vertex[2];
                Arr_cluster_params[3]	= (Float_t)N_close_vertex;
                Arr_cluster_params[4]	= (Float_t)dcaAB_min;
                Arr_cluster_params[5]	= (Float_t)TOFsignal_min;
                Arr_cluster_params[6]	= (Float_t)Track_length_min;
                Arr_cluster_params[7]	= (Float_t)TPCdEdx_min;
                Arr_cluster_params[8]	= (Float_t)dca_to_prim;
                Arr_cluster_params[9]	= (Float_t)pT_min;
                Arr_cluster_params[10]	= (Float_t)momentum_min;
                Arr_cluster_params[11]	= (Float_t)nuclev_bitmap;

                //NT_secondary_vertex_cluster->Fill(Arr_cluster_params);

                TRD_Nuclear_interaction = TRD_Self_Event -> createNucInteraction();

                TRD_Nuclear_interaction -> set_avg_sec_vertex((Float_t)TV3_avg_sec_vertex[0],(Float_t)TV3_avg_sec_vertex[1],(Float_t)TV3_avg_sec_vertex[2]);
                TRD_Nuclear_interaction -> set_N_close_vertex((Float_t)N_close_vertex);
                TRD_Nuclear_interaction -> set_dcaAB_min((Float_t)dcaAB_min);
                TRD_Nuclear_interaction -> set_dcaAB_min((Float_t)dcaAB_min);
                TRD_Nuclear_interaction -> set_TOFsignal_min((Float_t)TOFsignal_min);
                TRD_Nuclear_interaction -> set_Track_length_min((Float_t)Track_length_min);
                TRD_Nuclear_interaction -> set_TPCdEdx_min((Float_t)TPCdEdx_min);
                TRD_Nuclear_interaction -> set_pT_min((Float_t)pT_min);
                TRD_Nuclear_interaction -> set_momentum_min((Float_t)momentum_min);
                TRD_Nuclear_interaction -> set_nuclev_bitmap((Float_t)nuclev_bitmap);

                if (flag_TRD_TPC_tracks == 0) //Kalman
                {
                    for (Int_t i_track = 0; i_track < (Int_t)vec_idx_kalman_tracks_nuclev.size(); i_track++)
                    {
                        //printf("test1 \n");
                        TRD_Kalman_track = TRD_Nuclear_interaction -> createKalman_Track();
                        //printf("test2 \n");

                            TRD_Kalman_track ->setKalmanHelix_param(vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(0),vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(1),
                                vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(2),vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(3),
                                vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(4),vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(5)); 
                    }
                }
                if (flag_TRD_TPC_tracks == 1) //TPC
                {
                    for (Int_t i_track = 0; i_track < (Int_t)vec_idx_kalman_tracks_nuclev.size(); i_track++)
                    {
                        //printf("test1 \n");
                        TPC_track = TRD_Nuclear_interaction -> createTPC_Track();
                        //printf("test2 \n");

                            TPC_track ->setHelix(vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(0),vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(1),
                                vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(2),vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(3),
                                vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(4),vec_helices[vec_idx_kalman_tracks_nuclev[i_track]]->getHelix_param(5)); 
                    }
                }

                //-----------------------------------




#if defined(USEEVE)
                if(graphics)
                {
                    TEveP_nucl_int_vertices ->SetPoint(i_vertex_nucl_int,TV3_avg_sec_vertex[0],TV3_avg_sec_vertex[1],TV3_avg_sec_vertex[2]);

                    // Draw kalman tracks which contributed to nuclear interaction vertex
                    Draw_Kalman_Helix_Tracks(vec_idx_kalman_sec_vert[0][i_sec_vtx_A],kGreen+1,280.0,500.0);
                    Draw_Kalman_Helix_Tracks(vec_idx_kalman_sec_vert[1][i_sec_vtx_A],kGreen+1,280.0,500.0);

                    // Draw kalman tracklets
                    Draw_matched_Kalman_Tracklets(vec_idx_kalman_sec_vert[0][i_sec_vtx_A]);
                    Draw_matched_Kalman_Tracklets(vec_idx_kalman_sec_vert[1][i_sec_vtx_A]);

                    for(Int_t idx_B = 0; idx_B < (Int_t)idx_acc_vtx_B.size(); idx_B++)
                    {
                        Draw_Kalman_Helix_Tracks(vec_idx_kalman_sec_vert[0][idx_acc_vtx_B[idx_B]],kMagenta+1,280.0,500.0);
                        Draw_Kalman_Helix_Tracks(vec_idx_kalman_sec_vert[1][idx_acc_vtx_B[idx_B]],kMagenta+1,280.0,500.0);

                        // Draw kalman tracklets
                        Draw_matched_Kalman_Tracklets(vec_idx_kalman_sec_vert[0][idx_acc_vtx_B[idx_B]]);
                        Draw_matched_Kalman_Tracklets(vec_idx_kalman_sec_vert[1][idx_acc_vtx_B[idx_B]]);
                    }
                    if(idx_close_TPC_track >= 0) Draw_TPC_track(idx_close_TPC_track,kAzure-2,3,1000.0);

                    //Float_t pathA_dca, dcaAB_dca;
                    //fHelixAtoPointdca(TV3_EventVertex,vec_helices[i_track],pathA_dca,dcaAB_dca); // new helix to point dca calculation
                    //Draw_TPC_track(idx_matched_TPC_track,kAzure-2,3,1000.0);
                }
#endif
                i_vertex_nucl_int++;
            }
            //if(N_close_vertex > 3) printf("%s ----> N_close_vertex: %d %s, event: %lld \n",KRED,N_close_vertex,KNRM,Global_Event);

        }
        //------------------------
        //cout << "End of nuclear interaction finder" << endl;


#if defined(USEEVE)
        if(graphics)
        {
            TEveP_sec_vertices  ->SetMarkerSize(3);
            TEveP_sec_vertices  ->SetMarkerStyle(20);
            TEveP_sec_vertices  ->SetMarkerColor(kRed);
            gEve->AddElement(TEveP_sec_vertices);

            TEveP_photon_vertices  ->SetMarkerSize(3);
            TEveP_photon_vertices  ->SetMarkerStyle(20);
            TEveP_photon_vertices  ->SetMarkerColor(kCyan+2);
            gEve->AddElement(TEveP_photon_vertices);

            TEveP_nucl_int_vertices  ->SetMarkerSize(5);
            TEveP_nucl_int_vertices  ->SetMarkerStyle(20);
            TEveP_nucl_int_vertices  ->SetMarkerColor(kBlue);
            gEve->AddElement(TEveP_nucl_int_vertices);

            TEveP_close_TPC_photon  ->SetMarkerSize(5);
            TEveP_close_TPC_photon  ->SetMarkerStyle(20);
            TEveP_close_TPC_photon  ->SetMarkerColor(kGreen);
            //gEve->AddElement(TEveP_close_TPC_photon);

            TEveP_first_point_helix  ->SetMarkerSize(3);
            TEveP_first_point_helix  ->SetMarkerStyle(20);
            TEveP_first_point_helix  ->SetMarkerColor(kBlue);
            //gEve->AddElement(TEveP_first_point_helix);

            TEveP_second_point_helix  ->SetMarkerSize(3);
            TEveP_second_point_helix  ->SetMarkerStyle(20);
            TEveP_second_point_helix  ->SetMarkerColor(kMagenta);
            //gEve->AddElement(TEveP_second_point_helix);

            for(Int_t i_mother = 0; i_mother < (Int_t)TEveLine_mother.size(); i_mother++)
            {
                TEveLine_mother[i_mother]    ->SetLineStyle(9);
                TEveLine_mother[i_mother]    ->SetLineWidth(3);
                TEveLine_mother[i_mother]    ->SetMainColor(kMagenta);
                TEveLine_mother[i_mother]    ->SetMainAlpha(1.0);
                gEve->AddElement(TEveLine_mother[i_mother]);
            }

            gEve->Redraw3D(kTRUE);
        }
#endif
    }

    Tree_TRD_Self_Event_out ->Fill();


    return flag_found_good_AP_vertex;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_Kalman_helix_params(vector<vector<Double_t>> mHelices_kalman_in)
{
    //printf("Ali_TRD_ST_Analyze::set_Kalman_helix_params");
    mHelices_kalman = mHelices_kalman_in;

    vec_helices_TRD.clear();
    vec_helices_TRD.resize((Int_t)mHelices_kalman.size());

    for(Int_t i_track = 0; i_track < (Int_t)mHelices_kalman.size(); i_track++)
    {
        vec_helices_TRD[i_track] = new Ali_Helix();
        vec_helices_TRD[i_track] ->setHelix(mHelices_kalman[i_track][0],mHelices_kalman[i_track][1],mHelices_kalman[i_track][2],mHelices_kalman[i_track][3],mHelices_kalman[i_track][4],mHelices_kalman[i_track][5]);
        //printf("i_KF_track: %d, par: {%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f} \n",i_track,mHelices_kalman[i_track][0],mHelices_kalman[i_track][1],mHelices_kalman[i_track][2],mHelices_kalman[i_track][3],mHelices_kalman[i_track][4],mHelices_kalman[i_track][5]);
    }

}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_Kalman_chi_2(vector<Double_t> mChi_2s_kalman_in) //mChi_2s_tracker_KF
{
    //printf("Ali_TRD_ST_Analyze::set_Kalman_helix_params");
    mChi_2s_kalman = mChi_2s_kalman_in;

}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::set_TPC_helix_params(Long64_t i_event)
{
    //printf("Ali_TRD_ST_Analyze::set_TPC_helix_params");
    vec_helices_TPC.clear();


    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks   = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
    Float_t  V0MEq       = TRD_ST_Event ->getcent_class_V0MEq();
    //--------------------------------------------------



    //--------------------------------------------------
    // TPC track loop
    mHelices_TPC.clear();
    PID_params_TPC.clear();
    mHelices_TPC.resize(NumTracks);
    PID_params_TPC.resize(NumTracks);
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

        Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
        Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
        UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
        UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
        UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
        Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
        Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
        Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();
        Float_t charge = 1.0;
        if(dca < 0.0) charge = -1.0;

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t pZ_track        = TLV_part.Pz();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();

        vec_helices_TPC.push_back(new Ali_Helix());
        vec_helices_TPC[(Int_t)vec_helices_TPC.size()-1] ->setHelix(TRD_ST_TPC_Track->getHelix_param(0),TRD_ST_TPC_Track->getHelix_param(1),TRD_ST_TPC_Track->getHelix_param(2),TRD_ST_TPC_Track->getHelix_param(3),TRD_ST_TPC_Track->getHelix_param(4),TRD_ST_TPC_Track->getHelix_param(5));
        mHelices_TPC[i_track].resize(9);
        PID_params_TPC[i_track].resize(3);
        PID_params_TPC[i_track][0] = TPCdEdx;
        PID_params_TPC[i_track][1] = dca;
        PID_params_TPC[i_track][2] = momentum;
        for(Int_t i_param = 0; i_param < 9; i_param++)
        {
            mHelices_TPC[i_track][i_param] = TRD_ST_TPC_Track->getHelix_param(i_param);
            //printf("i_track: %d, i_param: %d, val: %4.3f \n",i_track,i_param,mHelices_TPC[i_track][i_param]);
        }
        mHelices_TPC[i_track][6] = pT_track;
        mHelices_TPC[i_track][7] = pZ_track;
    }

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_Kalman_TRD_tracklets(vector< vector<Ali_TRD_ST_Tracklets*> > vec_kalman_TRD_trackets_in)
{
    vec_kalman_TRD_trackets = vec_kalman_TRD_trackets_in;
    //printf("size: %d \n",(Int_t)vec_kalman_TRD_trackets.size());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Calc_Kalman_efficiency()
{
    // Loop over TPC tracks
    Int_t N_TPC_tracks = (Int_t)matched_tracks.size();
	Int_t i_TPC_mom_list=0;
    for(Int_t i_TPC_track = 0; i_TPC_track < N_TPC_tracks; i_TPC_track++)
    {
        //--------------------
        // Get information from TPC track
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_TPC_mom_list);
		TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
		Float_t momentum        = TLV_part.P();
        while(momentum<0.3)
		{	
			i_TPC_mom_list++;
			TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_TPC_mom_list);
			TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
			momentum        = TLV_part.P();
        
		}	
		
        Float_t pT_track        = TLV_part.Pt();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        Float_t charge = 1.0;
        if(dca < 0.0) charge = -1.0;
        //--------------------


        Int_t N_layers_matched_all    = 0;
        Int_t N_good_TPC_TRD_tracklets = 0;
        vector<Int_t> vec_N_layers_matched;
        Int_t N_Kalman_tracks = (Int_t)vec_kalman_TRD_trackets.size();
        vec_N_layers_matched.resize(N_Kalman_tracks);
        vector<Int_t> index_kalman_track_matched;
        index_kalman_track_matched.clear();
        for(Int_t i_kalm_track = 0; i_kalm_track < N_Kalman_tracks; i_kalm_track++)
        {
            vec_N_layers_matched[i_kalm_track] = 0;
        }
        //printf("N_layers TPC: %d \n",(Int_t)matched_tracks[i_TPC_track].size());

        // Loop over TRD layers matched with TPC track
        for(Int_t i_layer_TPC = 0; i_layer_TPC < (Int_t)matched_tracks[i_TPC_track].size(); i_layer_TPC++)
        {
            //matched_tracks
            if(matched_tracks[i_TPC_track][i_layer_TPC] == NULL) continue;
            TVector3 TPC_TV3_offset = matched_tracks[i_TPC_track][i_layer_TPC]->get_TV3_offset();
            N_good_TPC_TRD_tracklets++;

            // Loop over kalman matched TRD tracklets
            for(Int_t i_kalm_track = 0; i_kalm_track < N_Kalman_tracks; i_kalm_track++)
            {
                Int_t flag_match = 0;
                Int_t N_good_kalman_tracklets = 0;
                //printf("N_layers Kalman: %d \n",(Int_t)vec_kalman_TRD_trackets[i_kalm_track].size());

                // Loop over TRD layers matched with Kalman track
                for(Int_t i_layer = i_layer_TPC; i_layer <= i_layer_TPC; i_layer++)
                {
                    if(vec_kalman_TRD_trackets[i_kalm_track][i_layer] == NULL) continue;
                    N_good_kalman_tracklets++;
                    TVector3 kalman_TV3_offset = vec_kalman_TRD_trackets[i_kalm_track][i_layer] ->get_TV3_offset();

                    TVector3 TV3_diff =  kalman_TV3_offset - TPC_TV3_offset;
                    if(TV3_diff.Mag() < 0.1)
                    {
                        N_layers_matched_all++;
                        vec_N_layers_matched[i_kalm_track]++;
                        flag_match = 1;
                    }
                }
                if(vec_N_layers_matched[i_kalm_track] == 1 && flag_match) index_kalman_track_matched.push_back(i_kalm_track);
            }
        }

        // Calculate actuall matching (only tracklets matched) and tracking (inlcuding momentum) efficieny
        if(N_good_TPC_TRD_tracklets > 3) // Use only TPC tracks with a certain amount of TRD matched tracklets
        {
            // Check which kalman track has the best matching
            Int_t index_kalman_track_best = -1;
            Double_t matches_percent_best = 0.0;
            for(Int_t index_kalman_track = 0; index_kalman_track < (Int_t)index_kalman_track_matched.size(); index_kalman_track++)
            {
                // Fraction of matched Kalman to TPC track TRD tracklets
                Double_t matches_percent = 100.0*(Double_t)vec_N_layers_matched[index_kalman_track_matched[index_kalman_track]]/(Double_t)N_good_TPC_TRD_tracklets;
                if(matches_percent > matches_percent_best)
                {
                    matches_percent_best    = matches_percent;
                    index_kalman_track_best = index_kalman_track;
                }
            }

            // No matching found
            if(index_kalman_track_best < 0)
            {
                tp_efficiency_all_vs_pT      ->Fill(pT_track,0.0);
                tp_efficiency_matching_vs_pT ->Fill(pT_track,0.0);
            }
            else // Good matching found
            {
                for(Int_t index_kalman_track = index_kalman_track_best; index_kalman_track <= index_kalman_track_best; index_kalman_track++)
                {
                    // Fraction of matched Kalman to TPC track TRD tracklets
                    Int_t N_layers_matched_kalman = vec_N_layers_matched[index_kalman_track_matched[index_kalman_track]];
                    Double_t matches_percent = 100.0*(Double_t)N_layers_matched_kalman/(Double_t)N_good_TPC_TRD_tracklets;

                    if(matches_percent >= 50.0)
                    {
                        Int_t i_kalm_track = index_kalman_track_matched[index_kalman_track];

                        Double_t pT_kalman   = mHelices_kalman[i_kalm_track][6]; // pT
                        Double_t pz_kalman   = mHelices_kalman[i_kalm_track][7]; // pz
                        Double_t curv_kalman = mHelices_kalman[i_kalm_track][4]; // curvature
                        Double_t delta_pT    = pT_track - pT_kalman;
                        if(pT_track != 0.0)
                        {
                            Double_t rel_delta_pT = fabs(delta_pT)/pT_track;

                            //printf("%s i_TPC_track: %d %s matched to i_kalm_track: %d, %4.2f%% tracklets matched, delta_pT: %4.3f, pT(TPC): %4.3f, pT(Kalman): %4.3f, N_TPC_tracklets: %d \n",KGRN,i_TPC_track,KNRM,i_kalm_track,matches_percent,delta_pT,pT_track,pT_kalman,N_good_TPC_TRD_tracklets);

                            if(rel_delta_pT < 0.5)
                            {
                                tp_efficiency_all_vs_pT ->Fill(pT_track,1.0);
                            }
                            else
                            {
                                tp_efficiency_all_vs_pT ->Fill(pT_track,0.0);
                            }
                        }
                        tp_efficiency_matching_vs_pT ->Fill(pT_track,1.0);
                        //printf("N_layers_matched_kalman: %d, q*pT: %4.3f, delta_pT: %4.3f \n",N_layers_matched_kalman,charge*pT_track,delta_pT);
                        if(pT_track > 0.0 && pT_track < 100.0 && delta_pT < 10.0)
                        {
                            vec_h2D_delta_pT_all_vs_pT[N_layers_matched_kalman-1] ->Fill(charge*pT_track,delta_pT);
                        }
                    }
                    else
                    {
                        tp_efficiency_all_vs_pT      ->Fill(pT_track,0.0);
                        tp_efficiency_matching_vs_pT ->Fill(pT_track,0.0);
                    }
                }
            }
        }
		i_TPC_mom_list++;
    }

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Match_kalman_tracks_to_TPC_tracks(Int_t graphics, Int_t draw_matched_TPC_track, Int_t draw_matched_TRD_track, Int_t color)
{
    //printf("size: %d \n",(Int_t)vec_kalman_TRD_trackets.size());
    // Loop over kalman matched TRD tracklets
    for(Int_t i_kalm_track = 0; i_kalm_track < (Int_t)vec_kalman_TRD_trackets.size(); i_kalm_track++)
    {
        vector<Int_t> vec_N_layers_matched;
        Int_t N_TPC_tracks = (Int_t)matched_tracks.size();
        vec_N_layers_matched.resize(N_TPC_tracks);
        Int_t N_layers_matched_all = 0;
        vector<Int_t> index_TPC_track_matched;
        for(Int_t i_TPC_track = 0; i_TPC_track < N_TPC_tracks; i_TPC_track++)
        {
            vec_N_layers_matched[i_TPC_track] = 0;
        }
        index_TPC_track_matched.clear();
        vec_N_layers_matched.clear();

        //printf("layers: %d \n",(Int_t)vec_kalman_TRD_trackets[i_kalm_track].size());
        Int_t N_good_kalman_tracklets = 0;
        for(Int_t i_layer = 0; i_layer < (Int_t)vec_kalman_TRD_trackets[i_kalm_track].size(); i_layer++)
        {
            if(vec_kalman_TRD_trackets[i_kalm_track][i_layer] == NULL) continue;
            N_good_kalman_tracklets++;
            TVector3 kalman_TV3_offset = vec_kalman_TRD_trackets[i_kalm_track][i_layer] ->get_TV3_offset();
            //printf("i_kalman_tracklet: %d, i_layer: %d, offset: {%4.3f, %4.3f, %4.3f} \n",i_kalm_track,i_layer,kalman_TV3_offset.X(),kalman_TV3_offset.Y(),kalman_TV3_offset.Z());

            // Loop over TPC tracks
            for(Int_t i_TPC_track = 0; i_TPC_track < N_TPC_tracks; i_TPC_track++)
            {
                Int_t flag_match = 0;
                for(Int_t i_layer_TPC = 0; i_layer_TPC < (Int_t)matched_tracks[i_TPC_track].size(); i_layer_TPC++)
                {
                    //matched_tracks
                    if(matched_tracks[i_TPC_track][i_layer_TPC] == NULL) continue;
                    TVector3 TPC_TV3_offset = matched_tracks[i_TPC_track][i_layer_TPC]->get_TV3_offset();
                    TVector3 TV3_diff =  kalman_TV3_offset - TPC_TV3_offset;
                    if(TV3_diff.Mag() < 0.1)
                    {
                        vec_N_layers_matched[i_TPC_track]++;
                        N_layers_matched_all++;
                        flag_match = 1;
                        //printf("i_kalm_track: %d, matched: %d \n",i_kalm_track,N_layers_matched_all);
                    }
                }
                if(vec_N_layers_matched[i_TPC_track] == 1 && flag_match) index_TPC_track_matched.push_back(i_TPC_track);
            }
        }

        if(N_layers_matched_all >= 2)
        {
            //printf("Kalman track %d with %d total matches \n",i_kalm_track,N_layers_matched_all);
            for(Int_t index_TPC_track = 0; index_TPC_track < (Int_t)index_TPC_track_matched.size(); index_TPC_track++)
            {
                Double_t matches_percent = 100.0*(Double_t)vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]]/(Double_t)N_good_kalman_tracklets;
                //printf("  --> TPC_track: %d, N_good_kalman_tracklets: %d N_matches: %d, %4.1f%% matches \n",index_TPC_track_matched[index_TPC_track],N_good_kalman_tracklets,vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]],matches_percent);

                //if(matches_percent >= 90.0 && vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]] >= 5)
                if(matches_percent >= 74.0 && vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]] >= 4)
                {
                    Int_t idx_matched_TPC_track = vec_idx_matched_TPC_track[index_TPC_track_matched[index_TPC_track]];

                    TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(idx_matched_TPC_track);

                    TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
                    Float_t momentum        = TLV_part.P();
                    Float_t pT_track        = TLV_part.Pt();
                    Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex

                    Double_t pT_kalman   = mHelices_kalman[i_kalm_track][6]; // pT
                    Double_t pz_kalman   = mHelices_kalman[i_kalm_track][7]; // pz
                    Double_t curv_kalman = mHelices_kalman[i_kalm_track][4]; // curvature

                    if(graphics)
                    {
                        if(draw_matched_TPC_track) Draw_TPC_track(idx_matched_TPC_track,kAzure-2,3,1000.0);
                        //if(draw_matched_TPC_track) Draw_TPC_track(idx_matched_TPC_track,kAzure-2,3,280.0,500.0,1000.0);
                        if(draw_matched_TRD_track) Draw_Kalman_Helix_Tracks(i_kalm_track,color,0.0,500.0);
                    }
                    //printf("      >>>>> %s idx (Kalman): %d %s, %s idx (TPC): %d %s, pT (TPC): %4.3f, pT (kalman): %4.3f \n",KRED,i_kalm_track,KNRM,KBLU,idx_matched_TPC_track,KNRM,pT_track,pT_kalman);

                    TH2D_pT_TPC_vs_Kalman                                       ->Fill(-TMath::Sign(1,curv_kalman)*pT_kalman,TMath::Sign(1,dca)*pT_track);
                    vec_TH2D_pT_TPC_vs_Kalman[N_good_kalman_tracklets]          ->Fill(-TMath::Sign(1,curv_kalman)*pT_kalman,TMath::Sign(1,dca)*pT_track);
                    if((-TMath::Sign(1,curv_kalman)*pT_kalman) != 0.0 && pT_track != 0.0) vec_TH2D_one_over_pT_TPC_vs_Kalman[N_good_kalman_tracklets] ->Fill(1.0/(-TMath::Sign(1,curv_kalman)*pT_kalman),1.0/(TMath::Sign(1,dca)*pT_track));
                    //printf("N_good_kalman_tracklets: %d \n",N_good_kalman_tracklets);
                }
            }
        }

    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_single_helix_params(vector<Double_t> vec_params)
{
    for(Int_t i_param = 0; i_param < 6; i_param++)
    {
        aliHelix_params[i_param] = vec_params[i_param];
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Evaluate(Double_t t,Double_t r[3])  //radius vector
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives at given phase
  //--------------------------------------------------------------------
  float phase=aliHelix_params[4]*t+aliHelix_params[2];
  Double_t sn=sinf(phase), cs=cosf(phase);
  //  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = aliHelix_params[5] + sn/aliHelix_params[4];
  r[1] = aliHelix_params[0] - cs/aliHelix_params[4];
  r[2] = aliHelix_params[1] + aliHelix_params[3]*t;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_Kalman_Helix_Tracks(Int_t n_track, Int_t color, Double_t low_R, Double_t high_R)
{
    // n_track = -1 -> all tracks drawn
    Int_t i_track_start = 0;
    Int_t i_track_stop  = (Int_t)mHelices_kalman.size();
    if(n_track != -1)
    {
        i_track_start = n_track;
        i_track_stop  = n_track + 1;
    }

    printf("Ali_TRD_ST_Analyze::Draw_Kalman_Helix_Tracks, track: %d, start: %d, stop: %d \n",n_track,i_track_start,i_track_stop);

#if 1 // for calib tracks

    TVector3 vec_TV3_tracklet_vectors, TV3_tracklet_off_vector, vec_dir_vec_circle;
    Float_t pathA_dca, dcaAB_dca;
    Double_t track_posA[3], track_posB[3];
    Int_t detector[6] = {0};
#endif


    for(Int_t i_track = i_track_start; i_track < i_track_stop; i_track++)
    {
        

#if 1 //for calib tracks
        Double_t arr_path_layer[6] = {-999.0};
        Double_t arr_dca_layer[6]  = {-999.0};
        Double_t N_good_tracklets = 0.0;
        Double_t average_dca = 0.0;

        set_single_helix_params(mHelices_kalman[i_track]);
        Double_t pT_KF = mHelices_kalman[i_track][6]; // pT

        //if(pT_KF < 1.2) continue;

        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            //printf("i_layer: %d \n",i_layer);

            if(vec_kalman_TRD_trackets[i_track][i_layer] == NULL ) continue;

            TV3_tracklet_off_vector      = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset();
            detector[i_layer]        = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TRD_det();

            fHelixAtoPointdca(TV3_tracklet_off_vector,vec_helices_TRD[i_track],pathA_dca,dcaAB_dca); // new helix to point dca calculation
            arr_path_layer[i_layer] = pathA_dca;
            arr_dca_layer[i_layer]  = dcaAB_dca;
            //printf("   ---> pos trkl: {%4.3f, %4.3f, %4.3f}, pathA_dca: %4.3f, dcaAB_dca: %4.3f \n",TV3_tracklet_off_vector.X(),TV3_tracklet_off_vector.Y(),TV3_tracklet_off_vector.Z(),pathA_dca,dcaAB_dca);
            average_dca += dcaAB_dca;
            N_good_tracklets += 1.0;
        }

        if(N_good_tracklets > 0.0) average_dca /= N_good_tracklets;
        //if(average_dca > 3.0) continue;


#endif


#if 1
        //-------old part
        set_single_helix_params(mHelices_kalman[i_track]);
        Double_t pT_kalman   = mHelices_kalman[i_track][6]; // pT
        //if(pT_kalman < 0.5) continue;
        //-----------------
#endif

        Double_t track_pos[3];

#if defined(USEEVE)
        vec_TPL3D_helix_kalman.resize(i_track+1);
        vec_TPL3D_helix_kalman_hull.resize(i_track+1);
        vec_TPL3D_helix_kalman_inner.resize(i_track+1);
        vec_TPL3D_helix_kalman[i_track] = new TEveLine();
        vec_TPL3D_helix_kalman_hull[i_track] = new TEveLine();
        vec_TPL3D_helix_kalman_inner[i_track] = new TEveLine();
#endif

        Float_t pathA_dca, dcaAB_dca;
        fHelixAtoPointdca(TV3_EventVertex,vec_helices_TRD[i_track],pathA_dca,dcaAB_dca); // new helix to point dca calculation
        //if(dcaAB_dca > 50.0) continue;

        Double_t radius_helix = 0.0;
        for(Double_t track_path = 75.0; track_path > -350; track_path -= 1.0)
        {
            Evaluate(track_path,track_pos);
            radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
            //printf("i_track: %d, track_path: %4.3f, pos: {%4.3f, %4.3f, %4.3f}, radius_helix: %4.3f \n",i_track,track_path,track_pos[0],track_pos[1],track_pos[2],radius_helix);
            if(isnan(radius_helix)) break;
            if(radius_helix > high_R) break; // 500.0
            if(radius_helix < low_R) break; // 280.0
            if(fabs(track_pos[2]) > 430.0) break;
            //if(radius_helix > 80.0)
            {
#if defined(USEEVE)
                vec_TPL3D_helix_kalman[i_track]        ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
                vec_TPL3D_helix_kalman_hull[i_track]   ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
#endif
            }
        }



#if defined(USEEVE)
        HistName = "track kalman ";
        HistName += i_track;
        vec_TPL3D_helix_kalman[i_track]    ->SetName(HistName.Data());
        vec_TPL3D_helix_kalman[i_track]    ->SetLineStyle(1);
        vec_TPL3D_helix_kalman[i_track]    ->SetLineWidth(3);
        vec_TPL3D_helix_kalman[i_track]    ->SetMainColor(color);
        vec_TPL3D_helix_kalman[i_track]    ->SetMainAlpha(1.0);

        HistName = "track kalman (h) ";
        HistName += i_track;
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetName(HistName.Data());
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetLineStyle(1);
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetLineWidth(8);
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetMainColor(kOrange);
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetMainAlpha(0.3);

        gEve->AddElement(vec_TPL3D_helix_kalman[i_track]);
        gEve->AddElement(vec_TPL3D_helix_kalman_hull[i_track]);
#endif
    }

#if defined(USEEVE)
    gEve->Redraw3D(kTRUE);
#endif
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Init_tree(TString SEList)
{
    printf("Ali_TRD_ST_Analyze::Init_tree \n");
    TString pinputdir = input_dir;

    in_list_name = SEList;
    SEList = input_dir_lists + SEList;

    TRD_ST_Tracklet   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event      = new Ali_TRD_ST_Event();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( TRD_ST_TREE.Data(), TRD_ST_TREE.Data() );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = pinputdir+addfile;
                    input_SE ->AddFile(addfile.Data(),-1, TRD_ST_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( TRD_ST_BRANCH, &TRD_ST_Event );
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
        }
    }

    file_entries_total = input_SE->GetEntries();
    N_Events = file_entries_total;
    cout << "Total number of events in tree: " << file_entries_total << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::Loop_event(Long64_t i_event, Int_t graphics)
{
    //printf("Ali_TRD_ST_Analyze::Loop_event \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event


    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
    Global_Event         = i_event;
    Global_RunID         = TRD_ST_Event ->getid();
    TV3_EventVertex.SetXYZ(EventVertexX,EventVertexY,EventVertexZ);
    Float_t  V0MEq                = TRD_ST_Event ->getcent_class_V0MEq();

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_primary_vertex ->SetPoint(0,EventVertexX,EventVertexY,EventVertexZ);
        TEveP_primary_vertex ->SetMarkerStyle(20);
        TEveP_primary_vertex ->SetMarkerSize(4);
        TEveP_primary_vertex ->SetMarkerColor(kYellow+1);
        //gEve ->AddElement(TEveP_primary_vertex);
    }
#endif
    //--------------------------------------------------

    //printf("Event: %lld, TPC tracks: %d, TRD tracklets: %d \n",i_event,NumTracks,NumTracklets);

    //--------------------------------------------------
    // TPC track loop
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

        Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
        Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
        UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
        UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
        UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
        Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
        Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
        Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();
    }
    //--------------------------------------------------



    //--------------------------------------------------
    // TRD tracklet loop
    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t    i_det;

    Tracklets = new Ali_TRD_ST_Tracklets*[NumTracklets];
    Number_Tracklets=NumTracklets;
	
    for(Int_t i_tracklet = 0; i_tracklet < NumTracklets; i_tracklet++)
    {
        TRD_ST_Tracklet = TRD_ST_Event    ->getTracklet(i_tracklet);
        TRD_ST_Tracklet->set_TRD_index(i_tracklet);
        Tracklets[i_tracklet]=TRD_ST_Tracklet;
        TV3_offset      = TRD_ST_Tracklet ->get_TV3_offset();
        TV3_dir         = TRD_ST_Tracklet ->get_TV3_dir();
        i_det           = TRD_ST_Tracklet ->get_TRD_det();
		
        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        Double_t radius = TMath::Sqrt( TMath::Power(TV3_offset[0],2) + TMath::Power(TV3_offset[1],2) );

        th1d_TRD_layer_radii ->Fill(radius);
        vec_th1d_TRD_layer_radii_det[i_det] ->Fill(radius);
    }
    //--------------------------------------------------


    //TCanvas* can_TRD_layer_radii = new TCanvas("can_TRD_layer_radii","can_TRD_layer_radii",10,10,500,500);
    //can_TRD_layer_radii ->cd();
    //th1d_TRD_layer_radii ->DrawCopy("h");

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Animate_beams(Double_t beam_path)
{
     TEveP_beamA ->SetPoint(0,0,0,-beam_path);
     TEveP_beamA ->SetMarkerSize(3);
     TEveP_beamA ->SetMarkerStyle(20);
     TEveP_beamA ->SetMarkerColor(kRed);

     TEveP_beamB ->SetPoint(0,0,0,beam_path);
     TEveP_beamB ->SetMarkerSize(3);
     TEveP_beamB ->SetMarkerStyle(20);
     TEveP_beamB ->SetMarkerColor(kGreen+1);

     gEve->AddElement(TEveP_beamA);
     gEve->AddElement(TEveP_beamB);

     if(beam_path <= 0.0) gEve ->AddElement(TEveP_primary_vertex);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_TPC_track(Int_t i_track, Int_t color, Double_t line_width, Double_t max_path)
{
    Double_t track_pos[3];
    Double_t radius_helix;

    TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

#if defined(USEEVE)
    vec_TPL3D_helix.resize(i_track+1);
    vec_TPL3D_helix_hull.resize(i_track+1);
    vec_TPL3D_helix_inner.resize(i_track+1);
    vec_TPL3D_helix[i_track] = new TEveLine();
    vec_TPL3D_helix_hull[i_track] = new TEveLine();
    vec_TPL3D_helix_inner[i_track] = new TEveLine();
#endif
    //for(Int_t i_param=0;i_param<6;i_param++)
    //    cout<<TRD_ST_TPC_Track ->getHelix_param(i_param)<<" ";
    //cout<<endl;
    for(Double_t track_path = 0.0; track_path < max_path; track_path += 1.0)
    {
        TRD_ST_TPC_Track ->Evaluate(track_path,track_pos);
        radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
        if(radius_helix > 370.0) break;
        if(fabs(track_pos[2]) > 360.0) break;
        if(radius_helix > 80.0)
        {
#if defined(USEEVE)
            vec_TPL3D_helix[i_track]        ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
            vec_TPL3D_helix_hull[i_track]   ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
#endif
        }
        if(radius_helix < 80.0)
        {
#if defined(USEEVE)
            vec_TPL3D_helix_inner[i_track] ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
#endif
        }


        //if(i_track == 0) printf("track_path: %4.3f, pos: {%4.2f, %4.2f, %4.2f} \n",track_path,track_pos[0],track_pos[1],track_pos[2]);
    }

#if defined(USEEVE)
    HistName = "track ";
    HistName += i_track;
    vec_TPL3D_helix[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix[i_track]    ->SetLineWidth(line_width);
    vec_TPL3D_helix[i_track]    ->SetMainColor(color);
    vec_TPL3D_helix[i_track]    ->SetMainAlpha(1.0);

    HistName = "track (h) ";
    HistName += i_track;
    vec_TPL3D_helix_hull[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix_hull[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix_hull[i_track]    ->SetLineWidth(8);
    vec_TPL3D_helix_hull[i_track]    ->SetMainColor(kWhite);
    vec_TPL3D_helix_hull[i_track]    ->SetMainAlpha(0.3);

    HistName = "track (h) ";
    HistName += i_track;
    vec_TPL3D_helix_inner[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix_inner[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix_inner[i_track]    ->SetLineWidth(2);
    vec_TPL3D_helix_inner[i_track]    ->SetMainColor(kGray);
    vec_TPL3D_helix_inner[i_track]    ->SetMainAlpha(0.8);
    //if(i_track == 3)
    {
        gEve->AddElement(vec_TPL3D_helix[i_track]);
        gEve->AddElement(vec_TPL3D_helix_hull[i_track]);
        gEve->AddElement(vec_TPL3D_helix_inner[i_track]);
    }
#endif

    //printf("%s Ali_TRD_ST_Analyze::Draw_TPC_track, i_track: %d %s \n",KRED,i_track,KNRM,1000.0);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::Draw_event(Long64_t i_event, Int_t graphics, Int_t draw_tracks, Int_t draw_tracklets, Double_t track_path)
{
    //printf("Ali_TRD_ST_Analyze::Draw_event \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event


    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
    //--------------------------------------------------


    //--------------------------------------------------
    // TPC track loop
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

        Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
        Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
        UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
        UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
        UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
        Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
        Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
        Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();
		
        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();

        //if(momentum < 0.3 || dca > 2.0) continue;
        if(momentum < 0.3) continue;


        if(graphics && draw_tracks) Draw_TPC_track(i_track,track_color,3,track_path);
    }

    //--------------------------------------------------



    //--------------------------------------------------
    // TRD tracklet loop
    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t    i_det;

    for(Int_t i_tracklet = 0; i_tracklet < NumTracklets; i_tracklet++)
    {
        TRD_ST_Tracklet = TRD_ST_Event    ->getTracklet(i_tracklet);
        TV3_offset      = TRD_ST_Tracklet ->get_TV3_offset();
        TV3_dir         = TRD_ST_Tracklet ->get_TV3_dir();
        i_det           = TRD_ST_Tracklet ->get_TRD_det();

        if(TV3_offset.Mag() > 1000.0) continue;

        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        Double_t radius = TMath::Sqrt( TMath::Power(TV3_offset[0],2) + TMath::Power(TV3_offset[1],2) );
      

#if defined(USEEVE)
        if(graphics)
        {
            vec_TEveLine_tracklets[i_layer].resize(N_tracklets_layers[i_layer]+1);
            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] = new TEveLine();

            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0],TV3_offset[1],TV3_offset[2]);
            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0] + scale_length_vec*TV3_dir[0],TV3_offset[1] + scale_length_vec*TV3_dir[1],TV3_offset[2] + scale_length_vec*TV3_dir[2]);


            //printf("i_tracklet: %d, radius: %4.3f, pos A: {%4.2f, %4.2f, %4.2f}, pos B: {%4.2f, %4.2f, %4.2f} \n",i_tracklet,radius,TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + scale_length_vec*TV3_dir[0],TV3_offset[1] + scale_length_vec*TV3_dir[1],TV3_offset[2] + scale_length_vec*TV3_dir[2]);

            HistName = "tracklet ";
            HistName += i_tracklet;
            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetName(HistName.Data());
            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetLineStyle(1);
            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetLineWidth(3);
            vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetMainColor(color_layer[i_layer]);
            //if(i_tracklet == 63 || i_tracklet == 67 || i_tracklet == 72 || i_tracklet == 75 || i_tracklet == 83 || i_tracklet == 88)
            {
             if(graphics && draw_tracklets) gEve->AddElement(vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]);
            }
        }
#endif
        N_tracklets_layers[i_layer]++;
    }

#if defined(USEEVE)
    if(graphics) gEve->Redraw3D(kTRUE);
#endif
    //--------------------------------------------------

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::Do_TPC_TRD_matching(Long64_t i_event, Double_t xy_matching_window, Double_t z_matching_window, Int_t graphics)
{
    //printf("Ali_TRD_ST_Analyze::Do_TPC_TRD_matching \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

    Int_t offset_point = 0;
    Int_t TPC_at_offset_point = 0;
    matched_tracks.clear();
    vec_idx_matched_TPC_track.clear();

    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
    Float_t  V0MEq                = TRD_ST_Event ->getcent_class_V0MEq();
    //--------------------------------------------------



    //--------------------------------------------------
    // Store all TRD tracklets in an array to keep it in memory
    // TRD tracklet loop
    vector< vector<TVector3> > vec_TV3_dir_tracklets;
    vector< vector<TVector3> > vec_TV3_offset_tracklets;
    vector< vector<Int_t> > TRD_index_tracklets;

    vector <vector < vector<Double_t>>> vec_ADC_val;
    vec_TV3_dir_tracklets.resize(540);
    vec_TV3_offset_tracklets.resize(540);
    TRD_index_tracklets.resize(540);
    //vec_ADC_val.resize(540);

    //for (Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    vec_ADC_val[i_det].resize(0);
    //}

    //Double_t ADC_val[24];

    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t    i_det;

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_offset_points        = new TEvePointSet();
        TEveP_TPC_at_offset_points = new TEvePointSet();
    }
#endif

    for(Int_t i_tracklet = 0; i_tracklet < NumTracklets; i_tracklet++)
    {
        TRD_ST_Tracklet = TRD_ST_Event    ->getTracklet(i_tracklet);
        TV3_offset      = TRD_ST_Tracklet ->get_TV3_offset();
        TV3_dir         = TRD_ST_Tracklet ->get_TV3_dir();
        i_det           = TRD_ST_Tracklet ->get_TRD_det();

        //printf("test 1 \n");

        //vec_ADC_val[i_det].resize((Int_t)vec_ADC_val[i_det].size()+1); 

        //vec_ADC_val[i_det][i_tracklet].resize(24);

        //for (Int_t i_time = 0; i_time < 24; i_time++) { vec_ADC_val[i_det][i_tracklet][i_time]        = TRD_ST_Tracklet ->get_ADC_val(i_time); }

        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        vec_TV3_dir_tracklets[i_det].push_back(TV3_dir);
        vec_TV3_offset_tracklets[i_det].push_back(TV3_offset);
        TRD_index_tracklets[i_det].push_back(i_tracklet);
        //printf("test 2 \n");

    }
    //--------------------------------------------------


    //--------------------------------------------------
    // TPC track loop
    Double_t track_pos[3];
    Double_t radius_helix;
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        //for(Int_t i_track = 3; i_track < 4; i_track++)
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

        Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
        Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
        UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
        UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
        UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
        Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
        Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
        Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();
        Float_t charge = 1.0;
        if(dca < 0.0) charge = -1.0;

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();


        if(momentum < 0.3) continue;

        vector<TVector3> vec_TV3_helix_points_at_TRD_layers;
        vec_TV3_helix_points_at_TRD_layers.resize(6);
        matched_tracks.resize(matched_tracks.size() + 1);
        Int_t tracks_size = matched_tracks.size() - 1;  // index of current track
        vec_idx_matched_TPC_track.push_back(i_track);

        //printf("tracks_size: %d, size_TPC_idx: %d \n",tracks_size,(Int_t)vec_idx_matched_TPC_track.size());

        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            Double_t radius_layer_center = 0.5*(TRD_layer_radii[i_layer][0] + TRD_layer_radii[i_layer][1]);

            // Find the helix path which touches the first TRD layer
            Double_t track_path_add      = 10.0;
            Double_t track_path_layer0   = 0.0;
            Double_t radius_helix_layer0 = 0.0;
            for(Double_t track_path = TRD_layer_radii[i_layer][0]; track_path < 1000; track_path += track_path_add)
            {
                TRD_ST_TPC_Track ->Evaluate(track_path,track_pos);
                radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
                if(radius_helix < radius_layer_center)
                {
                    track_path_layer0   = track_path;
                    radius_helix_layer0 = radius_helix;
                    //printf("radius_helix_layer0: %4.3f, track_path_add: %4.3f \n",radius_helix_layer0,track_path_add);
                }
                else
                {
                    track_path -= track_path_add;
                    track_path_add *= 0.5;
                }
                if(track_path_add < 1.0)
                {
                    vec_TV3_helix_points_at_TRD_layers[i_layer].SetXYZ(track_pos[0],track_pos[1],track_pos[2]);
                    break;
                }
            }
            //printf("   --> i_layer: %d, track_path_layer0: %4.3f, radius_helix_layer0: %4.3f \n",i_layer,track_path_layer0,radius_helix_layer0);


            TVector3 TV3_diff_vec;
            Double_t dist_min      = 10.0;
            Double_t dist          = 0.0;
            Int_t    det_best      = -1;
            Int_t    tracklet_best = -1;

            for(Int_t i_det = 0; i_det < 540; i_det++)
            {
                Int_t i_layer_from_det  = i_det%6;
                if(i_layer_from_det != i_layer) continue;

                for(Int_t i_tracklet = 0; i_tracklet < (Int_t)vec_TV3_offset_tracklets[i_det].size(); i_tracklet++)
                {
                    TV3_diff_vec = vec_TV3_offset_tracklets[i_det][i_tracklet] - vec_TV3_helix_points_at_TRD_layers[i_layer];
                    dist = TV3_diff_vec.Mag();
                    if(dist < dist_min)
                    {
                        dist_min      = dist;
                        det_best      = i_det;
                        tracklet_best = i_tracklet;
                    }
                }
            }

            matched_tracks[tracks_size].push_back(NULL);


            if(det_best < 0 || tracklet_best < 0) continue;

            Ali_TRD_ST_Tracklets* temp_tracklet = new Ali_TRD_ST_Tracklets();
            temp_tracklet->set_TRD_det(det_best);
            temp_tracklet->set_TV3_offset(vec_TV3_offset_tracklets[det_best][tracklet_best]);
            temp_tracklet->set_TV3_dir(vec_TV3_dir_tracklets[det_best][tracklet_best]);
            temp_tracklet->set_TRD_index(TRD_index_tracklets[det_best][tracklet_best]);
            //for (Int_t i_time = 0; i_time < 24; i_time++)
            //{
            //    temp_tracklet ->set_ADC_val(i_time,vec_ADC_val[det_best][tracklet_best][i_time]); 
            //}
            matched_tracks[tracks_size][matched_tracks[tracks_size].size() - 1] = temp_tracklet;


//printf("tracklet TRD index: %d, i_track: %d, i_layer: %d \n",TRD_index_tracklets[det_best][tracklet_best],tracks_size,matched_tracks[tracks_size].size() - 1);


            //------------------------------------------
            // Calculate exact position of TPC helix at TRD tracklet offset position
            TPC_single_helix ->setHelix(TRD_ST_TPC_Track->getHelix_param(0),TRD_ST_TPC_Track->getHelix_param(1),TRD_ST_TPC_Track->getHelix_param(2),TRD_ST_TPC_Track->getHelix_param(3),TRD_ST_TPC_Track->getHelix_param(4),TRD_ST_TPC_Track->getHelix_param(5));
            Float_t pathA_dca, dcaAB_dca;
            fHelixAtoPointdca(vec_TV3_offset_tracklets[det_best][tracklet_best],TPC_single_helix,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            Double_t track_pos_at_tracklet[3];
            TRD_ST_TPC_Track ->Evaluate(pathA_dca,track_pos_at_tracklet);
            //printf("%s i_layer: %d %s, i_track: %d, dcaAB_dca: %4.3f, offset_point: {%4.3f, %4.3f, %4.3f}, TPC point: {%4.3f, %4.3f, %4.3f} \n",KMAG,i_layer,KNRM,i_track,dcaAB_dca,vec_TV3_offset_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2],track_pos_at_tracklet[0],track_pos_at_tracklet[1],track_pos_at_tracklet[2]);
#if defined(USEEVE)
            if(graphics) TEveP_TPC_at_offset_points  ->SetPoint(TPC_at_offset_point,track_pos_at_tracklet[0],track_pos_at_tracklet[1],track_pos_at_tracklet[2]);
#endif
            TPC_at_offset_point++;

            Double_t residuals_xyz;
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                residuals_xyz = track_pos_at_tracklet[i_xyz] - vec_TV3_offset_tracklets[det_best][tracklet_best][i_xyz];
                vec_h2D_pT_vs_TPC_TRD_residuals[det_best][i_xyz] ->Fill(residuals_xyz,charge*pT_track);
            }
            //------------------------------------------

#if defined(USEEVE)
            Int_t size_tracklet = 0;
            if(graphics)  size_tracklet = (Int_t)vec_TEveLine_tracklets_match[i_layer].size();
#endif

            //printf("i_layer: %d, size_tracklet: %d \n",i_layer,size_tracklet);

#if defined(USEEVE)


            if(graphics)
            {
                vec_TEveLine_tracklets_match[i_layer].resize(size_tracklet+1);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet] = new TEveLine();
                vec_TEveLine_tracklets_match[i_layer][size_tracklet] ->SetNextPoint(vec_TV3_offset_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2]);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet] ->SetNextPoint(vec_TV3_offset_tracklets[det_best][tracklet_best][0] + scale_length_vec*vec_TV3_dir_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1] + scale_length_vec*vec_TV3_dir_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2] + scale_length_vec*vec_TV3_dir_tracklets[det_best][tracklet_best][2]);

                HistName = "tracklet (m) ";
                HistName += size_tracklet;
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetName(HistName.Data());
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetLineStyle(1);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetLineWidth(6);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetMainColor(color_layer_match[i_layer]);

                TEveP_offset_points  ->SetPoint(offset_point,vec_TV3_offset_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2]);
                TEveP_offset_points  ->SetMarkerSize(2);
                TEveP_offset_points  ->SetMarkerStyle(20);
                TEveP_offset_points  ->SetMarkerColor(kCyan+1);
                offset_point++;
                //if(i_tracklet == 63 || i_tracklet == 67 || i_tracklet == 72 || i_tracklet == 75 || i_tracklet == 83 || i_tracklet == 88)

                {
                    gEve->AddElement(vec_TEveLine_tracklets_match[i_layer][size_tracklet]);
                }
            }
#endif
        }


        //--------------------------------------------------
        Int_t N_matched_tracklets = 0;
        for(Int_t i_layer = 0; i_layer < (Int_t)matched_tracks[tracks_size].size(); i_layer++)
        {
            if(matched_tracks[tracks_size][i_layer] != NULL) N_matched_tracklets++;
        }
        //printf("%s TPC track: %d %s, N_matched_tracklets: %d \n",KYEL,vec_idx_matched_TPC_track[tracks_size],KNRM,N_matched_tracklets);
        //if(N_matched_tracklets >= 4 && vec_idx_matched_TPC_track[tracks_size] == 25) Draw_TPC_track(vec_idx_matched_TPC_track[tracks_size],kAzure-2,2,1000.0);
        //printf(" \n");
        //--------------------------------------------------


    } // end of track loop
    //--------------------------------------------------

#if defined(USEEVE)
    if(graphics)
    {
        gEve->AddElement(TEveP_offset_points);

        TEveP_TPC_at_offset_points  ->SetMarkerSize(3);
        TEveP_TPC_at_offset_points  ->SetMarkerStyle(20);
        TEveP_TPC_at_offset_points  ->SetMarkerColor(kYellow);
        //gEve->AddElement(TEveP_TPC_at_offset_points);
        gEve->Redraw3D(kTRUE);
    }
#endif
    //printf("Ali_TRD_ST_Analyze::Do_TPC_TRD_matching, N TPC Tracks: %d, matched_tracks: %d \n",NumTracks,(Int_t)matched_tracks.size());

    return 1;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
TProfile* Ali_TRD_ST_Analyze::Calibrate_gain(Long64_t i_event, Int_t Bethe_flag)
{

    if (!input_SE->GetEntry( i_event )) return 0;
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks();

    Int_t N_TPC_tracks = (Int_t)matched_tracks.size();

    /*printf("det 0 vec_TV3_TRD_center[0]: {%4.3f,%4.3f,%4.3f}, vec_TV3_TRD_center[1]: {%4.3f,%4.3f,%4.3f},vec_TV3_TRD_center[2]: {%4.3f,%4.3f,%4.3f} \n",
        (*vec_TV3_TRD_center[0][0]).X(),(*vec_TV3_TRD_center[0][0]).Y(),(*vec_TV3_TRD_center[0][0]).Z(),
        (*vec_TV3_TRD_center[0][1]).X(),(*vec_TV3_TRD_center[0][1]).Y(),(*vec_TV3_TRD_center[0][1]).Z(),
        (*vec_TV3_TRD_center[0][2]).X(),(*vec_TV3_TRD_center[0][2]).Y(),(*vec_TV3_TRD_center[0][2]).Z());
*/


    ADC_rel.resize(540);  //relative ADC for each chamber

    for (Int_t i_TPC_track = 0; i_TPC_track < N_TPC_tracks; i_TPC_track++)
    //TVector3 TPC_TV3_offset = matched_tracks[i_TPC_track][i_layer_TPC]->get_TV3_offset();
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(vec_idx_matched_TPC_track[i_TPC_track]);

        Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
        Double_t nsigma_TPC_K   = TRD_ST_TPC_Track ->getnsigma_K_TPC();

        //if (i_TPC_track%10 == 0) 
            //printf("\ni_track: %d, nsigma_TPC_e: %4.3f, nsigma_TPC_pi: %4.3f, nsigma_TPC_p: %4.3f \n, nsigma_TPC_K: %4.3f \n",i_TPC_track,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p,nsigma_TPC_K);
        
        Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
        Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
        UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
        UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
        UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
        Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
        Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
        Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();
        Float_t charge = 1.0;
        if(dca < 0.0) charge = -1.0;

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();

        if(momentum < 0.3 && momentum > 0.7) continue;
        //if(momentum < 0.3 ) continue;
        //if (!(matched_tracks[i_TPC_track])==NULL) continue;

        if ((abs(nsigma_TPC_pi) < 2.5 && abs(nsigma_TPC_p) > 2.5 && abs(nsigma_TPC_K) > 2.5) // Pi
                || (abs(nsigma_TPC_pi) > 2.5 && abs(nsigma_TPC_p) < 2.5 && abs(nsigma_TPC_K) > 2.5) //p
                    || (abs(nsigma_TPC_pi) > 2.5 && abs(nsigma_TPC_p) > 2.5 && abs(nsigma_TPC_K) < 2.5)) //K
        {
            for (Int_t i_layer = 0; i_layer < (Int_t)matched_tracks[i_TPC_track].size(); i_layer++)
            {
//printf("i_TPC_track: %d, i_layer: %d, (Int_t)matched_tracks[i_TPC_track].size(): %d \n",i_TPC_track,i_layer,(Int_t)matched_tracks[i_TPC_track].size());

                if(matched_tracks[i_TPC_track][i_layer] != NULL)
                {
                    Int_t det = matched_tracks[i_TPC_track][i_layer] ->get_TRD_det();
                    //ADC_rel[det].push_back(NULL);

                    Double_t ADC = 0.0;
                    Int_t tb_count = 0;

                    Ali_TRD_ST_Tracklets* TRD_ST_Tracklet_adc = TRD_ST_Event    ->getTracklet((Int_t)matched_tracks[i_TPC_track][i_layer]->get_TRD_index());
                    for (Int_t i_time = 0; i_time < 24; i_time++) //24 time bins!
                    {

                        //TRD_index_tracklets[i_det][i_trkl] = index
                        if (TRD_ST_Tracklet_adc ->get_ADC_val(i_time) < 10.0) continue;
                        ADC += (TRD_ST_Tracklet_adc ->get_ADC_val(i_time) - 10.0);
                        tb_count++;
                    }

                    if (tb_count == 23) 
                    {

                        TVector3 TV3_offset_matched = matched_tracks[i_TPC_track][i_layer] ->get_TV3_offset();
                        TVector3 TV3_dir_matched    = matched_tracks[i_TPC_track][i_layer] ->get_TV3_dir();

                        if (TV3_dir_matched.Dot((*vec_TV3_TRD_center[det][2])) == 0.0) continue;

                        Double_t scale_factor = 3.7*((*vec_TV3_TRD_center[det][2]).Mag())/((*vec_TV3_TRD_center[det][2]).Dot(TV3_dir_matched));
                        // 3.7 divided by (dot product of vectors/length of norm vector)

                        //Double_t tracklet_length = 3.7/TMath::Sin(((-1*(*vec_TV3_TRD_center[det][2])).Angle(TV3_dir_matched)));

                        //TMath::Sqrt( TMath::Power(3.0*TV3_dir_matched[0],2) + TMath::Power(3.0*TV3_dir_matched[1],2))
                        Double_t tracklet_length = TMath::Sqrt( TMath::Power(scale_factor*TV3_dir_matched[0],2) + TMath::Power(scale_factor*TV3_dir_matched[1],2) + TMath::Power(scale_factor*TV3_dir_matched[2],2));
                        //same as scale factor 

                        //vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0],TV3_offset[1],TV3_offset[2]);
                        //vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0] + scale_length_vec*TV3_dir[0],TV3_offset[1] + scale_length_vec*TV3_dir[1],TV3_offset[2] + scale_length_vec*TV3_dir[2]);
                        //if (tracklet_length > 10.0) continue;
        //printf("\ntracklet TV3_dir_matched[0]: %4.3f, TV3_dir_matched[1]: %4.3f, TV3_dir_matched[2]: %4.3f, dot: %4.3f, tracklet_length: %4.3f \n",TV3_dir_matched[0],TV3_dir_matched[1],TV3_dir_matched[2],TV3_dir_matched.Dot((*vec_TV3_TRD_center[det][2])),tracklet_length);
        //printf("normZ[0]: %4.3f, normZ[1]: %4.3f, normZ[2]: %4.3f \n",(*vec_TV3_TRD_center[det][2]).X(),(*vec_TV3_TRD_center[det][2]).Y(),(*vec_TV3_TRD_center[det][2]).Z());
        //printf("x2: %4.3f, y2: %4.3f,  z2: %4.3f, sum: %4.3f, sqrt: %4.3f \n",TMath::Power(3.0*TV3_dir_matched[0],2),TMath::Power(3.0*TV3_dir_matched[1],2),TMath::Power(3.0*TV3_dir_matched[2],2),
            //TMath::Power(3.0*TV3_dir_matched[0],2) + TMath::Power(3.0*TV3_dir_matched[1],2) + TMath::Power(3.0*TV3_dir_matched[2],2),TMath::Sqrt( TMath::Power(3.0*TV3_dir_matched[0],2) + TMath::Power(3.0*TV3_dir_matched[1],2) + TMath::Power(3.0*TV3_dir_matched[2],2)));

                        //---------------calculate local position here ----------------------

//vec_TV3_TRD_center
                        TVector3 local_little_vec = (*vec_TV3_TRD_offset[det]) - TV3_offset_matched; // pointing from center of the chamber to the datapoint?
                        Double_t x_local = local_little_vec.Dot((*vec_TV3_TRD_center[det][0]))/(*vec_TV3_TRD_center[det][0]).Mag();
                        Double_t z_local = local_little_vec.Dot((*vec_TV3_TRD_center[det][1]))/(*vec_TV3_TRD_center[det][1]).Mag();
                        /*
                        printf("i_det: %d, vec_TV3_TRD_offset: {%4.3f, %4.3f, %4.3f}, TV3_offset_matched: {%4.3f, %4.3f, %4.3f}, little vector: {%4.3f, %4.3f, %4.3f}, local coordinate (x,z): {%4.3f, %4.3f} \n",
                            det,(*vec_TV3_TRD_offset[det]).X(),(*vec_TV3_TRD_offset[det]).Y(),(*vec_TV3_TRD_offset[det]).Z(),
                            TV3_offset_matched.X(),TV3_offset_matched.Y(),TV3_offset_matched.Z(),
                            local_little_vec.X(),local_little_vec.Y(),local_little_vec.Z(),x_local,z_local);
                        */
                        //-------------------------------------------------------------------
                        if(Bethe_flag == 1) //do via Bethe-Bloch
                        {
                            //vector<TGraph*> tg_bethe_bloch; // [e mu pi K p]

                            Double_t BB_val = -1.0; 

                            if (abs(nsigma_TPC_pi) < 2.5 ) BB_val = tg_bethe_bloch[2] ->Eval(momentum); //pi
                            if (abs(nsigma_TPC_p) < 2.5 )  BB_val = tg_bethe_bloch[4] ->Eval(momentum); //p
                            if (abs(nsigma_TPC_K) < 2.5 )  BB_val = tg_bethe_bloch[3] ->Eval(momentum); //p

                            if (BB_val == 0 || BB_val == -1.0) continue;

                            ADC = ADC/(BB_val*abs(tracklet_length));
                        }
                        if(Bethe_flag == 0) // do via TPC dEdx
                        {
                            if(TPCdEdx == 0) continue;
                            ADC = ADC/(TPCdEdx*abs(tracklet_length));
                        }

                        ADC_rel[det].push_back(ADC); 
                        tp_gain->Fill(det,ADC);
                        tp_gain_uncor->Fill(det,ADC*abs(tracklet_length));

                        th2d_gain->Fill(det,ADC);

                        vec_tp2d_gain_vs_xz[det]->Fill(x_local,z_local,ADC);

                        //printf("\ndet: %d, TP content: %4.3f, ADC: %4.3f, TPCdEdx: %4.3f, tracklet_length: %4.3f \n",det,tp_gain->GetBinContent(det),ADC,TPCdEdx,tracklet_length);
                    }
                }
            }
        }
    }

#if 0
    if((Int_t)vec_tp_gain.size() == 0)
    {
        vec_tp_gain.resize(540);
        vec_TH2D_gain.resize(540);

        for (Int_t i_det = 0; i_det < 540; i_det++)
        {
            vec_tp_gain[i_det] = new TProfile(Form("vec_tp_gain_%d",i_det),Form("vec_tp_gain_%d",i_det),360,-360,360);
            vec_TH2D_gain[i_det] = new TH2D(Form("vec_tp_gain_%d",i_det),Form("vec_tp_gain_%d",i_det),10,70,110,50,-25,25);
        }
    }
#endif

//tp_gain ->Draw();
//th2d_gain->Draw();  

return tp_gain;

}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_Save_Gain_calib()
{
    #if 0

    TH1D* h_dummy_TP_gain = new TH1D("h_dummy_TP_gain","h_dummy_TP_gain",90,50,140);
    h_dummy_TP_gain->SetStats(1);
    h_dummy_TP_gain->SetTitle("");
    h_dummy_TP_gain->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_TP_gain->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_TP_gain->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_TP_gain->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_TP_gain->GetXaxis()->SetLabelSize(0.04);
    h_dummy_TP_gain->GetYaxis()->SetLabelSize(0.04);
    h_dummy_TP_gain->GetXaxis()->SetTitleSize(0.04);
    h_dummy_TP_gain->GetYaxis()->SetTitleSize(0.04);
    //h_dummy_TP_gain->GetXaxis()->SetNdivisions(540,'N');
    //h_dummy_TP_gain->GetYaxis()->SetNdivisions(5,'N');
    h_dummy_TP_gain->GetXaxis()->CenterTitle();
    h_dummy_TP_gain->GetYaxis()->CenterTitle();
    h_dummy_TP_gain->GetXaxis()->SetTitle("Chamber ID");
    h_dummy_TP_gain->GetYaxis()->SetTitle("<ADC/dEdx>");
    h_dummy_TP_gain->GetXaxis()->SetRangeUser(0,539);
    //h_dummy_TP_gain->GetYaxis()->SetRangeUser(-24,24);
    #endif

    TProfile* tp_gain_avg_sec = new TProfile("tp_gain_avg_sec","tp_gain_avg_sec",540,0,539);
    TProfile* tp_gain_avg_sec_uncor = new TProfile("tp_gain_avg_sec_uncor","tp_gain_avg_sec_uncor",540,0,539);

    for(Int_t bin = 1; bin <= 30; bin++)
    {
        for (Int_t i_sect = 0; i_sect < 18; i_sect++)
        {
            Double_t bin_cont = tp_gain ->GetBinContent(bin+30*i_sect);
            if (bin_cont > 2.8 && bin_cont < 4.2) tp_gain_avg_sec ->Fill(bin,bin_cont);

            Double_t bin_cont_uncor = tp_gain_uncor ->GetBinContent(bin+30*i_sect);
            if (bin_cont_uncor > 7.0 && bin_cont_uncor < 9.0) tp_gain_avg_sec_uncor ->Fill(bin,bin_cont_uncor);
        }
    }

    for (Int_t bin = 1; bin <= 30; bin++)
    {
        for (Int_t i_sect = 1; i_sect < 18; i_sect++)
        {
            Double_t bin_cont = tp_gain_avg_sec ->GetBinContent(bin);
            tp_gain_avg_sec ->Fill(bin+30*i_sect,bin_cont);

            Double_t bin_cont_uncor = tp_gain_avg_sec_uncor ->GetBinContent(bin);
            tp_gain_avg_sec_uncor ->Fill(bin+30*i_sect,bin_cont_uncor);

        }
    }

    //---------------------------------------

    TFile* inputfile_dEdx_length = TFile::Open("./TRD_dEdx_length.root");
    tp_dEdx_length_det = (TProfile*)inputfile_dEdx_length ->Get("tp_dEdx_length_det");
    h_dEdx_length_det_orig = (TH1D*)tp_dEdx_length_det ->ProjectionX("h_dEdx_length_det_orig");
    h_dEdx_length_det = (TH1D*)tp_dEdx_length_det ->ProjectionX("h_dEdx_length_det");

#if 0 //draw gain hists

    for(Int_t bin = 1; bin <= h_dEdx_length_det->GetNbinsX(); bin++)
    {
        Double_t bin_cont = h_dEdx_length_det ->GetBinContent(bin);
        Double_t bin_cont_mod = bin_cont*0.375+2.23;
        h_dEdx_length_det ->SetBinContent(bin,bin_cont_mod);
    }

    //----------------------------------------
    
    TCanvas *can_TP_gain = new TCanvas("can_TP_gain", "can_TP_gain",10,10,1500,500);
    can_TP_gain->cd();
    tp_gain->GetXaxis()->SetTitle("Chamber ID");
    tp_gain->GetYaxis()->SetTitle("<ADC/dEdx>");
    tp_gain->GetYaxis()->SetRangeUser(1.5,5.5);
    tp_gain->SetStats(0);
    tp_gain->SetTitle("");
    //h_dummy_TP_gain->Draw("h");
    //tp_gain ->LabelsOption("");
    tp_gain ->DrawCopy("");

    h_dEdx_length_det ->SetMarkerColor(kGreen+2);
    h_dEdx_length_det ->SetMarkerSize(0.6);
    h_dEdx_length_det ->SetMarkerStyle(20);
    h_dEdx_length_det ->DrawCopy("same P");

    tp_gain_avg_sec ->SetMarkerColor(kRed+2);
    tp_gain_avg_sec ->SetMarkerSize(0.4);
    tp_gain_avg_sec ->SetMarkerStyle(20);
    tp_gain_avg_sec ->DrawCopy("same P");

    //----------------------------------------

    for(Int_t bin = 1; bin <= h_dEdx_length_det->GetNbinsX(); bin++)
    {
        Double_t bin_cont = h_dEdx_length_det_orig ->GetBinContent(bin);
        Double_t bin_cont_mod = (bin_cont)*0.75 + 5.425;
        h_dEdx_length_det ->SetBinContent(bin,bin_cont_mod);
    }

    TCanvas *can_TP_gain_uncor = new TCanvas("can_TP_gain_uncor", "can_TP_gain_uncor",10,10,1500,500);
    can_TP_gain_uncor->cd();
    tp_gain_uncor->GetXaxis()->SetTitle("Chamber ID");
    tp_gain_uncor->GetYaxis()->SetTitle("<ADC/dEdx>");
    tp_gain_uncor->GetYaxis()->SetRangeUser(6.0,10.0);
    tp_gain_uncor->SetStats(0);
    tp_gain_uncor->SetTitle("");
    //h_dummy_TP_gain->Draw("h");
    //tp_gain ->LabelsOption("");
    tp_gain_uncor ->DrawCopy("");

    h_dEdx_length_det ->SetMarkerColor(kGreen+2);
    h_dEdx_length_det ->SetMarkerSize(0.6);
    h_dEdx_length_det ->SetMarkerStyle(20);
    h_dEdx_length_det ->DrawCopy("same P");


    tp_gain_avg_sec_uncor ->SetMarkerColor(kRed+2);
    tp_gain_avg_sec_uncor ->SetMarkerSize(0.4);
    tp_gain_avg_sec_uncor ->SetMarkerStyle(20);
    tp_gain_avg_sec_uncor ->DrawCopy("same P");
    //----------------------------------------

    TCanvas *can_th2d_gain = new TCanvas("can_th2d_gain", "can_th2d_gain",10,10,1500,500);
    can_th2d_gain->cd();

    th2d_gain->GetXaxis()->SetTitle("Chamber ID");
    th2d_gain->GetYaxis()->SetTitle("<ADC/dEdx>");
    th2d_gain->SetStats(0);
    th2d_gain->SetTitle("");
    th2d_gain ->DrawCopy("COLZ");

    vector <TCanvas*> can_vec_tp2d_gain_vs_xz;
    can_vec_tp2d_gain_vs_xz.resize(540);

    for (Int_t i_det = 0; i_det < 3; i_det++)
    {
        can_vec_tp2d_gain_vs_xz[i_det] = new TCanvas(Form("can_tp2d_gain_%d",i_det), Form("can_tp2d_gain_%d",i_det),10,10,1500,500);
        can_vec_tp2d_gain_vs_xz[i_det] ->cd();
        vec_tp2d_gain_vs_xz[i_det]->GetXaxis()->SetTitle("local x, cm");
        vec_tp2d_gain_vs_xz[i_det]->GetYaxis()->SetTitle("local z, cm");
        vec_tp2d_gain_vs_xz[i_det]->SetStats(0);
        vec_tp2d_gain_vs_xz[i_det]->SetTitle("");
        vec_tp2d_gain_vs_xz[i_det] ->DrawCopy("COLZ");

    }
    #endif
    
    out_gain->cd();
    tp_gain->Write();
    th2d_gain->Write();

    for (Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp2d_gain_vs_xz[i_det] ->Write();
    }

}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_Kalman_Tracklets(vector< vector<Ali_TRD_ST_Tracklets*> > found_tracks)
{
#if defined(USEEVE)
    vector< vector<TEveLine*> > TEveLine_Kalman_found;
    TEveLine_Kalman_found.resize((Int_t)found_tracks.size());
#endif
    for(Int_t i_Track=0;i_Track<(Int_t)found_tracks.size();i_Track++)
    {
#if defined(USEEVE)
        TEveLine_Kalman_found[i_Track].resize(6);
#endif

        for(Int_t i_lay=0;i_lay<6;i_lay++)
        {
            if(found_tracks[i_Track][i_lay]!=NULL)
            {
#if defined(USEEVE)
                TEveLine_Kalman_found[i_Track][i_lay]=new TEveLine();
                TEveLine_Kalman_found[i_Track][i_lay]->SetNextPoint(found_tracks[i_Track][i_lay]->get_TV3_offset()[0],found_tracks[i_Track][i_lay]->get_TV3_offset()[1],found_tracks[i_Track][i_lay]->get_TV3_offset()[2]);
                TEveLine_Kalman_found[i_Track][i_lay]->SetNextPoint(found_tracks[i_Track][i_lay]->get_TV3_offset()[0] + scale_length_vec*found_tracks[i_Track][i_lay]->get_TV3_dir()[0],found_tracks[i_Track][i_lay]->get_TV3_offset()[1] + scale_length_vec*found_tracks[i_Track][i_lay]->get_TV3_dir()[1],found_tracks[i_Track][i_lay]->get_TV3_offset()[2] + scale_length_vec*found_tracks[i_Track][i_lay]->get_TV3_dir()[2]);
#endif

                Int_t i_det           = found_tracks[i_Track][i_lay] ->get_TRD_det();
                Int_t i_sector = (Int_t)(i_det/30);
                Int_t i_stack  = (Int_t)(i_det%30/6);
                Int_t i_layer  = i_det%6;


#if defined(USEEVE)
                HistName = "tracklet (kal) ";
                HistName += i_Track;
                HistName += "_";
                HistName += i_sector;
                HistName += "_";
                HistName += i_stack;
                HistName += "_";
                HistName += i_lay;

                TEveLine_Kalman_found[i_Track][i_lay]	->SetName(HistName.Data());
                TEveLine_Kalman_found[i_Track][i_lay]   ->SetLineStyle(1);
                TEveLine_Kalman_found[i_Track][i_lay]  	->SetLineWidth(6);
                TEveLine_Kalman_found[i_Track][i_lay]   ->SetMainColor(kOrange);
                gEve->AddElement(TEveLine_Kalman_found[i_Track][i_lay]);
#endif

            }
        }
    }
#if defined(USEEVE)
    gEve->Redraw3D(kTRUE);
#endif
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_matched_Kalman_Tracklets(Int_t i_track_plot)
{
#if defined(USEEVE)
    vector< vector<TEveLine*> > TEveLine_Kalman_found;
    TEveLine_Kalman_found.resize((Int_t)vec_kalman_TRD_trackets.size());
#endif
    for(Int_t i_Track=i_track_plot;i_Track< (i_track_plot+1);i_Track++)
    {
#if defined(USEEVE)
        TEveLine_Kalman_found[i_Track].resize(6);
#endif

        for(Int_t i_lay=0;i_lay<6;i_lay++)
        {
            if(vec_kalman_TRD_trackets[i_Track][i_lay]!=NULL)
            {
#if defined(USEEVE)
                TEveLine_Kalman_found[i_Track][i_lay]=new TEveLine();
                TEveLine_Kalman_found[i_Track][i_lay]->SetNextPoint(vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_offset()[0],vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_offset()[1],vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_offset()[2]);
                TEveLine_Kalman_found[i_Track][i_lay]->SetNextPoint(vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_offset()[0] + scale_length_vec*vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_dir()[0],vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_offset()[1] + scale_length_vec*vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_dir()[1],vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_offset()[2] + scale_length_vec*vec_kalman_TRD_trackets[i_Track][i_lay]->get_TV3_dir()[2]);
#endif

                Int_t i_det           = vec_kalman_TRD_trackets[i_Track][i_lay] ->get_TRD_det();
                Int_t i_sector = (Int_t)(i_det/30);
                Int_t i_stack  = (Int_t)(i_det%30/6);
                Int_t i_layer  = i_det%6;


#if defined(USEEVE)
                HistName = "tracklet (kal) ";
                HistName += i_Track;
                HistName += "_";
                HistName += i_sector;
                HistName += "_";
                HistName += i_stack;
                HistName += "_";
                HistName += i_lay;

                TEveLine_Kalman_found[i_Track][i_lay]	->SetName(HistName.Data());
                TEveLine_Kalman_found[i_Track][i_lay]   ->SetLineStyle(1);
                TEveLine_Kalman_found[i_Track][i_lay]  	->SetLineWidth(6);
                TEveLine_Kalman_found[i_Track][i_lay]   ->SetMainColor(kOrange);
                gEve->AddElement(TEveLine_Kalman_found[i_Track][i_lay]);
#endif

            }
        }
    }
#if defined(USEEVE)
    gEve->Redraw3D(kTRUE);
#endif
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Plot_AP()
{
    TCanvas *can_AP_plot = new TCanvas("can_AP_plot", "can_AP_plot",10,10,500,500);
    can_AP_plot->cd();
    TH2D_AP_plot ->DrawCopy("colz");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Plot_pT_TPC_vs_Kalman()
{
    TCanvas *can_pT_TPC_vs_Kalman = new TCanvas("can_pT_TPC_vs_Kalman", "can_pT_TPC_vs_Kalman",10,10,500,500);
    can_pT_TPC_vs_Kalman->cd();
    TH2D_pT_TPC_vs_Kalman ->DrawCopy("colz");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_hist_TPC_tracklet_diffs()
{
    printf("Ali_TRD_ST_Analyze::Draw_hist_TPC_tracklet_diffs \n");

    th1d_offset_diff->GetXaxis()->SetTitle("Offset Difference (cm)");
    th1d_offset_diff->GetYaxis()->SetTitle("Number of Tracklets");
    th1d_offset_diff->GetXaxis()->CenterTitle();
    th1d_offset_diff->GetYaxis()->CenterTitle();

    TCanvas *th1d_offset_diff_can = new TCanvas("th1d_offset_diff_can", "Max");
    th1d_offset_diff_can->cd();
    th1d_offset_diff->Draw();

    th1d_angle_diff->GetXaxis()->SetTitle("Angle Difference (deg)");
    th1d_angle_diff->GetYaxis()->SetTitle("Number of Tracklets");
    th1d_angle_diff->GetXaxis()->CenterTitle();
    th1d_angle_diff->GetYaxis()->CenterTitle();

    TCanvas *th1d_angle_diff_can = new TCanvas("th1d_angle_diff_can", "Plateu");
    th1d_angle_diff_can->cd();
    th1d_angle_diff->Draw();
}
//----------------------------------------------------------------------------------------

TH1I* Ali_TRD_ST_Analyze::get_h_good_bad_TRD_chambers(){
	return h_good_bad_TRD_chambers;
}


//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Write()
{
    printf("Write data to file \n");
    outputfile ->cd();
    outputfile ->mkdir("ST_Physics");
    outputfile ->cd("ST_Physics");

    Tree_TRD_Self_Event_out -> Write();

    outputfile ->cd();


    //NT_secondary_vertices         ->AutoSave("SaveSelf");
    //NT_secondary_vertex_cluster   ->AutoSave("SaveSelf");
    TH2D_AP_plot          ->Write();
    outputfile ->mkdir("AP_radii");
    outputfile ->cd("AP_radii");
    for(Int_t i_radius = 0; i_radius < N_AP_radii; i_radius++)
    {
        vec_TH2D_AP_plot_radius[i_radius] ->Write();
    }
    outputfile ->cd();
    TH2D_pT_TPC_vs_Kalman ->Write();
    for(Int_t i_pt_res = 0; i_pt_res < N_pT_resolution; i_pt_res++)
    {
        vec_TH2D_pT_TPC_vs_Kalman[i_pt_res] ->Write();
        vec_TH2D_one_over_pT_TPC_vs_Kalman[i_pt_res] ->Write();
    }

    outputfile ->mkdir("layer_radii");
    outputfile ->cd("layer_radii");
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        vec_th1d_TRD_layer_radii_det[TRD_detector] ->Write();
    }
    outputfile ->cd();

    outputfile ->mkdir("TPC_TRD_residuals");
    outputfile ->cd("TPC_TRD_residuals");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_h2D_pT_vs_TPC_TRD_residuals[i_det][i_xyz] ->Write();
        }
    }
    outputfile ->cd();

    tp_efficiency_matching_vs_pT ->Write();
    tp_efficiency_all_vs_pT      ->Write();
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_h2D_delta_pT_all_vs_pT[i_layer] ->Write();
    }

    printf("All data written \n");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Calibrate(Int_t graphics)
{
    //printf("Ali_TRD_ST_Analyze::Calibrate() \n");
    if((Int_t)vec_tp_Delta_vs_impact.size() == 0 && (Int_t)vec_tp_Delta_vs_impact_circle.size() == 0)
    {
        vec_tp_Delta_vs_impact.resize(540);
        vec_TH2D_Delta_vs_impact.resize(540);
        vec_tp_Delta_vs_impact_circle.resize(540);
        vec_TH2D_Delta_vs_impact_circle.resize(540);

        for (Int_t i_det = 0; i_det < 540; i_det++)
        {
            vec_tp_Delta_vs_impact[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_%d",i_det),Form("vec_th1d_Delta_vs_impact_%d",i_det),360,-360,360);
            vec_TH2D_Delta_vs_impact[i_det] = new TH2D(Form("vec_th2d_Delta_vs_impact_%d",i_det),Form("vec_th2d_Delta_vs_impact_%d",i_det),10,70,110,50,-25,25);
            vec_tp_Delta_vs_impact_circle[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),360,-360,360);
            vec_TH2D_Delta_vs_impact_circle[i_det] = new TH2D(Form("vec_th2d_Delta_vs_impact_circle_%d",i_det),Form("vec_th2d_Delta_vs_impact_circle_%d",i_det),13,70,110,200,-100,100);


        }
    }
    TVector3 vec_TV3_tracklet_vectors, TV3_tracklet_off_vector, vec_dir_vec_circle;
    vector< TVector3 > vec_dir_vec_circle_circle;
    Float_t pathA_dca, dcaAB_dca;
    Double_t track_posA[3], track_posB[3];
    Int_t detector[6] = {0};
    Int_t counter = 0;
    Int_t counter_circ = 0;

    vector< vector<TEveLine*> > TEveLine_vec_dir_vec_circle;
    vector< vector<TEveLine*> > TEveLine_vec_dir_vec_circle_circle;
    vector<TEveLine*> TEveLine_circle;

#if defined(USEEVE)
    if(graphics)
    {
        TEveLine_vec_dir_vec_circle.resize((Int_t)mHelices_kalman.size());
        TEveLine_vec_dir_vec_circle_circle.resize((Int_t)mHelices_kalman.size());
        TEveLine_circle.resize((Int_t)mHelices_kalman.size());
    }
#endif

    //printf("RunID: %d \n",Global_RunID);

    for(Int_t i_track = 0; i_track < (Int_t)mHelices_kalman.size(); i_track++)
    //for(Int_t i_track = 0; i_track < 1; i_track++)
    {
        //printf("i_track: %d \n",i_track);
        Double_t arr_path_layer[6] = {-999.0};
        Double_t arr_dca_layer[6]  = {-999.0};
        Double_t N_good_tracklets = 0.0;
        Double_t average_dca = 0.0;

        set_single_helix_params(mHelices_kalman[i_track]);
        Double_t pT_KF = mHelices_kalman[i_track][6]; // pT

        if(pT_KF < 1.2) continue;



#if 1 //----stuff for circle fit-------------------
        vector< TVector2 > vec_TV2_points;
        vec_TV2_points.resize(6);
        vec_Dt_digit_pos_cluster.resize(6);
        for (Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            vec_TV2_points[i_layer].SetX(-999.9);
            vec_TV2_points[i_layer].SetY(-999.9);
            vec_Dt_digit_pos_cluster[i_layer].resize(2);
            for (Int_t i_xy = 0; i_xy < 2; i_xy++)
            {
                vec_Dt_digit_pos_cluster[i_layer][i_xy] = -999.9;
            }
        }
        //-------------------------------------------------
#endif


        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            //printf("i_layer: %d \n",i_layer);

            if(vec_kalman_TRD_trackets[i_track][i_layer] == NULL ) continue;

            TV3_tracklet_off_vector  	 = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset();
            detector[i_layer]		 = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TRD_det();

            fHelixAtoPointdca(TV3_tracklet_off_vector,vec_helices_TRD[i_track],pathA_dca,dcaAB_dca); // new helix to point dca calculation
            arr_path_layer[i_layer] = pathA_dca;
            arr_dca_layer[i_layer]  = dcaAB_dca;
            //printf("   ---> pos trkl: {%4.3f, %4.3f, %4.3f}, pathA_dca: %4.3f, dcaAB_dca: %4.3f \n",TV3_tracklet_off_vector.X(),TV3_tracklet_off_vector.Y(),TV3_tracklet_off_vector.Z(),pathA_dca,dcaAB_dca);
            average_dca += dcaAB_dca;
            N_good_tracklets += 1.0;
        }

        if(N_good_tracklets > 0.0) average_dca /= N_good_tracklets;
        if(average_dca > 3.0) continue;

        Float_t dist_to_zero = primary_vertex_dca(i_track); //select only tracks that started at primary vertex
        //if(dist_to_zero > 3.0) continue;

        counter++;
        //printf("counter = %d \n",counter);

#if defined(USEEVE) //Draw tracks used in calibration
        if(graphics)
        {
            Draw_Kalman_Helix_Tracks(i_track,kGreen,3.0,500.0);
        }
#endif

#if defined(USEEVE)
        if(graphics)
        {
            TEveLine_vec_dir_vec_circle[i_track].resize(6);
        }
#endif

        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            //printf("i_layer: %d, i_det: %d, dca: %4.3f \n",i_layer,detector[i_layer],arr_dca_layer[i_layer]);

            if(vec_kalman_TRD_trackets[i_track][i_layer] == NULL ) continue;

            vec_TV3_tracklet_vectors 	 = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_dir();
            vec_TV3_tracklet_vectors[2]  = 0.0;
            vec_TV3_tracklet_vectors    *= 1/ vec_TV3_tracklet_vectors.Mag();

            Evaluate(arr_path_layer[i_layer],track_posA);
            Evaluate(arr_path_layer[i_layer]+0.1,track_posB);


            vec_dir_vec_circle[0]   = track_posB[0] - track_posA[0];
            vec_dir_vec_circle[1]   = track_posB[1] - track_posA[1];
            vec_dir_vec_circle[2]   = 0.0;

            //printf("Kalman: vec_dir_vec_circle[0] = %4.3f, vec_dir_vec_circle[1] = %4.3f \n",vec_dir_vec_circle.X(),vec_dir_vec_circle.Y());

            vec_dir_vec_circle	   *= 1/vec_dir_vec_circle.Mag();

            //printf("Kalman: track_posAX = %4.3f, track_posBX = %4.3f, track_posAY = %4.3f, track_posBY = %4.3f\n",track_posA[0],track_posB[0],track_posA[1],track_posB[1]);

        //printf("Kalman: vec_dir_vec_circle[0] = %4.3f, vec_dir_vec_circle[1] = %4.3f \n",vec_dir_vec_circle.X(),vec_dir_vec_circle.Y());

            Double_t delta_x_local_global_circle 	= vec_dir_vec_circle.Dot((*vec_TV3_TRD_center[detector[i_layer]][0]));
            //printf("delta_x_local_global_circle = %4.3f \n",delta_x_local_global_circle);

            // vec_dir_vec_circle[i_layer] - this is a tangent to the circle track
            // vec_TV3_TRD_center this is a vector from {0,0,0} to the middle of a chamber - from TRD geometry
            // but all related to the circle track that you probably dont need

            Double_t sign_direction_impact_circle 	= TMath::Sign(1.0,delta_x_local_global_circle);
            Double_t impact_angle_circle		= vec_dir_vec_circle.Angle((*vec_TV3_TRD_center[detector[i_layer]][2]));  //i need something like that for circles

            if(impact_angle_circle > TMath::Pi()*0.5) impact_angle_circle -= TMath::Pi();
            impact_angle_circle = 0.5*TMath::Pi() - sign_direction_impact_circle*impact_angle_circle;


            Double_t delta_x_local_tracklet = vec_TV3_tracklet_vectors.Dot((*vec_TV3_TRD_center[detector[i_layer]][0]));

            Double_t sign_angle_circle = 1.0;
            if(delta_x_local_tracklet < delta_x_local_global_circle) sign_angle_circle = -1.0;

            Double_t Delta_angle_circle  = sign_angle_circle*vec_dir_vec_circle.Angle(vec_TV3_tracklet_vectors);
        //if(impact_angle_circle*TMath::RadToDeg() > 85.0 && impact_angle_circle*TMath::RadToDeg() < 90.0) printf("delta_x_local_tracklet = %4.3f, sign_angle_circle: %4.3f, Delta_angle_circle: %4.3f \n",delta_x_local_tracklet,sign_angle_circle,Delta_angle_circle*TMath::RadToDeg());


            //printf("   ---> vec_KF: {%4.3f, %4.3f, %4.3f}, vec_trkl: {%4.3f, %4.3f, %4.3f} \n",vec_dir_vec_circle.X(),vec_dir_vec_circle.Y(),vec_dir_vec_circle.Z(),vec_TV3_tracklet_vectors.X(),vec_TV3_tracklet_vectors.Y(),vec_TV3_tracklet_vectors.Z());

            if(Delta_angle_circle > TMath::Pi()*0.5)  Delta_angle_circle -= TMath::Pi();
            if(Delta_angle_circle < -TMath::Pi()*0.5) Delta_angle_circle += TMath::Pi();

            vec_TH2D_Delta_vs_impact[detector[i_layer]]   ->Fill(impact_angle_circle*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());
            vec_tp_Delta_vs_impact[detector[i_layer]] ->Fill(impact_angle_circle*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());
        //if(impact_angle_circle*TMath::RadToDeg() > 85.0 && impact_angle_circle*TMath::RadToDeg() < 90.0) printf("--> i_track: %d, i_layer: %d \nimpact angle: %4.3f, delta angle: %4.3f \n",i_track,i_layer,impact_angle_circle*TMath::RadToDeg(),Delta_angle_circle*TMath::RadToDeg());


            //draw here -----------------------------------------

#if defined(USEEVE)
            if(graphics)
            {
                TEveLine_vec_dir_vec_circle[i_track][i_layer]=new TEveLine();
                TEveLine_vec_dir_vec_circle[i_track][i_layer]->SetNextPoint(vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[0],vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[1],vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[2]);
                TEveLine_vec_dir_vec_circle[i_track][i_layer]->SetNextPoint(vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[0] + scale_length_vec*vec_dir_vec_circle.X(),vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[1] + scale_length_vec*vec_dir_vec_circle.Y(),vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[2]);
            }
#endif

            //detector[i_layer]        = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TRD_det();
            //Int_t i_sector = (Int_t)(detector/30);
            //Int_t i_stack  = (Int_t)(detector%30/6);
            //Int_t i_layer  = i_det%6;


#if defined(USEEVE)
            if(graphics)
            {
                HistName = "vec_dir_vec_circle_";
                HistName += i_track;
                HistName += "_";
                //HistName += i_sector;
                //HistName += "_";
                //HistName += i_stack;
                //HistName += "_";
                HistName += i_layer;

                TEveLine_vec_dir_vec_circle[i_track][i_layer]   ->SetName(HistName.Data());
                TEveLine_vec_dir_vec_circle[i_track][i_layer]   ->SetLineStyle(1);
                TEveLine_vec_dir_vec_circle[i_track][i_layer]   ->SetLineWidth(6);
                TEveLine_vec_dir_vec_circle[i_track][i_layer]   ->SetMainColor(kMagenta);
                gEve->AddElement(TEveLine_vec_dir_vec_circle[i_track][i_layer]);
            }
 
#endif

#if defined(USEEVE)
            if(graphics)
            {
                gEve->Redraw3D(kTRUE);
            }
#endif
            //--------------------------------------------------- 

#if 1 //----stuff for circle fit---------------------
            vec_TV2_points[i_layer].SetX(vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[0]);
            vec_TV2_points[i_layer].SetY(vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[1]);

            vec_Dt_digit_pos_cluster[i_layer][0] = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[0];
            vec_Dt_digit_pos_cluster[i_layer][1] = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[1];
            //---------------------------------------------------
#endif
        }


#if 1 //----stuff for circle fit---------------------

        Double_t par_circle[3] = {0.0,0.0,0.0};
        vector<Double_t> vec_phi;
        vec_phi.resize(6);

        for (Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            vec_phi[i_layer] = 0.0;
        }

        //n_layer_notempty = i_layer_notempty-1;  ?? check later

        if ((Int_t)vec_TV2_points.size() < 2)  //<2
        {
            //printf("!!! less than 3 points to fit circle - you got bad circle !!! \n");
            //continue;
        }

        Double_t p0[3] = {10,20,1};
        Double_t tracklets_min_circle = -1.0;

        //printf(" -------- MINOS Fit get_2D_global_circle_fit -------- \n");
        TVirtualFitter *min = TVirtualFitter::Fitter(0,3);

        Double_t arglist_A[1] = {-1};
        Double_t arglist_B[1] = {0};
        Double_t arglist_C[1] = {-1};
        //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
        //min->ExecuteCommand("SHOw FCNvalue",arglist_A,1);
        min->ExecuteCommand("SET PRINT",arglist_C,1);
        //min->ExecuteCommand("SET PRINT",arglist_B,1);
        min->ExecuteCommand("SET NOWarnings",arglist_B,1);

        //ROOT::Math::Functor f(&Ali_TRD_ST_Analyze::sum_distance_circ_point_2D,1);

        //min->SetFCN(f);
        min->SetFCN(sum_distance_circ_point_2D);

        //Double_t arglist_A[1] = {-1};
        //Double_t arglist_B[1] = {0};
        //Double_t arglist_C[1] = {-1};
        //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
        //min->ExecuteCommand("SHOw FCNvalue",arglist_A,1);
        //min->ExecuteCommand("SET NOWarnings",arglist_B,1);
        //min->ExecuteCommand("SET PRINT",arglist_C,1);

        Double_t arglist[10];
        arglist[0] = 3;

        Double_t pStart[3] = {0.0,0.0,0.0}; //= {1,1,1,1};

        // old circle fit possibly causing layer 1 and 4 bias and radius underestimation
        // Double_t x1 = vec_TV2_points[0].X();
        // Double_t y1 = vec_TV2_points[0].Y();
        // Double_t x2 = vec_TV2_points[1].X();
        // Double_t y2 = vec_TV2_points[1].Y();
        // Double_t x3 = vec_TV2_points[2].X();
        // Double_t y3 = vec_TV2_points[2].Y();
        vector <Int_t> notempty;
        for (Int_t i_layerB = 0; i_layerB < 6; i_layerB++)
        {
            if (vec_TV2_points[i_layerB].X()!=-999.9) notempty.push_back(i_layerB);
        }

        Double_t x1 = vec_TV2_points[notempty[0]].X();
        Double_t y1 = vec_TV2_points[notempty[0]].Y();
        Double_t x2 = vec_TV2_points[notempty[notempty.size()-2]].X();
        Double_t y2 = vec_TV2_points[notempty[notempty.size()-2]].Y();
        Double_t x3 = vec_TV2_points[notempty[notempty.size()-1]].X();
        Double_t y3 = vec_TV2_points[notempty[notempty.size()-1]].Y();

        notempty.clear();
        
        //maybe implement 3rd point = 0,0,0?


        Double_t A_help = x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2;
        Double_t B_help = (x1*x1 + y1*y1)*(y3 - y2) + (x2*x2 + y2*y2)*(y1 - y3) + (x3*x3 + y3*y3)*(y2 - y1);
        Double_t C_help = (x1*x1 + y1*y1)*(x2 - x3) + (x2*x2 + y2*y2)*(x3 - x1) + (x3*x3 + y3*y3)*(x1 - x2);
        Double_t D_help = (x1*x1 + y1*y1)*(x3*y2 - x2*y3) + (x2*x2 + y2*y2)*(x1*y3 - x3*y1) + (x3*x3 + y3*y3)*(x2*y1 - x1*y2);

        //printf("A_help: {%4.3f, %4.3f}, B_help: {%4.3f, %4.3f}, C_help: {%4.3f, %4.3f}  \n",A_help,B_help,C_help);

        Double_t a_param = 0.0;
        Double_t b_param = 0.0;
        Double_t R_param = 0.0;
      

        if(A_help != 0.0)
        {
            a_param = -B_help/(2*A_help);
            b_param = -C_help/(2*A_help);
            R_param = TMath::Sqrt(TMath::Power(x1 - a_param,2.0) + TMath::Power(y1 - b_param,2.0));
        }

        //printf("pointA: {%4.3f, %4.3f}, pointB: {%4.3f, %4.3f}, pointC: {%4.3f, %4.3f}, circle(a,b,R): {%4.3f, %4.3f, %4.3f} \n",x1,y1,x2,y2,x3,y3,a_param,b_param,R_param);

        pStart[0] = a_param;
        pStart[1] = b_param;
        pStart[2] = R_param;


        //cout << "pStart[0]" << pStart[0] << endl;
        //cout << "pStart[1]" << pStart[1] << endl;
        //cout << "pStart[2]" << pStart[2] << endl;

        min->SetParameter(0,"a_param",pStart[0],0.01,0,0);
        min->SetParameter(1,"b_param",pStart[1],0.01,0,0);
        min->SetParameter(2,"R_param",pStart[2],0.01,0,0);

        arglist[0] = 1000; // number of function calls
        arglist[1] = 0.001; // tolerance

        //min->ExecuteCommand("MIGRAD",arglist,2);


        //min->ExecuteCommand("MINOS",arglist,2);
        min->ExecuteCommand("MINIMIZE",arglist,2);

        //min->ExecuteCommand("SET PRIntout",arglist_A,1); // http://www.fresco.org.uk/minuit/cern/node18.html
        //min->ExecuteCommand("SET NOWarnings",arglist_B,1);
        //min->ExecuteCommand("SET PRINT",arglist_C,1);
        //printf(" -------- END MINOS Fit get_2D_global_circle_fit -------- \n");

        //if (minos) min->ExecuteCommand("MINOS",arglist,0);
        Int_t nvpar,nparx;
        Double_t amin,edm, errdef;
        min->GetStats(amin,edm,errdef,nvpar,nparx);
        //min->PrintResults(1,amin);

        Double_t parFit_circ[3];

        for(Int_t i_par = 0; i_par <3; ++i_par)
        {
            parFit_circ[i_par] = min->GetParameter(i_par);
            //parFit[i] = pStart[i];
        }

        for(Int_t i_par = 0; i_par <3; ++i_par)
        {
            par_circle[i_par] = parFit_circ[i_par];
        }

        //printf("amin: %4.3f \n",amin);
        tracklets_min_circle = amin;


        // Calculate direction vectors from 2D circle fit
        TVector2 dir_vector;
        TVector2 circle_center_vector;
        circle_center_vector.SetX(parFit_circ[0]);
        circle_center_vector.SetY(parFit_circ[1]);
        vector <Int_t> notempt;

        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            if(vec_kalman_TRD_trackets[i_track][i_layer] == NULL ) continue;
            dir_vector = vec_TV2_points[i_layer];

            //printf("dir_vector.X(): %4.3f, dir_vector.Y(): %4.3f \n",dir_vector.X(),dir_vector.Y());

            dir_vector -= circle_center_vector;
            Double_t phi_val = dir_vector.Phi();

            //printf("phi_val: %4.3f \n",phi_val);

            //vec_phi.push_back(phi_val);
            vec_phi[i_layer] = phi_val;
            notempt.push_back(i_layer);

            //printf("i_point: %d, vec_phi[i_point]: %4.3f, \n",i_point,vec_phi[i_point]);
        }

        Double_t delta_phi = -(vec_phi[notempt[1]] - vec_phi[notempt[0]])/1000.0;

        //printf("delta_phi: %4.3f, vec_phi[1]: %4.3f, vec_phi[0]: %4.3f  \n",delta_phi,vec_phi[1],vec_phi[0]);

        vec_dir_vec_circle_circle.resize(6);

        for (Int_t i_layer = 0; i_layer <6; ++i_layer)
        {
            vec_dir_vec_circle_circle[i_layer].SetXYZ(-999.0,-999.0,-999.0);   //because i maybe need it to have all 6 dimensions in Calibrate
        }

#if defined(USEEVE)
        if(graphics)
        {
            TEveLine_vec_dir_vec_circle_circle[i_track].resize(6);
        }
#endif

#if defined(USEEVE) // Draw circle track ----------------------
        if(graphics)
        {
            TEveLine_circle[i_track]=new TEveLine();

            Double_t delta_i_y = par_circle[2]/200.0;

            Double_t radius_circle = 0.0;

            for(Double_t i_sign = -1.0; i_sign <= +1.0; i_sign += 2.0)
            {
                for(Double_t i_y = ((par_circle[1] - par_circle[2]) + delta_i_y); i_y < ((par_circle[1] + par_circle[2]) - delta_i_y); (i_y += delta_i_y))
                {
                    Double_t i_x = i_sign*TMath::Sqrt(TMath::Power(par_circle[2],2.0) - TMath::Power(i_y - par_circle[1],2.0)) + par_circle[0];
                    //radius_circle = TMath::Sqrt( TMath::Power(i_x,2) + TMath::Power(i_y,2) );
                    //if(isnan(radius_circle)) continue;
                    //if(radius_circle > 500.0) continue; // 500.0
                    //if(radius_circle < 280.0) continue; // 280.0
                    if(i_x > 500.0 || i_y > 500.0) continue; // 500.0
                    if(i_x < -500.0 || i_y < -500.0) continue; // 500.0
                    //if(i_x < 280.0) continue; // 280.0
                    //if(i_y > 500.0) continue; // 500.0
                    //if(i_y < 280.0) continue; // 280.0
                    Int_t non_null = 0;
                    for (Int_t i_layerB = 0; i_layerB < 6; i_layerB++)
                    {
                        if (vec_kalman_TRD_trackets[i_track][i_layerB]!=NULL) 
                            {
                                non_null = i_layerB;
                                break;
                            }
                    }
                    TEveLine_circle[i_track] ->SetNextPoint(i_x,i_y,vec_kalman_TRD_trackets[i_track][non_null] ->get_TV3_offset()[2]);
                    //if(fabs(track_pos[2]) > 430.0) break;
                    //if(isnan(radius_circle)) break;
                    //if(radius_circle > 500.0) break; // 500.0
                    //if(radius_circle < 280.0) break; // 280.0
                    if(i_x > 500.0 || i_y > 500.0) break; // 500.0
                    if(i_x < -500.0 || i_y < -500.0) break; // 500.0
                    //if(i_x < 280.0) break; // 280.0
                    //if(i_y > 500.0) break; // 500.0
                    //if(i_y < 280.0) break; // 280.0
                    //printf("circle point: {%4.3f, %4.3f} \n",i_x,i_y);

                }
            }

            HistName = "track_circle_";
            HistName += i_track;
            HistName += "_";

            TEveLine_circle[i_track]   ->SetName(HistName.Data());
            TEveLine_circle[i_track]   ->SetLineStyle(1);
            // TEveLine_circle[i_track] ->SetLineWidth(3);
            TEveLine_circle[i_track]   ->SetLineWidth(1);
            //TEveLine_circle[i_track] ->SetLineColor(kTeal+2);
            TEveLine_circle[i_track]   ->SetLineColor(kOrange);
            //TEveLine_circle[i_track]   ->DrawClone("l");
            //printf("test 1 \n");

            gEve->AddElement(TEveLine_circle[i_track]);
            //printf("test 2 \n");
            gEve->Redraw3D(kTRUE);

        }


        counter_circ++;
        //printf("counter_circ = %d \n",counter_circ);
#endif //------------------------------------------------------
        
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        //for(Int_t i_point = 0; i_point < 6; i_point++)
        {
            if(vec_kalman_TRD_trackets[i_track][i_layer] == NULL ) continue;
            //printf("from here \n");

            vec_TV3_tracklet_vectors     = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_dir();
            vec_TV3_tracklet_vectors[2]  = 0.0;
            vec_TV3_tracklet_vectors    *= 1/ vec_TV3_tracklet_vectors.Mag();
            Double_t phi_val = vec_phi[i_layer];

            Double_t x_val   =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val);
            Double_t y_val   =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val);

            Double_t x_valB  =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val + delta_phi);
            Double_t y_valB  =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val + delta_phi);

            if ((y_valB*y_valB + x_valB*x_valB) < (y_val*y_val + x_val*x_val))
            {
                x_valB  =  parFit_circ[0] + parFit_circ[2] * TMath::Cos(phi_val - delta_phi);
                y_valB  =  parFit_circ[1] + parFit_circ[2] * TMath::Sin(phi_val - delta_phi);
            }
            
            //TVector3 dir_vec_circle;               //maybe i can use this one

            vec_dir_vec_circle_circle[i_layer].SetX(x_valB - x_val);
            vec_dir_vec_circle_circle[i_layer].SetY(y_valB - y_val);
            vec_dir_vec_circle_circle[i_layer].SetZ(0.0);

            //printf("vec_dir_vec_circle_notempty[i_point].X: %4.3f, vec_dir_vec_circle_notempty[i_point].Y: %4.3f, vec_dir_vec_circle_notempty[i_point].Z: %4.3f \n",vec_dir_vec_circle_notempty[i_point].X(),vec_dir_vec_circle_notempty[i_point].Y(),vec_dir_vec_circle_notempty[i_point].Z());
            //printf("y_valB*y_valB + x_valB*x_valB = %4.3f,  y_val*y_val + x_val*x_val = %4.3f \n",y_valB*y_valB + x_valB*x_valB,y_val*y_val + x_val*x_val);


            //printf("x_val: %4.3f, x_valB: %4.3f, y_val: %4.3f, y_valB: %4.3f \n",x_val,x_valB,y_val,y_valB); //  probably correct
            //printf("vec_dir_vec_circle_circle[0] = %4.3f, vec_dir_vec_circle_circle[1] = %4.3f \n",vec_dir_vec_circle_circle[i_layer].X(),vec_dir_vec_circle_circle[i_layer].Y());


            //if(vec_dir_vec_circle_circle[i_layer].Mag() > 0.0)
            //{
                vec_dir_vec_circle_circle[i_layer] *= 1.0/vec_dir_vec_circle_circle[i_layer].Mag();
            //}

            //printf("vec_dir_vec_circle_circle[0] = %4.3f, vec_dir_vec_circle_circle[1] = %4.3f \n",vec_dir_vec_circle_circle[i_layer].X(),vec_dir_vec_circle_circle[i_layer].Y());
            
            Double_t delta_x_local_global_circle_circle    = vec_dir_vec_circle_circle[i_layer].Dot((*vec_TV3_TRD_center[detector[i_layer]][0]));
            //printf("delta_x_local_global_circle_circle = %4.3f \n",delta_x_local_global_circle_circle);
            // vec_dir_vec_circle[i_layer] - this is a tangent to the circle track
            // vec_TV3_TRD_center this is a vector from {0,0,0} to the middle of a chamber - from TRD geometry
            // but all related to the circle track that you probably dont need

            Double_t sign_direction_impact_circle_circle   = TMath::Sign(1.0,delta_x_local_global_circle_circle);
            Double_t impact_angle_circle_circle        = vec_dir_vec_circle_circle[i_layer].Angle((*vec_TV3_TRD_center[detector[i_layer]][2]));  //i need something like that for circles

            if(impact_angle_circle_circle > TMath::Pi()*0.5) impact_angle_circle_circle -= TMath::Pi();
            impact_angle_circle_circle = 0.5*TMath::Pi() - sign_direction_impact_circle_circle*impact_angle_circle_circle;


            Double_t delta_x_local_tracklet_circle_circle = vec_TV3_tracklet_vectors.Dot((*vec_TV3_TRD_center[detector[i_layer]][0]));

            Double_t sign_angle_circle_circle = 1.0;
            if(delta_x_local_tracklet_circle_circle < delta_x_local_global_circle_circle) sign_angle_circle_circle = -1.0;

            Double_t Delta_angle_circle_circle  = sign_angle_circle_circle*vec_dir_vec_circle_circle[i_layer].Angle(vec_TV3_tracklet_vectors);

            //printf("   ---> vec_KF: {%4.3f, %4.3f, %4.3f}, vec_trkl: {%4.3f, %4.3f, %4.3f} \n",vec_dir_vec_circle.X(),vec_dir_vec_circle.Y(),vec_dir_vec_circle.Z(),vec_TV3_tracklet_vectors.X(),vec_TV3_tracklet_vectors.Y(),vec_TV3_tracklet_vectors.Z());
        //if(impact_angle_circle_circle*TMath::RadToDeg() > 85.0 && impact_angle_circle_circle*TMath::RadToDeg() < 90.0) printf("delta_x_local_tracklet_circle_circle = %4.3f, sign_angle_circle_circle: %4.3f, Delta_angle_circle_circle: %4.3f \n",delta_x_local_tracklet_circle_circle,sign_angle_circle_circle,Delta_angle_circle_circle*TMath::RadToDeg());

            if(Delta_angle_circle_circle > TMath::Pi()*0.5)  Delta_angle_circle_circle -= TMath::Pi();
            if(Delta_angle_circle_circle < -TMath::Pi()*0.5) Delta_angle_circle_circle += TMath::Pi();

        //if(impact_angle_circle_circle*TMath::RadToDeg() > 85.0 && impact_angle_circle_circle*TMath::RadToDeg() < 90.0) printf("i_layer: %d, impact angle circle: %4.3f, delta angle circle: %4.3f \n",i_layer,impact_angle_circle_circle*TMath::RadToDeg(),Delta_angle_circle_circle*TMath::RadToDeg());

            vec_TH2D_Delta_vs_impact_circle[detector[i_layer]]   ->Fill(impact_angle_circle_circle*TMath::RadToDeg(),Delta_angle_circle_circle*TMath::RadToDeg());
            vec_tp_Delta_vs_impact_circle[detector[i_layer]]     ->Fill(impact_angle_circle_circle*TMath::RadToDeg(),Delta_angle_circle_circle*TMath::RadToDeg());

            //draw here -----------------------------------------

#if defined(USEEVE)
            if(graphics)
            {                
                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]=new TEveLine();
                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]->SetNextPoint(vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[0],vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[1],vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[2]);                
                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]->SetNextPoint(vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[0] + scale_length_vec*vec_dir_vec_circle_circle[i_layer].X(),vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[1] + scale_length_vec*vec_dir_vec_circle_circle[i_layer].Y(),vec_kalman_TRD_trackets[i_track][i_layer] ->get_TV3_offset()[2]);
            }
#endif

            //detector[i_layer]        = vec_kalman_TRD_trackets[i_track][i_layer] ->get_TRD_det();
            //Int_t i_sector = (Int_t)(detector/30);
            //Int_t i_stack  = (Int_t)(detector%30/6);
            //Int_t i_layer  = i_det%6;


#if defined(USEEVE)
            if(graphics)
            {
                HistName = "vec_dir_vec_circle_circle_";
                HistName += i_track;
                HistName += "_";
                //HistName += i_sector;
                //HistName += "_";
                //HistName += i_stack;
                //HistName += "_";
                HistName += i_layer;

                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]   ->SetName(HistName.Data());
                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]   ->SetLineStyle(1);
                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]   ->SetLineWidth(6);
                TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]   ->SetMainColor(kCyan);
                gEve->AddElement(TEveLine_vec_dir_vec_circle_circle[i_track][i_layer]);
            }
 
#endif

#if defined(USEEVE)
            if(graphics)
            {
                gEve->Redraw3D(kTRUE);
            }
#endif


        }
        vec_TV2_points.clear();
        vec_Dt_digit_pos_cluster.clear();
        //---------------------------------------------------
#endif
    }
}

//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

void Ali_TRD_ST_Analyze::Draw_n_Save_Calibration(TString out_dir, TString out_file_name_calib)
{
   
    #if 1 //------------for Kalman hists----------------------------

    vector<TCanvas*> vec_can_Delta_vs_impact;
    char NoP[50];

    vec_can_Delta_vs_impact.resize(6); // 6 sector blocks with 3 sectors each (18)

    TH1D* h_dummy_Delta_vs_impact = new TH1D("h_dummy_Delta_vs_impact","h_dummy_Delta_vs_impact",90,50,140);
    h_dummy_Delta_vs_impact->SetStats(0);
    h_dummy_Delta_vs_impact->SetTitle("");
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_Delta_vs_impact->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_Delta_vs_impact->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_Delta_vs_impact->GetXaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact->GetYaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact->GetXaxis()->CenterTitle();
    h_dummy_Delta_vs_impact->GetYaxis()->CenterTitle();
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitle("impact angle");
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitle("#Delta #alpha");
    h_dummy_Delta_vs_impact->GetXaxis()->SetRangeUser(70,110);
    h_dummy_Delta_vs_impact->GetYaxis()->SetRangeUser(-24,24);

    #if 1
    Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    for(Int_t i_sec_block = 0; i_sec_block < 6; i_sec_block++) //<6
    {
        HistName = "h_dummy_Delta_vs_impact_";
        HistName += i_sec_block;
        vec_can_Delta_vs_impact[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_can_Delta_vs_impact[i_sec_block] ->Divide(5,3); // x = stack, y = sector

        for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
        {
            Int_t i_sector = i_sec_block + 6*i_sec_sub;
            for(Int_t i_stack = 0; i_stack < 5; i_stack++)
            {
                Int_t iPad = i_sec_sub*5 + i_stack + 1;
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTicks(1,1);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetGrid(0,0);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetFillColor(10);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
                vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad);
                h_dummy_Delta_vs_impact->Draw("h");

                for(Int_t i_layer = 0; i_layer < 6; i_layer++)
                {
                    Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                    // printf("detector: %d \n",i_detector);
                    vec_tp_Delta_vs_impact[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
                    vec_tp_Delta_vs_impact[i_detector] ->SetLineWidth(2);
                    vec_tp_Delta_vs_impact[i_detector] ->SetLineStyle(1);
                    vec_tp_Delta_vs_impact[i_detector] ->Draw("same hl");

                    HistName = "";
                    sprintf(NoP,"%4.0f",(Double_t)i_detector);
                    HistName += NoP;
                    //plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                }
            }
        }
    }


    #endif
    #endif //--------------------------------------------------

    #if 1 //------------for circle hists----------------------------

    vector<TCanvas*> vec_can_Delta_vs_impact_circle;
    //char NoP[50];

    vec_can_Delta_vs_impact_circle.resize(6); // 6 sector blocks with 3 sectors each (18)

    TH1D* h_dummy_Delta_vs_impact_circle = new TH1D("h_dummy_Delta_vs_impact_circle","h_dummy_Delta_vs_impact_circle",90,50,140);
    h_dummy_Delta_vs_impact_circle->SetStats(0);
    h_dummy_Delta_vs_impact_circle->SetTitle("");
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact_circle->GetXaxis()->CenterTitle();
    h_dummy_Delta_vs_impact_circle->GetYaxis()->CenterTitle();
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetTitle("impact angle");
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetTitle("#Delta #alpha");
    h_dummy_Delta_vs_impact_circle->GetXaxis()->SetRangeUser(70,110);
    h_dummy_Delta_vs_impact_circle->GetYaxis()->SetRangeUser(-24,24);

    #if 1
    //Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    for(Int_t i_sec_block = 0; i_sec_block < 1; i_sec_block++) //<6
    {
        HistName = "vec_can_Delta_vs_impact_circle_";
        HistName += i_sec_block;
        vec_can_Delta_vs_impact_circle[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

        vec_can_Delta_vs_impact_circle[i_sec_block] ->Divide(5,3); // x = stack, y = sector

        for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
        {
            Int_t i_sector = i_sec_block + 6*i_sec_sub;
            for(Int_t i_stack = 0; i_stack < 5; i_stack++)
            {
                Int_t iPad = i_sec_sub*5 + i_stack + 1;
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTicks(1,1);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetGrid(0,0);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetFillColor(10);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
                vec_can_Delta_vs_impact_circle[i_sec_block] ->cd(iPad);
                h_dummy_Delta_vs_impact_circle->Draw("h");

                for(Int_t i_layer = 0; i_layer < 6; i_layer++)
                {
                    Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                    // printf("detector: %d \n",i_detector);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineWidth(2);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->SetLineStyle(1);
                    vec_tp_Delta_vs_impact_circle[i_detector] ->Draw("same hl");

                    HistName = "";
                    sprintf(NoP,"%4.0f",(Double_t)i_detector);
                    HistName += NoP;
                    //plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
                }
            }
        }
    }
    #endif 
    #endif//--------------------------------------------------
    

    HistName = out_dir;
    HistName += "/";
    HistName += out_file_name_calib;
    TFile* outputfile_hist = new TFile(HistName.Data(),"RECREATE");

    printf("Write data to output file \n");
    //TFile* h_detector_hit_outputfile = new TFile("./h_detector_hit.root","RECREATE");
    // h_detector_hit_outputfile ->cd();
    //  h_detector_hit->Write();

    // THIS NEEDED
    outputfile_hist ->cd();
    outputfile_hist ->mkdir("Delta_impact");
    outputfile_hist ->cd("Delta_impact");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact[i_det]   ->Write();
        vec_TH2D_Delta_vs_impact[i_det] ->Write();
    }

    outputfile_hist ->cd();
    outputfile_hist ->mkdir("Delta_impact_circle");
    outputfile_hist ->cd("Delta_impact_circle");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_tp_Delta_vs_impact_circle[i_det]   ->Write();
        vec_TH2D_Delta_vs_impact_circle[i_det] ->Write();
    }

    printf("All data written \n");

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

Float_t Ali_TRD_ST_Analyze::primary_vertex_dca(Int_t i_track)
{
	TVector3 zero;
	Float_t pathA_dca, dcaAB_dca;
    fHelixAtoPointdca(zero,vec_helices_TRD[i_track],pathA_dca,dcaAB_dca); // new helix to point dca calculation
    return dcaAB_dca;
}	

//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

void Ali_TRD_ST_Analyze::hists_pv_dca()
{
    //printf("Ali_TRD_ST_Analyze::Calibrate() \n");
	//cout<<"the init variable is "<<inits_vec_tp_pvdca_vs_sector<<endl;
    /*if(inits_vec_tp_pvdca_vs_sector == 0)
    
	//if(vec_tp_pvdca_vs_sector == NULL)
    {	*/
    /*    vec_tp_pvdca_vs_sector.resize(90);
        vec_TH2D_pvdca_vs_sector.resize(90);
        
        for (Int_t i_det = 0; i_det < 90; i_det++)
        {
		vec_tp_pvdca_vs_sector = new TProfile("vec_th1d_pvdca_vs_sector","vec_th1d_pvdca_vs_sector_%d",90,0,89);
    	vec_TH2D_pvdca_vs_sector = new TH2D("vec_th2d_pvdca_vs_sector","vec_th2d_pvdca_vs_sector_%d",90,0,89,100,0,50);
          //#
		//cout<<"initialized"<<endl;
		inits_vec_tp_pvdca_vs_sector++;
        //}	
		
   }*/
	 	
	Int_t detector;
	Int_t sector;
	Float_t dca;	
	//cout<<"number of tracks "<<(Int_t)mHelices_kalman.size()<<endl;
		
	for(Int_t i_track = 0; i_track < (Int_t)mHelices_kalman.size(); i_track++)
    {
    
		dca 		= primary_vertex_dca(i_track);
		Int_t tracklet =0;
		while( vec_kalman_TRD_trackets[i_track][tracklet] ==NULL) tracklet++;
		//cout<<"we got here this many times "<<i_track<<endl;
		
		detector	= vec_kalman_TRD_trackets[i_track][tracklet] ->get_TRD_det();
		sector 		= (Int_t) detector/6 ;
		//cout<<"we got here this many times "<<i_track<<endl;
		//cout<<"sector "<<sector<<" dca "<<dca<<endl;
		vec_tp_pvdca_vs_sector ->Fill(sector,dca);
		vec_TH2D_pvdca_vs_sector   ->Fill(sector,dca);
		//cout<<"we got here this many times "<<i_track<<endl;
		
	}            
	//cout<<"we got here!!"<<endl;
}

//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------

void Ali_TRD_ST_Analyze::Draw_n_Save_hists_pv_dca(TString out_dir, TString out_file_name_calib)
{
	TCanvas* vec_can_pvdca_vs_sector;
    char NoP[50];

    //vec_can_Delta_vs_impact.resize(6); // 6 sector blocks with 3 sectors each (18)

  /*  TH1D* h_dummy_Delta_vs_impact = new TH1D("h_dummy_Delta_vs_impact","h_dummy_Delta_vs_impact",90,50,140);
    h_dummy_Delta_vs_impact->SetStats(0);
    h_dummy_Delta_vs_impact->SetTitle("");
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitleOffset(0.85);
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitleOffset(0.78);
    h_dummy_Delta_vs_impact->GetXaxis()->SetLabelOffset(0.0);
    h_dummy_Delta_vs_impact->GetYaxis()->SetLabelOffset(0.01);
    h_dummy_Delta_vs_impact->GetXaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact->GetYaxis()->SetLabelSize(0.08);
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitleSize(0.08);
    h_dummy_Delta_vs_impact->GetXaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact->GetYaxis()->SetNdivisions(505,'N');
    h_dummy_Delta_vs_impact->GetXaxis()->CenterTitle();
    h_dummy_Delta_vs_impact->GetYaxis()->CenterTitle();
    h_dummy_Delta_vs_impact->GetXaxis()->SetTitle("impact angle");
    h_dummy_Delta_vs_impact->GetYaxis()->SetTitle("#Delta #alpha");
    h_dummy_Delta_vs_impact->GetXaxis()->SetRangeUser(70,110);
    h_dummy_Delta_vs_impact->GetYaxis()->SetRangeUser(-24,24);*/

    //#if 1
    Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

    //for(Int_t i_sec_block = 0; i_sec_block < 1; i_sec_block++) //<6
    //{
        HistName = "h_dummy_pvdca_vs_sector";
      //  HistName += i_sec_block;
        vec_can_pvdca_vs_sector = new TCanvas(HistName.Data(),HistName.Data(),1600,1000);

        //vec_can_Delta_vs_impact[i_sec_block] ->Divide(5,3); // x = stack, y = sector

       /* for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
        {
            Int_t i_sector = i_sec_block + 6*i_sec_sub;
            for(Int_t i_stack = 0; i_stack < 5; i_stack++)
            {
                Int_t iPad = i_sec_sub*5 + i_stack + 1;*/
                vec_can_pvdca_vs_sector ->SetTicks(1,1);
                vec_can_pvdca_vs_sector ->SetGrid(0,0);
                vec_can_pvdca_vs_sector ->SetFillColor(10);
                vec_can_pvdca_vs_sector ->SetRightMargin(0.01);
                vec_can_pvdca_vs_sector ->SetTopMargin(0.01);
                vec_can_pvdca_vs_sector ->SetBottomMargin(0.2);
                vec_can_pvdca_vs_sector ->SetLeftMargin(0.2);
				
              //  h_dummy_Delta_vs_impact->Draw();

               // for(Int_t i_layer = 0; i_layer < 6; i_layer++)
               // {
                 //   Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                    // printf("detector: %d \n",i_detector);
                    vec_tp_pvdca_vs_sector ->SetLineColor(arr_color_layer[1]);
                    vec_tp_pvdca_vs_sector ->SetLineWidth(2);
                    vec_tp_pvdca_vs_sector ->SetLineStyle(1);
                    vec_tp_pvdca_vs_sector ->Draw();
					vec_tp_pvdca_vs_sector->GetXaxis()->SetTitle("Stack Nbr");
    				vec_tp_pvdca_vs_sector->GetYaxis()->SetTitle("dca to primary vertex (cm)");
    				vec_tp_pvdca_vs_sector->SetStats(0);
   			 		vec_tp_pvdca_vs_sector->SetTitle("");


             /*       HistName = "";
                    sprintf(NoP,"%4.0f",(Double_t)i_detector);
                    HistName += NoP;*/
                    //plotTopLegend((char*)HistName.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
               // }
            //}
        //}
    //}
	HistName = out_dir;
    HistName += "/";
    HistName += out_file_name_calib;
    TFile* outputfile_hist = new TFile(HistName.Data(),"RECREATE");

    printf("Write data to output file \n");
    //TFile* h_detector_hit_outputfile = new TFile("./h_detector_hit.root","RECREATE");
    // h_detector_hit_outputfile ->cd();
    //  h_detector_hit->Write();

    // THIS NEEDED
    outputfile_hist ->cd();
    outputfile_hist ->mkdir("PV_dca");
    outputfile_hist ->cd("PV_dca");
	vec_tp_pvdca_vs_sector  ->Write();
	//vec_TH2D_Delta_vs_impact[i_det] ->Write();
    
    printf("All data written \n");
}
