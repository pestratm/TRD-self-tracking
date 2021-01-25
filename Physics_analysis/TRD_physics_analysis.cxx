
#include "TRD_physics_analysis.h"

// implementation of the main code

// Something like with your new classes


//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    //gStyle->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
									 TVector3 &point)
{
  // calculates the minimum distance of a point to a straight given as parametric straight x = base + n * dir

  if (!(dir.Mag()>0))
    {
      return -1000000.;
    }
  
  TVector3 diff = base-point;

  TVector3 cross = dir.Cross(diff);
  
  return cross.Mag()/dir.Mag();
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TVector3 calculateDCA_vec_StraightToPoint(TVector3 &base, TVector3 &dir, TVector3 &point)
{
  // calculates the minimum distance vector of a point to a straight given as parametric straight x = base + n * dir

    TVector3 diff = base-point;
    TVector3 dir_norm = dir;
    dir_norm *= (1.0/dir.Mag());
    Double_t proj_val = diff.Dot(dir_norm);
    TVector3 proj_dir = dir_norm;
    proj_dir *= proj_val;

    TVector3 dist_vec = proj_dir - diff;

    return dist_vec;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
TCanvas* Draw_2D_histo_and_canvas(TH2D* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.22);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(max_val > min_val)
    {
        hist->GetZaxis()->SetRangeUser(min_val,max_val);
    }
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TCanvas* Draw_1D_histo_and_canvas(TH1D* hist, TString name, Int_t x_size, Int_t y_size,
                              Double_t min_val, Double_t max_val, TString option)
{
    TCanvas* canvas = new TCanvas(name.Data(),name.Data(),10,10,x_size,y_size);
    canvas->SetFillColor(10);
    canvas->SetTopMargin(0.1);
    canvas->SetBottomMargin(0.2);
    canvas->SetRightMargin(0.05);
    canvas->SetLeftMargin(0.2);
    canvas->SetTicks(1,1);
    canvas->SetGrid(0,0);

    hist->SetStats(0);
    hist->SetTitle("");
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(1.2);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetNdivisions(505,'N');
    hist->GetYaxis()->SetNdivisions(505,'N');
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    if(min_val != max_val) hist->GetYaxis()->SetRangeUser(min_val,max_val);
    hist->DrawCopy(option.Data());

    return canvas;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
Float_t correct_pT_according_to_correlation_plot(Float_t pT_raw, TVector3 par)
{
    Float_t pT_new = par[0] + par[1]*pT_raw + par[2]*pT_raw*pT_raw;

    return pT_new;
}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
Ali_TRD_physics_analysis::Ali_TRD_physics_analysis(TString out_dir, TString out_file_name, Int_t graphics)
{
    SetRootGraphicStyle();

    HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    printf("test printf: %s \n",HistName.Data());

    TH2_vertex_photon_XY     = new TH2D("TH2_vertex_photon_XY","TH2_vertex_photon_XY",1600,-400,400,1600,-400,400);
    vec_TH1_vertex_photon_radius.resize(5);
    for(Int_t i_hist = 0; i_hist < (Int_t)vec_TH1_vertex_photon_radius.size(); i_hist++)
    {
        HistName = "vec_TH1_vertex_photon_radius_";
        HistName += i_hist;
        vec_TH1_vertex_photon_radius[i_hist] = new TH1D(HistName.Data(),HistName.Data(),1200,0,400.0);
    }
    TH2D* h2D_dEdx_vs_mom    = new TH2D("h2D_dEdx_vs_mom","h2D_dEdx_vs_mom",200,0,5,200,0,150.0);

    vec_TH1_mass_pi0.resize(5);
    vec_TH1_mass_pi0_ME.resize(5);

    vec_TLV_photon_mixed_events.clear();
    vec_nsigma_electron_mixed_events.clear();

    for(Int_t i_hist = 0; i_hist < (Int_t)vec_TH1_mass_pi0_ME.size(); i_hist++)
    {
        HistName = "vec_TH1_mass_pi0_ME_";
        HistName += i_hist;
        vec_TH1_mass_pi0_ME[i_hist] = new TH1D(HistName.Data(),HistName.Data(),300,0,6.0);
    }

    for(Int_t i_hist = 0; i_hist < (Int_t)vec_TH1_mass_pi0.size(); i_hist++)
    {
        HistName = "vec_TH1_mass_pi0_";
        HistName += i_hist;
        vec_TH1_mass_pi0[i_hist] = new TH1D(HistName.Data(),HistName.Data(),300,0,6.0);
    }

    TH1_angle_between_photons = new TH1D("TH1_angle_between_photons","TH1_angle_between_photons",180,-180,180);
    TH2D_angle_photons_vs_inv_mass = new TH2D("TH2D_angle_photons_vs_inv_mass","TH2D_angle_photons_vs_inv_mass",1200,0,6.0,180,0,180);
    th1d_angle_between = new TH1D("th1d_angle_between","th1d_angle_between",180,-180,180);

    h_angle_AB_raw = new TH1D("h_angle_AB_raw","h_angle_AB_raw",180,-180.0,180.0);
    h_angle_AB_raw_cutA = new TH1D("h_angle_AB_raw_cutA","h_angle_AB_raw_cutA",180,-180.0,180.0);
    h_angle_AB_raw_cutB = new TH1D("h_angle_AB_raw_cutB","h_angle_AB_raw_cutB",180,-180.0,180.0);

    TFile* pt_corr      = TFile::Open("./pt_corr.root");
    tg_pol2_par = (TGraph*)pt_corr->Get("tg_pol2_params");



    for(Int_t i_par = 0; i_par < 3; i_par++)
    {
        par_pT_corr_neg[i_par] = tg_pol2_par->GetY()[i_par];
        par_pT_corr_pos[i_par] = tg_pol2_par->GetY()[i_par+3];
    }

    #if defined(USEEVE)
    if(graphics)
    {
        SetRootGraphicStyle();
        //--------------------------
        printf("Ana_sec_vertices started \n");

        //gStyle->SetPalette(kDarkBodyRadiator); // https://root.cern.ch/doc/master/classTColor.html
        //gStyle->SetPalette(kInvertedDarkBodyRadiator); // https://root.cern.ch/doc/master/classTColor.html
        gStyle->SetPalette(kTemperatureMap); // https://root.cern.ch/doc/master/classTColor.html
        //--------------------------



        //--------------------------
        // Load TRD geometry
        TEveManager::Create();

        TH1I* h_good_bad_TRD_chambers;
        h_good_bad_TRD_chambers = new TH1I("h_good_bad_TRD_chambers","h_good_bad_TRD_chambers",540,0,540);
        TFile* file_TRD_QA_flags = TFile::Open("../Data/chamber_QC_flags.root");
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

        vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
        TFile* file_TRD_geom = TFile::Open("../Data/TRD_Geom.root");
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

        vector<TEveBox*> vec_eve_TRD_detector_box;
        vec_eve_TRD_detector_box.resize(540);
        vector< vector<TPolyLine*> > vec_PL_TRD_det_2D;
        vec_PL_TRD_det_2D.resize(18); // sectors
        for(Int_t i_sector = 0; i_sector < 18; i_sector++)
        {
            vec_PL_TRD_det_2D[i_sector].resize(6); // layers
            for(Int_t i_layer = 0; i_layer < 6; i_layer++)
            {
                vec_PL_TRD_det_2D[i_sector][i_layer] = new TPolyLine();
            }
        }
        for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
        {
            Int_t i_sector = (Int_t)(TRD_detector/30);
            Int_t i_stack  = (Int_t)(TRD_detector%30/6);
            Int_t i_layer  = TRD_detector%6;
            //Int_t i_det = layer + 6*stack + 30*sector;


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

                if(i_stack == 0 && i_vertex < 4)
                {
                    vec_PL_TRD_det_2D[i_sector][i_layer] ->SetPoint(i_vertex,arr_pos_glb[0],arr_pos_glb[1]);
                    if(i_vertex == 0) vec_PL_TRD_det_2D[i_sector][i_layer] ->SetPoint(4,arr_pos_glb[0],arr_pos_glb[1]);
                }
            }

            gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
        }
        gEve->Redraw3D(kTRUE);
        //--------------------------

        //-------------------------
        Double_t TRD_layer_radii[6][2] =
        {
            {297.5,306.5},
            {310.0,320.0},
            {323.0,333.0},
            {336.0,345.5},
            {348.0,357.0},
            {361.0,371.0}
        };

        //Int_t color_layer_match[6] = {kRed,kGreen,kSpring-6,kYellow,kPink-3,kOrange+8};
        Int_t color_layer_match[6] = {kGray+1,kAzure-2,kGreen+2,kCyan+2,kOrange+2,kRed};

        //--------------------------
    }

    #endif
}

//----------------------------------------------------------------------------------------
void Ali_TRD_physics_analysis::Init_tree(TString SEList) // read data from files and create containers
{
    printf("Ali_TRD_physics_analysis::Init_tree \n");
    TString pinputdir = input_dir;

    TString in_list_name = SEList;
    SEList = input_dir_lists + SEList;

    Kalman_Track_photon  			= new Ali_Kalman_Track();
    TPC_Track_photon  	   			= new Ali_TPC_Track();
    Kalman_Track_interact  			= new Ali_Kalman_Track();
    TPC_Track_interact 	   			= new Ali_TPC_Track();
    TRD_Photon     			= new Ali_TRD_Photon();
    TRD_Nuclear_interaction = new Ali_TRD_Nuclear_interaction();
    TRD_Self_Event 			= new Ali_TRD_Self_Event();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( TRD_Self_TREE.Data(), TRD_Self_TREE.Data() );
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
                    input_SE ->AddFile(addfile.Data(),-1, TRD_Self_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( TRD_Self_BRANCH, &TRD_Self_Event );
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
Int_t Ali_TRD_physics_analysis::Loop_event(Long64_t i_event, Double_t dist_max, Int_t graphics, Int_t ME, Int_t TRD_photon) //get all info from each event
{
    //printf("Ali_TRD_ST_Analyze::Loop_event \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_Self_Event class definition)
    UShort_t NumPhotons         = TRD_Self_Event ->getNumPhotons(); // number of tracks in this event
    Int_t    NumNucInteractions = TRD_Self_Event ->getNumNucInteractions();
    EventVertexX         = TRD_Self_Event ->getx();
    EventVertexY         = TRD_Self_Event ->gety();
    EventVertexZ         = TRD_Self_Event ->getz();
    Global_Event         = i_event;
    Global_RunID         = TRD_Self_Event ->getid();
    TV3_EventVertex.SetXYZ(EventVertexX,EventVertexY,EventVertexZ);
    Float_t  V0MEq                = TRD_Self_Event ->getcent_class_V0MEq();

    TVector3 TV3_prim_vertex(EventVertexX,EventVertexY,EventVertexZ);

    //--------------------------------------------------

    //printf(" \n");
    if (i_event==950) printf("\n Event: %lld, Photon conversions: %d, Nuclear interactions: %d \n",i_event,NumPhotons,NumNucInteractions);

    //--------------------------------------------------
    // Photon loop

    vec_PhotonVertex.clear();

    vec_TLV_photon.clear();
    vec_TLV_photon_daughters.clear();
    vec_TLV_photon_daughters.resize(2);
    vec_opening_angle_photon.clear();
    vec_nsigma_electron.clear();

    vec_photon_kalman_chi2.clear();
    vec_photon_kalman_chi2.resize(NumPhotons);

    vec_photon_kalman_helices.clear();
    vec_photon_kalman_helices.resize(NumPhotons);

    vec_photon_tpc_helices.clear();
    vec_photon_tpc_helices.resize(NumPhotons);

    Int_t i_vertex_photon = 0;

    vector<TLorentzVector> vec_TLV_photons;
    vec_TLV_photons.clear();
    vector<TLorentzVector> vec_TLV_photons_cutA;
    vec_TLV_photons_cutA.clear();
    vector<TLorentzVector> vec_TLV_photons_cutB;
    vec_TLV_photons_cutB.clear();

    for(Int_t i_photon = 0; i_photon < NumPhotons; i_photon++)
    {
        TRD_Photon = TRD_Self_Event ->getPhoton(i_photon);

        Float_t PhotonVertexX   = TRD_Photon ->get_vertex_point(0);
        Float_t PhotonVertexY   = TRD_Photon ->get_vertex_point(1);
        Float_t PhotonVertexZ   = TRD_Photon ->get_vertex_point(2);

        Double_t Photon_radius = TMath::Sqrt(PhotonVertexX*PhotonVertexX + PhotonVertexY*PhotonVertexY);

        TV3_PhotonVertex.SetXYZ(PhotonVertexX,PhotonVertexY,PhotonVertexZ);

        TVector3 TV3_dir_prim_to_secondary = TV3_PhotonVertex - TV3_EventVertex;
        TV3_dir_prim_to_secondary *= 1.0/TV3_dir_prim_to_secondary.Mag();

        vec_PhotonVertex.push_back(TV3_PhotonVertex);

        Float_t ph_TRD_layer_shared       = TRD_Photon ->get_bit_TRD_layer_shared(); // -1 for TPC track, for TRD it stores the information about independent and shared TRD layers between the two electrons

        //printf("ph_TRD_layer_shared: %4.3f \n",ph_TRD_layer_shared);

        if (TRD_photon)
        {
            if(ph_TRD_layer_shared < 0) continue; //TRD photons
        }

        else
        {
            if(ph_TRD_layer_shared > 0) continue; //TPC photons
        }

    if (i_event==950 || i_event==3592 || i_event==3447 || i_event==3305 || //cut
    i_event==3270 || i_event==3115 || i_event==2930 || i_event==2904 || //cut
     i_event==3498 || i_event==3304 || i_event==3298 || i_event==2904 || i_event==2904 || i_event==3270 //not cut 
     || i_event==2930 || i_event==2373) printf("Event: %lld, Photon conversions: %d \n",i_event,NumPhotons); //not cut 



        Float_t ph_pT_AB_raw 	          = TRD_Photon ->get_pT_AB();
        Float_t ph_AP_pT   		  = TRD_Photon ->get_AP_pT();
        Float_t ph_AP_alpha   		  = TRD_Photon ->get_AP_alpha();
        Float_t ph_dca_min_xy            = TRD_Photon ->get_dca_min_xy();
        Float_t ph_dca_min_z   		  = TRD_Photon ->get_dca_min_z();
        Float_t ph_path_min   		  = TRD_Photon ->get_path_min();
        Float_t ph_Inv_mass_AB   	  = TRD_Photon ->get_Inv_mass_AB();
        Float_t ph_Eta_AB  		  = TRD_Photon ->get_Eta_AB();
        Float_t ph_Phi_AB   		  = TRD_Photon ->get_Phi_AB();
        Float_t ph_dot_product_dir_vertex = TRD_Photon ->get_dot_product_dir_vertex();
        Float_t ph_Inv_mass_AB_K0s   	  = TRD_Photon ->get_Inv_mass_AB_K0s();
        Float_t ph_dcaAB_xy          = TRD_Photon ->get_dcaAB_xy();
        Float_t ph_dcaAB_z   		  = TRD_Photon ->get_dcaAB_z();
        Float_t ph_Inv_mass_AB_Lambda     = TRD_Photon ->get_Inv_mass_AB_Lambda();
        Float_t ph_Inv_mass_AB_antiLambda = TRD_Photon ->get_Inv_mass_AB_antiLambda();
        TLorentzVector TLV_part_A         = TRD_Photon ->get_TLV_part_A();
        TLorentzVector TLV_part_B         = TRD_Photon ->get_TLV_part_B();
        UShort_t N_TPC_tracks             = TRD_Photon ->getNumTPC_Tracks(); //should be always 2 or 0
        vector<Double_t> nsigma_e_TPC;
        nsigma_e_TPC.resize(2);

        TLorentzVector TLV_photon_raw = TLV_part_A + TLV_part_B; // TRD photon
        vec_TLV_photons.push_back(TLV_photon_raw);

        for (Int_t i_track_TPC = 0; i_track_TPC < N_TPC_tracks; i_track_TPC++)
        {
            TPC_Track_photon = TRD_Photon ->getTPC_Track(i_track_TPC);
            nsigma_e_TPC[i_track_TPC] = TPC_Track_photon ->getnsigma_e_TPC();

            //printf("i_track_TPC: %d, nsigma_e_TPC: %4.3f \n",i_track_TPC,nsigma_e_TPC[i_track_TPC]);
        }

        //if(ph_TRD_layer_shared < 0 && (fabs(nsigma_e_TPC[0]) > 2.5  || fabs(nsigma_e_TPC[1]) > 2.5 )) continue;

        //printf("i_event: %lld, pTA/B: {%4.3f, %4.3f}, ph_dcaAB: %4.3f, ph_TRD_layer_shared: %4.3f \n",i_event,TLV_part_A.Pt(),TLV_part_B.Pt(),ph_dcaAB,ph_TRD_layer_shared);

        if(ph_dcaAB_xy > 2.0 && ph_dcaAB_xy > 15.0) continue;
        vec_TLV_photons_cutA.push_back(TLV_photon_raw);

        printf("original ph_pT_AB: %4.3f \n",ph_pT_AB_raw);

        TVector3 param;

        Float_t pT_A_corr = fabs(1.0*correct_pT_according_to_correlation_plot(-TLV_part_A.Pt(),par_pT_corr_neg));
        Float_t pT_B_corr = fabs(1.0*correct_pT_according_to_correlation_plot(-TLV_part_B.Pt(),par_pT_corr_neg));


        TVector3 vec_dir_A = TLV_part_A.Vect();
        TVector3 vec_dir_B = TLV_part_B.Vect();
        Double_t angle_between = vec_dir_A.Angle(vec_dir_B)*TMath::RadToDeg();
        //th1d_angle_between ->Fill(angle_between);

        //if(angle_between > 5.0) continue;


        //printf("pT_A: {%4.3f, %4.3f} \n",TLV_part_A.Pt(),pT_A_corr);

        TLorentzVector TLV_part_A_corr;
        TLorentzVector TLV_part_B_corr;

        //TLV_part_A_corr.SetXYZM(TV3_dir_prim_to_secondary.X()*pT_A_corr,TV3_dir_prim_to_secondary.Y()*pT_A_corr,TV3_dir_prim_to_secondary.Z()*pT_A_corr,0.0);
        //TLV_part_A_corr.SetXYZM(TV3_dir_prim_to_secondary.X()*pT_B_corr,TV3_dir_prim_to_secondary.Y()*pT_B_corr,TV3_dir_prim_to_secondary.Z()*pT_B_corr,0.0);

        TLV_part_A_corr.SetPtEtaPhiM(pT_A_corr,TLV_part_A.Eta(),TLV_part_A.Phi(),TLV_part_A.M());
        TLV_part_B_corr.SetPtEtaPhiM(pT_B_corr,TLV_part_B.Eta(),TLV_part_B.Phi(),TLV_part_B.M());
        //TLV_part_A_corr.SetPtEtaPhiM(pT_A_corr,TLV_part_A.Eta(),TLV_part_A.Phi(),0.0);
        //TLV_part_B_corr.SetPtEtaPhiM(pT_B_corr,TLV_part_B.Eta(),TLV_part_B.Phi(),0.0);


        TLorentzVector TLV_photon = TLV_part_A_corr + TLV_part_B_corr; // TRD photon
        //TLorentzVector TLV_photon = TLV_part_A + TLV_part_B; // TRD photon
        if(ph_TRD_layer_shared < 0) TLV_photon = TLV_part_A + TLV_part_B; // TPC photon


        //------------------------------------------
        // Armenteros-Podolanski cuts
        // AP cuts from https://www.physi.uni-heidelberg.de//Publications/PhDThesis_Leardini.pdf (eqn. 5.2)
        Double_t AP_alpha_max = 0.75; // 0.95
        Double_t AP_qT_max    = 0.03; // 0.05
        Double_t AP_cut_value = 1.0;
        Double_t AP_value = TMath::Power(ph_AP_alpha/AP_alpha_max,2.0) + TMath::Power(ph_AP_pT/AP_qT_max,2.0);

        //if(fabs(AP_alpha) < 0.2 && AP_pT > 0.0 && AP_pT < 0.02 && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.5 && pTB < 0.5 && dot_product_dir_vertex > 0.9) // TRD photon conversion
        //printf("AP_value: %4.3f, dot_product_dir_vertex: %4.3f, CA*CB: %4.3f, pTA: %4.3f, pTB: %4.3f \n",AP_value,dot_product_dir_vertex,CA*CB,pTA,pTB);
        //if(AP_value < AP_cut_value);
        //------------------------------------------


        //printf("new ph_pT_AB: %4.3f \n",ph_pT_AB);

        //------------------------------------------
        // Analyze the bit map
        Int_t independent_layer_A[6] = {0};
        Int_t independent_layer_B[6] = {0};
        Int_t shared_layer[6]        = {0};
        for(Int_t i_bit = 0; i_bit < 18; i_bit++)
        {
            if(((Int_t)ph_TRD_layer_shared >> i_bit) & 1)
            {
                if(i_bit < 6)                independent_layer_A[i_bit]   = 1;
                if(i_bit >= 6 && i_bit < 12) independent_layer_B[i_bit-6] = 1;
                if(i_bit >= 12)              shared_layer[i_bit-12]       = 1;
            }
        }
        //------------------------------------------

        vector<Bool_t> conds;
        Bool_t cond1=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 1 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3]) > 1 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3]) > 1 );
        conds.push_back(cond1);
        Bool_t cond2=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2] + shared_layer[3]) > 2 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3]) > 0 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3]) > 0 );
        conds.push_back(cond2);

        Bool_t cond3=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 0 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3] + independent_layer_A[2]) > 2 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3] + independent_layer_B[2]) > 2 );
        conds.push_back(cond3);

        Bool_t cond4=(
                      (shared_layer[0] + shared_layer[1]) > 0 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3] + independent_layer_A[2] + independent_layer_A[1]) > 2 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3] + independent_layer_B[2] + independent_layer_B[1]) > 2 );
        conds.push_back(cond4);

        Bool_t cond5=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2]) == 0 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3] + independent_layer_A[2] + independent_layer_A[1] + independent_layer_A[0]) > 3 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3] + independent_layer_B[2] + independent_layer_B[1] + independent_layer_B[0]) > 3 );
        conds.push_back(cond5);

        Bool_t cond6=!(
                       (shared_layer[3] + shared_layer[4] + shared_layer[5]) > 0 &&
                       (independent_layer_A[0] + independent_layer_A[1] + independent_layer_A[2] ) > 0 &&
                       (independent_layer_B[0] + independent_layer_B[1] + independent_layer_B[2] ) > 0 );
        conds.push_back(cond6);
printf("cut 1 \n");

        //if(
        //   ((cond1 || cond2 || cond3 ||cond4 || cond5) && cond6) // TRD photon
        //   || ph_TRD_layer_shared < 0.0 // TPC photon
        //  )
        {
            vec_TLV_photons_cutB.push_back(TLV_photon_raw);

            if(ph_dca_min_xy > 10.0 || (ph_dca_min_xy <= 10.0 && ph_path_min < 0.0)) // no close by TPC track
            {
                

printf("cut 2 \n");

                TVector3 TV3_dir;
                TV3_dir = TLV_photon.Vect();

                TVector3 TV3_sec_vertex(PhotonVertexX,PhotonVertexY,PhotonVertexZ);
                TH2_vertex_photon_XY     ->Fill(PhotonVertexX,PhotonVertexY);
                if(AP_value < 1.0) vec_TH1_vertex_photon_radius[0] ->Fill(Photon_radius);
                if(AP_value < 0.9) vec_TH1_vertex_photon_radius[1] ->Fill(Photon_radius);
                if(AP_value < 0.8) vec_TH1_vertex_photon_radius[2] ->Fill(Photon_radius);
                if(AP_value < 0.6) vec_TH1_vertex_photon_radius[3] ->Fill(Photon_radius);
                if(AP_value < 0.4) vec_TH1_vertex_photon_radius[4] ->Fill(Photon_radius);

                //TVector3 point;
                //point.SetXYZ(0.0,0.0,0.0);

                Double_t dist = calculateMinimumDistanceStraightToPoint(TV3_sec_vertex, TV3_dir, TV3_prim_vertex);

                if(dist < dist_max)
                {
printf("cut 3 \n");


                    UShort_t N_Kalman_tracks = TRD_Photon ->getNumKalman_Tracks(); //should be always 2 or 0

                    if(N_Kalman_tracks > 0)
                    {
                        for(Int_t i_track = 0; i_track < N_Kalman_tracks; i_track++)
                        {
                            Kalman_Track_photon = TRD_Photon ->getKalman_Track(i_track);

                            //Double_t Chi2 = Kalman_Track_photon ->get_Chi2();
                            //vec_photon_kalman_chi2[i_photon].push_back(Chi2);

                            vec_photon_kalman_helices[i_photon].push_back(new Ali_Helix_copy());

                            vec_photon_kalman_helices[i_photon][(Int_t)vec_photon_kalman_helices[i_photon].size()-1] ->setHelix(Kalman_Track_photon->getKalmanHelix_param(0),
                                                                                                                                Kalman_Track_photon->getKalmanHelix_param(1),Kalman_Track_photon->getKalmanHelix_param(2),
                                                                                                                                Kalman_Track_photon->getKalmanHelix_param(3),Kalman_Track_photon->getKalmanHelix_param(4),Kalman_Track_photon->getKalmanHelix_param(5));

                        }


                    }


                    if(N_TPC_tracks > 0)
                    {
                        for (Int_t i_track = 0; i_track < N_TPC_tracks; i_track++)
                        {
                            TPC_Track_photon = TRD_Photon ->getTPC_Track(i_track);

                            vec_photon_tpc_helices[i_photon].push_back(new Ali_Helix_copy());
                            vec_photon_tpc_helices[i_photon][(Int_t)vec_photon_tpc_helices[i_photon].size()-1] ->setHelix(TPC_Track_photon->getHelix_param(0),
                                                                                                                          TPC_Track_photon->getHelix_param(1),TPC_Track_photon->getHelix_param(2),
                                                                                                                          TPC_Track_photon->getHelix_param(3),TPC_Track_photon->getHelix_param(4),TPC_Track_photon->getHelix_param(5));
                        }
                    }
                    vec_TLV_photon.push_back(TLV_photon);
                    vec_TLV_photon_daughters[0].push_back(TLV_part_A);
                    vec_TLV_photon_daughters[1].push_back(TLV_part_B);

                    if(ME) vec_TLV_photon_mixed_events.push_back(TLV_photon); //add all photons from all events

                    vec_opening_angle_photon.push_back(angle_between);
                    vec_nsigma_electron.push_back(nsigma_e_TPC);

                    if(ME) vec_nsigma_electron_mixed_events.push_back(nsigma_e_TPC);
                    printf("vec_TLV_photon size: %d \n",(Int_t)vec_TLV_photon.size());

                    //#if defined (USEEVE)
                    if(graphics)
                    {
                        Double_t distance_to_prim_vertex = (TV3_prim_vertex - TV3_sec_vertex).Mag();

                        TEveLine_mother.resize(i_vertex_photon+1);
                        TEveLine_mother[i_vertex_photon] = new TEveLine();

                        TEveLine_mother[i_vertex_photon] ->SetNextPoint(PhotonVertexX,PhotonVertexY,PhotonVertexZ);
                        TEveLine_mother[i_vertex_photon] ->SetNextPoint(PhotonVertexX - distance_to_prim_vertex*TV3_dir[0],
                                                                        PhotonVertexY - distance_to_prim_vertex*TV3_dir[1],PhotonVertexZ - distance_to_prim_vertex*TV3_dir[2]);

                        //printf("photon vertex: {%4.3f, %4.3f, %4.3f}, second point: {%4.3f, %4.3f, %4.3f} \n",PhotonVertexX,PhotonVertexY,PhotonVertexZ,
                        //PhotonVertexX - distance_to_prim_vertex*TV3_dir[0],PhotonVertexY - distance_to_prim_vertex*TV3_dir[1],PhotonVertexZ - distance_to_prim_vertex*TV3_dir[2]);
                        i_vertex_photon++;
                        //printf("i_vertex_photon: %d \n",i_vertex_photon);

                    }
                    //#endif
                }
            }
        }

        // Add some condition if its a good photon
        // Fill TLorentzVector with photon.

        // AP cuts from https://www.physi.uni-heidelberg.de//Publications/PhDThesis_Leardini.pdf (eqn. 5.2)
        //Double_t AP_alpha_max = 0.95;
        //Double_t AP_qT_max    = 0.05;
        //Double_t AP_cut_value = 1.0;
        //Double_t AP_value = TMath::Power(AP_alpha/AP_alpha_max,2.0) + TMath::Power(AP_pT/AP_qT_max,2.0);

        //if(AP_value < AP_cut_value && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.8 && pTB < 0.8 && dot_product_dir_vertex > 0.9)
        // vec_TLV_photon ->SetPtEtaPhiM(ph_pT_AB);

        #if 0 //not needed right now



        #endif

        //printf("--> Photon conversion number %d: TRD tracks: %d, TPC tracks: %d \n",i_photon,N_Kalman_tracks,N_TPC_tracks);

    }
    //Photon loop done

    //printf(" \n");
    //printf("--------------------- \n");

    for(Int_t i_photoA = 0; i_photoA < (Int_t)vec_TLV_photons.size(); i_photoA++)
    {
        for(Int_t i_photoB = (i_photoA+1); i_photoB < (Int_t)vec_TLV_photons.size(); i_photoB++)
        {
            Double_t Angle_AB = vec_TLV_photons[i_photoA].Angle(vec_TLV_photons[i_photoB].Vect());
            h_angle_AB_raw ->Fill(Angle_AB*TMath::RadToDeg());

            //printf("i_photoA: %d, i_photoB: %d, angleAB: %4.3f, pAB: {%4.3f, %4.3f} \n",i_photoA, i_photoB, Angle_AB*TMath::RadToDeg(), vec_TLV_photons[i_photoA].P(), vec_TLV_photons[i_photoB].P());
        }
    }
    for(Int_t i_photoA = 0; i_photoA < (Int_t)vec_TLV_photons_cutA.size(); i_photoA++)
    {
        for(Int_t i_photoB = (i_photoA+1); i_photoB < (Int_t)vec_TLV_photons_cutA.size(); i_photoB++)
        {
            Double_t Angle_AB = vec_TLV_photons_cutA[i_photoA].Angle(vec_TLV_photons_cutA[i_photoB].Vect());
            h_angle_AB_raw_cutA ->Fill(Angle_AB*TMath::RadToDeg());
        }
    }
    for(Int_t i_photoA = 0; i_photoA < (Int_t)vec_TLV_photons_cutB.size(); i_photoA++)
    {
        for(Int_t i_photoB = (i_photoA+1); i_photoB < (Int_t)vec_TLV_photons_cutB.size(); i_photoB++)
        {
            Double_t Angle_AB = vec_TLV_photons_cutB[i_photoA].Angle(vec_TLV_photons_cutB[i_photoB].Vect());
            h_angle_AB_raw_cutB ->Fill(Angle_AB*TMath::RadToDeg());
        }
    }

    //printf("--------------------- \n");
    //printf(" \n");

    // Combine photons from vec_TLV_photon and calculate pi0
    // 3D straight line to primary vertex
    // calculateMinimumDistanceStraightToPoint
    // TLV_pi0 = vec_TLV_photon[0] + vec_TLV_photon[1];
    // pi_hist ->Fill(TLV_pi0.M());
    //--------------------------------------------------

    //--------------------------------------------------
    //Nuclear interactions loop

    #if 0 //also not needed right now

    vec_NIVertex.clear();

    vec_ni_kalman_chi2.clear();
    vec_ni_kalman_chi2.resize(NumNucInteractions);

    vec_ni_kalman_helices.clear();
    vec_ni_kalman_helices.resize(NumNucInteractions);

    vec_ni_tpc_helices.clear();
    vec_ni_tpc_helices.resize(NumPhotons);

    for(Int_t i_interaction = 0; i_interaction < NumNucInteractions; i_interaction++)
    {
        TRD_Nuclear_interaction = TRD_Self_Event ->getNucInteraction(i_interaction);

        Float_t NIVertexX   = TRD_Nuclear_interaction ->get_avg_sec_vertex(0);
        Float_t NIVertexY   = TRD_Nuclear_interaction ->get_avg_sec_vertex(1);
        Float_t NIVertexZ   = TRD_Nuclear_interaction ->get_avg_sec_vertex(2);

        TV3_NIVertex.SetXYZ(NIVertexX,NIVertexY,NIVertexZ);

        vec_NIVertex.push_back(TV3_NIVertex);

        Float_t N_close_vertex   = TRD_Nuclear_interaction ->get_N_close_vertex();
        Float_t ni_dcaAB_min 				  = TRD_Nuclear_interaction ->get_dcaAB_min();
        Float_t ni_TOFsignal_min   				  = TRD_Nuclear_interaction ->get_TOFsignal_min();
        Float_t ni_Track_length_min   			  = TRD_Nuclear_interaction ->get_Track_length_min();
        Float_t ni_TPCdEdx_min   			  = TRD_Nuclear_interaction ->get_TPCdEdx_min();
        Float_t ni_pT_min   		      = TRD_Nuclear_interaction ->get_pT_min();
        Float_t ni_momentum_min   		  = TRD_Nuclear_interaction ->get_momentum_min();
        Float_t ni_nuclev_bitmap  			      = TRD_Nuclear_interaction ->get_nuclev_bitmap();

        UShort_t N_Kalman_tracks = TRD_Nuclear_interaction ->getNumKalman_Tracks();

        if (N_Kalman_tracks > 0)
        {
            for (Int_t i_track = 0; i_track < N_Kalman_tracks; i_track++)
            {
                Kalman_Track_interact = TRD_Nuclear_interaction ->getKalman_Track(i_track);

                Double_t Chi2 = Kalman_Track_interact ->get_Chi2();
                vec_ni_kalman_chi2[i_interaction].push_back(Chi2);

                vec_ni_kalman_helices[i_interaction].push_back(new Ali_Helix_copy());

                vec_ni_kalman_helices[i_interaction][(Int_t)vec_ni_kalman_helices[i_interaction].size()-1] ->setHelix(Kalman_Track_interact->getKalmanHelix_param(0),
                                                                                                                      Kalman_Track_interact->getKalmanHelix_param(1),Kalman_Track_interact->getKalmanHelix_param(2),
                                                                                                                      Kalman_Track_interact->getKalmanHelix_param(3),Kalman_Track_interact->getKalmanHelix_param(4),Kalman_Track_interact->getKalmanHelix_param(5));
            }
        }

        UShort_t N_TPC_tracks = TRD_Nuclear_interaction ->getNumTPC_Tracks();

        if (N_TPC_tracks > 0)
        {
            for (Int_t i_track = 0; i_track < N_TPC_tracks; i_track++)
            {
                TPC_Track_interact = TRD_Nuclear_interaction ->getTPC_Track(i_track);

                vec_ni_tpc_helices[i_interaction].push_back(new Ali_Helix_copy());
                vec_ni_tpc_helices[i_interaction][(Int_t)vec_ni_tpc_helices[i_interaction].size()-1] ->setHelix(TPC_Track_interact->getHelix_param(0),
                                                                                                                TPC_Track_interact->getHelix_param(1),TPC_Track_interact->getHelix_param(2),
                                                                                                                TPC_Track_interact->getHelix_param(3),TPC_Track_interact->getHelix_param(4),TPC_Track_interact->getHelix_param(5));
            }
        }

        //printf("--> Nuclear interaction number %d: TRD tracks: %d, TPC tracks: %d \n",i_interaction,N_Kalman_tracks,N_TPC_tracks);

    }
    #endif
    // Nuclear interactions done
    //--------------------------------------------------

    //#if defined(USEEVE)

    if (graphics)
    {

        for(Int_t i_mother = 0; i_mother < (Int_t)TEveLine_mother.size(); i_mother++)
        {
            TEveLine_mother[i_mother]    ->SetLineStyle(9);
            TEveLine_mother[i_mother]    ->SetLineWidth(3);
            TEveLine_mother[i_mother]    ->SetMainColor(kMagenta);
            TEveLine_mother[i_mother]    ->SetMainAlpha(1.0);
            gEve->AddElement(TEveLine_mother[i_mother]);
        }
        gEve->Redraw3D(kTRUE);
        //#endif
    }


    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_physics_analysis::Calculate_pi0_mass_SE() //loop over all good photons
{
    //printf("Calculate_pi0_mass_SE started \n");

    for(Int_t i_photon_A = 0; i_photon_A < (Int_t)vec_TLV_photon.size(); i_photon_A++)
    {
        Double_t ptA1 = vec_TLV_photon_daughters[0][i_photon_A].Pt();
        Double_t ptA2 = vec_TLV_photon_daughters[1][i_photon_A].Pt();

        for(Int_t i_photon_B = (i_photon_A + 1); i_photon_B < (Int_t)vec_TLV_photon.size(); i_photon_B++)
        {
            Double_t ptB1 = vec_TLV_photon_daughters[0][i_photon_B].Pt();
            Double_t ptB2 = vec_TLV_photon_daughters[1][i_photon_B].Pt();

            printf("ptA1: %4.5f, ptA2: %4.5f, ptB1: %4.5f, ptB2: %4.5f \n",ptA1,ptA2,ptB1,ptB2);
            
            if(fabs(ptA1 - ptB1) < 0.001 || fabs(ptA1 - ptB2) < 0.001 || fabs(ptA2 - ptB1) < 0.001 || fabs(ptA2 - ptB2) < 0.001)
            {
                //printf("ptA1: %4.5f, ptA2: %4.5f, ptB1: %4.5f, ptB2: %4.5f \n",ptA1,ptA2,ptB1,ptB2);
                continue;
            }

            TLorentzVector TLV_pi0;

            TLV_pi0 = vec_TLV_photon[i_photon_A] + vec_TLV_photon[i_photon_B];
            Double_t pT_pi0 = TLV_pi0.Pt();
            Double_t Mass_pi0 = TLV_pi0.M();
            vec_TH1_mass_pi0[0] ->Fill(Mass_pi0);

            Double_t angle_between_photons = TMath::RadToDeg()*vec_TLV_photon[i_photon_A].Angle(vec_TLV_photon[i_photon_B].Vect());

            TH1_angle_between_photons      ->Fill(angle_between_photons);
            TH2D_angle_photons_vs_inv_mass ->Fill(Mass_pi0,angle_between_photons);

            if(fabs(vec_nsigma_electron[i_photon_A][0]) < 2.5 && fabs(vec_nsigma_electron[i_photon_A][1]) < 2.5)
            {
                vec_TH1_mass_pi0[1] ->Fill(Mass_pi0);
                if(fabs(vec_nsigma_electron[i_photon_B][0]) < 2.5 && fabs(vec_nsigma_electron[i_photon_B][1]) < 2.5)
                {
                    vec_TH1_mass_pi0[2] ->Fill(Mass_pi0);
                }
            }

            if(fabs(vec_nsigma_electron[i_photon_A][0]) > 2.5 && fabs(vec_nsigma_electron[i_photon_A][1]) > 2.5 &&
               fabs(vec_nsigma_electron[i_photon_B][0]) > 2.5 && fabs(vec_nsigma_electron[i_photon_B][1]) > 2.5
              )
            {
                vec_TH1_mass_pi0[3] ->Fill(Mass_pi0);
            }

            //if(pT_pi0 > 0.2 && pT_pi0 <= 0.5) vec_TH1_mass_pi0[1] ->Fill(Mass_pi0);
            //if(pT_pi0 > 0.5 && pT_pi0 <= 1.0) vec_TH1_mass_pi0[2] ->Fill(Mass_pi0);
            //if(pT_pi0 > 1.0 && pT_pi0 <= 1.5) vec_TH1_mass_pi0[3] ->Fill(Mass_pi0);
            //if(pT_pi0 > 1.5 && pT_pi0 <= 2.0) vec_TH1_mass_pi0[4] ->Fill(Mass_pi0);

            if(Mass_pi0 < 0.16)
            {
                th1d_angle_between->Fill(vec_opening_angle_photon[i_photon_A]);
                th1d_angle_between->Fill(vec_opening_angle_photon[i_photon_B]);
            }
            //printf("number of entries: %d \n",TH1_mass_pi0->GetEntries());
        }

    }

}
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void Ali_TRD_physics_analysis::Calculate_pi0_mass_SE_and_ME() //loop over all good photons
{
    //printf("Calculate_pi0_mass_SE_and_ME started \n");


    Int_t photon_max = (Int_t)vec_TLV_photon_mixed_events.size();

    printf("photon max: %d \n",photon_max);
    
    if (photon_max > 5000) photon_max = 5000;

    printf("photon max: %d \n",photon_max);


    for(Int_t i_photon_A = 0; i_photon_A < photon_max; i_photon_A++)
    {
        for(Int_t i_photon_B = (i_photon_A + 1); i_photon_B < photon_max; i_photon_B++)
        {
            TLorentzVector TLV_pi0;

            TLV_pi0 = vec_TLV_photon_mixed_events[i_photon_A] + vec_TLV_photon_mixed_events[i_photon_B];
            Double_t pT_pi0 = TLV_pi0.Pt();
            Double_t Mass_pi0 = TLV_pi0.M();
            vec_TH1_mass_pi0_ME[0] ->Fill(Mass_pi0);

            //Double_t angle_between_photons = TMath::RadToDeg()*vec_TLV_photon[i_photon_A].Angle(vec_TLV_photon[i_photon_B].Vect());

            //TH1_angle_between_photons      ->Fill(angle_between_photons);
            //TH2D_angle_photons_vs_inv_mass ->Fill(Mass_pi0,angle_between_photons);

            if(fabs(vec_nsigma_electron_mixed_events[i_photon_A][0]) < 2.5 && fabs(vec_nsigma_electron_mixed_events[i_photon_A][1]) < 2.5)
            {
                vec_TH1_mass_pi0_ME[1] ->Fill(Mass_pi0);
                if(fabs(vec_nsigma_electron_mixed_events[i_photon_B][0]) < 2.5 && fabs(vec_nsigma_electron_mixed_events[i_photon_B][1]) < 2.5)
                {
                    vec_TH1_mass_pi0_ME[2] ->Fill(Mass_pi0);
                }
            }

            if(fabs(vec_nsigma_electron_mixed_events[i_photon_A][0]) > 2.5 && fabs(vec_nsigma_electron_mixed_events[i_photon_A][1]) > 2.5 &&
               fabs(vec_nsigma_electron_mixed_events[i_photon_B][0]) > 2.5 && fabs(vec_nsigma_electron_mixed_events[i_photon_B][1]) > 2.5
              )
            {
                vec_TH1_mass_pi0_ME[3] ->Fill(Mass_pi0);
            }

            //if(pT_pi0 > 0.2 && pT_pi0 <= 0.5) vec_TH1_mass_pi0[1] ->Fill(Mass_pi0);
            //if(pT_pi0 > 0.5 && pT_pi0 <= 1.0) vec_TH1_mass_pi0[2] ->Fill(Mass_pi0);
            //if(pT_pi0 > 1.0 && pT_pi0 <= 1.5) vec_TH1_mass_pi0[3] ->Fill(Mass_pi0);
            //if(pT_pi0 > 1.5 && pT_pi0 <= 2.0) vec_TH1_mass_pi0[4] ->Fill(Mass_pi0);

            //if(Mass_pi0 < 0.16)
            //{
            //    th1d_angle_between->Fill(vec_opening_angle_photon[i_photon_A]);
            //    th1d_angle_between->Fill(vec_opening_angle_photon[i_photon_B]);
            //}
            //printf("number of entries: %d \n",TH1_mass_pi0->GetEntries());
        }

    }

}
//----------------------------------------------------------------------------------------

#if 1
//----------------------------------------------------------------------------------------
void Ali_TRD_physics_analysis::Draw(Int_t ME)
{
    printf("total: %d \n",(Int_t)vec_TH1_mass_pi0[0]->GetEntries());

    Double_t start_norm = 0.5;
    Double_t stop_norm  = 2.0;

    TCanvas* can_angle_between_photons = Draw_1D_histo_and_canvas(TH1_angle_between_photons,"can_angle_between_photons",800,650,0,0,""); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    TCanvas* can_angle_between         = Draw_1D_histo_and_canvas(th1d_angle_between,"can_angle_between",800,650,0,0,""); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option


    //Normalize hists with pi0 mass:

    for (Int_t i_hist = 0; i_hist < (Int_t)vec_TH1_mass_pi0.size(); i_hist++)
    {
        //vec_TH1_mass_pi0[i_hist] ->Scale(1.0/vec_TH1_mass_pi0[i_hist] ->GetEntries());
    }

    Double_t norm_SE = vec_TH1_mass_pi0[0]->Integral(vec_TH1_mass_pi0[0]->FindBin(start_norm),vec_TH1_mass_pi0[0]->FindBin(stop_norm));

    vec_TH1_mass_pi0[0] ->SetLineColor(kBlack);
    vec_TH1_mass_pi0[0] ->GetXaxis()->SetTitle("M_{#gamma,#gamma} (GeV/c^{2})");
    vec_TH1_mass_pi0[0] ->GetYaxis()->SetTitle("counts");
    TCanvas* can_TH1_mass_pi0 = Draw_1D_histo_and_canvas(vec_TH1_mass_pi0[0],"can_TH1_mass_pi0",800,650,0,0,""); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    can_TH1_mass_pi0 ->cd();

    vec_TH1_mass_pi0[1] ->SetLineColor(kRed);
    //vec_TH1_mass_pi0[1] ->DrawCopy("same");

    vec_TH1_mass_pi0[2] ->SetLineColor(kBlue); //2 TPC electons 
    //vec_TH1_mass_pi0[2] ->DrawCopy("same");

    vec_TH1_mass_pi0[3] ->SetLineColor(kOrange);
    //vec_TH1_mass_pi0[3] ->DrawCopy("same");

    //vec_TH1_mass_pi0[4] ->SetLineColor(kMagenta+2);
    //vec_TH1_mass_pi0[4] ->DrawCopy("same");



    if(ME)
    {
        Double_t norm_ME = vec_TH1_mass_pi0_ME[0]->Integral(vec_TH1_mass_pi0_ME[0]->FindBin(start_norm),vec_TH1_mass_pi0_ME[0]->FindBin(stop_norm));

        vec_TH1_mass_pi0_ME[0] ->Scale(norm_SE/norm_ME);
        vec_TH1_mass_pi0_ME[0] ->SetLineColor(kGray+2);
        vec_TH1_mass_pi0_ME[0] ->SetFillColor(kGray);
        vec_TH1_mass_pi0_ME[0] ->SetFillStyle(3001);
        vec_TH1_mass_pi0_ME[0] ->DrawCopy("same hist E");
    }

    TH1D* TH1_mass_pi0_diff = (TH1D*)vec_TH1_mass_pi0[0] ->Clone("TH1_mass_pi0_diff");
    TH1_mass_pi0_diff ->Add(vec_TH1_mass_pi0_ME[0],-1.0);
    TH1_mass_pi0_diff ->SetLineColor(kMagenta);
    TH1_mass_pi0_diff ->DrawCopy("same");


    //---------------------
    double x[55];
    double y[55];
    double dxlow[55];
    double dxup[55];
    double dylow[55];
    double dyup[55];
    int n = 0;
    x[n] = 0.0581395;	y[n] = 436.776;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.0930233;	y[n] = 635.768;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.116279;	y[n] = 771.788;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.156977;	y[n] = 754.156;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.180233;	y[n] = 696.222;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.209302;	y[n] = 701.259;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.25;	y[n] = 691.184;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.267442;	y[n] = 796.977;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.296512;	y[n] = 774.307;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.319767;	y[n] = 824.685;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.360465;	y[n] = 804.534;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.372093;	y[n] = 892.695;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.424419;	y[n] = 870.025;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.424419;	y[n] = 983.375;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.476744;	y[n] = 943.073;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.5;	y[n] = 1001.01;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.569767;	y[n] = 1086.65;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.593023;	y[n] = 1131.99;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.622093;	y[n] = 1202.52;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    Double_t y_val_scale_Martin = vec_TH1_mass_pi0[0] ->GetBinContent(vec_TH1_mass_pi0[0]->FindBin(x[n-1]));
    Double_t scale_factor_Martin = y_val_scale_Martin/y[n-1];

    x[n] = 0.680233;	y[n] = 1230.23;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.72093;	y[n] = 1262.97;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.796512;	y[n] = 1270.53;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.866279;	y[n] = 1257.93;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.918605;	y[n] = 1237.78;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.953488;	y[n] = 1184.89;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 0.982558;	y[n] = 1142.07;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.02326;	y[n] = 1104.28;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.05814;	y[n] = 1046.35;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.08721;	y[n] = 985.894;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.11628;	y[n] = 948.111;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.14535;	y[n] = 895.214;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.20349;	y[n] = 842.317;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.25;	y[n] = 784.383;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.29651;	y[n] = 734.005;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.35465;	y[n] = 688.665;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.4186;	y[n] = 638.287;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.47093;	y[n] = 580.353;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.52907;	y[n] = 524.937;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.56395;	y[n] = 487.154;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.61047;	y[n] = 456.927;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.68023;	y[n] = 426.7;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.76744;	y[n] = 383.879;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.82558;	y[n] = 346.096;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.87209;	y[n] = 323.426;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 1.93605;	y[n] = 288.161;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.02326;	y[n] = 260.453;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.11047;	y[n] = 240.302;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.20349;	y[n] = 210.076;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.30233;	y[n] = 179.849;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.4186;	y[n] = 157.179;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.53488;	y[n] = 139.547;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.63372;	y[n] = 124.433;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.73837;	y[n] = 104.282;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.85465;	y[n] = 94.2065;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;
    x[n] = 2.96512;	y[n] = 79.0932;	dxlow[n] = 0;	dxup[n] = 0;	dylow[n] = 0;	dyup[n] = 0;	n++;

    for(Int_t i_point = 0; i_point < 55; i_point++)
    {
        y[i_point] = y[i_point]*scale_factor_Martin;
        //printf("i_point: %d, point: %4.3f \n",i_point,y[i_point]);
    }

    TGraphAsymmErrors *xyscan = new TGraphAsymmErrors(n, x, y, dxlow, dxup, dylow, dyup);
    xyscan->SetMarkerStyle(20);
    xyscan->SetMarkerColor(1);
    xyscan->Draw("PE same");
    //---------------------

    TCanvas* can_TH2_vertex_photon_XY = Draw_2D_histo_and_canvas(TH2_vertex_photon_XY,"can_TH2_vertex_photon_XY",800,650,0,0,"COLZ"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    can_TH1_mass_pi0 ->SetLogz(1);

    TCanvas* can_vertex_photon_radius = Draw_1D_histo_and_canvas(vec_TH1_vertex_photon_radius[0],"can_vertex_photon_radius",800,650,0,0,""); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    can_vertex_photon_radius ->cd();
    vec_TH1_vertex_photon_radius[1] ->Draw("same");
    vec_TH1_vertex_photon_radius[2] ->Draw("same");
    vec_TH1_vertex_photon_radius[3] ->Draw("same");
    vec_TH1_vertex_photon_radius[4] ->Draw("same");



    TCanvas* can_angle_photons_vs_inv_mass = Draw_2D_histo_and_canvas(TH2D_angle_photons_vs_inv_mass,"can_angle_photons_vs_inv_mass",800,650,0,0,"COLZ"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    can_angle_photons_vs_inv_mass ->SetLogz(1);


    TCanvas* can_angle_AB_raw = Draw_1D_histo_and_canvas(h_angle_AB_raw,"can_angle_AB_raw",800,650,0,0,""); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    can_angle_AB_raw ->cd();
    h_angle_AB_raw_cutA ->SetLineColor(kRed);
    h_angle_AB_raw_cutA ->DrawCopy("same");
    h_angle_AB_raw_cutB ->SetLineColor(kGreen+2);
    h_angle_AB_raw_cutB ->DrawCopy("same");
}
//----------------------------------------------------------------------------------------
#endif
