
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
Ali_TRD_physics_analysis::Ali_TRD_physics_analysis(TString out_dir, TString out_file_name)
{
    SetRootGraphicStyle();

    HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    printf("test printf: %s \n",HistName.Data());

    TH2_vertex_photon_XY = new TH2D("TH2_vertex_photon_XY","TH2_vertex_photon_XY",400,-400,400,400,-400,400);
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
Int_t Ali_TRD_physics_analysis::Loop_event(Long64_t i_event) //get all info from each event
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

    //--------------------------------------------------

    //printf("\n Event: %lld, Photon conversions: %d, Nuclear interactions: %d \n",i_event,NumPhotons,NumNucInteractions);

    //--------------------------------------------------
    // Photon loop

    vec_PhotonVertex.clear();

    vec_photon_kalman_chi2.clear();
    vec_photon_kalman_chi2.resize(NumPhotons);

    vec_photon_kalman_helices.clear();
    vec_photon_kalman_helices.resize(NumPhotons);

    vec_photon_tpc_helices.clear();
    vec_photon_tpc_helices.resize(NumPhotons);

    // vec_TLV_photon

    for(Int_t i_photon = 0; i_photon < NumPhotons; i_photon++)
    {
        TRD_Photon = TRD_Self_Event ->getPhoton(i_photon);

        Float_t PhotonVertexX   = TRD_Photon ->get_vertex_point(0);
        Float_t PhotonVertexY   = TRD_Photon ->get_vertex_point(1);
        Float_t PhotonVertexZ   = TRD_Photon ->get_vertex_point(2);

        TV3_PhotonVertex.SetXYZ(PhotonVertexX,PhotonVertexY,PhotonVertexZ);
        TH2_vertex_photon_XY ->Fill(PhotonVertexX,PhotonVertexY);


        vec_PhotonVertex.push_back(TV3_PhotonVertex);

        Float_t ph_TRD_layer_shared       = TRD_Photon ->get_bit_TRD_layer_shared(); // -1 for TPC track, for TRD it stores the information about independent and shared TRD layers between the two electrons
        Float_t ph_pT_AB 		  = TRD_Photon ->get_pT_AB();
        Float_t ph_AP_pT   		  = TRD_Photon ->get_AP_pT();
        Float_t ph_AP_alpha   		  = TRD_Photon ->get_AP_alpha();
        Float_t ph_dca_min   		  = TRD_Photon ->get_dca_min();
        Float_t ph_path_min   		  = TRD_Photon ->get_path_min();
        Float_t ph_Inv_mass_AB   	  = TRD_Photon ->get_Inv_mass_AB();
        Float_t ph_Eta_AB  		  = TRD_Photon ->get_Eta_AB();
        Float_t ph_Phi_AB   		  = TRD_Photon ->get_Phi_AB();
        Float_t ph_dot_product_dir_vertex = TRD_Photon ->get_dot_product_dir_vertex();
        Float_t ph_Inv_mass_AB_K0s   	  = TRD_Photon ->get_Inv_mass_AB_K0s();
        Float_t ph_dcaAB   		  = TRD_Photon ->get_dcaAB();
        Float_t ph_Inv_mass_AB_Lambda     = TRD_Photon ->get_Inv_mass_AB_Lambda();
        Float_t ph_Inv_mass_AB_antiLambda = TRD_Photon ->get_Inv_mass_AB_antiLambda();


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

        // Add some condition if its a good photon
        // Fill TLorentzVector with photon.

        // AP cuts from https://www.physi.uni-heidelberg.de//Publications/PhDThesis_Leardini.pdf (eqn. 5.2)
        //Double_t AP_alpha_max = 0.95;
        //Double_t AP_qT_max    = 0.05;
        //Double_t AP_cut_value = 1.0;
        //Double_t AP_value = TMath::Power(AP_alpha/AP_alpha_max,2.0) + TMath::Power(AP_pT/AP_qT_max,2.0);

        //if(AP_value < AP_cut_value && CA*CB < 0.0 && pTA > 0.04 && pTB > 0.04 && pTA < 0.8 && pTB < 0.8 && dot_product_dir_vertex > 0.9)
        // vec_TLV_photon ->SetPtEtaPhiM(ph_pT_AB);

        UShort_t N_Kalman_tracks = TRD_Photon ->getNumKalman_Tracks(); //should be always 2 or 0

        if (N_Kalman_tracks > 0)
        {
            for (Int_t i_track = 0; i_track < N_Kalman_tracks; i_track++)
            {
                Kalman_Track_photon = TRD_Photon ->getKalman_Track(i_track);

                Double_t Chi2 = Kalman_Track_photon ->get_Chi2();
                vec_photon_kalman_chi2[i_photon].push_back(Chi2);

                vec_photon_kalman_helices[i_photon].push_back(new Ali_Helix_copy());

                vec_photon_kalman_helices[i_photon][(Int_t)vec_photon_kalman_helices[i_photon].size()-1] ->setHelix(Kalman_Track_photon->getKalmanHelix_param(0),
                                                                                                                    Kalman_Track_photon->getKalmanHelix_param(1),Kalman_Track_photon->getKalmanHelix_param(2),
                                                                                                                    Kalman_Track_photon->getKalmanHelix_param(3),Kalman_Track_photon->getKalmanHelix_param(4),Kalman_Track_photon->getKalmanHelix_param(5));
            }
        }

        UShort_t N_TPC_tracks = TRD_Photon ->getNumTPC_Tracks(); //should be always 2 or 0

        if (N_TPC_tracks > 0)
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

        //printf("--> Photon conversion number %d: TRD tracks: %d, TPC tracks: %d \n",i_photon,N_Kalman_tracks,N_TPC_tracks);

    }
    //Photon loop done


    // Combine photons from vec_TLV_photon and calculate pi0
    // 3D straight line to primary vertex
    // calculateMinimumDistanceStraightToPoint
    // TLV_pi0 = vec_TLV_photon[0] + vec_TLV_photon[1];
    // pi_hist ->Fill(TLV_pi0.M());
    //--------------------------------------------------

    //--------------------------------------------------
    //Nuclear interactions loop

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
    // Nuclear interactions done
    //--------------------------------------------------


    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_physics_analysis::Draw()
{
    TCanvas* can_vertex_photon_XY = Draw_2D_histo_and_canvas(TH2_vertex_photon_XY,"can_vertex_photon_XY",800,650,0,0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size,Double_t min_val, Double_t max_val, TString option
    can_vertex_photon_XY ->SetLogz(1);
}
//----------------------------------------------------------------------------------------

