
#include "TRD_physics_analysis.h"

// implementation of the main code

// Something like with your new classes


//----------------------------------------------------------------------------------------
Ali_TRD_physics_analysis::Ali_TRD_physics_analysis(TString out_dir, TString out_file_name)
{
	HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    printf("test printf: %s \n",HistName.Data());
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

    printf("\n Event: %lld, Photon conversions: %d, Nuclear interactions: %d \n",i_event,NumPhotons,NumNucInteractions);

    //--------------------------------------------------
    // Photon loop

    vec_PhotonVertex.clear();

    vec_photon_kalman_chi2.clear();
    vec_photon_kalman_chi2.resize(NumPhotons);

    vec_photon_kalman_helices.clear();
    vec_photon_kalman_helices.resize(NumPhotons);

    vec_photon_tpc_helices.clear();
    vec_photon_tpc_helices.resize(NumPhotons);

    for(Int_t i_photon = 0; i_photon < NumPhotons; i_photon++)
    {
        TRD_Photon = TRD_Self_Event ->getPhoton(i_photon);

        Float_t PhotonVertexX   = TRD_Photon ->get_vertex_point(0);
        Float_t PhotonVertexY   = TRD_Photon ->get_vertex_point(1);
        Float_t PhotonVertexZ   = TRD_Photon ->get_vertex_point(2);

        TV3_PhotonVertex.SetXYZ(PhotonVertexX,PhotonVertexY,PhotonVertexZ);

        vec_PhotonVertex.push_back(TV3_PhotonVertex);

        Float_t ph_TRD_layer_shared   = TRD_Photon ->get_bit_TRD_layer_shared();
        Float_t ph_pT_AB 				  = TRD_Photon ->get_pT_AB();
        Float_t ph_AP_pT   				  = TRD_Photon ->get_AP_pT();
        Float_t ph_AP_alpha   			  = TRD_Photon ->get_AP_alpha();
        Float_t ph_dca_min   			  = TRD_Photon ->get_dca_min();
        Float_t ph_path_min   		      = TRD_Photon ->get_path_min();
        Float_t ph_Inv_mass_AB   		  = TRD_Photon ->get_Inv_mass_AB();
        Float_t ph_Eta_AB  			      = TRD_Photon ->get_Eta_AB();
        Float_t ph_Phi_AB   		      = TRD_Photon ->get_Phi_AB();
        Float_t ph_dot_product_dir_vertex = TRD_Photon ->get_dot_product_dir_vertex();
        Float_t ph_Inv_mass_AB_K0s   	  = TRD_Photon ->get_Inv_mass_AB_K0s();
        Float_t ph_dcaAB   				  = TRD_Photon ->get_dcaAB();
        Float_t ph_Inv_mass_AB_Lambda     = TRD_Photon ->get_Inv_mass_AB_Lambda();
        Float_t ph_Inv_mass_AB_antiLambda = TRD_Photon ->get_Inv_mass_AB_antiLambda();

        UShort_t N_Kalman_tracks = TRD_Photon ->getNumKalman_Tracks(); //should be always 2 or 0

        if (N_Kalman_tracks > 0) 
        	{
        		for (Int_t i_track = 0; i_track < N_Kalman_tracks; i_track++)
		        {
		        	Kalman_Track_photon = TRD_Photon ->getKalman_Track(i_track);

		        	Double_t Chi2 = Kalman_Track_photon ->get_Chi2();
		        	vec_photon_kalman_chi2[i_photon].push_back(Chi2);

		        	vec_photon_kalman_helices[i_photon].push_back(new Ali_Helix());

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

		        	vec_photon_tpc_helices[i_photon].push_back(new Ali_Helix());
		        	vec_photon_tpc_helices[i_photon][(Int_t)vec_photon_tpc_helices[i_photon].size()-1] ->setHelix(TPC_Track_photon->getHelix_param(0),
		        		TPC_Track_photon->getHelix_param(1),TPC_Track_photon->getHelix_param(2),
		        		TPC_Track_photon->getHelix_param(3),TPC_Track_photon->getHelix_param(4),TPC_Track_photon->getHelix_param(5));
		        }
		    }

    printf("--> Photon conversion number %d: TRD tracks: %d, TPC tracks: %d \n",i_photon,N_Kalman_tracks,N_TPC_tracks);

    }
    //Photon loop done 
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

		        	vec_ni_kalman_helices[i_interaction].push_back(new Ali_Helix());

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

		        	vec_ni_tpc_helices[i_interaction].push_back(new Ali_Helix());
		        	vec_ni_tpc_helices[i_interaction][(Int_t)vec_ni_tpc_helices[i_interaction].size()-1] ->setHelix(TPC_Track_interact->getHelix_param(0),
		        		TPC_Track_interact->getHelix_param(1),TPC_Track_interact->getHelix_param(2),
		        		TPC_Track_interact->getHelix_param(3),TPC_Track_interact->getHelix_param(4),TPC_Track_interact->getHelix_param(5));
		        }
		    }

    printf("--> Nuclear interaction number %d: TRD tracks: %d, TPC tracks: %d \n",i_interaction,N_Kalman_tracks,N_TPC_tracks);

    }
 	// Nuclear interactions done
    //--------------------------------------------------


    return 1;
}
//----------------------------------------------------------------------------------------
