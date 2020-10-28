
#include "TTRD_ST_Make_Tracklets.h"
#include "TTRD_ST_Make_Tracklets_LinkDef.h"


ClassImp(TTRD_ST_Make_Tracklets);

//----------------------------------------------------------------------------------------
TTRD_ST_Make_Tracklets::TTRD_ST_Make_Tracklets()
{
    Init_QA();


    vec_TV3_Tracklet_pos.resize(540);
    vec_TV3_Tracklet_dir.resize(540);
    vec_h_diff_helix_line_impact_angle.resize(3); // [all,-,+]
    vec_h_diff_ref_loc_off_trkl.resize(3); // [all,-,+]
    for(Int_t i_charge = 0; i_charge < 3; i_charge++)
    {
        vec_h_diff_helix_line_impact_angle[i_charge].resize(547);
        vec_h_diff_ref_loc_off_trkl[i_charge].resize(547);
    }

    // Standard time bins
    vec_merge_time_bins.resize(24+1);
    for(Int_t i_time = 0; i_time < (24+1); i_time++)
    {
        vec_merge_time_bins[i_time] = i_time;
    }


    h_delta_angle_perp_impact = new TH1D("h_delta_angle_perp_impact","h_delta_angle_perp_impact",240,-30,30);
    h_detector_hit            = new TH1D("h_detector_hit","h_detector_hit",540,0,540);

    vec_layer_in_fit.resize(7);


    vec_digit_single_info.resize(14); // x,y,z,time,ADC,sector,stack,layer,row,column,dca,dca_x,dca_y,dca_z
    vec_track_single_info.resize(12); // dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p


    vec_tracklet_fit_points.resize(7); // layers 0-5, 6 = fit through first time cluster points of all layers
    for(Int_t i_layer = 0; i_layer < 7; i_layer++)
    {
        vec_tracklet_fit_points[i_layer].resize(2); // start, stop point
        for(Int_t i_start_stop = 0; i_start_stop < 2; i_start_stop++)
        {
            vec_tracklet_fit_points[i_layer][i_start_stop].resize(3); // x,y,z
        }
    }


    vec_TV3_local_pos.resize(8);

    for(Int_t i_charge = 0; i_charge < 3; i_charge++)
    {
        for(Int_t TRD_detector = 0; TRD_detector < 547; TRD_detector++)
        {
            HistName = "vec_h_diff_helix_line_impact_angle_c_";
            HistName += i_charge;
            HistName += "_det_";
            HistName += TRD_detector;
            vec_h_diff_helix_line_impact_angle[i_charge][TRD_detector] = new TH1D(HistName.Data(),HistName.Data(),100,-10.0,10.0);

            HistName = "vec_h_diff_ref_loc_off_trkl_";
            HistName += i_charge;
            HistName += "_det_";
            HistName += TRD_detector;
            vec_h_diff_ref_loc_off_trkl[i_charge][TRD_detector] = new TH1D(HistName.Data(),HistName.Data(),400,-2.0,2.0);
        }
    }

    // OK
    fGeo = new AliTRDgeometry;

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
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] = new TH1D(HistName.Data(),HistName.Data(),540,0,540);
        }

    }
    vec_TV3_TRD_center.resize(540);
    vec_TV3_TRD_center_offset.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_TRD_center[i_det].resize(3);
    }

    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        Int_t flag_not_installed_TRD = 0;
        for(Int_t i_not_installed = 0; i_not_installed < 19; i_not_installed++)
        {
            if(TRD_detector == Not_installed_TRD_detectors[i_not_installed])
            {
                flag_not_installed_TRD = 1;
                break;
            }
        }
        if(flag_not_installed_TRD) continue;


        Int_t TRD_color = kGray+1;
        Int_t flag_defect_TRD = 0;
        for(Int_t i_defect = 0; i_defect < 84; i_defect++)
        {
            if(TRD_detector == Defect_TRD_detectors[i_defect])
            {
                flag_defect_TRD = 1;
                break;
            }
        }
        if(flag_defect_TRD)
        {
            TRD_color = kGreen+1;
            //continue;
        }


        Int_t                TRD_sector         = fGeo         ->GetSector(TRD_detector);
        Int_t                TRD_stack          = fGeo         ->GetStack(TRD_detector);
        Int_t                TRD_layer          = fGeo         ->GetLayer(TRD_detector);
        Float_t              TRD_time0          = fGeo         ->GetTime0(TRD_layer);

        // OK
        Float_t              TRD_chamber_length = fGeo->GetChamberLength(TRD_layer,TRD_stack);
        Float_t              TRD_chamber_width  = fGeo->GetChamberWidth(TRD_layer);
        Float_t              TRD_chamber_height = 8.4;

        AliTRDpadPlane*      padplane           = fGeo         ->GetPadPlane(TRD_detector);
        Double_t             TRD_col_end        = padplane     ->GetColEnd();
        Double_t             TRD_row_end        = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
        Double_t             TRD_col_start      = padplane     ->GetCol0();
        Double_t             TRD_row_start      = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
        Double_t             TRD_row_end_ROC    = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
        Double_t             TRD_col_spacing    = padplane     ->GetColSpacing();
        Double_t             TRD_row_spacing    = padplane     ->GetRowSpacing();

        // OK

        Double_t Rotation_angle     = ((360.0/18.0)/2.0) + ((Double_t)TRD_sector)*(360.0/18.0);

        Double_t             loc[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glb[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,loc,glb);


        Double_t             locZ[3]           = {TRD_time0-50.0,0.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbZ[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locZ,glbZ);

        Double_t             locX[3]           = {TRD_time0,50.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbX[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locX,glbX);

        Double_t             locY[3]           = {TRD_time0,0.0,50.0+(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbY[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locY,glbY);


        combitrans[TRD_detector] = new TGeoCombiTrans();
        combitrans[TRD_detector] ->RotateZ(Rotation_angle + 90.0);
        combitrans[TRD_detector] ->SetTranslation(glb[0],glb[1],glb[2]);


        // Not OK
        vec_TV3_local_pos[0].SetXYZ(-TRD_chamber_width/2.0,-TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[1].SetXYZ(TRD_chamber_width/2.0,-TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[2].SetXYZ(TRD_chamber_width/2.0,TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[3].SetXYZ(-TRD_chamber_width/2.0,TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[4].SetXYZ(-TRD_chamber_width/2.0,-TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[5].SetXYZ(TRD_chamber_width/2.0,-TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[6].SetXYZ(TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[7].SetXYZ(-TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0);


        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            Double_t arr_pos_loc[3] = {vec_TV3_local_pos[i_vertex][0],vec_TV3_local_pos[i_vertex][1],vec_TV3_local_pos[i_vertex][2]};
            Double_t arr_pos_glb[3] = {0.0,0.0,0.0};
            combitrans[TRD_detector] ->LocalToMaster(arr_pos_loc,arr_pos_glb);

            vec_TH1D_TRD_geometry[0][i_vertex] ->SetBinContent(TRD_detector,arr_pos_glb[0]);
            vec_TH1D_TRD_geometry[1][i_vertex] ->SetBinContent(TRD_detector,arr_pos_glb[1]);
            vec_TH1D_TRD_geometry[2][i_vertex] ->SetBinContent(TRD_detector,arr_pos_glb[2]);
        }


        vec_TV3_TRD_center_offset[TRD_detector].SetXYZ(glb[0],glb[1],glb[2]);

        vec_TV3_TRD_center[TRD_detector][0].SetXYZ(glbX[0]-glb[0],glbX[1]-glb[1],glbX[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][1].SetXYZ(glbY[0]-glb[0],glbY[1]-glb[1],glbY[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][2].SetXYZ(glbZ[0]-glb[0],glbZ[1]-glb[1],glbZ[2]-glb[2]);
    }



    vec_self_tracklet_points_matched.resize(540);
    vec_self_tracklet_fit_points_matched.resize(540);
    vec_self_tracklet_points_background.resize(540);
    vec_self_tracklet_fit_points_background.resize(540);
    trkl_min.resize(540);
    trkl_min_background.resize(540);
    vec_self_tracklet_points_bckg.resize(540);
    vec_self_tracklet_fit_points_bckg.resize(540);
    trkl_min_bckg.resize(540);

    vec_ADC_val.clear();
    vec_ADC_val.resize(540);

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
TTRD_ST_Make_Tracklets::~TTRD_ST_Make_Tracklets()
{

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TTRD_ST_Make_Tracklets::Init_QA()
{
    printf("TTRD_ST_Make_Tracklets::Init_QA() \n");
    TFile* inputfile_QA = TFile::Open("/home/ceres/schmah/ALICE/TRD_Run3_calib/QA_out_year_2016_V2.root");
    tg_HV_drift_vs_det = (TGraph*)inputfile_QA->Get("tg_HV_drift_vs_det");
    tg_HV_anode_vs_det = (TGraph*)inputfile_QA->Get("tg_HV_anode_vs_det");;
    tg_vdrift_vs_det   = (TGraph*)inputfile_QA->Get("tg_vdrift_vs_det");;
}
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
// Driver function to sort the 2D vector 
// on basis of a particular column 
bool sortcol_first( const vector<Double_t>& v1,
             const vector<Double_t>& v2 )
{
 return v1[3] > v2[3];  // second column
} 
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t TTRD_ST_Make_Tracklets::calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
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
TVector3 TTRD_ST_Make_Tracklets::calculateDCA_vec_StraightToPoint(TVector3 &base, TVector3 &dir, TVector3 &point)
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
TVector3 TTRD_ST_Make_Tracklets::calculate_point_on_Straight_dca_to_Point(TVector3 &base, TVector3 &dir, TVector3 &point)
{
  // calculates the TVector3 on the straight line which is closest to point

    TVector3 diff = point - base;
    TVector3 dir_norm = dir;
    dir_norm *= (1.0/dir.Mag());
    Double_t proj_val = diff.Dot(dir_norm);
    TVector3 proj_dir = dir_norm;
    proj_dir *= proj_val;

    TVector3 dist_vec = proj_dir + base;

    return dist_vec;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TTRD_ST_Make_Tracklets::Reset()
{
    printf("TTRD_ST_Make_Tracklets::Reset() \n");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void TTRD_ST_Make_Tracklets::Init_tree(TString SEList)
{
    cout << "Initialize tree" << endl;
    //TString pinputdir = "/misc/alidata120/alice_u/schmah/TRD_offline_calib/Data/";
    //TString inlistdir = "/home/ceres/schmah/ALICE/TRD_self_tracking/Lists/";
    //TString pinputdir = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/Data/";
    //TString pinputdir = "/home/ceres/berdnikova/TRD-Run3-Calibration/";
	TString inlistdir = "/home/ceres/hoppner/ALICE/TRD_self_tracking/TRD-self-tracking/";
    TString pinputdir = "/home/ceres/hoppner/ALICE/TRD_self_tracking/TRD-self-tracking/Data/";
    
    TString in_list_name = SEList;
    SEList = inlistdir + SEList;

    AS_Event = new Ali_AS_Event();
    AS_Track = new Ali_AS_Track();
    AS_Tracklet = new Ali_AS_Tracklet();
    AS_Digit = new Ali_AS_TRD_digit();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( JPsi_TREE.Data(), JPsi_TREE.Data() );
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
                    input_SE ->AddFile(addfile.Data(),-1, JPsi_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( JPsi_BRANCH, &AS_Event );
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
        }
    }

    file_entries_total = input_SE->GetEntries();
    N_Events = file_entries_total;
    cout << "Total number of events in tree: " << file_entries_total << endl;


    //------------------------------------------------
    printf("Create output file \n");
    //TString outfile_name = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/Calib_tracklets/" + in_list_name + "_out_V2.root";
    TString outfile_name = "/home/ceres/hoppner/ALICE/TRD_self_tracking/TRD-self-tracking/Calib_tracklets/" + in_list_name + "_out_V2.root";
    //outputfile = new TFile("./Data/TRD_Calib_ADC_X1.root","RECREATE");
    outputfile = new TFile(outfile_name.Data(),"RECREATE");
    outputfile ->cd();
    // TRD self tracking output data containers
    TRD_ST_Tracklet   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event      = new Ali_TRD_ST_Event();

    Tree_TRD_ST_Event  = NULL;
    Tree_TRD_ST_Event  = new TTree("Tree_TRD_ST_Event" , "TRD_ST_Events" );
    Tree_TRD_ST_Event  ->Branch("Tree_TRD_ST_Event_branch"  , "TRD_ST_Event", TRD_ST_Event );
    Tree_TRD_ST_Event  ->SetAutoSave( 5000000 );
    //------------------------------------------------
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Int_t TTRD_ST_Make_Tracklets::Loop_event(Long64_t event)
{
    printf("Loop event number: %lld \n",event);

    Event_active = event;

    if (!input_SE->GetEntry( event )) return 0; // take the event -> information is stored in event

    N_Digits = 0;


    //---------------------------------------------------------------------------
    UShort_t NumTracks            = AS_Event ->getNumTracks(); // number of tracks in this event
    NumTracks_event = NumTracks;
    Double_t EventVertexX         = AS_Event ->getx();
    Double_t EventVertexY         = AS_Event ->gety();
    Double_t EventVertexZ         = AS_Event ->getz();
    Int_t    N_tracks_event       = AS_Event ->getN_tracks();
    Int_t    N_TRD_tracklets      = AS_Event ->getN_TRD_tracklets();
    Int_t    N_TRD_tracklets_online = AS_Event ->getNumTracklets(); // online tracklet
    Float_t  V0MEq                = AS_Event ->getcent_class_V0MEq();

    printf("N_TRD_tracklets_online: %d \n",N_TRD_tracklets_online);

    N_Tracks = NumTracks;

    vec_track_info.clear();
    vec_digit_track_info.clear();
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_Tracklet_pos.clear();
        vec_TV3_Tracklet_dir.clear();
    }

    // Loop over all tracklets
    Double_t scale_factor_length = 2.0;
    for(UShort_t i_tracklet = 0; i_tracklet < N_TRD_tracklets_online; ++i_tracklet) // loop over all tracklets of the actual event
    {
        AS_Tracklet             = AS_Event    ->getTracklet( i_tracklet ); // take the track
        TVector3 TV3_offset     = AS_Tracklet ->get_TV3_offset(); // online tracklets
        TVector3 TV3_dir        = AS_Tracklet ->get_TV3_dir();    // online tracklets
        Short_t  i_det_tracklet = AS_Tracklet ->get_detector();

        vec_TV3_Tracklet_pos[i_det_tracklet].push_back(TV3_offset);
        vec_TV3_Tracklet_dir[i_det_tracklet].push_back(TV3_dir);

        Double_t impact_angle = TV3_dir.Angle(vec_TV3_TRD_center[i_det_tracklet][2]);
        if(impact_angle > TMath::Pi()*0.5) impact_angle -= TMath::Pi();

        Float_t points[6] =
        {
            (Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),
            (Float_t)(TV3_offset[0] + scale_factor_length*TV3_dir[0]),(Float_t)(TV3_offset[1] + scale_factor_length*TV3_dir[1]),(Float_t)(TV3_offset[2] + scale_factor_length*TV3_dir[2])
        };
        //printf("i_tracklet: %d, out of %d, impact_angle: %4.3f, offset: {%4.3f, %4.3f, %4.3f}, end: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,N_TRD_tracklets_online,impact_angle*TMath::RadToDeg(),TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + TV3_dir[0],TV3_offset[1] + TV3_dir[1],TV3_offset[2] + TV3_dir[2]);
    }

    // Loop over all tracks
    for(UShort_t i_track = 0; i_track < NumTracks; ++i_track) // loop over all tracks of the actual event
    {
        //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
        AS_Track      = AS_Event ->getTrack( i_track ); // take the track
        Double_t nsigma_TPC_e   = AS_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = AS_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = AS_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = AS_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = AS_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = AS_Track ->getTRDSignal();
        Double_t TRDsumADC      = AS_Track ->getTRDsumADC();
        Double_t dca            = AS_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = AS_Track ->get_TLV_part();
        UShort_t NTPCcls        = AS_Track ->getNTPCcls();
        UShort_t NTRDcls        = AS_Track ->getNTRDcls();
        UShort_t NITScls        = AS_Track ->getNITScls();
        Float_t TPCchi2         = AS_Track ->getTPCchi2();
        Float_t TPCdEdx         = AS_Track ->getTPCdEdx();
        Float_t TOFsignal       = AS_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = AS_Track ->getTrack_length();

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();

        vec_track_single_info[0]  = dca;
        vec_track_single_info[1]  = TPCdEdx;
        vec_track_single_info[2]  = momentum;
        vec_track_single_info[3]  = eta_track;
        vec_track_single_info[4]  = pT_track;
        vec_track_single_info[5]  = TOFsignal;
        vec_track_single_info[6]  = Track_length;
        vec_track_single_info[7]  = TRDsumADC;
        vec_track_single_info[8]  = TRD_signal;
        vec_track_single_info[9]  = nsigma_TPC_e;
        vec_track_single_info[10] = nsigma_TPC_pi;
        vec_track_single_info[11] = nsigma_TPC_p;
        vec_track_info.push_back(vec_track_single_info);

        //----------------------------------------------
        // TRD digit information
        UShort_t  fNumTRDdigits        = AS_Track ->getNumTRD_digits();
        UShort_t  fNumOfflineTracklets = AS_Track ->getNumOfflineTracklets();
        //--------------------------


        //--------------------------
        // Offline tracklet loop
        for(Int_t i_tracklet = 0; i_tracklet < fNumOfflineTracklets; i_tracklet++) // layers
        {
            AS_offline_Tracklet     = AS_Track            ->getOfflineTracklet( i_tracklet ); // take the track
            TVector3 TV3_offset     = AS_offline_Tracklet ->get_TV3_offset(); // offline tracklets
            TVector3 TV3_dir        = AS_offline_Tracklet ->get_TV3_dir();    // offline tracklets
            Short_t  i_det_tracklet = AS_offline_Tracklet ->get_detector();

            printf("offline, i_tracklet: %d, offset: {%4.3f, %4.3f, %4.3f}, dir: {%4.3f, %4.3f, %4.3f} \n",i_tracklet,(Float_t)(TV3_offset[0]),(Float_t)(TV3_offset[1]),(Float_t)(TV3_offset[2]),(Float_t)TV3_dir[0],(Float_t)TV3_dir[1],(Float_t)TV3_dir[2]);
        }
        //--------------------------


        printf("i_track: %d, pT: %4.3f, phi: %4.3f, eta: %4.3f, N_off_trkl: %d \n",i_track,pT_track,phi_track,eta_track,fNumOfflineTracklets);

        //printf("i_track: %d, fNumTRDdigits: %d \n",i_track,fNumTRDdigits);

        vec_digit_info.clear();
        for(UShort_t i_digits = 0; i_digits < fNumTRDdigits; i_digits++)
        {
            //cout << "i_digits: " << i_digits << ", of " << fNumTRDdigits << endl;
            AS_Digit              = AS_Track ->getTRD_digit(i_digits);
            Int_t    layer        = AS_Digit ->get_layer();
            Int_t    sector       = AS_Digit ->get_sector();
            Int_t    column       = AS_Digit ->get_column();
            Int_t    stack        = AS_Digit ->get_stack();
            Int_t    row          = AS_Digit ->get_row();
            Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);
            Float_t  dca_to_track = AS_Digit ->getdca_to_track();
            Float_t  dca_x        = AS_Digit ->getdca_x();
            Float_t  dca_y        = AS_Digit ->getdca_y();
            Float_t  dca_z        = AS_Digit ->getdca_z();
            Float_t  ImpactAngle  = AS_Digit ->getImpactAngle();

            for(Int_t i_time = 0; i_time < 24; i_time++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    digit_pos[i_xyz] = AS_Digit ->get_pos(i_time,i_xyz);
                }
                //printf("track: %d/%d, digit: %d/%d \n",i_track,NumTracks,i_digits,fNumTRDdigits);
                //printf("pos: {%4.3f, %4.3f, %4.3f} \n,",digit_pos[0],digit_pos[1],digit_pos[2]);

                Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time);


                // x,y,z,time,ADC,sector,stack,layer,row,column,dca
                vec_digit_single_info[0]  = digit_pos[0];
                vec_digit_single_info[1]  = digit_pos[1];
                vec_digit_single_info[2]  = digit_pos[2];
                vec_digit_single_info[3]  = i_time;
                vec_digit_single_info[4]  = ADC;
                vec_digit_single_info[5]  = sector;
                vec_digit_single_info[6]  = stack;
                vec_digit_single_info[7]  = layer;
                vec_digit_single_info[8]  = row;
                vec_digit_single_info[9]  = column;
                vec_digit_single_info[10] = dca_to_track;
                vec_digit_single_info[11] = dca_x;
                vec_digit_single_info[12] = dca_y;
                vec_digit_single_info[13] = dca_z;

                //printf("dca_full: %4.3f, dca: {%4.3f, %4.3f, %4.3f} \n",dca_to_track,dca_x,dca_y,dca_z);
                vec_digit_info.push_back(vec_digit_single_info);

                N_Digits ++;
            }

        } // end of digits loop

        vec_digit_track_info.push_back(vec_digit_info);

    } // end of track loop


    return 1;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void TTRD_ST_Make_Tracklets::Make_clusters_and_get_tracklets_fit(Double_t Delta_x, Double_t Delta_z, Double_t factor_layer, Double_t factor_missing)
{
    printf("TTRD_ST_Make_Tracklets::Make_clusters_and_get_tracklets_fit() \n");

    //Reset();
    //printf("test 1 \n");

    vector< vector< vector< vector<Double_t> > > > vec_all_TRD_digits;

    vector< vector< vector< vector<Double_t> > > > vec_all_TRD_digits_clusters;
	
	vector< vector< vector<Int_t> > >  vec_used_clusters;
    //vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_points;


    vec_all_TRD_digits.resize(540);
    vec_all_TRD_digits_clusters.resize(540);
    vec_used_clusters.resize(540);
        
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_all_TRD_digits[i_det].resize(24); // time bins
        vec_all_TRD_digits_clusters[i_det].resize(24); // time bins
    	vec_used_clusters[i_det].resize(24); // time bins
    }
    vector<Double_t> vec_digit_data; // x,y,z,ADC
    vector<Double_t> vec_digit_cluster_data; // x,y,z,ADC
    vec_digit_data.resize(4);
    vec_digit_cluster_data.resize(4);

    // Fill all the information in the hierachy of detectors and time bins
    Int_t    N_TRD_digits  = AS_Event ->getNumTRD_digits();
    for(Int_t i_digit = 0; i_digit < N_TRD_digits; i_digit++)
    {
        AS_Digit              = AS_Event ->getTRD_digit(i_digit);
        Int_t    layer        = AS_Digit ->get_layer();
        Int_t    sector       = AS_Digit ->get_sector();
        Int_t    column       = AS_Digit ->get_column();
        Int_t    stack        = AS_Digit ->get_stack();
        Int_t    row          = AS_Digit ->get_row();
        Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);

        Double_t radius_prev = 0.0;
        for(Int_t i_time = 0; i_time < 24; i_time++)
        {
            Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time);
            if(ADC < 0.0) continue;
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                digit_pos[i_xyz] = AS_Digit ->get_pos(i_time,i_xyz);
                vec_digit_data[i_xyz] = digit_pos[i_xyz];
            }
            vec_digit_data[3] = ADC;

            vec_all_TRD_digits[detector][i_time].push_back(vec_digit_data);

            Double_t radius = TMath::Sqrt(TMath::Power(digit_pos[0] ,2) + TMath::Power(digit_pos[1] ,2));
            //if(i_time > 0 && radius_prev < radius) printf("det: %d, layer: %d, sector: %d, stack: %d, i_time: %d, radius: %4.3f \n",detector,layer,sector,stack,i_time,radius);
            radius_prev = radius;
        }
    }



    //-------------------------------------------------------
    // Make clusters for each detector and time bin
    // Individial digits -> clusters/time bin
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        Int_t sector = (Int_t)(i_det/30);
        Int_t stack  = (Int_t)(i_det%30/6);
        Int_t layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        //printf("i_det: %d \n",i_det);
        for(Int_t i_time = 0; i_time < 24; i_time++)
        {
            //if(!(i_det == 0 && i_time == 0)) continue;
            //for(Int_t i_digit = 0; i_digit < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit++)
            //{
            //    printf("i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, ADC: %4.3f \n",i_digit,vec_all_TRD_digits[i_det][i_time][i_digit][0],vec_all_TRD_digits[i_det][i_time][i_digit][1],vec_all_TRD_digits[i_det][i_time][i_digit][2],vec_all_TRD_digits[i_det][i_time][i_digit][3]);
            //}

            //cout << "Sorting vector" << endl;
            // First order the vector from high to low ADC values

            // vec_all_TRD_digits // 540 chambers, 24 time bins, (x,y,z,ADC)
            std::sort(vec_all_TRD_digits[i_det][i_time].begin(),vec_all_TRD_digits[i_det][i_time].end(),sortcol_first); // large values to small values, last column sorted via function sortcol

            //for(Int_t i_digit = 0; i_digit < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit++)
            //{
            //    printf("i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, ADC: %4.3f \n",i_digit,vec_all_TRD_digits[i_det][i_time][i_digit][0],vec_all_TRD_digits[i_det][i_time][i_digit][1],vec_all_TRD_digits[i_det][i_time][i_digit][2],vec_all_TRD_digits[i_det][i_time][i_digit][3]);
            //}


            vector<Int_t> arr_used_digits;
            arr_used_digits.clear();
            arr_used_digits.resize((Int_t)vec_all_TRD_digits[i_det][i_time].size());
            for(Int_t i_digit_max = 0; i_digit_max < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit_max++)
            {
                arr_used_digits[i_digit_max] = 0;
            }

            // Start from the maximum ADC value(s)
            for(Int_t i_digit_max = 0; i_digit_max < ((Int_t)vec_all_TRD_digits[i_det][i_time].size() - 1); i_digit_max++)
            {
                if(arr_used_digits[i_digit_max]) continue;
                arr_used_digits[i_digit_max] = 1;

                Double_t pos_ADC_max[4] = {vec_all_TRD_digits[i_det][i_time][i_digit_max][0],vec_all_TRD_digits[i_det][i_time][i_digit_max][1],vec_all_TRD_digits[i_det][i_time][i_digit_max][2],vec_all_TRD_digits[i_det][i_time][i_digit_max][3]};
                Double_t N_digits_added = pos_ADC_max[3]; // ADC as weight
                Double_t pos_ADC_sum[4] = {pos_ADC_max[0]*pos_ADC_max[3],pos_ADC_max[1]*pos_ADC_max[3],pos_ADC_max[2]*pos_ADC_max[3],pos_ADC_max[3]}; // positions times ADC value [3]

                //printf("i_time: %d, i_digit_max: %d, ADC: %4.3f \n",i_time,i_digit_max,pos_ADC_max[3]);

                // Get all other digits within a certain radius
                for(Int_t i_digit_sub = (i_digit_max + 1); i_digit_sub < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit_sub++)
                {
                    if(arr_used_digits[i_digit_sub]) continue;

                    Double_t pos_ADC_sub[4] = {vec_all_TRD_digits[i_det][i_time][i_digit_sub][0],vec_all_TRD_digits[i_det][i_time][i_digit_sub][1],vec_all_TRD_digits[i_det][i_time][i_digit_sub][2],vec_all_TRD_digits[i_det][i_time][i_digit_sub][3]};
                    Double_t dist_digits_XY = TMath::Sqrt(TMath::Power(pos_ADC_max[0] - pos_ADC_sub[0],2) + TMath::Power(pos_ADC_max[1] - pos_ADC_sub[1],2));
                    Double_t dist_digits_Z  = fabs(pos_ADC_max[2] - pos_ADC_sub[2]);
                    if(dist_digits_XY > 2.4)  continue;
                    if(dist_digits_Z  > 15.0) continue;

                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        pos_ADC_sum[i_xyz] += pos_ADC_sub[3]*pos_ADC_sub[i_xyz];
                    }
                    pos_ADC_sum[3] += pos_ADC_sub[3];

                    arr_used_digits[i_digit_sub] = 1;
                    N_digits_added += pos_ADC_sub[3];
                }

                if(N_digits_added <= 0.0) continue;
                // Calculate average position
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    pos_ADC_sum[i_xyz] /= N_digits_added;
                }

                //TEve_clusters ->SetNextPoint(pos_ADC_sum[0],pos_ADC_sum[1],pos_ADC_sum[2]);
                for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
                {
                    vec_digit_cluster_data[i_xyzADC] = pos_ADC_sum[i_xyzADC];
                }

                //printf("i_time: %d, cluster pos: {%4.3f, %4.3f, %4.3f} \n",i_time,vec_digit_cluster_data[0],vec_digit_cluster_data[1],vec_digit_cluster_data[2]);

                vec_all_TRD_digits_clusters[i_det][i_time].push_back(vec_digit_cluster_data);
				vec_used_clusters[i_det][i_time].push_back(0);
            }
        } // end of timebin loop
    }
    printf("arrangement of clusters done \n");
    //-------------------------------------------------------


    //create "tracklets"
    //Int_t sector = (Int_t)(i_det/30);
    //Int_t stack  = (Int_t)(i_det%30/6);
    //Int_t layer  = i_det%6;
    //Int_t i_det = layer + 6*stack + 30*sector;


    //-------------------------------------------------------
    // Connect clusters within each chamber
    // Clusters/time bin -> connect them time bin wise
    vec_self_tracklet_points.clear();
    vec_self_tracklet_points.resize(540);
	Int_t min_nbr_cls = 15;
	
    for(Int_t i_det = 0; i_det < 540; i_det++) // is done chamber wise
    {
		    
        for(Int_t i_time = 0; i_time < 24 - min_nbr_cls; i_time++) // is done chamber wise
    	{
			
    		Int_t N_clusters = (Int_t)vec_all_TRD_digits_clusters[i_det][i_time].size();
			

			//printf("i_det: %d, N_clusters: %d \n",i_det,N_clusters);

			vec_self_tracklet_points[i_det].resize(N_clusters);

			for(Int_t i_cls = 0; i_cls < N_clusters; i_cls++) // loop over all clusters in one detector and for time bin 0
			{	
				if(vec_used_clusters[i_det][i_time][i_cls]) continue;
				Int_t n_clusters_attached = 0;
				Double_t pos_ADC_max[4] = {vec_all_TRD_digits_clusters[i_det][i_time][i_cls][0],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][1],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][2],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][3]};

				Double_t radius = TMath::Sqrt(TMath::Power(vec_all_TRD_digits_clusters[i_det][i_time][i_cls][0],2) + TMath::Power(vec_all_TRD_digits_clusters[i_det][i_time][i_cls][1],2));

				//printf("---> i_cls: %d, time 0, pos: {%4.3f, %4.3f, %4.3f}, radius: %4.3f \n",i_cls,vec_all_TRD_digits_clusters[i_det][i_time][i_cls][0],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][1],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][2],radius);

				vec_self_tracklet_points[i_det][i_cls].resize(24);

				for(Int_t i_timebin = 0; i_timebin < 24; i_timebin++)
				{
					vec_self_tracklet_points[i_det][i_cls][i_timebin].resize(4);
					for (Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
					{
						vec_self_tracklet_points[i_det][i_cls][i_timebin][i_xyzADC] = -999.0;
					}
				}

				for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
				{
					vec_self_tracklet_points[i_det][i_cls][i_time][i_xyzADC] = pos_ADC_max[i_xyzADC];
				}

				Int_t    i_time_start    = i_time + 1;
				Double_t scale_fac_add   = 1.0;
				Int_t    missed_time_bin = 0;
				Double_t radius_prev     = 0.0;

				for(Int_t i_time_sub = i_time_start; i_time_sub < 24; i_time_sub++)
				{
					Double_t scale_fac = 1.0*scale_fac_add;
					//printf("i_time_sub: %d, scale_fac: %4.3f \n",i_time_sub,scale_fac);
					//if(i_time_sub == 0) scale_fac = factor_layer*scale_fac_add;
					Int_t N_clusters_sub = (Int_t)vec_all_TRD_digits_clusters[i_det][i_time_sub].size();

					Int_t best_sub_cluster = -1;
					Double_t best_cluster_quality = 1000000.0;

					for(Int_t i_cls_sub = 0; i_cls_sub < N_clusters_sub; i_cls_sub++)
					{
						Double_t pos_ADC_sub[4] = {vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][0],vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][1],vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][2],vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][3]};

						Double_t dist_clusters_XY = TMath::Sqrt(TMath::Power(pos_ADC_max[0] - pos_ADC_sub[0],2) + TMath::Power(pos_ADC_max[1] - pos_ADC_sub[1],2));
						Double_t dist_clusters_Z  = fabs(pos_ADC_max[2] - pos_ADC_sub[2]);

						if(dist_clusters_XY > scale_fac*Delta_x)  continue;
						if(dist_clusters_Z  > Delta_z) continue;
						//if(dist_clusters_XY > scale_fac*3.0)  continue;
						//if(dist_clusters_Z  > 10.0) continue;

						// Matching quality - chi2 like
						Double_t sub_cluster_quality = TMath::Power(dist_clusters_XY/0.7,2.0) + TMath::Power(dist_clusters_Z/7.5,2.0);


						if(sub_cluster_quality < best_cluster_quality)
						{
							best_cluster_quality = sub_cluster_quality;
							best_sub_cluster     = i_cls_sub;
						}
					}


					//cout << "best_sub_cluster: " << best_sub_cluster << ", i_time_sub: " << i_time_sub << ", i_layer_sub: " << i_layer_sub << endl;
					if(best_sub_cluster < 0)
					{
						scale_fac_add *= factor_missing; // one time bin was missing, increase matching window
						missed_time_bin++;
						continue;
					}

					if(missed_time_bin > 3) break;
					scale_fac_add = 1.0; // reset additional matching window factor once a match was found

					// Define new pos_ADC_max
					for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
					{
						pos_ADC_max[i_xyzADC] = vec_all_TRD_digits_clusters[i_det][i_time_sub][best_sub_cluster][i_xyzADC];
						vec_self_tracklet_points[i_det][i_cls][i_time_sub][i_xyzADC] = pos_ADC_max[i_xyzADC];

					}
					vec_used_clusters[i_det][i_time_sub][best_sub_cluster] = 1;
					Double_t radius_sub = TMath::Sqrt(TMath::Power(vec_self_tracklet_points[i_det][i_cls][i_time_sub][0],2) + TMath::Power(vec_self_tracklet_points[i_det][i_cls][i_time_sub][1],2));


					// Already sometimes wrong radius-time ordering
					//if(i_time_sub > 1 && radius_prev < radius_sub) printf("---> i_det: %d, i_cls: %d, time %d, pos: {%4.3f, %4.3f, %4.3f}, radius_sub: %4.3f, radius_prev: %4.3f \n",i_det,i_cls,i_time_sub,vec_self_tracklet_points[i_det][i_cls][i_time_sub][0],vec_self_tracklet_points[i_det][i_cls][i_time_sub][1],vec_self_tracklet_points[i_det][i_cls][i_time_sub][2],radius_sub,radius_prev);
					radius_prev = radius_sub;

					//vec_TEveLine_cluster_tracks[(Int_t)vec_TEveLine_cluster_tracks.size()-1] ->SetNextPoint((Float_t)pos_ADC_max[0],(Float_t)pos_ADC_max[1],(Float_t)pos_ADC_max[2]);
					//vec_TPL_cluster_tracks[i_sector][i_stack][(Int_t)vec_TPL_cluster_tracks[i_sector][i_stack].size()-1] ->SetNextPoint((Float_t)pos_ADC_max[0],(Float_t)pos_ADC_max[1]);

					n_clusters_attached++;
				}
			}
		}	
    }

    printf("connection of clusters within detector done \n");
    printf("tracklets fit starts now \n");
    //-------------------------------------------------------



    //-------------------------------------------------------
    // Fit the time bin wise connected clusters

    //ready to fit vec_self_tracklet_points[i_det][i_cls][i_time_sub][i_xyzADC]

    // Is fitting the tracklets and doing the global fit through all first cluster points of all available layers
    //printf("TTRD_ST_Make_Tracklets::get_tracklets_fit((%d) \n",i_track);

    //fit merged digits with a straight line
    vec_self_tracklet_fit_points.clear();
    vec_self_tracklet_fit_points.resize(540);         //[i_det][i_trkl][i_start_stop][i_xyz]

    for(Int_t i_detector = 0; i_detector < 540; i_detector++)
    {
        vec_self_tracklet_fit_points[i_detector].resize((Int_t)vec_self_tracklet_points[i_detector].size());

        vec_ADC_val[i_detector].resize((Int_t)vec_self_tracklet_points[i_detector].size());

        for (Int_t i_trkl = 0; i_trkl < (Int_t)vec_self_tracklet_points[i_detector].size(); i_trkl++)
        {
            vec_self_tracklet_fit_points[i_detector][i_trkl].resize(2);
            for(Int_t i_start_stop = 0; i_start_stop < 2; i_start_stop++)
            {
                vec_self_tracklet_fit_points[i_detector][i_trkl][i_start_stop].resize(3);
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    vec_self_tracklet_fit_points[i_detector][i_trkl][i_start_stop][i_xyz] = -999.0;
                }
            }
        }
    }


    Double_t p0[4] = {10,20,1,2};
    //Int_t i_layer_notempty = 0;

    //declare things that will be needed later
    Double_t arglist[10];
    Double_t a0[3] = {0,0,0};
    Double_t a1[3] = {0,0,0};

    TVector3 vec_a0;
    TVector3 vec_a1;
    TVector3 vec_u;
    TVector3 vec_x0;
    TVector3 vec_u_perp;
    TVector3 vec_x1;

    Double_t pStart[4]; //= {1,1,1,1};

    int nvpar,nparx;
    Double_t amin,edm, errdef;

    Double_t parFit[4];

    Int_t    n;
    Double_t t0;
    Double_t dt;
    TVector3 TV3_line_point;
    Double_t i_point;


    //loop over layers
    Double_t layer_dist_min_max[6][2] =
    {
        {280.0,290.0},
        {290.0,300.0},
        {300.0,310.0},
        {310.0,320.0},
        {320.0,330.0},
        {330.0,340.0},
    };

    //Double_t delta_layer = 12.5;

    self_tracklets_min.resize(540);


    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        //if(i_det == 0) printf("fitting detector %d  \n",i_det);
        //if(i_det % 40 == 0) printf("fitting detector %d  \n",i_det);

        self_tracklets_min[i_det].resize(vec_self_tracklet_points[i_det].size());

        //printf("i_det: %d, N_tracklets: %d \n",i_det,(Int_t)vec_self_tracklet_points[i_det].size());

        for(Int_t i_trkl = 0; i_trkl < (Int_t)vec_self_tracklet_points[i_det].size(); i_trkl++)
        {
            self_tracklets_min[i_det][i_trkl] = -1.0;

            vec_ADC_val[i_det][i_trkl].resize((Int_t)vec_self_tracklet_points[i_det][i_trkl].size());

#if 0
            for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_self_tracklet_points[i_det][i_trkl].size(); i_time_merge++)
            {
                //printf("layer: %d, i_time_merge: %d, point: {%4.3f, %4.3f, %4.3f} \n",i_layer,i_time_merge,vec_Dt_digit_pos_cluster[i_layer][i_time_merge][0],vec_Dt_digit_pos_cluster[i_layer][i_time_merge][1],vec_Dt_digit_pos_cluster[i_layer][i_time_merge][2]);
            }
#endif

            Int_t i_time_merge_AB[2] = {-1,-1};
			Int_t number_ok_clusters = 0;
            for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_self_tracklet_points[i_det][i_trkl].size(); i_time_merge++)
            {
                if(vec_self_tracklet_points[i_det][i_trkl][i_time_merge][0] == -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][1] == -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][2] == -999.0) continue;
                else 
				{
					i_time_merge_AB[0] = i_time_merge;
					number_ok_clusters++;
				}
            }
			/*for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_self_tracklet_points[i_det][i_trkl].size(); i_time_merge++)
            {
                if(vec_self_tracklet_points[i_det][i_trkl][i_time_merge][0] != -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][1] != -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][2] != -999.0)
                {
					i_time_merge_AB[0] = i_time_merge; 
				 	break;
				}
            }*/

            if(i_time_merge_AB[0] == -1) continue; // all values are 0

            for(Int_t i_time_merge = ((Int_t)vec_self_tracklet_points[i_det][i_trkl].size() - 1); i_time_merge >= 0; i_time_merge--)
            {
                if(vec_self_tracklet_points[i_det][i_trkl][i_time_merge][0] == -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][1] == -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][2] == -999.0) continue;
                else i_time_merge_AB[1] = i_time_merge;
            }
			/*for(Int_t i_time_merge = ((Int_t)vec_self_tracklet_points[i_det][i_trkl].size() - 1); i_time_merge >= 0; i_time_merge--)
            {
                if(vec_self_tracklet_points[i_det][i_trkl][i_time_merge][0] == -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][1] == -999.0 && vec_self_tracklet_points[i_det][i_trkl][i_time_merge][2] == -999.0)
				{
					i_time_merge_AB[1] = i_time_merge;
					break;
				}	
            }*/

            if(i_time_merge_AB[0] == i_time_merge_AB[1]) continue; // no fit possible with just one point
			if(number_ok_clusters < min_nbr_cls) continue;
            //printf("TTRD_ST_Make_Tracklets::get_tracklets_fit(%d), i_layer: %d \n",i_track,i_layer);

            TVirtualFitter *min = TVirtualFitter::Fitter(0,4);
            //min->SetObjectFit(tracklets_gr);

            arglist[0] = 3;
            //min->ExecuteCommand("SET PRINT",arglist,1);
            Double_t arglist_B[1] = {-1};
            min->ExecuteCommand("SET PRIntout",arglist_B,1);
            min->ExecuteCommand("SET NOWarnings",arglist_B,1);

            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                a0[i_xyz] = vec_self_tracklet_points[i_det][i_trkl][i_time_merge_AB[0]][i_xyz];
                a1[i_xyz] = vec_self_tracklet_points[i_det][i_trkl][i_time_merge_AB[1]][i_xyz];
            }

            Det = i_det;
            Trkl = i_trkl;

            //printf("point start: {%4.3f, %4.3f, %4.3f} \n",a0[0],a0[1],a0[2]);
            //printf("point end: {%4.3f, %4.3f, %4.3f} \n",a1[0],a1[1],a1[2]);

            Int_t flag_XZ = 1;
            if(fabs(a0[0] - a1[0]) > fabs(a0[2] - a1[2]))
            {
                flag_XZ = 0; // x is used
            }

            vec_a0.SetXYZ(a0[0],a0[1],a0[2]);
            vec_a1.SetXYZ(a1[0],a1[1],a1[2]);
            vec_u = vec_a1 - vec_a0;
            if(flag_XZ == 0)
            {
                min->SetFCN(SumDistance2_self_X_tr);
                vec_x0 = vec_a0 - vec_u*(vec_a0.X()/vec_u.X());
                vec_u_perp.SetXYZ(vec_u[0],vec_u[1],vec_u[2]);
                vec_u_perp *= 1.0/vec_u[0];
                //printf("X_tr \n");
            }
            else
            {
                min->SetFCN(SumDistance2_self_tr);
                vec_x0 = vec_a0 - vec_u*(vec_a0.Z()/vec_u.Z());
                vec_u_perp.SetXYZ(vec_u[0],vec_u[1],vec_u[2]);
                vec_u_perp *= 1.0/vec_u[2];
                //printf("tr \n");
            }
            vec_x1 = vec_x0 + vec_u_perp;

            //TVector3 u = a1-a0;
            //x0 = a0  (a0z/uz)*u;
            //x1 = x0 + u / uz;

            if(flag_XZ == 0)
            {
                pStart[0] = vec_x0.Z();
                pStart[1] = vec_x1.Z() - pStart[0];
                pStart[2] = vec_x0.Y();
                pStart[3] = vec_x1.Y() - pStart[2];
            }
            else
            {
                pStart[0] = vec_x0.X();
                pStart[1] = vec_x1.X() - pStart[0];
                pStart[2] = vec_x0.Y();
                pStart[3] = vec_x1.Y() - pStart[2];
            }

            //cout << "pStart[0]" << pStart[0] << endl;
            //cout << "pStart[1]" << pStart[1] << endl;
            //cout << "pStart[2]" << pStart[2] << endl;
            //cout << "pStart[3]" << pStart[3] << endl;

            min->SetParameter(0,"x0",pStart[0],0.01,0,0);
            min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
            min->SetParameter(2,"y0",pStart[2],0.01,0,0);
            min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

            //printf("test 5.7 \n");

            arglist[0] = 1000; // number of function calls
            arglist[1] = 0.001; // tolerance

            //printf("test 5.71 \n");
            //printf("     ------------------ MIGRAD ------------------ \n");
            min->ExecuteCommand("MIGRAD",arglist,2);
            //printf("     ------------------ END ------------------ \n");
            //printf("-------------------------> test 5.72 \n");

            //if (minos) min->ExecuteCommand("MINOS",arglist,0);

            min->GetStats(amin,edm,errdef,nvpar,nparx);
            //min->PrintResults(1,amin);

            // get fit parameters
            for(int i = 0; i < 4; ++i)
            {
                parFit[i] = min->GetParameter(i);
                //parFit[i] = pStart[i];
            }

            self_tracklets_min[i_det][i_trkl] = amin;


            //-------------------------------------------------------
            // Calculate tracklet base and direction vectors

            // Arbitrary space point on fitted line
            Double_t x_A, y_A, z_A;
            if(flag_XZ == 0) line_X(0.0,parFit,x_A,y_A,z_A);
            else line(0.0,parFit,x_A,y_A,z_A);
            TVector3 TV3_base_fit(x_A,y_A,z_A);

            // normalized direction vector of fitted line
            if(flag_XZ == 0) line_X(1.0,parFit,x_A,y_A,z_A);
            else line(1.0,parFit,x_A,y_A,z_A);
            TVector3 TV3_dir_fit(x_A,y_A,z_A);

            TV3_dir_fit -= TV3_base_fit;
            Double_t dir_length = TV3_dir_fit.Mag();
            if(dir_length > 0.0) TV3_dir_fit *= 1.0/dir_length;

            // First space point of fitted clusters
			
            TVector3 TV3_t0_point(vec_self_tracklet_points[i_det][i_trkl][0][0],vec_self_tracklet_points[i_det][i_trkl][0][1],vec_self_tracklet_points[i_det][i_trkl][0][2]);

            // Space point on straight line which is closes to first space point of fitted clusters
            //TVector3 TV3_base_fit_t0 = calculate_point_on_Straight_dca_to_Point(TV3_base_fit,TV3_dir_fit,TV3_t0_point);
			TVector3 TV3_base_plane = vec_TV3_TRD_center_offset[i_det];
			TVector3 TV3_norm_plane = vec_TV3_TRD_center[i_det][0];
			TVector3 TV3_base_fit_t0 = intersect_line_plane(TV3_base_fit,TV3_dir_fit,TV3_base_plane,TV3_norm_plane);

            TVector3 vec_AB[2];
            vec_AB[0] = TV3_base_fit_t0;
            vec_AB[1] = TV3_base_fit_t0 + TV3_dir_fit;
            //if(vec_AB[1].Mag() > vec_AB[0].Mag())  // changed sign -> correct physical direction, pointing inwards
            if(vec_AB[1].Mag() < vec_AB[0].Mag())  // changed sign -> "uncorrect" direction, pointing outwards
            {
                TV3_dir_fit *= -1.0;
                vec_AB[1] = TV3_base_fit_t0 + TV3_dir_fit;
            }

            for(Int_t i_AB = 0; i_AB < 2; i_AB++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    vec_self_tracklet_fit_points[i_det][i_trkl][i_AB][i_xyz] = vec_AB[i_AB][i_xyz];
                }
            }

            // fill ADC val vector
            Int_t tbn_max = (Int_t)vec_self_tracklet_points[i_det][i_trkl].size();
            //Double_t radius_in  = TMath::Sqrt( TMath::Power(vec_self_tracklet_points[i_det][i_trkl][tbn_max-1][0],2) + TMath::Power(vec_self_tracklet_points[i_det][i_trkl][tbn_max-1][1],2) );
            //Double_t radius_out = TMath::Sqrt( TMath::Power(vec_self_tracklet_points[i_det][i_trkl][0][0],2) + TMath::Power(vec_self_tracklet_points[i_det][i_trkl][0][1],2) );

            for(Int_t i_timebin = 0; i_timebin < tbn_max; i_timebin++)
            {
                vec_ADC_val[i_det][i_trkl][i_timebin] = vec_self_tracklet_points[i_det][i_trkl][i_timebin][3];
            }

            //-------------------------------------------------------

            Double_t radius = TMath::Sqrt( TMath::Power(TV3_base_fit_t0[0],2) + TMath::Power(TV3_base_fit_t0[1],2) );
            //printf("amin: %4.3f, par: {%4.3f, %4.3f, %4.3f, %4.3f} \n",amin,parFit[0],parFit[1],parFit[2],parFit[3]);
            //printf("   --> radius: %4.3f, point first cluster: {%4.3f, %4.3f, %4.3f}, point line: {%4.3f, %4.3f, %4.3f} \n",radius,TV3_t0_point[0],TV3_t0_point[1],TV3_t0_point[2],TV3_base_fit_t0[0],TV3_base_fit_t0[1],TV3_base_fit_t0[2]);

            delete min;
        }
    }
    //-------------------------------------------------------

}

//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Int_t TTRD_ST_Make_Tracklets::Calibrate(Double_t Delta_x, Double_t Delta_z, Double_t factor_layer, Double_t factor_missing)
{
    printf("TTRD_ST_Make_Tracklets::Loop_full_self_tracking() \n");

    TVector3 TV3_trkl_offset;
    TVector3 TV3_trkl_dir;
    Float_t  helix_par[9];
    Double_t ADC_val[24];

    for(Long64_t i_event = 0; i_event < file_entries_total; i_event++)
    //for(Long64_t i_event = 0; i_event < 10; i_event++)
    {
        if(i_event % 5 == 0) printf("i_event: %lld out of %lld \n",i_event,file_entries_total);
        if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event


        //------------------------------------------
        // Fill event information
        TRD_ST_Event ->clearTrackList();
        TRD_ST_Event ->clearTrackletList();

        TRD_ST_Event ->setTriggerWord( AS_Event->getTriggerWord() );
        TRD_ST_Event ->setx( AS_Event->getx() );
        TRD_ST_Event ->sety( AS_Event->gety() );
        TRD_ST_Event ->setz( AS_Event->getz() );
        TRD_ST_Event ->setid( AS_Event->getid() );
        TRD_ST_Event ->setN_tracks( AS_Event->getN_tracks() );
        TRD_ST_Event ->setN_TRD_tracklets( AS_Event->getN_TRD_tracklets() );
        TRD_ST_Event ->setBeamIntAA( AS_Event->getBeamIntAA() );
        TRD_ST_Event ->setT0zVertex( AS_Event->getT0zVertex() );
        TRD_ST_Event ->setcent_class_ZNA( AS_Event->getcent_class_ZNA() );
        TRD_ST_Event ->setcent_class_ZNC( AS_Event->getcent_class_ZNC() );
        TRD_ST_Event ->setcent_class_V0A( AS_Event->getcent_class_V0A() );
        TRD_ST_Event ->setcent_class_V0C( AS_Event->getcent_class_V0C() );
        TRD_ST_Event ->setcent_class_V0M( AS_Event->getcent_class_V0M() );
        TRD_ST_Event ->setcent_class_CL0( AS_Event->getcent_class_CL0() );
        TRD_ST_Event ->setcent_class_CL1( AS_Event->getcent_class_CL1() );
        TRD_ST_Event ->setcent_class_SPD( AS_Event->getcent_class_SPD() );
        TRD_ST_Event ->setcent_class_V0MEq( AS_Event->getcent_class_V0MEq() );
        TRD_ST_Event ->setcent_class_V0AEq( AS_Event->getcent_class_V0AEq() );
        TRD_ST_Event ->setcent_class_V0CEq( AS_Event->getcent_class_V0CEq() );

        //printf("Event information filled \n");
        //------------------------------------------



        //------------------------------------------
        // Loop over all tracks
        UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event
        for(UShort_t i_track = 0; i_track < NumTracks; ++i_track) // loop over all tracks of the actual event
        {
            TRD_ST_TPC_Track = TRD_ST_Event ->createTrack(); // TPC track
            //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
            AS_Track      = AS_Event ->getTrack( i_track ); // take the track
            TRD_ST_TPC_Track ->setnsigma_e_TPC( AS_Track ->getnsigma_e_TPC() );
            TRD_ST_TPC_Track ->setnsigma_pi_TPC( AS_Track ->getnsigma_pi_TPC() );
            TRD_ST_TPC_Track ->setnsigma_p_TPC( AS_Track ->getnsigma_p_TPC() );
            TRD_ST_TPC_Track ->setnsigma_e_TOF( AS_Track ->getnsigma_e_TOF() );
            TRD_ST_TPC_Track ->setnsigma_pi_TOF( AS_Track ->getnsigma_pi_TOF() );
            TRD_ST_TPC_Track ->setTRDSignal( AS_Track ->getTRDSignal() );
            TRD_ST_TPC_Track ->setTRDsumADC( AS_Track ->getTRDsumADC() );
            TRD_ST_TPC_Track ->setdca( AS_Track ->getdca() );  // charge * distance of closest approach to the primary vertex
            TRD_ST_TPC_Track ->set_TLV_part( AS_Track ->get_TLV_part() );
            TRD_ST_TPC_Track ->setNTPCcls( AS_Track ->getNTPCcls() );
            TRD_ST_TPC_Track ->setNTRDcls( AS_Track ->getNTRDcls() );
            TRD_ST_TPC_Track ->setNITScls( AS_Track ->getNITScls() );
            TRD_ST_TPC_Track ->setTPCchi2( AS_Track ->getTPCchi2() );
            TRD_ST_TPC_Track ->setTPCdEdx( AS_Track ->getTPCdEdx() );
            TRD_ST_TPC_Track ->setTOFsignal( AS_Track ->getTOFsignal() ); // in ps (1E-12 s)
            TRD_ST_TPC_Track ->setTrack_length( AS_Track ->getTrack_length() );

            for(Int_t i_helix_par = 0; i_helix_par < 9; i_helix_par++)
            {
                helix_par[i_helix_par] =  AS_Track ->getHelix_param(i_helix_par);
            }

            TRD_ST_TPC_Track ->setHelix(helix_par[0],helix_par[1],helix_par[2],helix_par[3],helix_par[4],helix_par[5],helix_par[6],helix_par[7],helix_par[8]);
        }

        //printf("Track information filled \n");
        //------------------------------------------




        //------------------------------------------
        // Loop over all TRD tracklets
        Make_clusters_and_get_tracklets_fit(Delta_x,Delta_z,factor_layer,factor_missing);
        for(Int_t i_det = 0; i_det < 540; i_det++)
        {
            //printf("i_det: %d \n",i_det);
            for(Int_t i_trkl = 0; i_trkl < (Int_t)vec_self_tracklet_fit_points[i_det].size(); i_trkl++)
            {
                for (Int_t i_tbn = 0; i_tbn < 24; i_tbn++)
                {
                    ADC_val[i_tbn] = {-999.0};
                }

                //printf("i_trkl: %d \n",i_trkl);
                for(Int_t i_AB = 0; i_AB < 2; i_AB++)
                {
                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        //printf("i_AB: %d, i_xyz: %d \n",i_AB,i_xyz);
                        if(i_AB == 0) TV3_trkl_offset[i_xyz] = vec_self_tracklet_fit_points[i_det][i_trkl][i_AB][i_xyz];
                        if(i_AB == 1) TV3_trkl_dir[i_xyz]    = vec_self_tracklet_fit_points[i_det][i_trkl][i_AB][i_xyz];
                    }
                }

                TV3_trkl_dir -= TV3_trkl_offset;

                TRD_ST_Tracklet = TRD_ST_Event ->createTracklet(); // online tracklet
                TRD_ST_Tracklet  ->set_TRD_det(i_det);
                TRD_ST_Tracklet  ->set_TV3_offset(TV3_trkl_offset);
                TRD_ST_Tracklet  ->set_TV3_dir(TV3_trkl_dir);
                for(Int_t i_timebin = 0; i_timebin < (Int_t)vec_ADC_val[i_det][i_trkl].size(); i_timebin++)
                {
                    ADC_val[i_timebin] = vec_ADC_val[i_det][i_trkl][i_timebin];
                    TRD_ST_Tracklet  ->set_ADC_val(i_timebin,ADC_val[i_timebin]);
                }

            } // end tracklet loop
        } // end detector loop

        //printf("Tracklet information filled \n");
        //------------------------------------------


        //printf("Fill tree \n");
        Tree_TRD_ST_Event ->Fill();
        printf("Tree filled for event %lld \n",i_event);
    }

    printf("Write data to file \n");
    outputfile ->cd();
    Tree_TRD_ST_Event ->Write();

#if 0
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] ->Write();
        }

    }
#endif

    printf("All data written \n");

    return 1;
}
//----------------------------------------------------------------------------------------

