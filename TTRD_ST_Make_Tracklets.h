#ifndef __TTRD_MAKE_TRACKLETS_H__
#define __TTRD_MAKE_TRACKLETS_H__


#define USEEVE

using namespace std;
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TChain.h"

#include "TObject.h"

#include "Ali_AS_Event.h"
#include "Ali_AS_EventLinkDef.h"
#include "Ali_TRD_ST.h"
#include "Ali_TRD_ST_LinkDef.h"
#include "Ana_Digits_functions.h"
#include "TDecompSVD.h"
#include "TMatrixD.h"

#if defined(USEEVE)
#include "TEveBox.h"
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEvePointSet.h"
#endif


ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_Event)
ClassImp(Ali_TRD_ST_Tracklets)
ClassImp(Ali_TRD_ST_TPC_Track)
ClassImp(Ali_TRD_ST_Event)



//----------------------------------------------------------------------------------------
class TTRD_ST_Make_Tracklets
{
private:


    Ali_AS_Event*     AS_Event;
    Ali_AS_Track*     AS_Track;
    Ali_AS_Tracklet*  AS_Tracklet;
    Ali_AS_offline_Tracklet*  AS_offline_Tracklet;
    Ali_AS_TRD_digit* AS_Digit;


    Ali_TRD_ST_Tracklets* TRD_ST_Tracklet;
    Ali_TRD_ST_TPC_Track* TRD_ST_TPC_Track;
    Ali_TRD_ST_Event*     TRD_ST_Event;
    TTree* Tree_TRD_ST_Event;

    TFile* outputfile;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;

    TGraph* tg_HV_drift_vs_det;
    TGraph* tg_HV_anode_vs_det;
    TGraph* tg_vdrift_vs_det;

    Long64_t Event_active = 0;
    Double_t tracklets_min[7] = {-1.0};
    vector< vector<Double_t> > self_tracklets_min;
    vector< vector<Double_t> > trkl_min; //for matched self tracklets
    vector< vector<Double_t> > trkl_min_background; //for background self tracklets
    vector< vector<Double_t> > trkl_min_bckg; //for background self tracklets


    TChain* input_SE;
    TString JPsi_TREE   = "Tree_AS_Event";
    TString JPsi_BRANCH = "Tree_AS_Event_branch";
    Long64_t file_entries_total;

    Float_t digit_pos[3];

    vector<Float_t> vec_digit_single_info;
    vector< vector<Float_t> > vec_digit_info;
    vector< vector< vector<Float_t> > > vec_digit_track_info;
    vector<Float_t> vec_track_single_info;
    vector< vector<Float_t> > vec_track_info;
    vector<TVector3> vec_TV3_local_pos;



    TString HistName;
    char NoP[50];


    // TVector3 vec_datapoints;
    vector <TCanvas*> ADC_vs_time;

    TH1I* h_good_bad_TRD_chambers;

    vector<TVector3> vec_TV3_TRD_center_offset; // 540 chambers
    vector< vector<TVector3> >     vec_TV3_TRD_center; // 540 chambers, 3 axes

    Int_t arr_layer_detector[6];

    TH1D* h_delta_angle_perp_impact;
    TH1D* h_detector_hit;
    vector< vector<TH1D*> > vec_h_diff_helix_line_impact_angle; // [all,-,+][540]
    vector< vector<TH1D*> > vec_h_diff_ref_loc_off_trkl; // [all,-,+][540]


    // TRD offline Tracklets
    vector< vector<TVector3> > vec_TV3_Tracklet_pos;
    vector< vector<TVector3> > vec_TV3_Tracklet_dir;



    TGeoManager     *geom;
    TGeoMaterial    *vacuum;
    TGeoMedium      *Air;
    TGeoVolume      *top;
    TGeoMaterial *Fe;
    TGeoMaterial *M_outer_tube;
    TGeoMaterial *Material_TRD_box;
    TGeoMedium *Iron;
    TGeoMedium *Me_outer_tube;
    TGeoMedium *Medium_TRD_box;

    TGeoVolume *inner_field_tube;
    TGeoVolume *outer_field_tube;


    TGeoCombiTrans* combitrans[540];
    TGeoVolume *TRD_boxes[540];

    UShort_t NumTracks_event = 0;

    Int_t Not_installed_TRD_detectors[19] = {402,403,404,405,406,407,432,433,434,435,436,437,462,463,464,465,466,467,538};
    Int_t Defect_TRD_detectors[84]        = {2,5,17,24,26,30,31,32,36,40,41,43,49,50,59,62,64,78,88,92,93,107,110,111,113,116,119,131,161,
    165,182,184,188,190,191,215,219,221,223,226,227,233,236,239,241,249,255,265,277,287,302,308,310,311,318,319,320,326,328,335,348,354,368,377,380,
    386,389,452,455,456,470,474,476,483,484,485,490,491,493,494,500,502,504,506};

    Double_t max_dca_z_to_track = 8.0; // in cm
    Double_t max_dca_r_to_track = 1.0; // in cm
    vector<Int_t> vec_merge_time_bins;

    vector< vector<TVector3> > vec_TV3_digit_pos_cluster;    // layer, merged time bin
    vector< TVector3> vec_TV3_digit_pos_cluster_t0; // layer, x, y, z
    vector< vector<TH1F*> > th1f_ADC_vs_time;
    //Int_t color_layer[7] = {kOrange+2,kGreen,kBlue,kMagenta,kCyan,kPink+5,kYellow};
    Int_t color_layer[7] = {kGray,kGreen,kBlue,kMagenta,kCyan,kYellow,kOrange};
    Int_t line_width_layer[7] = {3,3,3,3,3,3,2};
    vector<Int_t> vec_layer_in_fit;
    vector< vector< vector<Double_t> > > vec_tracklet_fit_points;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_fit_points;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_points_matched;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_fit_points_matched;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_points_background;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_fit_points_background;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_points_bckg;
    vector< vector< vector< vector<Double_t> > > > vec_self_tracklet_fit_points_bckg;

    vector< vector< vector<Double_t> > > vec_ADC_val; //[i_det][i_trkl][i_timebin]

    vector<TProfile*> vec_tp_Delta_vs_impact;
    vector<TH2D*> vec_TH2D_Delta_vs_impact;

    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
    TH1D* radii_digits_initial;
    TH1D* radii_tracklets_final;

    TFile* layer_radii_file;
    TH1D* h_layer_radii_det;
    TFile* file_TRD_geometry;
    TH2D* h2D_TRD_det_coordinates;

    TVector3 TV3_SVD_tracklet_offset;
    TVector3 TV3_SVD_tracklet_dir;
    //Double_t SVD_chi2;


#if defined(USEEVE)
    TEveLine* TEveLine_beam_axis = NULL;
    TEveLine* TPL3D_helix = NULL;
    vector<TEveLine*> vec_TPL3D_helix;
    vector<TEveLine*> vec_TPL3D_helix_inner;
    vector<TEveLine*> vec_TPL3D_helix_hull;

    vector<TEveLine*> vec_TPL3D_helix_kalman;
    vector<TEveLine*> vec_TPL3D_helix_kalman_inner;
    vector<TEveLine*> vec_TPL3D_helix_kalman_hull;

    vector< vector<TEveLine*> > vec_TEveLine_tracklets;
    vector< vector<TEveLine*> > vec_TEveLine_tracklets_match;
    vector< vector<TEveLine*> > vec_TEveLine_self_matched_tracklets;
    TEvePointSet* TEveP_TRD_det_origin;
    TEvePointSet* TEveP_offset_points;
    TEvePointSet* TEveP_TPC_at_offset_points;
    vector<TEveBox*> vec_eve_TRD_detector_box;
    TEvePointSet* TEveP_sec_vertices;
    TEvePointSet* TEveP_close_TPC_photon;
    TEvePointSet* TEveP_nucl_int_vertices;
    TEvePointSet* TEveP_primary_vertex;
    TEvePointSet* TEveP_first_point_helix;
    TEvePointSet* TEveP_second_point_helix;
    vector<TEveLine*> TEveLine_mother;

    vector<TEvePointSet*> TEveP_digits;
    Double_t arr_ADC_color_digits[8] = {0,10.0,20.0,30.0,40.0,50.0,60.0,200.0};
    Int_t arr_color_ADC_digits[7] = {kGray,kGreen,kCyan,kBlue,kYellow,kOrange,kRed};

    TEvePointSet* TEve_clusters;
    TEvePointSet* TEve_connected_clusters;
    vector< vector<TEveLine*> > TEveLine_fitted_tracklets;
    Double_t scale_length_vec = -10.0;

#endif



public:
    TTRD_ST_Make_Tracklets(Int_t graphics);
    ~TTRD_ST_Make_Tracklets();


    void Init_QA();


    //bool sortcol_first( const vector<Double_t>& v1, const vector<Double_t>& v2 );
    Double_t calculateMinimumDistanceStraightToPoint(TVector3 &base, TVector3 &dir,
                                                     TVector3 &point);
    TVector3 calculateDCA_vec_StraightToPoint(TVector3 &base, TVector3 &dir, TVector3 &point);
    TVector3 calculate_point_on_Straight_dca_to_Point(TVector3 &base, TVector3 &dir, TVector3 &point);
    TVector3 calculate_point_on_Straight_dca_to_Point_2D(TVector3 &base, TVector3 &dir, TVector3 &point);

    void Reset();

    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t event);
    //void Calc_SVD_tracklet(Int_t i_det, Int_t i_trkl);
    Double_t Calc_SVD_tracklet(Int_t i_det, Int_t i_trkl);
    void Make_clusters_and_get_tracklets_fit(Double_t Delta_x, Double_t Delta_z, Double_t factor_missing, Int_t graphics);
    Int_t Calibrate(Double_t Delta_x, Double_t Delta_z, Double_t factor_missing, Int_t graphics);
	
    void plot_dem_histos1();
    void plot_dem_histos2();


    ClassDef(TTRD_ST_Make_Tracklets, 1)
};
//----------------------------------------------------------------------------------------


#endif // __TTRD_MAKE_TRACKLETS_H__