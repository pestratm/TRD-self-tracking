
#define USEEVE

using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TString.h"

#include "TObject.h"

#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "Math/Functor.h"

#include<TMath.h>

//for generator
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include <Math/BinaryOperators.h>
#include "Math/MConfig.h"
#include <iosfwd>
#include "Math/Expression.h"
#include "Math/MatrixRepresentationsStatic.h"
#include "Math/SMatrix.h"
#include "Math/MatrixFunctions.h"
#include "TProfile.h"
#include "TProfile2D.h"

#if defined(USEEVE)
#include "TEveBox.h"
#include "TEveArrow.h"
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEvePointSet.h"
#endif

#include "Ali_TRD_ST.h"
#include "Ali_TRD_ST_LinkDef.h"
#include "Ali_TRD_Self_Event.h" 
#include "Ali_TRD_Self_EventLinkDef.h" 


ClassImp(Ali_TRD_ST_Tracklets)
ClassImp(Ali_TRD_ST_TPC_Track)
ClassImp(Ali_TRD_ST_Event)


ClassImp(Ali_Kalman_Track)
ClassImp(Ali_TPC_Track)
ClassImp(Ali_TRD_Photon)
ClassImp(Ali_TRD_Nuclear_interaction)
ClassImp(Ali_TRD_Self_Event)



#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

    static vector< vector<Double_t> > vec_Dt_digit_pos_cluster;

//----------------------------------------------------------------------------------------
class Ali_TRD_ST_Analyze
{
private:
    TChain* input_SE;

    TString in_list_name;

    TString TRD_ST_TREE   = "Tree_TRD_ST_Event";
    TString TRD_ST_BRANCH = "Tree_TRD_ST_Event_branch";
    Long64_t file_entries_total;
    Long64_t N_Events;
    Ali_TRD_ST_Tracklets* TRD_ST_Tracklet;
    Ali_TRD_ST_TPC_Track* TRD_ST_TPC_Track;
    Ali_TRD_ST_Event*     TRD_ST_Event;

    Ali_TRD_ST_Tracklets* TRD_ST_Tracklet_out;
    Ali_TRD_ST_TPC_Track* TRD_ST_TPC_Track_out;
    Ali_TRD_ST_Event*     TRD_ST_Event_out;
    TTree* Tree_TRD_ST_Event_out;

    //stuff for new classes 
    Ali_Kalman_Track* TRD_Kalman_track;
    Ali_TPC_Track* TPC_track;
    Ali_TRD_Self_Event*   TRD_Self_Event;
    Ali_TRD_Photon* TRD_Photon;
    Ali_TRD_Nuclear_interaction* TRD_Nuclear_interaction;

    //Ali_Kalman_Track* TRD_Kalman_track_out;
    //Ali_TPC_Track* TPC_track_out;
    Ali_TRD_Self_Event*   TRD_Self_Event_out;
    Ali_TRD_Photon* TRD_Photon_out;
    Ali_TRD_Nuclear_interaction* TRD_Nuclear_interaction_out;

    
    TTree* Tree_TRD_Self_Event_out;


    TH2D* TH2D_AP_plot;
    const Int_t N_AP_radii = 20;
    const Int_t N_pT_resolution = 7;
    const Double_t Delta_AP_radius = 18.0;
    vector<TH2D*> vec_TH2D_AP_plot_radius;
    TH2D* TH2D_pT_TPC_vs_Kalman;
    vector<TH2D*> vec_TH2D_pT_TPC_vs_Kalman;
    vector<TH2D*> vec_TH2D_one_over_pT_TPC_vs_Kalman;

    TFile* outputfile;
    TFile* out_gain;
    Double_t test;

    TString HistName;
    TH1D* th1d_TRD_layer_radii;
    vector<TH1D*> vec_th1d_TRD_layer_radii_det;
    TH1D* th1d_offset_diff;
    TH1D* th1d_angle_diff;

    Double_t EventVertexX = -999.0;
    Double_t EventVertexY = -999.0;
    Double_t EventVertexZ = -999.0;
    Long64_t Global_Event = -999;
    Int_t    Global_RunID = -999;
    TVector3 TV3_EventVertex;

    vector<TVector3> vec_TV3_secondary_vertices;
    TNtuple* NT_secondary_vertices;
    TNtuple* NT_secondary_vertex_cluster;



    //static vector< vector< vector<Double_t> > > vec_Dt_digit_pos_cluster;    // layer, merged time bin. xyzADC for circle fits in Calibrate()
    //vector< vector<Double_t> > vec_Dt_digit_pos_cluster;


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
    TEvePointSet* TEveP_photon_vertices;
    TEvePointSet* TEveP_close_TPC_photon;
    TEvePointSet* TEveP_nucl_int_vertices;
    TEvePointSet* TEveP_primary_vertex;
    TEvePointSet* TEveP_first_point_helix;
    TEvePointSet* TEveP_second_point_helix;
    vector<TEveLine*> TEveLine_mother;
    TEvePointSet* TEveP_beamA;
    TEvePointSet* TEveP_beamB;

    vector< vector<TEveLine*> > TEveLine_vec_dir_vec_circle;
    vector< vector<TEveLine*> > TEveLine_vec_dir_vec_circle_circle;
    vector<TEveLine*> TEveLine_circle;

#endif

    Int_t N_tracklets_layers[6] = {0};
    Double_t scale_length_vec = -10.0;
    Int_t track_color    = kAzure-2;
    Int_t color_layer_match[6] = {kRed,kGreen,kCyan,kYellow,kPink-3,kOrange+8};
    Int_t color_layer[6] = {kGray,kGray,kGray,kGray,kGray,kGray};
    Double_t TRD_layer_radii[6][2] =
    {
        {297.5,306.5},
        {310.0,320.0},
        {323.0,333.0},
        {336.0,345.5},
        {348.0,357.0},
        {361.0,371.0}
    };

    Int_t Not_installed_TRD_detectors[19] = {402,403,404,405,406,407,432,433,434,435,436,437,462,463,464,465,466,467,538};
    TH1I* h_good_bad_TRD_chambers;

    // TRD 3D graphics
    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector

    vector<vector<Double_t>> mHelices_kalman; // Kalman helix parameters, based on AliHelix
    vector<Double_t> mChi_2s_kalman; // Kalman helix parameters, based on AliHelix
    vector<vector<Double_t>> mHelices_TPC; // Kalman helix parameters, based on AliHelix
    vector<vector<Double_t>> PID_params_TPC;
    Double_t aliHelix_params[6];
    vector<Ali_Helix*> vec_helices;
    vector<Ali_Helix*> vec_helices_TRD;
    vector<Ali_Helix*> vec_helices_TPC;
    Ali_Helix* TPC_single_helix;
    vector< vector<Ali_TRD_ST_Tracklets*> > vec_kalman_TRD_trackets;

    TFile* layer_radii_file;
    TH1D* h_layer_radii_det;

    TFile* file_TRD_geometry;
    TH2D* h2D_TRD_det_coordinates;

    vector< vector<TH2D*> > vec_h2D_pT_vs_TPC_TRD_residuals;
    TString input_dir;
    TString input_dir_lists;
    TProfile* tp_efficiency_matching_vs_pT;
    TProfile* tp_efficiency_all_vs_pT;
    vector<TH2D*> vec_h2D_delta_pT_all_vs_pT;

    vector<TProfile*> vec_tp_Delta_vs_impact;
    vector<TH2D*> vec_TH2D_Delta_vs_impact;
    vector<TProfile*> vec_tp_Delta_vs_impact_circle;
    vector<TH2D*> vec_TH2D_Delta_vs_impact_circle;
    vector< vector<TVector3*> >   vec_TV3_TRD_center;
    vector <TVector3*> vec_TV3_TRD_offset;

    TProfile* vec_tp_pvdca_vs_sector;
    TH2D* vec_TH2D_pvdca_vs_sector;
    //Int_t inits_vec_tp_pvdca_vs_sector = 0;

    vector<vector<Double_t>> ADC_rel;
    vector<TGraph*> tg_bethe_bloch; // [e mu pi K p]
    TProfile* tp_gain;
    TProfile* tp_gain_uncor;
    TH2D* th2d_gain;
    TProfile* tp_dEdx_length_det;
    TH1D* h_dEdx_length_det_orig;
    TH1D* h_dEdx_length_det;
    vector<TProfile2D*> vec_tp2d_gain_vs_xz;

    vector <vector <Double_t>> vec_trd_TRD_pp_geom;

    //TFile* out_gain;

public:
    Ali_TRD_ST_Analyze(TString out_dir, TString out_file_name, Int_t graphics);
    //~Ali_TRD_ST_Analyze();

    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t i_event, Int_t graphics);
    Int_t Draw_event(Long64_t i_event, Int_t graphics, Int_t draw_tracks, Int_t draw_tracklets, Double_t track_path);
    void Animate_beams(Double_t beam_path);
    Int_t Do_TPC_TRD_matching(Long64_t i_event, Double_t xy_matching_window, Double_t z_matching_window, Int_t graphics);
    void Draw_hist_TPC_tracklet_diffs();
    TH1I* get_h_good_bad_TRD_chambers();

    void set_self_event_info();

    Ali_TRD_ST_Tracklets** Tracklets;
    vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks;
    vector<Int_t> vec_idx_matched_TPC_track;
    Int_t Number_Tracklets;
    void Draw_Kalman_Tracklets(vector< vector<Ali_TRD_ST_Tracklets*> > found_tracks);
    void Draw_matched_Kalman_Tracklets(Int_t i_track_plot);
    void set_Kalman_helix_params(vector<vector<Double_t>> mHelices_kalman_in);
    void set_Kalman_chi_2(vector<Double_t> mChi_2s_in);
    Int_t set_TPC_helix_params(Long64_t i_event);
    void set_Kalman_TRD_tracklets(vector< vector<Ali_TRD_ST_Tracklets*> > vec_kalman_TRD_trackets_in);
    void Calc_Kalman_efficiency();
    void Match_kalman_tracks_to_TPC_tracks(Int_t graphics, Int_t draw_matched_TPC_track, Int_t draw_matched_TRD_track, Int_t color);
    void Draw_Kalman_Helix_Tracks(Int_t n_track, Int_t color, Double_t low_R, Double_t high_R);
    void set_single_helix_params(vector<Double_t> vec_params);
    void Evaluate(Double_t t, // helix evaluation, taken from AliHelix
                  Double_t r[3]);  //radius vector
    void fHelixAtoPointdca(TVector3 space_vec, Ali_Helix* helixA, Float_t &pathA, Float_t &dcaAB);
    void fHelixABdca(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB,Float_t pathA_in, Float_t pathB_in);
    Int_t fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c);
    TVector3 fDCA_Helix_Estimate(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB);
    Int_t Calculate_secondary_vertices(Int_t graphics, Int_t flag_TRD_TPC_tracks, Int_t flag_fill_tree);
    pair<Double_t,Double_t>fpathLength(Double_t r,Ali_Helix* helixA) const;
    Int_t fCircle_Interception(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c);
    void Plot_AP();
    void Plot_pT_TPC_vs_Kalman();
    void Draw_TPC_track(Int_t i_track, Int_t color, Double_t line_width, Double_t max_path);
    TH1D* get_layer_radii_hist() {return h_layer_radii_det;}
    Long64_t get_N_Events() {return N_Events;}
    void set_input_lists(TString input_dir_lists_in) {input_dir_lists = input_dir_lists_in;}
    void set_input_dir(TString input_dir_in) {input_dir = input_dir_in;}
    Float_t Calc_nuclev_bitmap(vector<Int_t> vec_idx_kalman_tracks_nuclev_in);
    void Write();
    void Calibrate(Int_t graphics);
    void Draw_n_Save_Calibration(TString out_dir, TString out_file_name_calib);
    void hists_pv_dca();
    void Draw_n_Save_hists_pv_dca(TString out_dir, TString out_file_name_calib);
    Float_t primary_vertex_dca(Int_t i_track);
    //static Double_t distance_circ_point_2D(Double_t x,Double_t y,Double_t *p);
    //static void sum_distance_circ_point_2D(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t );
    TProfile* Calibrate_gain(Long64_t i_event, Int_t Bethe_flag);
    void Draw_Save_Gain_calib();

    ClassDef(Ali_TRD_ST_Analyze, 1)
};
//----------------------------------------------------------------------------------------


