
//#define USEEVE

using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TString.h"

#include "TObject.h"


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

#if defined(USEEVE)
#include "TEveBox.h"
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEvePointSet.h"
#endif

#include "Ali_TRD_ST.h"
#include "Ali_TRD_ST_LinkDef.h"

ClassImp(Ali_TRD_ST_Tracklets)
ClassImp(Ali_TRD_ST_TPC_Track)
ClassImp(Ali_TRD_ST_Event)


#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

//----------------------------------------------------------------------------------------
class Ali_TRD_ST_Analyze
{
private:
    TChain* input_SE;

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

    TH2D* TH2D_AP_plot;
    TH2D* TH2D_pT_TPC_vs_Kalman;

    TFile* outputfile;
    Double_t test;

    TString HistName;
    TH1D* th1d_TRD_layer_radii;
    vector<TH1D*> vec_th1d_TRD_layer_radii_det;
    TH1D* th1d_offset_diff;
    TH1D* th1d_angle_diff;

    Double_t EventVertexX = -999.0;
    Double_t EventVertexY = -999.0;
    Double_t EventVertexZ = -999.0;
    TVector3 TV3_EventVertex;

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
    TEvePointSet* TEveP_offset_points;
    TEvePointSet* TEveP_TPC_at_offset_points;
    vector<TEveBox*> vec_eve_TRD_detector_box;
    TEvePointSet* TEveP_sec_vertices;
    TEvePointSet* TEveP_primary_vertex;
    TEvePointSet* TEveP_first_point_helix;
    TEvePointSet* TEveP_second_point_helix;
    vector<TEveLine*> TEveLine_mother;
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
    Double_t aliHelix_params[6];
    vector<Ali_Helix*> vec_helices;
    Ali_Helix* TPC_single_helix;
    vector< vector<Ali_TRD_ST_Tracklets*> > vec_kalman_TRD_trackets;

    TFile* layer_radii_file;
    TH1D* h_layer_radii_det;

    TFile* file_TRD_geometry;
    TH2D* h2D_TRD_det_coordinates;

    vector< vector<TH2D*> > vec_h2D_pT_vs_TPC_TRD_residuals;
    TString input_dir;


public:
    Ali_TRD_ST_Analyze(TString out_dir, TString out_file_name, Int_t graphics);
    //~Ali_TRD_ST_Analyze();

    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t i_event, Int_t graphics);
    Int_t Draw_event(Long64_t i_event);
    Int_t Do_TPC_TRD_matching(Long64_t i_event, Double_t xy_matching_window, Double_t z_matching_window, Int_t graphics);
    void Draw_hist_TPC_tracklet_diffs();
    TH1I* get_h_good_bad_TRD_chambers();

    Ali_TRD_ST_Tracklets** Tracklets;
    vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks;
    vector<Int_t> vec_idx_matched_TPC_track;
    Int_t Number_Tracklets;
    void Draw_Kalman_Tracks(vector< vector<Ali_TRD_ST_Tracklets*> > found_tracks);
    void set_Kalman_helix_params(vector<vector<Double_t>> mHelices_kalman_in);
    void set_Kalman_TRD_tracklets(vector< vector<Ali_TRD_ST_Tracklets*> > vec_kalman_TRD_trackets_in);
    void Match_kalman_tracks_to_TPC_tracks(Int_t graphics);
    void Draw_Kalman_Helix_Tracks(Int_t n_track, Int_t color);
    void set_single_helix_params(vector<Double_t> vec_params);
    void Evaluate(Double_t t, // helix evaluation, taken from AliHelix
                  Double_t r[3]);  //radius vector
    void fHelixAtoPointdca(TVector3 space_vec, Ali_Helix* helixA, Float_t &pathA, Float_t &dcaAB);
    void fHelixABdca(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB,Float_t pathA_in, Float_t pathB_in);
    Int_t fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c);
    Int_t fDCA_Helix_Estimate(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB);
    void Calculate_secondary_vertices(Int_t graphics);
    pair<Double_t,Double_t>fpathLength(Double_t r,Ali_Helix* helixA) const;
    Int_t fCircle_Interception(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c);
    void Plot_AP();
    void Plot_pT_TPC_vs_Kalman();
    void Draw_TPC_track(Int_t i_track, Int_t color, Double_t line_width);
    TH1D* get_layer_radii_hist() {return h_layer_radii_det;}
    Long64_t get_N_Events() {return N_Events;}
    void create_output_file(TString out_dir, TString out_file_name);
    void set_input_dir(TString input_dir_in) {input_dir = input_dir_in;}
    void Write();


    ClassDef(Ali_TRD_ST_Analyze, 1)
};
//----------------------------------------------------------------------------------------


