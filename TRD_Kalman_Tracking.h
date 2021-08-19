#ifndef TRD_KALMAN_TRACKING
#define TRD_KALMAN_TRACKING

#define USEEVE

#if defined(USEEVE)
#include "TEveBox.h"
#include "TEveArrow.h"
#include <TEveManager.h>
#include "TEveLine.h"
#include "TEvePointSet.h"
#endif

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TString.h"

#include "TObject.h"


#include<TMath.h>

//for generator
#include "TRandom.h"
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
#include "TH3D.h"
#include <Math/BinaryOperators.h>
#include "Math/MConfig.h"
#include <iosfwd>
#include "Math/Expression.h"
#include "Math/MatrixRepresentationsStatic.h"
#include "Math/SMatrix.h"
#include "Math/MatrixFunctions.h"
#include "Ali_TRD_ST.h"
#include "Ali_TRD_ST_LinkDef.h"

ClassImp(Ali_TRD_ST_Tracklets)

//#include "TRD_Kalman_Tracking_Source.cxx"



class TRD_Kalman_Trackfinder
{
 private:
  ROOT::Math::SMatrix<double, 5, 5> mCov;      //Covariance Matrix
  ROOT::Math::SMatrix<double, 4, 5> mObs = ROOT::Math::SMatrixIdentity();      //Observation Model Matrix
  ROOT::Math::SMatrix<double, 5, 4> mObs_t{ROOT::Math::Transpose(mObs)};      //Observation Model Matrix
  ROOT::Math::SMatrix<double, 5, 5> mTau;      //Propagation uncertainty Matrix
  ROOT::Math::SMatrix<double, 4, 4> mSig;      //Measure uncertainty Matrix
  ROOT::Math::SMatrix<double, 5, 4> mKal_Gain; //Kalman Gain Matrix
  ROOT::Math::SMatrix<double, 5, 5> mKal_Gain_RTS; //Kalman Gain Matrix for RTS smoother
  ROOT::Math::SMatrix<double, 4, 4> mCov_Res_Inv;
  ROOT::Math::SMatrix<double, 5, 5> A; //transportation matrix
  ROOT::Math::SMatrix<double, 5, 5> A_t; //transportation matrix transposed
  //Double_t curv;							//curvature
  ROOT::Math::SVector<double, 5> mMu;                       	//current estimate
  std::vector<ROOT::Math::SVector<double, 5>> mEstimate;   		//estimates of this track
  std::vector<std::vector<ROOT::Math::SVector<double, 5>>> mEstimates; 	//estimates of all tracks
  std::vector<ROOT::Math::SMatrix<double, 5, 5>> mCovs; //covariances of this track
  std::vector<ROOT::Math::SMatrix<double, 5, 5>> mCovsTrans; //covariances P_k|k-1 of this track
  std::vector<ROOT::Math::SMatrix<double, 5, 5>> mTrans; //transportation matrices for track 
  std::vector<ROOT::Math::SMatrix<double, 5, 5>> mTrans_t; //transportation matrices for track transposes
  ROOT::Math::SVector<double, 4>* mMeasurements;
  Double_t mChi_2;
  Double_t mDist;
  Int_t mCurrent_Det;
  Int_t mPrimVertex;
  ROOT::Math::SVector<double, 4> mUnc;
  ROOT::Math::SVector<double, 4> mMu_red;
  ROOT::Math::SVector<double, 4> mRes;
  
  TFile *correctfile;
  TTree *correcttree;
  TVector3 track_pos_before_tree;
  TVector3 track_pos_after_tree;
  TVector3 offset_tree;
  TVector3 dir_tree;
  Int_t layer_tree;
  Double_t momentum_tree;
  Double_t chi_sq_tree;
  
  TH1D* h_layer_radii;

  Double_t b_field	=	0.5;
  Double_t mchi_2_pen	=	18.5;

  TRandom r;

  std::vector<std::vector<Double_t>> mHelices; // Kalman helix parameters, based on AliHelix
  vector<Double_t> mChi_2s;          //chi_2 of every Kalman track
  
  
  //needed for intersection points method
  ROOT::Math::SVector<double, 4>* mMeasurements_old;

  Double_t mTRD_layer_radii[6][3] = {
    { 0,297.5, 306.5 },
    { 0,310.0, 320.0 },
    { 0,323.0, 333.0 },
    { 0,336.0, 345.5 },
    { 0,348.0, 357.0 },
    { 0,361.0, 371.0 }
  };

  Double_t mTRD_layer_radii_all[540];
  
  std::vector<std::vector<TVector3>> detector_geometry; //0: offset, 1,2,3: x,y,z direction
  std::vector<std::vector<std::vector<Double_t>>> pad_geometry; //detector, row, pos/size

  std::vector<std::vector<Ali_TRD_ST_Tracklets*>> mSeed;
  std::vector< std::vector<Bool_t> >  mVisited;
  std::vector<std::vector<Ali_TRD_ST_Tracklets*>> mBins; //Bins for all Tracklets corresponding to each Module and Layer
                                               //how many entries there are per bin
  std::vector<Ali_TRD_ST_Tracklets*> mTrack;        //Array with the Tracklets of the current Track
  Int_t mNbr_tracklets;
  std::vector<std::vector<Ali_TRD_ST_Tracklets*>> mFound_tracks; //Array with the Tracklets of all found Tracks
  std::vector<std::vector<Ali_TRD_ST_Tracklets*>> mFound_tracks_corr; //Array with the tracklets of all found tracks, with correct z offset

  Bool_t mShow;
  Bool_t mSearch_tracklets;
  ROOT::Math::SVector<double, 4> measure(Ali_TRD_ST_Tracklets* tracklet);
  Bool_t fitting(Ali_TRD_ST_Tracklets* a, Ali_TRD_ST_Tracklets* b);
  void get_seed(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets);
  Bool_t prediction(Double_t dist);
  void correction(ROOT::Math::SVector<double, 4> measure);
  //Bool_t fits(ROOT::Math::SVector<double, 4> measure);
  void Kalman(std::vector<Ali_TRD_ST_Tracklets*> start);
  Bool_t transport(ROOT::Math::SVector<double, 4> &meas, ROOT::Math::SVector<double, 5> mu, Double_t dist);
  Bool_t transport(ROOT::Math::SVector<double, 5> &mu, Double_t dist);
  std::vector<double> find_Intersect(vector<double> pt1, vector<double> dir1, double alpha1, vector<double> pt2, vector<double> dir2, double alpha2);
  TVector2 find_Intersect(TVector2 pt1, TVector2 dir1, TVector2 pt2, TVector2 dir2);
  //void Chi2_Kalman_tracklet(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t );
  Int_t GetPadRowNumber(Int_t det, Double_t z);

 public:
  void Initialize_Tree();
  void Write();
  //std::vector< std::vector<Ali_TRD_ST_Tracklets*> > Kalman_Trackfit(std::vector< std::vector<Ali_TRD_ST_Tracklets*> > tracks,Int_t prim_vertex);
  std::pair<std::vector< std::vector<Ali_TRD_ST_Tracklets*> >,std::vector< std::vector<Ali_TRD_ST_Tracklets*> >> Kalman_Trackfit(std::vector< std::vector<Ali_TRD_ST_Tracklets*> > tracks,Int_t prim_vertex);
  //std::vector<std::vector<Ali_TRD_ST_Tracklets*>> Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets, Int_t prim_vertex);
  std::pair<std::vector<std::vector<Ali_TRD_ST_Tracklets*>>,std::vector<std::vector<Ali_TRD_ST_Tracklets*>>> Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets, Int_t prim_vertex);
  std::vector<std::vector<Double_t>> get_Kalman_helix_params();
  vector<Double_t> get_Kalman_chi_2();
  void set_layer_radii_hist(TH1D* h_layer_radii_in)
  {
      h_layer_radii = h_layer_radii_in;
      for(Int_t i_det = 0; i_det < 540; i_det++)
      {
          mTRD_layer_radii_all[i_det] = h_layer_radii ->GetBinContent(i_det+1);
		  //cout<<"Det:"<<i_det<<"Det_lay: "<<i_det%6<<" radii: "<<mTRD_layer_radii_all[i_det]<<endl;
      }
  }
  
    void set_detector_geom()
  {
		TFile* file_TRD_geom_new_align = TFile::Open("./TRD_geometry_full_new_align.root");

		detector_geometry.resize(4);
		//detector_offsets.resize(540);
		//detector_x_axis.resize(540);
		//detector_y_axis.resize(540);
		//detector_z_axis.resize(540);

		for(Int_t i_vec = 0; i_vec < 4; i_vec++){
			detector_geometry[i_vec].resize(540);
			for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++){
				//HistName = "vec_TH1D_TV3_TRD_coordinates_";
				//HistName += i_vec;
				//HistName += "_V";
				//HistName += i_xyz;
				TH1D *hist_temp;
				if(i_vec==1){
					hist_temp = (TH1D*)file_TRD_geom_new_align->Get(Form("vec_TH1D_TV3_TRD_coordinates_%i_V%i", 2, i_xyz));
				}else if(i_vec==2){
					hist_temp = (TH1D*)file_TRD_geom_new_align->Get(Form("vec_TH1D_TV3_TRD_coordinates_%i_V%i", 1, i_xyz));
				}else{
					hist_temp = (TH1D*)file_TRD_geom_new_align->Get(Form("vec_TH1D_TV3_TRD_coordinates_%i_V%i", i_vec, i_xyz));
				}
				
				for(Int_t i_det = 0; i_det < 540; i_det++){
						detector_geometry[i_vec][i_det][i_xyz] = hist_temp->GetBinContent(i_det+1);
					/*if(i_vec==0){
						detector_offsets[i_det][i_xyz] = hist_temp->GetBinContent(i_det+1);
					}else if(i_vec==1){
						//detector_x_axis[i_det][i_xyz] = hist_temp->GetBinContent(i_det+1);
					}else if(i_vec==2){
						detector_y_axis[i_det][i_xyz] = hist_temp->GetBinContent(i_det+1);
					}else if(i_vec==3){
						detector_z_axis[i_det][i_xyz] = hist_temp->GetBinContent(i_det+1);
					}*/
				}
				
			}

		}

	};
	
	TVector3 get_detector_geom(Int_t i_det, Int_t i_vec)
	{
		return detector_geometry[i_vec][i_det];
	}
	
	void set_pad_geom(){
		TFile* file_pad_rows = new TFile("./PadGeometry.root");
		TH3D* hist_rowsize = (TH3D*)file_pad_rows->Get("row_size");
		TH3D* hist_rowpos = (TH3D*)file_pad_rows->Get("pad_centre_z");
		
		pad_geometry.resize(540);
		for(Int_t i_det=0; i_det < 540; i_det++){
			pad_geometry[i_det].resize(16);
			
			for(Int_t i_row = 0; i_row < 16; i_row++){
				pad_geometry[i_det][i_row].resize(2);
				pad_geometry[i_det][i_row][0] = hist_rowpos->GetBinContent(i_det + 1, 1, i_row + 1);
				pad_geometry[i_det][i_row][1] = hist_rowsize->GetBinContent(i_det + 1 , 1, i_row + 1);
			}
		}
	}
	
	Double_t get_pad_geom(Int_t i_det, Int_t i_row, Int_t i_vec){
		return pad_geometry[i_det][i_row][i_vec];
	}

};

#endif
