#ifndef TRD_KALMAN_TRACKING
#define TRD_KALMAN_TRACKING

//#include "TRD_Kalman_Tracking_Source.cxx"

class TRD_Kalman_Trackfinder
{
 private:
  ROOT::Math::SMatrix<double, 5, 5> mCov;      //Covariance Matrix
  ROOT::Math::SMatrix<double, 4, 5> mObs;      //Observation Model Matrix
  ROOT::Math::SMatrix<double, 5, 5> mTau;      //Propagation uncertainty Matrix
  ROOT::Math::SMatrix<double, 4, 4> mSig;      //Measure uncertainty Matrix
  ROOT::Math::SMatrix<double, 5, 4> mKal_Gain; //Kalman Gain Matrix
  ROOT::Math::SMatrix<double, 4, 4> mCov_Res_Inv;   
  //Double_t curv;							//curvature
  ROOT::Math::SVector<double, 5> mMu;                       	//current estimate
  vector<ROOT::Math::SVector<double, 5>> mEstimate;   		//estimates of this track
  vector<vector<ROOT::Math::SVector<double, 5>>> mEstimates; 	//estimates of all tracks
  ROOT::Math::SVector<double, 4>* mMeasurements;
  Double_t mChi_2;
  Double_t mDist;
  Int_t mCurrent_Det;
  ROOT::Math::SVector<double, 4> mUnc;
  ROOT::Math::SVector<double, 4> mMu_red;
  ROOT::Math::SVector<double, 4> mRes;

  Double_t b_field	=	0.1;

  vector<vector<Double_t>> mHelices; // Kalman helix parameters, based on AliHelix

  Double_t mTRD_layer_radii[6][2] = {
    { 297.5, 306.5 },
    { 310.0, 320.0 },
    { 323.0, 333.0 },
    { 336.0, 345.5 },
    { 348.0, 357.0 },
    { 361.0, 371.0 }
  };

  vector<vector<Ali_TRD_ST_Tracklets*>> mSeed;
  vector< vector<Bool_t> >  mVisited;
  vector<vector<Ali_TRD_ST_Tracklets*>> mBins; //Bins for all Tracklets corresponding to each Module and Layer
                                               //how many entries there are per bin
  vector<Ali_TRD_ST_Tracklets*> mTrack;        //Array with the Tracklets of the current Track
  Int_t mNbr_tracklets;
  vector<vector<Ali_TRD_ST_Tracklets*>> mFound_tracks; //Array with the Tracklets of all found Tracks

  Bool_t mShow;

  ROOT::Math::SVector<double, 4> measure(Ali_TRD_ST_Tracklets* tracklet);
  Bool_t fitting(Ali_TRD_ST_Tracklets* a, Ali_TRD_ST_Tracklets* b);
  void get_seed(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets);
  void prediction(Double_t dist);
  void correction(ROOT::Math::SVector<double, 4> measure);
  //Bool_t fits(ROOT::Math::SVector<double, 4> measure);
  void Kalman(vector<Ali_TRD_ST_Tracklets*> start);

 public:
  vector<vector<Ali_TRD_ST_Tracklets*>> Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets);
  vector<vector<Double_t>> get_Kalman_helix_params();

};
	
#endif
