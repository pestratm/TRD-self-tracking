#include "TRD_Kalman_Tracking.h"
#include "TVector3.h"
#include "TVirtualFitter.h"
#include "Math/MinimizerOptions.h"
#include "TMinuitMinimizer.h"
#include <chrono>
using namespace std;
using namespace ROOT::Math;

double T_calc; //Variable for the calculation time
double T_wall;

// Colors for printf statements
#define KNRM "\x1B[0m"
#define KRED "\x1B[31m"
#define KGRN "\x1B[32m"
#define KYEL "\x1B[33m"
#define KBLU "\x1B[34m"
#define KMAG "\x1B[35m"
#define KCYN "\x1B[36m"
#define KWHT "\x1B[37m"

static vector<Ali_TRD_ST_Tracklets*> track;
static vector<vector<TVector3>> det_geom;
void Chi2_Kalman_tracklet(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t );

/*static Ali_Helix *kalmanHelix;
static TVector3 tracklet_offset;
static TVector3 det_dir;
void dist_helix_tracklet(Int_t &, Double_t *, Double_t & dist, Double_t * par, Int_t );*/

vector<vector<Double_t>> TRD_Kalman_Trackfinder::get_Kalman_helix_params()
{
  return mHelices;
}
vector<Double_t> TRD_Kalman_Trackfinder::get_Kalman_chi_2()
{
  return mChi_2s;
}

//gives the 4 measurements of a tracklet (y,z,sin(phi),tan(lambda)) back
SVector<double, 4> TRD_Kalman_Trackfinder::measure(Ali_TRD_ST_Tracklets* tracklet)
{
  SVector<double, 4> measurement;
  TVector3 offset = tracklet->get_TV3_offset();
  TVector3 dir = tracklet->get_TV3_dir();
  Double_t inv_hypo = 1. / (TMath::Sqrt(dir[1] * dir[1] + dir[0] * dir[0]));
  measurement[0] = offset[1];
  measurement[1] = offset[2];
  measurement[2] = dir[1] * inv_hypo;
  measurement[3] = dir[2] * inv_hypo;
  return measurement;
}

//checks if 2 tracklets are more or less in a line
Bool_t TRD_Kalman_Trackfinder::fitting(Ali_TRD_ST_Tracklets* a, Ali_TRD_ST_Tracklets* b)
{
  //direction is not ok
  //if ((a->get_TRD_index()==1400)&&(b->get_TRD_index()==389)) cout<<"a2";

  TVector3 dir_a = a->get_TV3_dir();
  TVector3 dir_b = b->get_TV3_dir();
  TVector3 a_angle = (TVector3){ dir_a[0], dir_a[1], 0 };
  TVector3 b_angle = (TVector3){ dir_b[0], dir_b[1], 0 };
  /*
	if ((a->get_TRD_index()==1400)&&(b->get_TRD_index()==389)){
        */
  /*
        cout<<"a ";
		dir_a.Print();
		cout<<" b ";
		dir_b.Print();
                cout<<abs(a_angle.Angle(b_angle)*TMath::RadToDeg())<<" "<<abs(dir_a.Angle(dir_b)*TMath::RadToDeg() )<<endl;
                */
  /*}
	*/
  if (abs(a_angle.Angle(b_angle) * TMath::RadToDeg()) > 15)
    return 0;
  if (abs(dir_a.Angle(dir_b) * TMath::RadToDeg()) > 15)
    return 0;

  //position is not ok
  TVector3 off_b = b->get_TV3_offset();
  TVector3 off_a = a->get_TV3_offset();
  TVector3 z = off_a + dir_a * ((off_b[0] - off_a[0]) / dir_a[0]);
  /*
	if ((a->get_TRD_index()==1400)&&(b->get_TRD_index()==389)){
		cout<<"a ";
		z.Print();
		cout<<" b ";
		off_b.Print();
	}*/

  if (abs((z - off_b)[1]) > 7.)
    return 0;
  //if ((a->get_TRD_index()==1648) && (b->get_TRD_index()==1643)) cout<<"a"<<endl;

  if (abs((z - off_b)[2]) > 20)
    return 0;
  return 1;
}

void TRD_Kalman_Trackfinder::get_seed(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets)
{

  Int_t nbr_sectors = 18;
  Int_t nbr_stack = 5;
  Int_t nbr_layers = 6;
  mBins.clear();
  mBins.resize(nbr_sectors * nbr_stack * nbr_layers); //bins for every detector
  mVisited.clear();
  mVisited.resize(nbr_sectors * nbr_stack * nbr_layers);

  //sort then in bins
  for (Int_t i_tracklet = 0; i_tracklet < Num_Tracklets; i_tracklet++) {
    if (Tracklets[i_tracklet]->get_TV3_offset().Mag() > 1000.0)
      continue;

    Int_t i_det = Tracklets[i_tracklet]->get_TRD_det();
    Int_t i_sector = (Int_t)(i_det / 30);

    //rotate the tracklets to the local coordiante system
    TVector3 temp1 = Tracklets[i_tracklet]->get_TV3_offset();
    temp1.RotateZ((Double_t)(-1 * (2 * i_sector + 1) * TMath::Pi() / 18));
    Tracklets[i_tracklet]->set_TV3_offset(temp1);

    TVector3 temp2 = Tracklets[i_tracklet]->get_TV3_dir();
    temp2.RotateZ((Double_t)(-1 * (2 * i_sector + 1) * TMath::Pi() / 18));
    Tracklets[i_tracklet]->set_TV3_dir(temp2);
    Tracklets[i_tracklet]->set_TRD_index(i_tracklet);

    //sort them in
    mBins[i_det].push_back(Tracklets[i_tracklet]);
    mVisited[i_det].push_back(0);
  }

  //check for each tracklet in the 6. and 5. layer whether there are matching tracklets in the 4. or 5. layer
  Bool_t done = 0;
  for (int i_sec = 0; i_sec < nbr_sectors; i_sec++)
    for (int i_stck = 0; i_stck < nbr_stack; i_stck++)
      for (int i_lay = nbr_layers - 1; i_lay > 1; i_lay--) {
        Int_t i_det = i_sec * nbr_stack * nbr_layers + i_stck * nbr_layers + i_lay;
        for (int i_tracklet = 0; i_tracklet < (Int_t)mBins[i_det].size(); i_tracklet++)
          if (!(mVisited[i_det][i_tracklet])) {
            mVisited[i_det][i_tracklet] = 1;
            for (int j_lay = i_lay - 1; j_lay > i_lay - 3; j_lay--) {
              done = 0;
              Int_t j_det = i_det + j_lay - i_lay;
              for (int j_tracklet = 0; j_tracklet < (Int_t)mBins[j_det].size(); j_tracklet++)
                if (fitting(mBins[i_det][i_tracklet], mBins[j_det][j_tracklet])) {
                  mVisited[j_det][j_tracklet] = 1;
                  vector<Ali_TRD_ST_Tracklets*> temp_seed;
                  temp_seed.resize(2);
                  temp_seed[0] = mBins[i_det][i_tracklet];
                  temp_seed[1] = mBins[j_det][j_tracklet];
                  mSeed.push_back(temp_seed);
                  done = 1;
                  chrono::high_resolution_clock::time_point start_calc_time = chrono::high_resolution_clock::now();

                  Kalman(mSeed[(Int_t)mSeed.size() - 1]);
                  chrono::high_resolution_clock::time_point end_calc_time = chrono::high_resolution_clock::now();
                  T_calc += chrono::duration_cast<chrono::microseconds>(end_calc_time - start_calc_time).count();
                }
              if (done)
                break;
            }
          }
      }
  //for(int i=0;i<mBins.size();i++){
  //cout<<i<<": "<<mBins[i].size()<<endl;
  //}
}

//predict next location
Bool_t TRD_Kalman_Trackfinder::prediction(Double_t dist)
{
  if (abs(dist) < 1e-7)
    return true;

  //calculate new estimate
  Double_t f1 = mMu[2];
  Double_t b_fak = b_field * 3. / 1000.; //3/1000 so that units are correct GV/c *T = 1000/3 cm

  Double_t curvature = (mMu[4] * b_fak);
  Double_t f2 = f1 + dist * curvature;
  Double_t r1 = TMath::Sqrt((1. - f1) * (1. + f1));
  Double_t r2 = TMath::Sqrt((1. - f2) * (1. + f2));
  Double_t dy2dx = (f1 + f2) / (r1 + r2);
  if ((f2 >= 1) || (f2 <= -1))
    return 0;
  if ((f1 >= 1) || (f1 <= -1))
    return 0;
  if (r1 <= 0)
    return 0;
  if (r2 <= 0)
    return 0;

  mMu[0] += dist * dy2dx;
  //if (dist<-300) printf("values: curv: %f f1: %f f2: %f r1 %f r2 %f dy2dx: %f y: %f\n",curvature,f1,f2,r1,r2,dy2dx,mMu[0]);
  if (TMath::Abs(dist * curvature) < 0.05) {
    mMu[1] += dist * (r2 + f2 * dy2dx) * mMu[3];
  } else {
    Double_t rot = TMath::ASin(r1 * f2 - r2 * f1);
    if (f1 * f1 + f2 * f2 > 1 && f1 * f2 < 0) { // special cases of large rotations or large abs angles
      if (f2 > 0) {
        rot = TMath::Pi() - rot;
      } else {
        rot = -TMath::Pi() - rot;
      }
    }
    mMu[1] += mMu[3] / curvature * rot;
  }
  mMu[2] = f2;

  //Build transport matrix
  //SMatrix<double, 5, 5> A; //Transport Matrix
  A = SMatrixIdentity();

  Double_t dr = dist / (r1 * r1 * r1);
  A[0][2] = dr;
  A[0][4] = dist * dr * b_fak * 0.5;
  A[1][2] = A[0][2] * mMu[3] * f1;
  A[1][3] = dist / r1;
  A[1][4] = A[1][2] * dist * b_fak * 0.5;
  A[2][4] = dist * b_fak;
  //SMatrix<double, 5, 5> A_t = Transpose(A);
  A_t = Transpose(A);

  if (mShow) {
    //cout<<"Unc:"<<mUnc<<endl;
    cout << "dist: " << dist << endl;
    cout << "curv: " << curvature << endl;
    cout << "r1: " << r1 << endl;
    cout << "(1.- f1)*(1. + f1): " << (1. - f1) * (1. + f1) << endl;
    cout << "f1: " << f1 << endl;
    cout << "f2: " << f2 << endl;
    cout << "transport:" << A << endl;
    cout << "f_mult:" << A * mCov << endl;
  }

  //calculate new covariance
  mCov = A * mCov * A_t;

  mMu_red = mObs * mMu;       //made once and saved for later use
  for (int i = 0; i < 4; i++) //calculate the uncertainty (useful to know)
    mUnc[i] = TMath::Sqrt(mCov[i][i] + mSig[i][i]);
  return 1;
}

//----------------------------------------------------------------------------------------

//predict next locataion 
//similar to prediction, however here the vector that is transported is given to the function explicitely
Bool_t TRD_Kalman_Trackfinder::transport(ROOT::Math::SVector<double, 5> &mu, Double_t dist)
{
	
  /*Double_t b_fak = b_field * 3. / 1000.; //kappa*B
  Double_t r1 = TMath::Sqrt(1 - TMath::Power(mu[2], 2));
  Double_t F02 = dist / TMath::Power(r1, 3);
  Double_t F04 = TMath::Power(dist, 2) / (2 * TMath::Power(r1, 3)) * b_fak;
  Double_t F12 = dist*mu[3]*mu[2]/TMath::Power(r1, 3);
  Double_t F13 = dist/r1;
  Double_t F14 = TMath::Power(dist, 2) * mu[3] * mu[2] * b_fak / (2 * TMath::Power(r1, 3));
  Double_t F24 = dist * b_fak;
  
  Double_t f1 = mu[2];
  Double_t curvature = (mu[4] * b_fak);
  Double_t f2 = f1 + dist * curvature;
  //Double_t r1 = TMath::Sqrt((1. - f1) * (1. + f1));
  Double_t r2 = TMath::Sqrt((1. - f2) * (1. + f2));
  if ((f2 >= 1) || (f2 <= -1))
    return 0;
  if ((f1 >= 1) || (f1 <= -1))
    return 0;
  if (r1 <= 0)
    return 0;
  if (r2 <= 0)
    return 0;
  
  mu[0] += F02 * mu[2] + F04 * mu[4];
  mu[1] += F12 * mu[2] + F13 * mu[3] + F14 * mu[4];
  mu[2] += F24 * mu[4];
  return 1;*/
	
  //calculate new estimate
  Double_t f1 = mu[2];
  Double_t b_fak = b_field * 3. / 1000.; //3/1000 so that units are correct GV/c *T = 1000/3 cm

  Double_t curvature = (mu[4] * b_fak);
  Double_t f2 = f1 + dist * curvature;
  Double_t r1 = TMath::Sqrt((1. - f1) * (1. + f1));
  Double_t r2 = TMath::Sqrt((1. - f2) * (1. + f2));
  Double_t dy2dx = (f1 + f2) / (r1 + r2);
  if ((f2 >= 1) || (f2 <= -1))
    return 0;
  if ((f1 >= 1) || (f1 <= -1))
    return 0;
  if (r1 <= 0)
    return 0;
  if (r2 <= 0)
    return 0;

  mu[0] += dist * dy2dx;
  //if (dist<-300) printf("values: curv: %f f1: %f f2: %f r1 %f r2 %f dy2dx: %f y: %f\n",curvature,f1,f2,r1,r2,dy2dx,mu[0]);
  if (TMath::Abs(dist * curvature) < 0.05) {
    mu[1] += dist * (r2 + f2 * dy2dx) * mu[3];
  } else {
    Double_t rot = TMath::ASin(r1 * f2 - r2 * f1);
    if (f1 * f1 + f2 * f2 > 1 && f1 * f2 < 0) { // special cases of large rotations or large abs angles
      if (f2 > 0) {
        rot = TMath::Pi() - rot;
      } else {
        rot = -TMath::Pi() - rot;
      }
    }
    mu[1] += mu[3] / curvature * rot;
  }
  mu[2] += dist * curvature;
  //return mu;
  return 1;

}

//----------------------------------------------------------------------------------------


//predict next location
//here, a vector meas is transported to another laying using the vector mu
Bool_t TRD_Kalman_Trackfinder::transport(ROOT::Math::SVector<double, 4> &meas, ROOT::Math::SVector<double, 5> mu, Double_t dist)
{
	
  /*Double_t b_fak = b_field * 3. / 1000.; //kappa*B
  Double_t r1 = TMath::Sqrt(1 - TMath::Power(mu[2], 2));
  Double_t F02 = dist / TMath::Power(r1, 3);
  Double_t F04 = TMath::Power(dist, 2) / (2 * TMath::Power(r1, 3)) * b_fak;
  Double_t F12 = dist*mu[3]*mu[2]/TMath::Power(r1, 3);
  Double_t F13 = dist/r1;
  Double_t F14 = TMath::Power(dist, 2) * mu[3] * mu[2] * b_fak / (2 * TMath::Power(r1, 3));
  Double_t F24 = dist * b_fak;
  
  Double_t f1 = mu[2];
  Double_t curvature = (mu[4] * b_fak);
  Double_t f2 = f1 + dist * curvature;
  //Double_t r1 = TMath::Sqrt((1. - f1) * (1. + f1));
  Double_t r2 = TMath::Sqrt((1. - f2) * (1. + f2));
  if ((f2 >= 1) || (f2 <= -1))
    return 0;
  if ((f1 >= 1) || (f1 <= -1))
    return 0;
  if (r1 <= 0)
    return 0;
  if (r2 <= 0)
    return 0;
  
  mu[0] += F02 * mu[2] + F04 * mu[4];
  mu[1] += F12 * mu[2] + F13 * mu[3] + F14 * mu[4];
  mu[2] += F24 * mu[4];
  return 1;*/
	
  //calculate new estimate
  Double_t f1 = mu[2];
  Double_t b_fak = b_field * 3. / 1000.; //3/1000 so that units are correct GV/c *T = 1000/3 cm

  Double_t curvature = (mu[4] * b_fak);
  Double_t f2 = f1 + dist * curvature;
  Double_t r1 = TMath::Sqrt((1. - f1) * (1. + f1));
  Double_t r2 = TMath::Sqrt((1. - f2) * (1. + f2));
  Double_t dy2dx = (f1 + f2) / (r1 + r2);
  if ((f2 >= 1) || (f2 <= -1))
    return 0;
  if ((f1 >= 1) || (f1 <= -1))
    return 0;
  if (r1 <= 0)
    return 0;
  if (r2 <= 0)
    return 0;

  meas[0] += dist * dy2dx;
  //if (dist<-300) printf("values: curv: %f f1: %f f2: %f r1 %f r2 %f dy2dx: %f y: %f\n",curvature,f1,f2,r1,r2,dy2dx,mu[0]);
  if (TMath::Abs(dist * curvature) < 0.05) {
    meas[1] += dist * (r2 + f2 * dy2dx) * mu[3];
  } else {
    Double_t rot = TMath::ASin(r1 * f2 - r2 * f1);
    if (f1 * f1 + f2 * f2 > 1 && f1 * f2 < 0) { // special cases of large rotations or large abs angles
      if (f2 > 0) {
        rot = TMath::Pi() - rot;
      } else {
        rot = -TMath::Pi() - rot;
      }
    }
    meas[1] += mu[3] / curvature * rot;
  }
  meas[2] += dist * curvature;
  //return mu;
  return 1;

}

//----------------------------------------------------------------------------------------


//Find the intersection point of two 2d vectors
TVector2 TRD_Kalman_Trackfinder::find_Intersect(TVector2 pt1, TVector2 dir1, TVector2 pt2, TVector2 dir2){
	TVector2 delta = pt2 - pt1;
	Double_t t2 = (delta.X() * dir1.Y() - delta.Y() * dir1.X())/(dir2.Y() * dir1.X() - dir2.X() * dir1.Y());
	TVector2 intersection_point = pt2 + t2 * dir2;
	return intersection_point;

}

//----------------------------------------------------------------------------------------

//correction as used in every lecture
//see equation 1.9 to 1.14 in Bachelor thesis Sven Hoppner
void TRD_Kalman_Trackfinder::correction(SVector<double, 4> measure)
{

  SMatrix<double, 5, 5> Eye = SMatrixIdentity();

  mKal_Gain = mCov * mObs_t * mCov_Res_Inv;

  SVector<double, 4> res = measure - mObs * mMu;
  mMu += mKal_Gain * (res);

  mCov = (Eye - mKal_Gain * mObs) * mCov;
  res = measure - mObs * mMu;
  mChi_2 += Dot(res, (mCov_Res_Inv * res));
}
	

//----------------------------------------------------------------------------------------
void TRD_Kalman_Trackfinder::Kalman(vector<Ali_TRD_ST_Tracklets*> seed)
{	
	Int_t correct_offsets = 1; //Correct for pad tilting: 0 = don't correct, 1 = find intersection points using original Kalman fit, 2 = fit chi**2 and use to find intersection pts
	Bool_t use_RTS = false; //use Rauch-Tung-Striebel smoother (not fully functional yet)

  { //init
	if(mShow){
		  cout << "--------------------------------------" << endl;
		  cout << "New track" << endl;
		  cout << "--------------------------------------" << endl;
	}
    SVector<double, 4> mes = measure(seed[0]);
    for (int i = 0; i < 4; i++)
      mMu[i] = mes[i];

    mMu[4] = 0.0 / 1.0; // q/pT  0.0/1.0

    mEstimate.resize(0);
    mEstimate.resize(6);
    mTrack.resize(0);
    mTrack.resize(6);
    mCovs.resize(0);
    mCovs.resize(6);
    mCovsTrans.resize(0);
    mCovsTrans.resize(6);
    mTrans.resize(0);
    mTrans.resize(6);
    mTrans_t.resize(0);
    mTrans_t.resize(6);

    mCurrent_Det = seed[0]->get_TRD_det();
    //mDist=mTRD_layer_radii[5][1]-mTRD_layer_radii[mCurrent_Det%6][1];
    mEstimate[mCurrent_Det % 6] = mMu;

    // 0,0 = y, 1,1=z , 2,2=sin phi , 3,3=tan lambda

    Double_t dy2 = 0.2;         // 0.2  0.4
    Double_t dz2 = 4.0;         // 4.0  4.0
    Double_t dsin_phi = 7.0;    // 7.0  10.0
    Double_t dsin_theta = 18.0; // 18.0  25.0
    Double_t dpT = 10.0;        // 10.0  10.0

    mSig[0][0] = dy2;                                                           // 0.2
    mSig[1][1] = dz2;                                                           // 4.0
    mSig[2][2] = TMath::Power(TMath::Sin(dsin_phi * TMath::Pi() / 180.0), 2);   // 7.0
    mSig[3][3] = TMath::Power(TMath::Tan(dsin_theta * TMath::Pi() / 180.0), 2); // 18.0

    mCov = SMatrixIdentity();
    mCov[0][0] = dy2;                                                           // 0.2
    mCov[1][1] = dz2;                                                           // 4.0
    mCov[2][2] = TMath::Power(TMath::Sin(dsin_phi * TMath::Pi() / 180.0), 2);   // 7.0
    mCov[3][3] = TMath::Power(TMath::Tan(dsin_theta * TMath::Pi() / 180.0), 2); // 20.0
                                                                                //mCov[4][4]	=	0.09; // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15
    mCov[4][4] = dpT * dpT;                                                     // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15

    mChi_2 = 0;

    mMeasurements = new SVector<double, 4>[6]; //eventuell fixes array
    mMeasurements_old = new SVector<double, 4>[6]; //only needed for correct_offsets==1

    for (const auto i_seed : seed) { //store measurements from seed and safe tracklets in global coordinate system

      Int_t lay = i_seed->get_TRD_det() % 6;
      mMeasurements[lay] = measure(i_seed);
      TVector3 temp_vec;
      Int_t det;
      mTrack[lay] = new Ali_TRD_ST_Tracklets();

      temp_vec = i_seed->get_TV3_offset();
      det = i_seed->get_TRD_det();
      temp_vec.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18);

      mTrack[lay]->set_TRD_det(det);
      mTrack[lay]->set_TV3_offset(temp_vec);
      temp_vec = i_seed->get_TV3_dir();
      temp_vec.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18);

      mTrack[lay]->set_TV3_dir(temp_vec);
      mTrack[lay]->set_TRD_index(i_seed->get_TRD_index());

      if (mShow) {
        cout << "NEW TRACK" << endl;
        cout << "Layer:" << lay << endl;
        cout << "cov:" << mCov << endl;
        cout << "Mu:" << mMu << endl;
        cout << "mes:" << mMeasurements[lay] << endl;
        //cout<<"unc"<<mUnc<<endl;
      }
    }

    mNbr_tracklets = (Int_t)seed.size();
  }

  vector<SMatrix<double, 5, 5>> cov_per_layer(6);
  cov_per_layer[mCurrent_Det % 6] = mCov;
  mCovs[mCurrent_Det % 6] = mCov;

  for (Int_t i_seed = 1; i_seed < (Int_t)seed.size(); i_seed++) {
    Int_t i_det = seed[i_seed]->get_TRD_det();
    Int_t i_layer = i_det % 6;
    //mDist	=	mTRD_layer_radii[i_layer-1][1]-mTRD_layer_radii[i_layer][1];
    mDist = mTRD_layer_radii_all[i_det] - mTRD_layer_radii_all[mCurrent_Det];

    if (!(prediction(mDist)))
      break;
    mCovsTrans[mCurrent_Det%6] = mCov;
    mTrans[mCurrent_Det%6] = A;
    mTrans_t[mCurrent_Det%6] = A_t;
    mCurrent_Det = i_det;
    mCov_Res_Inv = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
    mCov_Res_Inv.Invert();

    if (mShow) {
      cout << "Layer:" << i_layer << endl;
      cout << "cov:" << mCov << endl;
      cout << "Mu:" << mMu << endl;
      cout << "unc" << mUnc << endl;
      cout << "mes:" << mMeasurements[i_layer] << endl;
    }
    correction(mMeasurements[i_layer]);
    mCovs[i_layer] = mCov;
    cov_per_layer[i_layer] = mCov;
    mEstimate[i_layer] = mMu;
  }
  { //go back to first seeding tracklet
    Int_t i_layer = seed[0]->get_TRD_det() % 6;
    mMu = mEstimate[i_layer];
    mCurrent_Det = mTrack[i_layer]->get_TRD_det();
    mDist = mTRD_layer_radii_all[mCurrent_Det - i_layer + 5] - mTRD_layer_radii_all[mCurrent_Det];
    mCov = cov_per_layer[i_layer];
  }
  //for(Int_t i_layer = 5;i_layer>=0;i_layer--) // Alex: changes where the kalman tracker seeds
  
  //cout << "--------------------------------" << endl;
  //cout << "First Kalman filter" << endl;
  for (Int_t i_layer = mCurrent_Det % 6; i_layer >= 0; i_layer--) // Alex: changes where the kalman tracker seeds
  {                                                               //seeding is done start tracking loop
	//cout << "\nLayer: " << i_layer << endl;

    if (mEstimate[i_layer] != 0) {
      //if we already have an estimate for this layer, update layer and skip calculations
      mMu = mEstimate[i_layer];
      mCurrent_Det = mTrack[i_layer]->get_TRD_det();
      mDist = mTRD_layer_radii_all[mCurrent_Det - 1] - mTRD_layer_radii_all[mCurrent_Det];
      mCov = cov_per_layer[i_layer];
      //cout << "Measurement: " << mMeasurements[i_layer] << endl;
      //cout << "Estimate: " << mEstimate[i_layer] << endl;
      continue;
    }

    if (!(prediction(mDist)))
      break;
    
    //cout << "Fmu prediction: " << A*mEstimate[i_layer+1] << endl;
      
    //cout << "Prediction: " << mMu << endl;
      
    mTrans[i_layer+1] = A;
    mTrans_t[i_layer+1] = A_t;  
    mCovsTrans[i_layer+1] = mCov;

    
    if (mShow) {
      cout << "Layer:" << i_layer << endl;
      cout << "cov:" << mCov << endl;
      cout << "Mu:" << mMu << endl;
      cout << "unc" << mUnc << endl;
    }

    mCurrent_Det--;
    mDist = mTRD_layer_radii_all[mCurrent_Det - 1] - mTRD_layer_radii_all[mCurrent_Det];
    SMatrix<double, 4, 4> Cov_res = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
    Cov_res.Invert();

    if (mMeasurements[i_layer] != 0) {
      printf("ERROR:: should never happen\n");
      mCov_Res_Inv = Cov_res;
      correction(mMeasurements[i_layer]);
    } else if (mSearch_tracklets) {
      vector<Int_t> Dets;
      if (0) {
      } //out of border todo

      else if (0) {
      } //out of border todo

      else //in border
        Dets.push_back(mCurrent_Det);

      for (const auto i_det : Dets) {

        Int_t i_sector = (Int_t)(i_det / 30);
        Int_t i_stack = (Int_t)(i_det % 30 / 6);

        if (i_sector != (Int_t)(mCurrent_Det / 30)) {
        } //change coord sys todo

        //vector<Ali_TRD_ST_Tracklets*> 			found_tracklets;
        //vector<Int_t> 							found_in;
        //vector<Double_t> 						chis;

        Ali_TRD_ST_Tracklets* found_tracklet = nullptr;
        Int_t found_in = -1;
        Int_t ind_nbr = -1;
        Double_t min_chi = mchi_2_pen;

        //vector<SVector<double,4>> 	resses;
        //vector<SVector<double,4>> 	meas_list;

        // TODO SOFORT!!!!!! MIN_CHI UPDATE!!!!!
        for (Int_t i_tracklet = 0; i_tracklet < (Int_t)mBins[i_det].size(); i_tracklet++) {
          SVector<double, 4> measurement = measure(mBins[i_det][i_tracklet]);
          SVector<double, 4> res = measurement - mMu_red;
          Double_t chi_2 = Dot(res, (Cov_res * res));

          if (mShow) {
            cout << "meas" << measurement << ", " << mBins[i_det][i_tracklet]->get_TV3_offset()[0] << ", " << mBins[i_det][i_tracklet]->get_TRD_index() << endl;
            cout << "chi:" << chi_2 << endl;
          }

          if (chi_2 < min_chi) {
            min_chi = chi_2;
            found_tracklet = mBins[i_det][i_tracklet];
            found_in = i_det;
            ind_nbr = i_tracklet;
          }
          //cout<<"dist to rad new:"<<mTRD_layer_radii_all[mCurrent_Det]-mBins[i_det][i_tracklet]->get_TV3_offset()[0]<<" dtor old:"<<mTRD_layer_radii[mCurrent_Det%6][1]-mBins[i_det][i_tracklet]->get_TV3_offset()[0]<<endl;
          //found_tracklets.push_back(mBins[i_det][i_tracklet]);
          //found_in.push_back(i_det);
          //chis.push_back(chi_2);
          //resses.push_back(res);
          //meas_list.push_back(measurement);
          //ind_nbr.push_back(i_tracklet);
        }
        /*
				for(Int_t i=0 ;i<(Int_t)found_tracklets.size();i++)
				{
					if (min_chi>chis[i])
					{
						min_ind=i;
						min_chi=chis[i];
					}
				}*/

		  bool tracklet_in_proximity = false;
		  for (Int_t i_tracklet = 0; i_tracklet < (Int_t)mBins[i_det].size(); i_tracklet++) {
			  if (!(min_chi < mchi_2_pen)) continue;
			  if (found_tracklet == mBins[i_det][i_tracklet]) continue;
			  TVector3 offset1 = found_tracklet->get_TV3_offset();
			  TVector3 offset2 = mBins[i_det][i_tracklet]->get_TV3_offset();
			  if ((offset2 - offset1).Mag() < 2.) tracklet_in_proximity = true;
			  //cout << "dist " << (offset2-offset1).Mag() << endl;
		  }


        if ((min_chi < mchi_2_pen) && !tracklet_in_proximity) { //new tracklet was found
			  
		  //chi_sq_tree = min_chi;
		  //correcttree->Fill();
          if ((Int_t)(found_in / 30) != (Int_t)(mCurrent_Det / 30)) {
          }
          //change_coord_sys(found_in_det[min_ind]); //MuST BE IMPLEMENTET !!!!!!!!

          //save tracklet info
          mVisited[found_in][ind_nbr] = 1;
          mMeasurements[i_layer] = measure(found_tracklet);
          mCov_Res_Inv = Cov_res;

          correction(mMeasurements[i_layer]);

          Ali_TRD_ST_Tracklets* temp_tracklet = new Ali_TRD_ST_Tracklets();
          mTrack[i_layer] = temp_tracklet;

          TVector3 temp_vec;
          Int_t det = found_in;
          mTrack[i_layer]->set_TRD_det(det);

          temp_vec = found_tracklet->get_TV3_offset();
          temp_vec.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18);
          mTrack[i_layer]->set_TV3_offset(temp_vec);

          temp_vec = found_tracklet->get_TV3_dir();
          temp_vec.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18);
          mTrack[i_layer]->set_TV3_dir(temp_vec);

          mTrack[i_layer]->set_TRD_index(found_tracklet->get_TRD_index());
          mNbr_tracklets++;
        } else {
          mTrack[i_layer] = nullptr;
        }
      }
    }
    mEstimate[i_layer] = mMu;
    mCovs[i_layer] = mCov;
    
    //if(mMeasurements[i_layer]!=0) cout << "Measurement: " << mMeasurements[i_layer] << endl;
    //cout << "Estimate: " << mEstimate[i_layer] << endl;
  }
  //if Track
  /*	Double_t Chi_2_list[7]={0.,0.,0.,21.,26.,35.,40.};
    if( (mNbr_tracklets>2)	&&	(mChi_2<=Chi_2_list[mNbr_tracklets]) ){ // Changed from 2, Alex: 20.07.2020
		//save Track
		cout<<"Nbr tracklets: "<<mNbr_tracklets<<" mChi_2: "<<mChi_2<<endl;
*/
  if ((mNbr_tracklets > 2) && (mChi_2 < 100)) { // Changed from 2, Alex: 20.07.2020
    //save Track

	//Have two arrays, containing the tracklets belonging to Kalman tracks, one corrected and one before correction
	//If no correction is performed or only the track but not the tracklets are corrected, those two are identical
    mFound_tracks.push_back(mTrack);
    vector<Ali_TRD_ST_Tracklets*> mTrack_corr;
    if((correct_offsets==1) || (correct_offsets==2)){
		mTrack_corr.resize(mTrack.size());
		mFound_tracks_corr.push_back(mTrack_corr);
	}else{
		mFound_tracks_corr.push_back(mTrack);
	}
    mEstimates.push_back(mEstimate);


    if (mPrimVertex == 1) {
      mDist = -mTRD_layer_radii_all[mCurrent_Det];
      //cout<<"\nTrack NBR: "<< mFound_tracks.size()-1<<endl;
      //cout<<"mDist after pred: "<<mDist<<endl;
      //cout<<"mMu before pred: "<<mMu<<endl;

      if (prediction(mDist)) //TODO new mSig for prim vertex where errors on direction are HUGE
      {
        //cout<<"mMu after pred: "<<mMu<<endl;
        SVector<double, 4> prim_vert_measurement;
        prim_vert_measurement[0] = 0;
        prim_vert_measurement[1] = 0;
        prim_vert_measurement[2] = mMu[2];
        prim_vert_measurement[3] = mMu[3];
        mCov_Res_Inv = mObs * mCov * Transpose(mObs) + 1.0 * mSig; //Measure uncertainty Matrix
        mCov_Res_Inv.Invert();

        correction(prim_vert_measurement);
        Double_t tempqpt = mMu[4];
        mMu = mEstimate[0];
        mMu[4] = tempqpt;
      }
    }

    //RTS Smoother as seen in https://jwmi.github.io/ASM/6-KalmanFilter.pdf
    //For the case of mTau=0, this is equivalent to transporting the estimate of the last layer back to every other layer without predicting.
	if(use_RTS){
		for (Int_t i_layer = mCurrent_Det % 6 +1; i_layer < 6; i_layer++){
			ROOT::Math::SMatrix<double, 5, 5> mCovsTrans_inv = mCovsTrans[i_layer];
			mCovsTrans_inv.Invert();
			mKal_Gain_RTS = mCovs[i_layer]*(mTrans_t[i_layer]*mCovsTrans_inv);
			//mCovsTrans[i_layer] = mTrans[i_layer] * (mCovs[i_layer] * Transpose(mTrans[i_layer]));
			//cout << "i_layer" << i_layer << endl;
			//cout << "old " <<  mEstimate[i_layer] << endl;
			mMu = mEstimate[i_layer];
			double dist = mTRD_layer_radii[i_layer-1][1]-mTRD_layer_radii[i_layer][1];
			ROOT::Math::SMatrix<double, 5, 5> mTrans_inv = mTrans[i_layer];
			mTrans_inv.Invert();
			ROOT::Math::SVector<double, 5> mEstimateTrans = mEstimate[i_layer];
			transport(mEstimateTrans, dist);
			mEstimate[i_layer] += mKal_Gain_RTS*(mEstimate[i_layer-1] - mEstimateTrans);
			//mEstimate[i_layer] = mEstimate[i_layer] + mKal_Gain_RTS*(mEstimate[i_layer-1] - mTrans[i_layer]*mEstimate[i_layer]);
			cout << "new " << mEstimate[i_layer] << endl;
			mCovs[i_layer] = mCovs[i_layer] + mKal_Gain_RTS*((mCovs[i_layer-1] - mCovsTrans[i_layer])*Transpose(mKal_Gain_RTS));
		}
	}else{
    
		//Loop again for better fit
		//cout << "----------------------" << endl;
		//cout << "Second Kalman filter" << endl;
		for (Int_t i_layer = mCurrent_Det % 6 + 1; i_layer < 6; i_layer++) {
			
		 //cout << "\nLayer: " << i_layer << endl;	
		  if (mShow){
			cout << "currently 1st looplayer: " << i_layer << endl;
		   }
		  mDist = mTRD_layer_radii_all[mCurrent_Det + 1] - mTRD_layer_radii_all[mCurrent_Det];
		  //printf("%s ---> i_layer: %d %s, mCurrent_Det+1: %d, mCurrent_Det: %d \n",KGRN,i_layer,KNRM,mCurrent_Det+1,mCurrent_Det);
		  if (!(prediction(mDist))){
			break;
		   }

		  mCurrent_Det++;
		  if (mMeasurements[i_layer] != 0) {
			if ((Int_t)(mTrack[i_layer]->get_TRD_det() / 30) != (Int_t)(mCurrent_Det / 30)) {
			}
			mCov_Res_Inv = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
			mCov_Res_Inv.Invert();
			
			mEstimate[i_layer] = mMu;
			//cout<<"mMu bef corri: "<<mMu<<endl;

			correction(mMeasurements[i_layer]);
			//cout<<"mMu after corri: "<<mMu<<endl;
		  }else{
			  mEstimate[i_layer] = mMu;
		  }
		  //cout << "Estimate: " << mEstimate[i_layer] << endl;
		}
		if (mShow) {
		  //	cout<<"Layer:"<< i_layer<<endl;
		  cout << "cov:" << mCov << endl;
		  cout << "Mu:" << mMu << endl;
		  cout << "unc" << mUnc << endl;
		  cout << "mes:" << mMeasurements[5] << endl;
		}
		
		//Loop again for better fit
		for (Int_t i_layer = mCurrent_Det % 6 - 1; i_layer >= 0; i_layer--) {
		  if (mShow){
			cout << "current 2. looplayer: " << i_layer << endl;
		}
		  mDist = mTRD_layer_radii_all[mCurrent_Det - 1] - mTRD_layer_radii_all[mCurrent_Det];

		  if (!(prediction(mDist)))
			break;

		  mCurrent_Det--;
		  if (mMeasurements[i_layer] != 0) {
			if ((Int_t)(mTrack[i_layer]->get_TRD_det() / 30) != (Int_t)(mCurrent_Det / 30)) {
			}
			mCov_Res_Inv = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
			//if(mShow) cout<<"mcov_res: "<<mCov_Res_Inv<<endl;
			mCov_Res_Inv.Invert();
			//if(mShow) cout<<"mcov_res_inv: "<<mCov_Res_Inv<<endl;
			correction(mMeasurements[i_layer]);
		  }
		}
	}
	//cout << mMu << endl;
    if (mPrimVertex == 1) {
      SVector<double, 5> mutemp = mMu;
      mDist = -mTRD_layer_radii_all[mCurrent_Det];
      //cout<<"mDist after pred: "<<mDist<<endl;
      //cout<<"mMu before pred: "<<mMu<<endl;

      Bool_t continueloop = prediction(mDist);
      if (continueloop) {
        //cout<<"mMu after pred: "<<mMu<<endl;
        SVector<double, 4> prim_vert_measurement;
        prim_vert_measurement[0] = 0;
        prim_vert_measurement[1] = 0;
        prim_vert_measurement[2] = mMu[2];
        prim_vert_measurement[3] = mMu[3];

        mCov_Res_Inv = mObs * mCov * Transpose(mObs) + 1.0 * mSig; //Measure uncertainty Matrix
        mCov_Res_Inv.Invert();

        correction(prim_vert_measurement);
        Double_t tempqpt = mMu[4];
        mMu = mutemp;
        mMu[4] = tempqpt;
      }
    }
    
    //Correct the tracklet offsets by calculating intersections with other layers where the pad tilt goes in the other direction and let Kalman filter run again.
    if(correct_offsets==1){
		
		//save old measurements for second correction round
		for(Int_t i_layer = 0; i_layer < 6; i_layer++){
			mMeasurements_old[i_layer] = mMeasurements[i_layer];
		}
		
		//cout << "----------------------------" << endl;
		//cout << "Correct offsets" << endl;
		vector<vector<vector<double>>> intersection_pts;	// intersection pts [y,z] of offset in one layer in all three other layers with opposite tilt
		intersection_pts.resize(6);
		
		//here use either all layers with opposite orientation or just adjacent layers
		//for(auto const& i_layer:{0,2,4}){ //over all even layer
		for(Int_t i_layer = 0; i_layer < 6; i_layer++){ //over adjacent layers
			if(mMeasurements[i_layer] == 0) continue;
			if(mEstimate[i_layer] == ROOT::Math::SVector<double, 5>()) continue;
			
			//for(auto const& j_layer:{1,3,5}){ //over all odd layers
			for(auto const& j_layer:{i_layer - 1, i_layer + 1}){
				if(j_layer < 0) continue; //only important if only iterating over adjacent layers
				if(j_layer > 5) continue; //--------------------------"--------------------------
				if(mMeasurements[j_layer] == 0) continue;
				if(mEstimate[j_layer] == ROOT::Math::SVector<double, 5>()) continue;
				
				//Get offset and z direction in layer i
				Int_t i_det = mTrack[i_layer]->get_TRD_det();
				ROOT::Math::SVector<double, 4> meas1 = mMeasurements[i_layer]; //measurement in layer i
				ROOT::Math::SVector<double, 5> mu1 = mEstimate[i_layer]; //estimate in layer i
				TVector3 offset1 = mTrack[i_layer]->get_TV3_offset();
				offset1.RotateZ(- (Double_t)(2 * (Int_t)(i_det / 30) + 1) * TMath::Pi() / 18.);
				Int_t i_pad = GetPadRowNumber(i_det, meas1[1]);
				if(i_pad==-1)continue;
				Double_t i_pad_centre = pad_geometry[i_det][i_pad][0];
				meas1[1] = i_pad_centre;
				
				//Get offset and z direction in layer j
				Int_t j_det = mTrack[j_layer]->get_TRD_det();
				ROOT::Math::SVector<double, 4> meas2 = mMeasurements[j_layer]; //estimate in layer j
				ROOT::Math::SVector<double, 5> mu2 = mEstimate[j_layer]; //estimate in layer j
				TVector3 offset2 = mTrack[j_layer]->get_TV3_offset();
				offset2.RotateZ(- (Double_t)(2 * (Int_t)(j_det / 30) + 1) * TMath::Pi() / 18.);
				Int_t j_pad = GetPadRowNumber(j_det, meas2[1]);
				if(j_pad==-1)continue;
				Double_t j_pad_centre = pad_geometry[j_det][j_pad][0];
				meas2[1] = j_pad_centre;
				
				//Transport offset j to layer i (or, alternatively, offset i to layer j)
				double dist = offset2[0] - offset1[0];
				//if(!transport(meas1, mu1, dist)) continue;
				//if(!transport(mu1, dist)) continue;
				if(!transport(meas2, mu2, -dist)) continue;
				if(!transport(mu2, -dist)) continue;

				//2d offset vector i in layer i (or, alternatively, layer j)
				TVector2 pt1{meas1[0], meas1[1]};

				//z direction vector i
				Double_t alpha1 = 2. * TMath::Power(-1, i_layer) * TMath::DegToRad();
				TVector3 tv3_dir1 = detector_geometry[3][i_det] * TMath::Cos(alpha1) - detector_geometry[2][i_det] * TMath::Sin(alpha1);
				tv3_dir1.RotateZ(- (Double_t)(2 * (Int_t)(i_det / 30) + 1) * TMath::Pi() / 18.);
				TVector2 tv2_dir1{tv3_dir1[1], tv3_dir1[2]};
				
				//2d offset vector j in layer i (or, alternatively, layer j)			
				TVector2 pt2{meas2[0], meas2[1]};
				
				//z direction vector j
				Double_t alpha2 = 2. * TMath::Power(-1, j_layer) * TMath::DegToRad();
				TVector3 tv3_dir2 = detector_geometry[3][j_det] * TMath::Cos(alpha2) - detector_geometry[2][j_det] * TMath::Sin(alpha2);
				tv3_dir2.RotateZ(- (Double_t)(2 * (Int_t)(j_det / 30) + 1) * TMath::Pi() / 18.);
				TVector2 tv2_dir2{tv3_dir2[1], tv3_dir2[2]};
				
				//Calculate intersection pts
				TVector2 pt_intersect = find_Intersect(pt1, tv2_dir1, pt2, tv2_dir2);
				
				//If intersection point is outside of pad
				//Take point at extreme end of pad
				//Or don't correct at all
				Double_t i_pad_up = pad_geometry[i_det][i_pad][0] + 0.5 * pad_geometry[i_det][i_pad][1];
				Double_t i_pad_low = pad_geometry[i_det][i_pad][0] - 0.5 * pad_geometry[i_det][i_pad][1];
							
				Double_t j_pad_up = pad_geometry[j_det][j_pad][0] + 0.5 * pad_geometry[i_det][i_pad][1];
				Double_t j_pad_low = pad_geometry[j_det][j_pad][0] - 0.5 * pad_geometry[j_det][j_pad][1];
				
				ROOT::Math::SVector<double, 4> meas_intersect = ROOT::Math::SVector<double, 4> (meas2);
				meas_intersect[0] = pt_intersect.X();
				meas_intersect[1] = pt_intersect.Y();
				transport(meas_intersect, mu2, dist);
				TVector2 pt_intersect2{meas_intersect[0], meas_intersect[1]};
				transport(meas2, mu2, dist);
				transport(mu2, dist);
				pt2 = {meas2[0], meas2[1]};
				
				//cout << "Intersection pts" << endl;
				//cout << "i" << endl;
				//cout << meas1[1] << " " << pt_intersect.Y() << " " << i_pad_low << " " << i_pad_up << endl;
				//cout << "j" << endl;
				//cout << meas2[1] << " " << pt_intersect2.Y() << " " << j_pad_low << " " << j_pad_up << endl;
				
				if((pt_intersect.Y() > i_pad_up) && (pt_intersect2.Y() < j_pad_low)){
					//Case 1: intersection point above pad 1, below pad 2
					//continue;
					//cout << "case 1" << endl;
					Double_t t1 = (i_pad_up - pt1.Y())/tv2_dir1.Y();
					Double_t t2 = (j_pad_low - pt2.Y())/tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = i_pad_up;
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = j_pad_low;
				}else if ((pt_intersect.Y() < i_pad_low) && (pt_intersect2.Y() > j_pad_up)){
					//Case 2: intersection point below pad 1, above pad 2
					//continue;
					//cout << "case 2" << endl;
					Double_t t1 = (i_pad_low - pt1.Y())/tv2_dir1.Y();
					Double_t t2 = (j_pad_up - pt2.Y())/tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = i_pad_low;
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = j_pad_up;
				}else if((pt_intersect.Y() > i_pad_up) || (pt_intersect2.Y() > j_pad_up)){
					//Case 3: intersection point above either pad 1 or 2
					//Use the lower value of the upper ends of the pads
					//continue;	
					//cout << "case 3" << endl;			
					Double_t t1 = pt_intersect.Y() - i_pad_up;
					Double_t t2 = pt_intersect2.Y() - j_pad_up;
					Double_t t = TMath::Max(t1, t2);
					t1 = t / tv2_dir1.Y();
					t2 = t / tv2_dir2.Y();
					
					meas1[0] = (pt_intersect - t1 * tv2_dir1).X();
					meas1[1] = (pt_intersect - t1 * tv2_dir1).Y();
					meas2[0] = (pt_intersect2 - t2 * tv2_dir2).X();
					meas2[1] = (pt_intersect2 - t2 * tv2_dir2).Y();
				}else if((pt_intersect.Y() < i_pad_low) || (pt_intersect2.Y() < j_pad_low)){
					//Case 4: intersection point below either pad 1 or pad 2
					//Use the higher value of the lower ends of the pads
					//continue;
					//cout << "case 4" << endl;
					Double_t t1 = pt_intersect.Y() - i_pad_low;
					Double_t t2 = pt_intersect2.Y() - j_pad_low;
					Double_t t = TMath::Min(t1, t2);
					t1 = t / tv2_dir1.Y();
					t2 = t / tv2_dir2.Y();
					
					meas1[0] = (pt_intersect - t1 * tv2_dir1).X();
					meas1[1] = (pt_intersect - t1 * tv2_dir1).Y();
					meas2[0] = (pt_intersect2 - t2 * tv2_dir2).X();
					meas2[1] = (pt_intersect2 - t2 * tv2_dir2).Y();
				}else{
					//Case 5: intersection point within pads
					//cout << "case 5" << endl;
					meas1[0] = pt_intersect.X();
					meas1[1] = pt_intersect.Y();
					meas2[0] = pt_intersect2.X();
					meas2[1] = pt_intersect2.Y();
				}
				//cout << "" << endl;
				//cout << "i" << endl;
				//cout << meas1[1] << " " << pt_intersect.Y() << " " << i_pad_low << " " << i_pad_up << endl;
				//cout << "j" << endl;
				//cout << meas2[1] << " " << pt_intersect2.Y() << " " << j_pad_low << " " << j_pad_up << endl;
				//cout << "\n" << endl;
				 
				/*meas1[0] = pt_intersect.X();
				meas1[1] = pt_intersect.Y();
				meas2[0] = pt_intersect.X();
				meas2[1] = pt_intersect.Y();*/			 
				 
				/*Double_t delta_z1 = pt_intersect.Y() - pt1.Y();
				Double_t delta_z2 = pt_intersect.Y() - pt2.Y();
				 
				Double_t i_pad_size = pad_geometry[i_det][i_pad][1];
				Double_t j_pad_size = pad_geometry[j_det][j_pad][1];
				
				if((delta_z1 > i_pad_size/2.) && (delta_z2 < -j_pad_size/2.)){
					//Case 1: intersection point above pad 1, below pad 2
					//continue; //disregard point
					
					Double_t t1 = i_pad_size / 2. / tv2_dir1.Y();
					Double_t t2 = -j_pad_size / 2. / tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else if ((delta_z1 < -i_pad_size/2.) && (delta_z2 > j_pad_size/2.)){
					//Case 2: intersection point below pad 1, above pad 2
					//continue; //disregard point
					
					Double_t t1 = -i_pad_size / 2. / tv2_dir1.Y();
					Double_t t2 = j_pad_size / 2. / tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else if((delta_z1 > i_pad_size/2.) || (delta_z2 > j_pad_size/2.)){
					//continue;
					//Case 3: intersection point above either pad 1 or 2
					//Use the lower value of the upper ends of the pads
					continue; //disregard point
					
					Double_t t1 = TMath::Min(i_pad_size / 2. / tv2_dir1.Y(), (i_pad_size / 2. - (pt1.Y() - pt2.Y()))/tv2_dir1.Y());
					Double_t t2 = TMath::Min(j_pad_size / 2. / tv2_dir2.Y(), (j_pad_size / 2. - (pt2.Y() - pt1.Y()))/tv2_dir2.Y());
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else if((delta_z1 < -i_pad_size/2.) || (delta_z2 < -j_pad_size/2.)){
					//continue;
					//Case 4: intersection point below either pad 1 or pad 2
					//Use the higher value of the lower ends of the pads
					continue; //disregard point
					
					Double_t t1 = TMath::Max(-i_pad_size / 2. / tv2_dir1.Y(), (-i_pad_size/2. + (pt1.Y() - pt2.Y()))/tv2_dir1.Y());
					Double_t t2 = TMath::Max(-j_pad_size / 2. / tv2_dir2.Y(), (-j_pad_size/2. + (pt2.Y() - pt1.Y()))/tv2_dir2.Y());
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else{
					//Case 5: intersection point within pads
					meas1[0] = pt_intersect.X();
					meas1[1] = pt_intersect.Y();
					meas2[0] = pt_intersect.X();
					meas2[1] = pt_intersect.Y();
				}
				
				//Transport intersection point back to layer j (or, alternatively, layer i)
				//transport(meas1, mu1, -dist); 
				//transport(mu1, -dist);
				transport(meas2, mu2, dist);
				transport(mu2, dist);*/
				
				intersection_pts[i_layer].push_back({meas1[0], meas1[1]});
				//intersection_pts[j_layer].push_back({meas2[0], meas2[1]}); //Uncomment this if iteration is done over odd and even layers instead of adjacent layers
			}
		}
		
		for(int i_layer=0; i_layer<6; i_layer++){ 
			//Take mean of intersection pts in every layer
			if(!mTrack[i_layer]) continue;
			if(intersection_pts[i_layer].size()==0){
				Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*mTrack[i_layer]);
				mFound_tracks_corr.back()[i_layer]=tracklet;
				//delete tracklet;
				continue;
			}
			mMeasurements[i_layer][0] = 0;
			mMeasurements[i_layer][1] = 0;
			for(unsigned int j=0; j<intersection_pts[i_layer].size();j++){
				mMeasurements[i_layer][0] += intersection_pts[i_layer][j][0];
				mMeasurements[i_layer][1] += intersection_pts[i_layer][j][1];
			}
			mMeasurements[i_layer][0] /= intersection_pts[i_layer].size();
			mMeasurements[i_layer][1] /= intersection_pts[i_layer].size();
			//cout << "New measurement mean " << i_layer << ": " << mMeasurements[i_layer]<< endl;
			
			//Change the offset of the track to the corrected one
			Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*mTrack[i_layer]);
			TVector3 offset = tracklet->get_TV3_offset();
			Int_t det = mFound_tracks.back()[i_layer]->get_TRD_det();
			offset.RotateZ(- (Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18.);
			offset[1] = mMeasurements[i_layer][0];
			offset[2] = mMeasurements[i_layer][1];
			offset.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18.);
			tracklet->set_TV3_offset(offset);
			mFound_tracks_corr.back()[i_layer]=tracklet;
			//delete tracklet;	
		}
		
		//Run Kalman filter over corrected measurements
		//cout << "---------------------------" << endl;
		//cout << "Third Kalman filter" << endl;
		for(int i_layer = 0; i_layer<6; i_layer++){
			if(mMeasurements[i_layer]==0){
				continue;
			}else{
				for (int i_param = 0; i_param < 4; i_param++) mMu[i_param] = mMeasurements[i_layer][i_param];
				mCurrent_Det = mTrack[i_layer]->get_TRD_det();
				break;
			}
		}
		
		int first_det = mCurrent_Det;
		mEstimate[mCurrent_Det % 6] = mMu;
		mMu[4] = 0.0 / 1.0;

		// 0,0 = y, 1,1=z , 2,2=sin phi , 3,3=tan lambda

		Double_t dy2 = 0.2;         // 0.2  0.4
		Double_t dz2 = 4.0;         // 4.0  4.0
		Double_t dsin_phi = 7.0;    // 7.0  10.0
		Double_t dsin_theta = 18.0; // 18.0  25.0
		Double_t dpT = 10.0;        // 10.0  10.0

		mSig[0][0] = dy2;                                                           // 0.2
		mSig[1][1] = dz2;                                                           // 4.0
		mSig[2][2] = TMath::Power(TMath::Sin(dsin_phi * TMath::Pi() / 180.0), 2);   // 7.0
		mSig[3][3] = TMath::Power(TMath::Tan(dsin_theta * TMath::Pi() / 180.0), 2); // 18.0

		mCov = SMatrixIdentity();
		mCov[0][0] = dy2;                                                           // 0.2
		mCov[1][1] = dz2;                                                           // 4.0
		mCov[2][2] = TMath::Power(TMath::Sin(dsin_phi * TMath::Pi() / 180.0), 2);   // 7.0
		mCov[3][3] = TMath::Power(TMath::Tan(dsin_theta * TMath::Pi() / 180.0), 2); // 20.0
																					//mCov[4][4]	=	0.09; // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15
		mCov[4][4] = dpT * dpT;                                                     // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15

		mCovs[mCurrent_Det % 6] = mCov;

		mCurrent_Det++;
		//cout << mMu << endl;
		for (Int_t i_layer = mCurrent_Det % 6 ; i_layer < 6; i_layer++) {
			//cout << "\nLayer: " << i_layer << endl;

			mDist = mTRD_layer_radii_all[mCurrent_Det] - mTRD_layer_radii_all[mCurrent_Det - 1];
			if(!(prediction(mDist))) break;
			//cout << "after prediction: " << (i_layer - 1) << " (" << mCurrent_Det%6 << ")" << endl;
			mCovsTrans[i_layer - 1] = mCov;
			mTrans[i_layer - 1] = A;
			mTrans_t[i_layer - 1] = A_t; 
			
			//cout << "Prediction: " << mMu << endl;
			
			mCurrent_Det++;
			
			if (mMeasurements[i_layer] != 0) {
				if ((Int_t)(mTrack[i_layer]->get_TRD_det() / 30) != (Int_t)(mCurrent_Det / 30)) {
				}
				mCov_Res_Inv = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
				mCov_Res_Inv.Invert();
				//cout<<"mMu bef corri: "<<mMu<<endl;

				correction(mMeasurements[i_layer]);
				//cout << "Measurement: " << mMeasurements[i_layer] << endl;
				//cout<<"mMu after corri: "<<mMu<<endl;
			}
			//cout << "after correction: " << i_layer << endl;
			mCovs [i_layer] = mCov; //Check if here or in if 
			mEstimate[i_layer] = mMu;
			//cout << "Estimate: " << mEstimate[i_layer] << endl;
		}

		mCurrent_Det --;
		
		if(use_RTS){
			//Second RTS Smoother
			//cout << mCurrent_Det << endl;
			for (Int_t i_layer = mCurrent_Det % 6 - 1; i_layer >= 0; i_layer--){
				//cout << "RTS smoother" << endl;
				ROOT::Math::SMatrix<double, 5, 5> mCovsTrans_inv = mCovsTrans[i_layer];
				mCovsTrans_inv.Invert();
				mKal_Gain_RTS = mCovs[i_layer]*(mTrans_t[i_layer]*mCovsTrans_inv);
				double dist = mTRD_layer_radii[i_layer+1][1]-mTRD_layer_radii[i_layer][1];
				ROOT::Math::SMatrix<double, 5, 5> mTrans_inv = mTrans[i_layer];
				mTrans_inv.Invert();
				ROOT::Math::SVector<double, 5> mEstimateTrans = mEstimate[i_layer];
				transport(mEstimateTrans, dist);
				//cout << "i_layer" << i_layer << endl;
				//cout << "old " <<  mEstimate[i_layer] << endl;
				//cout << "transported: " << mTrans[i_layer]*mEstimate[i_layer] << endl;
				mEstimate[i_layer] += mKal_Gain_RTS*(mEstimate[i_layer+1] - mEstimateTrans);
				//mEstimate[i_layer] = mEstimate[i_layer] + mKal_Gain_RTS*(mEstimate[i_layer+1] - mTrans[i_layer]*mEstimate[i_layer]);
				//cout << "newEstimate" << endl;
				//cout << "new " << mEstimate[i_layer] << endl;
				mCovs[i_layer] = mCovs[i_layer] + mKal_Gain_RTS*((mCovs[i_layer+1] - mCovsTrans[i_layer])*Transpose(mKal_Gain_RTS));
			}
	
		}else{
		
			//cout << "----------------------" << endl;
			//cout << "Fourth Kalman filter" << endl;
			for (Int_t i_layer = mCurrent_Det % 6 - 1; i_layer >= 0; i_layer--) {
			  
			  //cout << "\nLayer" << i_layer << endl;
			  
			  if (mShow)
				cout << "current 2. looplayer: " << i_layer << endl;
			  mDist = mTRD_layer_radii_all[mCurrent_Det - 1] - mTRD_layer_radii_all[mCurrent_Det];

			  if (!(prediction(mDist)))
				break;
			  //cout << "Prediction: " << mMu << endl;

			  mCurrent_Det--;
			  //if (mMeasurements[i_layer] != 0) {
			  //  if ((Int_t)(mTrack[i_layer]->get_TRD_det() / 30) != (Int_t)(mCurrent_Det / 30)) {
			  //  }
			  //  mCov_Res_Inv = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
			  //  mCov_Res_Inv.Invert();
			  //}
			  mEstimate[i_layer] = mMu;
			  //cout << i_layer << ": " << mMu << endl;
			}
			
			//cout << "Estimate first" << endl;
			//cout << mEstimates.back()[0] << " " << mEstimates.back()[1] << " " << mEstimates.back()[2] << " " << mEstimates.back()[3] << " " << mEstimates.back()[4] << " " << endl;
			mEstimates.back() = mEstimate;
			//cout << "Estimate then" << endl;
			//cout << mEstimates.back()[0] << " " << mEstimates.back()[1] << " " << mEstimates.back()[2] << " " << mEstimates.back()[3] << " " << mEstimates.back()[4] << " " << endl;
			
			//cout << mMu << endl;
		}
		
	//}

    //Do the same as before but use the corrected Kalman fit
	
		for(Int_t i_layer = 0; i_layer < 6; i_layer++){
			intersection_pts[i_layer] =  vector<vector<double>>();
		}

		//for(auto const& i_layer:{0,2,4}){ //over all even layer
		for(Int_t i_layer = 0; i_layer < 6; i_layer++){
			if(mMeasurements_old[i_layer] == 0) continue;
			if(mEstimate[i_layer] == ROOT::Math::SVector<double, 5>()) continue;
			
			//for(auto const& j_layer:{1,3,5}){ //over all odd layers
			for(auto const& j_layer:{i_layer - 1, i_layer + 1}){
				if(j_layer < 0) continue;
				if(j_layer > 5) continue;
				if(mMeasurements_old[j_layer] == 0) continue;
				//if(mEstimate[j_layer] == ROOT::Math::SVector<double, 5>()) continue;
				
				//Get offset and z direction in layer i
				Int_t i_det = mTrack[i_layer]->get_TRD_det();
				ROOT::Math::SVector<double, 4> meas1 = mMeasurements_old[i_layer]; //measurement in layer i
				ROOT::Math::SVector<double, 5> mu1 = mEstimate[i_layer]; //estimate in layer i
				TVector3 offset1 = mTrack[i_layer]->get_TV3_offset();
				offset1.RotateZ(- (Double_t)(2 * (Int_t)(i_det / 30) + 1) * TMath::Pi() / 18.);
				Int_t i_pad = GetPadRowNumber(i_det, meas1[1]);
				if(i_pad==-1)continue;
				Double_t i_pad_centre = pad_geometry[i_det][i_pad][0];
				meas1[1] = i_pad_centre;
				
				//Get offset and z direction in layer j
				Int_t j_det = mTrack[j_layer]->get_TRD_det();
				ROOT::Math::SVector<double, 4> meas2 = mMeasurements_old[j_layer]; //estimate in layer j
				ROOT::Math::SVector<double, 5> mu2 = mEstimate[j_layer]; //estimate in layer j
				TVector3 offset2 = mTrack[j_layer]->get_TV3_offset();
				offset2.RotateZ(- (Double_t)(2 * (Int_t)(j_det / 30) + 1) * TMath::Pi() / 18.);
				Int_t j_pad = GetPadRowNumber(j_det, meas2[1]);
				if(j_pad==-1)continue;
				Double_t j_pad_centre = pad_geometry[j_det][j_pad][0];
				meas2[1] = j_pad_centre;
				
				double dist = offset2[0] - offset1[0];
				//if(!transport(meas1, mu1, dist)) continue;
				//if(!transport(mu1, dist)) continue;
				if(!transport(meas2, mu2, -dist)) continue;
				if(!transport(mu2, -dist)) continue;
				
				//ROOT::Math::SVector<double, 5> mu1_var{meas1[0], meas1[1], meas1[2], meas1[3], mMu[4]};
				//if(!transport(meas1, mu1_var, dist)) continue;
				

				//Offset vector 1
				TVector2 pt1{meas1[0], meas1[1]};

				//Direction vector 1
				Double_t alpha1 = 2. * TMath::Power(-1, i_layer) * TMath::DegToRad();
				TVector3 tv3_dir1 = detector_geometry[3][i_det] * TMath::Cos(alpha1) - detector_geometry[2][i_det] * TMath::Sin(alpha1);
				tv3_dir1.RotateZ(- (Double_t)(2 * (Int_t)(i_det / 30) + 1) * TMath::Pi() / 18.);
				TVector2 tv2_dir1{tv3_dir1[1], tv3_dir1[2]};
				
				//Offset vector 2
				TVector2 pt2{meas2[0], meas2[1]};
				
				//Direction vector 2
				Double_t alpha2 = 2. * TMath::Power(-1, j_layer) * TMath::DegToRad();
				TVector3 tv3_dir2 = detector_geometry[3][j_det] * TMath::Cos(alpha2) - detector_geometry[2][j_det] * TMath::Sin(alpha2);
				tv3_dir2.RotateZ(- (Double_t)(2 * (Int_t)(j_det / 30) + 1) * TMath::Pi() / 18.);
				TVector2 tv2_dir2{tv3_dir2[1], tv3_dir2[2]};
				
				//Calculate intersection pts
				TVector2 pt_intersect = find_Intersect(pt1, tv2_dir1, pt2, tv2_dir2);
				
				//If intersection point is outside of pad
				//Take point at extreme end of pad
				//Or don't correct at all
				Double_t i_pad_up = pad_geometry[i_det][i_pad][0] + 0.5 * pad_geometry[i_det][i_pad][1];
				Double_t i_pad_low = pad_geometry[i_det][i_pad][0] - 0.5 * pad_geometry[i_det][i_pad][1];
						
				Double_t j_pad_up = pad_geometry[j_det][j_pad][0] + 0.5 * pad_geometry[i_det][i_pad][1];
				Double_t j_pad_low = pad_geometry[j_det][j_pad][0] - 0.5 * pad_geometry[j_det][j_pad][1];
				
				ROOT::Math::SVector<double, 4> meas_intersect = ROOT::Math::SVector<double, 4> (meas2);
				meas_intersect[0] = pt_intersect.X();
				meas_intersect[1] = pt_intersect.Y();
				transport(meas_intersect, mu2, dist);
				TVector2 pt_intersect2{meas_intersect[0], meas_intersect[1]};
				transport(meas2, mu2, dist);
				transport(mu2, dist);
				pt2 = {meas2[0], meas2[1]};
				if((pt_intersect.Y() > i_pad_up) && (pt_intersect2.Y() < j_pad_low)){
					//Case 1: intersection point above pad 1, below pad 2
					//continue;
					Double_t t1 = (i_pad_up - pt1.Y())/tv2_dir1.Y();
					Double_t t2 = (j_pad_low - pt2.Y())/tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = i_pad_up;
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = j_pad_low;
				}else if ((pt_intersect.Y() < i_pad_low) && (pt_intersect2.Y() > j_pad_up)){
					//Case 2: intersection point below pad 1, above pad 2
					//continue;
					Double_t t1 = (i_pad_low - pt1.Y())/tv2_dir1.Y();
					Double_t t2 = (j_pad_up - pt2.Y())/tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = i_pad_low;
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = j_pad_up;
				}else if((pt_intersect.Y() > i_pad_up) || (pt_intersect2.Y() > j_pad_up)){
					//Case 3: intersection point above either pad 1 or 2
					//Use the lower value of the upper ends of the pads
					//continue;				
					Double_t t1 = pt_intersect.Y() - i_pad_up;
					Double_t t2 = pt_intersect2.Y() - j_pad_up;
					Double_t t = TMath::Max(t1, t2);
					t1 = t / tv2_dir1.Y();
					t2 = t / tv2_dir2.Y();
					
					meas1[0] = (pt_intersect - t1 * tv2_dir1).X();
					meas1[1] = (pt_intersect - t1 * tv2_dir1).Y();
					meas2[0] = (pt_intersect2 - t2 * tv2_dir2).X();
					meas2[1] = (pt_intersect2 - t2 * tv2_dir2).Y();
				}else if((pt_intersect.Y() < i_pad_low) || (pt_intersect2.Y() < j_pad_low)){
					//Case 4: intersection point below either pad 1 or pad 2
					//Use the higher value of the lower ends of the pads
					//continue;
					Double_t t1 = pt_intersect.Y() - i_pad_low;
					Double_t t2 = pt_intersect2.Y() - j_pad_low;
					Double_t t = TMath::Min(t1, t2);
					t1 = t / tv2_dir1.Y();
					t2 = t / tv2_dir2.Y();
					
					meas1[0] = (pt_intersect - t1 * tv2_dir1).X();
					meas1[1] = (pt_intersect - t1 * tv2_dir1).Y();
					meas2[0] = (pt_intersect2 - t2 * tv2_dir2).X();
					meas2[1] = (pt_intersect2 - t2 * tv2_dir2).Y();
				}else{
					//Case 5: intersection point within pads
					meas1[0] = pt_intersect.X();
					meas1[1] = pt_intersect.Y();
					meas2[0] = pt_intersect2.X();
					meas2[1] = pt_intersect2.Y();
				}
				
				/*meas1[0] = pt_intersect.X();
				meas1[1] = pt_intersect.Y();
				meas2[0] = pt_intersect.X();
				meas2[1] = pt_intersect.Y();*/
				
				/*Double_t delta_z1 = pt_intersect.Y() - pt1.Y();
				Double_t delta_z2 = pt_intersect.Y() - pt2.Y();
				 
				Double_t i_pad_size = pad_geometry[i_det][i_pad][1];
				Double_t j_pad_size = pad_geometry[j_det][j_pad][1];
				
				if((delta_z1 > i_pad_size/2.) && (delta_z2 < -j_pad_size/2.)){
					//Case 1: intersection point above pad 1, below pad 2
					//continue; //disregard point
					
					Double_t t1 = i_pad_size / 2. / tv2_dir1.Y();
					Double_t t2 = -j_pad_size / 2. / tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else if ((delta_z1 < -i_pad_size/2.) && (delta_z2 > j_pad_size/2.)){
					//Case 2: intersection point below pad 1, above pad 2
					//continue; //disregard point
					
					Double_t t1 = -i_pad_size / 2. / tv2_dir1.Y();
					Double_t t2 = j_pad_size / 2. / tv2_dir2.Y();
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else if((delta_z1 > i_pad_size/2.) || (delta_z2 > j_pad_size/2.)){
					//continue;
					//Case 3: intersection point above either pad 1 or 2
					//Use the lower value of the upper ends of the pads
					continue; //disregard point
					
					Double_t t1 = TMath::Min(i_pad_size / 2. / tv2_dir1.Y(), (i_pad_size / 2. - (pt1.Y() - pt2.Y()))/tv2_dir1.Y());
					Double_t t2 = TMath::Min(j_pad_size / 2. / tv2_dir2.Y(), (j_pad_size / 2. - (pt2.Y() - pt1.Y()))/tv2_dir2.Y());
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else if((delta_z1 < -i_pad_size/2.) || (delta_z2 < -j_pad_size/2.)){
					//continue;
					//Case 4: intersection point below either pad 1 or pad 2
					//Use the higher value of the lower ends of the pads
					continue; //disregard point
					
					Double_t t1 = TMath::Max(-i_pad_size / 2. / tv2_dir1.Y(), (-i_pad_size/2. + (pt1.Y() - pt2.Y()))/tv2_dir1.Y());
					Double_t t2 = TMath::Max(-j_pad_size / 2. / tv2_dir2.Y(), (-j_pad_size/2. + (pt2.Y() - pt1.Y()))/tv2_dir2.Y());
					
					meas1[0] = (pt1 + t1 * tv2_dir1).X();
					meas1[1] = (pt1 + t1 * tv2_dir1).Y();
					meas2[0] = (pt2 + t2 * tv2_dir2).X();
					meas2[1] = (pt2 + t2 * tv2_dir2).Y();
				}else{
					//Case 5: intersection point within pads
					meas1[0] = pt_intersect.X();
					meas1[1] = pt_intersect.Y();
					meas2[0] = pt_intersect.X();
					meas2[1] = pt_intersect.Y();
				}
				
				//Transport intersection point back to layer j (or, alternatively, layer i)
				//transport(meas1, mu1, -dist); 
				//transport(mu1, -dist);
				transport(meas2, mu2, dist);
				transport(mu2, dist);*/
				
				intersection_pts[i_layer].push_back({meas1[0], meas1[1]});
				//intersection_pts[j_layer].push_back({meas2[0], meas2[1]});
			}
		}

		
		for(int i_layer=0; i_layer<6; i_layer++){ //Take mean of intersection pts in every layer
			if(!mTrack[i_layer]) continue;
			if(intersection_pts[i_layer].size()==0){
				Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*mTrack[i_layer]);
				mFound_tracks_corr.back()[i_layer]=tracklet;
				continue;
			}
			mMeasurements[i_layer][0] = 0;
			mMeasurements[i_layer][1] = 0;
			for(unsigned int j=0; j<intersection_pts[i_layer].size();j++){
				mMeasurements[i_layer][0] += intersection_pts[i_layer][j][0];
				mMeasurements[i_layer][1] += intersection_pts[i_layer][j][1];
			}
			mMeasurements[i_layer][0] /= intersection_pts[i_layer].size();
			mMeasurements[i_layer][1] /= intersection_pts[i_layer].size();
			//cout << "New measurement mean " << i_layer << ": " << mMeasurements[i_layer]<< endl;
			
			
			Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*mTrack[i_layer]);
			TVector3 offset = tracklet->get_TV3_offset();
			Int_t det = mFound_tracks.back()[i_layer]->get_TRD_det();
			offset.RotateZ(- (Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18.);
			offset[1] = mMeasurements[i_layer][0];
			offset[2] = mMeasurements[i_layer][1];
			offset.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18.);
			tracklet->set_TV3_offset(offset);
			//cout << "i_layer: " << offset[0] << " " << offset[1] << " " << offset[2] << endl;
			mFound_tracks_corr.back()[i_layer]=tracklet;
				
		}	
		
		//Run Kalman filter over corrected measurements
		//cout << "---------------------------" << endl;
		//cout << "Third Kalman filter" << endl;
		for(int i_layer = 0; i_layer<6; i_layer++){
			if(mMeasurements[i_layer]==0){
				continue;
			}else{
				for (int i_param = 0; i_param < 4; i_param++) mMu[i_param] = mMeasurements[i_layer][i_param];
				mCurrent_Det = mTrack[i_layer]->get_TRD_det();
				break;
			}
		}
		
		first_det = mCurrent_Det;
		mEstimate[mCurrent_Det % 6] = mMu;
		mMu[4] = 0.0 / 1.0;

		// 0,0 = y, 1,1=z , 2,2=sin phi , 3,3=tan lambda

		dy2 = 0.2;         // 0.2  0.4
		dz2 = 4.0;         // 4.0  4.0
		dsin_phi = 7.0;    // 7.0  10.0
		dsin_theta = 18.0; // 18.0  25.0
		dpT = 10.0;        // 10.0  10.0

		mSig[0][0] = dy2;                                                           // 0.2
		mSig[1][1] = dz2;                                                           // 4.0
		mSig[2][2] = TMath::Power(TMath::Sin(dsin_phi * TMath::Pi() / 180.0), 2);   // 7.0
		mSig[3][3] = TMath::Power(TMath::Tan(dsin_theta * TMath::Pi() / 180.0), 2); // 18.0

		mCov = SMatrixIdentity();
		mCov[0][0] = dy2;                                                           // 0.2
		mCov[1][1] = dz2;                                                           // 4.0
		mCov[2][2] = TMath::Power(TMath::Sin(dsin_phi * TMath::Pi() / 180.0), 2);   // 7.0
		mCov[3][3] = TMath::Power(TMath::Tan(dsin_theta * TMath::Pi() / 180.0), 2); // 20.0
																					//mCov[4][4]	=	0.09; // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15
		mCov[4][4] = dpT * dpT;                                                     // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15

		mCovs[mCurrent_Det % 6] = mCov;

		mCurrent_Det++;
		//cout << mMu << endl;
		for (Int_t i_layer = mCurrent_Det % 6 ; i_layer < 6; i_layer++) {
			//cout << "\nLayer: " << i_layer << endl;

			mDist = mTRD_layer_radii_all[mCurrent_Det] - mTRD_layer_radii_all[mCurrent_Det - 1];
			if(!(prediction(mDist))) break;
			//cout << "after prediction: " << (i_layer - 1) << " (" << mCurrent_Det%6 << ")" << endl;
			mCovsTrans[i_layer - 1] = mCov;
			mTrans[i_layer - 1] = A;
			mTrans_t[i_layer - 1] = A_t; 
			
			//cout << "Prediction: " << mMu << endl;
			
			mCurrent_Det++;
			
			if (mMeasurements[i_layer] != 0) {
				if ((Int_t)(mTrack[i_layer]->get_TRD_det() / 30) != (Int_t)(mCurrent_Det / 30)) {
				}
				mCov_Res_Inv = mObs * mCov * Transpose(mObs) + mSig; //Measure uncertainty Matrix
				mCov_Res_Inv.Invert();
				//cout<<"mMu bef corri: "<<mMu<<endl;

				correction(mMeasurements[i_layer]);
				//cout << "Measurement: " << mMeasurements[i_layer] << endl;
				//cout<<"mMu after corri: "<<mMu<<endl;
			}
			//cout << "after correction: " << i_layer << endl;
			mCovs [i_layer] = mCov; //Check if here or in if 
			mEstimate[i_layer] = mMu;
			//cout << "Estimate: " << mEstimate[i_layer] << endl;
		}

		mCurrent_Det --;
		
		if(use_RTS){
			//Second RTS Smoother
			//cout << mCurrent_Det << endl;
			for (Int_t i_layer = mCurrent_Det % 6 - 1; i_layer >= 0; i_layer--){
				//cout << "RTS smoother" << endl;
				ROOT::Math::SMatrix<double, 5, 5> mCovsTrans_inv = mCovsTrans[i_layer];
				mCovsTrans_inv.Invert();
				mKal_Gain_RTS = mCovs[i_layer]*(mTrans_t[i_layer]*mCovsTrans_inv);
				double dist = mTRD_layer_radii[i_layer+1][1]-mTRD_layer_radii[i_layer][1];
				ROOT::Math::SMatrix<double, 5, 5> mTrans_inv = mTrans[i_layer];
				mTrans_inv.Invert();
				ROOT::Math::SVector<double, 5> mEstimateTrans = mEstimate[i_layer];
				transport(mEstimateTrans, dist);
				//cout << "i_layer" << i_layer << endl;
				//cout << "old " <<  mEstimate[i_layer] << endl;
				//cout << "transported: " << mTrans[i_layer]*mEstimate[i_layer] << endl;
				mEstimate[i_layer] += mKal_Gain_RTS*(mEstimate[i_layer+1] - mEstimateTrans);
				//mEstimate[i_layer] = mEstimate[i_layer] + mKal_Gain_RTS*(mEstimate[i_layer+1] - mTrans[i_layer]*mEstimate[i_layer]);
				//cout << "newEstimate" << endl;
				//cout << "new " << mEstimate[i_layer] << endl;
				mCovs[i_layer] = mCovs[i_layer] + mKal_Gain_RTS*((mCovs[i_layer+1] - mCovsTrans[i_layer])*Transpose(mKal_Gain_RTS));
			}
	
		}else{
		
			//cout << "----------------------" << endl;
			//cout << "Fourth Kalman filter" << endl;
			for (Int_t i_layer = mCurrent_Det % 6 - 1; i_layer >= 0; i_layer--) {
			  
			  //cout << "\nLayer" << i_layer << endl;
			  
			  if (mShow)
				cout << "current 2. looplayer: " << i_layer << endl;
			  mDist = mTRD_layer_radii_all[mCurrent_Det - 1] - mTRD_layer_radii_all[mCurrent_Det];

			  if (!(prediction(mDist)))
				break;

			  mCurrent_Det--;
			  mEstimate[i_layer] = mMu;
			}
			
			mEstimates.back() = mEstimate;

		}
		
	}


    //calculate Helix param
    //mMu = mEstimate[0];
    //if(mTrack[0]){
	//	mCurrent_Det = mTrack[0]->get_TRD_det();
	//}else{
	//	mCurrent_Det = 1;}
    
    Double_t charge = 1.;
    if (mMu[4] < 0)
      charge = -1.;
    Double_t lam = TMath::ATan(mMu[3]);
    Double_t pxy = charge / mMu[4];

    TVector3 x_vek;
    TVector3 p_vek;
    x_vek[0] = mTRD_layer_radii_all[mCurrent_Det];
    /*if(mPrimVertex==1)
			x_vek[0]=0;*/
    x_vek[1] = mMu[0];
    x_vek[2] = mMu[1];
    p_vek[0] = TMath::Cos(TMath::ASin(mMu[2])) * pxy;
    p_vek[1] = mMu[2] * pxy;
    p_vek[2] = mMu[3] * pxy;
    x_vek.RotateZ((Double_t)(2 * (Int_t)(mCurrent_Det / 30) + 1) * TMath::Pi() / 18);
    p_vek.RotateZ((Double_t)(2 * (Int_t)(mCurrent_Det / 30) + 1) * TMath::Pi() / 18);
    Double_t x[3];
    Double_t p[3];
    x[0] = x_vek[0];
    x[1] = x_vek[1];
    x[2] = x_vek[2];
    p[0] = p_vek[0];
    p[1] = p_vek[1];
    p[2] = p_vek[2];
    //calculation of Helixparameter taken from http://alidoc.cern.ch/AliRoot/v5-09-36/_ali_helix_8cxx_source.html
    //AliHelix::AliHelix(Double_t x[3], Double_t p[3], Double_t charge, Double_t conversion)
    vector<Double_t> fHelix;
    fHelix.resize(8);
    Double_t pt = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]);
    //

    Double_t b_fak = b_field * 3. / 1000.;

    Double_t curvature = (mMu[4] * b_fak);

    fHelix[4] = curvature; // C
    fHelix[3] = p[2] / pt; // tgl
                           //
    Double_t xc, yc, rc;
    rc = 1 / fHelix[4];
    xc = x[0] - rc * p[1] / pt;
    yc = x[1] + rc * p[0] / pt;
    //
    fHelix[5] = x[0]; // x0
    fHelix[0] = x[1]; // y0
    fHelix[1] = x[2]; // z0
      //
    //fHelix[6] = xc;
    //fHelix[7] = yc;
    //fHelix[8] = TMath::Abs(rc);
    //
    fHelix[5] = xc;
    fHelix[0] = yc;

    fHelix[6] = pt;
    fHelix[7] = p[2];
    //
    if (TMath::Abs(p[1]) < TMath::Abs(p[0])) {
      fHelix[2] = TMath::ASin(p[1] / pt);
      //Helix[2]=asinf(p[1]/pt);
      if (charge * yc < charge * x[1])
        fHelix[2] = TMath::Pi() - fHelix[2];
    } else {
      fHelix[2] = TMath::ACos(p[0] / pt);
      //fHelix[2]=acosf(p[0]/pt);
      if (charge * xc > charge * x[0])
        fHelix[2] = -fHelix[2];
    }

	//Double_t delta_diff[6] = {0.,0.,0.,0.,0.,0.};
	//Double_t dist_diff[6] = {0.,0.,0.,0.,0.,0.};

	if(correct_offsets==2){
		det_geom = detector_geometry;
		track = mTrack;
				
		//Start fitter
		TVirtualFitter *min = TVirtualFitter::Fitter(0,2);
        min->SetFCN(Chi2_Kalman_tracklet);
        Double_t pStart[6] = {fHelix[0], fHelix[1], fHelix[2], fHelix[3], fHelix[4], fHelix[5]}; //Use Kalman helix as initial parameters

		//cout << "pars: " << pStart << endl;

        min->SetParameter(0,"xc",pStart[0],0.01,0,0);
        min->SetParameter(1,"yc",pStart[1],0.01,0,0);
        min->SetParameter(2,"curvature",pStart[2],0.01,0,0);
        min->SetParameter(3,"sin_phi",pStart[3],0.01,0,0);
        min->SetParameter(4,"tan_lambda",pStart[4],0.01,0,0);
        min->SetParameter(5,"zc",pStart[5],0.01,0,0);
		//min->SetParameter(5, "zc", pStart[5], 0.01, pStart[5] - 8.25/2., pStart[5] + 8.25/2.); //Limit z to pad length

        Double_t arglist[2];
        arglist[0] = 1000; // number of function calls
        arglist[1] = 0.001; // tolerance
        //arglist[0] = -1.0;
        //min->ExecuteCommand("SET PRINT", arglist, 1);

        Double_t arglist_A[1] = {-1};
        Double_t arglist_B[1] = {0};
        Double_t arglist_C[1] = {-1};

        min->ExecuteCommand("SET PRIntout",arglist_A,1);
        min->ExecuteCommand("SET NOWarnings",arglist_B,1);

        min->ExecuteCommand("MIGRAD",arglist,0);

		Ali_Helix *kHelix_old = new Ali_Helix();
		kHelix_old->setHelix(fHelix[0], fHelix[1], fHelix[2], fHelix[3], fHelix[4], fHelix[5]); 

        for(int i = 0; i < 6; ++i)
        {
			//cout << i << ": " << fHelix[i] << endl;
            fHelix[i] = min->GetParameter(i);	//Set helix with fitted parameters
			//cout << i << ": " << fHelix[i] << endl;
        }
  
		{
			//Use the fitted helix to find the intersection points of the offsets
			
			std::vector<ROOT::Math::SVector<double, 5>> mEstimate_corr;
			mEstimate_corr.resize(6);
			
			//Get the estimates from the helix in the first(?) layer
			TVector3 p_retransform{TMath::Cos(fHelix[2]), TMath::Sin(fHelix[2]), fHelix[3]};	
			TVector3 x_retransform{fHelix[5] + p_retransform[1]/fHelix[4], fHelix[0] - p_retransform[0]/fHelix[4], fHelix[1]};
			x_vek.Print();
			x_retransform.Print();
			p_retransform.RotateZ(-(Double_t)(2 * (Int_t)(mCurrent_Det / 30) + 1) * TMath::Pi() / 18);
			x_retransform.RotateZ(-(Double_t)(2 * (Int_t)(mCurrent_Det / 30) + 1) * TMath::Pi() / 18);
			ROOT::Math::SVector<double, 5> mu_retransform{x_retransform[1], x_retransform[2], p_retransform[1], p_retransform[2], fHelix[4]/b_fak};
			mEstimate_corr[mCurrent_Det%6] = mu_retransform;
			mCurrent_Det++;
			
			//Get the estimates in the other layers
			for(Int_t i_layer = mCurrent_Det%6; i_layer < 6; i_layer++){
				Double_t dist = mTRD_layer_radii_all[mCurrent_Det] - mTRD_layer_radii_all[mCurrent_Det - 1];
				transport(mu_retransform, dist);
				mEstimate_corr[mCurrent_Det%6] = mu_retransform;
				mCurrent_Det++;
			}
	  
			vector<vector<vector<double>>> intersection_pts;	// intersection pts [y,z] of offset in one layer in all three other layers with opposite tilt
			intersection_pts.resize(6);
			Int_t nozerotracklets=0;
			
			for(int i_layer=0; i_layer<6; i_layer++){
				if(mMeasurements[i_layer][0] == 0.) nozerotracklets++;
			}
			
			//cout << "number of zero tracklets: " << nozerotracklets << endl;
			
			//here use either all layers with opposite orientation or just adjacent layers
			//for(auto const& i_layer:{0,2,4}){ //over all even layer
			for(Int_t i_layer = 0; i_layer < 6; i_layer++){
				if(mMeasurements[i_layer] == 0) continue;
				if(mEstimate_corr[i_layer] == ROOT::Math::SVector<double, 5>()) continue;
				
				//for(auto const& j_layer:{1,3,5}){ //over all odd layers
				for(auto const& j_layer:{i_layer - 1, i_layer + 1}){
					if(j_layer < 0) continue;
					if(j_layer > 5) continue;
					if(mMeasurements[j_layer] == 0) continue;
					//if(mEstimate[j_layer] == ROOT::Math::SVector<double, 5>()) continue;
					
					Int_t i_det = mTrack[i_layer]->get_TRD_det();
					ROOT::Math::SVector<double, 4> meas1 = mMeasurements[i_layer]; //measurement in layer i
					ROOT::Math::SVector<double, 5> mu1 = mEstimate_corr[i_layer]; //estimate in layer i
					TVector3 offset1 = mTrack[i_layer]->get_TV3_offset();
					offset1.RotateZ(- (Double_t)(2 * (Int_t)(i_det / 30) + 1) * TMath::Pi() / 18.);
					Int_t i_pad = GetPadRowNumber(i_det, meas1[1]);
					
					Int_t j_det = mTrack[j_layer]->get_TRD_det();
					ROOT::Math::SVector<double, 4> meas2 = mMeasurements[j_layer]; //estimate in layer j
					ROOT::Math::SVector<double, 5> mu2 = mEstimate_corr[j_layer]; //estimate in layer j
					TVector3 offset2 = mTrack[j_layer]->get_TV3_offset();
					offset2.RotateZ(- (Double_t)(2 * (Int_t)(j_det / 30) + 1) * TMath::Pi() / 18.);
					Int_t j_pad = GetPadRowNumber(j_det, meas2[1]);
					
					double dist = offset2[0] - offset1[0];
					//if(!transport(meas1, mu1, dist)) continue;
					//if(!transport(mu1, dist)) continue;
					if(!transport(meas2, mu2, -dist)) continue;
					if(!transport(mu2, -dist)) continue;

					//Offset vector 1
					TVector2 pt1{meas1[0], meas1[1]};

					//Direction vector 1
					Double_t alpha1 = 2. * TMath::Power(-1, i_layer) * TMath::DegToRad();
					TVector3 tv3_dir1 = detector_geometry[3][i_det] * TMath::Cos(alpha1) - detector_geometry[2][i_det] * TMath::Sin(alpha1);
					tv3_dir1.RotateZ(- (Double_t)(2 * (Int_t)(i_det / 30) + 1) * TMath::Pi() / 18.);
					TVector2 tv2_dir1{tv3_dir1[1], tv3_dir1[2]};
					
					//Offset vector 2
					TVector2 pt2{meas2[0], meas2[1]};
					
					//Direction vector 2
					Double_t alpha2 = 2. * TMath::Power(-1, j_layer) * TMath::DegToRad();
					TVector3 tv3_dir2 = detector_geometry[3][j_det] * TMath::Cos(alpha2) - detector_geometry[2][j_det] * TMath::Sin(alpha2);
					tv3_dir2.RotateZ(- (Double_t)(2 * (Int_t)(j_det / 30) + 1) * TMath::Pi() / 18.);
					TVector2 tv2_dir2{tv3_dir2[1], tv3_dir2[2]};
					
					//Calculate intersection pts
					TVector2 pt_intersect = find_Intersect(pt1, tv2_dir1, pt2, tv2_dir2);
					
					//If intersection point is outside of pad
					//Take point at extreme end of pad 
					/*ROOT::Math::SVector<double, 4> meas_intersect = ROOT::Math::SVector<double, 4> (meas2);
					meas_intersect[0] = pt_intersect.X();
					meas_intersect[1] = pt_intersect.Y();
					transport(meas_intersect, mu2, dist);
					TVector2 pt_intersect2{meas_intersect[0], meas_intersect[1]};
					transport(meas2, mu2, dist);
					transport(mu2, dist);
					pt2 = {meas2[0], meas2[1]};
					
					Int_t i_pad = GetPadRowNumber(i_det, pt1.Y());
					Double_t i_pad_up = pad_geometry[i_det][i_pad][0];
					Double_t i_pad_low = pad_geometry[i_det][i_pad][0] - pad_geometry[i_det][i_pad][1];
					
					Int_t j_pad = GetPadRowNumber(j_det, pt2.Y());
					Double_t j_pad_up = pad_geometry[j_det][j_pad][0];
					Double_t j_pad_low = pad_geometry[j_det][j_pad][0] - pad_geometry[j_det][j_pad][1];
					
					if((pt_intersect.Y() > i_pad_up) && (pt_intersect2.Y() < j_pad_low)){
						//Case 1: intersection point above pad 1, below pad 2
						Double_t t1 = (i_pad_up - pt1.Y())/tv2_dir1.Y();
						Double_t t2 = (j_pad_low - pt2.Y())/tv2_dir2.Y();
						
						meas1[0] = (pt1 + t1 * tv2_dir1).X();
						meas1[1] = i_pad_up;
						meas2[0] = (pt2 + t2 * tv2_dir2).X();
						meas2[1] = j_pad_low;
					}else if ((pt_intersect.Y() < i_pad_low) && (pt_intersect2.Y() > j_pad_up)){
						//Case 2: intersection point below pad 1, above pad 2
						Double_t t1 = (i_pad_low - pt1.Y())/tv2_dir1.Y();
						Double_t t2 = (j_pad_up - pt2.Y())/tv2_dir2.Y();
						
						meas1[0] = (pt1 + t1 * tv2_dir1).X();
						meas1[1] = i_pad_low;
						meas2[0] = (pt2 + t2 * tv2_dir2).X();
						meas2[1] = j_pad_up;
					}else if((pt_intersect.Y() > i_pad_up) || (pt_intersect2.Y() > j_pad_up)){
						//Case 3: intersection point above either pad 1 or 2
						//Use the lower value of the upper ends of the pads
						//continue;				
						Double_t t1 = pt_intersect.Y() - i_pad_up;
						Double_t t2 = pt_intersect2.Y() - j_pad_up;
						Double_t t = TMath::Max(t1, t2);
						t1 = t / tv2_dir1.Y();
						t2 = t / tv2_dir2.Y();
						
						meas1[0] = (pt_intersect - t1 * tv2_dir1).X();
						meas1[1] = (pt_intersect - t1 * tv2_dir1).Y();
						meas2[0] = (pt_intersect2 - t2 * tv2_dir2).X();
						meas2[1] = (pt_intersect2 - t2 * tv2_dir2).Y();
					}else if((pt_intersect.Y() < i_pad_low) || (pt_intersect2.Y() < j_pad_low)){
						//Case 4: intersection point below either pad 1 or pad 2
						//Use the higher value of the lower ends of the pads
						//continue;
						Double_t t1 = pt_intersect.Y() - i_pad_low;
						Double_t t2 = pt_intersect2.Y() - j_pad_up;
						Double_t t = TMath::Min(t1, t2);
						t1 = t / tv2_dir1.Y();
						t2 = t / tv2_dir2.Y();
						
						meas1[0] = (pt_intersect - t1 * tv2_dir1).X();
						meas1[1] = (pt_intersect - t1 * tv2_dir1).Y();
						meas2[0] = (pt_intersect2 - t2 * tv2_dir2).X();
						meas2[1] = (pt_intersect2 - t2 * tv2_dir2).Y();
					}else{
						//Case 5: intersection point within pads
						meas1[0] = pt_intersect.X();
						meas1[1] = pt_intersect.Y();
						meas2[0] = pt_intersect2.X();
						meas2[1] = pt_intersect2.Y();
					}*/
					
					
					Double_t delta_z1 = pt_intersect.Y() - pt1.Y();
					Double_t delta_z2 = pt_intersect.Y() - pt2.Y();
					 
					Double_t i_pad_size = pad_geometry[i_det][i_pad][1];
					Double_t j_pad_size = pad_geometry[j_det][j_pad][1];
					
					if((delta_z1 > i_pad_size/2.) && (delta_z2 < -j_pad_size/2.)){
						//Case 1: intersection point above pad 1, below pad 2
						continue; //disregard point
						
						Double_t t1 = i_pad_size/2.;
						Double_t t2 = -j_pad_size/2.;
						
						meas1[0] = (pt1 + t1 * tv2_dir1).X();
						meas1[1] = (pt1 + t1 * tv2_dir1).Y();
						meas2[0] = (pt2 + t2 * tv2_dir2).X();
						meas2[1] = (pt2 + t2 * tv2_dir2).Y();
					}else if ((delta_z1 < -i_pad_size/2.) && (delta_z2 > j_pad_size/2.)){
						//Case 2: intersection point below pad 1, above pad 2
						continue; //disregard point
						
						Double_t t1 = -i_pad_size/2.;
						Double_t t2 = j_pad_size/2.;
						
						meas1[0] = (pt1 + t1 * tv2_dir1).X();
						meas1[1] = (pt1 + t1 * tv2_dir1).Y();
						meas2[0] = (pt2 + t2 * tv2_dir2).X();
						meas2[1] = (pt2 + t2 * tv2_dir2).Y();
					}else if((delta_z1 > i_pad_size/2.) || (delta_z2 > j_pad_size/2.)){
						//Case 3: intersection point above either pad 1 or 2
						//Use the lower value of the upper ends of the pads
						continue; //disregard point
						
						Double_t t1 = TMath::Min(i_pad_size/2., i_pad_size/2. - (pt1.Y() - pt2.Y()));
						Double_t t2 = TMath::Min(j_pad_size/2., j_pad_size/2. - (pt2.Y() - pt1.Y()));
						
						meas1[0] = (pt1 + t1 * tv2_dir1).X();
						meas1[1] = (pt1 + t1 * tv2_dir1).Y();
						meas2[0] = (pt2 + t2 * tv2_dir2).X();
						meas2[1] = (pt2 + t2 * tv2_dir2).Y();
					}else if((delta_z1 < -i_pad_size/2.) || (delta_z2 < -j_pad_size/2.)){
						//Case 4: intersection point below either pad 1 or pad 2
						//Use the higher value of the lower ends of the pads
						continue; //disregard point
						
						Double_t t1 = TMath::Max(-i_pad_size/2., -i_pad_size/2. + (pt1.Y() - pt2.Y()));
						Double_t t2 = TMath::Max(-j_pad_size/2., -j_pad_size/2. + (pt2.Y() - pt1.Y()));
						
						meas1[0] = (pt1 + t1 * tv2_dir1).X();
						meas1[1] = (pt1 + t1 * tv2_dir1).Y();
						meas2[0] = (pt2 + t2 * tv2_dir2).X();
						meas2[1] = (pt2 + t2 * tv2_dir2).Y();
					}else{
						//Case 5: intersection point within pads
						meas1[0] = pt_intersect.X();
						meas1[1] = pt_intersect.Y();
						meas2[0] = pt_intersect.X();
						meas2[1] = pt_intersect.Y();
					}
					
					//Transport intersection point back to layer j (or, alternatively, layer i)
					//transport(meas1, mu1, -dist); 
					//transport(mu1, -dist);
					transport(meas2, mu2, dist);
					transport(mu2, dist);
					
					intersection_pts[i_layer].push_back({meas1[0], meas1[1]});
					intersection_pts[j_layer].push_back({meas2[0], meas2[1]});
				}
			}
			
			
			for(int i_layer=0; i_layer<6; i_layer++){ //Take mean of intersection pts in every layer
				if(!mTrack[i_layer]) continue;
				if(intersection_pts[i_layer].size()==0){
					Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*mTrack[i_layer]);
					mFound_tracks_corr.back()[i_layer]=tracklet;
					//delete tracklet;
					continue;
				}
				mMeasurements[i_layer][0] = 0;
				mMeasurements[i_layer][1] = 0;
				for(unsigned int j=0; j<intersection_pts[i_layer].size();j++){
					mMeasurements[i_layer][0] += intersection_pts[i_layer][j][0];
					mMeasurements[i_layer][1] += intersection_pts[i_layer][j][1];
				}
				mMeasurements[i_layer][0] /= intersection_pts[i_layer].size();
				mMeasurements[i_layer][1] /= intersection_pts[i_layer].size();
				//cout << "New measurement mean " << i_layer << ": " << mMeasurements[i_layer]<< endl;
				
				Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*mTrack[i_layer]);
				TVector3 offset = tracklet->get_TV3_offset();
				Int_t det = mFound_tracks.back()[i_layer]->get_TRD_det();
				offset.RotateZ(- (Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18.);
				offset[1] = mMeasurements[i_layer][0];
				offset[2] = mMeasurements[i_layer][1];
				offset.RotateZ((Double_t)(2 * (Int_t)(det / 30) + 1) * TMath::Pi() / 18.);
				cout << "offset" << endl;
				tracklet->set_TV3_offset(offset);
				//cout << "i_layer: " << offset[0] << " " << offset[1] << " " << offset[2] << endl;
				mFound_tracks_corr.back()[i_layer]=tracklet;
				}
		}
  
		delete min;	
  		
	}
	
    mHelices.push_back(fHelix);
    mChi_2s.push_back(mChi_2);

    if (mShow) {
      cout << "pT: " << pt << endl;
      cout << "mMu[4]: " << mMu[4] << endl;
      cout << "x[1]: " << x[1] << endl;
    }
    
  }
}

void TRD_Kalman_Trackfinder::Write(){
	correctfile->cd();
	correcttree->Write();
	delete correctfile;
}
//vector<vector<Ali_TRD_ST_Tracklets*>> TRD_Kalman_Trackfinder::Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets, Int_t prim_vertex)
//Changed to return uncorrected and correct tracklet offsets (if offsets are not corrected, they are equal)
std::pair<vector<vector<Ali_TRD_ST_Tracklets*>>, vector<vector<Ali_TRD_ST_Tracklets*>>> TRD_Kalman_Trackfinder::Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets, Int_t prim_vertex)
{
  mShow = 0;
  mSearch_tracklets = 1;
  mHelices.clear();
  mFound_tracks.clear();
  mFound_tracks_corr.clear();
  mEstimates.clear();
  mPrimVertex = prim_vertex;
  //cout<<"TRD_Kalman_Trackfinder::Kalman_Trackfind"<<endl;
  T_calc = 0.;
  chrono::high_resolution_clock::time_point start_wall_time = chrono::high_resolution_clock::now();

  get_seed(Tracklets, Num_Tracklets);

  chrono::high_resolution_clock::time_point end_wall_time = chrono::high_resolution_clock::now();
  T_wall = chrono::duration_cast<chrono::microseconds>(end_wall_time - start_wall_time).count();
  //printf("T_wall is %f us\n",T_wall);
  //printf("T_calc is %f us\n",T_calc);

  /*for (int i=0; i<mSeed.size();i++){
     if (i == 55 ) mShow=1;
     Kalman(mSeed[i]);
     mShow=0;
     }*/
  //cout<<mFound_tracks[0][4]->get_TRD_index()<<endl;
  //return mFound_tracks;
  return std::pair<vector<vector<Ali_TRD_ST_Tracklets*>>, vector<vector<Ali_TRD_ST_Tracklets*>>>(mFound_tracks, mFound_tracks_corr);
  //returns both the uncorrected tracks and the corrected tracks (if they are corrected)
  //if no correction is performed, mFound_tracks and mFound_tracks_corr are identical
}

//vector<vector<Ali_TRD_ST_Tracklets*>> TRD_Kalman_Trackfinder::Kalman_Trackfit(vector<vector<Ali_TRD_ST_Tracklets*>> tracks, Int_t prim_vertex)
//Changed to return uncorrected and correct tracklet offsets (if offsets are not corrected, they are equal)
std::pair<vector<vector<Ali_TRD_ST_Tracklets*>>, vector<vector<Ali_TRD_ST_Tracklets*>>> TRD_Kalman_Trackfinder::Kalman_Trackfit(vector<vector<Ali_TRD_ST_Tracklets*>> tracks, Int_t prim_vertex)
{
  mShow = 0;
  mSearch_tracklets = 0;
  mHelices.clear();
  mFound_tracks.clear();
  mFound_tracks_corr.clear();
  mEstimates.clear();
  mPrimVertex = prim_vertex;
  //cout<<"TRD_Kalman_Trackfinder::Kalman_Trackfind"<<endl;
  cout << "entries total: " << tracks.size() << endl;
  //get_seed(Tracklets,Num_Tracklets);
  for (Int_t i_vec = 0; i_vec < (Int_t)tracks.size(); i_vec++) {
    //cout<<"entries now: "<< tracks[i_vec].size()<<endl;
    //tracks[i_vec][0]->get_TV3_offset().Print();
    Kalman(tracks[i_vec]);
  }
  //printf("T_wall is %f us\n",T_wall);
  //printf("T_calc is %f us\n",T_calc);

  /*for (int i=0; i<mSeed.size();i++){
     if (i == 55 ) mShow=1;
     Kalman(mSeed[i]);
     mShow=0;
     }*/
  //cout<<mFound_tracks[0][4]->get_TRD_index()<<endl;
  //return mFound_tracks;
  return std::pair<vector<vector<Ali_TRD_ST_Tracklets*>>, vector<vector<Ali_TRD_ST_Tracklets*>>>(mFound_tracks, mFound_tracks_corr);
   //returns both the uncorrected tracks and the corrected tracks (if they are corrected)
  //if no correction is performed, mFound_tracks and mFound_tracks_corr are identical
}

Int_t TRD_Kalman_Trackfinder::GetPadRowNumber(Int_t det, Double_t z){
  //
  // Finds the pad row number for a given z-position in local supermodule system
  // Taken from AliTRDpadPlane
  //

  Int_t row    = 0;
  Int_t nabove = 0;
  Int_t nbelow = 0;
  Int_t middle = 0;
  Int_t rownum = 0;
  if(pad_geometry[det][12][1] < 0){
	  rownum = 12;
  }else{
	  rownum = 16;
  }

  Double_t Row0 = pad_geometry[det][0][0] + 0.5 * pad_geometry[det][0][1];
  Double_t RowEnd = pad_geometry[det][rownum-1][0] - 0.5 * pad_geometry[det][rownum-1][1];

  if ((z > Row0) || (z < RowEnd)){

    row = -1;

  }else{

	nabove = rownum;
    nbelow = 0;
    while (nabove - nbelow > 0) {
      middle = (nabove + nbelow) / 2;
      if (z == pad_geometry[det][middle][0]) {
        row    = middle;
      }
      if (z  > pad_geometry[det][middle][0] - 0.5 * pad_geometry[det][middle][1]) {
        nabove = middle;
      }
      else {
        nbelow = middle + 1;
      }
    }
    row = nbelow;

  }

  return row;

}

void Chi2_Kalman_tracklet(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t ){
	//Distance between helix track and offset lines
	sum = 0; //sum of distances square to be minimized
	
	//calculate position of helix at layer
	//alternatively: use same t as before (small changes only?)
	//alternatively: calculate minimum distance between helix and tracklet line
	
	Ali_Helix *kalmanHelix = new Ali_Helix();
	kalmanHelix->setHelix(par[0], par[1], par[2], par[3], par[4], par[5]);
	
	//Calculate distance for each layer
	for(Int_t i_layer=0; i_layer < 6; i_layer++){
		//if(i_layer >= 2)continue;
	
		if(!track[i_layer]) continue;
		
		Ali_TRD_ST_Tracklets *tracklet = new Ali_TRD_ST_Tracklets(*track[i_layer]);
		TVector3 offset = tracklet->get_TV3_offset();
		Int_t i_det = tracklet->get_TRD_det();
		Double_t alpha = TMath::Power(-1, i_layer) * 2. * TMath::DegToRad();
		TVector3 dir = det_geom[3][i_det] * TMath::Cos(alpha) - det_geom[2][i_det] * TMath::Sin(alpha);

		
		//Double_t radius_layer_center = 0.5*(TRD_layer_radii[i_layer][1] + TRD_layer_radii[i_layer][2]);
		Double_t radius_layer_center = offset.Perp();

		Double_t track_pos[3];
		TVector3 tv3_track_pos;
		Double_t radius_helix;

		// Find the helix path which touches the first TRD layer
		Double_t track_path_add      = 10.0;
		Double_t track_path_layer0   = 0.0;
		Double_t radius_helix_layer0 = 0.0;
		//for(Double_t track_path = TRD_layer_radii[i_layer][0]; track_path < 1000; track_path += track_path_add)
		
		kalmanHelix->Evaluate(0, track_pos);
		
		//Find helix position at layer
		for(Double_t track_path = 290.-TMath::Sqrt(TMath::Power(track_pos[0], 2)+TMath::Power(track_pos[1], 2)); track_path < 300; track_path += track_path_add)
		{
			kalmanHelix ->Evaluate(track_path,track_pos);; 
			//track_path, track_pos);
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
				break;
			}
		}
		
		tv3_track_pos = (TVector3) track_pos;
	
		TVector3 delta = tv3_track_pos - offset; //distance between helix and offset point
		Double_t dist = TMath::Sqrt(delta*delta - (delta*dir)*(delta*dir)/(dir*dir)); //minimal distance between helix and line through offset point
		
		//cout << i_layer << " dist w/o penalty: " << dist << endl;
		//delta.Print();
		//Use penalty function if outside of pad, i.e. distance in z direction larger than half the pad length
		Double_t l_z = 8.25;
		if(TMath::Abs(delta[2])>l_z/2.){
			dist += TMath::Exp(TMath::Power((delta[2] - l_z/2.)/.2, 2.))-1;
			dist += (delta[2] - l_z/2.)*(delta[2] - l_z/2.)*10;
		}
		sum+= dist*dist;
	
		delete tracklet;
	}
	delete kalmanHelix;
}
