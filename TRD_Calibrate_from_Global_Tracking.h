// some useful includes

#include "AliTRDgeometry.h" // TRD geometry
#include "AliExternalTrackParam.h" // ALICE track model
#include "AliTrackerBase.h" // track propagation

#include "TTree.h"
#include "TFile.h"
#include "TVectorF.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGrid.h"

// database access for geometry, magnetic field etc.
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliGRPManager.h"
#include "functions.h"

// c++ stuff
#include <vector>
#include <iostream>
#include <math.h>
//


struct TRDTracklet {
  float x;
  float y;
  float z;
  float dy;
  //float tanphi = dy/3.;
  float sinphi;
  float sigy;
  float sigz;
  float sigyz;
  int det;
};


void init()
{
  // load the calibration database for the run which you are analysing (on the HD servers we have cvmfs installed, so no grid connection is needed)
  auto man = AliCDBManager::Instance();
  man->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2015/OCDB");
  //TGrid::Connect("alien");
  //man->SetDefaultStorage("raw://");
  man->SetRun(245353);
  // load the geometry and take alignment into account
  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF");
  // load magnetic field
  AliGRPManager grpMan;
  grpMan.ReadGRPEntry();
  grpMan.SetMagField();
}

bool loadInput(std::vector<AliExternalTrackParam>& tracks, std::vector<TRDTracklet>& tracklets, std::vector<int>& indices, const char* filename = "TRDhlt.root")
{
  // open file TRDhlt.root and load the tracks and tracklets into the vectors
  TTree* tree = 0x0;
  TVectorF* trackletX = 0x0;
  TVectorF* trackletY = 0x0;
  TVectorF* trackletZ = 0x0;
  TVectorF* trackletDy = 0x0;
  TVectorF* trackSec = 0x0;
  TVectorF* trackletYerr = 0x0;
  TVectorF* trackletZerr = 0x0;
  TVectorF* trackletYZerr = 0x0;
  TVectorF* update = 0x0;
  TVectorF* trackletDet = 0x0;
  TVectorF* trackPhi = 0x0;
  TVectorF* trackLambda = 0x0;
  TVectorF* trackQPt = 0x0;

  Int_t nTracklets = 0x0;
  Int_t nFake = 0x0;
  Int_t trackID = 0x0;

  TFile* f = 0x0;
  f = new TFile(filename, "open");
  if (!f) {
    printf("File could not be opened\n");
    return false;
  }
  tree = (TTree*)f->Get("tracksFinal");
  if (!tree) {
    printf("Tree could not be opened\n");
    return false;
  }

  tree->SetBranchAddress("trackletX.", &trackletX);
  tree->SetBranchAddress("trackletY.", &trackletY);
  tree->SetBranchAddress("trackletDy.", &trackletDy);
  tree->SetBranchAddress("trackletZ.", &trackletZ);
  tree->SetBranchAddress("trackSec.", &trackSec);
  //tree->SetBranchAddress("nFake", &nFake);
  tree->SetBranchAddress("trackID", &trackID);
  tree->SetBranchAddress("nTracklets", &nTracklets);
  tree->SetBranchAddress("trackletYerr.", &trackletYerr);
  tree->SetBranchAddress("trackletZerr.", &trackletZerr);
  tree->SetBranchAddress("trackletYZerr.", &trackletYZerr);
  tree->SetBranchAddress("update.", &update);
  tree->SetBranchAddress("trackletDet.", &trackletDet);
  tree->SetBranchAddress("trackPhi.", &trackPhi);
  tree->SetBranchAddress("trackLambda.", &trackLambda);
  tree->SetBranchAddress("trackQPt.", &trackQPt);


int index = 0;
  for (Int_t iTrk = 0; iTrk < tree->GetEntriesFast(); ++iTrk) {
  //for (Int_t iTrk = 0; iTrk < 1000; ++iTrk) {
    // loop over all tracks
    tree->GetEntry(iTrk);
    if (trackID < 0) {
      // no valid MC ID for this track
      continue;
    }
    if (nTracklets < 4) {
      // skip tracks without any TRD tracklets
      continue;
    }
    if (nFake > 0 ) {
      // skip tracks which do not have 6 tracklets with MC truth match
      continue;
    }
    indices.push_back(index);

    // for example create a track
    // take the variables from the input and set them for the current track
    float param[5] = { 0. }; // put here the track parameters y, z, sin(phi), tan(lambda), q/pt
    float x = 0; // current radial position of the track
    float alpha = 0; // angle to convert from sector coordinates into global coordinates
    float covariance[15] = { 0. }; // covariance matrix, set these to very high values
    covariance[0]=100*100;
    covariance[2]=100*100;
    covariance[5]=1*1;
    covariance[9]=1*1;
    covariance[14]=100*100;

    float dy = 0;
    int det_temp = 0;


    for (Int_t iLy = 0; iLy < 6; iLy++) {
      // loop over all TRD layers
      if ((*update)[iLy] < 1) { //no tracklet!
        // no tracklet available in this layer
        continue;
      }
      TRDTracklet trackletIn;
      det_temp = (*trackletDet)[iLy];
      alpha = (float)((2. * ((int)(det_temp / 30)) + 1) * 3.14159265359 / 18.);
      trackletIn.x      = (*trackletX)[iLy];
      trackletIn.y      = (*trackletY)[iLy];
      trackletIn.z      = (*trackletZ)[iLy];
      param[0]          = (*trackletY)[iLy];
      param[1]          = (*trackletZ)[iLy];
      param[2]          = (*trackPhi)[iLy];
      param[3]          = (*trackLambda)[iLy];
      param[4]          = (*trackQPt)[iLy];
      x                 = (*trackletX)[iLy];
      dy                = (*trackletDy)[iLy];
      trackletIn.dy     = dy;
      trackletIn.sigy   = (*trackletYerr)[iLy];
      trackletIn.sigz   = (*trackletZerr)[iLy];
      trackletIn.sigyz  = (*trackletYZerr)[iLy];
      trackletIn.det    = det_temp;
      //trackletIn.tanphi = dy/3.;
      trackletIn.sinphi = dy/(dy*dy + 9);


      tracklets.emplace_back(std::move(trackletIn));
      index++;

    }
    AliExternalTrackParam trackIn;

    trackIn.Set(x, alpha, param, covariance);
    // move the track into the vector
    tracks.emplace_back(std::move(trackIn));



  }
  indices.push_back(index);

  // you will need an additional helper for keeping track of the number of tracklets which belong to each track

  return true;
}


void mainFunction()
{
  // do initialization
  init();

  auto geo = new AliTRDgeometry();

  AliTrackerBase propagator;

  std::vector<AliExternalTrackParam> tracksIn; // input tracks
  std::vector<TRDTracklet> trackletsIn; // input tracklets
  std::vector<int> indices;

  vector<TProfile*> vec_tp_Delta_vs_impact(540);
  for (int i_det = 0; i_det < 540; i_det++)
    vec_tp_Delta_vs_impact[i_det] = new TProfile(Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),Form("vec_th1d_Delta_vs_impact_circle_%d",i_det),360,-360,360);



  if (!loadInput(tracksIn, trackletsIn,indices)) {
    printf("Error: no input found\n");
    return;
  }

  printf("nbr tracks %d,nbr indizes %d\n",tracksIn.size(),indices.size() );

  int nbr_angles=0;
  //for (auto& track : tracksIn) {
  for (int i_track = 0; i_track < tracksIn.size(); i_track++  ) {
    auto& track = tracksIn[i_track];
    bool successful_inward_prop = 1;
    bool successful_outward_prop = 1;

    // track loop

    //refit each track using the TRD tracklet information only
    //for (int iLayer = 5; iLayer >= 0; --iLayer) {
      //float xLayer = geo->GetTime0(iLayer);
    for (int i_tracklet = indices.at(i_track+1)-1; i_tracklet >= indices.at(i_track); --i_tracklet  ) {
      // propagate the track to the radius of xLayer, assume pion mass 0.139 GeV, maximum step size 2cm and don't rotate the track into the tracklet frame
      // this takes into acount the inhomogenic magnetic field, energy loss and MS contribution
      //if (!propagator.PropagateTrackToBxByBz(&track, trackletsIn[i_tracklet].x, 0.139, 2., kFALSE)) {
      int ret_val = propagator.PropagateTrackTo2(&track, trackletsIn[i_tracklet].x, 0.139, 2., kFALSE);
      if(ret_val!=1){
        printf("Warning: track failed inwards propagation step:ret_val = %d\n",ret_val);
        printf("x_val_tracklet = %f;x_val_track = %f\n",trackletsIn[i_tracklet].x,track.GetX());
        printf("y_val_tracklet = %f;y_val_track = %f\n",trackletsIn[i_tracklet].y,track.GetY());
        printf("z_val_tracklet = %f;z_val_track = %f\n",trackletsIn[i_tracklet].z,track.GetZ());
        printf("Layer %d of %d  \n",i_tracklet-indices[i_track],indices[i_track+1]-indices[i_track]);
        printf("index %d between %d and %d \n",i_tracklet,indices[i_track],indices[i_track+1]);
        successful_inward_prop = 0;
        //i_track++;
        break;
      }
      // TODO check if there is a tracklet in this layer, if not then skip it!
      // if there is probably better to propagate directly to the x of the tracklet..
      Double_t yz[2] = { trackletsIn[i_tracklet].y, trackletsIn[i_tracklet].z }; // put here the y and z position of the tracklet
      Double_t cov[3] = { trackletsIn[i_tracklet].sigy  ,trackletsIn[i_tracklet].sigyz ,trackletsIn[i_tracklet].sigz }; // and here the covariance 0:sigy^2 1:sigyz 2:sigz^2
      if (!track.Update(yz, cov)) {
        printf("Warning: track failed inwards update step\n");
        printf("x_val_tracklet = %f;x_val_track = %f\n",trackletsIn[i_tracklet].x,track.GetX());
        printf("y_val_tracklet = %f;y_val_track = %f\n",trackletsIn[i_tracklet].y,track.GetY());
        printf("z_val_tracklet = %f;z_val_track = %f\n",trackletsIn[i_tracklet].z,track.GetZ());
        printf("Layer %d of %d  \n",i_tracklet-indices[i_track],indices[i_track+1]-indices[i_track]);
        printf("index %d between %d and %d \n",i_tracklet,indices[i_track],indices[i_track+1]);
        successful_inward_prop = 0;
        //i_track++;
        break;
      }
    }

    if(!successful_inward_prop){
      //i_track++;
      continue;
    }
    // after the refit we can compare the the track to the tracklet angle at each layer
    //for (int iLayer = 0; iLayer <= 5; ++iLayer) {
    for (int i_tracklet = indices.at(i_track); i_tracklet < indices.at(i_track+1); i_tracklet++  ) {
      // propagate to x of layer and compare tracklet slope with track.GetSinPhi()
      //if (!propagator.PropagateTrackToBxByBz(&track, trackletsIn[i_tracklet].x, 0.139, 2., kFALSE)) {
      int ret_val = propagator.PropagateTrackTo2(&track, trackletsIn[i_tracklet].x, 0.139, 2., kFALSE);
      if(ret_val!=1){
        printf("Warning: track failed outwards propagation step:ret_val = %d\n",ret_val);
        printf("x_val_tracklet = %f;x_val_track = %f\n",trackletsIn[i_tracklet].x,track.GetX());
        successful_outward_prop=0;
        //i_track++;
        break;
      }

      Double_t yz[2] = { trackletsIn[i_tracklet].y, trackletsIn[i_tracklet].z }; // put here the y and z position of the tracklet
      Double_t cov[3] = { trackletsIn[i_tracklet].sigy  ,trackletsIn[i_tracklet].sigyz ,trackletsIn[i_tracklet].sigz }; // and here the covariance 0:sigy^2 1:sigyz 2:sigz^2
      if (!track.Update(yz, cov)) {
        printf("Warning: track failed outwards update step\n");
        successful_outward_prop=0;
        //i_track++;
        break;
      }

      // fill histograms
      float impact_angle = std::asin(track.GetSnp())*180/3.14159265359;
      float angle_diff = (std::asin(trackletsIn[i_tracklet].sinphi)*180/3.14159265359 - impact_angle);
      impact_angle +=90;
      vec_tp_Delta_vs_impact[trackletsIn[i_tracklet].det]     ->Fill(impact_angle,angle_diff);
      nbr_angles++;

      //printf("detector %d; impact_ange %f; angle_diff %f\n",(int)trackletsIn[i_tracklet].det,impact_angle, angle_diff);



    }
    //if(successful_inward_prop)
      //printf("\nSUCCESS\n\n");

    //i_track++;
  }
  //printf("nbr angles %d, tracks %d\n",nbr_angles,i_track );

  // check the histograms
  vector<TCanvas*> vec_can_Delta_vs_impact;
  char NoP[50];
  Int_t arr_color_layer[6] = {kBlack,kRed,kBlue,kGreen,kMagenta,kCyan};

  vec_can_Delta_vs_impact.resize(6); // 6 sector blocks with 3 sectors each (18)

  TH1D* h_dummy_Delta_vs_impact = new TH1D("h_dummy_Delta_vs_impact","h_dummy_Delta_vs_impact",90,50,140);
  h_dummy_Delta_vs_impact->SetStats(0);
  h_dummy_Delta_vs_impact->SetTitle("");
  h_dummy_Delta_vs_impact->GetXaxis()->SetTitleOffset(0.85);
  h_dummy_Delta_vs_impact->GetYaxis()->SetTitleOffset(0.78);
  h_dummy_Delta_vs_impact->GetXaxis()->SetLabelOffset(0.0);
  h_dummy_Delta_vs_impact->GetYaxis()->SetLabelOffset(0.01);
  h_dummy_Delta_vs_impact->GetXaxis()->SetLabelSize(0.08);
  h_dummy_Delta_vs_impact->GetYaxis()->SetLabelSize(0.08);
  h_dummy_Delta_vs_impact->GetXaxis()->SetTitleSize(0.08);
  h_dummy_Delta_vs_impact->GetYaxis()->SetTitleSize(0.08);
  h_dummy_Delta_vs_impact->GetXaxis()->SetNdivisions(505,'N');
  h_dummy_Delta_vs_impact->GetYaxis()->SetNdivisions(505,'N');
  h_dummy_Delta_vs_impact->GetXaxis()->CenterTitle();
  h_dummy_Delta_vs_impact->GetYaxis()->CenterTitle();
  h_dummy_Delta_vs_impact->GetXaxis()->SetTitle("impact angle");
  h_dummy_Delta_vs_impact->GetYaxis()->SetTitle("#Delta #alpha");
  h_dummy_Delta_vs_impact->GetXaxis()->SetRangeUser(70,110);
  h_dummy_Delta_vs_impact->GetYaxis()->SetRangeUser(-24,24);

  for(Int_t i_sec_block = 0; i_sec_block < 6; i_sec_block++) //<6
  {
      TString HistName = "vec_can_Delta_vs_impact_";
      HistName += i_sec_block;
      vec_can_Delta_vs_impact[i_sec_block] = new TCanvas(HistName.Data(),HistName.Data(),10,10,1600,1000);

      vec_can_Delta_vs_impact[i_sec_block] ->Divide(5,3); // x = stack, y = sector

      for(Int_t i_sec_sub = 0; i_sec_sub < 3; i_sec_sub++)
      {
          Int_t i_sector = i_sec_block + 6*i_sec_sub;
          for(Int_t i_stack = 0; i_stack < 5; i_stack++)
          {
              Int_t iPad = i_sec_sub*5 + i_stack + 1;
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTicks(1,1);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetGrid(0,0);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetFillColor(10);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetRightMargin(0.01);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetTopMargin(0.01);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetBottomMargin(0.2);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad)->SetLeftMargin(0.2);
              vec_can_Delta_vs_impact[i_sec_block] ->cd(iPad);
              h_dummy_Delta_vs_impact->Draw("h");

              for(Int_t i_layer = 0; i_layer < 6; i_layer++)
              {
                  Int_t i_detector = i_layer + 6*i_stack + 30*i_sector;
                  // printf("detector: %d \n",i_detector);
                  vec_tp_Delta_vs_impact[i_detector] ->SetLineColor(arr_color_layer[i_layer]);
                  vec_tp_Delta_vs_impact[i_detector] ->SetLineWidth(2);
                  vec_tp_Delta_vs_impact[i_detector] ->SetLineStyle(1);
                  vec_tp_Delta_vs_impact[i_detector] ->Draw("same hl");

                  TString HistName2 = "";
                  sprintf(NoP,"%4.0f",(Double_t)i_detector);
                  HistName2 += NoP;
                  plotTopLegend((char*)HistName2.Data(),0.24,0.89-i_layer*0.07,0.045,arr_color_layer[i_layer],0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
              }
          }
      }
      HistName += "_sim";
      HistName += ".png";

       vec_can_Delta_vs_impact[i_sec_block]->SaveAs(HistName.Data());
  }

  //TFile* outputfile_hist = new TFile("vec_tp_Delta_vs_impact","RECREATE");

  printf("Write data to output file \n");
  TFile* outputfile_hist = new TFile("./angle_detector_hit_sim.root","RECREATE");
  // h_detector_hit_outputfile ->cd();
  //  h_detector_hit->Write();

  // THIS NEEDED
  outputfile_hist ->cd();
  //outputfile_hist ->mkdir("Delta_impact");
  //outputfile_hist ->cd("Delta_impact");
  for(Int_t i_det = 0; i_det < 540; i_det++)
      vec_tp_Delta_vs_impact[i_det]   ->Write();

  delete geo;
  printf("Done\n");
}

void TRD_Calibrate_from_Global_Tracking()
{
  mainFunction();
}
