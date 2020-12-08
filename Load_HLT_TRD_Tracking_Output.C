#if !defined(__CLING__) || defined(__ROOTCLING__)
using namespace std;
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "TString.h"

#include "TLorentzVector.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TMath.h"
#include "TCanvas.h"
//#include "TStatToolkit.h"
#include "TString.h"
#include "TVector3.h"
#include "Ali_TRD_ST.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#endif


TFile* f = 0x0;    // input file (TRDhlt.root)
TTree* tree = 0x0; // input tree (tracksFinal)

// branches
Int_t event = 0x0;
Int_t nTracklets = 0x0;
Float_t chi2Total = 0x0;
Int_t nFake = 0x0;
Int_t trackID = 0x0;
TVectorF* update = 0x0;
TVectorF* trackX = 0x0;
TVectorF* trackY = 0x0;
TVectorF* trackYerr = 0x0;
TVectorF* trackPhi = 0x0;
TVectorF* trackLambda = 0x0;
TVectorF* trackQPt = 0x0;
TVectorF* trackZ = 0x0;
TVectorF* trackZerr = 0x0;
TVectorF* trackSec = 0x0;
TVectorF* trackletX = 0x0;
TVectorF* trackletY = 0x0;
TVectorF* trackletZ = 0x0;
TVectorF* trackletYerr = 0x0;
TVectorF* trackletYZerr = 0x0;
TVectorF* trackletZerr = 0x0;
TVectorF* trackletYRaw = 0x0;
TVectorF* trackletZRaw = 0x0;
TVectorF* trackletDy = 0x0;
TVectorF* trackletDet = 0x0;
//

// functions
Bool_t InitAnalysis(const char* filename = "TRDhlt.root", Bool_t isMC = kTRUE);
void Reset();

Bool_t InitTree(const char* filename)
{
  f = new TFile(filename, "open");
  if (!f) {
    printf("File could not be opened\n");
    return kFALSE;
  }
  tree = (TTree*)f->Get("tracksFinal");
  if (!tree) {
    printf("Tree could not be opened\n");
    return kFALSE;
  }
  return kTRUE;
}

void InitBranches()
{
  if (!tree) {
    return;
  }

  // for explanations of the variables see GPUTRDTrackerDebug.h
  tree->SetBranchAddress("update.", &update);
  tree->SetBranchAddress("trackX.", &trackX);
  tree->SetBranchAddress("trackY.", &trackY);
  tree->SetBranchAddress("trackYerr.", &trackYerr);
  tree->SetBranchAddress("trackPhi.", &trackPhi);
  tree->SetBranchAddress("trackLambda.", &trackLambda);
  tree->SetBranchAddress("trackQPt.", &trackQPt);
  tree->SetBranchAddress("trackZ.", &trackZ);
  tree->SetBranchAddress("trackZerr.", &trackZerr);
  tree->SetBranchAddress("trackSec.", &trackSec);

  tree->SetBranchAddress("trackletX.", &trackletX);
  tree->SetBranchAddress("trackletY.", &trackletY);
  tree->SetBranchAddress("trackletDy.", &trackletDy);
  tree->SetBranchAddress("trackletZ.", &trackletZ);
  tree->SetBranchAddress("trackletYerr.", &trackletYerr);
  tree->SetBranchAddress("trackletZerr.", &trackletZerr);
  tree->SetBranchAddress("trackletYZerr.", &trackletYZerr);
  tree->SetBranchAddress("trackletYRaw.", &trackletYRaw);
  tree->SetBranchAddress("trackletZRaw.", &trackletZRaw);
  tree->SetBranchAddress("trackletDet.", &trackletDet);

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("nTracklets", &nTracklets);
  tree->SetBranchAddress("chi2Total", &chi2Total);
  tree->SetBranchAddress("nFake", &nFake);
  tree->SetBranchAddress("trackID", &trackID);

}

//Speicher Aufräumen!
vector< vector<Ali_TRD_ST_Tracklets*> > get_found_tracklets_format()
{
	vector< vector<Ali_TRD_ST_Tracklets*> > all_tracks;
	Int_t index=0;
	for (Int_t iTrk = 0; iTrk < tree->GetEntriesFast(); ++iTrk) {
  		tree->GetEntry(iTrk);
		if (trackID < 0) {
		  // no valid MC ID for this track
		  continue;
		}
		if (nTracklets < 1) {
		  // skip tracks without any TRD tracklets
		  continue;
		}
		if ((nFake > 0) || (nTracklets < 6)) {
		  // skip tracks which do not have 6 tracklets with MC truth match
		  continue;
		}
		vector<Ali_TRD_ST_Tracklets*> current_track;
		//current_track.resize(6):	
		for (Int_t iLy = 5; iLy >= 0 ; iLy--) {
		  // loop over all TRD layers
		  	if ((*update)[iLy] < 1) {
			// no tracklet available in this layer
				continue;
		 	}
		  //float trkY = (*trackY)[iLy];
		  //float trkltY = (*trackletY)[iLy];
		 // hResY->Fill(trkY - trkltY);
			TVector3 offset = (TVector3){(*trackletX)[iLy], (*trackletY)[iLy], (*trackletZ)[iLy]};
			
			TVector3 direction = (TVector3){3, (*trackletDy)[iLy], (*trackLambda)[iLy] *3};
			
			Ali_TRD_ST_Tracklets* temp_tracklet	=new Ali_TRD_ST_Tracklets();
					
			Int_t det							=(*trackletDet)[iLy];
			temp_tracklet						->set_TRD_det(det);
			temp_tracklet						->set_TV3_offset(offset);
			temp_tracklet						->set_TV3_dir(direction);
			temp_tracklet						->set_TRD_index(index);
			index++;
			current_track.push_back(temp_tracklet);
			
		}
		all_tracks.push_back(current_track);
		
	}
	return all_tracks;
}	

void PlotResidualsY()
{
  TH1F* hResY = new TH1F("resY", ";#Delta y (cm);counts", 100, -3, 3);
  for (Int_t iTrk = 0; iTrk < tree->GetEntriesFast(); ++iTrk) {
    // loop over all tracks
    tree->GetEntry(iTrk);
    if (trackID < 0) {
      // no valid MC ID for this track
      continue;
    }
    if (nTracklets < 1) {
      // skip tracks without any TRD tracklets
      continue;
    }
    if (nFake > 0 && nTracklets < 6) {
      // skip tracks which do not have 6 tracklets with MC truth match
      continue;
    }
    for (Int_t iLy = 0; iLy < 6; iLy++) {
      // loop over all TRD layers
      if ((*update)[iLy] < 1) {
        // no tracklet available in this layer
        continue;
      }
      float trkY = (*trackY)[iLy];
      float trkltY = (*trackletY)[iLy];
      hResY->Fill(trkY - trkltY);
    }
  }
  TCanvas *c1 = new TCanvas("c1", "c1");
  hResY->Draw();
}


Bool_t InitAnalysis(const char* filename, Bool_t isMC)
{
  Reset();
  if (!InitTree(filename)) {
    return kFALSE;
  }
  //
  InitBranches();
  //
  return kTRUE;
}

void Reset()
{
  if (f) {
    delete f;
    f = 0x0;
  }
  tree = 0x0;
}

void checkDbgOutput()
{
  printf("Basic usage:\n");
  printf("InitAnalysis(\"TRDhlt.root\", 1);\n");
}
