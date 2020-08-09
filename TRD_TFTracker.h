#ifndef TRD_TFTracker
#define TRD_TFTracker

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
#include "Ali_TRD_ST.h"
#include "Ali_TRD_ST_LinkDef.h"


//#include <TPython.h>



//#include "TRD_Kalman_Tracking_Source.cxx"

class TRD_TFTrackMaker
{
public:
    TRD_TFTrackMaker();
    Int_t  Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets);
    vector<vector<Double_t>> get_Kalman_helix_params();
    vector<vector<Ali_TRD_ST_Tracklets*>> get_Tracklets();

private:
    vector<vector<Ali_TRD_ST_Tracklets*>> mFound_tracklets;
    vector<vector<Double_t>> mHelices;
};

#endif
