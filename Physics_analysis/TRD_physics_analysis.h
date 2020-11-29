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

#include "../Ali_TRD_Self_Event.h" 
#include "../Ali_TRD_Self_EventLinkDef.h"

ClassImp(Ali_Kalman_Track)
ClassImp(Ali_TPC_Track)
ClassImp(Ali_TRD_Photon)
ClassImp(Ali_TRD_Nuclear_interaction)
ClassImp(Ali_TRD_Self_Event)
ClassImp(Ali_Helix)

class Ali_TRD_physics_analysis
{
private:
    TChain* input_SE;

    TString TRD_Self_TREE   = "ST_Physics/Tree_TRD_Self_Event";
    TString TRD_Self_BRANCH = "Tree_TRD_Self_Event_branch";
    
 	Long64_t file_entries_total;
    Long64_t N_Events;

    //stuff for new classes 
    Ali_Kalman_Track* Kalman_Track_photon;
    Ali_TPC_Track* TPC_Track_photon;    
    Ali_Kalman_Track* Kalman_Track_interact;
    Ali_TPC_Track* TPC_Track_interact;
    Ali_TRD_Self_Event*   TRD_Self_Event;
    Ali_TRD_Photon* TRD_Photon;
    Ali_TRD_Nuclear_interaction* TRD_Nuclear_interaction;

    
    TFile* outputfile;
   

    TString HistName;
    TH1D* th1d_TRD_layer_radii;
    
    TString input_dir;
    TString input_dir_lists;

    Double_t EventVertexX = -999.0;
    Double_t EventVertexY = -999.0;
    Double_t EventVertexZ = -999.0;
    Long64_t Global_Event = -999;
    Int_t    Global_RunID = -999;
    TVector3 TV3_EventVertex;

    //things for Photon convertions
    TVector3 TV3_PhotonVertex; //temporary thing
    vector< TVector3 > vec_PhotonVertex; //[i_photon]
    
    vector< vector< Double_t>> vec_photon_kalman_chi2; //[i_photon][i_track]
    vector< vector< Ali_Helix* >> vec_photon_kalman_helices; //[i_photon][i_track]
    vector< vector< Ali_Helix* >> vec_photon_tpc_helices; //[i_photon][i_track]

    //things for Nuclear interactions
    TVector3 TV3_NIVertex;
    vector< TVector3 > vec_NIVertex;

    vector< vector< Double_t>> vec_ni_kalman_chi2; //[i_interaction][i_track]
    vector< vector< Ali_Helix* >> vec_ni_kalman_helices; //[i_interaction][i_track]
    vector< vector< Ali_Helix* >> vec_ni_tpc_helices; //[i_interaction][i_track]

    //things for Nuclear ingeractions
    
public:
    Ali_TRD_physics_analysis(TString out_dir, TString out_file_name);
    //~Ali_TRD_ST_Analyze();

    void Init_tree(TString SEList);
    Int_t Loop_event(Long64_t i_event);


    void set_input_lists(TString input_dir_lists_in) {input_dir_lists = input_dir_lists_in;}
    void set_input_dir(TString input_dir_in) {input_dir = input_dir_in;}

    Long64_t get_N_Events() {return N_Events;}
    
    ClassDef(Ali_TRD_physics_analysis, 1)
};
//----------------------------------------------------------------------------------------

