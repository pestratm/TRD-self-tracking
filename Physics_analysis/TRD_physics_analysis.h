
#include "../Ali_TRD_Self_Event.h" 
#include "../Ali_TRD_Self_EventLinkDef.h"

ClassImp(Ali_Kalman_Track)
ClassImp(Ali_TPC_Track)
ClassImp(Ali_TRD_Photon)
ClassImp(Ali_TRD_Nuclear_interaction)
ClassImp(Ali_TRD_Self_Event)

class Ali_TRD_physics_analysis
{
private:
    TChain* input_SE;

    TString TRD_ST_TREE   = "Tree_TRD_Self_Event";
    TString TRD_ST_BRANCH = "Tree_TRD_Self_Event_branch";
    
 	Long64_t file_entries_total;
    Long64_t N_Events;

    //stuff for new classes 
    Ali_Kalman_Track* TRD_Kalman_track;
    Ali_TPC_Track* TPC_track;
    Ali_TRD_Self_Event*   TRD_Self_Event;
    Ali_TRD_Photon* TRD_Photon;
    Ali_TRD_Nuclear_interaction* TRD_Nuclear_interaction;

    
    TFile* outputfile;
   

    TString HistName;
    TH1D* th1d_TRD_layer_radii;
    
    TString input_dir;
    TString input_dir_lists;
    

public:
    Ali_TRD_physics_analysis(TString out_dir, TString out_file_name);
    //~Ali_TRD_ST_Analyze();

    void Init_tree(TString SEList);

    void set_input_lists(TString input_dir_lists_in) {input_dir_lists = input_dir_lists_in;}
    void set_input_dir(TString input_dir_in) {input_dir = input_dir_in;}
    
    ClassDef(Ali_TRD_physics_analysis, 1)
};
//----------------------------------------------------------------------------------------

