#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

//------------------------
#include "AliHelix.h"
#include "TLorentzVector.h"
#include "TSystem.h"
//------------------------

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliKalmanTrack.h"

#include "AliTRDpadPlane.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"
#include "AliESDfriend.h"

#include "AliTRDdigitsManager.h"
#include "AliTRDarrayADC.h"

#include "AliPIDResponse.h"
#include "AliPID.h"

#include "AliESDtrackCuts.h"

#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliESDRun.h"

#include "AliMultSelection.h"

#include "AliCDBEntry.h"
#include "TClonesArray.h"
#include "TGeoMatrix.h"
#include "AliAlignObjParams.h"

#include "AliTRDdigitsParam.h"
#include "AliRunTag.h"
#include "TObjString.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"
#include "AliTRDCalOnlineGainTable.h"
#include "AliTRDCalROC.h"
#include "TPolyMarker.h"

#include "AliTRDCommonParam.h"

#include "AliESDTrdTracklet.h"
#include "TProfile.h"
#include "TGraph.h"

//------------------------
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
//------------------------

#include "Ali_make_tracklets_from_digits.h"


#include <iostream>
#include <iomanip>
using namespace std;

static Int_t flag_plot_event = 0;
static TString HistName;

static TFile* dfile;
static TFile* TRD_alignment_file;
static TFile* TRD_calibration_file_AA;
static TGraph* tg_v_fit_vs_det;
static TGraph* tg_LA_factor_fit_vs_det;
static TH1D* h_v_fit_vs_det;
static TH1D* h_LA_factor_fit_vs_det;
static TGeoHMatrix TM_TRD_rotation_det[540];
static TGeoHMatrix TM_TRD_rotation_sector[18];
static TVector3    TV3_TRD_translation[540];
static const Double_t TRD_lorentz_angle = TMath::DegToRad()*8.8;
static AliTRDCalDet *ChamberVdrift;
static AliTRDCalDet *ChamberT0;
static AliTRDCalDet *ChamberExB;
static AliTRDCalDet *chambergain;
static AliTRDCalOnlineGainTable *KryptoGain;
static AliTRDCalPad *PadNoise;
static AliTRDCalPad *LocalT0_pad;
static AliTRDCalROC* CalROC;
static AliTRDCalROC  *LocalT0;            //  Pad wise T0 calibration object
static const Int_t N_pT_bins = 5;
static const Double_t pT_ranges[N_pT_bins+1] = {0.2,0.5,1.0,2.0,3.0,5.0};
static const Double_t TRD_Impact_distance_in_drift = 3.35;

static AliTRDCommonParam* fParam;
static AliESDfriend *esdFr = NULL;
static AliESDInputHandler *esdH = NULL;
static TString esdFriendTreeFName;
//static Class_peak_finder my_class_peak_finder;

ClassImp(Ali_AS_Event)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_make_tracklets_from_digits)



//----------------------------------------------------------------------------------------
TVector3 intersect_line_plane(TVector3 TV3_base_line, TVector3 TV3_dir_line, TVector3 TV3_base_plane, TVector3 TV3_norm_plane)
{

    TVector3 TV3_base_diff 		= TV3_base_plane - TV3_base_line;
    Double_t TV3_norm_dir_mult 	= TV3_dir_line.Dot(TV3_norm_plane);
    Double_t TV3_norm_diff_mult = TV3_base_diff.Dot(TV3_norm_plane);

    Double_t x 					= TV3_norm_diff_mult / TV3_norm_dir_mult;

    TVector3 intersect_point 	= TV3_base_line + x * TV3_dir_line;

    return intersect_point;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// Driver function to sort the 2D vector 
// on basis of a particular column 
bool sortcol_first( const vector<Double_t>& v1,
             const vector<Double_t>& v2 )
{
 return v1[3] > v2[3];  // second column
} 
//----------------------------------------------------------------------------------------


//________________________________________________________________________
Ali_make_tracklets_from_digits::Ali_make_tracklets_from_digits(const char *name)
: AliAnalysisTaskSE(name),
fLocalMode(kFALSE),
fDigitsInputFileName("TRD.FltDigits.root"), fDigitsInputFile(0),
fDigitsOutputFileName(""), fDigitsOutputFile(0),
fDigMan(0),fGeo(0),AS_Event(0),AS_Track(0),AS_Tracklet(0),AS_offline_Tracklet(0),AS_Digit(0),Tree_AS_Event(0),TRD_ST_Tracklet(0),TRD_ST_TPC_Track(0),TRD_ST_Event(0),Tree_TRD_ST_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
fListOfHistos(0x0),fTree(0x0), fPIDResponse(0), EsdTrackCuts(0),aliHelix(),TV3_SVD_tracklet_offset(),TV3_SVD_tracklet_dir(),
vec_self_tracklet_fit_points(),vec_ADC_val(),vec_TV3_TRD_center_offset(),vec_TV3_TRD_center(),TV3_trkl_offset(),TV3_trkl_dir()
{
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());

    // Output slot #0 id reserved by the base class for AOD
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}

//_______________________________________________________________________
TFile* Ali_make_tracklets_from_digits::OpenDigitsFile(TString inputfile,
				   TString digfile,
				   TString opt)
{
    // we should check if we are reading ESDs or AODs - for now, only
    // ESDs are supported

    cout << "" << endl;
    cout << "In OpenDigitsFile" << endl;
    //cout << "Digits file name: " << digfile.Data() << endl;

    if(digfile == "")
    {
	cout << "WARNING: No TRD digits file available" << endl;
	return NULL;
    }

    // TGrid::Connect("alien")
    // construct the name of the digits file from the input file
    inputfile.ReplaceAll("AliESDs.root", digfile);
    //TString inputfile_LF = "alien:///alice/data/2016/LHC16q/000265525/pass1_CENT_wSDD/16000265525037.6203/TRD.FltDigits.root";
    TString inputfile_LF = inputfile;

    // open the file
    AliInfo( "opening digits file " + inputfile_LF + " with option \"" + opt + "\"");

    cout << "inputfile: " << inputfile_LF.Data() << endl;
    //TFile* dfile = new TFile(inputfile_LF, opt);
    //if(dfile) delete dfile;
    dfile = TFile::Open(inputfile_LF);
    cout << "After TRD digits file" << endl;
    cout << "" << endl;

    if(!dfile)
    {
	AliWarning("digits file '" + inputfile + "' cannot be opened");
    }

    return dfile;
}


//_______________________________________________________________________
Bool_t Ali_make_tracklets_from_digits::UserNotify()
{
    cout << "" << endl;
    cout << "In UserNotify" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;


    if(!EsdTrackCuts) EsdTrackCuts = new AliESDtrackCuts();
    if(!fGeo) fGeo = new AliTRDgeometry;
    if(!fDigMan)
    {
        fDigMan = new AliTRDdigitsManager;
        fDigMan->CreateArrays();
    }
    if(!h_v_fit_vs_det)         h_v_fit_vs_det         = new TH1D("h_v_fit_vs_det","h_v_fit_vs_det",540,0,540);
    if(!h_LA_factor_fit_vs_det) h_LA_factor_fit_vs_det = new TH1D("h_LA_factor_fit_vs_det","h_LA_factor_fit_vs_det",540,0,540);



    cout << "Connected to GRID" << endl;

    fParam = AliTRDCommonParam::Instance();

    delete fDigitsInputFile;
    delete fDigitsOutputFile;

    cout << "Digits file pointers deleted" << endl;

    cout << "All pointers deleted" << endl;

    //AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>
    //    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    if ( ! esdH ) return kFALSE;
    if ( ! esdH->GetTree() ) return kFALSE;
    if ( ! esdH->GetTree()->GetCurrentFile() ) return kFALSE;


    //-----------------------------------
    TList* list = esdH->GetUserInfo();
    //list->Print();
    //list->Dump();

    //cout << "List: " << endl;
    //list->ls();
    TList* list_cdblist = (TList*)list->FindObject("cdbList");  // list contains paths of calibration data
    Int_t list_size = list_cdblist->Capacity();
    //cout << "cdblist, with size: " << list_size <<  endl;
    //list_cdblist->ls();
    TObjString* obj_string = (TObjString*)list_cdblist->At(0);
    TString string = obj_string->GetString();
    //cout << "String: " << string.Data() << endl;
    Int_t index_A = string.Index("[");
    string.Remove(0,index_A+1);
    //cout << "String: " << string.Data() << endl;
    Int_t index_B = string.Index(",");
    string.Remove(index_B,string.Sizeof());
    Int_t run_number_from_list = string.Atoi();
    cout << "String: " << string.Data() << ", run_number_from_list: " << run_number_from_list << endl;
    //-----------------------------------


    AliCDBManager* CDBman = AliCDBManager::Instance();
    if (!CDBman->IsDefaultStorageSet()) {
      if (fLocalMode) {
        cout << "Running in local mode" << endl;
        CDBman->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB");
      } else {
        cout << "Open connection to GRID" << endl;
        TGrid::Connect("alien");
        CDBman->SetDefaultStorage("raw://");
      }
      CDBman->SetRun(run_number_from_list);
    }


    // AliCDBEntry->GetObject()->IsA()->GetName()
    //-----------------------------------
    // Pad noise
    cout << "Open pad noise calibration file from database" << endl;
    AliCDBEntry *entryB = CDBman->Get("TRD/Calib/PadNoise",run_number_from_list); // new
    PadNoise = (AliTRDCalPad*)entryB->GetObject();
    cout << "Calibration data opened" << endl;
    //-----------------------------------


    //-----------------------------------
    // ChamberVdrift
    cout << "Open ChamberVdrift calibration file from database" << endl;
    AliCDBEntry *entryC = CDBman->Get("TRD/Calib/ChamberVdrift",run_number_from_list); // new
    ChamberVdrift = (AliTRDCalDet*)entryC->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", chamber vdrift: " << ChamberVdrift->GetValue(i_det) << endl;
    //}
    //-----------------------------------


    //-----------------------------------
    // ChamberT0
    cout << "Open ChamberT0 calibration file from database" << endl;
    AliCDBEntry *entryC1 = CDBman->Get("TRD/Calib/ChamberT0",run_number_from_list);
    ChamberT0 = (AliTRDCalDet*)entryC1->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", chamber t0: " << ChamberT0->GetValue(i_det) << endl; // in the order of "-1.36613"
    //}
    //-----------------------------------


    //-----------------------------------
    // LocalT0
    //cout << "Open LocalT0 calibration file from database" << endl;
    //AliCDBEntry *entryC2 = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/LocalT0",run_number_from_list);
    //LocalT0_pad = (AliTRDCalPad*)entryC2->GetObject();
    //cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    Int_t col = 2;
    //    Int_t row = 2;
    //    cout << "i_det: " << i_det << ", local t0: " << LocalT0_pad->GetCalROC(i_det)->GetValue(col,row) << endl; // returns 0
    //}
    //-----------------------------------


#if 0
    //-----------------------------------
    // LocalVdrift
    // Values are all 0
    cout << "Open LocalVdrift calibration file from database" << endl;
    AliCDBEntry *entryD = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/LocalVdrift",run_number_from_list);
    AliTRDCalDet *LocalVdrift = (AliTRDCalDet*)entryD->GetObject();
    cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", local vdrift: " << LocalVdrift->GetValue(i_det) << endl;
    //}
    //-----------------------------------
#endif


    //-----------------------------------
    // ChamberExB
    // Values are all 0
    //cout << "Open ChamberExB calibration file from database" << endl;
    //AliCDBEntry *entryE = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/ChamberExB",run_number_from_list);
    //ChamberExB = (AliTRDCalDet*)entryE->GetObject();
    //cout << "Calibration data opened" << endl;
    //for(Int_t i_det = 0; i_det < 540; i_det++)
    //{
    //    cout << "i_det: " << i_det << ", ExB: " << ChamberExB->GetValue(i_det) << endl;
    //}
    //-----------------------------------



    //----------------------------------------------------------------------
    // Krypton calibration -> pad gain factors
    //cout << "Open Krypto calibration file from database" << endl;
    //AliCDBEntry *entryF = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/Krypton_2015-02",run_number_from_list);
    //KryptoGain = (AliTRDCalOnlineGainTable*)entryF->GetObject();
    //cout << "Calibration data opened" << endl;
    //Float_t GainFactor = KryptoGain ->GetGainCorrectionFactor(i_det,i_row,i_column);
    //----------------------------------------------------------------------



    //----------------------------------------------------------------------
    // Chamber gain factors
    //cout << "Open chamber gain file from database" << endl;
    //AliCDBEntry *entryG = AliCDBManager::Instance()->GetStorage(pathdatabase)->Get("TRD/Calib/ChamberGainFactor",run_number_from_list);
    //chambergain = (AliTRDCalDet*)entryG->GetObject();
    //Float_t valuegainiteration = chambergain->GetValue(det);
    //----------------------------------------------------------------------



    //TObjString* ts = (TObjString*)list->FindObject("TRD/Calib/ChamberGainFactor");
    //TString string = ts->GetString();
    //cout << "string: " << string.Data() << endl;

    TString fname = esdH->GetTree()->GetCurrentFile()->GetName();
    TString Tree_name = esdH->GetTree()->GetName();
    FileStat_t file_stat;
    Int_t PathInfo = gSystem->GetPathInfo(fname.Data(),file_stat);
    cout << "PathInfo: " << PathInfo << ", fname: " << fname << endl;
    //TFile* file = TFile::Open(fname.Data());
    //cout << "Zombie: " << file->IsZombie() << ", header size: " << file->Sizeof() << ", FileBytesRead: " << file->GetFileBytesRead() << endl;

    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!inputHandler)
    {
	printf("WARNING: Inputhandler not available \n");
    }
    else
    {
	printf("Inputhandler available \n");

	fPIDResponse = inputHandler->GetPIDResponse();

        cout << "Got PID response" << endl;
    }

    //AliAnalysisManager *manB = AliAnalysisManager::GetAnalysisManager();
    //Int_t run_id = manB->GetRunFromPath();
    //cout << "fname: " << fname << ", tree name: " << Tree_name << ", run_id: " << run_id << endl;


    fEventNoInFile = -1;
    //N_good_events  = 0;


    fDigitsInputFile = OpenDigitsFile(fname,fDigitsInputFileName,""); // <-


#if 0
    // Connect the friends
    //fESD ->SetESDfriend(esdFr);

    esdFriendTreeFName = fname;
    TString basename = gSystem->BaseName(esdFriendTreeFName);
    int index = basename.Index("#")+1;
    basename.Remove(index);
    basename += "AliESDfriends.root";
    TString dirname = gSystem->DirName(esdFriendTreeFName);
    dirname += "/";
    esdFriendTreeFName = dirname + basename;
    cout << "Friend name: " << esdFriendTreeFName.Data() << endl;
#endif




    cout << "Add EsdTrackCuts" << endl;
    if(EsdTrackCuts) cout << "EsdTrackCuts exists" << endl;
    EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
    EsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.52);
    EsdTrackCuts->AliESDtrackCuts::SetMinNClustersTPC(50); // 60, Automatically requires TPC refitted tracks?
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXY(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(10.0);
    EsdTrackCuts->AliESDtrackCuts::SetPtRange(0.15,200.0); // 0.15, 200.0
    EsdTrackCuts->AliESDtrackCuts::SetEtaRange(-1.0,1.0); // 0.85

    // create the digits manager
    cout << "" << endl;
    cout << "________________________________________________________________________" << endl;
    cout << "Created AliTRDdigitsManager" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;
    cout << "" << endl;

    // create a TRD geometry, needed for matching digits to tracks
    if(!fGeo)
    {
	AliFatal("cannot create geometry ");
    }

    //if(fDigMan) delete fDigMan;

    //for(Int_t i_det = 0; i_det < 5; i_det++)
    //{
    //    Int_t N_columns   = fDigMan->GetDigits(i_det)->GetNcol();
    //    cout << "i_det: " << i_det << ", N_columns: " << N_columns << endl;
    //}


    //---------------------------------------
    // Load TRD geometry, created by Create_TRD_geometry_files.cc
    TFile* file_TRD_geom = fLocalMode ? TFile::Open("Data/TRD_geometry_full.root") : TFile::Open("alien::///alice/cern.ch/user/a/aschmah/Data/TRD_geometry_full.root");

    vec_TV3_TRD_center.resize(540);
    vec_TV3_TRD_center_offset.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_TRD_center[i_det].resize(3);
    }

    vector< vector<TH1D*> > vec_TH1D_TV3_TRD_center;
    vec_TH1D_TV3_TRD_center.resize(3); // x,y,z axes
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TV3_TRD_center[i_xyz].resize(3); // vector direction
        for(Int_t i_dir_component = 0; i_dir_component < 3; i_dir_component++)
        {
            HistName = "vec_TH1D_TV3_TRD_center_";
            HistName += i_xyz;
            HistName += "_V";
            HistName += i_dir_component;
            vec_TH1D_TV3_TRD_center[i_xyz][i_dir_component] = (TH1D*)file_TRD_geom->Get(HistName.Data());

            for(Int_t i_det = 0; i_det < 540; i_det++)
            {
                Double_t value = vec_TH1D_TV3_TRD_center[i_xyz][i_dir_component] ->GetBinContent(i_det+1);
                vec_TV3_TRD_center[i_det][i_xyz][i_dir_component] = value;
                //printf("i_xyz: %d, i_dir_component: %d, i_det: %d, value: %4.3f \n",i_xyz,i_dir_component,i_det,value);
                //vec_TV3_TRD_center[i_det][i_xyz].Print();
            }
        }
    }

    vector<TH1D*> vec_TH1D_TV3_TRD_center_offset;
    vec_TH1D_TV3_TRD_center_offset.resize(3); // x,y,z axes
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        HistName = "vec_TH1D_TV3_TRD_center_offset_";
        HistName += i_xyz;
        vec_TH1D_TV3_TRD_center_offset[i_xyz] = (TH1D*)file_TRD_geom->Get(HistName.Data());

        for(Int_t i_det = 0; i_det < 540; i_det++)
        {
            vec_TV3_TRD_center_offset[i_det][i_xyz] = vec_TH1D_TV3_TRD_center_offset[i_xyz] ->GetBinContent(i_det+1);
        }
    }

    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
    vec_TH1D_TRD_geometry.resize(3); // x,y,z
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TRD_geometry[i_xyz].resize(8); // 8 vertices
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            HistName = "vec_TH1D_TRD_geometry_xyz_";
            HistName += i_xyz;
            HistName += "_V";
            HistName += i_vertex;
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] = (TH1D*)file_TRD_geom->Get(HistName.Data());
        }
    }
    //---------------------------------------

    if(fname.Contains("/home/"))
    {
        cout << "Load local alignment file" << endl;
        TRD_alignment_file = TFile::Open("/home/ceres/schmah/ALICE/Database/TRD_Align_2016.root");
	cout << "Local alignment file loaded" << endl;
    }
    else
    {
        cout << "Load alignment file" << endl;
        TRD_alignment_file = fLocalMode ? TFile::Open("/cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB/TRD/Align/Data/Run0_999999999_v1_s0.root") : TFile::Open("alien:///alice/data/2016/OCDB/TRD/Align/Data/Run0_999999999_v1_s0.root");
        //TRD_alignment_file = TFile::Open("/alice/data/2016/OCDB/TRD/Align/Data/Run0_999999999_v1_s0.root");
        cout << "Alignment file from database loaded" << endl;
    }

    cout << "Open calibration file" << endl;
    //TRD_calibration_file_AA = TFile::Open("alien::///alice/cern.ch/user/a/aschmah/Data/TRD_Calib_vDfit_and_LAfit_3456.root");
    TRD_calibration_file_AA = fLocalMode ? TFile::Open("Data/TRD_Calib_vDfit_and_LAfit_23_11_2020.root ") : TFile::Open("alien::///alice/cern.ch/user/a/aschmah/Data/TRD_Calib_vDfit_and_LAfit_23_11_2020.root ");
    cout << "Calibration file opened" << endl;
    tg_v_fit_vs_det         = (TGraph*)TRD_calibration_file_AA ->Get("tg_v_fit_vs_det");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        h_v_fit_vs_det ->SetBinContent(i_det+1,1.05);
    }
    for(Int_t i_point = 0; i_point < tg_v_fit_vs_det->GetN(); i_point++)
    {
        Double_t det, vD;
        tg_v_fit_vs_det->GetPoint(i_point,det,vD);
        h_v_fit_vs_det ->SetBinContent(det+1,vD);
    }
    tg_LA_factor_fit_vs_det = (TGraph*)TRD_calibration_file_AA ->Get("tg_LA_factor_fit_vs_det");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
         h_LA_factor_fit_vs_det->SetBinContent(i_det+1,-0.14);
    }
    for(Int_t i_point = 0; i_point < tg_LA_factor_fit_vs_det->GetN(); i_point++)
    {
        Double_t det, LA;
        tg_LA_factor_fit_vs_det->GetPoint(i_point,det,LA);
        h_LA_factor_fit_vs_det ->SetBinContent(det+1,LA);
    }
    cout << "Calibration file opened" << endl;

    AliCDBEntry* align_cdb = (AliCDBEntry*)TRD_alignment_file->Get("AliCDBEntry");
    TClonesArray* TRD_align_array = (TClonesArray*)align_cdb->GetObject();
    // Alignment is stored for every INSTALLED chamber. Array does NOT got to 540!
    // First do the alignment per sector, then per chamber
    for(Int_t i_entry = 0; i_entry < TRD_align_array->GetEntries(); i_entry++)
    {
	AliAlignObjParams* Align_test = (AliAlignObjParams*)TRD_align_array->At(i_entry);
	char* name = (char*)Align_test->GetSymName();
	TString sname(name);
        Int_t length = sname.Length();
	TString sector_name = sname(6,2);
	TString stack_name  = sname(11,1);
	TString layer_name  = sname(15,1);
	Int_t sector = sector_name.Atoi();
	Int_t stack  = stack_name.Atoi();
	Int_t layer  = layer_name.Atoi();
        Int_t detector = layer + stack * 6 + sector * 6 * 5;

	TGeoHMatrix TM_TRD_align;
        Align_test->GetMatrix(TM_TRD_align);

	//printf("i_entry: %d, name: %s, length: %d, sector: %s, stack: %s, layer: %s, detector: %d \n",i_entry,name,sname.Length(),sector_name.Data(),stack_name.Data(),layer_name.Data(),detector);
        //TM_TRD_align.Print();
	if(length > 8)
	{
	    Align_test->GetMatrix(TM_TRD_rotation_det[detector]);
	    Double_t* rot_matrix = TM_TRD_rotation_det[detector].GetRotationMatrix();
	    //for(Int_t i = 0; i < 9; i++)
	    //{
	    //    cout << "i: " << i << ", rot_matrix: " << rot_matrix[i] << endl;
	    //}
	    Double_t TRD_translation[3];
	    Align_test->GetTranslation(TRD_translation);
	    TV3_TRD_translation[detector].SetXYZ(TRD_translation[0],TRD_translation[1],TRD_translation[2]);
	    //cout << "i_entry: " << i_entry << ", trans = {" << TV3_TRD_translation[i_entry].X() << ", " << TV3_TRD_translation[i_entry].Y() << ", " << TV3_TRD_translation[i_entry].Z() << "}" << endl;
	    //TM_TRD_rotation_det[i_entry].Print();
	}
	else
	{
	    Align_test->GetMatrix(TM_TRD_rotation_sector[sector]);
	}
    }
    cout << "TRD_align_array entries: " << TRD_align_array->GetEntries() << endl;



    cout << "End of UserNotify" << endl;
    return kTRUE;
}


//________________________________________________________________________
void Ali_make_tracklets_from_digits::UserCreateOutputObjects()
{
    cout << "" << endl;
    cout << "In UserCreateOutputObjects" << endl;
    cout << "fDigitsInputFileName: " << fDigitsInputFileName.Data() << endl;


    OpenFile(1);
    cout << "File opened" << endl;

    fListOfHistos = new TList();
    fListOfHistos ->SetOwner();


    OpenFile(2);
    cout << "File opened" << endl;

    AS_Event       = new Ali_AS_Event();
    AS_Track       = new Ali_AS_Track();
    //AS_Tracklet    = new Ali_AS_Tracklet();
    //AS_offline_Tracklet    = new Ali_AS_offline_Tracklet();
    AS_Digit       = new Ali_AS_TRD_digit();
    Tree_AS_Event  = NULL;
    Tree_AS_Event  = new TTree("Tree_AS_Event" , "AS_Events" );
    Tree_AS_Event  ->Branch("Tree_AS_Event_branch"  , "AS_Event", AS_Event );


    TRD_ST_Tracklet   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event      = new Ali_TRD_ST_Event();

    Tree_TRD_ST_Event  = NULL;
    Tree_TRD_ST_Event  = new TTree("Tree_TRD_ST_Event" , "TRD_ST_Events" );
    Tree_TRD_ST_Event  ->Branch("Tree_TRD_ST_Event_branch"  , "TRD_ST_Event", TRD_ST_Event );


    PostData(1,fListOfHistos);
    //PostData(2,Tree_AS_Event);
    PostData(2,Tree_TRD_ST_Event);


    cout << "PostData called" << endl;

}



//________________________________________________________________________
Bool_t Ali_make_tracklets_from_digits::NextEvent(Bool_t preload)
{
    fEventNoInFile++;
    //cout << "fEventNoInFile: " << fEventNoInFile << endl;
    fDigitsLoadedFlag = kFALSE;


    if(preload)
    {
	//cout << "Preload"  << endl;
	return ReadDigits();
    }
    else
    {
	//cout << "No preload"  << endl;
	return kTRUE;
    }
}

//________________________________________________________________________
void Ali_make_tracklets_from_digits::UserExec(Option_t *)
{
    //cout << "" << endl;
    //cout << "Analysis started" << endl;
    //cout << "----------------------------------------------------------------------------------------------------------------------------------" << endl;

    Int_t flag_calibrated = 1; // 0 = standard fixed precalibration used, 1 = use pre calibration from root input file
    //-----------------------------------------------------------------
    // IMPORTANT: call NextEvent() for book-keeping
    NextEvent();
    //if(fEventNoInFile > 50) return;
    //-----------------------------------------------------------------



    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    // prepare event data structures
    AliESDEvent* fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD)
    {
	printf("ERROR: fESD not available\n");
	return;
    }

    //-----------------------------------------------------------------
    // Check if TRD digits (raw data) are available for this ESD event
    if(!ReadDigits()) return;
    //-----------------------------------------------------------------

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(man)
    {
        //Int_t run_id = man->GetRunFromPath(); // doesn't work
	//cout << "Got AliAnalysisManager, run_id: " << run_id << endl;
	AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
	if(inputHandler)
	{
	    //cout << "Got AliInputEventHandler" << endl;
	    fPIDResponse = inputHandler->GetPIDResponse();
	}
    }
    //cout << "cent: " << fPIDResponse->GetCurrentCentrality() << endl;

    Int_t 	       eventNumber      = fESD ->GetEventNumberInFile();

    Int_t          N_tracks         = fESD ->GetNumberOfTracks();
    Int_t          N_TRD_tracks     = fESD ->GetNumberOfTrdTracks();
    Int_t          N_TRD_tracklets  = fESD ->GetNumberOfTrdTracklets(); // online
    Float_t        magF             = fESD ->GetMagneticField();
    const AliESDVertex* PrimVertex  = fESD ->GetPrimaryVertex();
    Int_t          RunNum           = fESD ->GetRunNumber();
    Double_t       T0zVertex        = fESD ->GetT0zVertex();
    AliCentrality* Centrality       = fESD ->GetCentrality();
    Double_t       MeanBeamIntAA    = fESD ->GetESDRun()->GetMeanIntensity(0,0);

    //printf("RunNum: %d \n",RunNum);

    Double_t Sign_magnetic_field = (magF/fabs(magF));
    //cout << "Trigger: " <<  fESD->GetFiredTriggerClasses() << endl;

    Int_t digit_counter = 0;

    //------------------------------------------
    // For space points (digits) tree
    // Fill event information
    AS_Event ->clearTrackList();
    AS_Event ->clearTrackletList();
    AS_Event ->clearTRD_digit_list();
    AS_Event ->setTriggerWord(fESD->GetFiredTriggerClasses());
    AS_Event ->setx(PrimVertex->GetX());
    AS_Event ->sety(PrimVertex->GetY());
    AS_Event ->setz(PrimVertex->GetZ());
    AS_Event ->setid(RunNum);
    AS_Event ->setN_tracks(N_tracks);
    AS_Event ->setN_TRD_tracklets(N_TRD_tracklets); // online
    AS_Event ->setBeamIntAA(MeanBeamIntAA);
    AS_Event ->setT0zVertex(T0zVertex);
    AS_Event ->setEventNumber(eventNumber);

    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    if(MultSelection)
    {
	// V0MEq, V0AEq, V0CEq, SPDTracklets

	AS_Event ->setcent_class_ZNA(MultSelection->GetMultiplicityPercentile("ZNA"));
	AS_Event ->setcent_class_ZNC(MultSelection->GetMultiplicityPercentile("ZNC"));
	AS_Event ->setcent_class_V0A(MultSelection->GetMultiplicityPercentile("V0A"));
	AS_Event ->setcent_class_V0C(MultSelection->GetMultiplicityPercentile("V0C"));
	AS_Event ->setcent_class_V0M(MultSelection->GetMultiplicityPercentile("V0M"));
	AS_Event ->setcent_class_CL0(MultSelection->GetMultiplicityPercentile("CL0"));
	AS_Event ->setcent_class_CL1(MultSelection->GetMultiplicityPercentile("CL1"));
	AS_Event ->setcent_class_SPD(MultSelection->GetMultiplicityPercentile("SPDTracklets"));
	AS_Event ->setcent_class_V0MEq(MultSelection->GetMultiplicityPercentile("V0MEq"));
	AS_Event ->setcent_class_V0AEq(MultSelection->GetMultiplicityPercentile("V0AEq"));
	AS_Event ->setcent_class_V0CEq(MultSelection->GetMultiplicityPercentile("V0CEq"));
    }
    //------------------------------------------


    //------------------------------------------
    // Fill event information for tracklets tree
    TRD_ST_Event ->clearTrackList();
    TRD_ST_Event ->clearTrackletList();

    TRD_ST_Event ->setTriggerWord( AS_Event->getTriggerWord() );
    TRD_ST_Event ->setx( AS_Event->getx() );
    TRD_ST_Event ->sety( AS_Event->gety() );
    TRD_ST_Event ->setz( AS_Event->getz() );
    TRD_ST_Event ->setid( AS_Event->getid() );
    TRD_ST_Event ->setN_tracks( AS_Event->getN_tracks() );
    TRD_ST_Event ->setN_TRD_tracklets( AS_Event->getN_TRD_tracklets() );
    TRD_ST_Event ->setBeamIntAA( AS_Event->getBeamIntAA() );
    TRD_ST_Event ->setT0zVertex( AS_Event->getT0zVertex() );
    TRD_ST_Event ->setcent_class_ZNA( AS_Event->getcent_class_ZNA() );
    TRD_ST_Event ->setcent_class_ZNC( AS_Event->getcent_class_ZNC() );
    TRD_ST_Event ->setcent_class_V0A( AS_Event->getcent_class_V0A() );
    TRD_ST_Event ->setcent_class_V0C( AS_Event->getcent_class_V0C() );
    TRD_ST_Event ->setcent_class_V0M( AS_Event->getcent_class_V0M() );
    TRD_ST_Event ->setcent_class_CL0( AS_Event->getcent_class_CL0() );
    TRD_ST_Event ->setcent_class_CL1( AS_Event->getcent_class_CL1() );
    TRD_ST_Event ->setcent_class_SPD( AS_Event->getcent_class_SPD() );
    TRD_ST_Event ->setcent_class_V0MEq( AS_Event->getcent_class_V0MEq() );
    TRD_ST_Event ->setcent_class_V0AEq( AS_Event->getcent_class_V0AEq() );
    TRD_ST_Event ->setcent_class_V0CEq( AS_Event->getcent_class_V0CEq() );

    //printf("Event information filled \n");
    //------------------------------------------



    //-----------------------------------------------------------------
    //cout << "" << endl;
    //cout << "" << endl;
    //cout << "----------------------------------------------------------------------------------------" << endl;
    //cout << "Event number: " << fEventNoInFile << ", event number with TRD digits: " << N_good_events << endl;
    //printf("Event number: %d, N_tracks: %d, cent(V0M): %f , cent(CL0): %f \n",fEventNoInFile,N_tracks,MultSelection->GetMultiplicityPercentile("SPDTracklets"));
    //-----------------------------------------------------------------


    //-----------------------------------------------------------------
    Double_t TRD_time_per_bin        = 0.1;   // 100 ns = 0.1 mus
    Double_t TRD_drift_velocity      = 1.56;  // 1.56 cm/mus, this value is too high for p+Pb, database values are now used
    Double_t TPC_TRD_matching_window = 10.0;  // Matching window between TRD digits (pad position) and TPC track in cm
    Double_t TPC_min_radius_plot     = 114.0/2.0; // Minimum radius for TPC track to be plotted
    Double_t TPC_max_radius_plot     = 368.0; // Maximum radius for TPC track to be plotted
    Double_t TPC_radius_scan         = 368.0 - 30.0; // Radius for TPC track to be used for first match with TRD digits

    if(N_tracks == 0)
    {
	// Skip empty event
	return;
    }

    //printf("There are %d tracks in this event\n", N_tracks);
    //printf("There are %d TRD tracks in this event\n", N_TRD_tracks);
    //printf("There are %d TRD tracklets in this event\n", N_TRD_tracklets);


    TClonesArray* TC_tofHits = fESD->GetESDTOFHits();
    //printf(" \n");
    //printf("   ------------> TC_tofHits: \n");
    // AliCDBEntry->GetObject()->IsA()->GetName()   only root 5
    Int_t N_entries_TOF = TC_tofHits ->GetEntries();
    Double_t TOF_R = ((AliESDTOFHit*)TC_tofHits ->First())->GetR();
    //cout << "N_tracks: " << N_tracks <<  ", TC_tofHits: " << TC_tofHits << ", IsA: " << TC_tofHits ->First()->IsA()->GetName() << ", TOF_R: " << TOF_R << ", N_entries_TOF: " << N_entries_TOF << endl;
    for(Int_t i_tof = 0; i_tof < N_entries_TOF; i_tof++)
    {
        Double_t TOF_Z       = ((AliESDTOFHit*)TC_tofHits->At(i_tof))->GetZ();
        Int_t    TOF_channel = ((AliESDTOFHit*)TC_tofHits->At(i_tof))->GetTOFchannel();
        Double_t TOF_time    = ((AliESDTOFHit*)TC_tofHits->At(i_tof))->GetTime();
        //printf("i_tof: %d, TOF_Z: %4.3f, TOF_channel: %d, TOF_time: %d \n",i_tof,TOF_Z,TOF_channel,TOF_time);
    }
    //printf(" \n");


    //-----------------------------------------------------------------
    // TRD online tracklet loop
    for(Int_t iTracklet = 0; iTracklet < N_TRD_tracklets; iTracklet++)
    {
	AliESDTrdTracklet* tracklet = fESD->GetTrdTracklet(iTracklet);
        if(!tracklet)
        {
	    printf("ERROR: Could not receive tracklet %d\n", iTracklet);
	    continue;
        }
        Int_t det_tracklet   = tracklet ->GetDetector();
        Int_t BinY_tracklet  = tracklet ->GetBinY();
        Int_t BinDy_tracklet = tracklet ->GetBinDy();
        Int_t BinZ_tracklet  = tracklet ->GetBinZ(); // row

        // Table 10 in Jochen Kleins thesis:
        Double_t y_local  = ((Double_t)BinY_tracklet)*160.0*1E-04; // in cm
        Double_t dy_local = ((Double_t)BinDy_tracklet)*140.0*1E-04; // in cm over 3 cm
        //printf("iTracklet: %d, det_tracklet: %d \n",iTracklet,det_tracklet);

        AliTRDpadPlane* padplane = fGeo ->GetPadPlane(det_tracklet);
        Int_t i_layer  = fGeo ->GetLayer(det_tracklet);
        Int_t i_sector = fGeo ->GetSector(det_tracklet);
        Float_t TRD_time0 = fGeo ->GetTime0(i_layer); // in cm
        Double_t TRD_row_pos  = padplane ->GetRowPos(BinZ_tracklet);       // fPadRow[row] + fPadRowSMOffset;
        Double_t TRD_row_size = padplane ->GetRowSize(BinZ_tracklet);

        //------------------------
        // Calculate tracklet offset space point
        Double_t loc_tracklet[3] = {TRD_time0 - 0.0,y_local,TRD_row_pos - TRD_row_size/2.0};
        Double_t             glb_tracklet[3]              = {0.0,0.0,0.0};
        Double_t             glb_tracklet_align_sec[3]    = {0.0,0.0,0.0};
        Double_t             glb_tracklet_align[3]        = {0.0,0.0,0.0}; // offset space point for tracklet
        fGeo ->RotateBack(det_tracklet,loc_tracklet,glb_tracklet);
        TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_tracklet,glb_tracklet_align_sec);
        TM_TRD_rotation_det[det_tracklet].LocalToMaster(glb_tracklet_align_sec,glb_tracklet_align);
        //------------------------


        //------------------------
        // Calculate tracklet direction vector
        Double_t loc_tracklet_vec[3] = {TRD_time0 - 3.0,y_local - dy_local,TRD_row_pos - TRD_row_size/2.0};
        Double_t             glb_tracklet_vec[3]              = {0.0,0.0,0.0};
        Double_t             glb_tracklet_align_sec_vec[3]    = {0.0,0.0,0.0};
        Double_t             glb_tracklet_align_vec[3]        = {0.0,0.0,0.0}; // offset space point for tracklet
        fGeo ->RotateBack(det_tracklet,loc_tracklet_vec,glb_tracklet_vec);
        TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_tracklet_vec,glb_tracklet_align_sec_vec);
        TM_TRD_rotation_det[det_tracklet].LocalToMaster(glb_tracklet_align_sec_vec,glb_tracklet_align_vec);

        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            glb_tracklet_align_vec[i_xyz] -= glb_tracklet_align[i_xyz];
        }
        //------------------------

        TVector3 TV3_offset(glb_tracklet_align[0],glb_tracklet_align[1],glb_tracklet_align[2]);
        TVector3 TV3_dir(glb_tracklet_align_vec[0],glb_tracklet_align_vec[1],glb_tracklet_align_vec[2]);


        //printf("det: %d, sector: %d \n",det_tracklet,i_sector);
        //printf("Local  tracklet pos: {%4.3f, %4.3f, %4.3f} \n",TRD_time0 - 0.0,y_local,TRD_row_pos - TRD_row_size/2.0);
        //printf("Online tracklet pos: {%4.3f, %4.3f, %4.3f} \n",glb_tracklet_align[0],glb_tracklet_align[1],glb_tracklet_align[2]);
        //printf(" \n");

        AS_Tracklet  = AS_Event ->createTracklet(); // online tracklet
        AS_Tracklet  ->set_detector(det_tracklet);
        AS_Tracklet  ->set_TV3_offset(TV3_offset);
        AS_Tracklet  ->set_TV3_dir(TV3_dir);
        AS_Tracklet  ->set_online_dy(dy_local);
    }
    //-----------------------------------------------------------------



    Int_t N_good_tracks = 0;


    //-----------------------------------------------------------------
    // Loop over all TRD channels
    std::vector<TVector3> TV3_TRD_hits;
    std::vector<TVector3> TV3_TRD_hits_uncalib;
    std::vector<TVector3> TV3_TRD_hits_det_angle;
    std::vector<Double_t> vec_TRD_hits_ADC_value;
    std::vector<Double_t> vec_TRD_hits_time;
    std::vector< std::vector<Int_t> >  vec_TRD_det_lay_row_col;

    std::vector<TVector3> TV3_TRD_hits_middle;
    std::vector<TVector3> TV3_TRD_hits_det_angle_middle;
    std::vector< std::vector<Double_t> > vec_TRD_hits_ADC_values_time;
    std::vector< std::vector<TVector3> > vec_TRD_hits_points_time;
    std::vector< std::vector<TVector3> > vec_TRD_hits_points_time_uncalib;
    std::vector< std::vector<Int_t> >    vec_TRD_hits_det_lay_row_col;

    vec_TRD_det_lay_row_col.resize(4);

    Int_t max_N_rows    = 0;
    Int_t max_N_columns = 0;

    Double_t sum_full_ADC_digit_det[540];
    memset(sum_full_ADC_digit_det, 0, sizeof(sum_full_ADC_digit_det)); // for automatically-allocated arrays

    //cout << "Loop over all TRD detectors" << endl;
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
	Int_t N_columns   = fDigMan->GetDigits(i_det)->GetNcol();
	Int_t N_rows      = fDigMan->GetDigits(i_det)->GetNrow();
        Int_t N_times     = fDigMan->GetDigits(i_det)->GetNtime();
        //printf("N_times: %d \n",N_times);
	Int_t N_detectors = fDigMan->GetDigits(i_det)->GetNdet();

	Int_t i_sector = fGeo->GetSector(i_det);
	Int_t i_stack  = fGeo->GetStack(i_det);
        Int_t i_layer  = fGeo->GetLayer(i_det);

        //for uncalibrated digits: fixed vD and LA
        Double_t vD_calib = 1.546; // 1.546, 1.2
        Double_t LA_calib = -0.16133; // -0.16133

        //for closure test: calib parameters from our file
        if (flag_calibrated)
        {
            vD_calib = h_v_fit_vs_det         ->GetBinContent(i_det + 1);
            LA_calib = h_LA_factor_fit_vs_det ->GetBinContent(i_det + 1);
        }

        //printf("i_det(start from 0): %d, vD_calib: %4.5f, LA_calib: %4.5f \n",i_det,vD_calib,LA_calib);

	//printf("i_det: %d, N_columns: %d, N_rows: %d, N_times: %d, i_sector: %d, i_stack: %d, i_layer: %d \n",i_det,N_columns,N_rows,N_times,i_sector,i_stack,i_layer);

	if(N_columns > max_N_columns) max_N_columns = N_columns;
	if(N_rows    > max_N_rows)    max_N_rows    = N_rows;

	for(Int_t i_column = 0; i_column < N_columns; i_column++)
	{
	    for(Int_t i_row = 0; i_row < N_rows; i_row++)
            {
                Int_t x_TRD = i_column + i_sector*144;
                Int_t y_TRD = i_row    + i_stack*16 + i_layer*16*5;

		Double_t ADC_amplitude_sum_times = 0.0;
		std::vector<Double_t> vec_ADC_time_bins;
		vec_ADC_time_bins.resize(N_times);
                std::vector<TVector3> vec_points_time_bins;
                vec_points_time_bins.resize(N_times);
                std::vector<TVector3> vec_points_time_bins_uncalib;
                vec_points_time_bins_uncalib.resize(N_times);
		//cout << "i_column: " << i_column << ", i_row: " << i_row << ", N_times: " << N_times << endl;
                Double_t arr[30] = {0.0};

		for(Int_t i_time = 0; i_time < N_times; i_time++)
		{
		    Double_t ADC_amplitude = fDigMan->GetDigitAmp(i_row,i_column,i_time,i_det);
		    Int_t fBaseline = fDigMan->GetDigitsParam()->GetADCbaseline(i_det); // ADC baseline to be subtracted from digit ADC values, constant value of 10
                    ADC_amplitude_sum_times += ADC_amplitude;


                    //printf("ADC_amplitude: %4.3f \n",ADC_amplitude);

		    if(ADC_amplitude > 0.0)
		    {
                        sum_full_ADC_digit_det[i_det] += ADC_amplitude - fBaseline;

                        //if(i_time == 0) printf("sector: %d \n",i_sector);
                        //printf("i_time: %d, ADC: %f \n",i_time,ADC_amplitude);
			arr[i_time] = ADC_amplitude;
			//cout << "i_time: " << i_time << ", fill arr: " << arr[i_time] << endl;
                        //cout << "fBaseline: " << fBaseline << ", ADC: " << ADC_amplitude << endl;

			//----------------------
			// Calculate global position for fired TRD pad
			Float_t              TRD_time0                = fGeo         ->GetTime0(i_layer); // in cm
                        AliTRDpadPlane*      padplane                 = fGeo         ->GetPadPlane(i_det);
			Double_t             TRD_col_end              = padplane     ->GetColEnd();
			Double_t             TRD_row_end              = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
			Double_t             TRD_row_start            = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
			Double_t             TRD_row_end_ROC          = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
			Double_t             TRD_col_spacing          = padplane     ->GetColSpacing();
			Double_t             TRD_row_spacing          = padplane     ->GetRowSpacing();
			Double_t             TRD_col_pos              = padplane     ->GetColPos(i_column);
			Double_t             TRD_row_pos              = padplane     ->GetRowPos(i_row);       // fPadRow[row] + fPadRowSMOffset;
			Double_t             TRD_row_pos_ROC          = padplane     ->GetRowPosROC(i_row);    // fPadRow[row]; = pad border position
		        Double_t             TRD_col_size             = padplane     ->GetColSize(i_column);
			Double_t             TRD_row_size             = padplane     ->GetRowSize(i_row);
                        Float_t              TRD_loc_Y                = TRD_col_pos + 0.5*TRD_col_size;
                        Float_t              TRD_loc_Y_uncalib        = TRD_loc_Y;
                        //Double_t             TRD_drift_time           = ((Double_t)i_time)*TRD_time_per_bin*ChamberVdrift->GetValue(i_det); // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum
                        Double_t             TRD_drift_time           = ((Double_t)i_time)*TRD_time_per_bin*vD_calib; // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum
                        Double_t             TRD_drift_time_uncalib   = ((Double_t)i_time)*TRD_time_per_bin*1.56; // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum

                        //printf("TRD_tim0: %4.3f, vdrift: %4.3f \n",TRD_time0,ChamberVdrift->GetValue(i_det));

                        //Float_t lorentz_angle_corr_y = ChamberExB->GetValue(i_det)*TRD_drift_time;
                        Float_t lorentz_angle_corr_y = -TMath::Tan(LA_calib)*TRD_drift_time;
                        TRD_loc_Y -= Sign_magnetic_field*lorentz_angle_corr_y;

                        //printf("LA_calib: %4.3f, LA: %4.3f \n",LA_calib,ChamberExB->GetValue(i_det));

			//cout << "i_time: " << i_time << ", magF: " << magF << ", x: " << TRD_drift_time << ", lorentz_angle_corr_y: " << lorentz_angle_corr_y << endl;

			//Float_t              TRD_loc_Z        = (Float_t)i_row*8.5+8.5/2.0;
			//Float_t              TRD_loc_Y        = (Float_t)i_column*0.725 + 0.725/2.0 - (144.0*0.725)/2.0;
			//Double_t             loc[3]           = {TRD_time0,TRD_loc_Y,TRD_loc_Z+TRD_row_end};

			// Determine global position of TRD hit
                        //Double_t             loc[3]           = {TRD_time0 - TRD_drift_time,TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};
                        //Double_t             loc_uncalib[3]   = {TRD_time0 - TRD_drift_time_uncalib,TRD_loc_Y_uncalib,TRD_row_pos - TRD_row_size/2.0};
                        Double_t             loc_uncalib[3]  = {TRD_time0 - TRD_drift_time,TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};
                        Double_t             loc[3]          = {TRD_time0 - TRD_drift_time_uncalib,TRD_loc_Y_uncalib,TRD_row_pos - TRD_row_size/2.0};

                        Double_t             glb[3]           = {0.0,0.0,0.0};
                        Double_t             glb_uncalib[3]   = {0.0,0.0,0.0};
                        fGeo ->RotateBack(i_det,loc,glb);
                        fGeo ->RotateBack(i_det,loc_uncalib,glb_uncalib);

			// Apply alignment     
			Double_t glb_align_sec[3];
                        Double_t glb_align[3];
                        Double_t glb_align_sec_uncalib[3];
			Double_t glb_align_uncalib[3];

                        TM_TRD_rotation_sector[i_sector].LocalToMaster(glb,glb_align_sec);
                        TM_TRD_rotation_det[i_det].LocalToMaster(glb_align_sec,glb_align);
                        TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_uncalib,glb_align_sec_uncalib);
                        TM_TRD_rotation_det[i_det].LocalToMaster(glb_align_sec_uncalib,glb_align_uncalib);

			//glb_align[0] = glb[0];
                        //glb_align[1] = glb[1];
                        //glb_align[2] = glb[2];

			// Determine global pointing vector of TRD plane (vector perpendicular to plane in global system)  // here no changes
			Double_t             loc_vec[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
			Double_t             glb_vec[3]           = {0.0,0.0,0.0};
			fGeo ->RotateBack(i_det,loc_vec,glb_vec);


			// Apply alignment  // here no changes
			Double_t glb_vec_align_sec[3];
			Double_t glb_vec_align[3];
			TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_vec,glb_vec_align_sec);
			TM_TRD_rotation_det[i_det].LocalToMaster(glb_vec_align_sec,glb_vec_align);

                        //for(Int_t i = 0; i < 3; i++)
			//{
			//    glb_align[i] = glb[i];
                        //    glb_vec_align[i] = glb_vec[i];
			//}


                        Double_t             TRD_radius           = TMath::Sqrt(glb_align[0]*glb_align[0] + glb_align[1]*glb_align[1]);
                        Double_t             TRD_radius_uncalib   = TMath::Sqrt(glb_align_uncalib[0]*glb_align_uncalib[0] + glb_align_uncalib[1]*glb_align_uncalib[1]);

#if 0
			if(i_det >= 0 && i_time == 0 && i_column == 0 && i_row >= 0)
			{
			    cout << "i_det: " << i_det << ", i_row: " << i_row << ", N_rows: " << N_rows << ", row_end: " << TRD_row_end
				<< ", row_end_ROC: " << TRD_row_end_ROC << ", row_pos: " << TRD_row_pos
				<< ", row_pos_ROC: " << TRD_row_pos_ROC
				<< ", row_size: " << TRD_row_size << endl;
			}
#endif

#if 0
			cout << "ADC_amplitude: " << ADC_amplitude << ", i_row: " << i_row
			    << ", i_column: " << i_column << ", i_time: " << i_time << ", detector: " << i_det
			    << ", col_end: " << TRD_col_end << ", row_end: " << TRD_row_end
			    << ", col_spacing: " << TRD_col_spacing << ", row_spacing: " << TRD_row_spacing
			    << ", col_pos: " << TRD_col_pos << ", row_pos: " << TRD_row_pos
			    << ", col_size: " << TRD_col_size << ", row_size: " << TRD_row_size << endl;
#endif

			if(ADC_amplitude > 0.0)
                        {
                            //printf("i_det: %d, i_column: %d, i_row: %d, i_time: %d, radius: %4.3f, TRD_time0 - TRD_drift_time: %4.3f, vD_calib: %4.3f \n",i_det,i_column,i_row,i_time,TRD_radius_uncalib,TRD_time0 - TRD_drift_time,vD_calib);

			    TVector3 TV3_TRD_hit, TV3_TRD_hit_uncalib, TV3_TRD_hit_det_angle;
			    TV3_TRD_hit.SetXYZ(glb_align[0],glb_align[1],glb_align[2]);
                            TV3_TRD_hit_uncalib.SetXYZ(glb_align_uncalib[0],glb_align_uncalib[1],glb_align_uncalib[2]);
                            TV3_TRD_hit_det_angle.SetXYZ(glb_vec_align[0],glb_vec_align[1],0.0); // glb_vec is from {0,0,0} to center of detector, set z component to 0 to get vector perpendicular to detector plane
                            TV3_TRD_hits.push_back(TV3_TRD_hit);
                            TV3_TRD_hits_uncalib.push_back(TV3_TRD_hit_uncalib);
			    TV3_TRD_hits_det_angle.push_back(TV3_TRD_hit_det_angle);
			    vec_TRD_hits_ADC_value.push_back(ADC_amplitude);
			    vec_TRD_hits_time.push_back(((Double_t)i_time)*0.1); // time in mus
			    vec_TRD_det_lay_row_col[0].push_back(i_det);
			    vec_TRD_det_lay_row_col[1].push_back(i_layer);
			    vec_TRD_det_lay_row_col[2].push_back(i_row);
			    vec_TRD_det_lay_row_col[3].push_back(i_column);

			    vec_ADC_time_bins[i_time]            = ADC_amplitude;
                            vec_points_time_bins[i_time]         = TV3_TRD_hit;
                            vec_points_time_bins_uncalib[i_time] = TV3_TRD_hit_uncalib;
                            //cout << "vec calib: " << vec_points_time_bins[i_time] << endl;
                            //cout << "vec uncalib: " << vec_points_time_bins_uncalib[i_time] << endl;
                        }
			//----------------------
		    }


                } // end of time loop
		//--------------------------------------------------



                //--------------------------------------------------
		Int_t    max_time_bin = -1;
                Double_t max_ADC_val  = 0.0;
		if(ADC_amplitude_sum_times > 0.0)
		{
		    for(Int_t i_time = 0; i_time < N_times; i_time++)
		    {
			if(arr[i_time] > max_ADC_val)
			{
			    max_ADC_val  = arr[i_time];
                            max_time_bin = i_time;
			}
		    }
		}
                //--------------------------------------------------



		//------------------------------------------------------------
		// Determine reference position (time bin with highest ADC value) for one single pad
		// Fill vector with ADC values for all 24 time bins for this pad
		if(ADC_amplitude_sum_times > 0.0 && max_time_bin > 0)
		{
		    //----------------------
		    // Calculate global position for fired TRD pad
		    Float_t              TRD_time0        = fGeo         ->GetTime0(i_layer); // in cm
		    AliTRDpadPlane*      padplane         = fGeo         ->GetPadPlane(i_det);
		    Double_t             TRD_col_end      = padplane     ->GetColEnd();
		    Double_t             TRD_row_end      = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
		    Double_t             TRD_row_start    = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
		    Double_t             TRD_row_end_ROC  = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
		    Double_t             TRD_col_spacing  = padplane     ->GetColSpacing();
		    Double_t             TRD_row_spacing  = padplane     ->GetRowSpacing();
		    Double_t             TRD_col_pos      = padplane     ->GetColPos(i_column);
		    Double_t             TRD_row_pos      = padplane     ->GetRowPos(i_row);       // fPadRow[row] + fPadRowSMOffset;
		    Double_t             TRD_row_pos_ROC  = padplane     ->GetRowPosROC(i_row);    // fPadRow[row]; = pad border position
		    Double_t             TRD_col_size     = padplane     ->GetColSize(i_column);
		    Double_t             TRD_row_size     = padplane     ->GetRowSize(i_row);
		    Float_t              TRD_loc_Y        = TRD_col_pos + 0.5*TRD_col_size;
		    Double_t             TRD_drift_time   = ((Double_t)(max_time_bin + 0.5))*TRD_time_per_bin*ChamberVdrift->GetValue(i_det); // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum
                    //Double_t             TRD_drift_time   = ((Double_t)12)*TRD_time_per_bin*ChamberVdrift->GetValue(i_det); // 100 ns per time bin, 1.56 cm/mus drift velocity, 3 cm drift length at maximum

		    Double_t kX0shift     = AliTRDgeometry::AnodePos(); //[cm]
                    Double_t fZShiftIdeal = 0.5 * (TRD_row_start + TRD_row_end);
		    // Get the detector wise defined calibration values
		    Double_t fCalVdriftDetValue = ChamberVdrift ->GetValue(i_det);
		    Double_t fCalT0DetValue     = ChamberT0     ->GetValue(i_det);
		    //Double_t fCalExBDetValue    = ChamberExB    ->GetValue(i_det);

		    // Retrieve calibration values
		    // drift velocity
		    Double_t vd  = fCalVdriftDetValue * 1.0;
		    // t0
		    Double_t t0  = fCalT0DetValue     + 0.0;
                    Double_t fSamplingFrequency = fParam->GetSamplingFrequency();
		    t0 /= fSamplingFrequency;
		    // ExB correction
		    //Double_t exb = fCalExBDetValue;//AliTRDCommonParam::Instance()->GetOmegaTau(vd);

		    //Float_t lorentz_angle_corr_y = ChamberExB->GetValue(i_det)*TRD_drift_time;
		    //TRD_loc_Y -= Sign_magnetic_field*lorentz_angle_corr_y;

		    //printf("ExB correction: %f \n",Sign_magnetic_field*lorentz_angle_corr_y);

		    std::vector<Int_t> vec_det_info;
		    vec_det_info.resize(4);
		    vec_det_info[0] = i_det;
		    vec_det_info[1] = i_layer;
		    vec_det_info[2] = i_row;
		    vec_det_info[3] = i_column;

		    // Determine global position of TRD hit
		    Double_t             loc_ref[3]           = {TRD_time0 - TRD_drift_time,TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};
		    //Double_t             loc_ref[3]           = {kX0shift - (TRD_drift_time - TRD_time0),TRD_loc_Y,TRD_row_pos - TRD_row_size/2.0};

		    //printf("fZShiftIdeal: %f, old x: %f, new x: %f \n",fZShiftIdeal,TRD_time0 - TRD_drift_time,kX0shift - (TRD_drift_time - TRD_time0));

		    Double_t             glb_ref[3]           = {0.0,0.0,0.0};
		    fGeo ->RotateBack(i_det,loc_ref,glb_ref);

                    // Apply alignment
		    Double_t glb_ref_align_sec[3];
                    Double_t glb_ref_align[3];

		    TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_ref,glb_ref_align_sec);
                    TM_TRD_rotation_det[i_det].LocalToMaster(glb_ref_align_sec,glb_ref_align);

		    // Determine global pointing vector of TRD plane (vector perpendicular to plane in global system)
		    Double_t             loc_vec_ref[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
		    Double_t             glb_vec_ref[3]           = {0.0,0.0,0.0};
		    fGeo ->RotateBack(i_det,loc_vec_ref,glb_vec_ref);

		    // Apply alignment
                    Double_t glb_vec_ref_align_sec[3];
		    Double_t glb_vec_ref_align[3];
                    TM_TRD_rotation_sector[i_sector].LocalToMaster(glb_vec_ref,glb_vec_ref_align_sec);
		    TM_TRD_rotation_det[i_det].LocalToMaster(glb_vec_ref_align_sec,glb_vec_ref_align);

		    TVector3 TV3_TRD_hit_middle, TV3_TRD_hit_det_angle_middle;
		    TV3_TRD_hit_middle.SetXYZ(glb_ref_align[0],glb_ref_align[1],glb_ref_align[2]);
		    TV3_TRD_hit_det_angle_middle.SetXYZ(glb_vec_ref_align[0],glb_vec_ref_align[1],0.0); // glb_vec is from {0,0,0} to center of detector, set z component to 0 to get vector perpendicular to detector plane
		    TV3_TRD_hits_middle.push_back(TV3_TRD_hit_middle);
		    TV3_TRD_hits_det_angle_middle.push_back(TV3_TRD_hit_det_angle_middle);
		    vec_TRD_hits_det_lay_row_col.push_back(vec_det_info);

		    vec_TRD_hits_ADC_values_time.push_back(vec_ADC_time_bins);
                    vec_TRD_hits_points_time.push_back(vec_points_time_bins);
                    vec_TRD_hits_points_time_uncalib.push_back(vec_points_time_bins_uncalib);

                    //if(digit_counter > 66652) continue;
                    AS_Digit = AS_Event ->createTRD_digit();
                    AS_Digit ->sethit_ids(x_TRD,y_TRD);
                    AS_Digit ->setdca_to_track(0.0,0.0,0.0,0.0);
                    AS_Digit ->setImpactAngle(0.0);
                    for(Int_t i_time = 0; i_time < (Int_t)vec_ADC_time_bins.size(); i_time++)
                    //for(Int_t i_time = 0; i_time < 2; i_time++)
                    {
                        Double_t ADC_value  = vec_ADC_time_bins[i_time];
                        //printf("i_time: %d, ADC: %4.3f \n",i_time,ADC_value);
                        AS_Digit ->setADC_time_value(i_time,(Short_t)ADC_value);
                        AS_Digit ->set_pos(i_time,vec_points_time_bins_uncalib[i_time].X(),vec_points_time_bins_uncalib[i_time].Y(),vec_points_time_bins_uncalib[i_time].Z());
                        //if(i_time == 0)
                        //{
                            //printf("1st i_time: %d, counter: %d, i_sector: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_time,digit_counter,i_sector,AS_Digit->get_pos(i_time,0),AS_Digit->get_pos(i_time,1),AS_Digit->get_pos(i_time,2));
                            //AS_Digit              = AS_Event ->getTRD_digit(digit_counter);
                            //printf("2nd i_time: %d, counter: %d, i_sector: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_time,digit_counter,i_sector,AS_Digit->get_pos(i_time,0),AS_Digit->get_pos(i_time,1),AS_Digit->get_pos(i_time,2));
                        //}
                    }



                    //-------------
                    // Test
                    AS_Digit              = AS_Event ->getTRD_digit(digit_counter);
                    Int_t    bsector      = AS_Digit ->get_sector();

                    //for(Int_t i_time = 0; i_time < (Int_t)vec_ADC_time_bins.size(); i_time++)
                    //for(Int_t i_time = 0; i_time < 1; i_time++)
                    //{
                    //    printf("test i_time: %d, i_digit: %d, i_sector: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_time,digit_counter,bsector,AS_Digit->get_pos(i_time,0),AS_Digit->get_pos(i_time,1),AS_Digit->get_pos(i_time,2));
                    //}
                    //-------------
                    //printf(" \n");



                    digit_counter++;
		}
		//------------------------------------------------------------


	    } // end of row loop
	} // end of column loop
    } // end of detector loop


    //printf("digit_counter: %d \n",digit_counter);

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        AS_Event ->setADC_sum_det(i_det,sum_full_ADC_digit_det[i_det]);
    }
    //cout << "max_N_columns: " << max_N_columns << ", max_N_rows: " << max_N_rows << endl;
    //-----------------------------------------------------------------


    //-----------------------------------------------------------------
    // Track loop
    //cout << "" << endl;
    //cout << "-----------------------------------------------------------------" << endl;
    //cout << "Start matching " << N_tracks << " TPC tracks with " << TV3_TRD_hits_middle.size() << " TRD pads" << endl;
    N_good_tracks = 0;
    Int_t N_matched_TRD_hits_total = 0;
    for(Int_t iTracks = 0; iTracks < N_tracks; iTracks++)
    {
	//---------------------------------------------------------------
	// Gather track information

        //cout << "iTracks," << iTracks << endl;

	// We always want the ESD track
	AliESDtrack* track = fESD->GetTrack(iTracks);
	if(!track)
	{
	    printf("ERROR: Could not receive track %d\n", iTracks);
	    continue;
	}


        if(!EsdTrackCuts->AcceptTrack(track)) continue;


        Double_t TRD_signal   = track ->GetTRDsignal(); // truncated mean signal?
        Double_t Track_pT     = track ->Pt();
        Double_t Track_p      = track ->P();
        Double_t p_vec[3];
        track->GetPxPyPz(p_vec);
        Int_t    charge       = track ->Charge();
        Double_t Track_phi    = track ->Phi();
	Double_t Track_theta  = track ->Theta();
	Double_t Track_eta    = track ->Eta();
	Double_t TPC_chi2     = track ->GetTPCchi2();
	Double_t TPC_signal   = track ->GetTPCsignal(); // dE/dx?
	Double_t TOF_signal   = track ->GetTOFsignal(); // time-of-flight?
        Double_t Track_length = track ->GetIntegratedLength();
	UShort_t N_TPC_cls    = track ->GetTPCNcls();

        Int_t pT_bin;
	for(Int_t i_pT = 0; i_pT < N_pT_bins; i_pT++)
	{
	    if(Track_pT >= pT_ranges[i_pT] && Track_pT < pT_ranges[i_pT+1])
	    {
                pT_bin = i_pT;
                break;
	    }
        }


	TLorentzVector TLV_p_vec;
	Double_t p_vec_energy = TMath::Sqrt(p_vec[0]*p_vec[0] + p_vec[1]*p_vec[1] + p_vec[2]*p_vec[2] + 0.938*0.938);
	TLV_p_vec.SetPxPyPzE(p_vec[0],p_vec[1],p_vec[2],p_vec_energy);
	//cout << "TLV_p_vec.P: " << TLV_p_vec.P() << ", P: " << Track_p << ", TLV_p_vec.Theta: " << TLV_p_vec.Theta() << ", Theta: " << Track_theta
	//<< ", TLV_p_vec.Phi: " << TLV_p_vec.Phi() << ", phi: " << Track_phi  << endl;


	ULong_t status = track->GetStatus();
	Int_t ITS_refit = 0;
        Int_t TPC_refit = 0;
        Int_t track_status = 0;
	if(((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit))
	{
	    ITS_refit = 1;
	    track_status |= 1 << 0; // setting bit 0 to 1
	}
	if(((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit))
	{
	    TPC_refit = 1;
	    track_status |= 1 << 1; // setting bit 1 to 1
	}


          // e = 0, muon = 1, pion = 2, kaon = 3, proton = 4
	Double_t Track_PID[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};


        // nSigma TPC
	Track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron);
	Track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kMuon);
	Track_PID[2] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
	Track_PID[3] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
	Track_PID[4] = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);

        // nSigma TOF, -999 in case there is no TOF hit
	Track_PID[5] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kElectron);
	Track_PID[6] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kMuon);
	Track_PID[7] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
	Track_PID[8] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
	Track_PID[9] = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);


	Float_t track_xy_impact,track_z_impact;
	track->GetImpactParameters(track_xy_impact,track_z_impact);
	Double_t track_total_impact = TMath::Sqrt(track_xy_impact*track_xy_impact + track_z_impact*track_z_impact);


	Double_t TRD_ADC_bin_width = 100.0;

	//-------------------
	Int_t N_ITS_cls = 0;
	for(Int_t i_ITS_layer = 0; i_ITS_layer < 6; ++i_ITS_layer)
	{
	    if(track ->HasPointOnITSLayer(i_ITS_layer))
	    {
		N_ITS_cls |= 1 << i_ITS_layer; // setting bit i_ITS_layer to 1
	    }
	}
	//-------------------


	TLorentzVector TL_vec;
	TL_vec.SetPtEtaPhiM(Track_pT,Track_eta,Track_phi,0.1349766);
	AS_Track  = AS_Event ->createTrack();
        AS_Track  ->clearTRD_digit_list();
        AS_Track  ->clearOfflineTrackletList();
	AS_Track  ->set_TLV_part(TL_vec);
	AS_Track  ->setdca(((Double_t)charge)*track_total_impact);
	AS_Track  ->setnsigma_e_TPC(Track_PID[0]);
	AS_Track  ->setnsigma_e_TOF(Track_PID[5]);
	AS_Track  ->setnsigma_pi_TPC(Track_PID[2]);
	AS_Track  ->setnsigma_pi_TOF(Track_PID[7]);
	AS_Track  ->setnsigma_K_TPC(Track_PID[3]);
	AS_Track  ->setnsigma_K_TOF(Track_PID[8]);
	AS_Track  ->setnsigma_p_TPC(Track_PID[4]);
	AS_Track  ->setnsigma_p_TOF(Track_PID[9]);
	AS_Track  ->setTRDSignal(TRD_signal);
	AS_Track  ->setNTPCcls(N_TPC_cls);
	AS_Track  ->setNITScls(N_ITS_cls);
	AS_Track  ->setStatus(track_status);
	AS_Track  ->setTPCchi2(TPC_chi2);
	AS_Track  ->setTPCdEdx(TPC_signal);
	AS_Track  ->setTOFsignal(TOF_signal);
        AS_Track  ->setTrack_length(Track_length);


	//-------------------
	// Get TRD information
	Int_t N_TRD_cls = 0;
	Double_t TRD_sum_ADC = 0.0;
	Long64_t TRD_layer_info[6]; // each value stores all 8 time slices, Long64_t has 8 byte = 8*8 bit = 64 bit in total -> one byte = 8 bit = 256 per time bin
	memset(TRD_layer_info, 0, sizeof(TRD_layer_info)); // for automatically-allocated arrays


	//Helix
        FillHelix(track,magF);

        AS_Track ->setHelix(aliHelix.fHelix[0],aliHelix.fHelix[1],aliHelix.fHelix[2],aliHelix.fHelix[3],aliHelix.fHelix[4],aliHelix.fHelix[5],aliHelix.fHelix[6],aliHelix.fHelix[7],aliHelix.fHelix[8]);

	Int_t N_track_TRD_tracklets     = track     ->GetTRDntracklets();
	//printf("N_track_TRD_tracklets: %d \n",N_track_TRD_tracklets);


	N_good_tracks++;

    } // End of TPC track loop
    //cout << "Tracks matched" << endl;



    //-----------------------------------
    // Make tracklets
    Double_t Delta_x        = 1.8; // 1.8;
    Double_t Delta_z        = 10.0;
    Double_t factor_missing = 1.0;
    vec_ADC_val.clear();
    vec_ADC_val.resize(540);
    Make_clusters_and_get_tracklets_fit(Delta_x,Delta_z,factor_missing);
    //printf("Make_clusters_and_get_tracklets_fit done \n");

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        //printf("i_det: %d \n",i_det);
        for(Int_t i_trkl = 0; i_trkl < (Int_t)vec_self_tracklet_fit_points[i_det].size(); i_trkl++)
        {
            for (Int_t i_tbn = 0; i_tbn < 30; i_tbn++)
            {
                ADC_val[i_tbn] = {-999.0};
            }

            //printf("i_trkl: %d \n",i_trkl);
            for(Int_t i_AB = 0; i_AB < 2; i_AB++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    //printf("i_AB: %d, i_xyz: %d \n",i_AB,i_xyz);
                    if(i_AB == 0) TV3_trkl_offset[i_xyz] = vec_self_tracklet_fit_points[i_det][i_trkl][i_AB][i_xyz];
                    if(i_AB == 1) TV3_trkl_dir[i_xyz]    = vec_self_tracklet_fit_points[i_det][i_trkl][i_AB][i_xyz];
                }
            }

            TV3_trkl_dir -= TV3_trkl_offset;

            TRD_ST_Tracklet = TRD_ST_Event ->createTracklet(); // online tracklet
            TRD_ST_Tracklet  ->set_TRD_det(i_det);
            TRD_ST_Tracklet  ->set_TV3_offset(TV3_trkl_offset);
            TRD_ST_Tracklet  ->set_TV3_dir(TV3_trkl_dir);
            for(Int_t i_timebin = 0; i_timebin < (Int_t)vec_ADC_val[i_det][i_trkl].size(); i_timebin++) // ALEX
            {
                ADC_val[i_timebin] = vec_ADC_val[i_det][i_trkl][i_timebin];
                TRD_ST_Tracklet  ->set_ADC_val(i_timebin,ADC_val[i_timebin]);
            }

        } // end tracklet loop
    } // end detector loop

    printf("Tracklet information filled \n");
    //-----------------------------------



    //------------------------------------------
    // Loop over all tracks -> fill tree for tracklets tree
    UShort_t NumTracks = AS_Event ->getNumTracks(); // number of tracks in this event

    for(UShort_t i_track = 0; i_track < NumTracks; ++i_track) // loop over all tracks of the actual event
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->createTrack(); // TPC track
        //cout << "i_track: " << i_track << ", of " << NumTracks << endl;
        AS_Track      = AS_Event ->getTrack( i_track ); // take the track
        TRD_ST_TPC_Track ->setnsigma_e_TPC( AS_Track ->getnsigma_e_TPC() );
        TRD_ST_TPC_Track ->setnsigma_pi_TPC( AS_Track ->getnsigma_pi_TPC() );
        TRD_ST_TPC_Track ->setnsigma_p_TPC( AS_Track ->getnsigma_p_TPC() );
        TRD_ST_TPC_Track ->setnsigma_K_TPC( AS_Track ->getnsigma_K_TPC() );
        TRD_ST_TPC_Track ->setnsigma_e_TOF( AS_Track ->getnsigma_e_TOF() );
        TRD_ST_TPC_Track ->setnsigma_pi_TOF( AS_Track ->getnsigma_pi_TOF() );
        TRD_ST_TPC_Track ->setnsigma_K_TOF( AS_Track ->getnsigma_K_TOF() );
        TRD_ST_TPC_Track ->setTRDSignal( AS_Track ->getTRDSignal() );
        TRD_ST_TPC_Track ->setTRDsumADC( AS_Track ->getTRDsumADC() );
        TRD_ST_TPC_Track ->setdca( AS_Track ->getdca() );  // charge * distance of closest approach to the primary vertex
        TRD_ST_TPC_Track ->set_TLV_part( AS_Track ->get_TLV_part() );
        TRD_ST_TPC_Track ->setNTPCcls( AS_Track ->getNTPCcls() );
        TRD_ST_TPC_Track ->setNTRDcls( AS_Track ->getNTRDcls() );
        TRD_ST_TPC_Track ->setNITScls( AS_Track ->getNITScls() );
        TRD_ST_TPC_Track ->setTPCchi2( AS_Track ->getTPCchi2() );
        TRD_ST_TPC_Track ->setTPCdEdx( AS_Track ->getTPCdEdx() );
        TRD_ST_TPC_Track ->setTOFsignal( AS_Track ->getTOFsignal() ); // in ps (1E-12 s)
        TRD_ST_TPC_Track ->setTrack_length( AS_Track ->getTrack_length() );

        for(Int_t i_helix_par = 0; i_helix_par < 9; i_helix_par++)
        {
            helix_par[i_helix_par] =  AS_Track ->getHelix_param(i_helix_par);
        }

        TRD_ST_TPC_Track ->setHelix(helix_par[0],helix_par[1],helix_par[2],helix_par[3],helix_par[4],helix_par[5],helix_par[6],helix_par[7],helix_par[8]);
    }

    printf("Track information filled \n");
    //------------------------------------------



    //-----------------------------------
    //Tree_AS_Event ->Fill();
    Tree_TRD_ST_Event ->Fill(); // new tracklets tree to be filled
    //Long64_t size_of_tree = Tree_TRD_ST_Event ->GetEntries();
    //printf("Event: %d, tree filled, size of tree: %lld \n",N_good_events,size_of_tree);
    printf("Tree filled \n");
    //-----------------------------------


#if 0
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    AliInfoF("Processing event %i", fEventNoInFile);
    AliInfoF("Memory: RSS: %3ld VMEM: %3ld",procInfo.fMemResident/1024,procInfo.fMemVirtual/1024);
#endif
    printf("Event: %d \n",N_good_events);

    if (fLocalMode) {
      ProcInfo_t procInfo;
      gSystem->GetProcInfo(&procInfo);
      AliInfoF("Processing event %i", N_good_events);
      AliInfoF("Memory: RSS: %3ld VMEM: %3ld",procInfo.fMemResident/1024,procInfo.fMemVirtual/1024);
    }

    N_good_events++;

}


//----------------------------------------------------------------------------------------
TVector3 Ali_make_tracklets_from_digits::calculate_point_on_Straight_dca_to_Point_2D(TVector3 &base, TVector3 &dir, TVector3 &point)
{
  // calculates the TVector3 on the straight line which is closest to point

    
    TVector3 point2d;
    TVector3 base2d;
    TVector3 dir2d;

    point2d.SetXYZ(point[0],point[1],0.0);
    base2d.SetXYZ(base[0],base[1],0.0);
    dir2d.SetXYZ(dir[0],dir[1],0.0);

    TVector3 diff = point2d - base2d;
    TVector3 dir_norm = dir2d;
    dir_norm *= (1.0/dir2d.Mag());
    Double_t proj_val = diff.Dot(dir_norm);
    TVector3 proj_dir = dir_norm;
    proj_dir *= proj_val;

    TVector3 dist_vec = proj_dir + base2d;

    return dist_vec;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
Double_t Ali_make_tracklets_from_digits::Calc_SVD_tracklet(Int_t i_det, Int_t i_trkl)
{
    // https://www.codefull.net/2015/06/3d-line-fitting/ -> MATLAB code

    TVector3 TV3_point;
    TVector3 TV3_mean(0.0,0.0,0.0);
    vector<TVector3> arr_TV3_points;
    vector<TVector3> arr_TV3_points_mean;

    Double_t SVD_chi2 = 0.0;  //chi2 of the tracklet

    Int_t number_ok_clusters = 0;
    for(Int_t i_time = 0; i_time < (Int_t)vec_connected_clusters[i_det][i_trkl].size(); i_time++)
    {
        TV3_point.SetXYZ(vec_connected_clusters[i_det][i_trkl][i_time][0],vec_connected_clusters[i_det][i_trkl][i_time][1],vec_connected_clusters[i_det][i_trkl][i_time][2]);
        TV3_mean += TV3_point;
        arr_TV3_points.push_back(TV3_point);
        number_ok_clusters++;
    }

    //printf("number_ok_clusters: %d \n",number_ok_clusters);
    TArrayD arr_data_points(number_ok_clusters*3);

    // Calculate mean and subtract it
    if(number_ok_clusters > 0) TV3_mean *= 1.0/((Double_t)number_ok_clusters);

    //if(i_det == 244 && i_trkl == 4)
    //{
    //    printf(" --> number_ok_clusters: %d, mean vec: {%4.3f, %4.3f, %4.3f} \n",number_ok_clusters,TV3_mean.X(),TV3_mean.Y(),TV3_mean.Z());
    //}

    for(Int_t i_point = 0; i_point < number_ok_clusters; i_point++)
    {
        arr_TV3_points_mean.push_back(arr_TV3_points[i_point] - TV3_mean);
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            Int_t i_data = i_xyz + 3*i_point;
            arr_data_points[i_data] = arr_TV3_points_mean[i_point][i_xyz];
        }
    }

    TMatrixD A_Matrix(number_ok_clusters,3);
    A_Matrix.SetMatrixArray(arr_data_points.GetArray());
    TDecompSVD tSVD(A_Matrix);
    tSVD.Decompose();
    TMatrixD V_Matrix = tSVD.GetV();

    TVector3 TV3_fit_dir(V_Matrix[0][0],V_Matrix[1][0],V_Matrix[2][0]);
    TV3_fit_dir *= 1.0/TV3_fit_dir.Mag();

    TV3_SVD_tracklet_offset = TV3_mean;
    TV3_SVD_tracklet_dir    = TV3_fit_dir;

    for(Int_t i_point = 0; i_point < number_ok_clusters; i_point++)
    {
        TVector3 testpoint = calculate_point_on_Straight_dca_to_Point_2D(TV3_SVD_tracklet_offset, TV3_SVD_tracklet_dir, arr_TV3_points[i_point]);

        Double_t dist_DCA_XY = TMath::Sqrt(TMath::Power(arr_TV3_points[i_point][0] - testpoint[0],2) + TMath::Power(arr_TV3_points[i_point][1] - testpoint[1],2));
        //Double_t dist_DCA_Z  = fabs(arr_TV3_points[i_point][2] - testpoint[2]);
        Double_t dist_weighted = TMath::Sqrt(TMath::Power(dist_DCA_XY/0.7,2.0));// + TMath::Power(dist_DCA_Z/7.5,2.0));
        SVD_chi2 += abs(dist_weighted);

    }

    SVD_chi2 = SVD_chi2/number_ok_clusters;
    return SVD_chi2;
    //printf("SVD_chi2 = %4.3f \n",SVD_chi2);

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_make_tracklets_from_digits::Make_clusters_and_get_tracklets_fit(Double_t Delta_x, Double_t Delta_z, Double_t factor_missing)
{
    //printf("TTRD_ST_Make_Tracklets::Make_clusters_and_get_tracklets_fit() \n");

    //Reset();

    vector< vector< vector< vector<Double_t> > > > vec_all_TRD_digits;
    vector< vector< vector< vector<Double_t> > > > vec_all_TRD_digits_clusters;
    vector< vector< vector<Int_t> > >              vec_used_clusters;


    vec_all_TRD_digits.resize(540);
    vec_all_TRD_digits_clusters.resize(540);
    vec_used_clusters.resize(540);

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_all_TRD_digits[i_det].resize(30); // time bins
        vec_all_TRD_digits_clusters[i_det].resize(30); // time bins
    	vec_used_clusters[i_det].resize(30); // time bins
    }
    vector<Double_t> vec_digit_data; // x,y,z,ADC
    vector<Double_t> vec_digit_cluster_data; // x,y,z,ADC,N_digits
    vec_digit_data.resize(4);
    vec_digit_cluster_data.resize(5);

    // Fill all the information in the hierachy of detectors and time bins
    Int_t    N_TRD_digits  = AS_Event ->getNumTRD_digits();
    printf(" --> N_TRD_digits: %d \n",N_TRD_digits);

    for(Int_t i_digit = 0; i_digit < N_TRD_digits; i_digit++)
    {
        AS_Digit              = AS_Event ->getTRD_digit(i_digit);
        Int_t    layer        = AS_Digit ->get_layer();
        Int_t    sector       = AS_Digit ->get_sector();
        Int_t    column       = AS_Digit ->get_column();
        Int_t    stack        = AS_Digit ->get_stack();
        Int_t    row          = AS_Digit ->get_row();
        Int_t    detector     = AS_Digit ->get_detector(layer,stack,sector);

        //printf("i_digit: %d, out of %d, sector: %d \n",i_digit,N_TRD_digits,sector);


        Double_t radius_prev = 0.0;
        for(Int_t i_time = 0; i_time < 30; i_time++)
        {
            Float_t ADC = (Float_t)AS_Digit ->getADC_time_value(i_time);
            //if(i_time == 0) printf("read i_time: %d, counter: %d, i_sector: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_time,i_digit,sector,AS_Digit->get_pos(i_time,0),AS_Digit->get_pos(i_time,1),AS_Digit->get_pos(i_time,2));


            if(ADC < 0.0) continue;
            //if(i_time == 0 ) printf("i_digit: %d, i_sector: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_digit,sector,AS_Digit->get_pos(i_time,0),AS_Digit->get_pos(i_time,1),AS_Digit->get_pos(i_time,2));
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                digit_pos[i_xyz] = AS_Digit ->get_pos(i_time,i_xyz);
               
                vec_digit_data[i_xyz] = digit_pos[i_xyz];
            }
            vec_digit_data[3] = ADC;


            vec_all_TRD_digits[detector][i_time].push_back(vec_digit_data);

            Double_t radius = TMath::Sqrt(TMath::Power(digit_pos[0] ,2) + TMath::Power(digit_pos[1] ,2));
            //if(i_time > 0 && radius_prev < radius) printf("det: %d, layer: %d, sector: %d, stack: %d, i_time: %d, radius: %4.3f \n",detector,layer,sector,stack,i_time,radius);
            radius_prev = radius;

        }
    }


    //-------------------------------------------------------
    // Make clusters for each detector and time bin
    // Individial digits -> clusters/time bin
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        Int_t sector = (Int_t)(i_det/30);
        Int_t stack  = (Int_t)(i_det%30/6);
        Int_t layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        //printf("i_det: %d \n",i_det);
        for(Int_t i_time = 0; i_time < 30; i_time++)
        {
            //if(!(i_det == 0 && i_time == 0)) continue;
            //for(Int_t i_digit = 0; i_digit < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit++)
            //{
            //    printf("i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, ADC: %4.3f \n",i_digit,vec_all_TRD_digits[i_det][i_time][i_digit][0],vec_all_TRD_digits[i_det][i_time][i_digit][1],vec_all_TRD_digits[i_det][i_time][i_digit][2],vec_all_TRD_digits[i_det][i_time][i_digit][3]);
            //}

            //cout << "Sorting vector" << endl;
            // First order the vector from high to low ADC values

            // vec_all_TRD_digits // 540 chambers, 24 time bins, (x,y,z,ADC)
            std::sort(vec_all_TRD_digits[i_det][i_time].begin(),vec_all_TRD_digits[i_det][i_time].end(),sortcol_first); // large values to small values, last column sorted via function sortcol

            //for(Int_t i_digit = 0; i_digit < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit++)
            //{
            //    printf("i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, ADC: %4.3f \n",i_digit,vec_all_TRD_digits[i_det][i_time][i_digit][0],vec_all_TRD_digits[i_det][i_time][i_digit][1],vec_all_TRD_digits[i_det][i_time][i_digit][2],vec_all_TRD_digits[i_det][i_time][i_digit][3]);
            //}


            vector<Int_t> arr_used_digits;
            arr_used_digits.clear();
            arr_used_digits.resize((Int_t)vec_all_TRD_digits[i_det][i_time].size());
            for(Int_t i_digit_max = 0; i_digit_max < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit_max++)
            {
                arr_used_digits[i_digit_max] = 0;
            }


            //-------------------------------
            // Pre cleaning procedure -> flag digits which are alone, no other digits around in same time window

            // Start from the maximum ADC value(s)
            for(Int_t i_digit_max = 0; i_digit_max < ((Int_t)vec_all_TRD_digits[i_det][i_time].size()); i_digit_max++)
            {
                //if(arr_used_digits[i_digit_max]) continue;
                //arr_used_digits[i_digit_max] = 1;

                Double_t pos_ADC_max[4] = {vec_all_TRD_digits[i_det][i_time][i_digit_max][0],vec_all_TRD_digits[i_det][i_time][i_digit_max][1],vec_all_TRD_digits[i_det][i_time][i_digit_max][2],vec_all_TRD_digits[i_det][i_time][i_digit_max][3]};
                Int_t N_digits_added = 1;

                // Get all other digits within a certain radius
                for(Int_t i_digit_sub = 0; i_digit_sub < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit_sub++)
                {
                    if(i_digit_sub == i_digit_max) continue;
                    //if(arr_used_digits[i_digit_sub]) continue;

                    Double_t pos_ADC_sub[4] = {vec_all_TRD_digits[i_det][i_time][i_digit_sub][0],vec_all_TRD_digits[i_det][i_time][i_digit_sub][1],vec_all_TRD_digits[i_det][i_time][i_digit_sub][2],vec_all_TRD_digits[i_det][i_time][i_digit_sub][3]};
                    Double_t dist_digits_XY = TMath::Sqrt(TMath::Power(pos_ADC_max[0] - pos_ADC_sub[0],2) + TMath::Power(pos_ADC_max[1] - pos_ADC_sub[1],2));
                    Double_t dist_digits_Z  = fabs(pos_ADC_max[2] - pos_ADC_sub[2]);
                    if(dist_digits_XY < 1.0 && dist_digits_Z  < 5.0)
                    {
                        N_digits_added++;
                        break;
                    }

                }

                if(N_digits_added == 1) arr_used_digits[i_digit_max] = 1;

                //if(i_det == 390 && i_time == 0 && pos_ADC_max[0] < -44.0)  printf("i_digit_max: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_digit_max,pos_ADC_max[0],pos_ADC_max[1],pos_ADC_max[2]);

            } // end of digit max loop
            //-------------------------------


            // Start from the maximum ADC value(s)
            Double_t baseline = 10.0;
            for(Int_t i_digit_max = 0; i_digit_max < ((Int_t)vec_all_TRD_digits[i_det][i_time].size() - 1); i_digit_max++)
            {
                if(arr_used_digits[i_digit_max]) continue;
                arr_used_digits[i_digit_max] = 1;

                Double_t pos_ADC_max[4] = {vec_all_TRD_digits[i_det][i_time][i_digit_max][0],vec_all_TRD_digits[i_det][i_time][i_digit_max][1],vec_all_TRD_digits[i_det][i_time][i_digit_max][2],vec_all_TRD_digits[i_det][i_time][i_digit_max][3] - baseline};
                Int_t N_digits_added = 1;
                Double_t Sum_ADC_weight = pos_ADC_max[3]; // ADC as weight
                if(Sum_ADC_weight < 0.0) continue;
                Double_t pos_ADC_sum[4] = {pos_ADC_max[0]*pos_ADC_max[3],pos_ADC_max[1]*pos_ADC_max[3],pos_ADC_max[2]*pos_ADC_max[3],pos_ADC_max[3]}; // positions times ADC value [3]

                //printf("i_time: %d, i_digit_max: %d, ADC: %4.3f \n",i_time,i_digit_max,pos_ADC_max[3]);

                // Get all other digits within a certain radius
                for(Int_t i_digit_sub = (i_digit_max + 1); i_digit_sub < (Int_t)vec_all_TRD_digits[i_det][i_time].size(); i_digit_sub++)
                {
                    if(arr_used_digits[i_digit_sub]) continue;

                    Double_t pos_ADC_sub[4] = {vec_all_TRD_digits[i_det][i_time][i_digit_sub][0],vec_all_TRD_digits[i_det][i_time][i_digit_sub][1],vec_all_TRD_digits[i_det][i_time][i_digit_sub][2],vec_all_TRD_digits[i_det][i_time][i_digit_sub][3] - baseline};

                    Double_t ADC_sub = pos_ADC_sub[3];
                    if(ADC_sub < 0.0) continue;

                    Double_t dist_digits_XY = TMath::Sqrt(TMath::Power(pos_ADC_max[0] - pos_ADC_sub[0],2) + TMath::Power(pos_ADC_max[1] - pos_ADC_sub[1],2));
                    Double_t dist_digits_Z  = fabs(pos_ADC_max[2] - pos_ADC_sub[2]);
                    if(dist_digits_XY > 2.5)  continue; // 2.5
                    if(dist_digits_Z  > 10.0) continue; // 15.0

                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        pos_ADC_sum[i_xyz] += ADC_sub*pos_ADC_sub[i_xyz];
                    }
                    pos_ADC_sum[3] += pos_ADC_sub[3];

                    arr_used_digits[i_digit_sub] = 1;
                    Sum_ADC_weight += ADC_sub;
                    N_digits_added++;
                }

                if(Sum_ADC_weight <= 0.0) continue;
                // Calculate average position
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    pos_ADC_sum[i_xyz] /= Sum_ADC_weight;
                }





                for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
                {
                    vec_digit_cluster_data[i_xyzADC] = pos_ADC_sum[i_xyzADC];
                }
                vec_digit_cluster_data[4] = (Double_t)N_digits_added;

                //printf("i_time: %d, cluster pos: {%4.3f, %4.3f, %4.3f} \n",i_time,vec_digit_cluster_data[0],vec_digit_cluster_data[1],vec_digit_cluster_data[2]);

                vec_all_TRD_digits_clusters[i_det][i_time].push_back(vec_digit_cluster_data);
                vec_used_clusters[i_det][i_time].push_back(0);
            } // end of digit max loop
        } // end of timebin loop
    }
    printf("arrangement of clusters done \n");
    //-------------------------------------------------------



    //create "tracklets"
    //Int_t sector = (Int_t)(i_det/30);
    //Int_t stack  = (Int_t)(i_det%30/6);
    //Int_t layer  = i_det%6;
    //Int_t i_det = layer + 6*stack + 30*sector;


    //-------------------------------------------------------
    // Connect clusters within each chamber
    // Clusters/time bin -> connect them time bin wise

    vec_connected_clusters.clear(); // defined in Ana_Digits_functions.h as static
    vec_connected_clusters.resize(540); // i_det i_trkl i_point (up to 24 time bins, can be less) i_xyz

    vector< vector<Double_t> > vec_single_connected_clusters; // for one tracklet, i_point (up to 24 time bins, can be less) i_xyz
    vector<Double_t> vec_single_point; // x,y,z,ADC
    vec_single_point.resize(4); // x,y,z,ADC

    vector< vector<Int_t> > vec_N_clusters_self_tracklet_points;
    vec_N_clusters_self_tracklet_points.resize(540);

    Int_t min_nbr_cls = 10;


    for(Int_t i_det = 0; i_det < 540; i_det++) // is done chamber wise
    {

        for(Int_t i_time = 0; i_time < 30 - min_nbr_cls; i_time++) // is done chamber wise
        //for(Int_t i_time = 0; i_time < 1; i_time++) // is done chamber wise  ALEX
        {
            Int_t N_clusters = (Int_t)vec_all_TRD_digits_clusters[i_det][i_time].size();

            std::sort(vec_all_TRD_digits_clusters[i_det][i_time].begin(),vec_all_TRD_digits_clusters[i_det][i_time].end(),sortcol_first); // large values to small values, last column sorted via function sortcol


            //printf("i_det: %d, N_clusters: %d \n",i_det,N_clusters);

            vec_N_clusters_self_tracklet_points[i_det].resize(N_clusters);

            vector< vector<Int_t> > vec_cls_shared;
            vector<Int_t> vec_single_trkl_cls;

            for(Int_t i_cls = 0; i_cls < N_clusters; i_cls++) // loop over all clusters in one detector and for time bin 0
            {
                vec_single_connected_clusters.clear();

                vec_single_trkl_cls.push_back(i_cls);
                if(vec_used_clusters[i_det][i_time][i_cls]) continue;
                Int_t n_clusters_attached = 0;
                Double_t pos_ADC_max[4] = {vec_all_TRD_digits_clusters[i_det][i_time][i_cls][0],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][1],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][2],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][3]};

                Int_t N_digits_used = vec_all_TRD_digits_clusters[i_det][i_time][i_cls][4];
                if(N_digits_used <= 1) continue;


                Double_t radius = TMath::Sqrt(TMath::Power(vec_all_TRD_digits_clusters[i_det][i_time][i_cls][0],2) + TMath::Power(vec_all_TRD_digits_clusters[i_det][i_time][i_cls][1],2));

                //if(i_time == 0) printf("---> i_det: %d, i_cls: %d, time 0, pos: {%4.3f, %4.3f, %4.3f}, radius: %4.3f, ADC: %4.3f, N_digits_used: %d \n",i_det,i_cls,vec_all_TRD_digits_clusters[i_det][i_time][i_cls][0],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][1],vec_all_TRD_digits_clusters[i_det][i_time][i_cls][2],radius,vec_all_TRD_digits_clusters[i_det][i_time][i_cls][3],(Int_t)vec_all_TRD_digits_clusters[i_det][i_time][i_cls][4]);


                for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
                {
                    vec_single_point[i_xyzADC] = pos_ADC_max[i_xyzADC];
                }
                vec_single_connected_clusters.push_back(vec_single_point);


                Int_t    i_time_start    = i_time + 1;
                Double_t scale_fac_add   = 1.0;
                Int_t    missed_time_bin = 0;
                Double_t radius_prev     = 0.0;

                Double_t sum_cluster_quality = 0.0;
                for(Int_t i_time_sub = i_time_start; i_time_sub < 30; i_time_sub++)
                {
                    Double_t scale_fac = 1.0*scale_fac_add;
                    //printf("i_time_sub: %d, scale_fac: %4.3f \n",i_time_sub,scale_fac);
                    //if(i_time_sub == 0) scale_fac = factor_layer*scale_fac_add;
                    Int_t N_clusters_sub = (Int_t)vec_all_TRD_digits_clusters[i_det][i_time_sub].size();

                    Int_t best_sub_cluster = -1;
                    Double_t best_cluster_quality = 1000000.0;

                    for(Int_t i_cls_sub = 0; i_cls_sub < N_clusters_sub; i_cls_sub++)
                    {
                        if(vec_used_clusters[i_det][i_time_sub][i_cls_sub]) continue;

                        Double_t pos_ADC_sub[4] = {vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][0],vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][1],vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][2],vec_all_TRD_digits_clusters[i_det][i_time_sub][i_cls_sub][3]};

                        Double_t dist_clusters_XY = TMath::Sqrt(TMath::Power(pos_ADC_max[0] - pos_ADC_sub[0],2) + TMath::Power(pos_ADC_max[1] - pos_ADC_sub[1],2));
                        Double_t dist_clusters_Z  = fabs(pos_ADC_max[2] - pos_ADC_sub[2]);

                        if(dist_clusters_XY > scale_fac*Delta_x)  continue;
                        if(dist_clusters_Z  > Delta_z) continue;
                        //if(dist_clusters_XY > scale_fac*3.0)  continue;
                        //if(dist_clusters_Z  > 10.0) continue;

                        // Matching quality - chi2 like
                        Double_t sub_cluster_quality = TMath::Power(dist_clusters_XY/0.7,2.0) + TMath::Power(dist_clusters_Z/7.5,2.0);


                        if(sub_cluster_quality < best_cluster_quality)
                        {
                            best_cluster_quality = sub_cluster_quality;
                            best_sub_cluster     = i_cls_sub;
                        }
                    }

                    sum_cluster_quality += best_cluster_quality;
                    vec_single_trkl_cls.push_back(best_sub_cluster);

                    //cout << "best_sub_cluster: " << best_sub_cluster << ", i_time_sub: " << i_time_sub << ", i_layer_sub: " << i_layer_sub << endl;
                    if(best_sub_cluster < 0)
                    {
                        scale_fac_add *= factor_missing; // one time bin was missing, increase matching window
                        missed_time_bin++;
                        continue;
                    }

                    if(missed_time_bin > 3) break;
                    scale_fac_add = 1.0; // reset additional matching window factor once a match was found

                    //if(missed_time_bin > 3) break;
                    //scale_fac_add = 1.0; 

                    // Define new pos_ADC_max
                    for(Int_t i_xyzADC = 0; i_xyzADC < 4; i_xyzADC++)
                    {
                        pos_ADC_max[i_xyzADC] = vec_all_TRD_digits_clusters[i_det][i_time_sub][best_sub_cluster][i_xyzADC]; // The one connected is now the new "first" one for the next search step in time
                        vec_single_point[i_xyzADC] = pos_ADC_max[i_xyzADC];
                    }

                    vec_single_connected_clusters.push_back(vec_single_point);


                    vec_used_clusters[i_det][i_time_sub][best_sub_cluster] = 1;

                    n_clusters_attached++;
                } // end of time sub loop

                vec_N_clusters_self_tracklet_points[i_det][i_cls] = n_clusters_attached;


                vec_connected_clusters[i_det].push_back(vec_single_connected_clusters);
            } // end of cluster loop


        } // end of time loop
    }

    printf("connection of clusters within detector done \n");
    //-------------------------------------------------------



    //-------------------------------------------------------
    // Fit the time bin wise connected clusters


    // Is fitting the tracklets and doing the global fit through all first cluster points of all available layers
    //printf("TTRD_ST_Make_Tracklets::get_tracklets_fit((%d) \n",i_track);

    //fit merged digits with a straight line
    vec_self_tracklet_fit_points.clear();
    vec_self_tracklet_fit_points.resize(540);         //[i_det][i_trkl][i_start_stop][i_xyz]

    for(Int_t i_detector = 0; i_detector < 540; i_detector++)
    {
        vec_self_tracklet_fit_points[i_detector].resize((Int_t)vec_connected_clusters[i_detector].size());

        vec_ADC_val[i_detector].resize((Int_t)vec_connected_clusters[i_detector].size());

        for (Int_t i_trkl = 0; i_trkl < (Int_t)vec_connected_clusters[i_detector].size(); i_trkl++)
        {
            vec_self_tracklet_fit_points[i_detector][i_trkl].resize(2);
            for(Int_t i_start_stop = 0; i_start_stop < 2; i_start_stop++)
            {
                vec_self_tracklet_fit_points[i_detector][i_trkl][i_start_stop].resize(3);
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    vec_self_tracklet_fit_points[i_detector][i_trkl][i_start_stop][i_xyz] = -999.0;
                }
            }
        }
    }


    Int_t trkl_index_layer[6] = {0};
    Int_t trkl_index=0;
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        //if(i_det == 0) printf("fitting detector %d  \n",i_det);
        //if(i_det % 40 == 0) printf("fitting detector %d  \n",i_det);

        Int_t i_layer  = i_det%6;


        //printf("i_det: %d, N_tracklets: %d \n",i_det,(Int_t)vec_self_tracklet_points[i_det].size());

        //for(Int_t i_trkl = 0; i_trkl < (Int_t)vec_self_tracklet_points[i_det].size(); i_trkl++) // ALEX
        for(Int_t i_trkl = 0; i_trkl < (Int_t)vec_connected_clusters[i_det].size(); i_trkl++)
        {
            //printf("i_det: %d, i_trkl: %d, N_clusters: %d \n",i_det,i_trkl,vec_N_clusters_self_tracklet_points[i_det][i_trkl]);
            if((Int_t)vec_connected_clusters[i_det][i_trkl].size() < 10) continue; // ALEX

            //------------------------------
            Double_t SVD_chi2 = Calc_SVD_tracklet(i_det,i_trkl);

            if(SVD_chi2 > 0.5) continue; // 0.5 ALEX


            //TV3_SVD_tracklet_offset.Print();
            //TV3_SVD_tracklet_dir.Print();

            //vec_ADC_val[i_det][i_trkl].resize((Int_t)vec_self_tracklet_points[i_det][i_trkl].size());
            vec_ADC_val[i_det][i_trkl].resize((Int_t)vec_connected_clusters[i_det][i_trkl].size()); // ALEX


            //------------------------------

            //-------------------------------------------------------
            // Calculate tracklet base and direction vectors

            // Space point on straight line which is closes to first space point of fitted clusters
            TVector3 TV3_base_plane = vec_TV3_TRD_center_offset[i_det];
            TVector3 TV3_norm_plane = vec_TV3_TRD_center[i_det][2];
            TVector3 TV3_base_fit_t0 = intersect_line_plane(TV3_SVD_tracklet_offset,TV3_SVD_tracklet_dir,TV3_base_plane,TV3_norm_plane);

            Double_t TV3_base_fit_t0_radius = TMath::Sqrt(TMath::Power(TV3_base_fit_t0[0],2) + TMath::Power(TV3_base_fit_t0[1],2));

            //TV3_base_fit_t0.Print();
            //TV3_SVD_tracklet_dir.Print();


            TVector3 vec_AB[2];
            vec_AB[0] = TV3_base_fit_t0;
            vec_AB[1] = TV3_base_fit_t0 + TV3_SVD_tracklet_dir;
            //if(vec_AB[1].Mag() > vec_AB[0].Mag())  // changed sign -> correct physical direction, pointing inwards
            if(vec_AB[1].Perp() < vec_AB[0].Perp())  // changed sign -> "uncorrect" direction, pointing outwards
            {
                TV3_SVD_tracklet_dir *= -1.0;
                vec_AB[1] = TV3_base_fit_t0 + TV3_SVD_tracklet_dir;
            }


            for(Int_t i_AB = 0; i_AB < 2; i_AB++)
            {
                for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                {
                    vec_self_tracklet_fit_points[i_det][i_trkl][i_AB][i_xyz] = vec_AB[i_AB][i_xyz];
                }
            }


            // fill ADC val vector
            Int_t tbn_max = (Int_t)vec_connected_clusters[i_det][i_trkl].size();
            //Double_t radius_in  = TMath::Sqrt( TMath::Power(vec_connected_clusters[i_det][i_trkl][tbn_max-1][0],2) + TMath::Power(vec_connected_clusters[i_det][i_trkl][tbn_max-1][1],2) );
            //Double_t radius_out = TMath::Sqrt( TMath::Power(vec_connected_clusters[i_det][i_trkl][0][0],2) + TMath::Power(vec_connected_clusters[i_det][i_trkl][0][1],2) );

            for(Int_t i_timebin = 0; i_timebin < tbn_max; i_timebin++)
            {
                //printf("i_timebin: %d, value: %4.3f \n",i_timebin,vec_connected_clusters[i_det][i_trkl][i_timebin][3]);
                vec_ADC_val[i_det][i_trkl][i_timebin] = vec_connected_clusters[i_det][i_trkl][i_timebin][3]; // ALEX
            }

            //-------------------------------------------------------

            Double_t radius = TMath::Sqrt( TMath::Power(TV3_base_fit_t0[0],2) + TMath::Power(TV3_base_fit_t0[1],2) );
            //printf("amin: %4.3f, par: {%4.3f, %4.3f, %4.3f, %4.3f} \n",amin,parFit[0],parFit[1],parFit[2],parFit[3]);
            //printf("   --> radius: %4.3f, point first cluster: {%4.3f, %4.3f, %4.3f}, point line: {%4.3f, %4.3f, %4.3f} \n",radius,TV3_t0_point[0],TV3_t0_point[1],TV3_t0_point[2],TV3_base_fit_t0[0],TV3_base_fit_t0[1],TV3_base_fit_t0[2]);



            trkl_index_layer[i_layer]++;
            trkl_index++;
        }
    }

    //-------------------------------------------------------

}

//----------------------------------------------------------------------------------------



//________________________________________________________________________
void Ali_make_tracklets_from_digits::Terminate(Option_t *)
{
    cout << "In terminate" << endl;
}


//________________________________________________________________________
void Ali_make_tracklets_from_digits::FillHelix(AliESDtrack* track_in, Double_t magF_in)
{
    //-------------------
    // Get helix
    // Track parametrization:
    // https://www.physi.uni-heidelberg.de/~sma/alice/LukasLayer_bachelor.pdf
    Double_t alpha, alpha_alt, x_param, p_param[5]; // p are the track paramters, 0 = Y, 1 = Z, 2 = Snp, 3 = Tgl (p[2]/pt), 4 = Signed1Pt
    track_in->GetExternalParameters(x_param,p_param); // Those are the parameters used for vertexing
    alpha_alt=track_in->GetAlpha();

    //cout << "x_param: " << x_param << endl;
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    cout << "A i_param: " << i_param << ", param: " << p_param[i_param] << endl;
    //}

    //track_in->GetOuterExternalParameters(alpha,x_param,p_param); //
    //alpha_alt = alpha;
    Double_t fX     = track_in->GetX();

    //cout << "x_param: " << x_param << endl;
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    cout << "B i_param: " << i_param << ", param: " << p_param[i_param] << endl;
    //}

    //cout << "x_param: " << x_param << endl;

    //const AliExternalTrackParam* TPC_track_param = track_in->GetOuterParam();
    //const Double_t* p_param_b = TPC_track_param->GetParameter();
    //const Double_t* p_param_b = track_in->GetOuterParam()->GetParameter();
    //const Double_t* p_param_b = track_in->GetInnerParam()->GetParameter();
    //for(Int_t i_param = 0; i_param < 5; i_param++)
    //{
    //    p_param[i_param] = p_param_b[i_param];
    //}
    //x_param = fX;

    //-------------------
    // Correct way of filling aliHelix from
    // http://personalpages.to.infn.it/~puccio/htmldoc/src/AliHelix.cxx.html#PalT1E
    // line 52

    // CONCLUSION: AliTracker::GetBz() is not identical to GetC(magF_in), magnetic field is slightly different but it doesn't matter...

    Double_t fHelix_alt[9];
    Double_t x_alt,cs_alt,sn_alt;
    //track_in->GetExternalParameters(x_alt,fHelix_alt); // Those are the parameters used for vertexing
    x_alt = x_param;

    for(Int_t i_param = 0; i_param < 5; i_param++)
    {
	fHelix_alt[i_param] = p_param[i_param];
    }

    //cout << "alpha: " << alpha << ", alpha_alt: " << alpha_alt << endl;

    //
    //circle parameters
    //PH Sometimes fP4 and fHelix[4] are very big and the calculation
    //PH of the Sqrt cannot be done. To be investigated...

    // kB2C=-0.299792458e-3; // from /AliRoot/STEER/STEERBase/AliVParticle.h
    //Double_t kB2C_test =-0.299792458e-3;
    //Double_t GetC(Double_t b) const
    //{return fP[4]*b*kB2C;}

    //Double_t par4test = fHelix_alt[4]*magF_in*kB2C_test;
    fHelix_alt[4] = track_in->GetC(magF_in);
    //fHelix_alt[4] = par4test; // take the one with the magnetic field directly from the ESD file

    cs_alt = TMath::Cos(alpha_alt);
    sn_alt = TMath::Sin(alpha_alt);

    Double_t xc_alt, yc_alt, rc_alt;
    rc_alt  =  1/fHelix_alt[4];
    xc_alt  =  x_alt-fHelix_alt[2]*rc_alt;
    Double_t dummy = 1-(x_alt-xc_alt)*(x_alt-xc_alt)*fHelix_alt[4]*fHelix_alt[4];
    yc_alt  =  fHelix_alt[0]+TMath::Sqrt(dummy)/fHelix_alt[4];

    fHelix_alt[6] = xc_alt*cs_alt - yc_alt*sn_alt;
    fHelix_alt[7] = xc_alt*sn_alt + yc_alt*cs_alt;
    fHelix_alt[8] =  TMath::Abs(rc_alt);
    //
    //
    fHelix_alt[5]=x_alt*cs_alt - fHelix_alt[0]*sn_alt;            // x0
    fHelix_alt[0]=x_alt*sn_alt + fHelix_alt[0]*cs_alt;            // y0
    fHelix_alt[2]=TMath::ATan2(-(fHelix_alt[5]-fHelix_alt[6]),fHelix_alt[0]-fHelix_alt[7]); // phi0
    if (fHelix_alt[4]>0) fHelix_alt[2]-=TMath::Pi();
    fHelix_alt[5]   = fHelix_alt[6];
    fHelix_alt[0]   = fHelix_alt[7];

    for(Int_t i_param = 0; i_param < 9; i_param++)
    {
	aliHelix.fHelix[i_param] = fHelix_alt[i_param];
    }
    //-------------------
}



//________________________________________________________________________
void Ali_make_tracklets_from_digits::FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB)
{
    // V1.0
    Float_t pA[2] = {path_initA,path_initB}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    for(Int_t r = 0; r < 2; r++)
    {
	Double_t helix_point[3];
	helixA.Evaluate(pA[r],helix_point);
	testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[r]
	distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 30.0;
    while(fabs(scale_length) > 0.1 && loopcounter < 100) // stops when the length is too small
    {
	//cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
	//    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
	//    << ", d[1] = " << distarray[1] << ", flip = " << flip
	//    << ", scale_length = " << scale_length << endl;
	if(distarray[0] > distarray[1])
	{
	    if(loopcounter != 0)
	    {
		if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
		else scale = 0.7; // go on in this direction but only by the way * 0.7
	    }
	    scale_length = (pA[1]-pA[0])*scale; // the next length interval
	    pA[0]     = pA[1] + scale_length; // the new path

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[0],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[0]
	    distarray[0] = (testA-space_vec).Mag(); // new dca
	    flip = 1.0;
	}
	else
	{
	    if(loopcounter != 0)
	    {
		if(flip == -1.0) scale = 0.4;
		else scale = 0.7;
	    }
	    scale_length = (pA[0]-pA[1])*scale;
	    pA[1]     = pA[0] + scale_length;

	    Double_t helix_point[3];
	    helixA.Evaluate(pA[1],helix_point);
	    testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]); // 3D-vector of helixA at path pA[1]
	    distarray[1] = (testA-space_vec).Mag();
	    flip = -1.0;
	}
	loopcounter++;
    }

    if(loopcounter >= 100) cout << "WARNING: FindDCAHelixPoint exceeded maximum of 100 loops" << endl;

    if(distarray[0] < distarray[1])
    {
	pathA = pA[0];
	dcaAB = distarray[0];
    }
    else
    {
	pathA = pA[1];
	dcaAB = distarray[1];
    }
}


//________________________________________________________________________
Bool_t Ali_make_tracklets_from_digits::ReadDigits()
{
    //cout << "In ReadDigits" << endl;
    // don't do anything if the digits have already been loaded
    if (fDigitsLoadedFlag) return kTRUE;

    if(!fDigMan)
    {
	AliError("no digits manager");
	return kFALSE;
    }

    // reset digit arrays
    for(Int_t det = 0; det<  540; det++)
    {
	fDigMan->ClearArrays(det);
	fDigMan->ClearIndexes(det);
    }


    if(!fDigitsInputFile)
    {
	AliError("digits file not available");
	return kFALSE;
    }


    // read digits from file
    TTree* tr = (TTree*)fDigitsInputFile->Get(Form("Event%d/TreeD",fEventNoInFile));

    if(!tr)
    {
	//AliWarning(Form("digits tree for event %d not found", fEventNoInFile));
	return kFALSE;
    }

    //printf("  --> read");
    fDigMan->ReadDigits(tr);
    delete tr;

    //Int_t sum_dim_before = 0;
    //Int_t sum_dim_after  = 0;

    // expand digits for use in this task
    for(Int_t det = 0; det < 540; det++)
    {
	if(fDigMan->GetDigits(det))
        {
            //Int_t dim_before = fDigMan->GetDigits(det)->GetDim();
            fDigMan->GetDigits(det)->Expand();
            //Int_t dim_after = fDigMan->GetDigits(det)->GetDim();

            //sum_dim_before += dim_before;
            //sum_dim_after  += dim_after;
            //printf("det: %d, dim_before: %d, dim_after: %d \n",det,dim_before,dim_after);
	}
    }

    //cout << " --> size of fDigMan: " << sizeof(fDigMan) << endl;
    //printf("sum_dim_before: %d, sum_dim_after: %d \n",sum_dim_before,sum_dim_after);


    fDigitsLoadedFlag = kTRUE;
    return kTRUE;
}


