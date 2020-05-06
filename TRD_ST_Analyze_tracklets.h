
#include "Ali_TRD_ST.h"
#include "Ali_TRD_ST_LinkDef.h"



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

    TString HistName;

    TEveLine* TEveLine_beam_axis = NULL;
    TEveLine* TPL3D_helix = NULL;
    vector<TEveLine*> vec_TPL3D_helix;
    vector< vector<TEveLine*> > vec_TEveLine_tracklets;
    Int_t N_tracklets_layers[6] = {0};
    Int_t color_layer[6] = {kBlue,kGreen,kOrange,kCyan,kMagenta,kYellow};

public:
    Ali_TRD_ST_Analyze();
    ~Ali_TRD_ST_Analyze();
    void Init_tree(TString SEList);
    void Loop_event(Long64_t i_event);
    void Draw_event(Long64_t i_event);

    ClassDef(Ali_TRD_ST_Analyze, 1)
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Ali_TRD_ST_Analyze::Ali_TRD_ST_Analyze()
{
    // constructor
    TEveManager::Create();

    TPL3D_helix = new TEveLine();
    TEveLine_beam_axis = new TEveLine();
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,-300.0);
    TEveLine_beam_axis ->SetNextPoint(0.0,0.0,300.0);
    TEveLine_beam_axis ->SetName("beam axis");
    TEveLine_beam_axis ->SetLineStyle(1);
    TEveLine_beam_axis ->SetLineWidth(4);
    TEveLine_beam_axis ->SetMainColor(kBlue);
    gEve->AddElement(TEveLine_beam_axis);

    vec_TEveLine_tracklets.resize(6); // layers
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Init_tree(TString SEList)
{
    printf("Ali_TRD_ST_Analyze::Init_tree \n");
    TString pinputdir = "./Data/";

    TRD_ST_Tracklet   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event      = new Ali_TRD_ST_Event();

    // Same event input
    if (!SEList.IsNull())   // if input file is ok
    {
        cout << "Open same event file list " << SEList << endl;
        ifstream in(SEList);  // input stream
        if(in)
        {
            cout << "file list is ok" << endl;
            input_SE  = new TChain( TRD_ST_TREE.Data(), TRD_ST_TREE.Data() );
            char str[255];       // char array for each file name
            Long64_t entries_save = 0;
            while(in)
            {
                in.getline(str,255);  // take the lines of the file list
                if(str[0] != 0)
                {
                    TString addfile;
                    addfile = str;
                    addfile = pinputdir+addfile;
                    input_SE ->AddFile(addfile.Data(),-1, TRD_ST_TREE.Data() );
                    Long64_t file_entries = input_SE->GetEntries();
                    cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
                    entries_save = file_entries;
                }
            }
            input_SE  ->SetBranchAddress( TRD_ST_BRANCH, &TRD_ST_Event );
        }
        else
        {
            cout << "WARNING: SE file input is problemtic" << endl;
        }
    }

    file_entries_total = input_SE->GetEntries();
    N_Events = file_entries_total;
    cout << "Total number of events in tree: " << file_entries_total << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Loop_event(Long64_t i_event)
{
    printf("Ali_TRD_ST_Analyze::Loop_event \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event


    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    Double_t EventVertexX         = TRD_ST_Event ->getx();
    Double_t EventVertexY         = TRD_ST_Event ->gety();
    Double_t EventVertexZ         = TRD_ST_Event ->getz();
    Float_t  V0MEq                = TRD_ST_Event ->getcent_class_V0MEq();
    //--------------------------------------------------



    //--------------------------------------------------
    // TPC track loop
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

        Double_t nsigma_TPC_e   = TRD_ST_TPC_Track ->getnsigma_e_TPC();
        Double_t nsigma_TPC_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TPC();
        Double_t nsigma_TPC_p   = TRD_ST_TPC_Track ->getnsigma_p_TPC();
        Double_t nsigma_TOF_e   = TRD_ST_TPC_Track ->getnsigma_e_TOF();
        Double_t nsigma_TOF_pi  = TRD_ST_TPC_Track ->getnsigma_pi_TOF();
        Double_t TRD_signal     = TRD_ST_TPC_Track ->getTRDSignal();
        Double_t TRDsumADC      = TRD_ST_TPC_Track ->getTRDsumADC();
        Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex
        TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
        UShort_t NTPCcls        = TRD_ST_TPC_Track ->getNTPCcls();
        UShort_t NTRDcls        = TRD_ST_TPC_Track ->getNTRDcls();
        UShort_t NITScls        = TRD_ST_TPC_Track ->getNITScls();
        Float_t TPCchi2         = TRD_ST_TPC_Track ->getTPCchi2();
        Float_t TPCdEdx         = TRD_ST_TPC_Track ->getTPCdEdx();
        Float_t TOFsignal       = TRD_ST_TPC_Track ->getTOFsignal(); // in ps (1E-12 s)
        Float_t Track_length    = TRD_ST_TPC_Track ->getTrack_length();

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();
    }
    //--------------------------------------------------



    //--------------------------------------------------
    // TRD tracklet loop
    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t    i_det;

    for(Int_t i_tracklet = 0; i_tracklet < NumTracklets; i_tracklet++)
    {
        TRD_ST_Tracklet = TRD_ST_Event    ->getTracklet(i_tracklet);
        TV3_offset      = TRD_ST_Tracklet ->get_TV3_offset();
        TV3_dir         = TRD_ST_Tracklet ->get_TV3_dir();
        i_det           = TRD_ST_Tracklet ->get_TRD_det();

        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

    }
    //--------------------------------------------------



}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_event(Long64_t i_event)
{
    printf("Ali_TRD_ST_Analyze::Draw_event \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event


    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    Double_t EventVertexX         = TRD_ST_Event ->getx();
    Double_t EventVertexY         = TRD_ST_Event ->gety();
    Double_t EventVertexZ         = TRD_ST_Event ->getz();
    //--------------------------------------------------



    //--------------------------------------------------
    // TPC track loop
    Double_t track_pos[3];
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
    {
        TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

        vec_TPL3D_helix.resize(i_track+1);
        vec_TPL3D_helix[i_track] = new TEveLine();

        for(Double_t track_path = 0.0; track_path < 300; track_path += 1.0)
        {
            TRD_ST_TPC_Track ->Evaluate(track_path,track_pos);
            vec_TPL3D_helix[i_track] ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);

            //if(i_track == 0) printf("track_path: %4.3f, pos: {%4.2f, %4.2f, %4.2f} \n",track_path,track_pos[0],track_pos[1],track_pos[2]);
        }

        HistName = "track ";
        HistName += i_track;
        vec_TPL3D_helix[i_track]    ->SetName(HistName.Data());
        vec_TPL3D_helix[i_track]    ->SetLineStyle(1);
        vec_TPL3D_helix[i_track]    ->SetLineWidth(2);
        vec_TPL3D_helix[i_track]    ->SetMainColor(kRed);
        gEve->AddElement(vec_TPL3D_helix[i_track]);
    }

    //--------------------------------------------------



    //--------------------------------------------------
    // TRD tracklet loop
    Double_t scale_length = 3.0;

    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t    i_det;

    for(Int_t i_tracklet = 0; i_tracklet < NumTracklets; i_tracklet++)
    {
        TRD_ST_Tracklet = TRD_ST_Event    ->getTracklet(i_tracklet);
        TV3_offset      = TRD_ST_Tracklet ->get_TV3_offset();
        TV3_dir         = TRD_ST_Tracklet ->get_TV3_dir();
        i_det           = TRD_ST_Tracklet ->get_TRD_det();

        if(TV3_offset.Mag() > 1000.0) continue;

        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        vec_TEveLine_tracklets[i_layer].resize(N_tracklets_layers[i_layer]+1);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] = new TEveLine();

        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0],TV3_offset[1],TV3_offset[2]);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0] + scale_length*TV3_dir[0],TV3_offset[1] + scale_length*TV3_dir[1],TV3_offset[2] + scale_length*TV3_dir[2]);

        Double_t radius = TMath::Sqrt( TMath::Power(TV3_offset[0],2) + TMath::Power(TV3_offset[1],2) );
        printf("i_tracklet: %d, radius: %4.3f, pos A: {%4.2f, %4.2f, %4.2f}, pos B: {%4.2f, %4.2f, %4.2f} \n",i_tracklet,radius,TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + scale_length*TV3_dir[0],TV3_offset[1] + scale_length*TV3_dir[1],TV3_offset[2] + scale_length*TV3_dir[2]);

        HistName = "tracklet ";
        HistName += i_tracklet;
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetName(HistName.Data());
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetLineStyle(1);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetLineWidth(2);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetMainColor(color_layer[i_layer]);
        gEve->AddElement(vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]);

        N_tracklets_layers[i_layer]++;
    }

    gEve->Redraw3D(kTRUE);
    //--------------------------------------------------
}
//----------------------------------------------------------------------------------------

