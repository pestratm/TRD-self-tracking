


//#include "Ana_Digits_functions.h"

#include "Ali_AS_Event.h"
#include "TBase_TRD_Calib.h"
#include "Ali_AS_EventLinkDef.h"

#include "TBase_TRD_Calib.h"
#include "TBase_TRD_CalibLinkDef.h"
ClassImp(TBase_TRD_Calib)


ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Event)


//---------------------------------------------------------------------------------
//class TGUI_TRD_Calib : public TGMainFrame, public TBase_TRD_Calib
class TGUI_TRD_Calib : public TGMainFrame
{
private:
    TRootEmbeddedCanvas *fCanvas_HV_vs_time        = NULL;
    TGMainFrame* Frame_Main;
    TGHorizontalFrame *hframe_Main[6];
    TGVerticalFrame   *vframe_Main[6];
    TGVerticalFrame   *vframe_stat_Main[6];

    TGNumberEntry*     arr_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_ana_params[4];
    TGLabel*           arr_Label_NEntry_stat[3];

    TGTextButton *Button_exit;
    TGTextButton *Button_load;
    TGTextButton *Button_save;
    TGTextButton *Button_draw3D;
    TGTextButton *Button_draw3D_track;
    TGTextButton *Button_draw3D_online_tracklets;
    TGTextButton *Button_draw3D_offline_tracklets;
    TGTextButton *Button_Calibrate;
    TGTextButton *Button_Track_Tracklets;
    TGTextButton *Button_draw2D_track;
    TGTextButton *Button_draw_digits;
    TGTextButton *Button_draw_TRD_tracks;
    TGTextButton *Button_draw_all_tracks;

    TBase_TRD_Calib *Base_TRD_Calib;
    vector<TEvePointSet*> vec_TPM3D_digits;
    vector<TEveLine*>   vec_TPL3D_tracklets;
    vector<TEveLine*>   vec_TPL3D_tracklets_match;
    TEvePointSet* TPM3D_cluster;
    TPolyMarker*  TPM_cluster;

    Long64_t N_Events;
    Long64_t N_Tracks;
    Long64_t N_Digits;

    vector< vector< vector<Float_t> > > vec_digit_track_info; // [track number][digit number][info ->] x,y,z,time,ADC,sector,stack,layer,row,column,dca
    vector< vector<Float_t> >           vec_track_info; // [track number][info ->] dca,TPCdEdx,momentum,eta_track,pT_track,TOFsignal,Track_length,TRDsumADC,TRD_signal,nsigma_TPC_e,nsigma_TPC_pi,nsigma_TPC_p

    TCanvas* c_3D        = NULL;
    TCanvas* c_3D_track  = NULL;

    Int_t color_layer[7] = {kGray,kGreen,kBlue,kMagenta,kCyan,kYellow,kOrange};
    vector<TEvePointSet*> vec_TPM3D_single_track_digit_layer;
    vector<TEvePointSet*> vec_TPM3D_single_track_digit;
    vector<TPolyMarker*> vec_TPM_single_track_digit_layer;
    TEvePointSet* TPM3D_single;

public:
    TGUI_TRD_Calib();
    virtual ~TGUI_TRD_Calib();
    Int_t LoadData();
    Int_t Draw3D();
    Int_t Draw3D_track();
    Int_t Draw_2D_track();
    void  Draw_digits();
    void  Draw_TRD_tracks();
    Int_t Draw_online_tracklets();
    Int_t Draw_offline_tracklets();
    Int_t Draw_all_tracks();
    Int_t Calibrate();
    Int_t Track_Tracklets();
    ClassDef(TGUI_TRD_Calib, 0)
};
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
TGUI_TRD_Calib::TGUI_TRD_Calib() : TGMainFrame(gClient->GetRoot(), 100, 100)
{
    //-------------------------------------
    cout << "TGUI_TRD_Calib started" << endl;
    TGaxis::SetMaxDigits(3);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    SetCleanup(kDeepCleanup);
    TGGC myGC = *gClient->GetResourcePool()->GetFrameGC();
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
    //-------------------------------------



    //-------------------------------------
    vec_TPM3D_single_track_digit_layer.resize(6); // layer
    vec_TPM_single_track_digit_layer.resize(6); // layer
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_single_track_digit_layer[i_layer] = new TEvePointSet();
        vec_TPM_single_track_digit_layer[i_layer]   = new TPolyMarker();
    }
    TPM3D_single  = new TEvePointSet();
    TPM3D_cluster = new TEvePointSet();
    TPM_cluster   = new TPolyMarker();
    //-------------------------------------



    //-------------------------------------
    // Create horizontal splitter
    Frame_Main = new TGMainFrame(gClient->GetRoot(), 400, 100);
    Frame_Main ->SetWindowName("Buttons");

    //--------------
    // A horizontal frame
    hframe_Main[0]  = new TGHorizontalFrame(Frame_Main,200,100);

    // exit button
    Button_exit = new TGTextButton(hframe_Main[0], "&Exit ","gApplication->Terminate(0)");
    hframe_Main[0]->AddFrame(Button_exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // load button
    Button_load = new TGTextButton(hframe_Main[0], "&Load ",10);
    Button_load->Connect("Clicked()", "TGUI_TRD_Calib", this, "LoadData()");
    hframe_Main[0]->AddFrame(Button_load, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // save button
    Button_save = new TGTextButton(hframe_Main[0], "&Save ","gApplication->Terminate(0)");
    hframe_Main[0]->AddFrame(Button_save, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw3D = new TGTextButton(hframe_Main[0], "&Draw 3D ",10);
    Button_draw3D->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw3D()");
    hframe_Main[0]->AddFrame(Button_draw3D, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw3D_track = new TGTextButton(hframe_Main[0], "&Draw 3D track ",10);
    Button_draw3D_track->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw3D_track()");
    hframe_Main[0]->AddFrame(Button_draw3D_track, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw3D_online_tracklets = new TGTextButton(hframe_Main[0], "&Draw on. trk. ",10);
    Button_draw3D_online_tracklets->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw_online_tracklets()");
    hframe_Main[0]->AddFrame(Button_draw3D_online_tracklets, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw3D_offline_tracklets = new TGTextButton(hframe_Main[0], "&Draw off. trk. ",10);
    Button_draw3D_offline_tracklets->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw_offline_tracklets()");
    hframe_Main[0]->AddFrame(Button_draw3D_offline_tracklets, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw_digits = new TGTextButton(hframe_Main[0], "&Draw digits ",10);
    Button_draw_digits->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw_digits()");
    hframe_Main[0]->AddFrame(Button_draw_digits, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw 3D button
    Button_draw_TRD_tracks = new TGTextButton(hframe_Main[0], "&Draw TRD tracks ",10);
    Button_draw_TRD_tracks->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw_TRD_tracks()");
    hframe_Main[0]->AddFrame(Button_draw_TRD_tracks, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    Frame_Main ->AddFrame(hframe_Main[0], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------

    //--------------
    // A horizontal frame
    hframe_Main[4]  = new TGHorizontalFrame(Frame_Main,200,100);

    // draw 3D button
    Button_draw_all_tracks = new TGTextButton(hframe_Main[4], "&Draw all tracks ",10);
    Button_draw_all_tracks->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw_all_tracks()");
    hframe_Main[4]->AddFrame(Button_draw_all_tracks, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    Frame_Main ->AddFrame(hframe_Main[4], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------


    //--------------
    TString arr_label_stat[3] = {"# events:          ","# tracks:       ","# digits:           "};
    hframe_Main[1]  = new TGHorizontalFrame(Frame_Main,200,100);
    for(Int_t i_param = 0; i_param < 3; i_param++)
    {
        vframe_stat_Main[i_param] = new TGVerticalFrame(hframe_Main[1], 400,200);
        TString label_entry = arr_label_stat[i_param];
        arr_Label_NEntry_stat[i_param] = new TGLabel(vframe_stat_Main[i_param], label_entry.Data(), myGC(), myfont->GetFontStruct());
        vframe_stat_Main[i_param]->AddFrame(arr_Label_NEntry_stat[i_param], new TGLayoutHints(kLHintsCenterX,20,20,2,2));
        hframe_Main[1]->AddFrame(vframe_stat_Main[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    }
    Frame_Main ->AddFrame(hframe_Main[1], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    //--------------
    TString arr_label_params_ana[4] = {"Event","Track","min pT","max pT"};
    hframe_Main[2]  = new TGHorizontalFrame(Frame_Main,200,100);
    Double_t min_max_pT_start[2] = {0.5,5.0};
    for(Int_t i_param = 0; i_param < 4; i_param++)
    {
        vframe_Main[i_param] = new TGVerticalFrame(hframe_Main[2], 200,200);
        if(i_param < 2)
        {
            arr_NEntry_ana_params[i_param] = new TGNumberEntry(vframe_Main[i_param], 0.0, 12,(TGNumberFormat::EStyle) 0);
            arr_NEntry_ana_params[i_param] ->SetNumStyle( TGNumberFormat::kNESInteger); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
        }
        else
        {
            arr_NEntry_ana_params[i_param] = new TGNumberEntry(vframe_Main[i_param], min_max_pT_start[i_param - 2], 12,(TGNumberFormat::EStyle) 1);
            arr_NEntry_ana_params[i_param] ->SetNumStyle( TGNumberFormat::kNESRealTwo); // https://root.cern.ch/doc/master/classTGNumberFormat.html#a8a0f81aac8ac12d0461aef554c6271ad
        }
        vframe_Main[i_param]->AddFrame(arr_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));

        TString label_entry = arr_label_params_ana[i_param];
        arr_Label_NEntry_ana_params[i_param] = new TGLabel(vframe_Main[i_param], label_entry.Data(), myGC(), myfont->GetFontStruct());
        vframe_Main[i_param]->AddFrame(arr_Label_NEntry_ana_params[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
        hframe_Main[2]->AddFrame(vframe_Main[i_param], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    }
    Frame_Main ->AddFrame(hframe_Main[2], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------


    //--------------
    // A horizontal frame
    hframe_Main[3]  = new TGHorizontalFrame(Frame_Main,200,100);

    // draw button
    Button_Calibrate = new TGTextButton(hframe_Main[3], "&Calibrate ",10);
    Button_Calibrate->Connect("Clicked()", "TGUI_TRD_Calib", this, "Calibrate()");
    hframe_Main[3]->AddFrame(Button_Calibrate, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw button
    Button_Track_Tracklets = new TGTextButton(hframe_Main[3], "&TrackTracklets ",10);
    Button_Track_Tracklets->Connect("Clicked()", "TGUI_TRD_Calib", this, "Track_Tracklets()");
    hframe_Main[3]->AddFrame(Button_Track_Tracklets, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

    // draw button
    Button_draw2D_track = new TGTextButton(hframe_Main[3], "&Draw 2D ",10);
    Button_draw2D_track->Connect("Clicked()", "TGUI_TRD_Calib", this, "Draw_2D_track()");
    hframe_Main[3]->AddFrame(Button_draw2D_track, new TGLayoutHints(kLHintsCenterX,5,5,3,4));


    Frame_Main ->AddFrame(hframe_Main[3], new TGLayoutHints(kLHintsCenterX,2,2,2,2));
    //--------------



    Frame_Main ->Resize(750,300); // size of frame
    Frame_Main ->MapSubwindows();
    Frame_Main ->MapWindow();
    Frame_Main ->Move(1050,50); // position of frame
    //-------------------------------------

    Base_TRD_Calib = new TBase_TRD_Calib();
    //Base_TRD_Calib ->Init_tree("list_tree_all_digits.txt");
    Base_TRD_Calib ->Init_tree("list_calib.txt");
    Base_TRD_Calib ->Loop_event(0);
    vec_TPM3D_digits    = Base_TRD_Calib ->get_PM3D_digits();
    vec_TPL3D_tracklets = Base_TRD_Calib ->get_PL3D_tracklets();
    vector<Int_t> vec_merge_time_bins;
#if 0
    vec_merge_time_bins.resize(4);
    vec_merge_time_bins[0] = 0;
    vec_merge_time_bins[1] = 5;
    vec_merge_time_bins[2] = 12;
    vec_merge_time_bins[3] = 23;
#endif
    vec_merge_time_bins.resize(24);
    for(Int_t i_time = 0; i_time < 24; i_time++)
    {
        vec_merge_time_bins[i_time] = i_time;
    }

    Base_TRD_Calib ->set_merged_time_bins(vec_merge_time_bins);

    LoadData();

}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
TGUI_TRD_Calib::~TGUI_TRD_Calib()
{
    // Clean up
    Cleanup();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::LoadData()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);

    Base_TRD_Calib ->Reset();

    Int_t Event = arr_NEntry_ana_params[0]->GetNumberEntry()->GetNumber();
    Base_TRD_Calib ->Loop_event(Event);
    N_Events = Base_TRD_Calib ->get_N_Events();
    N_Tracks = Base_TRD_Calib ->get_N_Tracks();
    N_Digits = Base_TRD_Calib ->get_N_Digits();
    vec_TPM3D_digits    = Base_TRD_Calib ->get_PM3D_digits();
    vec_TPL3D_tracklets = Base_TRD_Calib ->get_PL3D_tracklets();

    arr_Label_NEntry_stat[0] ->SetText(Form("# events: %lld",N_Events));
    arr_Label_NEntry_stat[1] ->SetText(Form("# tracks: %lld",N_Tracks));
    arr_Label_NEntry_stat[2] ->SetText(Form("# digits: %lld",N_Digits));

    vec_track_info       = Base_TRD_Calib ->get_track_info();
    vec_digit_track_info = Base_TRD_Calib ->get_digit_track_info();

#if 0
    for(Int_t i_track = 0; i_track < (Int_t)vec_digit_track_info.size(); i_track++)
    {
        Int_t n_digits_track = (Int_t)vec_digit_track_info[i_track].size();

        printf("TGUI_TRD_Calib::LoadData(), i_track: %d, n_digits_track: %d \n",i_track,n_digits_track);
    }
#endif

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw3D()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw3D->ChangeBackground(green);

    //if(!c_3D) c_3D    = new TCanvas("c_3D","c_3D",10,10,800,800);

    TGLViewer* TGL_viewer = Base_TRD_Calib ->Draw_TRD();
    TGL_viewer ->UpdateScene(kFALSE);
    TGL_viewer ->SetSmartRefresh(kFALSE);

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_digits[i_layer] ->SetMarkerColor(color_layer[i_layer]);
        vec_TPM3D_digits[i_layer] ->SetMarkerSize(0.5); // 5.0
        vec_TPM3D_digits[i_layer] ->SetMarkerStyle(20);
        vec_TPM3D_digits[i_layer] ->DrawClone("");
    }

#if 0
    Int_t N_tracklets_offline = (Int_t)vec_TPL3D_tracklets.size();
    for(Int_t i_tracklet = 0; i_tracklet < N_tracklets_offline; i_tracklet++)
    //for(Int_t i_tracklet = 0; i_tracklet < 5; i_tracklet++)
    {
        printf("TGUI_TRD_Calib::Draw3D(), i_tracklet: %d, out of %d \n",i_tracklet,N_tracklets_offline);
        vec_TPL3D_tracklets[i_tracklet] ->SetLineWidth(4);
        vec_TPL3D_tracklets[i_tracklet] ->SetLineColor(kRed);
        vec_TPL3D_tracklets[i_tracklet] ->SetLineStyle(1);
        vec_TPL3D_tracklets[i_tracklet] ->DrawClone("ogl");
    }
#endif

#if 0
    Int_t N_tracklets_offline_match = (Int_t)vec_TPL3D_tracklets_match.size();
    for(Int_t i_tracklet = 0; i_tracklet < N_tracklets_offline_match; i_tracklet++)
    //for(Int_t i_tracklet = 0; i_tracklet < 5; i_tracklet++)
    {
        printf("TGUI_TRD_Calib::Draw3D(), i_tracklet_match: %d, out of %d \n",i_tracklet,N_tracklets_offline_match);
        vec_TPL3D_tracklets_match[i_tracklet] ->SetLineWidth(6);
        vec_TPL3D_tracklets_match[i_tracklet] ->SetLineColor(kWhite);
        vec_TPL3D_tracklets_match[i_tracklet] ->SetLineStyle(1);
        vec_TPL3D_tracklets_match[i_tracklet] ->DrawClone("ogl");
    }
#endif

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Calibrate()
{
    printf("TGUI_TRD_Calib::Calibrate() \n");
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_Calibrate->ChangeBackground(green);

    Base_TRD_Calib ->Calibrate();

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Track_Tracklets()
{
    printf("TGUI_TRD_Calib::Track_Tracklets() \n");
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_Track_Tracklets->ChangeBackground(green);

    Base_TRD_Calib ->Track_Tracklets();
    vec_TPL3D_tracklets_match = Base_TRD_Calib ->get_PL3D_tracklets_match();

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
void TGUI_TRD_Calib::Draw_digits()
{
    printf("TGUI_TRD_Calib::Draw_digits() \n");
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw_digits->ChangeBackground(green);

    Base_TRD_Calib ->Draw_digits();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
void TGUI_TRD_Calib::Draw_TRD_tracks()
{
    printf("TGUI_TRD_Calib::Draw_TRD_tracks() \n");
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw_TRD_tracks->ChangeBackground(green);

    Base_TRD_Calib ->Draw_TRD_tracks();
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw_2D_track()
{
    printf("TGUI_TRD_Calib::Draw_2D_track() \n");
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw2D_track->ChangeBackground(green);

    Int_t i_track = arr_NEntry_ana_params[1]->GetNumberEntry()->GetNumber();

    vector<TVector2> vec_TV2_clusters;
    vec_TV2_clusters.resize(6);
    Int_t N_clusters_circle_fit = 0;
    vector< vector<TVector3> > vec_TV3_digit_pos_cluster = Base_TRD_Calib ->make_clusters(i_track); // layer, merged time bin
    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_TV3_digit_pos_cluster[i_layer].size(); i_time_merge++)
        {
            if(vec_TV3_digit_pos_cluster[i_layer][i_time_merge][0] > -999.0)
            {
                TPM_cluster ->SetNextPoint(vec_TV3_digit_pos_cluster[i_layer][i_time_merge][0],vec_TV3_digit_pos_cluster[i_layer][i_time_merge][1]);
                if(i_time_merge == 3 && i_layer < 3)
                {
                    vec_TV2_clusters[i_layer].SetX(vec_TV3_digit_pos_cluster[i_layer][i_time_merge][0]);
                    vec_TV2_clusters[i_layer].SetY(vec_TV3_digit_pos_cluster[i_layer][i_time_merge][1]);
                    N_clusters_circle_fit++;
                    //printf("i_layer: %d, i_time_merge: %d, N_clusters_circle_fit: %d, point: {%4.3f, %4.3f} \n",i_layer,i_time_merge,N_clusters_circle_fit,vec_TV3_digit_pos_cluster[i_layer][i_time_merge][0],vec_TV3_digit_pos_cluster[i_layer][i_time_merge][1]);
                }
            }
        }
    }

    Base_TRD_Calib ->Draw_2D_track(i_track);
    //printf("N_clusters_circle_fit: %d \n",N_clusters_circle_fit);
    if(N_clusters_circle_fit == 3) Base_TRD_Calib ->Draw_2D_circle_3points(vec_TV2_clusters);


    Float_t dca            = vec_track_info[i_track][0];
    Float_t TPCdEdx        = vec_track_info[i_track][1];
    Float_t momentum       = vec_track_info[i_track][2];
    Float_t eta_track      = vec_track_info[i_track][3];
    Float_t pT_track       = vec_track_info[i_track][4];
    Float_t TOFsignal      = vec_track_info[i_track][5];
    Float_t Track_length   = vec_track_info[i_track][6];
    Float_t TRDsumADC      = vec_track_info[i_track][7];
    Float_t TRD_signal     = vec_track_info[i_track][8];
    Float_t nsigma_TPC_e   = vec_track_info[i_track][9];
    Float_t nsigma_TPC_pi  = vec_track_info[i_track][10];
    Float_t nsigma_TPC_p   = vec_track_info[i_track][11];

    Int_t n_digits_track = (Int_t)vec_digit_track_info[i_track].size();

    for(Int_t i_digit = 0; i_digit < n_digits_track; i_digit++)
    {
        Float_t x_pos   = vec_digit_track_info[i_track][i_digit][0];
        Float_t y_pos   = vec_digit_track_info[i_track][i_digit][1];
        Float_t z_pos   = vec_digit_track_info[i_track][i_digit][2];
        Int_t time_bin  = (Int_t)vec_digit_track_info[i_track][i_digit][3];
        Float_t raw_ADC = vec_digit_track_info[i_track][i_digit][4];
        Int_t sector    = (Int_t)vec_digit_track_info[i_track][i_digit][5];
        Int_t stack     = (Int_t)vec_digit_track_info[i_track][i_digit][6];
        Int_t layer     = (Int_t)vec_digit_track_info[i_track][i_digit][7];
        Int_t row       = (Int_t)vec_digit_track_info[i_track][i_digit][8];
        Int_t column    = (Int_t)vec_digit_track_info[i_track][i_digit][9];
        Float_t dca     = vec_digit_track_info[i_track][i_digit][10];
        Float_t dca_x   = vec_digit_track_info[i_track][i_digit][11];
        Float_t dca_y   = vec_digit_track_info[i_track][i_digit][12];
        Float_t dca_z   = vec_digit_track_info[i_track][i_digit][13];

        Float_t dca_phi = TMath::Sqrt(dca_x*dca_x + dca_y*dca_y);
        //printf("dca_full: %4.3f, dca: {%4.3f, %4.3f, %4.3f}, dca_phi: %4.3f \n",dca,dca_x,dca_y,dca_z,dca_phi);

        Double_t digit_size = 1.0;

        if(dca_phi > 30.0) continue; // 3.0
        if(raw_ADC < 50.0) digit_size = 1.0;
        if(raw_ADC >= 50.0) digit_size = 10.0;

        vec_TPM_single_track_digit_layer[layer]   ->SetNextPoint(x_pos,y_pos);

        //printf("  ->2D track i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, row: %d, column: %d, time_bin: %d, sector: %d, stack: %d, layer: %d, raw_ADC: %4.3f \n",i_digit,x_pos,y_pos,z_pos,row,column,time_bin,sector,stack,layer,raw_ADC);
    }

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM_single_track_digit_layer[i_layer] ->SetMarkerColor(color_layer[i_layer]);
        vec_TPM_single_track_digit_layer[i_layer] ->SetMarkerSize(0.7);
        vec_TPM_single_track_digit_layer[i_layer] ->SetMarkerStyle(24);
        vec_TPM_single_track_digit_layer[i_layer] ->Draw("");
    }

    TPM_cluster ->SetMarkerColor(kRed);
    TPM_cluster ->SetMarkerSize(0.9);
    TPM_cluster ->SetMarkerStyle(20);
    TPM_cluster ->Draw("");

    Base_TRD_Calib ->Draw_tracklets_line_2D(i_track);

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw_online_tracklets()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw3D_online_tracklets->ChangeBackground(green);
    Base_TRD_Calib ->Draw_online_tracklets();

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw_offline_tracklets()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw3D_offline_tracklets->ChangeBackground(green);
    Base_TRD_Calib ->Draw_offline_tracklets();

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw_all_tracks()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw_all_tracks->ChangeBackground(green);

    Double_t min_pT = arr_NEntry_ana_params[2] ->GetNumberEntry()->GetNumber();
    Double_t max_pT = arr_NEntry_ana_params[3] ->GetNumberEntry()->GetNumber();

    Base_TRD_Calib ->Draw_all_tracks(min_pT,max_pT);

    gEve->Redraw3D(kTRUE);

    return 1;
}
//---------------------------------------------------------------------------------



//---------------------------------------------------------------------------------
Int_t TGUI_TRD_Calib::Draw3D_track()
{
    Pixel_t green;
    gClient->GetColorByName("green", green);
    Button_draw3D_track->ChangeBackground(green);


    //if(!c_3D_track) c_3D_track    = new TCanvas("c_3D_track","c_3D_track",10,10,800,800);

    Int_t i_track = arr_NEntry_ana_params[1]->GetNumberEntry()->GetNumber();

    //Base_TRD_Calib ->Draw_TRD();
    Base_TRD_Calib ->Draw_track(i_track);
    Base_TRD_Calib ->Draw_neighbor_tracks(i_track);


    vector< vector<TVector3> > vec_TV3_digit_pos_cluster = Base_TRD_Calib ->make_clusters(i_track); // layer, merged time bin

    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        for(Int_t i_time_merge = 0; i_time_merge < (Int_t)vec_TV3_digit_pos_cluster[i_layer].size(); i_time_merge++)
        {
            TPM3D_cluster ->SetNextPoint(vec_TV3_digit_pos_cluster[i_layer][i_time_merge][0],vec_TV3_digit_pos_cluster[i_layer][i_time_merge][1],vec_TV3_digit_pos_cluster[i_layer][i_time_merge][2]);
        }
    }

    gEve->AddElement(TPM3D_cluster);


    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        vec_TPM3D_single_track_digit_layer[i_layer] ->Clear();
    }

    vec_TPM3D_single_track_digit.clear();


    Float_t dca            = vec_track_info[i_track][0];
    Float_t TPCdEdx        = vec_track_info[i_track][1];
    Float_t momentum       = vec_track_info[i_track][2];
    Float_t eta_track      = vec_track_info[i_track][3];
    Float_t pT_track       = vec_track_info[i_track][4];
    Float_t TOFsignal      = vec_track_info[i_track][5];
    Float_t Track_length   = vec_track_info[i_track][6];
    Float_t TRDsumADC      = vec_track_info[i_track][7];
    Float_t TRD_signal     = vec_track_info[i_track][8];
    Float_t nsigma_TPC_e   = vec_track_info[i_track][9];
    Float_t nsigma_TPC_pi  = vec_track_info[i_track][10];
    Float_t nsigma_TPC_p   = vec_track_info[i_track][11];

    Int_t n_digits_track = (Int_t)vec_digit_track_info[i_track].size();

    printf("TGUI_TRD_Calib::Draw3D_track(), i_track: %d, n_digits_track: %d \n",i_track,n_digits_track);


    for(Int_t i_digit = 0; i_digit < n_digits_track; i_digit++)
    {
        Float_t x_pos   = vec_digit_track_info[i_track][i_digit][0];
        Float_t y_pos   = vec_digit_track_info[i_track][i_digit][1];
        Float_t z_pos   = vec_digit_track_info[i_track][i_digit][2];
        Int_t time_bin  = (Int_t)vec_digit_track_info[i_track][i_digit][3];
        Float_t raw_ADC = vec_digit_track_info[i_track][i_digit][4];
        Int_t sector    = (Int_t)vec_digit_track_info[i_track][i_digit][5];
        Int_t stack     = (Int_t)vec_digit_track_info[i_track][i_digit][6];
        Int_t layer     = (Int_t)vec_digit_track_info[i_track][i_digit][7];
        Int_t row       = (Int_t)vec_digit_track_info[i_track][i_digit][8];
        Int_t column    = (Int_t)vec_digit_track_info[i_track][i_digit][9];
        Float_t dca     = vec_digit_track_info[i_track][i_digit][10];
        Float_t dca_x   = vec_digit_track_info[i_track][i_digit][11];
        Float_t dca_y   = vec_digit_track_info[i_track][i_digit][12];
        Float_t dca_z   = vec_digit_track_info[i_track][i_digit][13];

        Float_t dca_phi = TMath::Sqrt(dca_x*dca_x + dca_y*dca_y);
        //printf("dca_full: %4.3f, dca: {%4.3f, %4.3f, %4.3f}, dca_phi: %4.3f \n",dca,dca_x,dca_y,dca_z,dca_phi);

        Double_t digit_size = 1.0;

        //if(raw_ADC < 20.0) continue;
        if(dca_phi > 30.0) continue; // 3.0

        if(raw_ADC < 50.0) digit_size = 1.0;
        if(raw_ADC >= 50.0) digit_size = 10.0;
        //if(raw_ADC >= 50.0 && raw_ADC < 100.0) digit_size = 2.0;
        //if(raw_ADC >= 100.0 && raw_ADC < 150.0) digit_size = 3.0;
        //if(raw_ADC >= 150.0) digit_size = 4.0;

        //if(!(row == 14 && column == 120)) continue;
        TPM3D_single ->SetNextPoint(x_pos,y_pos,z_pos);
        vec_TPM3D_single_track_digit_layer[layer]   ->SetNextPoint(x_pos,y_pos,z_pos);
        vec_TPM3D_single_track_digit.push_back((TEvePointSet*)TPM3D_single ->Clone());
        Int_t N_digits = (Int_t)vec_TPM3D_single_track_digit.size();
        vec_TPM3D_single_track_digit[N_digits - 1] ->SetMarkerSize(digit_size);
        vec_TPM3D_single_track_digit[N_digits - 1] ->SetMarkerColor(kGreen+2);
        vec_TPM3D_single_track_digit[N_digits - 1] ->SetMarkerStyle(20);

        printf("  -> i_digit: %d, pos: {%4.3f, %4.3f, %4.3f}, row: %d, column: %d, time_bin: %d, sector: %d, stack: %d, layer: %d, raw_ADC: %4.3f \n",i_digit,x_pos,y_pos,z_pos,row,column,time_bin,sector,stack,layer,raw_ADC);
    }


    for(Int_t i_layer = 0; i_layer < 6; i_layer++)
    {
        //for (Int_t i_digit_per_layer = 0; i_digit_per_layer< vec_TPM3D_single_track_digit_layer[i_layer]->GetN(); i_digit_per_layer++)
        //{
        //    vec_TPM3D_single_track_digit_layer[i_layer][i_digit_per_layer].SetMarkerColor(color_layer[i_layer]-i_digit_per_layer);
        //}

        //cout << "i_layer : " << i_layer << endl;
        //cout << "vec_TPM3D_single_track_digit_layer[i_layer]->Size() : " << vec_TPM3D_single_track_digit_layer[i_layer]->Size() << endl;

        vec_TPM3D_single_track_digit_layer[i_layer] ->SetMarkerColor(color_layer[i_layer]);
        vec_TPM3D_single_track_digit_layer[i_layer] ->SetMarkerSize(1.0);
        vec_TPM3D_single_track_digit_layer[i_layer] ->SetMarkerStyle(20);
        //vec_TPM3D_single_track_digit_layer[i_layer] ->DrawClone("");
        gEve->AddElement(vec_TPM3D_single_track_digit_layer[i_layer]);
    }


    TPM3D_cluster ->SetMarkerColor(kRed);
    TPM3D_cluster ->SetMarkerSize(1.5);
    TPM3D_cluster ->SetMarkerStyle(20);
    gEve->AddElement(TPM3D_cluster);
    //TPM3D_cluster ->DrawClone("");


#if 0
    for(Int_t i_digit = 0; i_digit < (Int_t)vec_TPM3D_single_track_digit.size(); i_digit++)
    {
        printf("i_digit: %d \n",i_digit);
        vec_TPM3D_single_track_digit[i_digit] ->DrawClone("");
    }
#endif


    //Base_TRD_Calib ->Draw_line(i_track);
    Base_TRD_Calib ->Draw_tracklets_line(i_track);
    //Base_TRD_Calib ->make_plots_ADC(i_track);


    gEve->Redraw3D(kTRUE);
    return 1;
}
//---------------------------------------------------------------------------------


