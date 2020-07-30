
#include "./functions_BW.h"

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

Double_t ExpFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1;
    par0  = fabs(par[0]);
    par1  = par[1];
    x = x_val[0];
    y = par0*TMath::Exp(x*par1);
    return y;
}
//----------------------------------------------------------------------------------------



void Ana_sec_vertices()
{



    //--------------------------
    printf("Ana_sec_vertices started \n");

    //gStyle->SetPalette(kDarkBodyRadiator); // https://root.cern.ch/doc/master/classTColor.html
    //gStyle->SetPalette(kInvertedDarkBodyRadiator); // https://root.cern.ch/doc/master/classTColor.html
    gStyle->SetPalette(kTemperatureMap); // https://root.cern.ch/doc/master/classTColor.html
    //--------------------------



    //--------------------------
    // Load TRD geometry
    TEveManager::Create();

    TH1I* h_good_bad_TRD_chambers;
    h_good_bad_TRD_chambers = new TH1I("h_good_bad_TRD_chambers","h_good_bad_TRD_chambers",540,0,540);
    TFile* file_TRD_QA_flags = TFile::Open("./Data/chamber_QC_flags.root");
    vector<int> *t_flags;
    file_TRD_QA_flags ->GetObject("QC_flags", t_flags);

    // Its a 3 digit binary number. LSB is ADC good = 0 or bad = 1, next bit is anode HV good = 0, or bad = 1, and last bit is drift HV
    // so a 3 means that the ADC and the anode HV was bad, but the drift HV was okay

    // LSB = official QA, bit 1 = no fit, bit 2 = anode HV defect, bit 3 = drift HV defect, bit 4 = adc defect

    // number   adc defect   drift HV defect   anode HD defect    no fit   official QA
    //   0          0               0                0               0          0         --> all good
    //   1          0               0                0               0          1         --> official QA bad, rest good
    //  ...
    //   31         1               1                1               1          1         --> all bad

    Int_t i_chamber = 0;
    for(vector<int>::iterator it = t_flags->begin(); it != t_flags->end(); ++it)
    {
        //cout << "chamber: " << i_chamber << ", it: "  << *it << ", " << t_flags->at(i_chamber) << endl;
        h_good_bad_TRD_chambers ->SetBinContent(i_chamber+1,t_flags->at(i_chamber));
        i_chamber++;
    }

    Int_t color_flag_QC[32];
    for(Int_t i_QC_flag = 0; i_QC_flag < 32; i_QC_flag++)
    {
        color_flag_QC[i_QC_flag] = kCyan;

        Int_t k_bit = 1; // fit
        Int_t bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // no fit
        {
            color_flag_QC[i_QC_flag] = kPink;
        }

        k_bit = 4; // ADC value
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // ADC low
        {
            color_flag_QC[i_QC_flag] = kMagenta;
        }

        k_bit = 2; // anode HV
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // anode HV low
        {
            color_flag_QC[i_QC_flag] = kYellow;
        }

        k_bit = 3; // drift HV bit
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // drift HV defect
        {
            color_flag_QC[i_QC_flag] = kOrange;
        }

        k_bit = 0; // official QA
        bit_value = (i_QC_flag & ( 1 << k_bit )) >> k_bit;
        if(bit_value == 1) // official QA bad
        {
            color_flag_QC[i_QC_flag] = kRed;
        }
    }
    color_flag_QC[31] = kRed;

    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
    TFile* file_TRD_geom = TFile::Open("./Data/TRD_Geom.root");
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

    TH2D* h2D_dEdx_vs_mom = new TH2D("h2D_dEdx_vs_mom","h2D_dEdx_vs_mom",200,0,5,200,0,150.0);


    vector<TEveBox*> vec_eve_TRD_detector_box;
    vec_eve_TRD_detector_box.resize(540);
    vector< vector<TPolyLine*> > vec_PL_TRD_det_2D;
    vec_PL_TRD_det_2D.resize(18); // sectors
    for(Int_t i_sector = 0; i_sector < 18; i_sector++)
    {
        vec_PL_TRD_det_2D[i_sector].resize(6); // layers
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            vec_PL_TRD_det_2D[i_sector][i_layer] = new TPolyLine();
        }
    }
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        Int_t i_sector = (Int_t)(TRD_detector/30);
        Int_t i_stack  = (Int_t)(TRD_detector%30/6);
        Int_t i_layer  = TRD_detector%6;
        //Int_t i_det = layer + 6*stack + 30*sector;


        vec_eve_TRD_detector_box[TRD_detector] = new TEveBox;

        HistName = "TRD_box_";
        HistName += TRD_detector;
        vec_eve_TRD_detector_box[TRD_detector] ->SetName(HistName.Data());
        Int_t flag_QC = h_good_bad_TRD_chambers ->GetBinContent(TRD_detector+1);
        if(!flag_QC) // chamber is OK flagged by QA
        {
            vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(kCyan);
            vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(95); // the higher the value the more transparent
        }
        else // bad chamber
        {
            vec_eve_TRD_detector_box[TRD_detector]->SetMainColor(color_flag_QC[flag_QC]);
            vec_eve_TRD_detector_box[TRD_detector]->SetMainTransparency(85); // the higher the value the more transparent
        }

        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            Double_t arr_pos_glb[3] = {vec_TH1D_TRD_geometry[0][i_vertex]->GetBinContent(TRD_detector),vec_TH1D_TRD_geometry[1][i_vertex]->GetBinContent(TRD_detector),vec_TH1D_TRD_geometry[2][i_vertex]->GetBinContent(TRD_detector)};
            vec_eve_TRD_detector_box[TRD_detector]->SetVertex(i_vertex,arr_pos_glb[0],arr_pos_glb[1],arr_pos_glb[2]);

            if(i_stack == 0 && i_vertex < 4)
            {
                vec_PL_TRD_det_2D[i_sector][i_layer] ->SetPoint(i_vertex,arr_pos_glb[0],arr_pos_glb[1]);
                if(i_vertex == 0) vec_PL_TRD_det_2D[i_sector][i_layer] ->SetPoint(4,arr_pos_glb[0],arr_pos_glb[1]);
            }
        }

        gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
    }
    gEve->Redraw3D(kTRUE);
    //--------------------------

    //-------------------------
    Double_t TRD_layer_radii[6][2] =
    {
        {297.5,306.5},
        {310.0,320.0},
        {323.0,333.0},
        {336.0,345.5},
        {348.0,357.0},
        {361.0,371.0}
    };

    //Int_t color_layer_match[6] = {kRed,kGreen,kSpring-6,kYellow,kPink-3,kOrange+8};
    Int_t color_layer_match[6] = {kGray+1,kAzure-2,kGreen+2,kCyan+2,kOrange+2,kRed};

    //--------------------------

    //--------------------------
    // Open Ntuple
    //TFile* inputfile = TFile::Open("./ST_out/Merge_ST_hoppner_V2.root");
    TFile* inputfile = TFile::Open("./ST_out/Merge_ST_ABCDE_V9.root");
    TNtuple* NT_sec_vertices = (TNtuple*)inputfile->Get("NT_secondary_vertices");
    Float_t x_sec,y_sec,z_sec,bitmap_sec,pT_AB_sec;
    Float_t qpT_A_sec,qpT_B_sec,AP_pT_sec,AP_alpha_sec,dcaTPC_sec,pathTPC_sec,InvM_sec;
    NT_sec_vertices ->SetBranchAddress("x",&x_sec);
    NT_sec_vertices ->SetBranchAddress("y",&y_sec);
    NT_sec_vertices ->SetBranchAddress("z",&z_sec);
    NT_sec_vertices ->SetBranchAddress("ntracks",&bitmap_sec); // bitmap for shared clusters
    NT_sec_vertices ->SetBranchAddress("pT_AB",&pT_AB_sec);
    NT_sec_vertices ->SetBranchAddress("qpT_A",&qpT_A_sec);
    NT_sec_vertices ->SetBranchAddress("qpT_B",&qpT_B_sec);
    NT_sec_vertices ->SetBranchAddress("AP_pT",&AP_pT_sec);
    NT_sec_vertices ->SetBranchAddress("AP_alpha",&AP_alpha_sec);
    NT_sec_vertices ->SetBranchAddress("dcaTPC",&dcaTPC_sec);
    NT_sec_vertices ->SetBranchAddress("pathTPC",&pathTPC_sec);
    NT_sec_vertices ->SetBranchAddress("InvM",&InvM_sec);  // of photon candidates

    TNtuple* NT_sec_cluster = (TNtuple*)inputfile->Get("NT_secondary_vertex_cluster");
    Float_t x_clus,y_clus,z_clus,ntracks_clus,dcaTPC_clus,tof,trklength,dEdx,dcaprim,pT,mom;
    NT_sec_cluster ->SetBranchAddress("x",&x_clus);
    NT_sec_cluster ->SetBranchAddress("y",&y_clus);
    NT_sec_cluster ->SetBranchAddress("z",&z_clus);
    NT_sec_cluster ->SetBranchAddress("nvertices",&ntracks_clus);
    NT_sec_cluster ->SetBranchAddress("dcaTPC",&dcaTPC_clus);
    NT_sec_cluster ->SetBranchAddress("tof",&tof);
    NT_sec_cluster ->SetBranchAddress("trklength",&trklength);
    NT_sec_cluster ->SetBranchAddress("dEdx",&dEdx);
    NT_sec_cluster ->SetBranchAddress("dcaprim",&dcaprim);
    NT_sec_cluster ->SetBranchAddress("pT",&pT);
    NT_sec_cluster ->SetBranchAddress("mom",&mom);

    //--------------------------

    Int_t bin_nbr=200;
    TH2D* h2d_vertex_pos_xy = new TH2D("h2d_vertex_pos_xy","h2d_vertex_pos_xy",500,-400,400,500,-400,400);
    vector<TH1D*> h1d_vertex_pos_r;
    h1d_vertex_pos_r.resize(6);
    for (Int_t i_a=0; i_a<6; i_a++) h1d_vertex_pos_r[i_a] = new TH1D("h1d_vertex_pos_r","h1d_vertex_pos_r",bin_nbr,200,400);
    TH1D* h1d_vertex_mom_pT = new TH1D("h1d_vertex_mom_pT","h1d_vertex_mom_pT",bin_nbr,0,1);
    TH1D* h1d_vertex_mom_pT_exp = new TH1D("h1d_vertex_mom_pT_exp","h1d_vertex_mom_pT_exp",bin_nbr,0,1);
    TEvePointSet* TEveP_offset_points = new TEvePointSet();

    TH2D* h2d_cluster_pos_xy = new TH2D("h2d_cluster_pos_xy","h2d_cluster_pos_xy",500,-400,400,500,-400,400);
    TH1D* h1d_cluster_pos_r = new TH1D("h1d_cluster_pos_r","h1d_cluster_pos_r",200,250,450);
    TH1D** h1d_cluster_pos_r_array = new TH1D*[5];
    TString HistName;
    for(Int_t i_clus_arr=0; i_clus_arr<5;i_clus_arr++)
    {
        HistName="h1d_cluster_pos_r_";
        HistName+=i_clus_arr;
        h1d_cluster_pos_r_array[i_clus_arr] = new TH1D(HistName,HistName,200,270,380);
    }
    //--------------------------
    // Loop over photon conversion data
    Long64_t N_entries = NT_sec_vertices->GetEntries();
    for(Long64_t i_entry = 0; i_entry < N_entries; i_entry++)
        //for(Long64_t i_entry = 0; i_entry < 50; i_entry++)
    {
        if (i_entry != 0  &&  i_entry % 50 == 0)
            cout << "." << flush;
        if (i_entry != 0  &&  i_entry % 500 == 0)
        {
            printf("i_entry: %lld out of %lld, %4.2f%% total done \n",i_entry,N_entries,((Double_t)i_entry/(Double_t)N_entries)*100.0);
        }

        NT_sec_vertices ->GetEntry(i_entry);

        // Test the bit map
        Int_t independent_layer_A[6] = {0};
        Int_t independent_layer_B[6] = {0};
        Int_t shared_layer[6]        = {0};
        for(Int_t i_bit = 0; i_bit < 18; i_bit++)
        {
            if(((Int_t)bitmap_sec >> i_bit) & 1)
            {
                if(i_bit < 6)                independent_layer_A[i_bit]   = 1;
                if(i_bit >= 6 && i_bit < 12) independent_layer_B[i_bit-6] = 1;
                if(i_bit >= 12)              shared_layer[i_bit-12]       = 1;
            }
        }
        vector<Bool_t> conds;
        Bool_t cond1=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 1 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3]) > 1 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3]) > 1 );
        conds.push_back(cond1);
        Bool_t cond2=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2] + shared_layer[3]) > 2 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3]) > 0 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3]) > 0 );
        conds.push_back(cond2);

        Bool_t cond3=(
                      (shared_layer[0] + shared_layer[1] + shared_layer[2]) > 0 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3] + independent_layer_A[2]) > 2 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3] + independent_layer_B[2]) > 2 );
        conds.push_back(cond3);

        Bool_t cond4=(
                      (shared_layer[0] + shared_layer[1] ) > 0 &&
                      (independent_layer_A[5] + independent_layer_A[4] + independent_layer_A[3] + independent_layer_A[2]+ independent_layer_A[1]) > 2 &&
                      (independent_layer_B[5] + independent_layer_B[4] + independent_layer_B[3] + independent_layer_B[2]+ independent_layer_B[1]) > 2 );
        conds.push_back(cond4);

        Bool_t cond5=!(
                       (shared_layer[3] + shared_layer[4] + shared_layer[5]) > 0 &&
                       (independent_layer_A[0] + independent_layer_A[1] + independent_layer_A[2] ) > 0 &&
                       (independent_layer_B[0] + independent_layer_B[1] + independent_layer_B[2] ) > 0 );
        conds.push_back(cond5);

        if((cond1 || cond2 || cond3 ||cond4) && cond5)
        {
            if(dcaTPC_sec > 10.0 || (dcaTPC_sec <= 10.0 && pathTPC_sec < 0.0)) // no close by TPC track
            {
                h2d_vertex_pos_xy ->Fill(x_sec,y_sec);
                Double_t rad=TMath::Sqrt(x_sec*x_sec +y_sec*y_sec);
                h1d_vertex_pos_r[0]  ->Fill(rad);
                for(Int_t i_conds=0; i_conds<(Int_t)conds.size(); i_conds++)
                    if (conds[i_conds]) h1d_vertex_pos_r[i_conds+1]  ->Fill(rad);

                h1d_vertex_mom_pT ->Fill(pT_AB_sec);

                TEveP_offset_points  ->SetPoint(i_entry,x_sec,y_sec,z_sec);
            }
        }

        //printf("i_entry: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_entry,x_sec,y_sec,z_sec);
    }

    for(Int_t i_bin=0; i_bin<bin_nbr; i_bin++)
    {
        Double_t binContent = h1d_vertex_mom_pT->GetBinContent(i_bin);
        Double_t binCenter = h1d_vertex_mom_pT->GetBinCenter(i_bin);
        h1d_vertex_mom_pT_exp->SetBinContent(i_bin,binContent/(binCenter));

    }
    //--------------------------



    //--------------------------
    // Loop over nuclear interaction data
    N_entries = NT_sec_cluster->GetEntries();
    for(Long64_t i_entry = 0; i_entry < N_entries; i_entry++)
        //for(Long64_t i_entry = 0; i_entry < 50; i_entry++)
    {
        NT_sec_cluster	 ->GetEntry(i_entry);
        //if(ntracks_clus < 7) continue;
        if(dcaTPC_clus > 3.0) continue;
        h2d_cluster_pos_xy->Fill(x_clus,y_clus);
        Float_t clus_temp=TMath::Sqrt(x_clus*x_clus + y_clus*y_clus);
        h1d_cluster_pos_r ->Fill(clus_temp);

        h2D_dEdx_vs_mom ->Fill(mom,dEdx);

        for(Int_t i_clusnbr = 0; i_clusnbr < 5;i_clusnbr++)
        {
            //cout<<"ntracks_clus"<<ntracks_clus<<endl;
            if(ntracks_clus < 4+i_clusnbr) continue;
            h1d_cluster_pos_r_array[i_clusnbr]->Fill(clus_temp);
            //cout<<"fil:"<<clus_temp<<" in"<<i_clusnbr<<endl;
        }
    }
    //--------------------------



    //--------------------------
    h2D_dEdx_vs_mom ->GetYaxis()->SetTitle("TPC dE/dx (KeV/cm)");
    h2D_dEdx_vs_mom ->GetXaxis()->SetTitle("p (GeV/c)");
    h2D_dEdx_vs_mom->GetZaxis()->SetTitle("counts");
    TCanvas* can_dEdx_vs_mom = Draw_2D_histo_and_canvas(h2D_dEdx_vs_mom,"can_dEdx_vs_mom",1010,820,0.0,-1.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_dEdx_vs_mom->cd()->SetRightMargin(0.1);
    can_dEdx_vs_mom->cd()->SetTopMargin(0.05);
    can_dEdx_vs_mom->cd()->SetLogz(1);
    can_dEdx_vs_mom->cd();
    HistName = "p-Pb, #sqrt{s_{NN}}=5.02 TeV";
    plotTopLegend((char*)HistName.Data(),0.26,0.95,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    //--------------------------



    //--------------------------
    gEve->AddElement(TEveP_offset_points);
    TEveP_offset_points  ->SetMarkerSize(1);
    TEveP_offset_points  ->SetMarkerStyle(20);
    TEveP_offset_points  ->SetMarkerColor(kBlue+1);
    gEve->Redraw3D(kTRUE);
    //--------------------------




    //--------------------------
    // Plot y vs x
    h2d_vertex_pos_xy ->GetXaxis()->SetTitle("x (cm)");
    h2d_vertex_pos_xy ->GetYaxis()->SetTitle("y (cm)");
    TCanvas* can_vertex_pos_xy = Draw_2D_histo_and_canvas(h2d_vertex_pos_xy,"can_vertex_pos_xy",1010,820,0.0,-1.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_pos_xy->cd()->SetRightMargin(0.20);
    can_vertex_pos_xy->cd()->SetTopMargin(0.08);
    can_vertex_pos_xy->cd()->SetLogz(0);
    can_vertex_pos_xy->cd();

    for(Int_t i_sector = 0; i_sector < 18; i_sector++)
    {
        vec_PL_TRD_det_2D[i_sector].resize(6); // layers
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            vec_PL_TRD_det_2D[i_sector][i_layer] ->SetLineColor(kBlack);
            vec_PL_TRD_det_2D[i_sector][i_layer] ->SetLineWidth(1);
            vec_PL_TRD_det_2D[i_sector][i_layer] ->SetLineStyle(1);
            vec_PL_TRD_det_2D[i_sector][i_layer] ->Draw("l");
        }
    }
    //--------------------------

    //--------------------------
    // Plot y vs x for cluster
    h2d_cluster_pos_xy ->GetXaxis()->SetTitle("x (cm)");
    h2d_cluster_pos_xy ->GetYaxis()->SetTitle("y (cm)");
    TCanvas* can_cluster_pos_xy = Draw_2D_histo_and_canvas(h2d_cluster_pos_xy,"can_cluster_pos_xy",1010,820,0.0,-1.0,"colz"); // TH2D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_cluster_pos_xy->cd()->SetRightMargin(0.20);
    can_cluster_pos_xy->cd()->SetTopMargin(0.08);
    can_cluster_pos_xy->cd()->SetLogz(0);
    can_cluster_pos_xy->cd();

    for(Int_t i_sector = 0; i_sector < 18; i_sector++)
    {
        vec_PL_TRD_det_2D[i_sector].resize(6); // layers
        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            vec_PL_TRD_det_2D[i_sector][i_layer] ->SetLineColor(kBlack);
            vec_PL_TRD_det_2D[i_sector][i_layer] ->SetLineWidth(1);
            vec_PL_TRD_det_2D[i_sector][i_layer] ->SetLineStyle(1);
            vec_PL_TRD_det_2D[i_sector][i_layer] ->Draw("l");
        }
    }
    //--------------------------
    //--------------------------
    // Plot photons and nuclear interactions together
    TH1D* h1d_vertex_pos_r_norm = (TH1D*)h1d_vertex_pos_r[0] ->Clone("h1d_vertex_pos_r_norm");
    h1d_vertex_pos_r_norm ->Scale(1.0/h1d_vertex_pos_r_norm->GetBinContent(h1d_vertex_pos_r_norm->GetMaximumBin()));
    TH1D* h1d_cluster_pos_r_norm = (TH1D*)h1d_cluster_pos_r ->Clone("h1d_cluster_pos_r_norm");
    h1d_cluster_pos_r_norm ->Scale(1.0/h1d_cluster_pos_r_norm->GetBinContent(h1d_cluster_pos_r_norm->GetMaximumBin()));
    TCanvas* can_vertex_photons_and_nucl = Draw_1D_histo_and_canvas(h1d_vertex_pos_r_norm,"can_vertex_photons_and_nucl",1010,820,0.0,0.0,"h"); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_photons_and_nucl->cd()->SetRightMargin(0.20);
    can_vertex_photons_and_nucl->cd()->SetTopMargin(0.08);
    //can_vertex_photons_and_nucl->cd()->SetLogz(0);
    can_vertex_photons_and_nucl->cd();
    h1d_vertex_pos_r_norm ->SetLineColor(kBlue+2);
    h1d_vertex_pos_r_norm ->SetLineWidth(2);
    h1d_vertex_pos_r_norm ->SetFillColor(kBlue);
    h1d_vertex_pos_r_norm ->SetFillStyle(3022);
    h1d_vertex_pos_r_norm ->DrawCopy("same h");
    h1d_cluster_pos_r_norm ->SetLineColor(kRed+2);
    h1d_cluster_pos_r_norm ->SetLineWidth(2);
    h1d_cluster_pos_r_norm ->SetFillColor(kRed);
    h1d_cluster_pos_r_norm ->SetFillStyle(3022);
    h1d_cluster_pos_r_norm ->DrawCopy("same h");
    for(Int_t i_line = 0;i_line < 6; i_line++)
    {
        //cout<<"im in the loop "<<endl;
        Double_t temp_r=(TRD_layer_radii[i_line][0]+TRD_layer_radii[i_line][1])*0.5;
        if (!(i_line) || (i_line==5))
            PlotArrowHist(h1d_cluster_pos_r_norm,temp_r-1.0,.05,0.02,kBlack,3,1,45.0); //  hist,x_val,arrow_length,arrow_offset,Line_Col,LineWidth,LineStyle,angle
        else
            PlotArrowHist(h1d_vertex_pos_r_norm,temp_r-1.0,.05,0.02,kBlack,3,1,45.0); //  hist,x_val,arrow_length,arrow_offset,Line_Col,LineWidth,LineStyle,angle

    }

    /*
     for (Int_t i_line =0;i_line<6 ;i_line++)
     {
     //Double_t temp_r=(TRD_layer_radii[i_line][0]+TRD_layer_radii[i_line][1])*0.5;
     Double_t temp_r=TRD_layer_radii[i_line][0];
     //PlotLine(temp_r,temp_r,0.0,0.9,kBlue,2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
     PlotLine(temp_r,temp_r,0.0,0.9,color_layer_match[i_line],2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
     temp_r=TRD_layer_radii[i_line][1];
     PlotLine(temp_r,temp_r,0.0,0.9,color_layer_match[i_line],2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

     }*/
    TLegend* legend_p_n = new TLegend(0.22,0.70,0.48,0.8);
    legend_p_n->AddEntry(h1d_vertex_pos_r_norm,"photon vertices","f");
    legend_p_n->AddEntry(h1d_cluster_pos_r_norm,"nuclear vertices","f");
    legend_p_n->SetLineColor(10);
    legend_p_n->Draw();

    can_vertex_photons_and_nucl ->SaveAs("can_vertex_photons_and_nucl.png");
    //----------------------------



    //--------------------------
    // Plot r cluster
    /*  h1d_cluster_pos_r ->GetXaxis()->SetTitle("r (cm)");
     h1d_cluster_pos_r ->GetYaxis()->SetTitle("counts");
     TCanvas* can_cluster_pos_r = Draw_1D_histo_and_canvas(h1d_cluster_pos_r,"can_cluster_pos_r",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
     can_cluster_pos_r->cd()->SetRightMargin(0.20);
     can_cluster_pos_r->cd()->SetTopMargin(0.08);
     can_cluster_pos_r->cd()->SetLogy(1);
     can_cluster_pos_r->cd();
     //----------------------------
     */
    h1d_cluster_pos_r_array[0] ->GetXaxis()->SetTitle("r (cm)");
    h1d_cluster_pos_r_array[0] ->GetYaxis()->SetTitle("counts");
    TCanvas* can_cluster_pos_r = Draw_1D_histo_and_canvas(h1d_cluster_pos_r_array[0],"can_cluster_pos_r",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_cluster_pos_r->cd()->SetRightMargin(0.20);
    can_cluster_pos_r->cd()->SetTopMargin(0.08);
    can_cluster_pos_r->cd()->SetLogy(1);
    can_cluster_pos_r->cd();

    //--------------------------
    // Plot r vertex
    h1d_vertex_pos_r[0] ->GetXaxis()->SetTitle("r (cm)");
    h1d_vertex_pos_r[0] ->GetYaxis()->SetTitle("counts");
    h1d_vertex_pos_r[0]->GetXaxis()->SetRangeUser(260,390);
    TCanvas* can_vertex_pos_r = Draw_1D_histo_and_canvas(h1d_vertex_pos_r[0],"can_vertex_pos_r",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_pos_r->cd()->SetRightMargin(0.20);
    can_vertex_pos_r->cd()->SetTopMargin(0.08);
    can_vertex_pos_r->cd()->SetLogy(0);

    can_vertex_pos_r->cd();
    for(Int_t i_clusnbr = 0; i_clusnbr < 5; i_clusnbr++)
    {
        h1d_vertex_pos_r[i_clusnbr] ->SetLineColor(color_layer_match[i_clusnbr]+2);
        h1d_vertex_pos_r[i_clusnbr] ->SetLineWidth(2);
        h1d_vertex_pos_r[i_clusnbr] ->SetFillColor(color_layer_match[i_clusnbr]);
        h1d_vertex_pos_r[i_clusnbr] ->SetFillStyle(3022);
        h1d_vertex_pos_r[i_clusnbr]->GetXaxis()->SetRangeUser(260,360);

        h1d_vertex_pos_r[i_clusnbr] ->DrawCopy("same h");
    }
    for(Int_t i_line = 0;i_line < 6; i_line++)
    {
        Double_t temp_r=(TRD_layer_radii[i_line][0]+TRD_layer_radii[i_line][1])*0.5;
        PlotArrowHist(h1d_vertex_pos_r[0],temp_r-1.0,200.0,100.0,kBlack,3,1,45); //  hist,x_val,arrow_length,arrow_offset,Line_Col,LineWidth,LineStyle,angle
    }
    TLegend* legend_vert = new TLegend(0.57,0.7,0.8,0.9);
    legend_vert ->SetBorderSize(0);
    legend_vert ->SetFillColor(10);
    legend_vert ->SetTextSize(0.03);
    //legend_vert->SetHeader("The legend_vert Title","C"); // option "C" allows to center the header
    for(Int_t i_clusnbr = 0; i_clusnbr < 5; i_clusnbr++)
    {
        if (!(i_clusnbr))
            HistName="All cuts allowed";
        else
        {
            HistName = "Cut ";
            HistName += i_clusnbr;
        }
        legend_vert->AddEntry(h1d_vertex_pos_r[i_clusnbr],HistName,"f");
    }
    legend_vert->Draw();
    HistName = "p-Pb, #sqrt{s_{NN}}=5.02 TeV, #gamma candidates";
    plotTopLegend((char*)HistName.Data(),0.26,0.95,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    //----------------------------

    //------------------------------------
    //Plot all r for different cuts
    TCanvas* can_cluster_pos_r_comp = Draw_1D_histo_and_canvas(h1d_cluster_pos_r_array[0],"can_vertex_nucl_compare",1010,820,0.0,0.0,"h"); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_cluster_pos_r_comp->cd()->SetRightMargin(0.02);
    can_cluster_pos_r_comp->cd()->SetTopMargin(0.02);
    //can_vertex_photons_and_nucl->cd()->SetLogz(0);
    can_cluster_pos_r_comp->cd();
    for(Int_t i_clusnbr = 0; i_clusnbr < 5; i_clusnbr++)
    {
        h1d_cluster_pos_r_array[i_clusnbr] ->SetLineColor(color_layer_match[i_clusnbr]+2);
        //h1d_cluster_pos_r_array[i_clusnbr] ->SetLineColor(kRed);
        h1d_cluster_pos_r_array[i_clusnbr] ->SetLineWidth(2);
        h1d_cluster_pos_r_array[i_clusnbr] ->SetFillColor(color_layer_match[i_clusnbr]);
        //h1d_cluster_pos_r_array[i_clusnbr] ->SetFillColor(kRed);
        h1d_cluster_pos_r_array[i_clusnbr] ->SetFillStyle(3022);
        h1d_cluster_pos_r_array[i_clusnbr] ->DrawCopy("same h");
    }
    for(Int_t i_line = 0;i_line < 6; i_line++)
    {
        Double_t temp_r=(TRD_layer_radii[i_line][0]+TRD_layer_radii[i_line][1])*0.5;
        //Double_t temp_r=TRD_layer_radii[i_line][0];
        //PlotLine(temp_r,temp_r,0.0,600,kBlue,2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
        //PlotLine(temp_r,temp_r,0.0,0.9,color_layer_match[i_line],2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
        //temp_r=TRD_layer_radii[i_line][1];
        //PlotLine(temp_r,temp_r,0.0,0.9,color_layer_match[i_line],2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)

        PlotArrowHist(h1d_cluster_pos_r_array[0],temp_r-1.0,40.0,10.0,kBlack,3,1,45); //  hist,x_val,arrow_length,arrow_offset,Line_Col,LineWidth,LineStyle,angle

        /*
         TArrow *ar1 = new TArrow(x_val,hist_val+arrow_length+arrow_offset,x_val,hist_val+arrow_offset,0.01,"|>"); // x1,y1,x2,y2
         ar1->SetAngle(30.0);
         ar1->SetLineWidth(2);
         ar1->SetLineColor(vec_circle_color[N_acc_ring]);
         //ar1->SetFillStyle(3008);
         ar1->SetFillColor(vec_circle_color[N_acc_ring]);
         ar1->Draw();
         */
    }
    TLegend* legend = new TLegend(0.75,0.7,0.9,0.92);
    legend ->SetBorderSize(0);
    legend ->SetFillColor(10);
    legend ->SetTextSize(0.03);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    for(Int_t i_clusnbr = 0; i_clusnbr < 5; i_clusnbr++)
    {
        HistName = "Cut n_tracks > ";
        HistName += i_clusnbr +3;
        legend->AddEntry(h1d_cluster_pos_r_array[i_clusnbr],HistName,"f");
    }
    legend->Draw();
    //--------------------------------------




    //-------------------------------
    // Plot pT
    h1d_vertex_mom_pT ->GetXaxis()->SetTitle("pT (GeV/c)");
    h1d_vertex_mom_pT ->GetYaxis()->SetTitle("counts");
    TCanvas* can_vertex_mom_pT = Draw_1D_histo_and_canvas(h1d_vertex_mom_pT,"can_vertex_mom_pT",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_mom_pT->cd()->SetRightMargin(0.20);
    can_vertex_mom_pT->cd()->SetTopMargin(0.08);
    //can_vertex_pos_r->cd()->SetLogz(0);
    can_vertex_mom_pT->cd();
    //-------------------

    //-----------------------------
    // Plot pT_exp
    h1d_vertex_mom_pT_exp ->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1d_vertex_mom_pT_exp ->GetYaxis()->SetTitle("counts");
    TCanvas* can_vertex_mom_pT_exp = Draw_1D_histo_and_canvas(h1d_vertex_mom_pT_exp,"can_vertex_mom_pT_exp",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_mom_pT_exp->cd()->SetRightMargin(0.20);
    can_vertex_mom_pT_exp->cd()->SetTopMargin(0.08);
    can_vertex_mom_pT_exp->cd()->SetLogy();
    can_vertex_mom_pT_exp->cd();
    //-----------------------------------

    //fit

    TF1* func_Exp_fit = new TF1("func_Exp_fit",ExpFitFunc,0.0,1.0,2);

    for(Int_t i = 0; i < 2; i++)
    {
        func_Exp_fit ->ReleaseParameter(i);
        func_Exp_fit ->SetParameter(i,0.0);
        func_Exp_fit ->SetParError(i,0.0);
    }

    func_Exp_fit ->SetParameter(0,1);
    func_Exp_fit ->SetParameter(1,-1);
    func_Exp_fit ->SetRange(0.3,0.8);

    h1d_vertex_mom_pT_exp->Fit("func_Exp_fit","QMN","",0.3,0.8);
    Double_t amplitude = func_Exp_fit ->GetParameter(0);
    Double_t slope     = func_Exp_fit ->GetParameter(1);

    func_Exp_fit ->SetLineStyle(1);
    func_Exp_fit ->SetLineColor(kRed);
    func_Exp_fit ->SetLineWidth(3);
    func_Exp_fit ->DrawCopy("same l");

    cout<<"slope:"<<slope<<endl;
    HistName = "1/slope = ";
    sprintf(NoP,"%4.2f MeV",1000/slope);
    HistName += NoP;
    //can_vertex_mom_pT_exp->cd()->SetLogy(0);

    plotTopLegend((char*)HistName.Data(),0.45,0.83,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


}