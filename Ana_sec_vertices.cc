
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



    //--------------------------
    // Open Ntuple
	TFile* inputfile = TFile::Open("./ST_out/Merge_ST_hoppner_V2.root");
    TNtuple* NT_sec_vertices = (TNtuple*)inputfile->Get("NT_secondary_vertices");
    Float_t x_sec,y_sec,z_sec,ntracks_sec,pT_AB_sec;
    NT_sec_vertices ->SetBranchAddress("x",&x_sec);
    NT_sec_vertices ->SetBranchAddress("y",&y_sec);
    NT_sec_vertices ->SetBranchAddress("z",&z_sec);
    NT_sec_vertices ->SetBranchAddress("ntracks",&ntracks_sec);
	NT_sec_vertices ->SetBranchAddress("pT_AB",&pT_AB_sec);
	
	TNtuple* NT_sec_cluster = (TNtuple*)inputfile->Get("NT_secondary_vertex_cluster");
    Float_t x_clus,y_clus,z_clus,ntracks_clus;
    NT_sec_cluster ->SetBranchAddress("x",&x_clus);
    NT_sec_cluster ->SetBranchAddress("y",&y_clus);
    NT_sec_cluster ->SetBranchAddress("z",&z_clus);
    NT_sec_cluster ->SetBranchAddress("nvertices",&ntracks_clus);
	
    //--------------------------

	Int_t bin_nbr=100;
    TH2D* h2d_vertex_pos_xy = new TH2D("h2d_vertex_pos_xy","h2d_vertex_pos_xy",500,-400,400,500,-400,400);
    TH1D* h1d_vertex_pos_r = new TH1D("h1d_vertex_pos_r","h1d_vertex_pos_r",bin_nbr,200,400);
    TH1D* h1d_vertex_mom_pT = new TH1D("h1d_vertex_mom_pT","h1d_vertex_mom_pT",bin_nbr,0,1);
    TH1D* h1d_vertex_mom_pT_exp = new TH1D("h1d_vertex_mom_pT_exp","h1d_vertex_mom_pT_exp",bin_nbr,0,1);
    TEvePointSet* TEveP_offset_points = new TEvePointSet();

	TH2D* h2d_cluster_pos_xy = new TH2D("h2d_cluster_pos_xy","h2d_cluster_pos_xy",500,-400,400,500,-400,400);
    TH1D* h1d_cluster_pos_r = new TH1D("h1d_cluster_pos_r","h1d_cluster_pos_r",bin_nbr,270,370);
    
    //--------------------------
    // Loop over data
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
		h2d_vertex_pos_xy ->Fill(x_sec,y_sec);
		h1d_vertex_pos_r  ->Fill(TMath::Sqrt(x_sec*x_sec +y_sec*y_sec));
		h1d_vertex_mom_pT ->Fill(pT_AB_sec); 
		
		TEveP_offset_points  ->SetPoint(i_entry,x_sec,y_sec,z_sec);

        //printf("i_entry: %d, pos: {%4.3f, %4.3f, %4.3f} \n",i_entry,x_sec,y_sec,z_sec);
    }
    //--------------------------
	
	N_entries = NT_sec_cluster->GetEntries();
    for(Long64_t i_entry = 0; i_entry < N_entries; i_entry++)
        //for(Long64_t i_entry = 0; i_entry < 50; i_entry++)
    {
    	NT_sec_cluster	 ->GetEntry(i_entry);
        h2d_cluster_pos_xy->Fill(x_clus,y_clus);
		h1d_cluster_pos_r ->Fill(TMath::Sqrt(x_clus*x_clus +y_clus*y_clus));
		
	}
	
	
	for(Int_t i_bin=0; i_bin<bin_nbr; i_bin++)
	{
		Double_t binContent = h1d_vertex_mom_pT->GetBinContent(i_bin);
		Double_t binCenter = h1d_vertex_mom_pT->GetBinCenter(i_bin);
		h1d_vertex_mom_pT_exp->SetBinContent(i_bin,binContent/(binCenter));
		
	}


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
    // Plot r cluster
    h1d_cluster_pos_r ->GetXaxis()->SetTitle("r (cm)");
    h1d_cluster_pos_r ->GetYaxis()->SetTitle("counts");
    TCanvas* can_cluster_pos_r = Draw_1D_histo_and_canvas(h1d_cluster_pos_r,"can_cluster_pos_r",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_cluster_pos_r->cd()->SetRightMargin(0.20);
    can_cluster_pos_r->cd()->SetTopMargin(0.08);
    //can_vertex_pos_r->cd()->SetLogz(0);
    can_cluster_pos_r->cd();
	//----------------------------
	
	//--------------------------
    // Plot r 
    h1d_vertex_pos_r ->GetXaxis()->SetTitle("r (cm)");
    h1d_vertex_pos_r ->GetYaxis()->SetTitle("counts");
    TCanvas* can_vertex_pos_r = Draw_1D_histo_and_canvas(h1d_vertex_pos_r,"can_vertex_pos_r",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_pos_r->cd()->SetRightMargin(0.20);
    can_vertex_pos_r->cd()->SetTopMargin(0.08);
    //can_vertex_pos_r->cd()->SetLogz(0);
    can_vertex_pos_r->cd();
	//----------------------------
	
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
    h1d_vertex_mom_pT_exp ->GetXaxis()->SetTitle("pT (GeV/c)");
    h1d_vertex_mom_pT_exp ->GetYaxis()->SetTitle("counts");
    TCanvas* can_vertex_mom_pT_exp = Draw_1D_histo_and_canvas(h1d_vertex_mom_pT_exp,"can_vertex_mom_pT_exp",1010,820,0.0,0.0,""); // TH1D* hist, TString name, Int_t x_size, Int_t y_size, Double_t min_val, Double_t max_val, TString option
    can_vertex_mom_pT_exp->cd()->SetRightMargin(0.20);
    can_vertex_mom_pT_exp->cd()->SetTopMargin(0.08);
    //can_vertex_pos_r->cd()->SetLogz(0);
    can_vertex_mom_pT_exp->cd();
	//-----------------------------------
	
	//fit
	
	TF1* func_Exp_fit            = new TF1("func_Exp_fit",ExpFitFunc,0.,10000,3);
    func_Exp_fit ->SetParameter(0,1);
    func_Exp_fit ->SetParameter(1,-1);
	
	h1d_vertex_mom_pT_exp->Fit("func_Exp_fit","QMN","",0.3,0.8);
	Double_t amplitude = func_Exp_fit ->GetParameter(0);
    Double_t slope = func_Exp_fit ->GetParameter(1);
        
    cout<<"slope:"<<slope<<endl;    
	 HistName = "#chi_{MC}^{2}/ndf = ";
    sprintf(NoP,"%4.2f",slope);
    HistName += NoP;
     plotTopLegend((char*)HistName.Data(),0.18,0.86,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
	

}