
static TString HistName;
static char NoP[50];

//----------------------------------------------------------------------------------------
TLine* PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
    TLine* Zero_line = new TLine();
    Zero_line -> SetX1(x1_val);
    Zero_line -> SetX2(x2_val);
    Zero_line -> SetY1(y1_val);
    Zero_line -> SetY2(y2_val);
    Zero_line -> SetLineWidth(LineWidth);
    Zero_line -> SetLineStyle(LineStyle);
    Zero_line -> SetLineColor(Line_Col);
    Zero_line -> Draw();
    //delete Zero_line;

    return Zero_line;
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    // Text alignment: https://root.cern.ch/doc/master/classTAttText.html#T1
    // align = 10*HorizontalAlign + VerticalAlign
    // horizontal: 1=left adjusted, 2=centered, 3=right adjusted
    // vertical: 1=bottom adjusted, 2=centered, 3=top adjusted


    if((x<0||y<0) && NDC == 1)
    {   // defaults
      x=gPad->GetLeftMargin()*1.15;
      y=(1-gPad->GetTopMargin())*1.04;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->SetTextAlign(align);
    text->Draw();
    return text;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Double_t GaussFitFunc(Double_t* x_val, Double_t* par)
{
    Double_t x, y, par0, par1, par2;
    par0  = fabs(par[0]);
    par1  = par[1];
    par2  = fabs(par[2]);
    x = x_val[0];
    y = par0*TMath::Gaus(x,par1,par2,0);
    return y;
}
//----------------------------------------------------------------------------------------


void pT_resolution()
{
    const Int_t N_pT_resolution = 7;

    TF1* func_Gauss_fit            = new TF1("func_Gauss_fit",GaussFitFunc,0.,10000,3);
    vector<TH2D*> TH2D_pT_TPC_vs_Kalman;
    vector<vector<TH2D*>> vec_TH2D_pT_TPC_vs_Kalman;
    vec_TH2D_pT_TPC_vs_Kalman.resize(2);
    vec_TH2D_pT_TPC_vs_Kalman[0].resize(N_pT_resolution);
	vec_TH2D_pT_TPC_vs_Kalman[1].resize(N_pT_resolution);

    vector<TFile*> input_file;
	input_file.push_back(TFile::Open("./ST_out/Merge_ST_ABCDE_V9.root"));
	input_file.push_back(TFile::Open("./ST_out/Merge_ST_ABCDE_V9_prim_vertex.root"));
	//TFile* input_file = TFile::Open("./ST_out/Merge_ST_ABCDE_V9.root");
	//TFile* input_file_b = TFile::Open("./ST_out/Merge_ST_ABCDE_V9_prim_vertex.root");
	TH2D_pT_TPC_vs_Kalman.push_back( (TH2D*)(input_file[0]->Get("TH2D_pT_TPC_vs_Kalman")));
	TH2D_pT_TPC_vs_Kalman.push_back( (TH2D*)(input_file[1]->Get("TH2D_pT_TPC_vs_Kalman")));

	for(Int_t i_pt_res = 0; i_pt_res < N_pT_resolution; i_pt_res++)
	{
		HistName = "vec_TH2D_pT_TPC_vs_Kalman_";
		HistName += i_pt_res;
		vec_TH2D_pT_TPC_vs_Kalman[0][i_pt_res]=(TH2D*)(input_file[0]->Get(HistName.Data()));
		vec_TH2D_pT_TPC_vs_Kalman[1][i_pt_res]=(TH2D*)(input_file[1]->Get(HistName.Data()));
	}

	//TH2D_pT_TPC_vs_Kalman = vec_TH2D_pT_TPC_vs_Kalman[6]; // set number of tracklets for Kalman track, 4..6

	TH2D_pT_TPC_vs_Kalman[0] = vec_TH2D_pT_TPC_vs_Kalman[0][6];
    TH2D_pT_TPC_vs_Kalman[1] = vec_TH2D_pT_TPC_vs_Kalman[1][6];
    TCanvas* can_TPC_vs_Kalman = new TCanvas("can_TPC_vs_Kalman","can_TPC_vs_Kalman",10,10,800,800);
    can_TPC_vs_Kalman->cd();
    can_TPC_vs_Kalman->SetLeftMargin(0.22);
    can_TPC_vs_Kalman->SetRightMargin(0.2);
    can_TPC_vs_Kalman->SetBottomMargin(0.15);
    can_TPC_vs_Kalman->SetTopMargin(0.05);
    can_TPC_vs_Kalman->SetLogz(1);
	
    Double_t Label_size_mean_pT = 0.05;
    TH2D_pT_TPC_vs_Kalman[0]->SetTitle("");
    TH2D_pT_TPC_vs_Kalman[0]->SetStats(0);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetLabelOffset(0.0);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetLabelOffset(0.01);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetLabelSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetLabelSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetTitleSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetTitleSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetNdivisions(505,'N');
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetNdivisions(505,'N');
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetTitle("p_{T}^{Kalman} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetTitle("p_{T}^{TPC} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman[0]->DrawCopy("colz");
    PlotLine(-4.1,4.1,-4.1,4.1,kRed,2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
 
	TCanvas* can_TPC_vs_Kalman_prim_vert = new TCanvas("can_TPC_vs_Kalman_prim_vert","can_TPC_vs_Kalman_prim_vert",10,10,800,800);
    can_TPC_vs_Kalman_prim_vert->cd();
    can_TPC_vs_Kalman_prim_vert->SetLeftMargin(0.22);
    can_TPC_vs_Kalman_prim_vert->SetRightMargin(0.2);
    can_TPC_vs_Kalman_prim_vert->SetBottomMargin(0.15);
    can_TPC_vs_Kalman_prim_vert->SetTopMargin(0.05);
    can_TPC_vs_Kalman_prim_vert->SetLogz(1);

    //Double_t Label_size_mean_pT = 0.05;
    TH2D_pT_TPC_vs_Kalman[1]->SetTitle("");
    TH2D_pT_TPC_vs_Kalman[1]->SetStats(0);
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetLabelOffset(0.0);
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetLabelOffset(0.01);
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetLabelSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetLabelSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetTitleSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetTitleSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetNdivisions(505,'N');
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetNdivisions(505,'N');
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetTitle("p_{T}^{Kalman} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetTitle("p_{T}^{TPC} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman[1]->GetXaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman[1]->DrawCopy("colz");

    PlotLine(-4.1,4.1,-4.1,4.1,kRed,2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
   
    //----------------------------------------------------
    // Do projections to get resolution
    //TGraphErrors* tge_pT_resolution_RMS = new TGraphErrors();
    Int_t i_point = 0;
    vector<vector<TH1D*>> vec_h_proj_x;
	vec_h_proj_x.resize(6);
    vector<vector<Double_t>> vec_pT_TPC;
	vec_pT_TPC.resize(6);
	//TH2D_pT_TPC_vs_Kalman[0] = vec_TH2D_pT_TPC_vs_Kalman[0][6];
    //TH2D_pT_TPC_vs_Kalman[1] = vec_TH2D_pT_TPC_vs_Kalman[1][6];
    
	for(Int_t i_all = 5; i_all >=0; i_all--)
	{
    	for(Int_t biny = 1; biny <= vec_TH2D_pT_TPC_vs_Kalman[1-(Int_t)i_all/3][(i_all+1) +3*(1-(Int_t)i_all/3)]->GetNbinsY(); biny++)
    	{
			Double_t pT_TPC = vec_TH2D_pT_TPC_vs_Kalman[1-(Int_t)i_all/3][i_all+1+3*(1-(Int_t)i_all/3)]->GetYaxis()->GetBinCenter(biny);
			if(fabs(pT_TPC) > 4.0 || fabs(pT_TPC) < 0.35) continue;
			printf("i_point: %d, pT_TPC: %4.3f \n",i_point,pT_TPC);
			TH1D* h_proj_x = vec_TH2D_pT_TPC_vs_Kalman[1-(Int_t)i_all/3][i_all+1+3*(1-(Int_t)i_all/3)] ->ProjectionX("blubb",biny,biny);
			vec_h_proj_x[5-i_all].push_back((TH1D*)h_proj_x->Clone());
			vec_pT_TPC[5-i_all].push_back(pT_TPC);
			//Double_t pT_res = h_proj_x ->GetRMS()/fabs(pT_TPC);
			//tge_pT_resolution_RMS ->SetPoint(i_point,pT_TPC,pT_res);
			i_point++;
		}
	}
	/*for(Int_t biny = 1; biny <= TH2D_pT_TPC_vs_Kalman[1]->GetNbinsY(); biny++)
    {
        Double_t pT_TPC = TH2D_pT_TPC_vs_Kalman[1]->GetYaxis()->GetBinCenter(biny);
        if(fabs(pT_TPC) > 4.0 || fabs(pT_TPC) < 0.35) continue;
        printf("i_point: %d, pT_TPC: %4.3f \n",i_point,pT_TPC);
        TH1D* h_proj_x = TH2D_pT_TPC_vs_Kalman[1] ->ProjectionX("blubb",biny,biny);
        vec_h_proj_x[3*0 +0].push_back((TH1D*)h_proj_x->Clone());
        vec_pT_TPC[3*1 +0].push_back(pT_TPC);
        //Double_t pT_res = h_proj_x ->GetRMS()/fabs(pT_TPC);
        //tge_pT_resolution_RMS ->SetPoint(i_point,pT_TPC,pT_res);
        i_point++;
    }*/
    printf("Number of projections: %d \n",i_point);

    TCanvas* can_pT_res = new TCanvas("can_pT_res","can_pT_res",10,10,850,800);
    can_pT_res->cd();
    can_pT_res->SetLeftMargin(0.22);
    can_pT_res->SetRightMargin(0.05);
    can_pT_res->SetBottomMargin(0.15);
    can_pT_res->SetTopMargin(0.05);
    can_pT_res->cd()->SetTicks(1,1);
    can_pT_res->SetLogz(1);
/*
    tge_pT_resolution_RMS->SetTitle("");
    tge_pT_resolution_RMS->GetXaxis()->SetTitleOffset(1.2);
    tge_pT_resolution_RMS->GetYaxis()->SetTitleOffset(1.5);
    tge_pT_resolution_RMS->GetXaxis()->SetLabelOffset(0.0);
    tge_pT_resolution_RMS->GetYaxis()->SetLabelOffset(0.01);
    tge_pT_resolution_RMS->GetXaxis()->SetLabelSize(Label_size_mean_pT);
    tge_pT_resolution_RMS->GetYaxis()->SetLabelSize(Label_size_mean_pT);
    tge_pT_resolution_RMS->GetXaxis()->SetTitleSize(Label_size_mean_pT);
    tge_pT_resolution_RMS->GetYaxis()->SetTitleSize(Label_size_mean_pT);
    tge_pT_resolution_RMS->GetXaxis()->SetNdivisions(505,'N');
    tge_pT_resolution_RMS->GetYaxis()->SetNdivisions(505,'N');
    tge_pT_resolution_RMS->GetXaxis()->CenterTitle();
    tge_pT_resolution_RMS->GetYaxis()->CenterTitle();
    tge_pT_resolution_RMS->GetXaxis()->SetTitle("q*p_{T}^{TPC} (GeV/c)");
    tge_pT_resolution_RMS->GetYaxis()->SetTitle("#sigma(p_{T}^{Kalm})/p_{T}^{TPC} (%)");
    tge_pT_resolution_RMS->SetMarkerStyle(20);
    tge_pT_resolution_RMS->SetMarkerSize(1.1);
    tge_pT_resolution_RMS->SetMarkerColor(kBlack);
    tge_pT_resolution_RMS->GetYaxis()->SetRangeUser(0.0,1.1);
    //tge_pT_resolution_RMS->Draw("AP");
*/

    vector<TGraphErrors*> tge_pT_resolution_Gauss; 
	tge_pT_resolution_Gauss.resize(6);
	for(Int_t i_all = 0; i_all <=5; i_all++)
	{
    
	tge_pT_resolution_Gauss[i_all]= new TGraphErrors();
    for(Int_t i_vec = 0; i_vec < (Int_t)vec_pT_TPC[i_all].size(); i_vec++)
    {
        // Do a Gaussian fit
        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        Double_t amplitude = vec_h_proj_x[i_all][i_vec]->GetBinContent(vec_h_proj_x[i_all][i_vec]->GetMaximumBin());
        Double_t mean      = vec_h_proj_x[i_all][i_vec]->GetBinCenter(vec_h_proj_x[i_all][i_vec]->GetMaximumBin());
        Double_t sigma     = vec_h_proj_x[i_all][i_vec]->GetRMS()*1.0;
        //Double_t sigma     = 0.02;
        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->FixParameter(1,mean);
        func_Gauss_fit ->FixParameter(2,sigma);


        //Fit
        vec_h_proj_x[i_all][i_vec]->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));


        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);

        //Fit
        vec_h_proj_x[i_all][i_vec] ->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

        // Get parameters
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));

        func_Gauss_fit ->SetRange(mean-2.0*sigma,mean+2.0*sigma);

        Double_t pT_TPC = vec_pT_TPC[0][i_vec];
        Double_t pT_res = sigma/fabs(pT_TPC);

        tge_pT_resolution_Gauss[i_all] ->SetPoint(i_vec,pT_TPC,pT_res*100.0);
    }
	}
    can_pT_res ->cd();

	Int_t color_layer_match[6] = {kGray+1,kAzure-2,kGreen+2,kCyan+2,kOrange+2,kRed};
	
    tge_pT_resolution_Gauss[0]->SetTitle("");
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetTitleOffset(1.2);
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetTitleOffset(1.5);
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetLabelOffset(0.0);
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetLabelOffset(0.01);
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetLabelSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetLabelSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetTitleSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetTitleSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetNdivisions(505,'N');
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetNdivisions(505,'N');
    tge_pT_resolution_Gauss[0]->GetXaxis()->CenterTitle();
    tge_pT_resolution_Gauss[0]->GetYaxis()->CenterTitle();
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetTitle("q*p_{T}^{TPC} (GeV/c)");
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetTitle("#sigma^{Kalm-TRD}(p_{T})/p_{T}^{TPC} (%)");
    tge_pT_resolution_Gauss[0]->SetMarkerStyle(20);
    tge_pT_resolution_Gauss[0]->SetMarkerSize(1.1);
    tge_pT_resolution_Gauss[0]->SetMarkerColor(kBlack);
    tge_pT_resolution_Gauss[0]->GetXaxis()->SetRangeUser(-3.52,3.52);
    tge_pT_resolution_Gauss[0]->GetYaxis()->SetRangeUser(0.0,35.0);
    tge_pT_resolution_Gauss[0]->Draw("AP");
	tge_pT_resolution_Gauss[0]->SetMarkerColor(color_layer_match[0]);
    

/*    tge_pT_resolution_Gauss->SetMarkerStyle(20);
    tge_pT_resolution_Gauss->SetMarkerSize(1.1);
    tge_pT_resolution_Gauss->SetMarkerColor(kRed);
    tge_pT_resolution_Gauss->Draw("same P");
*/	
	
	TLegend* legend = new TLegend(0.23,0.7,0.4,0.92);
    legend ->SetBorderSize(0);
    legend ->SetFillColor(10);
    legend ->SetTextSize(0.03);
    HistName = "no prim vertex w. ";
	HistName += 6;
	HistName += "tracklets";
	legend->AddEntry(tge_pT_resolution_Gauss[0],HistName,"f");	
	
	for(Int_t i_all = 1; i_all <6; i_all++)
	{
	if (i_all<3)	
	HistName = "no prim vertex w. ";
	else	
   	HistName = "prim vertex w. ";
	HistName += 6 -i_all +3*((Int_t)i_all/3);
	HistName += "tracklets";
		
    legend->AddEntry(tge_pT_resolution_Gauss[i_all],HistName,"f");	
		
    tge_pT_resolution_Gauss[i_all]->SetMarkerStyle(20);
    tge_pT_resolution_Gauss[i_all]->SetMarkerSize(1.1);
    tge_pT_resolution_Gauss[i_all]->SetMarkerColor(color_layer_match[i_all%3]);
    tge_pT_resolution_Gauss[i_all]->Draw("same P");
	}
	legend->Draw();
    
	
	
	
	
	
	

    TCanvas* can_proj = new TCanvas("can_proj","can_proj",10,10,1200,1000);
    can_proj->Divide(3,2);
    //Int_t arr_proj_sel[6] = {67,74,82,63,109,119};
    Int_t arr_proj_sel[6] = {370,74,82,63,109,119};
    for(Int_t i_proj_plot = 0; i_proj_plot < 6; i_proj_plot++)
    {
        Int_t i_proj_sel = arr_proj_sel[i_proj_plot];

        can_proj->cd(i_proj_plot+1);
        can_proj->cd(i_proj_plot+1)->SetLeftMargin(0.22);
        can_proj->cd(i_proj_plot+1)->SetRightMargin(0.05);
        can_proj->cd(i_proj_plot+1)->SetBottomMargin(0.15);
        can_proj->cd(i_proj_plot+1)->SetTopMargin(0.05);
        can_proj->cd(i_proj_plot+1)->SetTicks(1,1);
        can_proj->cd(i_proj_plot+1)->SetLogy(0);

        vec_h_proj_x[3*0 +0][i_proj_sel]->SetTitle("");
        vec_h_proj_x[3*0 +0][i_proj_sel]->SetStats(0);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->SetTitleOffset(1.2);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetTitleOffset(1.5);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->SetLabelOffset(0.0);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetLabelOffset(0.01);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->SetLabelSize(Label_size_mean_pT);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetLabelSize(Label_size_mean_pT);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->SetTitleSize(Label_size_mean_pT);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetTitleSize(Label_size_mean_pT);
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->SetNdivisions(505,'N');
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetNdivisions(505,'N');
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->CenterTitle();
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->CenterTitle();
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetXaxis()->SetTitle("p_{T}^{Kalman} (GeV/c)");
        vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetTitle("counts");
        vec_h_proj_x[3*0 +0][i_proj_sel]->SetLineColor(kBlack);
        //vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetRangeUser(0.0,1.1);
        vec_h_proj_x[3*0 +0][i_proj_sel]->Draw("h");

        // Do a Gaussian fit
        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        Double_t amplitude = vec_h_proj_x[3*0 +0][i_proj_sel]->GetBinContent(vec_h_proj_x[3*0 +0][i_proj_sel]->GetMaximumBin());
        Double_t mean      = vec_h_proj_x[3*0 +0][i_proj_sel]->GetBinCenter(vec_h_proj_x[3*0 +0][i_proj_sel]->GetMaximumBin());
        Double_t sigma     = vec_h_proj_x[3*0 +0][i_proj_sel]->GetRMS()*1.0;
        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);

        printf("mean: %4.3f, sigma: %4.3f, amplitude: %4.3f \n",mean,sigma,amplitude);

        //Fit
        vec_h_proj_x[3*0 +0][i_proj_sel]->Fit("func_Gauss_fit","QMN","",mean-4.0*sigma,mean+4.0*sigma);

        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));


        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);

        //Fit
        vec_h_proj_x[3*0 +0][i_proj_sel] ->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

        // Get parameters
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));

        func_Gauss_fit ->SetRange(mean-2.0*sigma,mean+2.0*sigma);

        func_Gauss_fit ->SetLineColor(kRed);
        func_Gauss_fit ->SetLineWidth(3);
        func_Gauss_fit ->DrawCopy("same");

        Double_t pT_TPC = vec_pT_TPC[3*0 +0][i_proj_sel];
        Double_t pT_res = sigma/fabs(pT_TPC);

        HistName = "p_{T}^{TPC} = ";
        sprintf(NoP,"%4.1f",pT_TPC);
        HistName += NoP;
        HistName += " GeV/c";
        plotTopLegend((char*)HistName.Data(),0.3,0.86,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        HistName = "#sigma = ";
        sprintf(NoP,"%4.3f",sigma);
        HistName += NoP;
        plotTopLegend((char*)HistName.Data(),0.3,0.8,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        HistName = "R = ";
        sprintf(NoP,"%4.1f",pT_res*100.0);
        HistName += NoP;
        HistName += "%";
        plotTopLegend((char*)HistName.Data(),0.3,0.74,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
 
    }
    //----------------------------------------------------



}