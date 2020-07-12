
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


void pT_resolution(Int_t input_file = 0)
{
    TF1* func_Gauss_fit            = new TF1("func_Gauss_fit",GaussFitFunc,0.,10000,3);
    TH2D* TH2D_pT_TPC_vs_Kalman;

    if(input_file == 0) // cov 1.5
    {
        TFile* input_file = TFile::Open("pT_TPC_vs_Kalman.root");
        TCanvas* can_pT_TPC_vs_Kalman = (TCanvas*)input_file->Get("can_pT_TPC_vs_Kalman");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(can_pT_TPC_vs_Kalman->FindObject("TH2D_pT_TPC_vs_Kalman_copy"));
    }
    if(input_file == 1) // cov 1.2
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov1_2.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 2) // cov 2.0
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov2_0.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 3) // cov 4.0
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov4_0.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 4) // cov 10.0   -> best
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov10_0.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 5) // cov 20.0
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov20_0.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 6) // cov 10.0, q/pT = 1.0/1.0
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov10_0_qpT1_0.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 7) // cov 10.0, q/pT = 0.0/1.0, sig2 = 3, sig3 = 9
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov10_0_qpT0_sig2_3_sig3_9.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 8) // cov 10.0, q/pT = 0.0/1.0, sig2 = 9, sig3 = 22
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov10_0_qpT0_sig2_9_sig3_22.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }
    if(input_file == 9) // cov 10.0, q/pT = 0.0/1.0, sig2 = 7, sig3 = 18, 18k events
    {
        TFile* input_file = TFile::Open("TRD_Calib_matched_cov10_0_qpT0_sig2_7_sig3_18_18kevents.root");
        TH2D_pT_TPC_vs_Kalman = (TH2D*)(input_file->Get("TH2D_pT_TPC_vs_Kalman"));
    }



    TCanvas* can_TPC_vs_Kalman = new TCanvas("can_TPC_vs_Kalman","can_TPC_vs_Kalman",10,10,800,800);
    can_TPC_vs_Kalman->cd();
    can_TPC_vs_Kalman->SetLeftMargin(0.22);
    can_TPC_vs_Kalman->SetRightMargin(0.2);
    can_TPC_vs_Kalman->SetBottomMargin(0.15);
    can_TPC_vs_Kalman->SetTopMargin(0.05);
    can_TPC_vs_Kalman->SetLogz(1);

    Double_t Label_size_mean_pT = 0.05;
    TH2D_pT_TPC_vs_Kalman->SetTitle("");
    TH2D_pT_TPC_vs_Kalman->SetStats(0);
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetLabelOffset(0.0);
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetLabelOffset(0.01);
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetLabelSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetLabelSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetTitleSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetTitleSize(Label_size_mean_pT);
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetNdivisions(505,'N');
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetNdivisions(505,'N');
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetTitle("p_{T}^{Kalman} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetTitle("p_{T}^{TPC} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman->GetXaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman->GetYaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman->DrawCopy("colz");

    PlotLine(-4.1,4.1,-4.1,4.1,kRed,2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
   



    //----------------------------------------------------
    // Do projections to get resolution
    TGraphErrors* tge_pT_resolution_RMS = new TGraphErrors();
    Int_t i_point = 0;
    vector<TH1D*> vec_h_proj_x;
    vector<Double_t> vec_pT_TPC;
    for(Int_t biny = 1; biny <= TH2D_pT_TPC_vs_Kalman->GetNbinsY(); biny++)
    {
        Double_t pT_TPC = TH2D_pT_TPC_vs_Kalman->GetYaxis()->GetBinCenter(biny);
        if(fabs(pT_TPC) > 3.0 || fabs(pT_TPC) < 0.35) continue;
        printf("i_point: %d, pT_TPC: %4.3f \n",i_point,pT_TPC);
        TH1D* h_proj_x = TH2D_pT_TPC_vs_Kalman ->ProjectionX("blubb",biny,biny);
        vec_h_proj_x.push_back((TH1D*)h_proj_x->Clone());
        vec_pT_TPC.push_back(pT_TPC);
        Double_t pT_res = h_proj_x ->GetRMS()/fabs(pT_TPC);
        tge_pT_resolution_RMS ->SetPoint(i_point,pT_TPC,pT_res);
        i_point++;
    }
    printf("Number of projections: %d \n",i_point);

    TCanvas* can_pT_res = new TCanvas("can_pT_res","can_pT_res",10,10,850,800);
    can_pT_res->cd();
    can_pT_res->SetLeftMargin(0.22);
    can_pT_res->SetRightMargin(0.05);
    can_pT_res->SetBottomMargin(0.15);
    can_pT_res->SetTopMargin(0.05);
    can_pT_res->cd()->SetTicks(1,1);
    can_pT_res->SetLogz(1);

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
    tge_pT_resolution_RMS->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T} (%)");
    tge_pT_resolution_RMS->SetMarkerStyle(20);
    tge_pT_resolution_RMS->SetMarkerSize(1.1);
    tge_pT_resolution_RMS->SetMarkerColor(kBlack);
    tge_pT_resolution_RMS->GetYaxis()->SetRangeUser(0.0,1.1);
    //tge_pT_resolution_RMS->Draw("AP");


    TGraphErrors* tge_pT_resolution_Gauss = new TGraphErrors();
    for(Int_t i_vec = 0; i_vec < (Int_t)vec_pT_TPC.size(); i_vec++)
    {
        // Do a Gaussian fit
        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        Double_t amplitude = vec_h_proj_x[i_vec]->GetBinContent(vec_h_proj_x[i_vec]->GetMaximumBin());
        Double_t mean      = vec_h_proj_x[i_vec]->GetBinCenter(vec_h_proj_x[i_vec]->GetMaximumBin());
        Double_t sigma     = vec_h_proj_x[i_vec]->GetRMS()*0.3;
        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);


        //Fit
        vec_h_proj_x[i_vec]->Fit("func_Gauss_fit","QMN","",mean-4.0*sigma,mean+4.0*sigma);

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
        vec_h_proj_x[i_vec] ->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

        // Get parameters
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));

        func_Gauss_fit ->SetRange(mean-2.0*sigma,mean+2.0*sigma);

        Double_t pT_TPC = vec_pT_TPC[i_vec];
        Double_t pT_res = sigma/fabs(pT_TPC);

        tge_pT_resolution_Gauss ->SetPoint(i_vec,pT_TPC,pT_res*100.0);
    }

    can_pT_res ->cd();

    tge_pT_resolution_Gauss->SetTitle("");
    tge_pT_resolution_Gauss->GetXaxis()->SetTitleOffset(1.2);
    tge_pT_resolution_Gauss->GetYaxis()->SetTitleOffset(1.5);
    tge_pT_resolution_Gauss->GetXaxis()->SetLabelOffset(0.0);
    tge_pT_resolution_Gauss->GetYaxis()->SetLabelOffset(0.01);
    tge_pT_resolution_Gauss->GetXaxis()->SetLabelSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss->GetYaxis()->SetLabelSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss->GetXaxis()->SetTitleSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss->GetYaxis()->SetTitleSize(Label_size_mean_pT);
    tge_pT_resolution_Gauss->GetXaxis()->SetNdivisions(505,'N');
    tge_pT_resolution_Gauss->GetYaxis()->SetNdivisions(505,'N');
    tge_pT_resolution_Gauss->GetXaxis()->CenterTitle();
    tge_pT_resolution_Gauss->GetYaxis()->CenterTitle();
    tge_pT_resolution_Gauss->GetXaxis()->SetTitle("q*p_{T}^{TPC} (GeV/c)");
    tge_pT_resolution_Gauss->GetYaxis()->SetTitle("TRD-Kalman #sigma(p_{T})/p_{T} (%)");
    tge_pT_resolution_Gauss->SetMarkerStyle(20);
    tge_pT_resolution_Gauss->SetMarkerSize(1.1);
    tge_pT_resolution_Gauss->SetMarkerColor(kBlack);
    tge_pT_resolution_Gauss->GetXaxis()->SetRangeUser(-1.9,1.9);
    tge_pT_resolution_Gauss->GetYaxis()->SetRangeUser(0.0,35.0);
    tge_pT_resolution_Gauss->Draw("AP");

    tge_pT_resolution_Gauss->SetMarkerStyle(20);
    tge_pT_resolution_Gauss->SetMarkerSize(1.1);
    tge_pT_resolution_Gauss->SetMarkerColor(kRed);
    tge_pT_resolution_Gauss->Draw("same P");

    TCanvas* can_proj = new TCanvas("can_proj","can_proj",10,10,1200,1200);
    can_proj->Divide(2,2);
    Int_t arr_proj_sel[4] = {67,74,82,63};
    for(Int_t i_proj_plot = 0; i_proj_plot < 4; i_proj_plot++)
    {
        Int_t i_proj_sel = arr_proj_sel[i_proj_plot];

        can_proj->cd(i_proj_plot+1);
        can_proj->cd(i_proj_plot+1)->SetLeftMargin(0.22);
        can_proj->cd(i_proj_plot+1)->SetRightMargin(0.05);
        can_proj->cd(i_proj_plot+1)->SetBottomMargin(0.15);
        can_proj->cd(i_proj_plot+1)->SetTopMargin(0.05);
        can_proj->cd(i_proj_plot+1)->SetTicks(1,1);
        can_proj->cd(i_proj_plot+1)->SetLogy(0);

        vec_h_proj_x[i_proj_sel]->SetTitle("");
        vec_h_proj_x[i_proj_sel]->SetStats(0);
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetTitleOffset(1.2);
        vec_h_proj_x[i_proj_sel]->GetYaxis()->SetTitleOffset(1.5);
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetLabelOffset(0.0);
        vec_h_proj_x[i_proj_sel]->GetYaxis()->SetLabelOffset(0.01);
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetLabelSize(Label_size_mean_pT);
        vec_h_proj_x[i_proj_sel]->GetYaxis()->SetLabelSize(Label_size_mean_pT);
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetTitleSize(Label_size_mean_pT);
        vec_h_proj_x[i_proj_sel]->GetYaxis()->SetTitleSize(Label_size_mean_pT);
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetNdivisions(505,'N');
        vec_h_proj_x[i_proj_sel]->GetYaxis()->SetNdivisions(505,'N');
        vec_h_proj_x[i_proj_sel]->GetXaxis()->CenterTitle();
        vec_h_proj_x[i_proj_sel]->GetYaxis()->CenterTitle();
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetTitle("p_{T}^{Kalman} (GeV/c)");
        vec_h_proj_x[i_proj_sel]->GetYaxis()->SetTitle("counts");
        vec_h_proj_x[i_proj_sel]->SetLineColor(kBlack);
        //vec_h_proj_x[i_proj_sel]->GetYaxis()->SetRangeUser(0.0,1.1);
        vec_h_proj_x[i_proj_sel]->Draw("h");

        // Do a Gaussian fit
        for(Int_t i = 0; i < 3; i++)
        {
            func_Gauss_fit ->ReleaseParameter(i);
            func_Gauss_fit ->SetParameter(i,0.0);
            func_Gauss_fit ->SetParError(i,0.0);
        }

        Double_t amplitude = vec_h_proj_x[i_proj_sel]->GetBinContent(vec_h_proj_x[i_proj_sel]->GetMaximumBin());
        Double_t mean      = vec_h_proj_x[i_proj_sel]->GetBinCenter(vec_h_proj_x[i_proj_sel]->GetMaximumBin());
        Double_t sigma     = vec_h_proj_x[i_proj_sel]->GetRMS()*0.3;
        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);


        //Fit
        vec_h_proj_x[i_proj_sel]->Fit("func_Gauss_fit","QMN","",mean-4.0*sigma,mean+4.0*sigma);

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
        vec_h_proj_x[i_proj_sel] ->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

        // Get parameters
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));

        func_Gauss_fit ->SetRange(mean-2.0*sigma,mean+2.0*sigma);

        func_Gauss_fit ->SetLineColor(kRed);
        func_Gauss_fit ->SetLineWidth(3);
        func_Gauss_fit ->DrawCopy("same");

        Double_t pT_TPC = vec_pT_TPC[i_proj_sel];
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