
static TString HistName;
static char NoP[50];

//for generator
#include "TStyle.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TFile.h"

//----------------------------------------------------------------------------------------
void SetRootGraphicStyle()
{
    cout << "Set basic ROOT graphics style" << endl;
    //gStyle->Reset();
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetFrameFillColor(10);
    //gStyle->SetFillColor(4);
    TGaxis::SetMaxDigits(4);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadRightMargin(0.14);
    gStyle->SetPadLeftMargin(0.18);
    gStyle->SetLabelSize(0.07,"X");
    gStyle->SetLabelSize(0.07,"Y");
    gStyle->SetTitleSize(0.07,"X");
    gStyle->SetTitleSize(0.07,"Y");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,greens, blues, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetEndErrorSize(3);
    TRandom3 r3b;
    r3b.SetSeed(0); // seed for random number generator changes every second
    gRandom->SetSeed(0);
}
//----------------------------------------------------------------------------------------



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


void pT_resolution_physics()
{

    SetRootGraphicStyle();
    const Int_t N_pT_resolution = 7;

    TF1* func_Gauss_fit            = new TF1("func_Gauss_fit",GaussFitFunc,0.,10000,3);
    vector<TH2D*> TH2D_pT_TPC_vs_Kalman;
    vector<vector<TH2D*>> vec_TH2D_pT_TPC_vs_Kalman;
    vec_TH2D_pT_TPC_vs_Kalman.resize(2);
    vec_TH2D_pT_TPC_vs_Kalman[0].resize(N_pT_resolution);
    vec_TH2D_pT_TPC_vs_Kalman[1].resize(N_pT_resolution);

    vector<vector<TH2D*>> vec_TH2D_one_over_pT_TPC_vs_Kalman;
    vec_TH2D_one_over_pT_TPC_vs_Kalman.resize(2);
    vec_TH2D_one_over_pT_TPC_vs_Kalman[0].resize(N_pT_resolution);
    vec_TH2D_one_over_pT_TPC_vs_Kalman[1].resize(N_pT_resolution);

    vector<TFile*> input_file;
    //input_file.push_back(TFile::Open("./ST_out/Merge_ST_ABCDE_V9.root"));
    //input_file.push_back(TFile::Open("./ST_out/Merge_ST_ABCDE_V9_prim_vertex.root"));
    input_file.push_back(TFile::Open("../ST_out/Merge_photons_nucl_A_V0_bug_fixed_from_NASTIA_and_TLV.root"));
    //input_file.push_back(TFile::Open("./ST_out/Merge_ST_ABCDE_V13_PV.root"));
    //TFile* input_file = TFile::Open("./ST_out/Merge_ST_ABCDE_V9.root");
    //TFile* input_file_b = TFile::Open("./ST_out/Merge_ST_ABCDE_V9_prim_vertex.root");

    //TH2D_pT_TPC_vs_Kalman.push_back( (TH2D*)(input_file[0]->Get("TH2D_pT_TPC_vs_Kalman"))); /normally
    TH2D_pT_TPC_vs_Kalman.push_back( (TH2D*)(input_file[0]->Get("vec_TH2D_pT_TPC_vs_Kalman_6")));  //test

    //TH2D_pT_TPC_vs_Kalman.push_back( (TH2D*)(input_file[1]->Get("TH2D_pT_TPC_vs_Kalman")));

    for(Int_t i_pt_res = 0; i_pt_res < N_pT_resolution; i_pt_res++)
    {
        HistName = "vec_TH2D_pT_TPC_vs_Kalman_";
        HistName += i_pt_res;
        vec_TH2D_pT_TPC_vs_Kalman[0][i_pt_res]=(TH2D*)(input_file[0]->Get(HistName.Data()));
        //vec_TH2D_pT_TPC_vs_Kalman[1][i_pt_res]=(TH2D*)(input_file[1]->Get(HistName.Data()));

        HistName = "vec_TH2D_one_over_pT_TPC_vs_Kalman_";
        HistName += i_pt_res;
        vec_TH2D_one_over_pT_TPC_vs_Kalman[0][i_pt_res]=(TH2D*)(input_file[0]->Get(HistName.Data()));
        //vec_TH2D_one_over_pT_TPC_vs_Kalman[1][i_pt_res]=(TH2D*)(input_file[1]->Get(HistName.Data()));
    }

    //TH2D_pT_TPC_vs_Kalman = vec_TH2D_pT_TPC_vs_Kalman[6]; // set number of tracklets for Kalman track, 4..6




    //----------------------------------------------------
    // pT correlation


    //TH2D_pT_TPC_vs_Kalman[0] = vec_TH2D_pT_TPC_vs_Kalman[0][6]; //removed for test 


    //TH2D_pT_TPC_vs_Kalman[1] = vec_TH2D_pT_TPC_vs_Kalman[1][6];
    TCanvas* can_TPC_vs_Kalman = new TCanvas("can_TPC_vs_Kalman","can_TPC_vs_Kalman",10,10,800,800);
    can_TPC_vs_Kalman->cd();
    can_TPC_vs_Kalman->SetLeftMargin(0.22);
    can_TPC_vs_Kalman->SetRightMargin(0.2);
    can_TPC_vs_Kalman->SetBottomMargin(0.15);
    can_TPC_vs_Kalman->SetTopMargin(0.07);
    can_TPC_vs_Kalman->SetLogz(1);

    Double_t Label_size_mean_pT = 0.05;
    TH2D_pT_TPC_vs_Kalman[0]->SetTitle("");
    TH2D_pT_TPC_vs_Kalman[0]->SetStats(0);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetTitleOffset(1.2);
    TH2D_pT_TPC_vs_Kalman[0]->GetZaxis()->SetTitleOffset(1.2);
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
    TH2D_pT_TPC_vs_Kalman[0]->GetZaxis()->CenterTitle();
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetTitleFont(42);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetTitleFont(42);
    TH2D_pT_TPC_vs_Kalman[0]->GetZaxis()->SetTitleFont(42);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetLabelFont(42);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetLabelFont(42);
    TH2D_pT_TPC_vs_Kalman[0]->GetZaxis()->SetLabelFont(42);
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetTitle("q*p_{T}^{Kalman} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetTitle("q*p_{T}^{TPC} (GeV/c)");
    TH2D_pT_TPC_vs_Kalman[0]->GetZaxis()->SetTitle("counts");
    TH2D_pT_TPC_vs_Kalman[0]->GetXaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->SetRangeUser(-4.2,4.2);
    TH2D_pT_TPC_vs_Kalman[0]->DrawCopy("colz");
    PlotLine(-4.1,4.1,-4.1,4.1,kGray,2,9); // (Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
    HistName = "p-Pb, #sqrt{s_{NN}}=5.02 TeV, N_{trk}= 6";
    plotTopLegend((char*)HistName.Data(),0.26,0.95,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    HistName = "TRD only";
    plotTopLegend((char*)HistName.Data(),0.26,0.88,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1


    //----------------------------------------------------
    // Do projections to get resolution
    //TGraphErrors* tge_pT_resolution_RMS = new TGraphErrors();
    Int_t i_point = 0;
    vector<TH1D*> vec_h_proj_x;
    //vec_h_proj_x.resize(6);
    vector<Double_t> vec_pT_TPC;
    //vec_pT_TPC.resize(6);
    //TH2D_pT_TPC_vs_Kalman[0] = vec_TH2D_pT_TPC_vs_Kalman[0][6];
    //TH2D_pT_TPC_vs_Kalman[1] = vec_TH2D_pT_TPC_vs_Kalman[1][6];

    for(Int_t biny = 1; biny <= TH2D_pT_TPC_vs_Kalman[0]->GetNbinsY(); biny++)
    {
        Double_t pT_TPC = TH2D_pT_TPC_vs_Kalman[0]->GetYaxis()->GetBinCenter(biny);
        if(fabs(pT_TPC) > 4.0 || fabs(pT_TPC) < 0.35) continue;
        //printf("i_point: %d, pT_TPC: %4.3f \n",i_point,pT_TPC);
        TH1D* h_proj_x;
        if (fabs(pT_TPC) < 1.2 )
            h_proj_x = TH2D_pT_TPC_vs_Kalman[0] ->ProjectionX("blubb",biny,biny);
        else
        {
            h_proj_x = TH2D_pT_TPC_vs_Kalman[0] ->ProjectionX("blubb",biny,biny+2);
            biny+=2;
        }
        vec_h_proj_x.push_back((TH1D*)h_proj_x->Clone());
        vec_pT_TPC.push_back(pT_TPC);
            //Double_t pT_res = h_proj_x ->GetRMS()/fabs(pT_TPC);
            //tge_pT_resolution_RMS ->SetPoint(i_point,pT_TPC,pT_res);
        i_point++;
    }
    
    printf("Number of projections: %d \n",i_point);



    //------------------------------------------------------
    TCanvas* can_pT_res = new TCanvas("can_pT_res","can_pT_res",10,10,850,800);
    can_pT_res->cd();
    can_pT_res->SetLeftMargin(0.22);
    can_pT_res->SetRightMargin(0.05);
    can_pT_res->SetBottomMargin(0.15);
    can_pT_res->SetTopMargin(0.07);
    can_pT_res->cd()->SetTicks(1,1);
    can_pT_res->SetLogz(1);

    Label_size_mean_pT = 0.05;

    TGraphErrors* tge_pT_resolution_Gauss;
    TGraphErrors* tge_pT_mean_Gauss;

    //tge_pT_resolution_Gauss.resize(6);
    //for(Int_t i_all = 0; i_all <= 5; i_all++)
    //{

    tge_pT_resolution_Gauss = new TGraphErrors();
    tge_pT_mean_Gauss       = new TGraphErrors();

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
        Double_t sigma     = vec_h_proj_x[i_vec]->GetRMS()*1.0;
            //Double_t sigma     = 0.02;
        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->FixParameter(1,mean);
        func_Gauss_fit ->FixParameter(2,sigma);


            //Fit
        vec_h_proj_x[i_vec]->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

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
        vec_h_proj_x[i_vec] ->Fit("func_Gauss_fit","QMN","",mean-1.2*sigma,mean+1.2*sigma);

            // Get parameters
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));

        func_Gauss_fit ->SetRange(mean-1.2*sigma,mean+1.2*sigma);

        Double_t pT_TPC  = vec_pT_TPC[i_vec];
        Double_t pT_res  = sigma/fabs(pT_TPC);
        Double_t pT_mean = mean;

        tge_pT_resolution_Gauss ->SetPoint(i_vec,pT_TPC,pT_res*100.0);
        tge_pT_mean_Gauss       ->SetPoint(i_vec,mean,pT_TPC);
    }

    can_pT_res ->cd();

    Int_t color_layer_match[6] = {kRed+1,kGray+1,kAzure-2,kCyan+2,kOrange+2,kRed};
    Int_t marker_layer_match[6] = {24,26,30,20,22,29};
    Double_t marker_size_match[6] = {0.8,0.9,1.0,0.8,0.9,1.0};

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
    tge_pT_resolution_Gauss->GetYaxis()->SetTitle("#sigma^{Kalm-TRD}(p_{T})/p_{T}^{TPC} (%)");
    tge_pT_resolution_Gauss->GetXaxis()->SetRangeUser(-3.52,3.52);
    tge_pT_resolution_Gauss->GetYaxis()->SetRangeUser(0.0,35.0);


    TLegend* legend_TRD_only = new TLegend(0.72,0.56,0.92,0.68);
    legend_TRD_only ->SetBorderSize(0);
    legend_TRD_only ->SetFillColor(10);
    legend_TRD_only ->SetTextSize(0.04);

    HistName = "N_{trk}=6";
        

    tge_pT_resolution_Gauss->SetMarkerStyle(marker_layer_match[0]);
    tge_pT_resolution_Gauss->SetMarkerSize(marker_size_match[0]);
    tge_pT_resolution_Gauss->SetMarkerColor(color_layer_match[0]);
    tge_pT_resolution_Gauss->Draw("AP");
        //else tge_pT_resolution_Gauss[i_all]->Draw("same P");

    legend_TRD_only->AddEntry(tge_pT_resolution_Gauss,HistName,"p");


    legend_TRD_only->Draw();
    plotTopLegend((char*)"TRD only",0.74,0.685,0.04,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    //plotTopLegend((char*)"TRD+PV",0.32,0.485,0.04,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    plotTopLegend((char*)"p-Pb, #sqrt{s_{NN}}=5.02 TeV",0.28,0.95,0.04,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    
    //------------------------------------------------------
    //-- Add TRD pT on top of 2D hist ----------------------
    can_TPC_vs_Kalman   ->cd();
    tge_pT_mean_Gauss   ->SetMarkerStyle(20);
    tge_pT_mean_Gauss   ->SetMarkerSize(0.5);
    tge_pT_mean_Gauss   ->SetMarkerColor(kRed);
    //tge_pT_mean_Gauss   ->Draw("same P");

    Double_t x_pos_neg;
    Double_t y_pos_neg;

    TGraphErrors* tge_pT_mean_Gauss_pos = new TGraphErrors();
    TGraphErrors* tge_pT_mean_Gauss_neg = new TGraphErrors();
    Int_t i_points_neg = 0;
    Int_t i_points_pos = 0;

    for (Int_t i_point = 0; i_point < tge_pT_mean_Gauss->GetN(); i_point++)
    {
        tge_pT_mean_Gauss ->GetPoint(i_point,x_pos_neg,y_pos_neg);

        if (x_pos_neg > -1.2 && x_pos_neg < -0.3 && y_pos_neg > -2.0 && y_pos_neg < 0.0) 
            {tge_pT_mean_Gauss_neg ->SetPoint(i_points_neg,x_pos_neg,y_pos_neg); i_points_neg++;}
        if (x_pos_neg > 0.3 && x_pos_neg < 1.2 && y_pos_neg > 0.0 && y_pos_neg < 2.0) 
            {tge_pT_mean_Gauss_pos ->SetPoint(i_points_pos,x_pos_neg,y_pos_neg); i_points_pos++;}

    }

    printf("i_points_neg: %d, i_points_pos: %d \n",i_points_neg,i_points_pos);

    //----negative half------------
    
    tge_pT_mean_Gauss_neg ->Fit("pol2","F","",-1.2,-0.3);

    TF1 *myfunc_neg = tge_pT_mean_Gauss_neg->GetFunction("pol2");
    Double_t par[6];
    for (Int_t i_par = 0; i_par < 3; i_par++)
    {
        par[i_par] = myfunc_neg->GetParameter(i_par); //3 parameters for negative half
    }

    //----positve half------------
    tge_pT_mean_Gauss_pos ->Fit("pol2","F","",0.3,1.2);
    TF1 *myfunc_pos = tge_pT_mean_Gauss_pos->GetFunction("pol2");
    
    for (Int_t i_par = 0; i_par < 3; i_par++)
    {
        par[i_par+3] = myfunc_pos->GetParameter(i_par); //3 parameters for positive half
    }
    
    can_TPC_vs_Kalman   ->cd();
    tge_pT_mean_Gauss_neg   ->SetMarkerStyle(9);
    tge_pT_mean_Gauss_neg   ->SetMarkerSize(1.8);
    tge_pT_mean_Gauss_neg   ->SetMarkerColor(kBlue);
    tge_pT_mean_Gauss_neg   ->SetLineColor(kBlue);
    tge_pT_mean_Gauss_neg   ->Draw("same P");

    tge_pT_mean_Gauss_pos   ->SetMarkerStyle(9);
    tge_pT_mean_Gauss_pos   ->SetMarkerSize(1.8);
    tge_pT_mean_Gauss_pos   ->SetMarkerColor(kRed);
    tge_pT_mean_Gauss_pos   ->SetLineColor(kRed);
    tge_pT_mean_Gauss_pos   ->Draw("same");

    //----save parameters------

    TGraph* tg_pol2_params = new TGraph();

    tg_pol2_params ->SetNameTitle("tg_pol2_params","tg_pol2_params");

    for (Int_t i_par = 0; i_par < 6; i_par++)
    {
        tg_pol2_params ->SetPoint(i_par,i_par,par[i_par]);
        printf("par[%d]: %4.3f \n",i_par,par[i_par]);
    }

    TFile* pt_corr = new TFile("pt_corr.root","RECREATE");

    pt_corr ->cd();
    tg_pol2_params ->Write();

    //------------------------------------------------------



    TCanvas* can_proj = new TCanvas("can_proj","can_proj",10,10,1200,500);
    can_proj->Divide(3,1);
    //Int_t arr_proj_sel[6] = {67,74,82,63,109,119};
    Int_t arr_proj_sel[6] = {182,200,289,267,277,287};
    for(Int_t i_proj_plot = 0; i_proj_plot < 3; i_proj_plot++)
    {
        Int_t i_proj_sel = arr_proj_sel[i_proj_plot];

        can_proj->cd(i_proj_plot+1);
        can_proj->cd(i_proj_plot+1)->SetLeftMargin(0.22);
        can_proj->cd(i_proj_plot+1)->SetRightMargin(0.05);
        can_proj->cd(i_proj_plot+1)->SetBottomMargin(0.15);
        can_proj->cd(i_proj_plot+1)->SetTopMargin(0.08);
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
        //vec_h_proj_x[3*0 +0][i_proj_sel]->GetYaxis()->SetRangeUser(0.0,1.1);
        if(i_proj_plot == 1) vec_h_proj_x[i_proj_sel]->Rebin(4);
        if(i_proj_plot == 2) vec_h_proj_x[i_proj_sel]->Rebin(8);
        vec_h_proj_x[i_proj_sel]->GetXaxis()->SetRangeUser(-0.7,3.6);
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
        Double_t sigma     = vec_h_proj_x[i_proj_sel]->GetRMS()*1.0;
        func_Gauss_fit ->SetParameter(0,amplitude);
        func_Gauss_fit ->SetParameter(1,mean);
        func_Gauss_fit ->SetParameter(2,sigma);

        //printf("mean: %4.3f, sigma: %4.3f, amplitude: %4.3f \n",mean,sigma,amplitude);

        //Fit
        vec_h_proj_x[i_proj_sel]->Fit("func_Gauss_fit","QMN","",mean-2.0*sigma,mean+2.0*sigma);

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
        vec_h_proj_x[i_proj_sel] ->Fit("func_Gauss_fit","QMN","",mean-1.2*sigma,mean+1.2*sigma);

        // Get parameters
        amplitude = func_Gauss_fit ->GetParameter(0);
        mean      = func_Gauss_fit ->GetParameter(1);
        sigma     = fabs(func_Gauss_fit ->GetParameter(2));

        func_Gauss_fit ->SetRange(mean-1.2*sigma,mean+1.2*sigma);

        func_Gauss_fit ->SetLineColor(kRed);
        func_Gauss_fit ->SetLineWidth(3);
        //	func_Gauss_fit ->SetRange(0.,2.);

        func_Gauss_fit ->DrawCopy("same");

        Double_t pT_TPC = vec_pT_TPC[i_proj_sel];
        Double_t pT_res = sigma/fabs(pT_TPC);

        Double_t x_offset = 0.0;
        if(i_proj_plot == 0) x_offset = 0.2;
        if(i_proj_plot == 1) x_offset = 0.24;
        if(i_proj_plot == 2) x_offset = -0.05;

        HistName = "p_{T}^{TPC} = ";
        sprintf(NoP,"%4.1f",pT_TPC);
        HistName += NoP;
        HistName += " GeV/c";
        plotTopLegend((char*)HistName.Data(),0.3+x_offset,0.86,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        HistName = "#sigma = ";
        sprintf(NoP,"%4.3f",sigma);
        HistName += NoP;
        HistName += " GeV";
        plotTopLegend((char*)HistName.Data(),0.3+x_offset,0.8,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

        HistName = "R = ";
        sprintf(NoP,"%4.1f",pT_res*100.0);
        HistName += NoP;
        HistName += "%";
        plotTopLegend((char*)HistName.Data(),0.3+x_offset,0.74,0.055,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1

    }

    can_proj->cd(1);
    HistName = "p-Pb, #sqrt{s_{NN}}=5.02 TeV, N_{trk}= 6";
    plotTopLegend((char*)HistName.Data(),0.26,0.95,0.045,kBlack,0.0,42,1,1); // char* label,Float_t x=-1,Float_t y=-1, Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1, Int_t align = 1
    //----------------------------------------------------


}