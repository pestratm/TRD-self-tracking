#include "TRD_TFTracker.h"
#include <TPython.h>


vector<vector<Double_t>> TRD_TFTrackMaker::get_Kalman_helix_params()
{
    return mHelices;
}


vector<vector<Ali_TRD_ST_Tracklets*>> TRD_TFTrackMaker::get_Tracklets()
{
    return mFound_tracklets;

}

Int_t TRD_TFTrackMaker::Trackfind(Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets){
    TTree tree;
    Int_t dets;
    Int_t tag;
    Float_t x1;
    Float_t x2;
    Float_t x3;
    Float_t v1;
    Float_t v2;
    Float_t v3;
    tree.Branch("x1",&x1);
    tree.Branch("x2",&x2);
    tree.Branch("x3",&x3);
    tree.Branch("v1",&v1);
    tree.Branch("v2",&v2);
    tree.Branch("v3",&v3);
    tree.Branch("det",&dets);
    tree.Branch("tag",&tag);

    for(Int_t i_tracklet=0;i_tracklet<Num_Tracklets;i_tracklet++) {
        if (Tracklets[i_tracklet]->get_TV3_offset().Mag() > 1000.0) continue;
        Int_t det =  Tracklets[i_tracklet] ->get_TRD_det();
        TVector3 offset = Tracklets[i_tracklet]->get_TV3_offset();
        TVector3 dir = Tracklets[i_tracklet]->get_TV3_dir();
        x1 = offset[0];
        x2 = offset[1];
        x3 = offset[2];
        v1 = dir[0];
        v2 = dir[1];
        v3 = dir[2];
        dets = det;
        tag = i_tracklet;
        tree.Fill();
    }
    
    vector<Double_t> fHelix;
    fHelix.resize(8);
    TVectorD fHelixS;
    TVectorD trcklts;
    TPython::Bind(&fHelixS,"fHelixS");
    TPython::Bind(&trcklts,"trcklts");
    TPython::Bind(&tree,"tree");
    TPython::Exec("inputs=tree.AsMatrix()");
    TPython::Exec("cls, results=worker.getdata(inputs)");
    int lr = TPython::Eval("len(results)");
    if(lr>0){
        TPython::Exec("fHelixS.ResizeTo(len(results))");
        TPython::Exec("fHelixS.SetElements(results)");
        TPython::Exec("trcklts.ResizeTo(len(cls))");
        TPython::Exec("cls = cls.astype(float)");
        TPython::Exec("trcklts.SetElements(cls)");
        Int_t nrows = fHelixS.GetNrows();
        if(nrows>0){
            Int_t ntracks=nrows/6;
            for(Int_t idx=0;idx<ntracks;idx++){
                Double_t x0  = fHelixS[6*idx+0];
                Double_t y0  = fHelixS[6*idx+1];
                Double_t z0  = fHelixS[6*idx+2];
                Double_t tgl = fHelixS[6*idx+3];
                Double_t C   = fHelixS[6*idx+4];
                Double_t phi0= fHelixS[6*idx+5];
                Double_t bfield=0.5;
                Double_t b_fak = bfield*3./2000.;
                Double_t pt  = b_fak/C;
                Double_t pz  = pt*tgl;
                fHelix[0] = y0;
                fHelix[1] = z0;
                fHelix[2] = phi0;
                fHelix[3] = tgl;
                fHelix[4] = C;
                fHelix[5] = x0;
                fHelix[6] = pt;
                fHelix[7] = pz;
                mHelices.push_back(fHelix);
            }


            Int_t NCls = trcklts.GetNrows();
            for(Int_t i = 0; i<NCls/6; i++){
                vector<Ali_TRD_ST_Tracklets*> track;
                for(Int_t layeri = 0; layeri <6; layeri++){
                    Int_t clsN =  (Int_t)trcklts[6*i+layeri];
                    track.push_back(Tracklets[clsN]);
                }
                mFound_tracklets.push_back(track);
            }
        }
    }
    return 1;
}
TRD_TFTrackMaker::TRD_TFTrackMaker(){
	cout << "loading Python interface" << endl;
	TPython::Exec("import numpy as np");
    	TPython::Exec("import tensorflow as tf");
	TPython::Exec("import  TFTracker");
	TPython::Exec("worker=TFTracker.trdpyinterface(\"tf_trd_tracker\")");
	TPython::Exec("worker.worker()");




}
