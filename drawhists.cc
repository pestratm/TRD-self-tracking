#include "TRD_ST_Analyze_tracklets.h"
#include "TRD_Kalman_Tracking.h"
#include "TRD_Kalman_Tracking.cxx"
Bool_t fitting_track(Ali_TRD_ST_Tracklets* a,Ali_TRD_ST_Tracklets* b){
    //direction is not ok

    TVector3 a_angle = (TVector3){a->get_TV3_dir()[0], a->get_TV3_dir()[1], 0};
    TVector3 b_angle = (TVector3){b->get_TV3_dir()[0], b->get_TV3_dir()[1], 0};

    if (abs(a_angle.Angle(b_angle)*TMath::RadToDeg() )>7)
        return 0;
    if (abs(a->get_TV3_dir().Angle(b->get_TV3_dir())*TMath::RadToDeg() )>30)
        return 0;
    //position is not ok
    Double_t cylind_rad_dist=sqrt(b->get_TV3_offset()[0]*b->get_TV3_offset()[0] +b->get_TV3_offset()[1]*b->get_TV3_offset()[1]) -sqrt(a->get_TV3_offset()[0]*a->get_TV3_offset()[0] +a->get_TV3_offset()[1]*a->get_TV3_offset()[1]);
    Double_t cylind_rad_len=sqrt(a->get_TV3_dir()[0]*a->get_TV3_dir()[0] +a->get_TV3_dir()[1]*a->get_TV3_dir()[1]);
    TVector3 z=a->get_TV3_offset()+a->get_TV3_dir()*(cylind_rad_dist/cylind_rad_len );
    /*if (abs((a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] ) -b->get_TV3_offset())[1])>1)
     return 0;
     if (abs((a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] ) -b->get_TV3_offset())[2])>4)
     return 0;*/

    if ((z-b->get_TV3_offset())[0]>2.)
        return 0;

    if ((z-b->get_TV3_offset())[1]>2.)
        return 0;
    if ((z-b->get_TV3_offset())[2]>10.)
        return 0;

    return 1;
}
void drawhists()
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(11);
    //gStyle->SetOptFit(1111);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFillColor(0);
    gStyle->SetPalette(27);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);



    TH1F *histo = new TH1F("histogram","efficiency Kalman Trackfinder",20,0,1.2);

    TH1F *histo2 = new TH1F("histogram","purity Kalman Trackfinder",20,0,1.2);


    printf("TRD_ST_Analyze_tracklets started \n");
    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze();
    TRD_ST_Analyze ->Init_tree("List_data_ADC.txt");
    //Long64_t event = 10;

    TRD_Kalman_Trackfinder kalid;
    for (Long64_t event = 1; event < 2 ; event++)
    {

        TRD_ST_Analyze ->Loop_event(event);
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;

        //TRD_ST_Analyze ->Draw_event(event);
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;
        TRD_ST_Analyze ->Do_TPC_TRD_matching(event,3.0,10.0);
        //TRD_ST_Analyze ->Do_TPC_TRD_matching_allEvents(3.0,10.0);
        vector< vector<Ali_TRD_ST_Tracklets*> > kalman_found_tracks = kalid.Kalman_Trackfind(TRD_ST_Analyze->Tracklets,TRD_ST_Analyze->Number_Tracklets);
        TRD_ST_Analyze ->Draw_Kalman_Tracks(kalman_found_tracks);
        vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks=TRD_ST_Analyze->matched_tracks;
        vector< vector<Ali_TRD_ST_Tracklets*> > matched_beautiful_tracks;

        vector<vector<Double_t>> mHelices_kalman = kalid.get_Kalman_helix_params();
        printf("size of mHelices_kalman: %d \n",(Int_t)mHelices_kalman.size());
        TRD_ST_Analyze ->set_Kalman_helix_params(mHelices_kalman);
        TRD_ST_Analyze ->Draw_Kalman_Helix_Tracks();
        TRD_ST_Analyze ->Calculate_secondary_vertices(1); // 0 = no graphics
        //vector< vector<Ali_TRD_ST_Tracklets*> > kalman_found_tracks=kalid.found_tracks;
        //vector<Double_t> track_accuracy;

        //Ali_TRD_ST_Tracklets* last_tracklet;
        /*
         for(Int_t i_track=0;i_track< matched_tracks.size();i_track++){
         Bool_t beautiful=0;
         last_tracklet=NULL;

         for(Int_t i_layer=0;i_layer< matched_tracks[i_track].size();i_layer++){
         if 	(matched_tracks[i_track][i_layer]==NULL)
         continue;
         if (last_tracklet==NULL) last_tracklet=matched_tracks[i_track][i_layer];
         else{
         beautiful=fitting_track(last_tracklet,matched_tracks[i_track][i_layer]);
         if(!(beautiful))break;
         last_tracklet=matched_tracks[i_track][i_layer];
         }
         }
         if(beautiful) matched_beautiful_tracks.push_back(matched_tracks[i_track]);
         }
         cout<<"len_beatutiful:"<<matched_beautiful_tracks.size()<<endl;
         //TRD_ST_Analyze ->Draw_Kalman_Tracks(matched_tracks);
         */

        /*
        matched_beautiful_tracks=matched_tracks;
        TH1I* h_good_bad_TRD_chambers=TRD_ST_Analyze ->get_h_good_bad_TRD_chambers();;


        for(Int_t i_track=0;i_track< kalman_found_tracks.size();i_track++){

            Int_t number_to_find_real=0;
            Int_t number_found_max=0;
            Int_t number_of_noise_real=0;
            Int_t len_found_track_real=0;
            for(Int_t i_track_match=0;i_track_match< matched_beautiful_tracks.size();i_track_match++){

                Int_t number_to_find=0;
                Int_t number_found=0;
                Int_t number_of_noise=0;
                Int_t len_found_track=kalman_found_tracks[i_track].size();

                for (Int_t i_layer=0;i_layer< len_found_track;i_layer++){
                    if(kalman_found_tracks[i_track][i_layer]!=NULL)
                        number_of_noise++;
                    if (matched_beautiful_tracks[i_track_match][i_layer]==NULL) continue;
                    Int_t TRD_detector=matched_beautiful_tracks[i_track_match][i_layer]->get_TRD_det();
                    if (h_good_bad_TRD_chambers ->GetBinContent(TRD_detector+1)){
                        if(kalman_found_tracks[i_track][i_layer]!=NULL)
                            if(kalman_found_tracks[i_track][i_layer]->get_TRD_index()==matched_beautiful_tracks[i_track_match][i_layer]->get_TRD_index()){
                                number_to_find++;
                                number_found++;
                                number_of_noise--;
                            }
                        continue;
                    }
                    //cout<<"det:"<<(h_good_bad_TRD_chambers ->GetBinContent(TRD_detector+1))<<endl;
                    number_to_find++;
                    if(kalman_found_tracks[i_track][i_layer]==NULL) continue;

                    if(kalman_found_tracks[i_track][i_layer]->get_TRD_index()!=matched_beautiful_tracks[i_track_match][i_layer]->get_TRD_index())
                        continue;
                    //Double_t rel_off_diff=(kalman_found_tracks[i_track][i_layer]->get_TV3_offset()-matched_beautiful_tracks[i_track_match][i_layer]->get_TV3_offset()).Mag() /kalman_found_tracks[i_track][i_layer]->get_TV3_offset().Mag();
                    //if(rel_off_diff >0.001) continue;
                    //Double_t rel_dir_diff=(kalman_found_tracks[i_track][i_layer]->get_TV3_dir()-matched_beautiful_tracks[i_track_match][i_layer]->get_TV3_dir()).Mag() /kalman_found_tracks[i_track][i_layer]->get_TV3_dir().Mag();
                    //if(rel_dir_diff >0.001) continue;
                    number_found++;
                    number_of_noise--;

                }
                if(number_found >number_found_max){
                    number_found_max=number_found;
                    number_to_find_real=number_to_find;
                    number_of_noise_real=number_of_noise;
                    len_found_track_real=len_found_track;
                }
            }
            if(number_to_find_real>0) {
                histo->Fill((Double_t)number_found_max/number_to_find_real);
                histo2->Fill((1. -(Double_t)number_of_noise_real/len_found_track_real));
                cout<<(Double_t)number_found_max/number_to_find_real<<" "<<number_to_find_real<<" "<<i_track<<endl;}
            else
                cout<<"found_new_track"<<endl;
        }
        */
    }
    //gRandom->Rndm();


    TCanvas * c1= new TCanvas("c1", "fitted data",5,5,800,600);

    TColor *c3 = new TColor(9003,0,1,0);
    histo->SetLineColor(9003);

    histo->SetTitle("Efficiency Of Tracklet Detection");
    histo->SetXTitle("efficiency");
    histo->SetYTitle("count");

    histo->Draw("histo");
    TCanvas * c2= new TCanvas("c2", "fitted data",5,5,800,600);

    //TColor *c3 = new TColor(9003,0,1,0);


    histo2->SetLineColor(9003);

    histo2->SetTitle("Purity Of Tracklet Detection");
    histo2->SetXTitle("purity");
    histo2->SetYTitle("count");

    histo2->Draw("histo");


    //TRD_ST_Analyze ->Draw_Kalman_Tracks(kalid.found_tracks,kalid.nbr_tracks);

}