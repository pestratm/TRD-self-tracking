#ifndef ALI_AS_ANALYSIS_TRD_DIGITS_H
#define ALI_AS_ANALYSIS_TRD_DIGITS_H

class AliTRDdigitsManager;

#include "AliAnalysisTaskSE.h"

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "Ali_AS_Event.h"
#include "Ali_AS_EventLinkDef.h"

ClassImp(Ali_AS_TRD_digit)
ClassImp(Ali_AS_Track)
ClassImp(Ali_AS_Tracklet)
ClassImp(Ali_AS_offline_Tracklet)
ClassImp(Ali_AS_Event)

/*
class Ali_AS_analysis_TRD_digits : public AliAnalysisTaskSE
{
public:
    Ali_AS_analysis_TRD_digits()
	: AliAnalysisTaskSE(),
	fDigitsInputFileName("TRD.FltDigits.root"), fDigitsInputFile(0),
	fDigitsOutputFileName(""), fDigitsOutputFile(0),
	fDigMan(0),fGeo(0),AS_Event(0),AS_Track(0),AS_Tracklet(0),AS_Digit(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
	h_dca(0x0),h_dca_xyz(0x0),h2D_TPC_dEdx_vs_momentum(0x0),h_ADC_tracklet(0x0)
    {
	cout << "" << endl;
	cout << "***************************************************************************************" << endl;
	cout << "In Ali_AS_analysis_TRD_digits.h constructor" << endl;
	cout << "fDigitsInputFileName: " << fDigitsInputFileName << endl;
	cout << "***************************************************************************************" << endl;
	//AS_Event       = new Ali_AS_Event();
	//AS_Track       = new Ali_AS_Track();
	cout << "" << endl;
    }
	Ali_AS_analysis_TRD_digits(const char *name);
	//virtual ~Ali_AS_analysis_TRD_digits() {}

	virtual void   UserCreateOutputObjects();
	virtual Bool_t UserNotify();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
        void FillHelix(AliESDtrack* track_in, Double_t magF_in);
	void FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB);

	void SetDigitsInputFilename(TString x)
	{
	    fDigitsInputFileName=x;
	    cout << "" << endl;
	    cout << "***************************************************************************************" << endl;
	    cout << "fDigitsInputFileName: " << fDigitsInputFileName << endl;
	    cout << "***************************************************************************************" << endl;
	}
	void SetDigitsOutputFilename(TString x) {fDigitsOutputFileName=x;}

	AliHelix aliHelix;

    protected:

	Bool_t NextEvent(Bool_t preload=kFALSE);
	Bool_t ReadDigits();
	Bool_t WriteDigits();
        void   func_tail_cancellation(Short_t *arr, Int_t nexp);


	Int_t FindDigitsTrkl(const AliTRDtrackV1* trdTrack, Int_t layer,
			     Int_t* det, Int_t* row, Int_t* col,
			     Float_t* x, Float_t* y, Float_t* z);

	Int_t FindDigits(const AliExternalTrackParam* param,
			 Float_t bfield, Int_t layer,
			 Int_t* det, Int_t* row, Int_t* col);

	TList           *fListOfHistos;       //! list of output histograms
	TTree           *fTree;               //! output tree
	AliPIDResponse  *fPIDResponse;        //! PID handling
	AliESDtrackCuts *EsdTrackCuts;        //!

    private:

	TFile* OpenDigitsFile(TString inputfile, TString digfile, TString opt);

	TString fDigitsInputFileName;         //! Name of digits file for reading
	TFile*  fDigitsInputFile;             //! Digits file for reading
	TString fDigitsOutputFileName;        //! Name of digits file for writing
	TFile*  fDigitsOutputFile;            //! Digits file for writing

	AliTRDdigitsManager* fDigMan; //! digits manager
	AliTRDgeometry* fGeo; //! TRD geometry
	Ali_AS_Event* AS_Event;
        Ali_AS_Track* AS_Track;
        Ali_AS_Tracklet* AS_Tracklet;
        Ali_AS_TRD_digit* AS_Digit;
	TTree       *Tree_AS_Event;

	Int_t fEventNoInFile;
	Int_t N_good_events;
	Int_t fDigitsLoadedFlag;

	std::vector<TH1D*> h_dca;
	std::vector< std::vector<TH1D*> > h_dca_xyz;
        TH2D* h2D_TPC_dEdx_vs_momentum;
        vector<TH1D*> h_ADC_tracklet;

	Ali_AS_analysis_TRD_digits(const Ali_AS_analysis_TRD_digits&); // not implemented
	Ali_AS_analysis_TRD_digits& operator=(const Ali_AS_analysis_TRD_digits&); // not implemented

	ClassDef(Ali_AS_analysis_TRD_digits, 1);
};
*/


class Ali_AS_analysis_TRD_digits : public AliAnalysisTaskSE
{
public:
    Ali_AS_analysis_TRD_digits()
	: AliAnalysisTaskSE(),
	fDigitsInputFileName("TRD.FltDigits.root"), fDigitsInputFile(0),
	fDigitsOutputFileName(""), fDigitsOutputFile(0),
	fDigMan(0),fGeo(0),AS_Event(0),AS_Track(0),AS_Tracklet(0),AS_offline_Tracklet(0),AS_Digit(0),Tree_AS_Event(0), fEventNoInFile(-2), N_good_events(0), fDigitsLoadedFlag(kFALSE),
	h_dca(0x0),h_dca_xyz(0x0),h2D_TPC_dEdx_vs_momentum(0x0),h_ADC_tracklet(0x0),h_ADC_vs_time(0x0)
    {
	cout << "" << endl;
	cout << "***************************************************************************************" << endl;
	cout << "In Ali_AS_analysis_TRD_digits.h constructor" << endl;
	cout << "fDigitsInputFileName: " << fDigitsInputFileName << endl;
	cout << "***************************************************************************************" << endl;
	//AS_Event       = new Ali_AS_Event();
	//AS_Track       = new Ali_AS_Track();
	cout << "" << endl;
    }
	Ali_AS_analysis_TRD_digits(const char *name);
	//virtual ~Ali_AS_analysis_TRD_digits() {}

	virtual void   UserCreateOutputObjects();
	virtual Bool_t UserNotify();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
        void FillHelix(AliESDtrack* track_in, Double_t magF_in);
	void FindDCAHelixPoint(TVector3 space_vec, AliHelix helixA, Float_t path_initA, Float_t path_initB, Float_t &pathA, Float_t &dcaAB);

	void SetDigitsInputFilename(TString x)
	{
	    fDigitsInputFileName=x;
	    cout << "" << endl;
	    cout << "***************************************************************************************" << endl;
	    cout << "fDigitsInputFileName: " << fDigitsInputFileName << endl;
	    cout << "***************************************************************************************" << endl;
	}
	void SetDigitsOutputFilename(TString x) {fDigitsOutputFileName=x;}

	AliHelix aliHelix;

    protected:

	Bool_t NextEvent(Bool_t preload=kFALSE);
	Bool_t ReadDigits();
	Bool_t WriteDigits();
        void   func_tail_cancellation(Short_t *arr, Int_t nexp);


	Int_t FindDigitsTrkl(const AliTRDtrackV1* trdTrack, Int_t layer,
			     Int_t* det, Int_t* row, Int_t* col,
			     Float_t* x, Float_t* y, Float_t* z);

	Int_t FindDigits(const AliExternalTrackParam* param,
			 Float_t bfield, Int_t layer,
			 Int_t* det, Int_t* row, Int_t* col);

	TList           *fListOfHistos;       //! list of output histograms
	TTree           *fTree;               //! output tree
	AliPIDResponse  *fPIDResponse;        //! PID handling
	AliESDtrackCuts *EsdTrackCuts;        //!

    private:

	TFile* OpenDigitsFile(TString inputfile, TString digfile, TString opt);

	TString fDigitsInputFileName;         //! Name of digits file for reading
	TFile*  fDigitsInputFile;             //! Digits file for reading
	TString fDigitsOutputFileName;        //! Name of digits file for writing
	TFile*  fDigitsOutputFile;            //! Digits file for writing

	AliTRDdigitsManager* fDigMan; //! digits manager
	AliTRDgeometry* fGeo; //! TRD geometry
	Ali_AS_Event* AS_Event;
        Ali_AS_Track* AS_Track;
        Ali_AS_Tracklet* AS_Tracklet;
        Ali_AS_offline_Tracklet* AS_offline_Tracklet;
        Ali_AS_TRD_digit* AS_Digit;
	TTree       *Tree_AS_Event;

	Int_t fEventNoInFile;
	Int_t N_good_events;
	Int_t fDigitsLoadedFlag;

	std::vector<TH1D*> h_dca;
	std::vector< std::vector<TH1D*> > h_dca_xyz;
        TH2D* h2D_TPC_dEdx_vs_momentum;
        vector<TH1D*> h_ADC_tracklet;
        vector<TProfile*> h_ADC_vs_time;

	Ali_AS_analysis_TRD_digits(const Ali_AS_analysis_TRD_digits&); // not implemented
	Ali_AS_analysis_TRD_digits& operator=(const Ali_AS_analysis_TRD_digits&); // not implemented

	ClassDef(Ali_AS_analysis_TRD_digits, 1);
};



//----------------------------------------------------------------------------------------
class Class_peak_finder
{
    TH1D* TH1D_in;
    Double_t search_width;
    Double_t search_min_diff; // minimum difference between signal and background
    Double_t exclude_width; // range to be excluded between two peaks
    Double_t bin_width;
    Int_t    N_maxima;
    Double_t min_radius;
    Double_t max_width_scan;
    Int_t    flag_peak_rim; // 0 = peak finding, 1 = rim finding
    Int_t    flag_debug = 0; // 1 = return debug output
    TPolyMarker* pm_peaks;
    std::vector< std::vector<Double_t> > vec_peak_positions;
public:
    void  add_TH1D(TH1D* TH1D_in_a);
    void  set_finding_parameters(Int_t flag_peak_rim_in, Double_t search_width_in,
                                 Double_t search_min_diff_in, Double_t exclude_width_in,
                                 Int_t N_maxima_in, Double_t min_radius_in, Double_t max_width_scan_in);
    void  find_peaks();
    TPolyMarker* get_polymarker();
    std::vector< std::vector<Double_t> > get_peak_positions();
    void set_debug(Int_t flag_debug_in);
    void  clear();
};

void Class_peak_finder::set_finding_parameters(Int_t flag_peak_rim_in, Double_t search_width_in,
                                               Double_t search_min_diff_in, Double_t exclude_width_in,
                                               Int_t N_maxima_in, Double_t min_radius_in, Double_t max_width_scan_in)
{
    flag_peak_rim   = flag_peak_rim_in;
    search_width    = search_width_in;
    search_min_diff = search_min_diff_in;
    exclude_width   = exclude_width_in;
    N_maxima        = N_maxima_in;
    min_radius      = min_radius_in;
    max_width_scan  = max_width_scan_in;
}

void Class_peak_finder::add_TH1D(TH1D* TH1D_in_a)
{
    TH1D_in = (TH1D*)TH1D_in_a->Clone("TH1D_in");
    bin_width = TH1D_in->GetBinWidth(1);
    pm_peaks = new TPolyMarker();
}

void Class_peak_finder::find_peaks()
{
    Int_t max_bin = TH1D_in->GetNbinsX();
    Double_t min_val_x = TH1D_in->GetBinCenter(1);
    Double_t max_val_x = TH1D_in->GetBinCenter(max_bin);

    std::vector< std::vector<Double_t> > vec_diff;
    vec_diff.resize(8);

    vec_peak_positions.resize(4); // radial distance, height, density, width

    // Loop over all bins
    for(Int_t i_bin = 1; i_bin <= max_bin; i_bin++)
    {
        Double_t center_val = TH1D_in->GetBinCenter(i_bin);
        Double_t signal_height = TH1D_in->GetBinContent(i_bin);
        Double_t hist_start = TH1D_in->GetBinCenter(1);

        if(center_val < min_radius) continue;

        // Determine search width with maximum density = number of points weighted with Gaussian / Gaussian area
        Double_t max_signal         = 0.0;
        Double_t max_width          = 0.0;
        Double_t max_signal_density = 0.0;
        Double_t max_veto           = 1.0;
        Double_t max_bkgr_left      = 0.0;
        Double_t max_bkgr_right     = 0.0;
        for(Int_t i_search = 0; i_search < 20; i_search++)
        {
            Double_t search_width_scan          = search_width + search_width*i_search;
            Int_t    search_width_scan_bin_half = TH1D_in->FindBin(hist_start + search_width_scan/2.0);

            if(search_width_scan > max_width_scan)
            {
                //printf("i_bin: %d, i_search: %d, search_width_scan: %f, max_width_scan: %f \n",i_bin,i_search,search_width_scan,max_width_scan);
                break;
            }

            Int_t bin_left_bkgr    = i_bin - 2*search_width_scan_bin_half;
            Int_t bin_left_signal  = i_bin - search_width_scan_bin_half;
            Int_t bin_right_signal = i_bin + search_width_scan_bin_half;
            Int_t bin_right_bkgr   = i_bin + 2*search_width_scan_bin_half;

            //printf("i_bin: %d, i_search: %d, bin_width: %d, bins: {%d,%d,%d,%d} \n",i_bin,i_search,search_width_scan_bin_half,bin_left_bkgr,bin_left_signal,bin_right_signal,bin_right_bkgr);

            if(bin_left_bkgr < 1) break;
            if(bin_right_bkgr > max_bin) break;

            // Find minimum value in search range
            Double_t min_val_scan_range = 0.0;
            for(Int_t i_bin_search_min = bin_left_bkgr; i_bin_search_min <= bin_right_bkgr; i_bin_search_min++)
            {
                Double_t bin_cont = TH1D_in->GetBinContent(i_bin_search_min);
                if(bin_cont < min_val_scan_range) min_val_scan_range = bin_cont;
            }

            Double_t sum_signal     = 0.0;
            Double_t sum_bkgr       = 0.0;
            Double_t sum_bkgr_left  = 0.0;
            Double_t sum_bkgr_right = 0.0;
            for(Int_t i_bin_search_min = bin_left_bkgr; i_bin_search_min <= bin_right_bkgr; i_bin_search_min++)
            {
                Double_t bin_cont   = TH1D_in->GetBinContent(i_bin_search_min) - min_val_scan_range;
                if(i_bin_search_min >= bin_left_signal && i_bin_search_min <= bin_right_signal) sum_signal     += bin_cont;
                if(i_bin_search_min >= bin_left_bkgr && i_bin_search_min < bin_left_signal)     sum_bkgr_left  += bin_cont;
                if(i_bin_search_min > bin_right_signal && i_bin_search_min <= bin_right_bkgr)   sum_bkgr_right += bin_cont;
            }
            sum_bkgr = sum_bkgr_left + sum_bkgr_right;

            if(sum_bkgr > sum_signal) continue;
            if(sum_bkgr_left > sum_signal/2.0) continue;
            if(sum_bkgr_right > sum_signal/2.0) continue;
            if(sum_bkgr == 0.0) continue;

            Double_t signal = sum_signal - sum_bkgr;
            Double_t signal_density = signal/search_width_scan;

            //printf("i_bin: %d, i_search: %d, center_val: %f, signal: %f, signal_density: %f,  width: %f \n",i_bin,i_search,center_val,signal,signal_density,search_width_scan);
     

            if(signal_density > max_signal_density)
            {
                max_signal         = signal;
                max_signal_density = signal_density;
                max_width          = search_width_scan;
                max_veto           = 0.0;
                max_bkgr_left      = sum_bkgr_left;
                max_bkgr_right     = sum_bkgr_right;
                //printf("i_bin: %d, i_search: %d, center_val: %f, max_signal: %f, max_signal_density: %f,  max_width: %f \n",i_bin,i_search,center_val,max_signal,max_signal_density,max_width);
            }
        }

        if(max_width <= 0.0) continue;

        vec_diff[0].push_back(center_val);
        vec_diff[1].push_back(signal_height);
        vec_diff[2].push_back(max_signal_density);
        vec_diff[3].push_back(max_signal);
        vec_diff[4].push_back(max_veto);
        vec_diff[5].push_back(max_bkgr_left);
        vec_diff[6].push_back(max_bkgr_right);
        vec_diff[7].push_back(max_width);

        //if(flag_debug) printf("i_bin: %d, center_val: %f, S/B: %f, peak_total: %f, bkgr_left: %f, bkgr_right: %f \n",i_bin,center_val,diff_peak_bkgr,peak_total,bkgr_left,bkgr_right);
    } // end of histogram bin loop

    // Find all maxima, exclude ranges within exclude_width after each loop
    Int_t max_i_diff = 0;
    while(max_i_diff >= 0 && vec_diff[0].size() > 0)
    {
        Double_t best_center_val     = 0.0;
        Double_t best_signal_height  = 0.0;
        Double_t best_signal_density = 0.0;
        Double_t best_signal         = 0.0;
        Double_t best_bkgr_left      = 0.0;
        Double_t best_bkgr_right     = 0.0;
        Double_t best_width          = 0.0;

        max_i_diff = -1;
        if(flag_peak_rim == 0) // peak finding
        {
            // Find maximum signal to background above minimum defined
            for(Int_t i_diff = 0; i_diff < vec_diff[0].size(); i_diff++)
            {
                //printf("   i_diff: %d, S/B: %f, veto: %f \n",i_diff,vec_diff[1][i_diff],vec_diff[3][i_diff]);
                if(vec_diff[2][i_diff] > best_signal_density
                   //&& vec_diff[1][i_diff] >= search_min_diff
                   && vec_diff[4][i_diff] < 1.0)
                {
                    // Check if this point with its width is already in range of a previously found point
                    // Avoid peaks on top of peaks
                    Int_t flag_in_range_of_accepted_point = 0;
                    for(Int_t i_already_accepted_point = 0; i_already_accepted_point < vec_peak_positions[0].size(); i_already_accepted_point++)
                    {
                        Double_t acc_center       = vec_peak_positions[0][i_already_accepted_point];
                        Double_t this_point_width = vec_diff[7][i_diff];

                        if(fabs(vec_diff[0][i_diff] - acc_center) < this_point_width)
                        {
                            flag_in_range_of_accepted_point = 1;
                            vec_diff[4][i_diff] = 1.0; // veto this point
                        }
                    }

                    if(!flag_in_range_of_accepted_point)
                    {
                        best_center_val     = vec_diff[0][i_diff];
                        best_signal_height  = vec_diff[1][i_diff];
                        best_signal_density = vec_diff[2][i_diff];
                        best_signal         = vec_diff[3][i_diff];
                        best_bkgr_left      = vec_diff[5][i_diff];
                        best_bkgr_right     = vec_diff[6][i_diff];
                        best_width          = vec_diff[7][i_diff];
                        max_i_diff          = i_diff;
                        //if(flag_debug) printf("    -> peak found at x_pos: %f",max_diff_x);
                    }
                }
            }
        }
        else // rim finding
        {
            /*
            // Find maximum signal to background above minimum defined
            for(Int_t i_diff = 0; i_diff < vec_diff[0].size(); i_diff++)
            {
                Double_t slope_left  = vec_diff[2][i_diff] - vec_diff[4][i_diff];
                Double_t slope_right = vec_diff[5][i_diff] - vec_diff[2][i_diff];

                Double_t slope_change = fabs(slope_right - slope_left);
                if(1)
                {
                    //printf("max_i_diff: %d, heights: {%f,%f,%f}, slopes: {%f,%f}, slope_change: %f \n",max_i_diff,vec_diff[4][i_diff],vec_diff[2][i_diff],vec_diff[5][i_diff],slope_left,slope_right,slope_change);
                    if(slope_change > max_diff &&  slope_change > search_min_diff)
                    {
                        max_diff_x = vec_diff[0][i_diff];
                        max_diff   = slope_change;
                        max_height = vec_diff[2][i_diff];
                        max_i_diff = i_diff;
                    }
                }
            }
            */
        }

        if(best_width > exclude_width) exclude_width = best_width;

        // Find range within exclude_width to be erased
        Int_t i_diff_left_erase = 0;
        for(Int_t i_diff = max_i_diff; i_diff >= 0; i_diff--)
        {
            if(fabs(best_center_val - vec_diff[0][i_diff]) > exclude_width)
            {
                i_diff_left_erase = i_diff;
                break;
            }
        }

        Int_t i_diff_right_erase = vec_diff[0].size();
        for(Int_t i_diff = max_i_diff; i_diff <= vec_diff[0].size(); i_diff++)
        {
            if(fabs(best_center_val - vec_diff[0][i_diff]) > exclude_width)
            {
                i_diff_right_erase = i_diff;
                break;
            }
        }


        //cout << "start erase: "<<  i_diff_left_erase << ", stop erase: " << i_diff_right_erase << endl;
      
        if(max_i_diff >= 0)
        {
            //printf("polymarker set at: %f, height: %f \n",best_center_val,best_signal_height);
            pm_peaks ->SetNextPoint(best_center_val,best_signal_height);
            vec_peak_positions[0].push_back(best_center_val);
            vec_peak_positions[1].push_back(best_signal_height);
            vec_peak_positions[2].push_back(best_signal_density);
            vec_peak_positions[3].push_back(best_width);

            if(flag_debug) printf("max_i_diff: %d, x,height: {%f,%f}, best_density: %f, best_width: %f, S,Lbkgr,Rbkgr: {%f,%f,%f} \n",max_i_diff,best_center_val,best_signal_height,best_signal_density,best_width,best_signal,best_bkgr_left,best_bkgr_right);

            //cout << "max_diff_x: " << max_diff_x << ", max_diff: " << max_diff << ", max_i_diff: " << max_i_diff << ", size: " << vec_diff[0].size()
            //    << ", peak_total: " << vec_diff[2][max_i_diff] << ", bkgr_left: " << vec_diff[4][max_i_diff] << ", bkgr_right: " << vec_diff[5][max_i_diff] << endl;

        }

        if(vec_peak_positions[0].size() >= N_maxima) break;

        vec_diff[0].erase(vec_diff[0].begin() + i_diff_left_erase, vec_diff[0].begin() + i_diff_right_erase);
        vec_diff[1].erase(vec_diff[1].begin() + i_diff_left_erase, vec_diff[1].begin() + i_diff_right_erase);
        vec_diff[2].erase(vec_diff[2].begin() + i_diff_left_erase, vec_diff[2].begin() + i_diff_right_erase);
        vec_diff[3].erase(vec_diff[3].begin() + i_diff_left_erase, vec_diff[3].begin() + i_diff_right_erase);
        vec_diff[4].erase(vec_diff[4].begin() + i_diff_left_erase, vec_diff[4].begin() + i_diff_right_erase);
        vec_diff[5].erase(vec_diff[5].begin() + i_diff_left_erase, vec_diff[5].begin() + i_diff_right_erase);
        vec_diff[6].erase(vec_diff[6].begin() + i_diff_left_erase, vec_diff[6].begin() + i_diff_right_erase);
        vec_diff[7].erase(vec_diff[7].begin() + i_diff_left_erase, vec_diff[7].begin() + i_diff_right_erase);
    }
}

TPolyMarker* Class_peak_finder::get_polymarker()
{
    return pm_peaks;
}

std::vector< std::vector<Double_t> > Class_peak_finder::get_peak_positions()
{
    return vec_peak_positions;
}

void Class_peak_finder::set_debug(Int_t flag_debug_in)
{
    flag_debug = flag_debug_in;
}

void Class_peak_finder::clear()
{
    delete pm_peaks;
    vec_peak_positions.clear();
    delete TH1D_in;
}
//----------------------------------------------------------------------------------------




#endif
