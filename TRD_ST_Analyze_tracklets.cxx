
#include "TRD_ST_Analyze_tracklets.h"


//----------------------------------------------------------------------------------------
Ali_TRD_ST_Analyze::Ali_TRD_ST_Analyze(TString out_dir, TString out_file_name, Int_t graphics)
{
    //layer_radii_file = TFile::Open("./TRD_layer_radii.root");
    //h_layer_radii_det = (TH1D*)layer_radii_file ->Get("h_layer_radii_det");

    TPC_single_helix = new Ali_Helix();
    vec_h2D_pT_vs_TPC_TRD_residuals.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_h2D_pT_vs_TPC_TRD_residuals[i_det].resize(3); // x,y,z
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            HistName = "vec_h2D_pT_vs_TPC_TRD_residuals_det";
            HistName += i_det;
            HistName += "_xyz";
            HistName += i_xyz;
            vec_h2D_pT_vs_TPC_TRD_residuals[i_det][i_xyz] = new TH2D(HistName.Data(),HistName.Data(),200,-10.0,10.0,10,-5,5);
        }
    }

    TFile* file_TRD_geometry = TFile::Open("./TRD_Geometry.root");
    TList* list_geom = (TList*)file_TRD_geometry->Get("TRD_Digits_output");
    h2D_TRD_det_coordinates = (TH2D*)list_geom ->FindObject("h2D_TRD_det_coordinates");
    h_layer_radii_det = new TH1D("h_layer_radii_det","h_layer_radii_det",540,0,540);

    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        Double_t x_val = h2D_TRD_det_coordinates ->GetBinContent(i_det+1,1);
        Double_t y_val = h2D_TRD_det_coordinates ->GetBinContent(i_det+1,2);
        Double_t z_val = h2D_TRD_det_coordinates ->GetBinContent(i_det+1,3);

        Double_t radius_align = TMath::Sqrt(x_val*x_val + y_val*y_val);
        h_layer_radii_det ->SetBinContent(i_det+1,radius_align);

    }


    HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    outputfile = new TFile(HistName.Data(),"RECREATE");
    //------------------------------------------------
    outputfile ->cd();
    // TRD self tracking output data containers
    TRD_ST_Tracklet_out   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track_out  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event_out      = new Ali_TRD_ST_Event();

    Tree_TRD_ST_Event_out  = NULL;
    Tree_TRD_ST_Event_out  = new TTree("Tree_TRD_ST_Event_out" , "TRD_ST_Events_out" );
    Tree_TRD_ST_Event_out  ->Branch("Tree_TRD_ST_Event_branch_out"  , "TRD_ST_Event_out", TRD_ST_Event_out );
    Tree_TRD_ST_Event_out  ->SetAutoSave( 5000000 );
    //------------------------------------------------


    // constructor
#if defined(USEEVE)
    if(graphics)
    {
        TEveManager::Create();

        TEveP_primary_vertex = new TEvePointSet();

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
        vec_TEveLine_tracklets_match.resize(6); // layers
    }
#endif

    th1d_TRD_layer_radii = new TH1D("th1d_TRD_layer_radii","th1d_TRD_layer_radii",900,250,400.0);
    vec_th1d_TRD_layer_radii_det.resize(540);
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        HistName = "vec_th1d_TRD_layer_radii_det_";
        HistName += TRD_detector;
        vec_th1d_TRD_layer_radii_det[TRD_detector] = new TH1D(HistName.Data(),HistName.Data(),900,250,400.0);
    }




    //--------------------------
    // Open histogram which defines good and bad chambers
    //TFile* file_TRD_QA = TFile::Open("./Data/chamber_QC.root");
    //h_good_bad_TRD_chambers = (TH1D*)file_TRD_QA ->Get("all_defects_hist");

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
    //--------------------------

    //--------------------------
    // Load TRD geometry
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
#if defined(USEEVE)
    vec_eve_TRD_detector_box.resize(540);
#endif
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
    //= {kCyan,kPink,kMagenta,kMagenta+2,kOrange,kOrange+2,kRed,kRed+2};

#if defined(USEEVE)
    if(graphics)
    {
        for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
        {
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
            }

            gEve->AddElement(vec_eve_TRD_detector_box[TRD_detector]);
        }
    }
#endif
    //--------------------------


    TH2D_AP_plot          = new TH2D("TH2D_AP_plot","TH2D_AP_plot",200,-2.0,2.0,400,-0.1,4.0);
    TH2D_pT_TPC_vs_Kalman = new TH2D("TH2D_pT_TPC_vs_Kalman","TH2D_pT_TPC_vs_Kalman",2000,-10.0,10.0,2000,-10.0,10.0);

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::create_output_file(TString out_dir, TString out_file_name)
{
    HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    outputfile = new TFile(HistName.Data(),"RECREATE");
    //------------------------------------------------
    outputfile ->cd();
    // TRD self tracking output data containers
    TRD_ST_Tracklet_out   = new Ali_TRD_ST_Tracklets();
    TRD_ST_TPC_Track_out  = new Ali_TRD_ST_TPC_Track();
    TRD_ST_Event_out      = new Ali_TRD_ST_Event();

    Tree_TRD_ST_Event_out  = NULL;
    Tree_TRD_ST_Event_out  = new TTree("Tree_TRD_ST_Event_out" , "TRD_ST_Events_out" );
    Tree_TRD_ST_Event_out  ->Branch("Tree_TRD_ST_Event_branch_out"  , "TRD_ST_Event_out", TRD_ST_Event_out );
    Tree_TRD_ST_Event_out  ->SetAutoSave( 5000000 );
    //------------------------------------------------
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::fHelixAtoPointdca(TVector3 space_vec, Ali_Helix* helixA, Float_t &pathA, Float_t &dcaAB)
{
    // V1.1
    Float_t pA[2] = {100.0,-100.0}; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    Float_t distarray[2];
    TVector3 testA;
    Double_t helix_point[3];
    for(Int_t r = 0; r < 2; r++)
    {
        helixA ->Evaluate(pA[r],helix_point);  // 3D-vector of helixB point at path pB[r]
        testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        distarray[r] = (testA-space_vec).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.01 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pA[0] = " << pA[0]
        //    << ", pA[1] = " << pA[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pA[1]-pA[0])*scale; // the next length interval
            pA[0]     = pA[1] + scale_length; // the new path

            helixA ->Evaluate(pA[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            distarray[0] = (testA-space_vec).Mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pA[0]-pA[1])*scale;
            pA[1]     = pA[0] + scale_length;
            helixA ->Evaluate(pA[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[1] = (testA-space_vec).Mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathA = pA[1];
        dcaAB = distarray[1];
    }
    //cout << "pathA = " << pathA << ", dcaAB = " << dcaAB << endl;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,
                           Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
    // (x1,y1) -> center of circle 1, r1 -> radius of circle 1
    // (x2,y2) -> center of circle 2, r2 -> radius of circle 2
    // (x1_c,y1_c) -> crossing point 1
    // (x2_c,y2_c) -> crossing point 2
    // Solution see Wikipedia

    Double_t s = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));  // distance between the center of the circles

    x1_c = 0;
    y1_c = 0;
    x2_c = 0;
    y2_c = 0;

    //printf("s: %4.3f, r1+r2: %4.3f, pos1: {%4.3f, %4.3f}, pos2: {%4.3f, %4.3f}, r1: %4.3f, r2: %4.3f \n",s,r1+r2,x1,y1,x2,y2,r1,r2);

    if(x1 != x2 && y1 != y2 && s < (r1+r2))
    {
        Double_t m  = (x1-x2)/(y2-y1);
        Double_t n  = (-r2*r2 + r1*r1 + y2*y2 - y1*y1 + x2*x2 - x1*x1)/(2.0*(y2-y1));
        Double_t p  = (2.0*(-x1 + m*(n-y1)))/(1.0 + m*m);
        Double_t q  = (x1*x1 + (n-y1)*(n-y1) -r1*r1)/(1.0 + m*m);
        Double_t p2 = (p/2.0)*(p/2.0);


        if(p2 >= q)
        {
            x1_c = (-p/2.0) + sqrt(p2 - q);
            x2_c = (-p/2.0) - sqrt(p2 - q);
            y1_c = m*x1_c + n;
            y2_c = m*x2_c + n;
            return 1;
        }
        else return 0;
    }
    else return 0;

} 
//----------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::fCircle_Interception(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2,
                                               Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c)
{
	Double_t dif_x = -(x1 - x2);
	Double_t dif_y = -(y1 - y2);
	Double_t dist=TMath::Sqrt( (dif_x*dif_x) + (dif_y*dif_y) );
        if(dist==0) // circles overlap -> most likely identical helices
        {
            x1_c = 0.0;
            y1_c = 0.0;
            x2_c = 0.0;
            y2_c = 0.0;
            return 0;
        };
        Double_t dist_inv =1/dist;
	
        if(dist<(r1+r2))
        {
	//2 intersections
		Double_t x_loc 		=((dist*dist +r1*r1 -r2*r2)*(0.5*dist_inv));
		Double_t y_loc		=TMath::Sqrt(r1*r1 -x_loc*x_loc);
		
                x1_c	=(x1*dist+x_loc*dif_x -y_loc*dif_y)*dist_inv;
                y1_c	=(y1*dist+x_loc*dif_y +y_loc*dif_x)*dist_inv;
		
                x2_c	=(x1*dist+x_loc*dif_x +y_loc*dif_y)*dist_inv;
                y2_c	=(y1*dist+x_loc*dif_y -y_loc*dif_x)*dist_inv;

                return 1;
	}
        else
        {
	//no intersection (maybe 1)
            Double_t normrad	=(r1+dist-r2)*dist_inv*0.5;
            
            //Double_t normrad	= (r1 + 0.5*(dist - r1 - r2))*dist_inv;
                x1_c		=x1 + normrad*dif_x;
                y1_c		=y1 + normrad*dif_y;
		
                x2_c		=x1_c;
                y2_c		=y1_c;
        
            	return 2;
	}	
	
	return 3;
}	
//------------------------------------------------------------------------------------------------------------



//---------------------------------------------------------------------------------------------------------
pair<Double_t,Double_t> Ali_TRD_ST_Analyze::fpathLength(Double_t r,Ali_Helix* helixA) const
{
	//taken from https://www.star.bnl.gov/webdata/dox/html/StHelix_8cc_source.html
	
	pair<Double_t,Double_t> value;
	pair<Double_t,Double_t> VALUE(999999999.,999999999.);
	Double_t curvature		= fabs(helixA->getHelix_param(4));
	Double_t radius			= fabs(1/curvature);
	Double_t x0				= helixA->getHelix_param(5) +radius*TMath::Sin(helixA->getHelix_param(2));
	Double_t y0				= helixA->getHelix_param(0) -radius*TMath::Cos(helixA->getHelix_param(2));
	Double_t z0				= helixA->getHelix_param(1);

	Double_t phase			= helixA->getHelix_param(2) -TMath::Pi()/2;
	//Double_t phase=atan2f(-(x0-helixA->getHelix_param(5)),(y0-helixA->getHelix_param(0)));
	
	while(phase<-TMath::Pi()) phase	+= 2.*TMath::Pi() ;
	while(phase>TMath::Pi()) phase	-= 2.*TMath::Pi() ;
	
	Double_t dipangle		= TMath::ATan(helixA->getHelix_param(3));
	Double_t cosdipangle	= TMath::Cos(dipangle);
	Double_t sindipangle	= TMath::Sin(dipangle);
	Double_t h				= TMath::Sign(1,helixA->getHelix_param(4));
	if(phase> TMath::Pi()) phase	-= 2.*TMath::Pi() ;
	Double_t cosphase		= TMath::Cos(phase);
	Double_t sinphase		= TMath::Sin(phase);

        //printf("h: %4.3f, phase: %4.3f, dipangle: %4.3f \n",h,phase,dipangle);

	Double_t t1 	= y0*curvature;
	Double_t t2 	= sinphase;
  	Double_t t3 	= curvature*curvature;
  	Double_t t4 	= y0*t2;
	Double_t t5 	= cosphase;
	Double_t t6 	= x0*t5;
	Double_t t8 	= x0*x0;
	Double_t t11	= y0*y0;
	Double_t t14 	= r*r;
	Double_t t15 	= t14*curvature;
	Double_t t17 	= t8*t8;
	Double_t t19 	= t11*t11;
	Double_t t21 	= t11*t3;
	Double_t t23 	= t5*t5;
	Double_t t32 	= t14*t14;
	Double_t t35 	= t14*t3;
	Double_t t38 	= 8.0*t4*t6 - 4.0*t1*t2*t8 - 4.0*t11*curvature*t6 +
				4.0*t15*t6 + t17*t3 + t19*t3 + 2.0*t21*t8 + 4.0*t8*t23 -
				4.0*t8*x0*curvature*t5 - 4.0*t11*t23 -
				4.0*t11*y0*curvature*t2 + 4.0*t11 - 4.0*t14 +
				t32*t3 + 4.0*t15*t4 - 2.0*t35*t11 - 2.0*t35*t8;
	//cout<<t38<<endl;
	Double_t t40 	= (-t3*t38);
	if (t40<0.) return VALUE;
        t40 = ::sqrt(t40);

        if (h<0) phase += TMath::Pi();
	
	if(phase>TMath::Pi()) phase	-= 2.*TMath::Pi();

	Double_t t43 	= x0*curvature;
	Double_t t45 	= 2.0*t5 - t35 + t21 + 2.0 - 2.0*t1*t2 -2.0*t43 - 2.0*t43*t5 + t8*t3;
	Double_t t46 	= h*cosdipangle*curvature;

	value.first 	= (-phase + 2.0*atan((-2.0*t1 + 2.0*t2 + t40)/t45))/t46;
	value.second 	= -(phase + 2.0*atan((2.0*t1 - 2.0*t2 + t40)/t45))/t46;

	//
	//   Solution can be off by +/- one period, select smallest
	//
	Double_t p 		= fabs(2*TMath::Pi()/(h*curvature*cosdipangle));
	if (! std::isnan(value.first)) {
		if (fabs(value.first-p) < fabs(value.first)) value.first = value.first-p;
	   	else if (fabs(value.first+p) < fabs(value.first)) value.first = value.first+p;
	}
	if (! std::isnan(value.second)) {
	   	if (fabs(value.second-p) < fabs(value.second)) value.second = value.second-p;
	   	else if (fabs(value.second+p) < fabs(value.second)) value.second = value.second+p;
	}
	
	if (value.first > value.second)
	swap(value.first,value.second);
	return(value);

	
}
//----------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::fDCA_Helix_Estimate(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB)
{

    // Calculates the 2D crossing point, calculates the corresponding 3D point and returns pathA and pathB

    Double_t helix_point[3];

    Double_t x1 = helixA->getHelix_param(5);
    Double_t y1 = helixA->getHelix_param(0);
    Double_t x2 = helixB->getHelix_param(5);
    Double_t y2 = helixB->getHelix_param(0);
    Double_t c1 = helixA->getHelix_param(4);
    Double_t c2 = helixB->getHelix_param(4);
    Double_t r1 = 0.0;
    Double_t r2 = 0.0;
    if(c1 != 0 && c2 != 0)
    {
        r1 = fabs(1.0/c1);
        r2 = fabs(1.0/c2);
    } else return 0;

    Double_t x1_c = 0.0;
    Double_t y1_c = 0.0;
    Double_t x2_c = 0.0;
    Double_t y2_c = 0.0;

    //Int_t bCross_points = fCross_points_Circles(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);
    //printf("2D circle cross points (Alex): {%4.3f, %4.3f}, {%4.3f, %4.3f} \n",x1_c,y1_c,x2_c,y2_c);

    pathA = 0.0;
    pathB = 0.0;
    dcaAB = 0.0;


    Int_t cCross_points = fCircle_Interception(x1,y1,r1,x2,y2,r2,x1_c,y1_c,x2_c,y2_c);
    //printf("2D circle cross points (Sven): {%4.3f, %4.3f}, {%4.3f, %4.3f}, return: %d \n",x1_c,y1_c,x2_c,y2_c,cCross_points);

    Double_t radiusA = sqrt(x1_c*x1_c+y1_c*y1_c);
    Double_t radiusB = sqrt(x2_c*x2_c+y2_c*y2_c);
    //printf("radiusA: %4.3f, radiusB: %4.3f \n",radiusA,radiusB);

    //cout << "bCross_points = " << bCross_points << ", xyr(1) = {" << x1 << ", " << y1 << ", " << r1
    //    << "}, xyr(2) = {"  << x2 << ", " << y2 << ", " << r2 << "}, p1 = {" << x1_c << ", " << y1_c << "}, p2 = {" << x2_c << ", " << y2_c << "}" << endl;

    //if(bCross_points == 0) return 0;

    TVector3 pointA,pointB,pointA1,pointB1,pointA2,pointB2;

    Double_t path_lengthA_c1,path_lengthA_c2,path_lengthB_c1,path_lengthB_c2;

    // first crossing point for helix A
    pair< double, double > path_lengthA = fpathLength(radiusA,helixA);
    Double_t path_lengthA1 = path_lengthA.first;
    Double_t path_lengthA2 = path_lengthA.second;

    helixA ->Evaluate(path_lengthA1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x1_c-pointA1.x())*(x1_c-pointA1.x()) + (y1_c-pointA1.y())*(y1_c-pointA1.y())) <
       ((x1_c-pointA2.x())*(x1_c-pointA2.x()) + (y1_c-pointA2.y())*(y1_c-pointA2.y())))
    {
        path_lengthA_c1 = path_lengthA1;
    }
    else
    {
        path_lengthA_c1 = path_lengthA2;
    }

    // second crossing point for helix A
    path_lengthA = fpathLength(radiusB,helixA);
    path_lengthA1 = path_lengthA.first;
    path_lengthA2 = path_lengthA.second;

    helixA ->Evaluate(path_lengthA1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x2_c-pointA1.x())*(x2_c-pointA1.x()) + (y2_c-pointA1.y())*(y2_c-pointA1.y())) <
       ((x2_c-pointA2.x())*(x2_c-pointA2.x()) + (y2_c-pointA2.y())*(y2_c-pointA2.y())))
    {
        path_lengthA_c2 = path_lengthA1;
    }
    else
    {
        path_lengthA_c2 = path_lengthA2;
    }

    // first crossing point for helix B
    pair< double, double > path_lengthB = fpathLength(radiusA,helixB);
    Double_t path_lengthB1 = path_lengthB.first;
    Double_t path_lengthB2 = path_lengthB.second;

    helixB ->Evaluate(path_lengthB1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x1_c-pointB1.x())*(x1_c-pointB1.x()) + (y1_c-pointB1.y())*(y1_c-pointB1.y())) <
       ((x1_c-pointB2.x())*(x1_c-pointB2.x()) + (y1_c-pointB2.y())*(y1_c-pointB2.y())))
    {
        path_lengthB_c1 = path_lengthB1;
    }
    else
    {
        path_lengthB_c1 = path_lengthB2;
    }

    // second crossing point for helix B
    path_lengthB = fpathLength(radiusB,helixB);
    path_lengthB1 = path_lengthB.first;
    path_lengthB2 = path_lengthB.second;

    helixB ->Evaluate(path_lengthB1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    if( ((x2_c-pointB1.x())*(x2_c-pointB1.x()) + (y2_c-pointB1.y())*(y2_c-pointB1.y())) <
       ((x2_c-pointB2.x())*(x2_c-pointB2.x()) + (y2_c-pointB2.y())*(y2_c-pointB2.y())))
    {
        path_lengthB_c2 = path_lengthB1;
    }
    else
    {
        path_lengthB_c2 = path_lengthB2;
    }

    helixA ->Evaluate(path_lengthA_c1,helix_point);
    pointA1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixA ->Evaluate(path_lengthA_c2,helix_point);
    pointA2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

    helixB ->Evaluate(path_lengthB_c1,helix_point);
    pointB1.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
    helixB ->Evaluate(path_lengthB_c2,helix_point);
    pointB2.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);


    Double_t dcaAB1 = (pointA1-pointB1).Mag();
    Double_t dcaAB2 = (pointA2-pointB2).Mag();
    Double_t dcaAB3 = (pointA1-pointB2).Mag();
    Double_t dcaAB4 = (pointA2-pointB1).Mag();

#if 0
    printf("pointA1: {%4.3f, %4.3f, %4.3f} \n",pointA1.X(),pointA1.Y(),pointA1.Z());
    printf("pointA2: {%4.3f, %4.3f, %4.3f} \n",pointA2.X(),pointA2.Y(),pointA2.Z());
    printf("pointB1: {%4.3f, %4.3f, %4.3f} \n",pointB1.X(),pointB1.Y(),pointB1.Z());
    printf("pointB2: {%4.3f, %4.3f, %4.3f} \n",pointB2.X(),pointB2.Y(),pointB2.Z());
    printf("dcaAB1: %4.3f, dcaAB2: %4.3f, dcaAB3: %4.3f, dcaAB4: %4.3f \n",dcaAB1,dcaAB2,dcaAB3,dcaAB4);
#endif

    if(dcaAB1 < dcaAB2)
    {
        pathA = path_lengthA_c1;
        pathB = path_lengthB_c1;
        dcaAB = dcaAB1;
    }
    else
    {
        pathA = path_lengthA_c2;
        pathB = path_lengthB_c2;
        dcaAB = dcaAB2;
    }


    return 1;

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
// V1.0
void Ali_TRD_ST_Analyze::fHelixABdca(Ali_Helix* helixA, Ali_Helix* helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB,
                                    Float_t pathA_in, Float_t pathB_in)
{
    //cout << "Standard fHelixABdca called..." << endl;
    Float_t pA[2];
    Float_t pB[2]; // the two start values for pathB, 0.0 is the origin of the helix at the first measured point
    if(pathA_in < -9990.0 && pathB_in < -9990.0)
    {
        pA[0] = 0.0;
        pA[1] = 0.0;
        pB[0] = 0.0;
        pB[1] = -70.0;
    }
    else
    {
        pA[0] = pathA_in+5.0;
        pA[1] = pathA_in-5.0;
        pB[0] = pathB_in+5.0;
        pB[1] = pathB_in-5.0;
    }
    Float_t distarray[2];
    TVector3 testA, testB;
    Double_t helix_point[3];
    for(Int_t r = 0; r < 2; r++)
    {
        helixB ->Evaluate(pB[r],helix_point);  // 3D-vector of helixB point at path pB[r]
        testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        Float_t pathA_dca = -999.0;
        Float_t dcaAB_dca = -999.0;
        fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation

        helixA ->Evaluate(pathA_dca,helix_point);  // 3D-vector of helixB point at path pB[r]
        testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

        distarray[r] = (testA-testB).Mag(); // dca between helixA and helixB
    }
    Int_t loopcounter = 0;
    Float_t scale = 1.0;
    Float_t flip  = 1.0; // checks if the minimization direction changed
    Float_t scale_length = 100.0;
    while(fabs(scale_length) > 0.05 && loopcounter < 100) // stops when the length is too small
    {
        //cout << "n = " << loopcounter << ", pB[0] = " << pB[0]
        //    << ", pB[1] = " << pB[1] << ", d[0] = " << distarray[0]
        //    << ", d[1] = " << distarray[1] << ", flip = " << flip
        //    << ", scale_length = " << scale_length << endl;
        if(distarray[0] > distarray[1])
        {
            if(loopcounter != 0)
            {
                if(flip == 1.0) scale = 0.4; // if minimization direction changes -> go back, but only the way * 0.4
                else scale = 0.7; // go on in this direction but only by the way * 0.7
            }
            scale_length = (pB[1]-pB[0])*scale; // the next length interval
            pB[0]     = pB[1] + scale_length; // the new path
            helixB ->Evaluate(pB[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[0] = pathA_dca;
            //pA[0]     = helixA.pathLength(testB); // pathA at dca to helixB

            helixA ->Evaluate(pA[0],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[0] = (testA-testB).Mag(); // new dca
            flip = 1.0;
        }
        else
        {
            if(loopcounter != 0)
            {
                if(flip == -1.0) scale = 0.4;
                else scale = 0.7;
            }
            scale_length = (pB[0]-pB[1])*scale;
            pB[1]     = pB[0] + scale_length;
            helixB ->Evaluate(pB[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testB.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);

            Float_t pathA_dca = -999.0;
            Float_t dcaAB_dca = -999.0;
            fHelixAtoPointdca(testB,helixA,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            pA[1] = pathA_dca;
            //pA[1]     = helixA.pathLength(testB); // pathA at dca to helixB

            helixA ->Evaluate(pA[1],helix_point);  // 3D-vector of helixB point at path pB[r]
            testA.SetXYZ(helix_point[0],helix_point[1],helix_point[2]);
            distarray[1] = (testA-testB).Mag();
            flip = -1.0;
        }
        loopcounter++;
    }
    if(distarray[0] < distarray[1])
    {
        pathB = pB[0];
        pathA = pA[0];
        dcaAB = distarray[0];
    }
    else
    {
        pathB = pB[1];
        pathA = pA[1];
        dcaAB = distarray[1];
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Calculate_secondary_vertices(Int_t graphics)
{
    Double_t helix_pointA[3];
    Double_t helix_pointB[3];
    Double_t helix_pointAs[3];
    Double_t helix_pointBs[3];
    Double_t vertex_point[3];
    Double_t helix_pointA_est[3];
    Double_t helix_pointB_est[3];
    Double_t vertex_point_est[3];
    TVector3 TV3_dirA;
    TVector3 TV3_dirB;
    TVector3 TV3_dirAB;
    TLorentzVector TLV_A;
    TLorentzVector TLV_B;
    TLorentzVector TLV_AB;

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_sec_vertices        = new TEvePointSet();
        TEveP_first_point_helix   = new TEvePointSet();
        TEveP_second_point_helix  = new TEvePointSet();
    }
#endif

    Int_t i_vertex     = 0;
    Int_t i_comb       = 0;
    Int_t i_vertex_acc = 0;
    //printf("N_helices: %d \n",(Int_t)vec_helices.size());
    for(Int_t i_track_A = 0; i_track_A < ((Int_t)vec_helices.size() - 1); i_track_A++)
    {
        for(Int_t i_track_B = (i_track_A+1); i_track_B < (Int_t)vec_helices.size(); i_track_B++)
        {
            Float_t pathA_est, pathB_est, dcaAB_est;
            //printf("i_track_A: %d, i_track_B: %d, i_comb: %d \n",i_track_A,i_track_B,i_comb);
            //printf("3D cross point: {%4.3f, %4.3f, %4.3f} \n",vertex_point[0],vertex_point[1],vertex_point[2]);
            Int_t est_return = fDCA_Helix_Estimate(vec_helices[i_track_A],vec_helices[i_track_B],pathA_est,pathB_est,dcaAB_est);

            //printf(" \n");

            if(dcaAB_est > 15.0) continue;

            vec_helices[i_track_A] ->Evaluate(pathA_est,helix_pointA_est);
            vec_helices[i_track_B] ->Evaluate(pathB_est,helix_pointB_est);
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                vertex_point_est[i_xyz] = (helix_pointA_est[i_xyz] + helix_pointB_est[i_xyz])/2.0;
            }
            Double_t radius_est = TMath::Sqrt(vertex_point_est[0]*vertex_point_est[0] + vertex_point_est[1]*vertex_point_est[1]);
            if(radius_est > 400.0) continue;


            Float_t pathA, pathB, dcaAB;
            if(est_return == 0)
            {
                pathA_est = -9999.0;
                pathB_est = -9999.0;
            }
            fHelixABdca(vec_helices[i_track_A],vec_helices[i_track_B],pathA,pathB,dcaAB,pathA_est,pathB_est);
            //printf("track A,B: {%d, %d}, dcaAB: %4.3f \n",i_track_A,i_track_B,dcaAB);
            if(dcaAB < 10.0)
            {
                //if(i_comb == 38)
                {
                    vec_helices[i_track_A] ->Evaluate(pathA,helix_pointA);
                    vec_helices[i_track_B] ->Evaluate(pathB,helix_pointB);
                    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
                    {
                        vertex_point[i_xyz] = (helix_pointA[i_xyz] + helix_pointB[i_xyz])/2.0;
                    }

                    if(graphics)
                    {
                        //Draw_Kalman_Helix_Tracks(i_track_A,kGreen);
                        //Draw_Kalman_Helix_Tracks(i_track_B,kRed);
                    }

                    Double_t radius_vertex = TMath::Sqrt(vertex_point[0]*vertex_point[0] + vertex_point[1]*vertex_point[1]);
                    //if(radius_vertex > 240.0 && radius_vertex < 360.0)
#if defined(USEEVE)
                    if(graphics) TEveP_sec_vertices ->SetPoint(i_vertex,vertex_point[0],vertex_point[1],vertex_point[2]);
#endif

#if 1
                    //-------------------------------------------------
                    // Armenteros-Podolanski
                    if(radius_vertex < 250.0)
                    {
                        Double_t pTA = mHelices_kalman[i_track_A][6]; // pT
                        Double_t pzA = mHelices_kalman[i_track_A][7]; // pz
                        Double_t pA  = TMath::Sqrt(pTA*pTA + pzA*pzA); // p
                        Double_t pTB = mHelices_kalman[i_track_B][6]; // pT
                        Double_t pzB = mHelices_kalman[i_track_B][7]; // pz
                        Double_t pB  = TMath::Sqrt(pTB*pTB + pzB*pzB); // p
                        Double_t CA  = mHelices_kalman[i_track_A][4]; // curvature
                        Double_t CB  = mHelices_kalman[i_track_B][4]; // curvature

                        vec_helices[i_track_A] ->Evaluate(pathA+0.1,helix_pointAs);
                        vec_helices[i_track_B] ->Evaluate(pathB+0.1,helix_pointBs);
                        TV3_dirA.SetXYZ(helix_pointAs[0] - helix_pointA[0],helix_pointAs[1] - helix_pointA[1],helix_pointAs[2] - helix_pointA[2]);
                        TV3_dirB.SetXYZ(helix_pointBs[0] - helix_pointB[0],helix_pointBs[1] - helix_pointB[1],helix_pointBs[2] - helix_pointB[2]);

                        TV3_dirA *= pA/TV3_dirA.Mag();
                        TV3_dirB *= pB/TV3_dirB.Mag();
                        TV3_dirAB = TV3_dirA + TV3_dirB;
                        TV3_dirAB *= 1.0/TV3_dirAB.Mag();
                        Double_t projA = TV3_dirA.Dot(TV3_dirAB);
                        Double_t projB = TV3_dirB.Dot(TV3_dirAB);
                        Double_t AP_pT = TMath::Sqrt(TMath::Power(TV3_dirA.Mag(),2) - projA*projA);
                        Double_t AP_pTB = TMath::Sqrt(TMath::Power(TV3_dirB.Mag(),2) - projB*projB);
                        Double_t AP_alpha = (projA - projB)/(projA + projB);
                        if(CA < CB) AP_alpha *= -1.0;


                        //printf("i_comb: %d, A p,pT,pz: {%4.3f, %4.3f, %4.3f}, B p,pT,pz: {%4.3f, %4.3f, %4.3f}, AP pT, alpha: {%4.3f, %4.3f}, TV3_dirAB: {%4.3f, %4.3f, %4.3f} \n",i_comb,pA,pTA,pzA,pB,pTB,pzB,AP_pT,AP_alpha,TV3_dirAB[0],TV3_dirAB[1],TV3_dirAB[2]);

#if defined(USEEVE)
                        if(graphics)
                        {
                            TEveP_first_point_helix  ->SetPoint(i_vertex_acc,helix_pointA[0],helix_pointA[1],helix_pointA[2]);
                            TEveP_second_point_helix ->SetPoint(i_vertex_acc,helix_pointAs[0],helix_pointAs[1],helix_pointAs[2]);
                            TEveP_first_point_helix  ->SetPoint(i_vertex_acc+2,helix_pointB[0],helix_pointB[1],helix_pointB[2]);
                            TEveP_second_point_helix ->SetPoint(i_vertex_acc+2,helix_pointBs[0],helix_pointBs[1],helix_pointBs[2]);

                            TEveLine_mother.resize(i_vertex_acc+1);
                            TEveLine_mother[i_vertex_acc] = new TEveLine();

                            TEveLine_mother[i_vertex_acc] ->SetNextPoint(vertex_point[0],vertex_point[1],vertex_point[2]);
                            TEveLine_mother[i_vertex_acc] ->SetNextPoint(vertex_point[0] + 30.0*TV3_dirAB[0],vertex_point[1] + 30.0*TV3_dirAB[1],vertex_point[2] + 30.0*TV3_dirAB[2]);
                        }
#endif

                        if(pTA > 0.3 && pTB > 0.3 && CA*CB < 0.0 && (pzA + pzB) != 0.0)
                        {
                            //printf("CA: %4.3f, CB: %4.3f \n",CA,CB);


                            //printf("   ----> AP_alpha: %4.3f, AP_pT: %4.3f, AP_pTB: %4.3f, projA: %4.3f, projB: %4.3f \n",AP_alpha,AP_pT,AP_pTB,projA,projB);
                            TH2D_AP_plot ->Fill(AP_alpha,AP_pT);
                        }
                        i_vertex_acc++;
                    }
                    //-------------------------------------------------
#endif

                    i_vertex++;
                }
                i_comb++;
            }
        }
    }

    //printf("Number of secondary vertices found: %d \n",i_vertex);

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_sec_vertices  ->SetMarkerSize(3);
        TEveP_sec_vertices  ->SetMarkerStyle(20);
        TEveP_sec_vertices  ->SetMarkerColor(kRed);
        gEve->AddElement(TEveP_sec_vertices);

        TEveP_first_point_helix  ->SetMarkerSize(3);
        TEveP_first_point_helix  ->SetMarkerStyle(20);
        TEveP_first_point_helix  ->SetMarkerColor(kBlue);
        //gEve->AddElement(TEveP_first_point_helix);

        TEveP_second_point_helix  ->SetMarkerSize(3);
        TEveP_second_point_helix  ->SetMarkerStyle(20);
        TEveP_second_point_helix  ->SetMarkerColor(kMagenta);
        //gEve->AddElement(TEveP_second_point_helix);

        for(Int_t i_mother = 0; i_mother < (Int_t)TEveLine_mother.size(); i_mother++)
        {
            TEveLine_mother[i_mother]    ->SetLineStyle(1);
            TEveLine_mother[i_mother]    ->SetLineWidth(3);
            TEveLine_mother[i_mother]    ->SetMainColor(kMagenta);
            TEveLine_mother[i_mother]    ->SetMainAlpha(1.0);
            //gEve->AddElement(TEveLine_mother[i_mother]);
        }

        gEve->Redraw3D(kTRUE);
    }
#endif
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_Kalman_helix_params(vector<vector<Double_t>> mHelices_kalman_in)
{
    //printf("Ali_TRD_ST_Analyze::set_Kalman_helix_params");
    mHelices_kalman = mHelices_kalman_in;

    vec_helices.clear();
    vec_helices.resize((Int_t)mHelices_kalman.size());

    for(Int_t i_track = 0; i_track < (Int_t)mHelices_kalman.size(); i_track++)
    {
        vec_helices[i_track] = new Ali_Helix();
        vec_helices[i_track] ->setHelix(mHelices_kalman[i_track][0],mHelices_kalman[i_track][1],mHelices_kalman[i_track][2],mHelices_kalman[i_track][3],mHelices_kalman[i_track][4],mHelices_kalman[i_track][5]);
    }

}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_Kalman_TRD_tracklets(vector< vector<Ali_TRD_ST_Tracklets*> > vec_kalman_TRD_trackets_in)
{
    vec_kalman_TRD_trackets = vec_kalman_TRD_trackets_in;
    //printf("size: %d \n",(Int_t)vec_kalman_TRD_trackets.size());
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Match_kalman_tracks_to_TPC_tracks(Int_t graphics)
{
    //printf("size: %d \n",(Int_t)vec_kalman_TRD_trackets.size());
    // Loop over kalman matched TRD tracklets
    for(Int_t i_kalm_track = 0; i_kalm_track < (Int_t)vec_kalman_TRD_trackets.size(); i_kalm_track++)
    {
        vector<Int_t> vec_N_layers_matched;
        Int_t N_TPC_tracks = (Int_t)matched_tracks.size();
        vec_N_layers_matched.resize(N_TPC_tracks);
        Int_t N_layers_matched_all = 0;
        vector<Int_t> index_TPC_track_matched;
        for(Int_t i_TPC_track = 0; i_TPC_track < N_TPC_tracks; i_TPC_track++)
        {
            vec_N_layers_matched[i_TPC_track] = 0;
        }
        index_TPC_track_matched.clear();
        vec_N_layers_matched.clear();

        //printf("layers: %d \n",(Int_t)vec_kalman_TRD_trackets[i_kalm_track].size());
        Int_t N_good_kalman_tracklets = 0;
        for(Int_t i_layer = 0; i_layer < (Int_t)vec_kalman_TRD_trackets[i_kalm_track].size(); i_layer++)
        {
            if(vec_kalman_TRD_trackets[i_kalm_track][i_layer] == NULL) continue;
            N_good_kalman_tracklets++;
            TVector3 kalman_TV3_offset = vec_kalman_TRD_trackets[i_kalm_track][i_layer] ->get_TV3_offset();
            //printf("i_kalman_tracklet: %d, i_layer: %d, offset: {%4.3f, %4.3f, %4.3f} \n",i_kalm_track,i_layer,kalman_TV3_offset.X(),kalman_TV3_offset.Y(),kalman_TV3_offset.Z());

            // Loop over TPC tracks
            for(Int_t i_TPC_track = 0; i_TPC_track < N_TPC_tracks; i_TPC_track++)
            {
                Int_t flag_match = 0;
                for(Int_t i_layer_TPC = 0; i_layer_TPC < (Int_t)matched_tracks[i_TPC_track].size(); i_layer_TPC++)
                {
                    //matched_tracks
                    if(matched_tracks[i_TPC_track][i_layer_TPC] == NULL) continue;
                    TVector3 TPC_TV3_offset = matched_tracks[i_TPC_track][i_layer_TPC]->get_TV3_offset();
                    TVector3 TV3_diff =  kalman_TV3_offset - TPC_TV3_offset;
                    if(TV3_diff.Mag() < 0.1)
                    {
                        vec_N_layers_matched[i_TPC_track]++;
                        N_layers_matched_all++;
                        flag_match = 1;
                        //printf("i_kalm_track: %d, matched: %d \n",i_kalm_track,N_layers_matched_all);
                    }
                }
                if(vec_N_layers_matched[i_TPC_track] == 1 && flag_match) index_TPC_track_matched.push_back(i_TPC_track);
            }
        }

        if(N_layers_matched_all >= 2)
        {
            //printf("Kalman track %d with %d total matches \n",i_kalm_track,N_layers_matched_all);
            for(Int_t index_TPC_track = 0; index_TPC_track < (Int_t)index_TPC_track_matched.size(); index_TPC_track++)
            {
                Double_t matches_percent = 100.0*(Double_t)vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]]/(Double_t)N_good_kalman_tracklets;
                //printf("  --> TPC_track: %d, N_good_kalman_tracklets: %d N_matches: %d, %4.1f%% matches \n",index_TPC_track_matched[index_TPC_track],N_good_kalman_tracklets,vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]],matches_percent);

                if(matches_percent >= 90.0 && vec_N_layers_matched[index_TPC_track_matched[index_TPC_track]] >= 5)
                {
                    Int_t idx_matched_TPC_track = vec_idx_matched_TPC_track[index_TPC_track_matched[index_TPC_track]];

                    TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(idx_matched_TPC_track);

                    TLorentzVector TLV_part = TRD_ST_TPC_Track ->get_TLV_part();
                    Float_t momentum        = TLV_part.P();
                    Float_t pT_track        = TLV_part.Pt();
                    Double_t dca            = TRD_ST_TPC_Track ->getdca();  // charge * distance of closest approach to the primary vertex

                    Double_t pT_kalman   = mHelices_kalman[i_kalm_track][6]; // pT
                    Double_t pz_kalman   = mHelices_kalman[i_kalm_track][7]; // pz
                    Double_t curv_kalman = mHelices_kalman[i_kalm_track][4]; // curvature

                    if(graphics)
                    {
                        Draw_TPC_track(idx_matched_TPC_track,kAzure-2,3);
                        Draw_Kalman_Helix_Tracks(i_kalm_track,kGreen+1);
                    }
                    //printf("      >>>>> %s idx (Kalman): %d %s, %s idx (TPC): %d %s, pT (TPC): %4.3f, pT (kalman): %4.3f \n",KRED,i_kalm_track,KNRM,KBLU,idx_matched_TPC_track,KNRM,pT_track,pT_kalman);

                    TH2D_pT_TPC_vs_Kalman ->Fill(-TMath::Sign(1,curv_kalman)*pT_kalman,TMath::Sign(1,dca)*pT_track);
                }
            }
        }

    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::set_single_helix_params(vector<Double_t> vec_params)
{
    for(Int_t i_param = 0; i_param < 6; i_param++)
    {
        aliHelix_params[i_param] = vec_params[i_param];
    }
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Evaluate(Double_t t,Double_t r[3])  //radius vector
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives at given phase
  //--------------------------------------------------------------------
  float phase=aliHelix_params[4]*t+aliHelix_params[2];
  Double_t sn=sinf(phase), cs=cosf(phase);
  //  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = aliHelix_params[5] + sn/aliHelix_params[4];
  r[1] = aliHelix_params[0] - cs/aliHelix_params[4];
  r[2] = aliHelix_params[1] + aliHelix_params[3]*t;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_Kalman_Helix_Tracks(Int_t n_track, Int_t color)
{
    // n_track = -1 -> all tracks drawn
    Int_t i_track_start = 0;
    Int_t i_track_stop  = (Int_t)mHelices_kalman.size();
    if(n_track != -1)
    {
        i_track_start = n_track;
        i_track_stop  = n_track + 1;
    }
    for(Int_t i_track = i_track_start; i_track < i_track_stop; i_track++)
    {
        set_single_helix_params(mHelices_kalman[i_track]);
        Double_t pT_kalman   = mHelices_kalman[i_track][6]; // pT
        //if(pT_kalman < 0.5) continue;
        Double_t track_pos[3];

#if defined(USEEVE)
        vec_TPL3D_helix_kalman.resize(i_track+1);
        vec_TPL3D_helix_kalman_hull.resize(i_track+1);
        vec_TPL3D_helix_kalman_inner.resize(i_track+1);
        vec_TPL3D_helix_kalman[i_track] = new TEveLine();
        vec_TPL3D_helix_kalman_hull[i_track] = new TEveLine();
        vec_TPL3D_helix_kalman_inner[i_track] = new TEveLine();
#endif

        Float_t pathA_dca, dcaAB_dca;
        fHelixAtoPointdca(TV3_EventVertex,vec_helices[i_track],pathA_dca,dcaAB_dca); // new helix to point dca calculation
        //if(dcaAB_dca > 50.0) continue;

        Double_t radius_helix = 0.0;
        for(Double_t track_path = 75.0; track_path > -350; track_path -= 1.0)
        {
            Evaluate(track_path,track_pos);
            radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
            //printf("i_track: %d, track_path: %4.3f, pos: {%4.3f, %4.3f, %4.3f}, radius_helix: %4.3f \n",i_track,track_path,track_pos[0],track_pos[1],track_pos[2],radius_helix);
            if(isnan(radius_helix)) break;
            if(radius_helix > 500.0) break;
            if(fabs(track_pos[2]) > 430.0) break;
            //if(radius_helix > 80.0)
            {
#if defined(USEEVE)
                vec_TPL3D_helix_kalman[i_track]        ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
                vec_TPL3D_helix_kalman_hull[i_track]   ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
#endif
            }
        }

#if defined(USEEVE)
        HistName = "track kalman ";
        HistName += i_track;
        vec_TPL3D_helix_kalman[i_track]    ->SetName(HistName.Data());
        vec_TPL3D_helix_kalman[i_track]    ->SetLineStyle(1);
        vec_TPL3D_helix_kalman[i_track]    ->SetLineWidth(3);
        vec_TPL3D_helix_kalman[i_track]    ->SetMainColor(color);
        vec_TPL3D_helix_kalman[i_track]    ->SetMainAlpha(1.0);

        HistName = "track kalman (h) ";
        HistName += i_track;
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetName(HistName.Data());
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetLineStyle(1);
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetLineWidth(8);
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetMainColor(kOrange);
        vec_TPL3D_helix_kalman_hull[i_track]    ->SetMainAlpha(0.3);

        gEve->AddElement(vec_TPL3D_helix_kalman[i_track]);
        gEve->AddElement(vec_TPL3D_helix_kalman_hull[i_track]);
#endif
    }

#if defined(USEEVE)
    gEve->Redraw3D(kTRUE);
#endif
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Init_tree(TString SEList)
{
    printf("Ali_TRD_ST_Analyze::Init_tree \n");
    TString pinputdir = input_dir;

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
Int_t Ali_TRD_ST_Analyze::Loop_event(Long64_t i_event, Int_t graphics)
{
    //printf("Ali_TRD_ST_Analyze::Loop_event \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event


    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
    TV3_EventVertex.SetXYZ(EventVertexX,EventVertexY,EventVertexZ);
    Float_t  V0MEq                = TRD_ST_Event ->getcent_class_V0MEq();

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_primary_vertex ->SetPoint(0,EventVertexX,EventVertexY,EventVertexZ);
        TEveP_primary_vertex ->SetMarkerStyle(20);
        TEveP_primary_vertex ->SetMarkerSize(4);
        TEveP_primary_vertex ->SetMarkerColor(kYellow+1);
        gEve ->AddElement(TEveP_primary_vertex);
    }
#endif
    //--------------------------------------------------

    //printf("Event: %lld, TPC tracks: %d, TRD tracklets: %d \n",i_event,NumTracks,NumTracklets);

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

    Tracklets = new Ali_TRD_ST_Tracklets*[NumTracklets];
    Number_Tracklets=NumTracklets;
	
    for(Int_t i_tracklet = 0; i_tracklet < NumTracklets; i_tracklet++)
    {
        TRD_ST_Tracklet = TRD_ST_Event    ->getTracklet(i_tracklet);
        TRD_ST_Tracklet->set_TRD_index(i_tracklet);
        Tracklets[i_tracklet]=TRD_ST_Tracklet;
        TV3_offset      = TRD_ST_Tracklet ->get_TV3_offset();
        TV3_dir         = TRD_ST_Tracklet ->get_TV3_dir();
        i_det           = TRD_ST_Tracklet ->get_TRD_det();
		
        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        Double_t radius = TMath::Sqrt( TMath::Power(TV3_offset[0],2) + TMath::Power(TV3_offset[1],2) );

        th1d_TRD_layer_radii ->Fill(radius);
        vec_th1d_TRD_layer_radii_det[i_det] ->Fill(radius);
    }
    //--------------------------------------------------


    //TCanvas* can_TRD_layer_radii = new TCanvas("can_TRD_layer_radii","can_TRD_layer_radii",10,10,500,500);
    //can_TRD_layer_radii ->cd();
    //th1d_TRD_layer_radii ->DrawCopy("h");

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_TPC_track(Int_t i_track, Int_t color, Double_t line_width)
{
    Double_t track_pos[3];
    Double_t radius_helix;

    TRD_ST_TPC_Track = TRD_ST_Event ->getTrack(i_track);

#if defined(USEEVE)
    vec_TPL3D_helix.resize(i_track+1);
    vec_TPL3D_helix_hull.resize(i_track+1);
    vec_TPL3D_helix_inner.resize(i_track+1);
    vec_TPL3D_helix[i_track] = new TEveLine();
    vec_TPL3D_helix_hull[i_track] = new TEveLine();
    vec_TPL3D_helix_inner[i_track] = new TEveLine();
#endif
    //for(Int_t i_param=0;i_param<6;i_param++)
    //    cout<<TRD_ST_TPC_Track ->getHelix_param(i_param)<<" ";
    //cout<<endl;
    for(Double_t track_path = 0.0; track_path < 1000; track_path += 1.0)
    {
        TRD_ST_TPC_Track ->Evaluate(track_path,track_pos);
        radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
        if(radius_helix > 370.0) break;
        if(fabs(track_pos[2]) > 360.0) break;
        if(radius_helix > 80.0)
        {
#if defined(USEEVE)
            vec_TPL3D_helix[i_track]        ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
            vec_TPL3D_helix_hull[i_track]   ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
#endif
        }
        if(radius_helix < 80.0)
        {
#if defined(USEEVE)
            vec_TPL3D_helix_inner[i_track] ->SetNextPoint(track_pos[0],track_pos[1],track_pos[2]);
#endif
        }


        //if(i_track == 0) printf("track_path: %4.3f, pos: {%4.2f, %4.2f, %4.2f} \n",track_path,track_pos[0],track_pos[1],track_pos[2]);
    }

#if defined(USEEVE)
    HistName = "track ";
    HistName += i_track;
    vec_TPL3D_helix[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix[i_track]    ->SetLineWidth(line_width);
    vec_TPL3D_helix[i_track]    ->SetMainColor(color);
    vec_TPL3D_helix[i_track]    ->SetMainAlpha(1.0);

    HistName = "track (h) ";
    HistName += i_track;
    vec_TPL3D_helix_hull[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix_hull[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix_hull[i_track]    ->SetLineWidth(8);
    vec_TPL3D_helix_hull[i_track]    ->SetMainColor(kWhite);
    vec_TPL3D_helix_hull[i_track]    ->SetMainAlpha(0.3);

    HistName = "track (h) ";
    HistName += i_track;
    vec_TPL3D_helix_inner[i_track]    ->SetName(HistName.Data());
    vec_TPL3D_helix_inner[i_track]    ->SetLineStyle(1);
    vec_TPL3D_helix_inner[i_track]    ->SetLineWidth(2);
    vec_TPL3D_helix_inner[i_track]    ->SetMainColor(kGray);
    vec_TPL3D_helix_inner[i_track]    ->SetMainAlpha(0.8);
    //if(i_track == 3)
    {
        gEve->AddElement(vec_TPL3D_helix[i_track]);
        gEve->AddElement(vec_TPL3D_helix_hull[i_track]);
        gEve->AddElement(vec_TPL3D_helix_inner[i_track]);
    }
#endif

    //printf("%s Ali_TRD_ST_Analyze::Draw_TPC_track, i_track: %d %s \n",KRED,i_track,KNRM);
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::Draw_event(Long64_t i_event)
{
    //printf("Ali_TRD_ST_Analyze::Draw_event \n");
	
    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

	
    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
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

        if(momentum < 0.3 || dca > 2.0) continue;


        Draw_TPC_track(i_track,track_color,3);
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

        if(TV3_offset.Mag() > 1000.0) continue;

        //create "tracklets"
        Int_t i_sector = (Int_t)(i_det/30);
        Int_t i_stack  = (Int_t)(i_det%30/6);
        Int_t i_layer  = i_det%6;
        //Int_t i_det = layer + 6*stack + 30*sector;

        Double_t radius = TMath::Sqrt( TMath::Power(TV3_offset[0],2) + TMath::Power(TV3_offset[1],2) );
      

#if defined(USEEVE)
        vec_TEveLine_tracklets[i_layer].resize(N_tracklets_layers[i_layer]+1);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] = new TEveLine();

        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0],TV3_offset[1],TV3_offset[2]);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]] ->SetNextPoint(TV3_offset[0] + scale_length_vec*TV3_dir[0],TV3_offset[1] + scale_length_vec*TV3_dir[1],TV3_offset[2] + scale_length_vec*TV3_dir[2]);


        //printf("i_tracklet: %d, radius: %4.3f, pos A: {%4.2f, %4.2f, %4.2f}, pos B: {%4.2f, %4.2f, %4.2f} \n",i_tracklet,radius,TV3_offset[0],TV3_offset[1],TV3_offset[2],TV3_offset[0] + scale_length_vec*TV3_dir[0],TV3_offset[1] + scale_length_vec*TV3_dir[1],TV3_offset[2] + scale_length_vec*TV3_dir[2]);

        HistName = "tracklet ";
        HistName += i_tracklet;
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetName(HistName.Data());
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetLineStyle(1);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetLineWidth(3);
        vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]    ->SetMainColor(color_layer[i_layer]);
        //if(i_tracklet == 63 || i_tracklet == 67 || i_tracklet == 72 || i_tracklet == 75 || i_tracklet == 83 || i_tracklet == 88)
        {
             gEve->AddElement(vec_TEveLine_tracklets[i_layer][N_tracklets_layers[i_layer]]);
        }
#endif

        N_tracklets_layers[i_layer]++;
    }

#if defined(USEEVE)
    gEve->Redraw3D(kTRUE);
#endif
    //--------------------------------------------------

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
Int_t Ali_TRD_ST_Analyze::Do_TPC_TRD_matching(Long64_t i_event, Double_t xy_matching_window, Double_t z_matching_window, Int_t graphics)
{
    //printf("Ali_TRD_ST_Analyze::Do_TPC_TRD_matching \n");

    if (!input_SE->GetEntry( i_event )) return 0; // take the event -> information is stored in event

    Int_t offset_point = 0;
    Int_t TPC_at_offset_point = 0;
    matched_tracks.clear();
    vec_idx_matched_TPC_track.clear();

    //--------------------------------------------------
    // Event information (more data members available, see Ali_TRD_ST_Event class definition)
    UShort_t NumTracks            = TRD_ST_Event ->getNumTracks(); // number of tracks in this event
    Int_t    NumTracklets         = TRD_ST_Event ->getNumTracklets();
    EventVertexX         = TRD_ST_Event ->getx();
    EventVertexY         = TRD_ST_Event ->gety();
    EventVertexZ         = TRD_ST_Event ->getz();
    Float_t  V0MEq                = TRD_ST_Event ->getcent_class_V0MEq();
    //--------------------------------------------------



    //--------------------------------------------------
    // Store all TRD tracklets in an array to keep it in memory
    // TRD tracklet loop
    vector< vector<TVector3> > vec_TV3_dir_tracklets;
    vector< vector<TVector3> > vec_TV3_offset_tracklets;
    vector< vector<Int_t> > TRD_index_tracklets;
    vec_TV3_dir_tracklets.resize(540);
    vec_TV3_offset_tracklets.resize(540);
    TRD_index_tracklets.resize(540);

    Double_t ADC_val[24];

    TVector3 TV3_offset;
    TVector3 TV3_dir;
    Int_t    i_det;

#if defined(USEEVE)
    if(graphics)
    {
        TEveP_offset_points        = new TEvePointSet();
        TEveP_TPC_at_offset_points = new TEvePointSet();
    }
#endif

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

        vec_TV3_dir_tracklets[i_det].push_back(TV3_dir);
        vec_TV3_offset_tracklets[i_det].push_back(TV3_offset);
        TRD_index_tracklets[i_det].push_back(i_tracklet);
    }
    //--------------------------------------------------

    //--------------------------------------------------
    // TPC track loop
    Double_t track_pos[3];
    Double_t radius_helix;
    for(Int_t i_track = 0; i_track < NumTracks; i_track++)
        //for(Int_t i_track = 3; i_track < 4; i_track++)
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
        Float_t charge = 1.0;
        if(dca < 0.0) charge = -1.0;

        Float_t momentum        = TLV_part.P();
        Float_t eta_track       = TLV_part.Eta();
        Float_t pT_track        = TLV_part.Pt();
        Float_t theta_track     = TLV_part.Theta();
        Float_t phi_track       = TLV_part.Phi();

        if(momentum < 0.3) continue;

        vector<TVector3> vec_TV3_helix_points_at_TRD_layers;
        vec_TV3_helix_points_at_TRD_layers.resize(6);
        matched_tracks.resize(matched_tracks.size() + 1);
        Int_t tracks_size = matched_tracks.size() - 1;  // index of current track
        vec_idx_matched_TPC_track.push_back(i_track);

        //printf("tracks_size: %d, size_TPC_idx: %d \n",tracks_size,(Int_t)vec_idx_matched_TPC_track.size());

        for(Int_t i_layer = 0; i_layer < 6; i_layer++)
        {
            Double_t radius_layer_center = 0.5*(TRD_layer_radii[i_layer][0] + TRD_layer_radii[i_layer][1]);

            // Find the helix path which touches the first TRD layer
            Double_t track_path_add      = 10.0;
            Double_t track_path_layer0   = 0.0;
            Double_t radius_helix_layer0 = 0.0;
            for(Double_t track_path = TRD_layer_radii[i_layer][0]; track_path < 1000; track_path += track_path_add)
            {
                TRD_ST_TPC_Track ->Evaluate(track_path,track_pos);
                radius_helix = TMath::Sqrt( TMath::Power(track_pos[0],2) + TMath::Power(track_pos[1],2) );
                if(radius_helix < radius_layer_center)
                {
                    track_path_layer0   = track_path;
                    radius_helix_layer0 = radius_helix;
                    //printf("radius_helix_layer0: %4.3f, track_path_add: %4.3f \n",radius_helix_layer0,track_path_add);
                }
                else
                {
                    track_path -= track_path_add;
                    track_path_add *= 0.5;
                }
                if(track_path_add < 1.0)
                {
                    vec_TV3_helix_points_at_TRD_layers[i_layer].SetXYZ(track_pos[0],track_pos[1],track_pos[2]);
                    break;
                }
            }
            //printf("   --> i_layer: %d, track_path_layer0: %4.3f, radius_helix_layer0: %4.3f \n",i_layer,track_path_layer0,radius_helix_layer0);


            TVector3 TV3_diff_vec;
            Double_t dist_min      = 10.0;
            Double_t dist          = 0.0;
            Int_t    det_best      = -1;
            Int_t    tracklet_best = -1;

            for(Int_t i_det = 0; i_det < 540; i_det++)
            {
                Int_t i_layer_from_det  = i_det%6;
                if(i_layer_from_det != i_layer) continue;

                for(Int_t i_tracklet = 0; i_tracklet < (Int_t)vec_TV3_offset_tracklets[i_det].size(); i_tracklet++)
                {
                    TV3_diff_vec = vec_TV3_offset_tracklets[i_det][i_tracklet] - vec_TV3_helix_points_at_TRD_layers[i_layer];
                    dist = TV3_diff_vec.Mag();
                    if(dist < dist_min)
                    {
                        dist_min      = dist;
                        det_best      = i_det;
                        tracklet_best = i_tracklet;
                    }
                }
            }

            matched_tracks[tracks_size].push_back(NULL);

            if(det_best < 0 || tracklet_best < 0) continue;

            Ali_TRD_ST_Tracklets* temp_tracklet = new Ali_TRD_ST_Tracklets();
            temp_tracklet->set_TRD_det(det_best);
            temp_tracklet->set_TV3_offset(vec_TV3_offset_tracklets[det_best][tracklet_best]);
            temp_tracklet->set_TV3_dir(vec_TV3_dir_tracklets[det_best][tracklet_best]);
            temp_tracklet->set_TRD_index(TRD_index_tracklets[det_best][tracklet_best]);
            matched_tracks[tracks_size][matched_tracks[tracks_size].size() - 1] = temp_tracklet;


            //------------------------------------------
            // Calculate exact position of TPC helix at TRD tracklet offset position
            TPC_single_helix ->setHelix(TRD_ST_TPC_Track->getHelix_param(0),TRD_ST_TPC_Track->getHelix_param(1),TRD_ST_TPC_Track->getHelix_param(2),TRD_ST_TPC_Track->getHelix_param(3),TRD_ST_TPC_Track->getHelix_param(4),TRD_ST_TPC_Track->getHelix_param(5));
            Float_t pathA_dca, dcaAB_dca;
            fHelixAtoPointdca(vec_TV3_offset_tracklets[det_best][tracklet_best],TPC_single_helix,pathA_dca,dcaAB_dca); // new helix to point dca calculation
            Double_t track_pos_at_tracklet[3];
            TRD_ST_TPC_Track ->Evaluate(pathA_dca,track_pos_at_tracklet);
            //printf("%s i_layer: %d %s, i_track: %d, dcaAB_dca: %4.3f, offset_point: {%4.3f, %4.3f, %4.3f}, TPC point: {%4.3f, %4.3f, %4.3f} \n",KMAG,i_layer,KNRM,i_track,dcaAB_dca,vec_TV3_offset_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2],track_pos_at_tracklet[0],track_pos_at_tracklet[1],track_pos_at_tracklet[2]);
#if defined(USEEVE)
            if(graphics) TEveP_TPC_at_offset_points  ->SetPoint(TPC_at_offset_point,track_pos_at_tracklet[0],track_pos_at_tracklet[1],track_pos_at_tracklet[2]);
#endif
            TPC_at_offset_point++;

            Double_t residuals_xyz;
            for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
            {
                residuals_xyz = track_pos_at_tracklet[i_xyz] - vec_TV3_offset_tracklets[det_best][tracklet_best][i_xyz];
                vec_h2D_pT_vs_TPC_TRD_residuals[det_best][i_xyz] ->Fill(residuals_xyz,charge*pT_track);
            }
            //------------------------------------------

#if defined(USEEVE)
            Int_t size_tracklet = (Int_t)vec_TEveLine_tracklets_match[i_layer].size();
#endif

            //printf("i_layer: %d, size_tracklet: %d \n",i_layer,size_tracklet);

#if defined(USEEVE)
            if(graphics)
            {
                vec_TEveLine_tracklets_match[i_layer].resize(size_tracklet+1);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet] = new TEveLine();
                vec_TEveLine_tracklets_match[i_layer][size_tracklet] ->SetNextPoint(vec_TV3_offset_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2]);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet] ->SetNextPoint(vec_TV3_offset_tracklets[det_best][tracklet_best][0] + scale_length_vec*vec_TV3_dir_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1] + scale_length_vec*vec_TV3_dir_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2] + scale_length_vec*vec_TV3_dir_tracklets[det_best][tracklet_best][2]);

                HistName = "tracklet (m) ";
                HistName += size_tracklet;
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetName(HistName.Data());
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetLineStyle(1);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetLineWidth(6);
                vec_TEveLine_tracklets_match[i_layer][size_tracklet]    ->SetMainColor(color_layer_match[i_layer]);

                TEveP_offset_points  ->SetPoint(offset_point,vec_TV3_offset_tracklets[det_best][tracklet_best][0],vec_TV3_offset_tracklets[det_best][tracklet_best][1],vec_TV3_offset_tracklets[det_best][tracklet_best][2]);
                TEveP_offset_points  ->SetMarkerSize(2);
                TEveP_offset_points  ->SetMarkerStyle(20);
                TEveP_offset_points  ->SetMarkerColor(kCyan+1);
                offset_point++;
                //if(i_tracklet == 63 || i_tracklet == 67 || i_tracklet == 72 || i_tracklet == 75 || i_tracklet == 83 || i_tracklet == 88)
                {
                    gEve->AddElement(vec_TEveLine_tracklets_match[i_layer][size_tracklet]);
                }
            }
#endif
        }


        //--------------------------------------------------
        Int_t N_matched_tracklets = 0;
        for(Int_t i_layer = 0; i_layer < (Int_t)matched_tracks[tracks_size].size(); i_layer++)
        {
            if(matched_tracks[tracks_size][i_layer] != NULL) N_matched_tracklets++;
        }
        //printf("%s TPC track: %d %s, N_matched_tracklets: %d \n",KYEL,vec_idx_matched_TPC_track[tracks_size],KNRM,N_matched_tracklets);
        //if(N_matched_tracklets >= 4 && vec_idx_matched_TPC_track[tracks_size] == 25) Draw_TPC_track(vec_idx_matched_TPC_track[tracks_size],kAzure-2,2);
        //printf(" \n");
        //--------------------------------------------------


    } // end of track loop
    //--------------------------------------------------

#if defined(USEEVE)
    if(graphics)
    {
        gEve->AddElement(TEveP_offset_points);

        TEveP_TPC_at_offset_points  ->SetMarkerSize(3);
        TEveP_TPC_at_offset_points  ->SetMarkerStyle(20);
        TEveP_TPC_at_offset_points  ->SetMarkerColor(kYellow);
        //gEve->AddElement(TEveP_TPC_at_offset_points);
        gEve->Redraw3D(kTRUE);
    }
#endif
    //printf("Ali_TRD_ST_Analyze::Do_TPC_TRD_matching, N TPC Tracks: %d, matched_tracks: %d \n",NumTracks,(Int_t)matched_tracks.size());

    return 1;
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_Kalman_Tracks(vector< vector<Ali_TRD_ST_Tracklets*> > found_tracks)
{
#if defined(USEEVE)
    vector< vector<TEveLine*> > TEveLine_Kalman_found;
    TEveLine_Kalman_found.resize((Int_t)found_tracks.size());
#endif
    for(Int_t i_Track=0;i_Track<(Int_t)found_tracks.size();i_Track++)
    {
#if defined(USEEVE)
        TEveLine_Kalman_found[i_Track].resize(6);
#endif

        for(Int_t i_lay=0;i_lay<6;i_lay++)
        {
            if(found_tracks[i_Track][i_lay]!=NULL)
            {
#if defined(USEEVE)
                TEveLine_Kalman_found[i_Track][i_lay]=new TEveLine();
                TEveLine_Kalman_found[i_Track][i_lay]->SetNextPoint(found_tracks[i_Track][i_lay]->get_TV3_offset()[0],found_tracks[i_Track][i_lay]->get_TV3_offset()[1],found_tracks[i_Track][i_lay]->get_TV3_offset()[2]);
                TEveLine_Kalman_found[i_Track][i_lay]->SetNextPoint(found_tracks[i_Track][i_lay]->get_TV3_offset()[0] + scale_length_vec*found_tracks[i_Track][i_lay]->get_TV3_dir()[0],found_tracks[i_Track][i_lay]->get_TV3_offset()[1] + scale_length_vec*found_tracks[i_Track][i_lay]->get_TV3_dir()[1],found_tracks[i_Track][i_lay]->get_TV3_offset()[2] + scale_length_vec*found_tracks[i_Track][i_lay]->get_TV3_dir()[2]);
#endif

                Int_t i_det           = found_tracks[i_Track][i_lay] ->get_TRD_det();
                Int_t i_sector = (Int_t)(i_det/30);
                Int_t i_stack  = (Int_t)(i_det%30/6);
                Int_t i_layer  = i_det%6;


#if defined(USEEVE)
                HistName = "tracklet (kal) ";
                HistName += i_Track;
                HistName += "_";
                HistName += i_sector;
                HistName += "_";
                HistName += i_stack;
                HistName += "_";
                HistName += i_lay;

                TEveLine_Kalman_found[i_Track][i_lay]	->SetName(HistName.Data());
                TEveLine_Kalman_found[i_Track][i_lay]   ->SetLineStyle(1);
                TEveLine_Kalman_found[i_Track][i_lay]  	->SetLineWidth(6);
                TEveLine_Kalman_found[i_Track][i_lay]   ->SetMainColor(kOrange);
                gEve->AddElement(TEveLine_Kalman_found[i_Track][i_lay]);
#endif

            }
        }
    }
#if defined(USEEVE)
    gEve->Redraw3D(kTRUE);
#endif
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Plot_AP()
{
    TCanvas *can_AP_plot = new TCanvas("can_AP_plot", "can_AP_plot",10,10,500,500);
    can_AP_plot->cd();
    TH2D_AP_plot ->DrawCopy("colz");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Plot_pT_TPC_vs_Kalman()
{
    TCanvas *can_pT_TPC_vs_Kalman = new TCanvas("can_pT_TPC_vs_Kalman", "can_pT_TPC_vs_Kalman",10,10,500,500);
    can_pT_TPC_vs_Kalman->cd();
    TH2D_pT_TPC_vs_Kalman ->DrawCopy("colz");
}
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Draw_hist_TPC_tracklet_diffs()
{
    printf("Ali_TRD_ST_Analyze::Draw_hist_TPC_tracklet_diffs \n");

    th1d_offset_diff->GetXaxis()->SetTitle("Offset Difference (cm)");
    th1d_offset_diff->GetYaxis()->SetTitle("Number of Tracklets");
    th1d_offset_diff->GetXaxis()->CenterTitle();
    th1d_offset_diff->GetYaxis()->CenterTitle();

    TCanvas *th1d_offset_diff_can = new TCanvas("th1d_offset_diff_can", "Max");
    th1d_offset_diff_can->cd();
    th1d_offset_diff->Draw();

    th1d_angle_diff->GetXaxis()->SetTitle("Angle Difference (deg)");
    th1d_angle_diff->GetYaxis()->SetTitle("Number of Tracklets");
    th1d_angle_diff->GetXaxis()->CenterTitle();
    th1d_angle_diff->GetYaxis()->CenterTitle();

    TCanvas *th1d_angle_diff_can = new TCanvas("th1d_angle_diff_can", "Plateu");
    th1d_angle_diff_can->cd();
    th1d_angle_diff->Draw();
}
//----------------------------------------------------------------------------------------

TH1I* Ali_TRD_ST_Analyze::get_h_good_bad_TRD_chambers(){
	return h_good_bad_TRD_chambers;
}


//----------------------------------------------------------------------------------------
void Ali_TRD_ST_Analyze::Write()
{
    printf("Write data to file \n");
    outputfile ->cd();
    TH2D_AP_plot          ->Write();
    TH2D_pT_TPC_vs_Kalman ->Write();
    outputfile ->mkdir("layer_radii");
    outputfile ->cd("layer_radii");
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        vec_th1d_TRD_layer_radii_det[TRD_detector] ->Write();
    }
    outputfile ->cd();

    outputfile ->mkdir("TPC_TRD_residuals");
    outputfile ->cd("TPC_TRD_residuals");
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
        {
            vec_h2D_pT_vs_TPC_TRD_residuals[i_det][i_xyz] ->Write();
        }
    }
    outputfile ->cd();

    printf("All data written \n");
}
//----------------------------------------------------------------------------------------


