#include "TRD_Kalman_Tracking.h"
#include <chrono>

double T_calc;                          //Variable for the calculation time
double T_wall;

// Colors for printf statements
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

vector<vector<Double_t>> TRD_Kalman_Trackfinder::get_Kalman_helix_params()
{
    return mHelices;
}

//gives the 4 measurements of a tracklet (y,z,sin(phi),tan(lambda)) back
ROOT::Math::SVector<double,4> TRD_Kalman_Trackfinder::measure(Ali_TRD_ST_Tracklets* tracklet){
	ROOT::Math::SVector<double,4> measurement;
	TVector3 offset=tracklet->get_TV3_offset();
   	TVector3 dir=tracklet->get_TV3_dir();
	Double_t hypo=1./(TMath::Sqrt(dir[1]*dir[1] +dir[0]*dir[0]));
	measurement[0]=offset[1];
   	measurement[1]=offset[2];
	measurement[2]=dir[1]*hypo;
	measurement[3]=dir[2]*hypo;
	return measurement;
}

//checks if 2 tracklets are more or less in a line
Bool_t TRD_Kalman_Trackfinder::fitting(Ali_TRD_ST_Tracklets* a,Ali_TRD_ST_Tracklets* b){
	//direction is not ok 
	//if ((a->get_TRD_index()==1400)&&(b->get_TRD_index()==389)) cout<<"a2";
	
	
	TVector3 dir_a		=	a->get_TV3_dir();
	TVector3 dir_b		=	b->get_TV3_dir();
	TVector3 a_angle 	= 	(TVector3){dir_a[0], dir_a[1], 0};
	TVector3 b_angle 	= 	(TVector3){dir_b[0], dir_b[1], 0};
	/*
	if ((a->get_TRD_index()==1400)&&(b->get_TRD_index()==389)){ 
        */
        /*
        cout<<"a ";
		dir_a.Print();
		cout<<" b ";
		dir_b.Print();
                cout<<abs(a_angle.Angle(b_angle)*TMath::RadToDeg())<<" "<<abs(dir_a.Angle(dir_b)*TMath::RadToDeg() )<<endl;
                */
	/*}
	*/
	if (abs(a_angle.Angle(b_angle)*TMath::RadToDeg() )>15)
		return 0;
	if (abs(dir_a.Angle(dir_b)*TMath::RadToDeg() )>15)
		return 0;	
	
		
	//position is not ok
	TVector3 off_b 		=	b->get_TV3_offset();
	TVector3 off_a 		=	a->get_TV3_offset();
	TVector3 z			=	off_a+dir_a*((off_b[0]- off_a[0])/dir_a[0] );
	/*
	if ((a->get_TRD_index()==1400)&&(b->get_TRD_index()==389)){ 
		cout<<"a ";
		z.Print();
		cout<<" b ";
		off_b.Print();
	}*/
	
	if (abs((z -off_b)[1])>7.)
		return 0;
	//if ((a->get_TRD_index()==1648) && (b->get_TRD_index()==1643)) cout<<"a"<<endl;
	
	if (abs((z - off_b)[2])>20)
		return 0;
	return 1;
}

void TRD_Kalman_Trackfinder::get_seed( Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets){
		
		
	Int_t nbr_sectors=18;
	Int_t nbr_stack=5;
	Int_t nbr_layers=6; 	
	mBins.clear();
	mBins.resize(nbr_sectors*nbr_stack*nbr_layers);//bins for every detector
	mVisited.clear();
	mVisited.resize(nbr_sectors*nbr_stack*nbr_layers);
	
	//sort then in bins
	for(Int_t i_tracklet=0;i_tracklet<Num_Tracklets;i_tracklet++){
		if(Tracklets[i_tracklet]->get_TV3_offset().Mag() > 1000.0) continue;

			
		Int_t i_det = Tracklets[i_tracklet] ->get_TRD_det();
		Int_t i_sector = (Int_t)(i_det/30);
	
		//rotate the tracklets to the local coordiante system
		TVector3 temp1=Tracklets[i_tracklet]->get_TV3_offset();
		temp1.RotateZ((Double_t)(-1*(2*i_sector+1)*TMath::Pi()/18));
		Tracklets[i_tracklet]->set_TV3_offset(temp1);
			
		TVector3 temp2=Tracklets[i_tracklet]->get_TV3_dir();
		temp2.RotateZ((Double_t)(-1*(2*i_sector+1)*TMath::Pi()/18));
		Tracklets[i_tracklet]->set_TV3_dir(temp2);
		Tracklets[i_tracklet]->set_TRD_index(i_tracklet);
		
		//sort them in
		mBins[i_det].push_back(Tracklets[i_tracklet]);			
		mVisited[i_det].push_back(0);
	}		
	
	//check for each tracklet in the 6. and 5. layer whether there are matching tracklets in the 4. or 5. layer
	Bool_t done=0;
	for(int i_sec=0;i_sec<nbr_sectors;i_sec++) 
			for(int i_stck=0;i_stck<nbr_stack;i_stck++)
				for(int i_lay=nbr_layers-1;i_lay>1;i_lay--){
					Int_t i_det=i_sec*nbr_stack*nbr_layers+i_stck*nbr_layers+i_lay;
					for(int i_tracklet=0;i_tracklet<(Int_t)mBins[i_det].size();i_tracklet++)
						if (!(mVisited[i_det][i_tracklet])){
							mVisited[i_det][i_tracklet]=1;
							for(int j_lay=i_lay-1; j_lay>i_lay-3;j_lay--){
								done=0;
								Int_t j_det=i_det+j_lay -i_lay;
								for(int j_tracklet=0;j_tracklet<(Int_t)mBins[j_det].size();j_tracklet++)
									if(fitting(mBins[i_det][i_tracklet],mBins[j_det][j_tracklet])){
										mVisited[j_det][j_tracklet]=1;
										vector<Ali_TRD_ST_Tracklets*> temp_seed;
										temp_seed.resize(2);
										temp_seed[0]=mBins[i_det][i_tracklet];
										temp_seed[1]=mBins[j_det][j_tracklet];
										mSeed.push_back(temp_seed);
										done=1;
										chrono::high_resolution_clock::time_point start_calc_time=chrono::high_resolution_clock::now();

										Kalman(mSeed[(Int_t)mSeed.size()-1]);
										chrono::high_resolution_clock::time_point end_calc_time=chrono::high_resolution_clock::now();
										T_calc += chrono::duration_cast<chrono::microseconds>(end_calc_time-start_calc_time).count();
	
									}			
								if (done) break;
							}
						}	
				}
	//for(int i=0;i<mBins.size();i++){
		//cout<<i<<": "<<mBins[i].size()<<endl;
	//}	
}

//redict next location
void TRD_Kalman_Trackfinder::prediction(Double_t dist){
	if (dist ==0)
		return;
	
	//calculate new estimate
	Double_t f1			= 	mMu[2];
        Double_t b_fak		=	b_field*3./2000.;
		
	Double_t curvature 	= 	(2.*mMu[4]*b_fak);
	Double_t f2			= 	f1 + dist*curvature;
	Double_t r1=TMath::Sqrt((1.- f1)*(1. + f1));
	Double_t r2=TMath::Sqrt((1.- f2)*(1. + f2));
	Double_t dy2dx 		=	(f1 +f2)/(r1 +r2);
		
	mMu[0]+=dist*dy2dx;
        //if (dist<-300) printf("values: curv: %f f1: %f f2: %f r1 %f r2 %f dy2dx: %f y: %f\n",curvature,f1,f2,r1,r2,dy2dx,mMu[0]);
        if(TMath::Abs(dist*curvature)<0.05)
        {
            mMu[1]+= dist*(r2 + f2*dy2dx)*mMu[3];
        }
        else
        {
            Double_t rot=TMath::ASin(r1*f2 -r2*f1);
            if(f1*f1+f2*f2>1 && f1*f2<0)
            {          // special cases of large rotations or large abs angles
                if(f2>0)
                {
                    rot =  TMath::Pi() - rot;
                }
                else
                {
                    rot = -TMath::Pi() - rot;
                }
            }
            mMu[1]	+=mMu[3]/curvature*rot ;
        }
        mMu[2]	= f2;

		
	//Build transport matrix	
	ROOT::Math::SMatrix<double,5,5> A; //Transport Matrix
	A=ROOT::Math::SMatrixIdentity();
	
	Double_t dr	=	dist/(r1*r1*r1);
	A[0][2]	=	dr;
   	A[0][4]	=	dist * dr *b_fak;
	A[1][2]	=	A[0][2]*mMu[3]*f1;
   	A[1][3]	=	dist/r1;
   	A[1][4]	=	A[1][2]*dist *b_fak;
	A[2][4]	=	dist *b_fak*2;
	ROOT::Math::SMatrix<double,5,5> A_t=ROOT::Math::Transpose(A);
	
	if(mShow){
		//cout<<"Unc:"<<mUnc<<endl;
		cout<<"transport:"<<A<<endl;
		cout<<"f_mult:" <<A * mCov<<endl;
	}	
	
	//calculate new covariance
	mCov= A * mCov * A_t +mTau;    
	
	
	mMu_red=mObs*mMu; //made once and saved for later use 
	for( int i=0;i<4;i++)//calculate the uncertainty (useful to know)
		mUnc[i]=TMath::Sqrt(mCov[i][i] + mSig[i][i]);	
	
}


//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------

//correction as used in every lecture
void TRD_Kalman_Trackfinder::correction(ROOT::Math::SVector<double,4> measure)
{
	
	ROOT::Math::SMatrix<double,4,5> C	= mObs;
	ROOT::Math::SMatrix<double,5,4> C_t	= ROOT::Math::Transpose(C);
	ROOT::Math::SMatrix<double,5,5> P	= mCov;
	ROOT::Math::SMatrix<double,5,5> Eye	= ROOT::Math::SMatrixIdentity();
		
	mKal_Gain							= P*C_t*mCov_Res_Inv;
	
	ROOT::Math::SVector<double,4> res	=measure - C * mMu;
	mMu+=mKal_Gain*(res);
		
	mCov=(Eye -mKal_Gain*C) *P;
	res	=measure - C * mMu;
	mChi_2+=ROOT::Math::Dot(res,(mCov_Res_Inv*res));
					
	
}

//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
void TRD_Kalman_Trackfinder::Kalman(vector<Ali_TRD_ST_Tracklets*> start)
{
	{//init
		
		ROOT::Math::SVector<double,4> mes=measure(start[0]);
		for( int i=0;i<4;i++)
			mMu[i]=mes[i];
			
                mMu[4]=0.0/1.0; // q/pT  0.0/1.0
		
		mEstimate.resize(0);
		mEstimate.resize(6);
		mTrack.resize(0);
		mTrack.resize(6);
	
		mCurrent_Det=start[0]->get_TRD_det();
		//mDist=mTRD_layer_radii[5][1]-mTRD_layer_radii[mCurrent_Det%6][1];
		mEstimate[mCurrent_Det%6]=mMu;

		
		mObs		=	ROOT::Math::SMatrixIdentity();

                // 0,0 = y, 1,1=z , 2,2=sin phi , 3,3=tan lambda

		Double_t dy         = 0.2; // 0.2  0.4
		Double_t dz         = 4.0; // 4.0  4.0
		Double_t dsin_phi   = 10.0; // 7.0  10.0
		Double_t dsin_theta = 25.0; // 18.0  25.0
		Double_t dpT        = 10.0; // 10.0  10.0

		mSig[0][0]	=	dy; // 0.2
		mSig[1][1]	=	dz; // 4.0
		mSig[2][2]	=	TMath::Power(TMath::Sin(dsin_phi*TMath::Pi()/180.0),2); // 7.0
		mSig[3][3]	=	TMath::Power(TMath::Tan(dsin_theta*TMath::Pi()/180.0),2); // 18.0

		mCov		=	ROOT::Math::SMatrixIdentity();
		mCov[0][0]	=	dy; // 0.2
		mCov[1][1]	=	dz; // 4.0
		mCov[2][2]	=	TMath::Power(TMath::Sin(dsin_phi*TMath::Pi()/180.0),2);   // 7.0
		mCov[3][3]	=	TMath::Power(TMath::Tan(dsin_theta*TMath::Pi()/180.0),2); // 20.0
		//mCov[4][4]	=	0.09; // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15
     	mCov[4][4]	=	dpT*dpT; // 0.3*0.3  B -> 2.0 B 0.3 -> 0.15
	        
		mChi_2		=	0;	
		
		
		mMeasurements =new ROOT::Math::SVector<double,4>[6];
		
		for(Int_t ind=0; ind<(Int_t)start.size();ind++)
		{
			
			Int_t lay			=	start[ind]->get_TRD_det() %6;
			mMeasurements[lay] 	=	measure(start[ind]);
			TVector3 temp_vec;
			Int_t det;
			mTrack[lay]			=	new Ali_TRD_ST_Tracklets();
			temp_vec=start[ind]	->	get_TV3_offset();
			det=start[ind]		->	get_TRD_det();
			temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
			
			mTrack[lay]			->	set_TRD_det(det);
			mTrack[lay]			->	set_TV3_offset(temp_vec);
			temp_vec			=	start[ind]->get_TV3_dir();
			temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
			mTrack[lay]			->	set_TV3_dir(temp_vec);
			mTrack[lay]			->	set_TRD_index(start[ind]->get_TRD_index());
			
			if(mShow && ind==0)
			{
				cout<<"NEW TRACK"<<endl;
				cout<<"Layer:"<< lay<<endl;
				cout<<"cov:"<< mCov<<endl;
				cout<<"Mu:"<< mMu<<endl;
				cout<<"mes:"<<mMeasurements[lay]<<endl;
				//cout<<"unc"<<mUnc<<endl;			
			}
		}
		
		
		mNbr_tracklets		=	(Int_t)start.size();
					
	}
	vector<ROOT::Math::SMatrix<double, 5, 5>> cov_per_layer;
	cov_per_layer.resize(6);
	cov_per_layer[mCurrent_Det%6]=mCov;
	
	
  	for(Int_t i_start=1;i_start<(Int_t)start.size();i_start++)
   	{
		Int_t i_det   = start[i_start]->get_TRD_det();
       	Int_t i_layer = i_det%6;
                //mDist	=	mTRD_layer_radii[i_layer-1][1]-mTRD_layer_radii[i_layer][1];
		mDist	=	mTRD_layer_radii_all[i_det] - mTRD_layer_radii_all[mCurrent_Det];
		
		prediction(mDist);
		mCurrent_Det=mCurrent_Det-(mCurrent_Det%6) + i_layer;
		mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
		mCov_Res_Inv.Invert();
		
		
		if(mShow)
		{
			cout<<"Layer:"<< i_layer<<endl;
			cout<<"cov:"<< mCov<<endl;
			cout<<"Mu:"<< mMu<<endl;
			cout<<"unc"<<mUnc<<endl;	
			cout<<"mes:"<<mMeasurements[i_layer]<<endl;
				
		}
		correction(mMeasurements[i_layer]);
		cov_per_layer[i_layer]=mCov;
		
	}
	{
		Int_t i_layer=start[0]->get_TRD_det() %6;
		mMu=mEstimate[i_layer];
		mCurrent_Det=mTrack[i_layer]->get_TRD_det();
		mDist	=	mTRD_layer_radii_all[mCurrent_Det-i_layer+5]-mTRD_layer_radii_all[mCurrent_Det];
		mCov=cov_per_layer[i_layer];
		
		
	}
	Double_t chi_2_pen	=	18.5;
        //for(Int_t i_layer = 5;i_layer>=0;i_layer--) // Alex: changes where the kalman tracker starts
        for(Int_t i_layer=mCurrent_Det%6;i_layer>=0;i_layer--) // Alex: changes where the kalman tracker starts
        {
		
		if (mEstimate[i_layer]!=0)
		{
			mMu				= mEstimate[i_layer];
			mCurrent_Det	= mTrack[i_layer]->get_TRD_det();
			mDist			= mTRD_layer_radii_all[mCurrent_Det-1]-mTRD_layer_radii_all[mCurrent_Det];
			mCov			= cov_per_layer[i_layer];
			continue;
		}
		
		prediction(mDist);
		
		if(mShow)
		{
			cout<<"Layer:"<< i_layer<<endl;
			cout<<"cov:"<< mCov<<endl;
			cout<<"Mu:"<< mMu<<endl;
			cout<<"unc"<<mUnc<<endl;			
		}
		
		mCurrent_Det	= mCurrent_Det-(mCurrent_Det%6) + i_layer;
		mDist			= mTRD_layer_radii_all[mCurrent_Det-1]-mTRD_layer_radii_all[mCurrent_Det];
			
		ROOT::Math::SMatrix<double,4,4> Cov_res		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
		Cov_res.Invert();
				
		if (mMeasurements[i_layer]!=0)
		{	
			mCov_Res_Inv	= Cov_res;
			correction(mMeasurements[i_layer]);	
		}
		else
		{
			vector<Int_t> Dets;
			if (0){} //out of border
			
			else if (0){}//out of border
			
			else //in border
				Dets.push_back(mCurrent_Det-(mCurrent_Det%6) + i_layer);
			
			for(Int_t i_list_det=0;i_list_det<(Int_t)Dets.size();i_list_det++)
			{

				Int_t i_det		=	Dets[i_list_det];
				Int_t i_sector 	= 	(Int_t)(i_det/30);
        		Int_t i_stack  	= 	(Int_t)(i_det%30/6);
        		        
				if (i_sector != (Int_t)(mCurrent_Det/30)){} //change coord sys
				
				vector<Ali_TRD_ST_Tracklets*> 			found_tracklets;
				vector<Int_t> 							found_in;
				vector<Int_t> 							ind_nbr;
				vector<Double_t> 						chis;
				
				//vector<ROOT::Math::SVector<double,4>> 	resses;
				//vector<ROOT::Math::SVector<double,4>> 	meas_list;
					
				
				for(Int_t i_tracklet=0; i_tracklet<(Int_t)mBins[i_det].size();i_tracklet++){
					ROOT::Math::SVector<double,4> 	measurement	= measure(mBins[i_det][i_tracklet]);
					ROOT::Math::SVector<double,4> 	abs			= ROOT::Math::fabs(mMu_red -measurement);
					ROOT::Math::SVector<double,4> 	res			= measurement- mMu_red;
					Double_t 						chi_2		= ROOT::Math::Dot(res,(Cov_res*res));
					
					if(mShow){
						cout<<"meas"<<measurement<<", "<<mBins[i_det][i_tracklet]->get_TV3_offset()[0]<<", "<<mBins[i_det][i_tracklet]->get_TRD_index()<<endl;
						cout<<"chi:"<<chi_2<<endl;}
					//cout<<"dist to rad new:"<<mTRD_layer_radii_all[mCurrent_Det]-mBins[i_det][i_tracklet]->get_TV3_offset()[0]<<" dtor old:"<<mTRD_layer_radii[mCurrent_Det%6][1]-mBins[i_det][i_tracklet]->get_TV3_offset()[0]<<endl;
					//Bool_t fitting= ((abs[0]<2*mUnc[0]) && (abs[1]<3*mUnc[1]) && (abs[2]<2*mUnc[2]));
					if (1)// (fitting)
					{
						found_tracklets.push_back(mBins[i_det][i_tracklet]);
						found_in.push_back(i_det);
						chis.push_back(chi_2);
						//resses.push_back(res);
						//meas_list.push_back(measurement);
						ind_nbr.push_back(i_tracklet);
						
					}
				}
				
				Int_t 		min_ind	=-1;
				Double_t 	min_chi	=chi_2_pen;
				for(Int_t i=0 ;i<(Int_t)found_tracklets.size();i++)
					if (min_chi>chis[i])
					{
						min_ind=i;
						min_chi=chis[i];
					}
				
				if(min_chi<chi_2_pen)
				{ 	//new tracklet was found 
					
					if ((Int_t)(found_in[min_ind]/30) !=(Int_t)(mCurrent_Det/30)){}
						//change_coord_sys(found_in_det[min_ind]); //MuST BE IMPLEMENTET !!!!!!!!
       				
					//save tracklet info
					mVisited[found_in[min_ind]][ind_nbr[min_ind]]=1;
					mMeasurements[i_layer]=measure(found_tracklets[min_ind]);
					mCov_Res_Inv=Cov_res;
					
					
					correction(mMeasurements[i_layer]);
					
					
					Ali_TRD_ST_Tracklets* temp_tracklet	=new Ali_TRD_ST_Tracklets();
					mTrack[i_layer]						=temp_tracklet;
					
					TVector3 temp_vec;
					Int_t det							=found_in[min_ind];
					mTrack[i_layer]						->set_TRD_det(det);
				
					temp_vec							=found_tracklets[min_ind]->get_TV3_offset();
					temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
					mTrack[i_layer]						->set_TV3_offset(temp_vec);
					
					temp_vec=found_tracklets[min_ind]	->get_TV3_dir();
					temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
					mTrack[i_layer]						->set_TV3_dir(temp_vec);
					
					mTrack[i_layer]						->set_TRD_index(found_tracklets[min_ind]->get_TRD_index());
					mNbr_tracklets++;
				}
				else
				{									
					mTrack[i_layer]	=NULL;	
				}
			}
			
		}
		mEstimate[i_layer]	=mMu;
			
	}
	//if Track
        if(mNbr_tracklets>2){ // Changed from 2, Alex: 20.07.2020
		//save Track
		mFound_tracks.push_back(mTrack);
		mEstimates.push_back(mEstimate);
		
		if(mPrimVertex==1)
		{
			mDist	=	-mTRD_layer_radii_all[mCurrent_Det];
			//cout<<"\nTrack NBR: "<< mFound_tracks.size()-1<<endl;
			//cout<<"mDist after pred: "<<mDist<<endl;
			//cout<<"mMu before pred: "<<mMu<<endl;
			
			prediction(mDist);
			//cout<<"mMu after pred: "<<mMu<<endl;
			ROOT::Math::SVector<double,4> prim_vert_measurement;
			prim_vert_measurement[0]=0;
			prim_vert_measurement[1]=0;
			prim_vert_measurement[2]=mMu[2];
			prim_vert_measurement[3]=mMu[3];
			mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) + 1.0*mSig;		//Measure uncertainty Matrix
			mCov_Res_Inv.Invert();
			
			correction(prim_vert_measurement);
			//cout<<"mMu after corr: "<<mMu<<endl;
			/*
			mDist	=	 mTRD_layer_radii_all[mCurrent_Det];
			prediction(mDist);
			mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
			mCov_Res_Inv.Invert();
			if (mMeasurements[0]!=0){
				if ((Int_t)(mTrack[0]->get_TRD_det() /30) !=(Int_t)(mCurrent_Det/30)){}
				mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
				mCov_Res_Inv.Invert();
		
				correction(mMeasurements[0]);
				*/
			Double_t tempqpt=mMu[4];
			mMu=mEstimate[0];
			mMu[4]=tempqpt;
			
			//cout<<"mMu after corr2: "<<mMu<<endl;
			
				
		}	
		//used to look at specific track
		if(mFound_tracks.size()==44)mShow=1;
		if(mFound_tracks.size()==45)mShow=0;
		
		//Loop again for better fit
		for(Int_t i_layer=1;i_layer<6;i_layer++){
			mDist	=	mTRD_layer_radii_all[mCurrent_Det+1]-mTRD_layer_radii_all[mCurrent_Det];
                        //printf("%s ---> i_layer: %d %s, mCurrent_Det+1: %d, mCurrent_Det: %d \n",KGRN,i_layer,KNRM,mCurrent_Det+1,mCurrent_Det);
                        prediction(mDist);
			mCurrent_Det=mCurrent_Det-(mCurrent_Det%6) + i_layer;
			if (mMeasurements[i_layer]!=0){
				if ((Int_t)(mTrack[i_layer]->get_TRD_det() /30) !=(Int_t)(mCurrent_Det/30)){}
				mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
				mCov_Res_Inv.Invert();
				//cout<<"mMu bef corri: "<<mMu<<endl;
			
				correction(mMeasurements[i_layer]);
				//cout<<"mMu after corri: "<<mMu<<endl;
			
			}
		}
		//cout<<"mMu after first loop: "<<mMu<<endl;
			
		//Loop again for better fit
		for(Int_t i_layer=4;i_layer>=0;i_layer--){
			mDist	=	mTRD_layer_radii_all[mCurrent_Det-1]-mTRD_layer_radii_all[mCurrent_Det];
			prediction(mDist);
			mCurrent_Det=mCurrent_Det-(mCurrent_Det%6) + i_layer;
			if (mMeasurements[i_layer]!=0){
				if ((Int_t)(mTrack[i_layer]->get_TRD_det() /30) !=(Int_t)(mCurrent_Det/30)){}
				mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
				mCov_Res_Inv.Invert();
		
				correction(mMeasurements[i_layer]);
			}
		}
		
		if(mPrimVertex==1)
		{
			ROOT::Math::SVector<double, 5> mutemp=mMu;
			mDist	=	-mTRD_layer_radii_all[mCurrent_Det];
			//cout<<"mDist after pred: "<<mDist<<endl;
			//cout<<"mMu before pred: "<<mMu<<endl;
			
			prediction(mDist);
			//cout<<"mMu after pred: "<<mMu<<endl;
			ROOT::Math::SVector<double,4> prim_vert_measurement;
			prim_vert_measurement[0]=0;
			prim_vert_measurement[1]=0;
			prim_vert_measurement[2]=0;
			prim_vert_measurement[3]=mMu[3];
			
			mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) + 1.0*mSig;		//Measure uncertainty Matrix
			mCov_Res_Inv.Invert();
			
			correction(prim_vert_measurement);	
			//cout<<"mMu after corr: "<<mMu<<endl;
			/*
			mDist	=	mTRD_layer_radii_all[mCurrent_Det];
			prediction(mDist);
			mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
			mCov_Res_Inv.Invert();
			if (mMeasurements[0]!=0){
				if ((Int_t)(mTrack[0]->get_TRD_det() /30) !=(Int_t)(mCurrent_Det/30)){}
				mCov_Res_Inv		=	mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
				mCov_Res_Inv.Invert();
		
				correction(mMeasurements[0]);
			}
			
			*/
			
			Double_t tempqpt=mMu[4];
			mMu=mutemp;
			mMu[4]=tempqpt;
			
				
		}	
		
		//calculate Helix param
		Double_t charge=1.;
		if (mMu[4]<0) charge=-1.;
		Double_t lam=TMath::ATan( mMu[3]);
		Double_t pxy=charge/mMu[4] ;
		
		TVector3 x_vek;
		TVector3 p_vek;
		x_vek[0]	=	mTRD_layer_radii_all[mCurrent_Det];
		/*if(mPrimVertex==1)
			x_vek[0]=0;*/
		x_vek[1]	=	mMu[0];
		x_vek[2]	=	mMu[1];
		p_vek[0]	=	TMath::Cos(TMath::ASin(mMu[2]))*pxy;
		p_vek[1]	=	mMu[2]*pxy;
		p_vek[2]	=	mMu[3]*pxy;
		x_vek.RotateZ((Double_t)(2*(Int_t)(mCurrent_Det/30)+1)*TMath::Pi()/18);
		p_vek.RotateZ((Double_t)(2*(Int_t)(mCurrent_Det/30)+1)*TMath::Pi()/18);
		Double_t x[3];
		Double_t p[3];
		x[0]		=	x_vek[0];
		x[1]		=	x_vek[1];
		x[2]		=	x_vek[2]; 
		p[0]		=	p_vek[0];
		p[1]		=	p_vek[1];
		p[2]		=	p_vek[2];
		
		//calculation of Helixparameter taken from http://alidoc.cern.ch/AliRoot/v5-09-36/_ali_helix_8cxx_source.html 
		//AliHelix::AliHelix(Double_t x[3], Double_t p[3], Double_t charge, Double_t conversion)
		vector<Double_t> fHelix;
		fHelix.resize(8);
		Double_t pt = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
  //  	
		
		Double_t b_fak		=	b_field*3./1000.;
		
		Double_t curvature 	= 	(mMu[4]*b_fak);
	
  		fHelix[4] 	= curvature; // C
  		fHelix[3] 	= p[2]/pt;    // tgl
  //  
  		Double_t xc, yc, rc;
  		rc			=  1/fHelix[4];
		xc  		=  x[0]  -rc*p[1]/pt;
		yc  		=  x[1]  +rc*p[0]/pt; 
		  //
		fHelix[5] 	= x[0];   // x0
		fHelix[0] 	= x[1];   // y0
		fHelix[1] 	= x[2];   // z0
		  //
		//fHelix[6] = xc;
		//fHelix[7] = yc;
		//fHelix[8] = TMath::Abs(rc);
		  //
		fHelix[5] 	= xc; 
                fHelix[0] 	= yc;

                fHelix[6]       = pt;
                fHelix[7]       = p[2];
		  //
		if (TMath::Abs(p[1])<TMath::Abs(p[0])){     
			fHelix[2]=TMath::ASin(p[1]/pt);
			//Helix[2]=asinf(p[1]/pt);
			if (charge*yc<charge*x[1])  fHelix[2] = TMath::Pi()-fHelix[2];
		}
		else{
			fHelix[2]=TMath::ACos(p[0]/pt);
			//fHelix[2]=acosf(p[0]/pt);
			if (charge*xc>charge*x[0])  fHelix[2] = -fHelix[2];
		}
		mHelices.push_back(fHelix);

                /*
		cout<<mHelices.size()<<": ";
		for(int i=0;i<fHelix.size();i++)
		cout<<fHelix[i]<<" ";
                cout<<endl;
                */
	}
			
}

vector< vector<Ali_TRD_ST_Tracklets*> > TRD_Kalman_Trackfinder::Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets,Int_t Num_Tracklets,Int_t prim_vertex)
{
    mShow=0;
    mHelices.clear();
    mFound_tracks.clear();
    mEstimates.clear();
	mPrimVertex=prim_vertex;
    //cout<<"TRD_Kalman_Trackfinder::Kalman_Trackfind"<<endl;
	T_calc=0.;
	chrono::high_resolution_clock::time_point start_wall_time=chrono::high_resolution_clock::now();
				
    get_seed(Tracklets,Num_Tracklets);
    
	chrono::high_resolution_clock::time_point end_wall_time=chrono::high_resolution_clock::now();
	T_wall= chrono::duration_cast<chrono::microseconds>(end_wall_time-start_wall_time).count();
	//printf("T_wall is %f us\n",T_wall);
	//printf("T_calc is %f us\n",T_calc);
	
	/*for (int i=0; i<mSeed.size();i++){
     if (i == 55 ) mShow=1;
     Kalman(mSeed[i]);
     mShow=0;
     }*/
    //cout<<mFound_tracks[0][4]->get_TRD_index()<<endl;
    return mFound_tracks;
}

