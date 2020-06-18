#include "TRD_Kalman_Tracking.h"

vector<vector<Double_t>> TRD_Kalman_Trackfinder::get_Kalman_helix_params()
{
    return mHelices;
}

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

Bool_t TRD_Kalman_Trackfinder::fitting(Ali_TRD_ST_Tracklets* a,Ali_TRD_ST_Tracklets* b){
		//direction is not ok 
	TVector3 a_angle = (TVector3){a->get_TV3_dir()[0], a->get_TV3_dir()[1], 0};
	TVector3 b_angle = (TVector3){b->get_TV3_dir()[0], b->get_TV3_dir()[1], 0};
	
	if (abs(a_angle.Angle(b_angle)*TMath::RadToDeg() )>5)
		return 0;
	if (abs(a->get_TV3_dir().Angle(b->get_TV3_dir())*TMath::RadToDeg() )>15)
		return 0;	
		
		//position is not ok
	TVector3 z=a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] );
	if (abs((a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] ) -b->get_TV3_offset())[1])>3)
		return 0;
	if (abs((a->get_TV3_offset()+a->get_TV3_dir()*((b->get_TV3_offset()[0]- a->get_TV3_offset()[0])/a->get_TV3_dir()[0] ) -b->get_TV3_offset())[2])>10)
		return 0;
	return 1;
}

void TRD_Kalman_Trackfinder::get_seed( Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets){
		
		
	Int_t nbr_sectors=18;
	Int_t nbr_stack=5;
	Int_t nbr_layers=6; 	
	
	mBins.resize(nbr_sectors*nbr_stack*nbr_layers);
	vector< vector<Bool_t> >  visited;
	visited.resize(nbr_sectors*nbr_stack*nbr_layers);
	
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
		visited[i_det].push_back(0);
	}		
	
	Bool_t done=0;
	for(int i_sec=0;i_sec<nbr_sectors;i_sec++) 
			for(int i_stck=0;i_stck<nbr_stack;i_stck++)
				for(int i_lay=nbr_layers-1;i_lay>nbr_layers-3;i_lay--){
					Int_t i_det=i_sec*nbr_stack*nbr_layers+i_stck*nbr_layers+i_lay;
					for(int i_tracklet=0;i_tracklet<mBins[i_det].size();i_tracklet++)
						if (!(visited[i_det][i_tracklet])){
							visited[i_det][i_tracklet]=1;
							for(int j_lay=i_lay-1; j_lay>nbr_layers-4;j_lay--){
								done=0;
								Int_t j_det=i_det+j_lay -i_lay;
								for(int j_tracklet=0;j_tracklet<mBins[j_det].size();j_tracklet++)
									if(fitting(mBins[i_det][i_tracklet],mBins[j_det][j_tracklet])){
										visited[j_det][j_tracklet]=1;
										vector<Ali_TRD_ST_Tracklets*> temp_seed;
										temp_seed.resize(2);
										temp_seed[0]=mBins[i_det][i_tracklet];
										temp_seed[1]=mBins[j_det][j_tracklet];
										mSeed.push_back(temp_seed);
										done=1;
									}			
								if (done) break;
							}
						}	
				}
	for(int i=0;i<mBins.size();i++){
		//cout<<i<<": "<<mBins[i].size()<<endl;
	}	
}


void TRD_Kalman_Trackfinder::prediction(Double_t dist){
	if (dist ==0)
		return;
	Double_t f1= mMu[2];
	Double_t b_field=0.5;
	Double_t b_fak=b_field*1000./6.;
		
	Double_t curvature = (2.*mMu[4]*b_fak);
	Double_t f2= f1 + dist*curvature;
	Double_t r1=TMath::Sqrt((1.- f1)*(1. + f2));
	Double_t r2=TMath::Sqrt((1.- f1)*(1. + f2));
	Double_t dy2dx =(f1 +f2)/(r1 +r2);
		
	mMu[0]+=dist*dy2dx;
	if (TMath::Abs(dist*curvature)<0.05)
		mMu[1]+= dist*(r2 + f2*dy2dx)*mMu[3];
	else{
		Double_t rot=TMath::ASin(r1*f2 -r2*f1);
		if(f1*f1+f2*f2>1 && f1*f2<0) {          // special cases of large rotations or large abs angles
    		if (f2>0) rot =  TMath::Pi() - rot;    
    		else      rot = -TMath::Pi() - rot;
    	}
    	mMu[1]+=mMu[3]/curvature*rot ;
	}	
	mMu[2]= f2;
		
		
		
	ROOT::Math::SMatrix<double,5,5> A; //Transport Matrix
	A=ROOT::Math::SMatrixIdentity();
	
	Double_t dr=dist/(r1*r1*r1);
	A[0][2]=dr;
   	A[0][4]=dist * dr *b_fak;
	A[1][2]=A[0][2]*mMu[3]*f1;
   	A[1][3]=dist/r1;
   	A[1][4]=A[1][2]*dist *b_fak;
	A[2][4]=dist *b_fak*2;
	ROOT::Math::SMatrix<double,5,5> A_t=ROOT::Math::Transpose(A);
	mCov= A * mCov * A_t +mTau;    
		
	mMu_red=mObs*mMu;
	for( int i=0;i<4;i++)
		mUnc[i]=TMath::Sqrt(mCov[i][i] + mSig[i][i]);	
	if(mShow) cout<<mUnc<<endl;
}

void TRD_Kalman_Trackfinder::correction(ROOT::Math::SVector<double,4> measure){
	ROOT::Math::SMatrix<double,4,5> C=mObs;
	ROOT::Math::SMatrix<double,5,4> C_t=ROOT::Math::Transpose(C);
	ROOT::Math::SMatrix<double,5,5> P=mCov;
	ROOT::Math::SMatrix<double,5,5> Eye=ROOT::Math::SMatrixIdentity();
		
		
	ROOT::Math::SMatrix<double,4,4> A= C*P*C_t +mSig;
	A.Invert();
		
	ROOT::Math::SMatrix<double,5,4> K=P*C_t*A;
	
	ROOT::Math::SVector<double,4> res=measure - C * mMu;
	mMu+=K*(res);
		
	mKal_Gain=K;
	mCov=(Eye -K*C) *P;	
	
	mChi_2+=ROOT::Math::Dot(res,(A*res));
	
}

Bool_t TRD_Kalman_Trackfinder::fits(ROOT::Math::SVector<double,4> measure){
	ROOT::Math::SVector<double,4> abs=ROOT::Math::fabs(mMu_red -measure);
	Bool_t bol=abs< 2*mUnc;
	return bol;
}


void TRD_Kalman_Trackfinder::Kalman(vector<Ali_TRD_ST_Tracklets*> start){
		
	{//init
		//cout<<"a"<<endl;
		ROOT::Math::SVector<double,4> mes=measure(start[0]);
		for( int i=0;i<4;i++)
			mMu[i]=mes[i];
			
		mMu[4]=0;
		mEstimate.resize(0);		
		mEstimate.resize(6);
		//mEstimate=new ROOT::Math::SVector<double,5>[6];
		mCurrent_Det=start[0]->get_TRD_det();
		Int_t ind1=	start[0]->get_TRD_det() %6;
		Int_t ind2=	start[1]->get_TRD_det() %6;
		mDist=mTRD_layer_radii[5][0]-mTRD_layer_radii[ind1][0];
		
		mEstimate[ ind1]=mMu;

		mMeasurements =new ROOT::Math::SVector<double,4>[6];
		mMeasurements[ind1] = mes;
		
		mMeasurements[ind2] = measure(start[1]);
		mTrack.resize(0);
		mTrack.resize(6);
		TVector3 temp_vec;
		Int_t det;
		mTrack[ind1]=new Ali_TRD_ST_Tracklets();
		temp_vec=start[0]->get_TV3_offset();
		det=start[0]->get_TRD_det();
		temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
			
		mTrack[ind1]->set_TRD_det(det);
		mTrack[ind1]->set_TV3_offset(temp_vec);
		temp_vec=start[0]->get_TV3_dir();
		temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
		mTrack[ind1]->set_TV3_dir(temp_vec);
		mTrack[ind1]->set_TRD_index(start[0]->get_TRD_index());
			
		mTrack[ind2]=new Ali_TRD_ST_Tracklets();
		temp_vec=start[1]->get_TV3_offset();
		det=start[1]->get_TRD_det();
		temp_vec.RotateZ((Double_t)((2*(Int_t)(det/30))+1)*TMath::Pi()/18);	
		mTrack[ind2]->set_TRD_det(det);
		mTrack[ind2]->set_TV3_offset(temp_vec);
		temp_vec=start[1]->get_TV3_dir();
		temp_vec.RotateZ((Double_t)((2*(Int_t)(det/30))+1)*TMath::Pi()/18);
		mTrack[ind2]->set_TV3_dir(temp_vec);
		mTrack[ind2]->set_TRD_index(start[1]->get_TRD_index());
			
		mNbr_tracklets=2;
			
		mObs=ROOT::Math::SMatrixIdentity();

		mSig[0][0]=0.2;
		mSig[1][1]=16;
		mSig[2][2]=TMath::Power(TMath::Sin(6*TMath::Pi()/180),2);
		mSig[3][3]=TMath::Power(TMath::Tan(7*TMath::Pi()/180),2);

		mCov=ROOT::Math::SMatrixIdentity();
		mCov[0][0]=2;
		mCov[1][1]=20;
		mCov[2][2]=TMath::Power(TMath::Sin(7*TMath::Pi()/180),2);
		mCov[3][3]=TMath::Power(TMath::Tan(7*TMath::Pi()/180),2);
		mCov[4][4]=0.09;
		
		mChi_2=0;	
		
	}
	Double_t chi_2_pen=18.5;
	for(Int_t i_layer=5;i_layer>=0;i_layer--){
		
		if (mEstimate[i_layer]!=0){
			mMu=mEstimate[i_layer];
			mDist=mTRD_layer_radii[i_layer-1][0]-mTRD_layer_radii[i_layer][0];
		
			continue;}
		prediction(mDist);
		if(mShow)
		{
			cout<<"Layer:"<< i_layer<<endl;
			cout<<"cov:"<< mCov<<endl;
			cout<<"Mu:"<< mMu<<endl;
			cout<<"abs"<<mUnc<<endl;			
		}	
		mDist=mTRD_layer_radii[i_layer-1][0]-mTRD_layer_radii[i_layer][0];
		if (mMeasurements[i_layer]!=0)
				correction(mMeasurements[i_layer]);
		else{
			vector<Int_t> Dets;
			if (0){} //out of border
			
			else if (0){}//out of border
			
			else //in border
				Dets.push_back(mCurrent_Det-(mCurrent_Det%6) + i_layer);
			
			for(Int_t i_list_det=0;i_list_det<Dets.size();i_list_det++){

				Int_t i_det=Dets[i_list_det];
				Int_t i_sector = (Int_t)(i_det/30);
        		Int_t i_stack  = (Int_t)(i_det%30/6);
        		        
				if (i_sector != (Int_t)(mCurrent_Det/30)){} //change coord sys
				
				vector<Ali_TRD_ST_Tracklets*> found_tracklets;
				vector<Int_t> found_in;
				vector<Double_t> chis;
				
				vector<ROOT::Math::SVector<double,4>> resses;
				vector<ROOT::Math::SMatrix<double,4,4> >Cov_resses;
				vector<ROOT::Math::SVector<double,4>> meas_list;
				
				for(Int_t i_tracklet=0; i_tracklet<mBins[i_det].size();i_tracklet++){
					ROOT::Math::SVector<double,4> measurement=measure(mBins[i_det][i_tracklet]);
					ROOT::Math::SVector<double,4> abs=ROOT::Math::fabs(mMu_red -measurement);
					if(mShow) 
					cout<<"meas"<<measurement<<", "<<mBins[i_det][i_tracklet]->get_TV3_offset()[0]<<endl;
					ROOT::Math::SVector<double,4> res= measurement- mMu_red;
					ROOT::Math::SMatrix<double,4,4>  Cov_res=mObs*mCov*ROOT::Math::Transpose(mObs) +mSig;		//Measure uncertainty Matrix
					Cov_res.Invert();
					Double_t chi_2=ROOT::Math::Dot(res,(Cov_res*res));
					if(mShow)
					cout<<"chi:"<<chi_2<<endl;
					
					
					
					
					//Bool_t fitting= ((abs[0]<2*mUnc[0]) && (abs[1]<3*mUnc[1]) && (abs[2]<2*mUnc[2]));
					if (1)// (fitting)
					{
						found_tracklets.push_back(mBins[i_det][i_tracklet]);
						found_in.push_back(i_det);
						//Double_t dist=0;
						//for (int i=0;i<4;i++)
						//	dist+=abs[i];
						chis.push_back(chi_2);
						resses.push_back(res);
						Cov_resses.push_back(Cov_res);
						meas_list.push_back(measurement);
						
					}
				}
				
				Int_t min_ind=-1;
				Double_t min_chi=chi_2_pen;
				for(Int_t i=0 ;i<found_tracklets.size();i++)
					if (min_chi>chis[i]){
						min_ind=i;
						min_chi=chis[i];
					}
				if(min_chi<chi_2_pen){
					
					if ((Int_t)(found_in[min_ind]/30) !=(Int_t)(mCurrent_Det/30)){}
						//change_coord_sys(found_in_det[min_ind]); //MuST BE IMPLEMENTET !!!!!!!!
       					
					mMeasurements[i_layer]=measure(found_tracklets[min_ind]);
					correction(mMeasurements[i_layer]);
					
					if(mShow)
					{
						cout<<"Layer:"<< i_layer<<endl;
						cout<<"cov:"<< mCov<<endl;
						cout<<"Mu:"<< mMu<<endl;
						cout<<"abs"<<mUnc<<endl;			
					}	

					
					TVector3 temp_vec;
					Int_t det;
					Ali_TRD_ST_Tracklets* temp_tracklet=new Ali_TRD_ST_Tracklets();
					mTrack[i_layer]=temp_tracklet;
					temp_vec=found_tracklets[min_ind]->get_TV3_offset();
					det=found_in[min_ind];
					temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
					mTrack[i_layer]->set_TRD_det(det);
					mTrack[i_layer]->set_TV3_offset(temp_vec);
					temp_vec=found_tracklets[min_ind]->get_TV3_dir();
					temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
					mTrack[i_layer]->set_TV3_dir(temp_vec);
					mTrack[i_layer]->set_TRD_index(found_tracklets[min_ind]->get_TRD_index());
					mNbr_tracklets++;
					mChi_2+=min_chi;
				}
				else{									
					mTrack[i_layer]=NULL;
					mChi_2+=chi_2_pen;	
				}
			}
			
		}
		mEstimate[i_layer]=mMu;
			
	}
	
	if(mNbr_tracklets>2){ 
		mFound_tracks.push_back(mTrack);
		mShow=0;
		//if (mFound_tracks.size() ==-4) mShow=1;
		mEstimates.push_back(mEstimate);
		
		
		TVector3 x_vek;
		x_vek[0]=mTRD_layer_radii[0][0];
		x_vek[1]=mEstimate[0][0];
		x_vek[2]=mEstimate[0][1];
		x_vek.RotateZ((Double_t)(2*(Int_t)(mCurrent_Det/30)+1)*TMath::Pi()/18);
		Double_t charge=1.;
		if (mEstimate[0][4]<0) charge=-1.;
		Double_t lam=TMath::ATan( mEstimate[0][3]);
		Double_t pxy=TMath::Cos(lam)/mEstimate[0][4] *charge;
		TVector3 p_vek;
		p_vek[0]=TMath::Cos(TMath::ASin(mEstimate[0][2]))*pxy;
		p_vek[1]=mEstimate[0][2]*pxy;
		p_vek[2]=TMath::Sin(lam)/mEstimate[0][4] *charge;
		p_vek.RotateZ((Double_t)(2*(Int_t)(mCurrent_Det/30)+1)*TMath::Pi()/18);
		Double_t x[3];
		x[0]=x_vek[0];
		x[1]=x_vek[1];
		x[2]=x_vek[2];
		
		Double_t p[3];
		p[0]=p_vek[0];
		p[1]=p_vek[1];
		p[2]=p_vek[2];
		
		
		Double_t conversion=1;
		vector<Double_t> fHelix;
		fHelix.resize(6);
		Double_t pt = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
  //  
  		fHelix[4] = charge/(conversion*pt); // C
  		fHelix[3] = p[2]/pt;    // tgl
  //  
  		Double_t xc, yc, rc;
  		rc  =  1/fHelix[4];
		xc  =  x[0]  -rc*p[1]/pt;
		yc  =  x[1]  +rc*p[0]/pt; 
		  //
		fHelix[5] = x[0];   // x0
		fHelix[0] = x[1];   // y0
		fHelix[1] = x[2];   // z0
		  //
		//fHelix[6] = xc;
		//fHelix[7] = yc;
		//fHelix[8] = TMath::Abs(rc);
		  //
		fHelix[5]=xc; 
		fHelix[0]=yc; 
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
		
	}
			
}

vector< vector<Ali_TRD_ST_Tracklets*> > TRD_Kalman_Trackfinder::Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets,Int_t Num_Tracklets){
	
	
	mShow=0;
	get_seed(Tracklets,Num_Tracklets);
	cout<<"Nbr:"<<mSeed[0][0]<<endl;
	//mShow=1;
	for (int i=0; i<mSeed.size();i++){
		//if (i == 3 || i==4) mShow=1;
		Kalman(mSeed[i]);
		//mShow=0;
	}
	cout<<mFound_tracks[0][4]->get_TRD_index()<<endl;
	return mFound_tracks;
}

