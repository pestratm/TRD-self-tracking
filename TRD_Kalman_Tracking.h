class Kalman_ID{
	public:
		ROOT::Math::SMatrix<double,5,5>  Cov;		//Covariance Matrix
		ROOT::Math::SMatrix<double,4,5>  Obs;		//Observation Model Matrix
		ROOT::Math::SMatrix<double,5,5>  Tau;		//Propagation uncertainty Matrix
		ROOT::Math::SMatrix<double,4,4>  Sig;		//Measure uncertainty Matrix
		ROOT::Math::SMatrix<double,5,4>   Kal_Gain; //Kalman Gain Matrix
		//Double_t curv;							//curvature
		ROOT::Math::SVector<double,5> mu; 			//current estimate
		vector< ROOT::Math::SVector<double,5> > estimate;	//estimates of this track
		vector< vector< ROOT::Math::SVector<double,5> > > estimates;	//estimates of all tracks
		ROOT::Math::SVector<double,4>* measurements;
		Double_t chi_2;
		Double_t dist;
		Int_t current_Det;
		ROOT::Math::SVector<double,4> unc;
		ROOT::Math::SVector<double,4> mu_red;
		
		Double_t TRD_layer_radii[6][2] =
    {
        {297.5,306.5},
        {310.0,320.0},
        {323.0,333.0},
        {336.0,345.5},
        {348.0,357.0},
        {361.0,371.0}
    };
		
		
		vector< vector<Ali_TRD_ST_Tracklets*> > Seed;
	
		
		vector< vector<Ali_TRD_ST_Tracklets*> >  Bins;				//Bins for all Tracklets corresponding to each Module and Layer
									//how many entries there are per bin
		vector<Ali_TRD_ST_Tracklets*> track; 				//Array with the Tracklets of the current Track
		Int_t nbr_tracklets;
		vector< vector<Ali_TRD_ST_Tracklets*> > found_tracks;		//Array with the Tracklets of all found Tracks	
		
		Bool_t show;
		
};
	
	
ROOT::Math::SVector<double,4> measure(Ali_TRD_ST_Tracklets* tracklet){
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

Bool_t fitting(Ali_TRD_ST_Tracklets* a,Ali_TRD_ST_Tracklets* b){
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

void get_seed(Kalman_ID& kalid, Ali_TRD_ST_Tracklets** Tracklets, Int_t Num_Tracklets){
		
		
	Int_t nbr_sectors=18;
	Int_t nbr_stack=5;
	Int_t nbr_layers=6; 	
	
	kalid.Bins.resize(nbr_sectors*nbr_stack*nbr_layers);
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
			
			//sort them in
		kalid.Bins[i_det].push_back(Tracklets[i_tracklet]);			
		visited[i_det].push_back(0);
	}		
	
	Bool_t done=0;
	for(int i_sec=0;i_sec<nbr_sectors;i_sec++) 
			for(int i_stck=0;i_stck<nbr_stack;i_stck++)
				for(int i_lay=nbr_layers-1;i_lay>nbr_layers-3;i_lay--){
					Int_t i_det=i_sec*nbr_stack*nbr_layers+i_stck*nbr_layers+i_lay;
					for(int i_tracklet=0;i_tracklet<kalid.Bins[i_det].size();i_tracklet++)
						if (!(visited[i_det][i_tracklet])){
							visited[i_det][i_tracklet]=1;
							for(int j_lay=i_lay-1; j_lay>nbr_layers-4;j_lay--){
								done=0;
								Int_t j_det=i_det+j_lay -i_lay;
								for(int j_tracklet=0;j_tracklet<kalid.Bins[j_det].size();j_tracklet++)
									if(fitting(kalid.Bins[i_det][i_tracklet],kalid.Bins[j_det][j_tracklet])){
										visited[j_det][j_tracklet]=1;
										vector<Ali_TRD_ST_Tracklets*> temp_seed;
										temp_seed.resize(2);
										temp_seed[0]=kalid.Bins[i_det][i_tracklet];
										temp_seed[1]=kalid.Bins[j_det][j_tracklet];
										kalid.Seed.push_back(temp_seed);
										done=1;
									}			
								if (done) break;
							}
						}	
				}
	for(int i=0;i<kalid.Bins.size();i++){
		//cout<<i<<": "<<kalid.Bins[i].size()<<endl;
	}	
}


void prediction(Kalman_ID& kalid,Double_t dist){
	if (dist ==0)
		return;
	Double_t f1= kalid.mu[2];
	Double_t b_field=0.5;
	Double_t b_fak=b_field*3. /2000.;
		
	Double_t curvature = (2.*kalid.mu[4]*b_fak);
	Double_t f2= f1 + dist*curvature;
	Double_t r1=TMath::Sqrt((1.- f1)*(1. + f2));
	Double_t r2=TMath::Sqrt((1.- f1)*(1. + f2));
	Double_t dy2dx =(f1 +f2)/(r1 +r2);
		
	kalid.mu[0]+=dist*dy2dx;
	if (TMath::Abs(dist*curvature)<0.05)
		kalid.mu[1]+= dist*(r2 + f2*dy2dx)*kalid.mu[3];
	else{
		Double_t rot=TMath::ASin(r1*f2 -r2*f1);
    	kalid.mu[1]+=kalid.mu[3]*rot /curvature;
	}	
	kalid.mu[2]= f2;
		
		
		
	ROOT::Math::SMatrix<double,5,5> A; //Transport Matrix
	A=ROOT::Math::SMatrixIdentity();
	
	Double_t dr=dist/(r1*r1*r1);
	A[0][2]=dr;
   	A[0][4]=dist * dr *b_fak;
	A[1][2]=A[0][2]*kalid.mu[3]*f1;
   	A[1][3]=dist/r1;
   	A[1][4]=A[1][2]*dist *b_fak;
	A[2][4]=dist *b_fak*2;
	ROOT::Math::SMatrix<double,5,5> A_t=ROOT::Math::Transpose(A);
	kalid.Cov= A * kalid.Cov * A_t +kalid.Tau;    
		
	kalid.mu_red=kalid.Obs*kalid.mu;
	for( int i=0;i<4;i++)
		kalid.unc[i]=TMath::Sqrt(kalid.Cov[i][i] + kalid.Sig[i][i]);	
	if(kalid.show) cout<<kalid.unc<<endl;
}

void correction(Kalman_ID& kalid,ROOT::Math::SVector<double,4> measure){
	ROOT::Math::SMatrix<double,4,5> C=kalid.Obs;
	ROOT::Math::SMatrix<double,5,4> C_t=ROOT::Math::Transpose(C);
	ROOT::Math::SMatrix<double,5,5> P=kalid.Cov;
	ROOT::Math::SMatrix<double,5,5> Eye=ROOT::Math::SMatrixIdentity();
		
		
	ROOT::Math::SMatrix<double,4,4> A= C*P*C_t +kalid.Sig;
	A.Invert();
		
	ROOT::Math::SMatrix<double,5,4> K=P*C_t*A;
	
	ROOT::Math::SVector<double,4> res=measure - C * kalid.mu;
	kalid.mu+=K*(res);
		
	kalid.Kal_Gain=K;
	kalid.Cov=(Eye -K*C) *P;	
	
	kalid.chi_2+=ROOT::Math::Dot(res,(A*res));
	
}

Bool_t fits(Kalman_ID& kalid,ROOT::Math::SVector<double,4> measure){
	ROOT::Math::SVector<double,4> abs=ROOT::Math::fabs(kalid.mu_red -measure);
	Bool_t bol=abs< 2*kalid.unc;
	return bol;
}


void Kalman(Kalman_ID& kalid,vector<Ali_TRD_ST_Tracklets*> start){
		
	{//init
		//cout<<"a"<<endl;
		ROOT::Math::SVector<double,4> mes=measure(start[0]);
		for( int i=0;i<4;i++)
			kalid.mu[i]=mes[i];
			
		kalid.mu[4]=0;
		kalid.estimate.resize(0);		
		kalid.estimate.resize(6);
		//kalid.estimate=new ROOT::Math::SVector<double,5>[6];
		kalid.current_Det=start[0]->get_TRD_det();
		Int_t ind1=	start[0]->get_TRD_det() %6;
		Int_t ind2=	start[1]->get_TRD_det() %6;
		kalid.dist=kalid.TRD_layer_radii[5][0]-kalid.TRD_layer_radii[ind1][0];
		
		kalid.estimate[ ind1]=kalid.mu;

		kalid.measurements =new ROOT::Math::SVector<double,4>[6];
		kalid.measurements[ind1] = mes;
		
		kalid.measurements[ind2] = measure(start[1]);
		kalid.track.resize(0);
		kalid.track.resize(6);
		TVector3 temp_vec;
		Int_t det;
		kalid.track[ind1]=new Ali_TRD_ST_Tracklets();
		temp_vec=start[0]->get_TV3_offset();
		det=start[0]->get_TRD_det();
		temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
			
		kalid.track[ind1]->set_TRD_det(det);
		kalid.track[ind1]->set_TV3_offset(temp_vec);
		temp_vec=start[0]->get_TV3_dir();
		temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
		kalid.track[ind1]->set_TV3_dir(temp_vec);
		kalid.track[ind1]->set_TRD_index(start[0]->get_TRD_index());
			
			
		kalid.track[ind2]=new Ali_TRD_ST_Tracklets();
		temp_vec=start[1]->get_TV3_offset();
		det=start[1]->get_TRD_det();
		temp_vec.RotateZ((Double_t)((2*(Int_t)(det/30))+1)*TMath::Pi()/18);	
		kalid.track[ind2]->set_TRD_det(det);
		kalid.track[ind2]->set_TV3_offset(temp_vec);
		temp_vec=start[1]->get_TV3_dir();
		temp_vec.RotateZ((Double_t)((2*(Int_t)(det/30))+1)*TMath::Pi()/18);
		kalid.track[ind2]->set_TV3_dir(temp_vec);
		kalid.track[ind2]->set_TRD_index(start[1]->get_TRD_index());
			
		kalid.nbr_tracklets=2;
			
		kalid.Obs=ROOT::Math::SMatrixIdentity();

		kalid.Sig[0][0]=0.09;
		kalid.Sig[1][1]=12;
		kalid.Sig[2][2]=TMath::Power(TMath::Sin(5*TMath::Pi()/180),2);
		kalid.Sig[3][3]=TMath::Power(TMath::Tan(7*TMath::Pi()/180),2);

		kalid.Cov=ROOT::Math::SMatrixIdentity();
		kalid.Cov[0][0]=2;
		kalid.Cov[1][1]=4;
		kalid.Cov[2][2]=TMath::Power(TMath::Sin(3*TMath::Pi()/180),2);
		kalid.Cov[3][3]=TMath::Power(TMath::Tan(3*TMath::Pi()/180),2);
		kalid.Cov[4][4]=0.09;
		
		kalid.chi_2=0;
	}
	Double_t chi_2_pen=17;
	for(Int_t i_layer=5;i_layer>=0;i_layer--){
		
		if (kalid.estimate[i_layer]!=0){
			kalid.mu=kalid.estimate[i_layer];
			kalid.dist=kalid.TRD_layer_radii[i_layer-1][0]-kalid.TRD_layer_radii[i_layer][0];
		
			continue;}
		prediction(kalid,kalid.dist);
		if(kalid.show)
		{
			cout<<"Layer:"<< i_layer<<endl;
			cout<<"cov:"<< kalid.Cov<<endl;
			cout<<"mu:"<< kalid.mu<<endl;
			cout<<"abs"<<kalid.unc<<endl;			
		}	
		kalid.dist=kalid.TRD_layer_radii[i_layer-1][0]-kalid.TRD_layer_radii[i_layer][0];
		if (kalid.measurements[i_layer]!=0)
				correction(kalid,kalid.measurements[i_layer]);
		else{
			vector<Int_t> Dets;
			if (0){} //out of border
			
			else if (0){}//out of border
			
			else //in border
				Dets.push_back(kalid.current_Det-(kalid.current_Det%6) + i_layer);
			
			for(Int_t i_list_det=0;i_list_det<Dets.size();i_list_det++){

				Int_t i_det=Dets[i_list_det];
				Int_t i_sector = (Int_t)(i_det/30);
        		Int_t i_stack  = (Int_t)(i_det%30/6);
        		        
				if (i_sector != (Int_t)(kalid.current_Det/30)){} //change coord sys
				
				vector<Ali_TRD_ST_Tracklets*> found_tracklets;
				vector<Int_t> found_in;
				vector<Double_t> chis;
				
				vector<ROOT::Math::SVector<double,4>> resses;
				vector<ROOT::Math::SMatrix<double,4,4> >Cov_resses;
				vector<ROOT::Math::SVector<double,4>> meas_list;
				
				for(Int_t i_tracklet=0; i_tracklet<kalid.Bins[i_det].size();i_tracklet++){
					ROOT::Math::SVector<double,4> measurement=measure(kalid.Bins[i_det][i_tracklet]);
					ROOT::Math::SVector<double,4> abs=ROOT::Math::fabs(kalid.mu_red -measurement);
					if(kalid.show) 
					cout<<"meas"<<measurement<<", "<<kalid.Bins[i_det][i_tracklet]->get_TV3_offset()[0]<<endl;
					ROOT::Math::SVector<double,4> res= measurement- kalid.mu_red;
					ROOT::Math::SMatrix<double,4,4>  Cov_res=kalid.Obs*kalid.Cov*ROOT::Math::Transpose(kalid.Obs) +kalid.Sig;		//Measure uncertainty Matrix
					Cov_res.Invert();
					Double_t chi_2=ROOT::Math::Dot(res,(Cov_res*res));
					if(kalid.show)
					cout<<"chi:"<<chi_2<<endl;
					
					
					
					
					//Bool_t fitting= ((abs[0]<2*kalid.unc[0]) && (abs[1]<3*kalid.unc[1]) && (abs[2]<2*kalid.unc[2]));
					if (1)// (fitting)
					{
						found_tracklets.push_back(kalid.Bins[i_det][i_tracklet]);
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
					
					if ((Int_t)(found_in[min_ind]/30) !=(Int_t)(kalid.current_Det/30)){}
						//change_coord_sys(kalid,found_in_det[min_ind]); //MUST BE IMPLEMENTET !!!!!!!!
       					
					kalid.measurements[i_layer]=measure(found_tracklets[min_ind]);
					correction(kalid,kalid.measurements[i_layer]);
						
					TVector3 temp_vec;
					Int_t det;
					Ali_TRD_ST_Tracklets* temp_tracklet=new Ali_TRD_ST_Tracklets();
					kalid.track[i_layer]=temp_tracklet;
					temp_vec=found_tracklets[min_ind]->get_TV3_offset();
					det=found_in[min_ind];
					temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
					kalid.track[i_layer]->set_TRD_det(det);
					kalid.track[i_layer]->set_TV3_offset(temp_vec);
					temp_vec=found_tracklets[min_ind]->get_TV3_dir();
					temp_vec.RotateZ((Double_t)(2*(Int_t)(det/30)+1)*TMath::Pi()/18);
					kalid.track[i_layer]->set_TV3_dir(temp_vec);
					kalid.track[i_layer]->set_TRD_index(found_tracklets[min_ind]->get_TRD_index());
					kalid.nbr_tracklets++;
					kalid.chi_2+=min_chi;
				}
				else{									
					kalid.track[i_layer]=NULL;
					kalid.chi_2+=chi_2_pen;	
				}
			}
			
		}
		kalid.estimate[i_layer]=kalid.mu;
			
	}
	
	if(kalid.nbr_tracklets>1){ 
		kalid.found_tracks.push_back(kalid.track);
		kalid.estimates.push_back(kalid.estimate);
	}
			
}

Kalman_ID Kalman_Trackfind(Ali_TRD_ST_Tracklets** Tracklets,Int_t Num_Tracklets){
	
	
	
	Kalman_ID kalid;
	kalid.show=0;
	get_seed(kalid,Tracklets,Num_Tracklets);
	cout<<"Nbr:"<<kalid.Seed[0][0]<<endl;
	for (int i=0; i<kalid.Seed.size();i++){
		if (i == 3 || i==4) kalid.show=1;
		Kalman(kalid,kalid.Seed[i]);
		kalid.show=0;
	}
	return kalid;
}


