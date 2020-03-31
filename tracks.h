//This programm constructs the MC simulation of the detector


//all parameters of the detector area
class Box_ID{
	public:
	Int_t layer_nbr; //number of detection layers
	Double_t box_height;
	Double_t box_length;
	Double_t box_dist;
	Double_t border_dist;
	
	Double_t layer_ineff;//percentile of inefficiency of each layer
	Double_t zero[2]; //the coordinates of the collision point in this coordinate system (left lower point of layer 0 is point 0)
	Double_t b_field; //stength of the b_field (direction in positive z- axis ) 
	TBox **box; //array containing the boxes
};

//all parameters for the track
class Particle_ID{
	public:
	Double_t track_deg_unc; //uncertainty in the direction
	Double_t track_pos_unc; //uncertainty in the position
	Double_t part_charge;  //charge of the particle in elementary charge
	Double_t min_p; //minimal momentum
	Double_t max_p;
	Int_t detected_tracks; //number of tracks that actually got detected
};	

//--------------------------------------------------------------------------------------------------

//this function just draws the detector area and returns the array containing these boxes 
TBox** draw_boxes(Box_ID& boxid){
	
	TBox **box = new TBox*[boxid.layer_nbr];
	
	for (int i = 0; i < boxid.layer_nbr; i++)
	{
		box[i]=new TBox(0,i*(boxid.box_height+boxid.box_dist),boxid.box_length,i*(boxid.box_height+boxid.box_dist)+boxid.box_height);
		box[i]->SetFillColor(7);
		box[i]->SetLineColor(1);
		box[i]->SetLineWidth(3);
		box[i]->Draw();	
		
		//for the numbering at the side
		TString nbr="";
		nbr+=i;
		TLatex *tex = new TLatex(-boxid.border_dist/2,i*(boxid.box_height+boxid.box_dist)+boxid.box_height*3/10,nbr);
   		tex->SetLineWidth(4);
   		tex->Draw();
	}
	//title
	TLatex *tex = new TLatex(boxid.box_length*(0.5 -1.5/10) ,boxid.layer_nbr*(boxid.box_height+boxid.box_dist)+boxid.border_dist*3/10,"Random Particles with Track");
   	tex->SetLineWidth(4);
   	tex->Draw();
	
	
	return box;
}

//--------------------------------------------------------------------------------------------------


//This function draws the noise tracklets and returns them in an array
//each tracklet gets assigned a layer, a x position and a direction "randomly" within the borders
TPolyLine** draw_lines(Int_t& lines_nbr,Box_ID& boxid,Double_t& rand_seed,Double_t& max_deg){
	
	TRandom3 *gRandom=new TRandom3(rand_seed);		
	
	TPolyLine** line =new TPolyLine*[lines_nbr];
	
	Double_t max_dist=TMath::Tan(max_deg)*boxid.box_height;
	
	for (int i = 0; i < lines_nbr; i++)
	{
		Double_t* x_pos=new Double_t[2];
		Double_t* y_pos=new Double_t[2];
			
		x_pos[0]= (gRandom->Rndm() )*boxid.box_length;
		Int_t layer =gRandom->Integer(boxid.layer_nbr);
		x_pos[1]= (gRandom->Rndm() )*2*max_dist +x_pos[0]-max_dist;
		
		y_pos[0]=layer*(boxid.box_height+boxid.box_dist) +0.4;
		y_pos[1]=y_pos[0]+boxid.box_height-0.8;
		
		line[i]=new TPolyLine(2,x_pos,y_pos);
		line[i]->SetLineColor(kGreen +3);
		line[i]->SetLineWidth(4);
		line[i]->Draw();	
	
	}
	return line;
}

//--------------------------------------------------------------------------------------------------
//this function creates a new track of (layer_nbr) tracklets
TPolyLine** new_track(Box_ID& boxid,Double_t rand_seed,Particle_ID& partid){
	
	TRandom3 *gRandom=new TRandom3(rand_seed);		
	
	TPolyLine** track =new TPolyLine*[boxid.layer_nbr];
	
	//random momentum within borders
	Double_t mom=(gRandom->Rndm()*(partid.max_p-partid.min_p)+partid.min_p);
	//random location within borders
	Double_t x0[2]={(gRandom->Rndm() )*boxid.box_length,0};
	
	Double_t rad_dist=mom/(boxid.b_field*partid.part_charge)*10000/3;//radial distance
	
	//just used to calculate circle_center
	Double_t diff[2]={x0[0]-boxid.zero[0],x0[1]-boxid.zero[1]};	
	Double_t diff_dist=TMath::Sqrt(diff[0]*diff[0] +diff[1]*diff[1]); 
	Double_t ortho_dist=TMath::Sqrt(rad_dist*rad_dist - diff_dist*diff_dist/4);
	Double_t ortho[2]={-diff[1]*ortho_dist/diff_dist *TMath::Sign(1,partid.part_charge),diff[0]*ortho_dist/diff_dist *TMath::Sign(1,partid.part_charge)};
	
	//center of the circle on which the particle is moving
	Double_t circle_center[2]={boxid.zero[0]+ortho[0]+diff[0]/2,boxid.zero[1]+ortho[1]+diff[1]/2};
	
	
	
	Double_t* x_pos=new Double_t[2];
	Double_t* y_pos=new Double_t[2];
	
	for (int i = 0; i < boxid.layer_nbr; i++)
	{
		
		Double_t height=i*(boxid.box_height +boxid.box_dist);
		
		
		
		y_pos[0]=height +0.4;
		y_pos[1]=height+boxid.box_height-0.8;
		
		//first (exact) point in this layer
		x_pos[0]= TMath::Sqrt(rad_dist*rad_dist -TMath::Power((y_pos[0]- circle_center[1]),2)) *TMath::Sign(1,partid.part_charge)+circle_center[0];
		// (exact) degree of second point in this layer
		Double_t x_pos2= TMath::Sqrt(rad_dist*rad_dist -TMath::Power((y_pos[1]- circle_center[1]),2))*TMath::Sign(1,partid.part_charge)+ circle_center[0];
		Double_t deg=TMath::ATan((x_pos2-x_pos[0]) /boxid.box_height);
		
		//add uncertainties randomly
		x_pos[0]+= (gRandom->Gaus(0, partid.track_pos_unc/3) );
		x_pos[1]= x_pos[0] +TMath::Tan(deg+(gRandom->Gaus(0, partid.track_deg_unc/3)))*boxid.box_height;
		
		
		
		track[i]=new TPolyLine(2,x_pos,y_pos);
		track[i]->SetLineColor(kRed );
		track[i]->SetLineWidth(4);
	}
	
	//repeat if out of bound
	if ((track[boxid.layer_nbr-1]->GetX())[track[boxid.layer_nbr-1]->GetLastPoint()]<0)
		track= new_track(boxid,gRandom->Rndm(),partid);
	if ((track[boxid.layer_nbr-1]->GetX())[track[boxid.layer_nbr-1]->GetLastPoint()]>boxid.box_length)
		track= new_track(boxid,gRandom->Rndm(),partid);
	
	partid.detected_tracks=boxid.layer_nbr;
	return track;

}

//--------------------------------------------------------------------------------------------------

//this function checks which tracklets of the track get detected
TPolyLine** detect_track(Double_t rand_seed,TPolyLine** track, Particle_ID& partid, Box_ID& boxid){
	TRandom3 *gRandom=new TRandom3(rand_seed);		
	
	TPolyLine** track2 =new TPolyLine*[boxid.layer_nbr];
	Int_t count=0;
	//write all detected tracklets in a new array
	for (int i = 0; i < partid.detected_tracks; i++){
		//check for detection
		if (gRandom->Rndm() < boxid.layer_ineff)
			continue;
		
		track2[count]=track[i];
		count++;
	}
	
	//write them in an array with the right size
	TPolyLine** track3 =new TPolyLine*[count];
	for (int i = 0; i < count; i++){
		track3[i]=track2[i];
	}
	
	//repeat if no track was detected
	if (count==0)
		track3=detect_track(gRandom->Rndm(), track , partid, boxid);
	else	
		partid.detected_tracks=count;
	
	return track3;
}

//--------------------------------------------------------------------------------------------------

//draw the track and return it
TPolyLine** draw_track(Box_ID& boxid,Double_t rand_seed,Particle_ID& partid){
	
	TPolyLine** track =new TPolyLine*[boxid.layer_nbr];
	//create track
	track= new_track(boxid,rand_seed,partid);
	//check for detection
	track=detect_track(rand_seed, track , partid, boxid);
	//draw
	for (int i = 0; i < partid.detected_tracks; i++){
		//cout<<track[i]<<endl;
		track[i]->Draw();
	}	
	return track;

}

//--------------------------------------------------------------------------------------------------
//main function
TPolyLine** tracks(Double_t* param,Box_ID& boxid){
	
	Particle_ID partid;
	
	//number of noise tracklets
	Int_t lines_nbr= static_cast<Int_t>(param[0]);
	
	boxid.layer_ineff=param[1];
	Double_t rand_seed=param[2];
	
	Double_t m_d=param[3];
	Double_t t_d=param[4];
	Double_t t_p=param[5];
	
	boxid.layer_nbr=6;
	
	//max degree of the noise tracklets
	Double_t max_deg=m_d *TMath::Pi()/180;
	
   	partid.track_deg_unc=t_d *TMath::Pi()/180;
	partid.track_pos_unc=t_p;
	
	
	partid.part_charge=param[6];//in e 
	partid.min_p=param[7]; //in GeV/c
	//partid.min_p=0.435; //in GeV/c
	partid.max_p=param[8]; //in GeV/c
	
	//1 Pixel = 1mm
	//values from ALICE data sheet
	boxid.box_height=33.5;
	boxid.box_length=1000;
	boxid.box_dist=30+3.5+47;
	boxid.border_dist=boxid.box_length/10;
	boxid.b_field=0.5; //in Tesla
	boxid.zero[0]=static_cast<Double_t>( boxid.box_length)/2;
	boxid.zero[1]=	-2947.;
	
	//make canvas
	TCanvas *c1 = new TCanvas("c1", "Tracks",1920,1080,1280,720);
   	c1->Range(-boxid.border_dist,-boxid.border_dist,boxid.box_length+ boxid.border_dist,boxid.layer_nbr*(boxid.box_height+boxid.box_dist) +boxid.border_dist);
   	
	//draw detector boxes
	boxid.box=draw_boxes(boxid);
	//draw noise and track
	TPolyLine** line= draw_lines(lines_nbr,boxid, rand_seed, max_deg);
	TPolyLine** track= draw_track(boxid,rand_seed,partid);
	
	
	//from here: create an array with all tracklets and shuffle the real track a bit under the noise tracklets
	param[0]=partid.detected_tracks+lines_nbr;
	
	TPolyLine** tracklets= new TPolyLine*[static_cast<Int_t>(param[0])];
	
	TRandom3 *gRandom=new TRandom3(rand_seed);		
	Int_t layer =gRandom->Integer(param[0]);
	Int_t rnd_tr=lines_nbr;
	Int_t tra_tr=partid.detected_tracks;
	Bool_t bol=1;
	for (int i = 0; i < param[0]; i++)
	{
		if (rnd_tr==0){
			tra_tr--;
			tracklets[i]=track[tra_tr];
			continue;
			
		}
		else if (tra_tr==0){
			rnd_tr--;
			tracklets[i]=line[rnd_tr];
			continue;
		}
		else if (i<layer){
			rnd_tr--;
			tracklets[i]=line[rnd_tr];
			continue;
		}
		else if (bol){
			tra_tr--;
			tracklets[i]=track[tra_tr];
			bol=0;
			continue;

		}
		else {
			rnd_tr--;
			tracklets[i]=line[rnd_tr];
			bol=1;
			continue;
		}	
	}
	//return it
	return tracklets;
	
}