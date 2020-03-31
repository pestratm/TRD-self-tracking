//class with a lot of variables
class Trackfind_ID{
	public:
	TPolyLine*** bins;//bins with all tracklets
	Int_t* bin_nbrs; //number of entries for each bin 
	Double_t** angles; //angle of each tracklet
	
	TPolyLine** track; //track which is looked at now
	TPolyLine** track_backup;
	
	Int_t longest_track;//number of tracklets which each track at least needs to have
	Int_t track_now; //number of tracklets the "track" has now
	Int_t now[2];//indices of the last verified tracklet
	
	Bool_t fitting(Int_t[],Int_t[]);//check wheter two tracks fit 
	//fitting criterias
	Double_t deg_unc;
	Double_t loc_unc;
	
	//list for found tracks
	TPolyLine*** tracklist;
	Int_t tracklist_pntr;
	Int_t tracklist_cap;
	TPolyLine** found_tracks;
	
	Bool_t** not_visited;//whether each tracklet has been visited already
	Int_t* track_lenghts;//length of each found track
	};	

//check whether two tracklets belong together based on similarity of (estimated) position and direction
Bool_t Trackfind_ID::fitting(Int_t prob[],Int_t ref[]){
	
	Double_t x_set=TMath::Tan(angles[ref[0]][ref[1]])*(bins[prob[0]][prob[1]]->GetY()[0] -bins[ref[0]][ref[1]]->GetY()[0])+ bins[ref[0]][ref[1]]->GetX()[0];
	Bool_t bol= ((TMath::Abs(x_set-bins[prob[0]][prob[1]]->GetX()[0])<loc_unc) &&(TMath::Abs(angles[ref[0]][ref[1]]-angles[prob[0]][prob[1]]) <deg_unc));	
	return bol; 
}

//sort each bin according to the x position of the tracklets
void sort_bins(TPolyLine** bin,Int_t bin_nbr){
	TPolyLine* temp;
	for (int i = 1; i < bin_nbr; i++){
		temp = bin[i];
      	Int_t j = i;
      	while ((j > 0) && (bin[j-1]->GetX()[0] > temp->GetX()[0])){
           bin[j] = bin[j - 1];
           j = j - 1;
      	}
      	bin[j] = temp;
	}
}

//sort the tracklets in (layer_nbr) bins according to their y position 
void sort_tracklets(TPolyLine** tracklets, Int_t nbr_tracklets, Box_ID& boxid,Trackfind_ID& findid){
	
	Double_t center;
	
	for (int j = 0; j < boxid.layer_nbr; j++)
		findid.bin_nbrs[j]=0;
	
	for (int i = 0; i < nbr_tracklets; i++){
		center= (tracklets[i]->GetY()[0] + tracklets[i]->GetY()[1])/2;
		for (int j = 0; j < boxid.layer_nbr; j++){
			if (center>boxid.box[j]->GetY1() && center<boxid.box[j]->GetY2()){
				findid.bin_nbrs[j]++;
				break;
			}
		}
	}
	
		
	Int_t* bin_nbrs2= new Int_t[boxid.layer_nbr];
	for (int j = 0; j < boxid.layer_nbr; j++){	
		findid.bins[j]=new TPolyLine*[findid.bin_nbrs[j]];
		bin_nbrs2[j]=findid.bin_nbrs[j];
	}	
	
	for (int i = 0; i < nbr_tracklets; i++){
		center= (tracklets[i]->GetY()[0] + tracklets[i]->GetY()[1])/2;
		for (int j = 0; j < boxid.layer_nbr; j++){
			if (center>boxid.box[j]->GetY1() && center<boxid.box[j]->GetY2()){
				bin_nbrs2[j]--;
				findid.bins[j][bin_nbrs2[j]]=tracklets[i];
				break;
			}
		}
	}
	
	for (int j = 0; j < boxid.layer_nbr; j++)
		sort_bins(findid.bins[j],findid.bin_nbrs[j]);
}


//calculate the angle of a tracklet by it's x and y positions
void calc_angles(Trackfind_ID& findid,Box_ID& boxid){
	findid.angles= new Double_t*[boxid.layer_nbr];
	for (int j = 0; j < boxid.layer_nbr; j++){
		findid.angles[j]=new Double_t[findid.bin_nbrs[j]];
		for (int i = 0; i < findid.bin_nbrs[j]; i++){
			findid.angles[j][i]=TMath::ATan((findid.bins[j][i]->GetX()[1]-findid.bins[j][i]->GetX()[0]) /(findid.bins[j][i]->GetY()[1]-findid.bins[j][i]->GetY()[0]));
		}
	}	
	
}



//gets a startpoint-tracklet and looks for all other members who fit this tracklet
void walk_track(Box_ID& boxid, Trackfind_ID& findid,Int_t start_at,Int_t delete_to){
	
	for (int j = start_at; j < boxid.layer_nbr; j++){//for every layer
		Int_t count=0;
					
		for (int i = 0; i < findid.bin_nbrs[j]; i++){//get the number of fitting tracklets
			
			Int_t now[2]={j,i};//the tracklet which is looked at
			if (findid.fitting(now,findid.now)){
				
				count++;
				findid.not_visited[j][i]=0;//set it to visited
			}	
		}
		//only one fitting tracklet was found
		if (count==1){
		
			for (int i = 0; i < findid.bin_nbrs[j]; i++){ 
				Int_t now[2]={j,i};
				if (findid.fitting(now,findid.now)){//find it
					//add it to the track and update the last found tracklet
					findid.track[j]=	findid.bins[j][i];
					findid.now[0]=j;
					findid.now[1]=i;
					findid.track_now++;
					break;
				}
			}	
		}
		//more than one tracklet has been found
		if (count>1){
			//find them and put them in a queue
			Int_t* indices= new Int_t[count];
			Int_t stack_cnt=0;
			for (int i = 0; i < findid.bin_nbrs[j]; i++){ 
				Int_t now[2]={j,i};
				if (findid.fitting(now,findid.now)){
					indices[stack_cnt]=i;
					stack_cnt++;
					
				}
			}
			//for each found tracklet we add it to the track and recall this function but starting at one layer above
			Int_t track_now=findid.track_now;
			for (int i = 0; i < count; i++){
				
				findid.track_now=track_now;
				findid.now[0]=j;
				findid.now[1]=indices[i];
				findid.track[j]=	findid.bins[j][indices[i]];
				findid.track_now++;
				if (i<count-1)
					walk_track(boxid, findid, j+1,j);
		
			}
		}
		
	}
	//if the track fulfills the length requirement and the list is not full write it into it

	if(findid.track_now >= findid.longest_track && findid.tracklist_pntr<findid.tracklist_cap){
		for (int j = 0; j < boxid.layer_nbr; j++){	
			//findid.track_backup[j]=findid.track[j];
			findid.tracklist[findid.tracklist_pntr][j]=findid.track[j];
		}
		findid.track_lenghts[findid.tracklist_pntr]=findid.track_now;
		findid.tracklist_pntr++;	
	}
	
	//reset the list for the up to date track	
	for (int j = delete_to; j < boxid.layer_nbr; j++){	
		findid.track[j]=0;
	}	
		
}

//just paint all tracks that have been found
void paint_track(Trackfind_ID& findid,Box_ID& boxid){
	findid.found_tracks=new TPolyLine*[findid.tracklist_pntr];
	for(int i=0; i<findid.tracklist_pntr;i++){
		
		Double_t *x_vals=new Double_t[findid.track_lenghts[i]];
		Double_t *y_vals=new Double_t[findid.track_lenghts[i]];
		Int_t k=0;	
		for (int j = 0; j < boxid.layer_nbr; j++){
			if (findid.tracklist[i][j]!=0){
			Double_t center_y= (findid.tracklist[i][j]->GetY()[0] + findid.tracklist[i][j]->GetY()[1])/2;
			Double_t center_x= (findid.tracklist[i][j]->GetX()[0] + findid.tracklist[i][j]->GetX()[1])/2;
			x_vals[j+k]=center_x;
			y_vals[j+k]=center_y;		
			//cout<<center_x<<" "<<center_y<<" "<<j+k<<findid.longest_track<<endl;	
			}
			else
				k--;
				
			
		}
		
		findid.found_tracks[i]=new TPolyLine(findid.track_lenghts[i],x_vals,y_vals);
		findid.found_tracks[i]->SetLineColor(kMagenta -4);
		findid.found_tracks[i]->SetLineWidth(2);
		findid.found_tracks[i]->Draw();	
			
		
	}
	
	
}
//function which finds all tracks
void det_track(Trackfind_ID& findid,Box_ID& boxid){
	//calculate the angle of each track
	calc_angles(findid,boxid);
	
	//create list for track and list for all found tracks
	findid.track =new TPolyLine*[boxid.layer_nbr];
	findid.track_backup =new TPolyLine*[boxid.layer_nbr];
	findid.tracklist_cap=20;
	findid.tracklist_pntr=0;
	findid.tracklist=new TPolyLine**[findid.tracklist_cap];
	findid.track_lenghts=new Int_t[findid.tracklist_cap];  
	findid.not_visited=new Bool_t*[boxid.layer_nbr];
	
	//initiate all values
	for (int j = 0; j < findid.tracklist_cap; j++){	
		findid.tracklist[j]=new TPolyLine*[boxid.layer_nbr];
		//findid.track[j]=0;
	}

	for (int j = 0; j < boxid.layer_nbr; j++){	
		findid.track_backup[j]=0;
		findid.track[j]=0;
		findid.not_visited[j]=new Bool_t[findid.bin_nbrs[j]];
		for (int i = 0; i < findid.bin_nbrs[j]; i++)
			findid.not_visited[j][i]=1;
	}

	//minimum lenght of track
	findid.longest_track=3;
	
	//actual lenght of track
	findid.track_now=0;
	
	//for each tracklet starting bottom left:
	for (int j = 0; j < boxid.layer_nbr-findid.longest_track; j++){
		for (int i = 0; i < findid.bin_nbrs[j]; i++){
			//if it hasn't been visited 
			if (findid.not_visited[j][i]){
				
				//set it as start of the track and walk along the track to find all other members
				findid.track_now=1;
				findid.now[0]=j;
				findid.now[1]=i;
				findid.track[j]=	findid.bins[j][i];					
				walk_track(boxid, findid, j+1,0);
			}
		}	
	}
	//paint it
	paint_track(findid,boxid);
		
}



//find all tracks
void find_track(TPolyLine** tracklets, Int_t nbr_tracklets, Box_ID& boxid,Double_t* param){
	Trackfind_ID findid;
	//create bins
	findid.bins= new TPolyLine**[boxid.layer_nbr];
	findid.bin_nbrs= new Int_t[boxid.layer_nbr];
	findid.deg_unc=param[0]*TMath::Pi();
	findid.loc_unc=param[1]*boxid.box_length;
	//sort the tracklets
	sort_tracklets( tracklets, nbr_tracklets, boxid, findid);
	//find the tracks
	det_track(findid,boxid);	
	
}






