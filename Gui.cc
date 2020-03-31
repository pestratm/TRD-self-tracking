# include "tracks.h"
# include "find_track.h"

//This GUI is for setting parameters and delivering data from one funtion to the other
//not worth looking through

class MyMainFrame {
   	RQ_OBJECT("MyMainFrame")
private:
	const Int_t nbr_param =9;
   	const Int_t nbr_param2 =2;
   	TGMainFrame          *fMain;
	TGHorizontalFrame 	 *fLeader;
	TGLayoutHints        *fLeaderL;
   	
	TGVerticalFrame      *fF1;
   	TGVerticalFrame      *fF2;
   	TGHorizontalFrame   **fF;
   	TGLayoutHints        *fL1;
   	TGLayoutHints        *fL2;
   	TGLabel             **fLabel;
   	TGNumberEntry       **fNumericEntries;
	TGVerticalFrame      *fFh1;
   	TGVerticalFrame      *fFh2;
   	TGHorizontalFrame   **fFh;
   	TGLayoutHints        *fLh1;
   	TGLayoutHints        *fLh2;
   	TGLabel             **fLabel2;
   	TGNumberEntry       **fNumericEntries2;
	TPolyLine			**tracklets;
	Int_t 				  num_tracklets;
	Box_ID 				  boxid;
public:
   	MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
   	virtual ~MyMainFrame();
	void Draw();
	void Find();
};

MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h) {
   // Create a main frame
	
   	fMain = new TGMainFrame(p,w,h);
	
	//Variable deklaration
	const Double_t numinit[] = {20, 0.1, 60,20.,2.,2.,1,0.5,1.5};
	const char *numlabel[] = {"lines_nbr","layer_ineff","rand_seed", "max_deg" ,
							  "max_track_deg","max_track_pos","part_charge","min_p","max_p"};
	
	
	fLeader=new TGHorizontalFrame(fMain,200,300);
	fLeaderL = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
   	fMain->AddFrame(fLeader, fLeaderL);
	
	fF=new TGHorizontalFrame*[nbr_param];
	fLabel=new TGLabel*[nbr_param];
	fNumericEntries=new TGNumberEntry*[nbr_param];
	fF1 = new TGVerticalFrame(fLeader, 200, 300);
   	fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
   	fLeader->AddFrame(fF1, fL1);
   	fL2 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2);
   

	for (Int_t i = 0; i < nbr_param; i++) {
		fF[i] = new TGHorizontalFrame(fF1, 200, 30);
      	fF1->AddFrame(fF[i], fL2);
		
		fNumericEntries[i] = new TGNumberEntry(fF[i], numinit[i], 12, i + 20,(TGNumberFormat::EStyle) 5);
      	
		fF[i]->AddFrame(fNumericEntries[i], fL2);
      	fLabel[i] = new TGLabel(fF[i], numlabel[i]);
      	fF[i]->AddFrame(fLabel[i], fL2);
   	}
	

	const Double_t numinit2[] = {0.025, 0.015};
	const char *numlabel2[] = {"deg_unc","log_unc"};
	
	
	fFh=new TGHorizontalFrame*[nbr_param2];
	fLabel2=new TGLabel*[nbr_param2];
	fNumericEntries2=new TGNumberEntry*[nbr_param2];
	
	fFh1 = new TGVerticalFrame(fLeader, 200, 500);
   	fLh1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
   	fLeader->AddFrame(fFh1, fLh1);
	fLh2 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2);
   
   	for (Int_t i = 0; i < nbr_param2; i++) {
		fFh[i] = new TGHorizontalFrame(fFh1, 200, 30);
      	fFh1->AddFrame(fFh[i], fLh2);
		
		fNumericEntries2[i] = new TGNumberEntry(fFh[i], numinit2[i], 12, i + 20,(TGNumberFormat::EStyle) 5);
      	
		fFh[i]->AddFrame(fNumericEntries2[i], fLh2);
      	fLabel2[i] = new TGLabel(fFh[i], numlabel2[i]);
      	fFh[i]->AddFrame(fLabel2[i], fLh2);
   	}
	
	
	//Buttons
	TGHorizontalFrame *hframe = new TGHorizontalFrame(fMain,200,40);
   	TGTextButton *draw = new TGTextButton(hframe,"&Draw_next");
   	draw->Connect("Clicked()","MyMainFrame",this,"Draw()");
   	hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
	
	TGTextButton *find = new TGTextButton(hframe,"&Find_track");
   	find->Connect("Clicked()","MyMainFrame",this,"Find()");
   	hframe->AddFrame(find, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
	
	TGTextButton *exit = new TGTextButton(hframe,"&Exit","gApplication->Terminate(0)");
   	hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
   	fMain->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

	fMain->MapSubwindows();
	fMain->Resize(fMain->GetDefaultSize());
	fMain->MapWindow();
	
	num_tracklets=0;
}

MyMainFrame::~MyMainFrame() {
   fMain->Cleanup();
   	delete fMain;
}

void MyMainFrame::Draw() {
	
	
	Double_t *param=new Double_t[nbr_param];
	for (Int_t i = 0; i < nbr_param; i++) 
		param[i]=fNumericEntries[i]->GetNumber();
	
	tracklets=tracks(param,boxid);
	fNumericEntries[2]->SetNumber(fNumericEntries[2]->GetNumber()+1);
	num_tracklets=static_cast<Int_t>(param[0]);

}

void MyMainFrame::Find() {
	Double_t *param=new Double_t[nbr_param2];
	for (Int_t i = 0; i < nbr_param2; i++) 
		param[i]=fNumericEntries2[i]->GetNumber();

	if (num_tracklets>0)
		find_track(tracklets,num_tracklets, boxid,param);
}

void Gui	() {
   // Popup the GUI...
   	new MyMainFrame(gClient->GetRoot(),200,800);
}

