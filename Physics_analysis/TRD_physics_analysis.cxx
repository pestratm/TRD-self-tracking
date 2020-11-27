
#include "TRD_physics_analysis.h"

// implementation of the main code

// Something like with your new classes


//----------------------------------------------------------------------------------------
Ali_TRD_physics_analysis::Ali_TRD_physics_analysis(TString out_dir, TString out_file_name)
{
	HistName = out_dir;
    HistName += "/";
    HistName += out_file_name;
    printf("test printf: %s \n",HistName.Data());
}

//----------------------------------------------------------------------------------------
void Ali_TRD_physics_analysis::Init_tree(TString SEList)
{
    printf("Ali_TRD_physics_analysis::Init_tree \n");
    TString pinputdir = input_dir;

    TString in_list_name = SEList;
    SEList = input_dir_lists + SEList;

    Kalman_Track  			= new Ali_Kalman_Track();
    TPC_Track  	   			= new Ali_TPC_Track();
    TRD_Photon     			= new Ali_TRD_Photon();
    TRD_Nuclear_interaction = new Ali_TRD_Nuclear_interaction();
    TRD_Self_Event 			= new Ali_TRD_Self_Event();

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
