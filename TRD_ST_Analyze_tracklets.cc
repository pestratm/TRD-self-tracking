#include "TRD_ST_Analyze_tracklets.h"

void TRD_ST_Analyze_tracklets()
{
    printf("TRD_ST_Analyze_tracklets started \n");

    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze();
    TRD_ST_Analyze ->Init_tree("List_data.txt");

    Long64_t event = 3;
    TRD_ST_Analyze ->Loop_event(event);
    TRD_ST_Analyze ->Draw_event(event);
    TRD_ST_Analyze ->Do_TPC_TRD_matching(event,3.0,10.0);
    // TRD_ST_Analyze ->Do_TPC_TRD_matching_allEvents(3.0,10.0);
<<<<<<< HEAD
    //TRD_ST_Analyze ->Do_TRD_self_matching(event,10.0,45.0);

     TRD_ST_Analyze ->Draw_hist_TPC_tracklet_diffs();
=======
    TRD_ST_Analyze ->Do_TRD_self_matching(event,10.0,45.0);

    // TRD_ST_Analyze ->Draw_hist_TPC_tracklet_diffs();
>>>>>>> ea2fbd66992f6ea7d5f13b33b6b4fcfaf5b22db8
}