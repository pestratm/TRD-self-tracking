

#include "TRD_ST_Analyze_tracklets.h"

void TRD_ST_Analyze_tracklets()
{
    printf("TRD_ST_Analyze_tracklets started \n");
    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze();
    TRD_ST_Analyze ->Init_tree("List_data.txt");
    TRD_ST_Analyze ->Loop_event(0);
    TRD_ST_Analyze ->Draw_event(0);
    TRD_ST_Analyze ->Do_TPC_TRD_matching(0,3.0,10.0);

}