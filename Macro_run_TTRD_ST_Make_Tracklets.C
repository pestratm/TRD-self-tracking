
R__LOAD_LIBRARY(TTRD_ST_Make_Tracklets_cxx.so);

void Macro_run_TTRD_ST_Make_Tracklets(TString In_list)
{
    // .L TTRD_ST_Make_Tracklets.cxx++

    //root Macro_run_TTRD_ST_Make_Tracklets.C\(\"Split_ST_digits_vD_1.546_LA_0.16133_191-195.txt\"\)
    gSystem ->Load("TTRD_ST_Make_Tracklets_cxx.so");
    //std::this_thread::sleep_for(std::chrono::milliseconds(20000));
    //Macro_run_TBase(In_list);
    Int_t graphics = 1;
    TTRD_ST_Make_Tracklets* ST_Make_Tracklets = new TTRD_ST_Make_Tracklets(graphics);
    ST_Make_Tracklets ->Init_tree(In_list.Data());

    ST_Make_Tracklets ->Loop_event(0);

    Double_t Delta_x        = 3.0;
    Double_t Delta_z        = 10.0;
    Double_t factor_layer   = 6.0;
    Double_t factor_missing = 1.0;
    ST_Make_Tracklets ->Calibrate(Delta_x,Delta_z,factor_layer,factor_missing);
    ST_Make_Tracklets -> plot_dem_histos1();
    ST_Make_Tracklets -> plot_dem_histos2();
}
