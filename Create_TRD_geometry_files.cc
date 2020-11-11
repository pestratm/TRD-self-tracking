

void Create_TRD_geometry_files()
{
    printf("Create_TRD_geometry_files started \n");

    TFile* outputfile = new TFile("TRD_geometry_full.root","RECREATE");

    TString HistName;

    AliTRDgeometry* fGeo = new AliTRDgeometry;
    TGeoCombiTrans* combitrans[540];
    TGeoVolume *TRD_boxes[540];
    vector< vector<TH1D*> > vec_TH1D_TRD_geometry; // store for all 540 chambers the 8 corner vertices per detector
    vec_TH1D_TRD_geometry.resize(3); // x,y,z

    vector< vector<TH1D*> > vec_TH1D_TV3_TRD_center;
    vec_TH1D_TV3_TRD_center.resize(3); // x,y,z axes
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TV3_TRD_center[i_xyz].resize(3); // vector direction
        for(Int_t i_vertex = 0; i_vertex < 3; i_vertex++)
        {
            HistName = "vec_TH1D_TV3_TRD_center_";
            HistName += i_xyz;
            HistName += "_V";
            HistName += i_vertex;
            vec_TH1D_TV3_TRD_center[i_xyz][i_vertex] = new TH1D(HistName.Data(),HistName.Data(),540,0,540);
        }
    }

    vector<TH1D*> vec_TH1D_TV3_TRD_center_offset;
    vec_TH1D_TV3_TRD_center_offset.resize(3); // x,y,z axes
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        HistName = "vec_TH1D_TV3_TRD_center_offset_";
        HistName += i_xyz;
        vec_TH1D_TV3_TRD_center_offset[i_xyz] = new TH1D(HistName.Data(),HistName.Data(),540,0,540);
    }

    vector<TVector3> vec_TV3_TRD_center_offset; // 540 chambers
    vector< vector<TVector3> >     vec_TV3_TRD_center; // 540 chambers, 3 axes
    vec_TV3_TRD_center.resize(540);
    vec_TV3_TRD_center_offset.resize(540);
    for(Int_t i_det = 0; i_det < 540; i_det++)
    {
        vec_TV3_TRD_center[i_det].resize(3);
    }
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TRD_geometry[i_xyz].resize(8); // 8 vertices
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            HistName = "vec_TH1D_TRD_geometry_xyz_";
            HistName += i_xyz;
            HistName += "_V";
            HistName += i_vertex;
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] = new TH1D(HistName.Data(),HistName.Data(),540,0,540);
        }
    }


    vector<TVector3> vec_TV3_local_pos;
    vec_TV3_local_pos.resize(8);
    for(Int_t TRD_detector = 0; TRD_detector < 540; TRD_detector++)
    {
        Int_t                TRD_sector         = fGeo         ->GetSector(TRD_detector);
        Int_t                TRD_stack          = fGeo         ->GetStack(TRD_detector);
        Int_t                TRD_layer          = fGeo         ->GetLayer(TRD_detector);
        Float_t              TRD_time0          = fGeo         ->GetTime0(TRD_layer);

        // OK
        Float_t              TRD_chamber_length = fGeo->GetChamberLength(TRD_layer,TRD_stack);
        Float_t              TRD_chamber_width  = fGeo->GetChamberWidth(TRD_layer);
        Float_t              TRD_chamber_height = 8.4;

        AliTRDpadPlane*      padplane           = fGeo         ->GetPadPlane(TRD_detector);
        Double_t             TRD_col_end        = padplane     ->GetColEnd();
        Double_t             TRD_row_end        = padplane     ->GetRowEnd();            // fPadRow[fNrows-1] - fLengthOPad + fPadRowSMOffset;
        Double_t             TRD_col_start      = padplane     ->GetCol0();
        Double_t             TRD_row_start      = padplane     ->GetRow0();              // fPadRow[0] + fPadRowSMOffset
        Double_t             TRD_row_end_ROC    = padplane     ->GetRowEndROC();         // fPadRow[fNrows-1] - fLengthOPad;
        Double_t             TRD_col_spacing    = padplane     ->GetColSpacing();
        Double_t             TRD_row_spacing    = padplane     ->GetRowSpacing();

        // OK

        Double_t Rotation_angle     = ((360.0/18.0)/2.0) + ((Double_t)TRD_sector)*(360.0/18.0);

        Double_t             loc[3]           = {TRD_time0,0.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glb[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,loc,glb);


        Double_t             locZ[3]           = {TRD_time0-50.0,0.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbZ[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locZ,glbZ);

        Double_t             locX[3]           = {TRD_time0,50.0,(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbX[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locX,glbX);

        Double_t             locY[3]           = {TRD_time0,0.0,50.0+(TRD_row_end + TRD_row_start)/2.0};
        Double_t             glbY[3]           = {0.0,0.0,0.0};
        fGeo ->RotateBack(TRD_detector,locY,glbY);


        combitrans[TRD_detector] = new TGeoCombiTrans();
        combitrans[TRD_detector] ->RotateZ(Rotation_angle + 90.0);
        combitrans[TRD_detector] ->SetTranslation(glb[0],glb[1],glb[2]);


        // Not OK
        vec_TV3_local_pos[0].SetXYZ(-TRD_chamber_width/2.0,-TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[1].SetXYZ(TRD_chamber_width/2.0,-TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[2].SetXYZ(TRD_chamber_width/2.0,TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[3].SetXYZ(-TRD_chamber_width/2.0,TRD_chamber_height/2.0,-TRD_chamber_length/2.0);
        vec_TV3_local_pos[4].SetXYZ(-TRD_chamber_width/2.0,-TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[5].SetXYZ(TRD_chamber_width/2.0,-TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[6].SetXYZ(TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0);
        vec_TV3_local_pos[7].SetXYZ(-TRD_chamber_width/2.0,TRD_chamber_height/2.0,TRD_chamber_length/2.0);


        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            Double_t arr_pos_loc[3] = {vec_TV3_local_pos[i_vertex][0],vec_TV3_local_pos[i_vertex][1],vec_TV3_local_pos[i_vertex][2]};
            Double_t arr_pos_glb[3] = {0.0,0.0,0.0};
            combitrans[TRD_detector] ->LocalToMaster(arr_pos_loc,arr_pos_glb);

            vec_TH1D_TRD_geometry[0][i_vertex] ->SetBinContent(TRD_detector+1,arr_pos_glb[0]);
            vec_TH1D_TRD_geometry[1][i_vertex] ->SetBinContent(TRD_detector+1,arr_pos_glb[1]);
            vec_TH1D_TRD_geometry[2][i_vertex] ->SetBinContent(TRD_detector+1,arr_pos_glb[2]);
        }


        vec_TV3_TRD_center_offset[TRD_detector].SetXYZ(glb[0],glb[1],glb[2]);

        for(Int_t i_axis = 0; i_axis < 3; i_axis++)
        {
            vec_TH1D_TV3_TRD_center_offset[i_axis] ->SetBinContent(TRD_detector+1,vec_TV3_TRD_center_offset[TRD_detector][i_axis]);
        }

        vec_TV3_TRD_center[TRD_detector][0].SetXYZ(glbX[0]-glb[0],glbX[1]-glb[1],glbX[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][1].SetXYZ(glbY[0]-glb[0],glbY[1]-glb[1],glbY[2]-glb[2]);
        vec_TV3_TRD_center[TRD_detector][2].SetXYZ(glbZ[0]-glb[0],glbZ[1]-glb[1],glbZ[2]-glb[2]);

        for(Int_t i_axis = 0; i_axis < 3; i_axis++)
        {
            for(Int_t i_dir_component = 0; i_dir_component < 3; i_dir_component++)
            {
                vec_TH1D_TV3_TRD_center[i_axis][i_dir_component] ->SetBinContent(TRD_detector+1,vec_TV3_TRD_center[TRD_detector][i_axis][i_dir_component]);
                //printf("TRD_detector: %d, i_axis: %d, i_dir_component: %d, value: %4.3f \n",TRD_detector,i_axis,i_dir_component,vec_TV3_TRD_center[TRD_detector][i_axis][i_dir_component]);
            }
        }

    }

    outputfile ->cd();
    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        for(Int_t i_vertex = 0; i_vertex < 3; i_vertex++)
        {
            vec_TH1D_TV3_TRD_center[i_xyz][i_vertex] ->Write();
        }
    }

    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        vec_TH1D_TV3_TRD_center_offset[i_xyz] ->Write();
    }

    for(Int_t i_xyz = 0; i_xyz < 3; i_xyz++)
    {
        for(Int_t i_vertex = 0; i_vertex < 8; i_vertex++)
        {
            vec_TH1D_TRD_geometry[i_xyz][i_vertex] ->Write();
        }
    }
    outputfile ->Close();


}