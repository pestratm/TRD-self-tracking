#ifndef __ALI_TRD_SELF_EVENT_H__
#define __ALI_TRD_SELF_EVENT_H__

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

//----------------------------------------------------------------------------------------
class Ali_Kalman_Track : public TObject
{
private:
    // Track properties
    
    TLorentzVector TLV_part; // Lorentz vector properties of this particle
    ULong64_t      TRD_ADC_time_layer[6];
    Float_t        aliHelix_params[6]; //need
    Double_t       Chi_2;
    //bitTRDlayer

public:
    Ali_Kalman_Track() :
    TLV_part(),TRD_ADC_time_layer(),aliHelix_params(),Chi_2(-1)
    {
        
    }
    ~Ali_Kalman_Track()
    {
           
    }

    // setters
    
    void set_TLV_part(TLorentzVector tlv)     { TLV_part = tlv; }
    
    void setTRD_layer(Int_t i_layer, ULong64_t l)  { TRD_ADC_time_layer[i_layer] = l;}
       
    void setKalmanHelix_param(Float_t a, Float_t b,Float_t c,Float_t d,Float_t e,Float_t f)
    {
        aliHelix_params[0] = a;
        aliHelix_params[1] = b;
        aliHelix_params[2] = c;
        aliHelix_params[3] = d;
        aliHelix_params[4] = e;
        aliHelix_params[5] = f;
    }

    void set_Chi2(Double_t chi2)  { Chi_2 = chi2;}

    // getters
   
    TLorentzVector get_TLV_part() const       { return TLV_part;   }
    
    ULong64_t getTRD_layer(Int_t i_layer) const   { return TRD_ADC_time_layer[i_layer]; }
        
    Float_t   getKalmanHelix_param(Int_t i_param) const              {return aliHelix_params[i_param]; }

    Double_t get_Chi2() const { return Chi_2;}
    /*
    Float_t   getTRD_ADC(Int_t i_layer, Int_t i_time_bin) const
    {
        if(i_layer < 0 || i_layer > 5 || i_time_bin < 0 || i_time_bin > 7) return -1; // out of range
        ULong64_t TRD_value = 0;
        for(Int_t i_bit = 0; i_bit < 8; i_bit++) // One single time slice 8 bit = 256
        {
        Int_t bitcheck = i_bit + 8*i_time_bin; // range: 0..63 = 64 bit = 8 byte = Long64_t
        Int_t bit_status = (TRD_ADC_time_layer[i_layer] >> bitcheck) & 1; // check bit bitcheck
        if(bit_status) TRD_value |= (ULong64_t)1 << i_bit; // setting bit i_bit to 1
        }
        Float_t TRD_value_decode = (Float_t)TRD_value * 100.0; // * TRD_ADC_bin_width
        return TRD_value_decode;
    }
    */



    ClassDef(Ali_Kalman_Track,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
class Ali_TPC_Track : public TObject
{
private:
    // Track properties
    Float_t        nsigma_e_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_e_TOF; // nsigma TOF of particle
    Float_t        nsigma_pi_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_pi_TOF; // nsigma TOF of particle
    Float_t        nsigma_K_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_K_TOF; // nsigma TOF of particle
    Float_t        nsigma_p_TPC; // nsigma dE/dx of particle
    Float_t        nsigma_p_TOF; // nsigma TOF of particle
    Float_t        dca; // distance of closest approach of particle A
    TLorentzVector TLV_part; // Lorentz vector properties of this particle
    UShort_t       NTPCcls; // Number of TPC clusters
    UShort_t       NITScls; // Number of TRD clusters
    UShort_t       status; // status of track: bit 0: ITS refit, bit1: TPC refit
    Float_t        TPCchi2; // TPC chi2
    Float_t        TPCdEdx; // Energy loss information of TPC
    Float_t        TOFsignal; // Time-of-flight
    Float_t        Track_length; // length of track
    Float_t        aliHelix_params[9];

public:
    Ali_TPC_Track() :
    nsigma_e_TPC(-1),nsigma_e_TOF(-1),nsigma_pi_TPC(-1),nsigma_pi_TOF(-1),nsigma_K_TPC(-1),nsigma_K_TOF(-1),nsigma_p_TPC(-1),nsigma_p_TOF(-1),dca(-1),TLV_part(),NTPCcls(-1),NITScls(-1),status(-1),TPCchi2(-1),
        TPCdEdx(-1),TOFsignal(-1),Track_length(-1),aliHelix_params()
    {}

    ~Ali_TPC_Track() 
    {}


    // setters
    void setnsigma_e_TPC(Float_t f)                     { nsigma_e_TPC = f;         }
    void setnsigma_e_TOF(Float_t f)                     { nsigma_e_TOF = f;         }
    void setnsigma_pi_TPC(Float_t f)                     { nsigma_pi_TPC = f;         }
    void setnsigma_pi_TOF(Float_t f)                     { nsigma_pi_TOF = f;         }
    void setnsigma_K_TPC(Float_t f)                     { nsigma_K_TPC = f;         }
    void setnsigma_K_TOF(Float_t f)                     { nsigma_K_TOF = f;         }
    void setnsigma_p_TPC(Float_t f)                     { nsigma_p_TPC = f;         }
    void setnsigma_p_TOF(Float_t f)                     { nsigma_p_TOF = f;         }
    void setdca(Float_t f)                    { dca = f;        }
    void set_TLV_part(TLorentzVector tlv)     { TLV_part = tlv; }
    void setNTPCcls(UShort_t s)               { NTPCcls = s;}
    void setNITScls(UShort_t s)               { NITScls = s;}
    void setStatus(UShort_t s)                { status = s;}
    void setTPCchi2(Float_t f)                { TPCchi2 = f;}
    void setTPCdEdx(Float_t f)                       {TPCdEdx = f;}
    void setTOFsignal(Float_t f)                     {TOFsignal = f;}
    void setTrack_length(Float_t f)                  {Track_length = f;}
    void setHelix(Float_t a, Float_t b,Float_t c,Float_t d,Float_t e,Float_t f)
        {
            aliHelix_params[0] = a;
            aliHelix_params[1] = b;
            aliHelix_params[2] = c;
            aliHelix_params[3] = d;
            aliHelix_params[4] = e;
            aliHelix_params[5] = f;

        }

    // getters
    Float_t getnsigma_e_TPC() const                     { return nsigma_e_TPC;         }
    Float_t getnsigma_e_TOF() const                     { return nsigma_e_TOF;         }
    Float_t getnsigma_pi_TPC() const                     { return nsigma_pi_TPC;         }
    Float_t getnsigma_pi_TOF() const                     { return nsigma_pi_TOF;         }
    Float_t getnsigma_K_TPC() const                     { return nsigma_K_TPC;         }
    Float_t getnsigma_K_TOF() const                     { return nsigma_K_TOF;         }
    Float_t getnsigma_p_TPC() const                     { return nsigma_p_TPC;         }
    Float_t getnsigma_p_TOF() const                     { return nsigma_p_TOF;         }
    Float_t getdca() const                    { return dca;        }
    TLorentzVector get_TLV_part() const       { return TLV_part;   }
    UShort_t getNTPCcls() const               { return NTPCcls;    }
    UShort_t getNITScls() const               { return NITScls;    }
    UShort_t getStatus() const               { return status;    }
    Float_t  getTPCchi2() const              { return TPCchi2; }
    Float_t   getTPCdEdx() const                { return TPCdEdx; }
    Float_t   getTOFsignal() const              { return TOFsignal; }
    Float_t   getTrack_length() const           { return Track_length; }
    Float_t   getHelix_param(Int_t i_param) const              {return aliHelix_params[i_param]; }


    //Float_t   getTRD_ADC(Int_t i_layer, Int_t i_time_bin) const
    //{
    //    if(i_layer < 0 || i_layer > 5 || i_time_bin < 0 || i_time_bin > 7) return -1; // out of range
    //    ULong64_t TRD_value = 0;
    //    for(Int_t i_bit = 0; i_bit < 8; i_bit++) // One single time slice 8 bit = 256
    //    {
    //  Int_t bitcheck = i_bit + 8*i_time_bin; // range: 0..63 = 64 bit = 8 byte = Long64_t
    //  Int_t bit_status = (TRD_ADC_time_layer[i_layer] >> bitcheck) & 1; // check bit bitcheck
    //  if(bit_status) TRD_value |= (ULong64_t)1 << i_bit; // setting bit i_bit to 1
    //    }
    //    Float_t TRD_value_decode = (Float_t)TRD_value * 100.0; // * TRD_ADC_bin_width
    //    return TRD_value_decode;
    //}

    Int_t     HasITShit_on_layer(Int_t ilayer) { return ((NITScls >> ilayer) & 1);}  // ITShit -> LOL

    ClassDef(Ali_TPC_Track,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
class Ali_TRD_Photon : public TObject
{
private:
                                
    Float_t vertex_point[3];
    Float_t bit_TRD_layer_shared;
    Float_t pT_AB;
    Float_t AP_pT;
    Float_t AP_alpha;
    Float_t dca_min;
    Float_t path_min;
    Float_t Inv_mass_AB;
    Float_t Eta_AB;
    Float_t Phi_AB;
    Float_t dot_product_dir_vertex;
    Float_t Inv_mass_AB_K0s;
    Float_t dcaAB;
    Float_t Inv_mass_AB_Lambda;
    Float_t Inv_mass_AB_antiLambda;

    UShort_t fNumTPC_Tracks;
    UShort_t fNumKalman_Tracks;

    TLorentzVector TLV_part_A; // Lorentz vector properties of this particle
    TLorentzVector TLV_part_B; // Lorentz vector properties of this particle

    TClonesArray* fTPC_Tracks;              //->
    TClonesArray* fKalman_Tracks;           //->

public:
    Ali_TRD_Photon() :
        vertex_point(),bit_TRD_layer_shared(-1),pT_AB(-1),AP_pT(-1),AP_alpha(-1),dca_min(-1),path_min(-1),Inv_mass_AB(-1),Eta_AB(-1),Phi_AB(-1),
        dot_product_dir_vertex(-1),Inv_mass_AB_K0s(-1),dcaAB(-1),Inv_mass_AB_Lambda(-1),Inv_mass_AB_antiLambda(-1),
        fNumTPC_Tracks(0),fNumKalman_Tracks(0),TLV_part_A(),TLV_part_B()
    {
        fTPC_Tracks             = new TClonesArray( "Ali_TPC_Track", 10 );
        fKalman_Tracks          = new TClonesArray( "Ali_Kalman_Track", 10 );
    }
	~Ali_TRD_Photon() 
    {
        delete fTPC_Tracks;
        fTPC_Tracks = NULL;
        delete fKalman_Tracks;
        fKalman_Tracks = NULL;
    }

        // setters
        void set_vertex_point(Float_t x, Float_t y, Float_t z)         {    vertex_point[0] = x; vertex_point[1] = y; vertex_point[2] = z;  }
        void set_bit_TRD_layer_shared(Float_t bit_layer_shared)        {    bit_TRD_layer_shared = bit_layer_shared;                        }
        void set_pT_AB(Float_t pt_ab)                                  {    pT_AB = pt_ab;                                                  }
        void set_AP_pT(Float_t ap_pt)                                  {    AP_pT = ap_pt;                                                  }
        void set_AP_alpha(Float_t ap_alpha)                            {    AP_alpha = ap_alpha;                                            }
        void set_dca_min(Float_t dca)                                  {    dca_min = dca;                                                  }
        void set_path_min(Float_t path)                                {    path_min = path;                                                }
        void set_Inv_mass_AB(Float_t mass_ab)                          {    Inv_mass_AB = mass_ab;                                          }
        void set_Eta_AB(Float_t eta_ab)                                {    Eta_AB = eta_ab;                                                }
        void set_Phi_AB(Float_t phi_ab)                                {    Phi_AB = phi_ab;                                                }
        void set_dot_product_dir_vertex(Float_t product)               {    dot_product_dir_vertex = product;                               }
        void set_Inv_mass_AB_K0s(Float_t mass_k0)                      {    Inv_mass_AB_K0s = mass_k0;                                      }
        void set_dcaAB(Float_t dca_ab)                                 {    dcaAB = dca_ab;                                                 }
        void set_Inv_mass_AB_Lambda(Float_t mass_l)                    {    Inv_mass_AB_Lambda = mass_l;                                    }
        void set_Inv_mass_AB_antiLambda(Float_t mass_al)               {    Inv_mass_AB_antiLambda = mass_al;                               }
        void set_TLV_part_A(TLorentzVector tlv)                        {    TLV_part_A = tlv;                                               }
        void set_TLV_part_B(TLorentzVector tlv)                        {    TLV_part_B = tlv;                                               }


        // getters
        Float_t get_vertex_point(Int_t i_xyz) const    {    return vertex_point[i_xyz];         }
        Float_t get_bit_TRD_layer_shared() const       {    return bit_TRD_layer_shared; }
        Float_t get_pT_AB() const                      {    return pT_AB;   }
        Float_t get_AP_pT() const                      {    return AP_pT; }
        Float_t get_AP_alpha() const                   {    return AP_alpha; }
        Float_t get_dca_min() const                                 {    return dca_min; }
        Float_t get_path_min() const                               {    return path_min; }
        Float_t get_Inv_mass_AB() const                         {    return Inv_mass_AB; }
        Float_t get_Eta_AB() const                               {    return Eta_AB; }
        Float_t get_Phi_AB() const                               {    return Phi_AB;  }
        Float_t get_dot_product_dir_vertex() const              {    return dot_product_dir_vertex; }
        Float_t get_Inv_mass_AB_K0s() const                     {    return Inv_mass_AB_K0s; }
        Float_t get_dcaAB() const                                {    return dcaAB; }
        Float_t get_Inv_mass_AB_Lambda() const                   {    return Inv_mass_AB_Lambda; }
        Float_t get_Inv_mass_AB_antiLambda() const              {    return Inv_mass_AB_antiLambda; }
        TLorentzVector get_TLV_part_A() const       { return TLV_part_A;   }
        TLorentzVector get_TLV_part_B() const       { return TLV_part_B;   }

        //-----------------------------------

        Ali_TPC_Track* createTPC_Track()
        {
            if (fNumTPC_Tracks == fTPC_Tracks->GetSize())
                fTPC_Tracks->Expand( fNumTPC_Tracks + 10 );
            if (fNumTPC_Tracks >= 10000)
            {
                Fatal( "Ali_TRD_photons::createTPC_Track()", "ERROR: Too many TPC tracks (>10000)!" );
                exit( 2 );
            }

            new((*fTPC_Tracks)[fNumTPC_Tracks++]) Ali_TPC_Track;
            return (Ali_TPC_Track*)((*fTPC_Tracks)[fNumTPC_Tracks - 1]);
        }
        void clearTPC_TrackList()
        {
            fNumTPC_Tracks   = 0;
            fTPC_Tracks      ->Clear();
        }
        UShort_t getNumTPC_Tracks() const
        {
            return fNumTPC_Tracks;
        }
        Ali_TPC_Track* getTPC_Track(UShort_t i) const
        {
            return i < fNumTPC_Tracks ? (Ali_TPC_Track*)((*fTPC_Tracks)[i]) : NULL;
        }

        //-----------------------------------


        //-----------------------------------

        Ali_Kalman_Track* createKalman_Track()
        {
            if (fNumKalman_Tracks == fKalman_Tracks->GetSize())
                fKalman_Tracks->Expand( fNumKalman_Tracks + 10 );
            if (fNumKalman_Tracks >= 10000)
            {
                Fatal( "Ali_TRD_photons::createKalman_Track()", "ERROR: Too many Kalman tracks (>10000)!" );
                exit( 2 );
            }

            new((*fKalman_Tracks)[fNumKalman_Tracks++]) Ali_Kalman_Track;
            return (Ali_Kalman_Track*)((*fKalman_Tracks)[fNumKalman_Tracks - 1]);
        }
        void clearKalman_TrackList()
        {
            fNumKalman_Tracks   = 0;
            fKalman_Tracks      ->Clear();
        }
        UShort_t getNumKalman_Tracks() const
        {
            return fNumKalman_Tracks;
        }
        Ali_Kalman_Track* getKalman_Track(UShort_t i) const
        {
            return i < fNumKalman_Tracks ? (Ali_Kalman_Track*)((*fKalman_Tracks)[i]) : NULL;
        }


        ClassDef(Ali_TRD_Photon,1);  //
};
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
class Ali_TRD_Nuclear_interaction : public TObject
{
private:

    Float_t avg_sec_vertex[3];
    Float_t N_close_vertex;
    Float_t dcaAB_min;
    Float_t TOFsignal_min;
    Float_t Track_length_min;
    Float_t TPCdEdx_min;
    //Float_t dca_to_prim;
    Float_t pT_min;
    Float_t momentum_min;
    Float_t nuclev_bitmap;

    UShort_t fNumTPC_Tracks;
    UShort_t fNumKalman_Tracks;

    TClonesArray* fTPC_Tracks;              //->
    TClonesArray* fKalman_Tracks;           //->

public:
    Ali_TRD_Nuclear_interaction() :
        avg_sec_vertex(),N_close_vertex(-1),dcaAB_min(-1),TOFsignal_min(-1),Track_length_min(-1),TPCdEdx_min(-1),pT_min(-1),momentum_min(-1),nuclev_bitmap(-1),fNumTPC_Tracks(0),fNumKalman_Tracks(0)
    {
        fTPC_Tracks             = new TClonesArray( "Ali_TPC_Track", 10 );
        fKalman_Tracks          = new TClonesArray( "Ali_Kalman_Track", 10 );
    }
        ~Ali_TRD_Nuclear_interaction()
        {
            delete fTPC_Tracks;
            fTPC_Tracks = NULL;
            delete fKalman_Tracks;
            fKalman_Tracks = NULL;
        }

        // setters
        void set_avg_sec_vertex(Float_t x, Float_t y, Float_t z)         {    avg_sec_vertex[0] = x; avg_sec_vertex[1] = y; avg_sec_vertex[2] = z;  }
        void set_N_close_vertex(Float_t n)                              {    N_close_vertex = n;                        }
        void set_dcaAB_min(Float_t dca)                                  {    dcaAB_min = dca;                                                  }
        void set_TOFsignal_min(Float_t tof)                                      {    TOFsignal_min = tof;                                         }
        void set_Track_length_min(Float_t length)                                      {    Track_length_min = length;                                    }
        void set_TPCdEdx_min(Float_t dedx)                                  {    TPCdEdx_min = dedx;                                                  }
        //void set_dca_to_prim(Float_t dcaprim)                            {    dca_to_prim = dcaprim;                                            }
        void set_pT_min(Float_t pt)                                  {    pT_min = pt;                                                  }
        void set_momentum_min(Float_t mom)                                {    momentum_min = mom;                                                }
        void set_nuclev_bitmap(Float_t bitmap)                          {    nuclev_bitmap = bitmap;                                          }

        // getters

        Float_t get_avg_sec_vertex(Int_t i_xyz) const    {    return avg_sec_vertex[i_xyz];         }
        Float_t get_N_close_vertex() const       {    return N_close_vertex; }
        Float_t get_dcaAB_min() const                      {    return dcaAB_min;   }
        Float_t get_TOFsignal_min() const                        {    return TOFsignal_min; }
        Float_t get_Track_length_min() const                        {    return Track_length_min; }
        Float_t get_TPCdEdx_min() const                      {    return TPCdEdx_min; }
        //Float_t get_dca_to_prim() const                   {    return dca_to_prim; }
        Float_t get_pT_min() const                                 {    return pT_min; }
        Float_t get_momentum_min() const                               {    return momentum_min; }
        Float_t get_nuclev_bitmap() const                         {    return nuclev_bitmap; }

        //-----------------------------------

        Ali_TPC_Track* createTPC_Track()
        {
            if (fNumTPC_Tracks == fTPC_Tracks->GetSize())
                fTPC_Tracks->Expand( fNumTPC_Tracks + 10 );
            if (fNumTPC_Tracks >= 10000)
        {
        Fatal( "Ali_TRD_Nuclear_interaction::createTPC_Track()", "ERROR: Too many TPC tracks (>10000)!" );
        exit( 2 );
        }

        new((*fTPC_Tracks)[fNumTPC_Tracks++]) Ali_TPC_Track;
        return (Ali_TPC_Track*)((*fTPC_Tracks)[fNumTPC_Tracks - 1]);
    }
    void clearTPC_TrackList()
    {
        fNumTPC_Tracks   = 0;
        fTPC_Tracks      ->Clear();
    }
    UShort_t getNumTPC_Tracks() const
    {
        return fNumTPC_Tracks;
    }
    Ali_TPC_Track* getTPC_Track(UShort_t i) const
    {
        return i < fNumTPC_Tracks ? (Ali_TPC_Track*)((*fTPC_Tracks)[i]) : NULL;
    }
     
    //-----------------------------------


    //-----------------------------------   

    Ali_Kalman_Track* createKalman_Track()
    {
        if (fNumKalman_Tracks == fKalman_Tracks->GetSize())
        fKalman_Tracks->Expand( fNumKalman_Tracks + 10 );
        if (fNumKalman_Tracks >= 10000)
        {
        Fatal( "Ali_TRD_Nuclear_interaction::createKalman_Track()", "ERROR: Too many Kalman tracks (>10000)!" );
        exit( 2 );
        }

        new((*fKalman_Tracks)[fNumKalman_Tracks++]) Ali_Kalman_Track;
        return (Ali_Kalman_Track*)((*fKalman_Tracks)[fNumKalman_Tracks - 1]);
    }
    void clearKalman_TrackList()
    {
        fNumKalman_Tracks   = 0;
        fKalman_Tracks      ->Clear();
    }
    UShort_t getNumKalman_Tracks() const
    {
        return fNumKalman_Tracks;
    }
    Ali_Kalman_Track* getKalman_Track(UShort_t i) const
    {
        return i < fNumKalman_Tracks ? (Ali_Kalman_Track*)((*fKalman_Tracks)[i]) : NULL;
    }
    

    ClassDef(Ali_TRD_Nuclear_interaction,1);  //
};
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
class Ali_TRD_Self_Event : public TObject
{
private:
    Int_t eventNumber;
    Float_t x; // Event vertex x
    Float_t y; // Event vertex y
    Float_t z; // Event vertex z
    Int_t   id; // Run id
    Int_t   N_TPC_tracks; // total number of tracks
    Int_t   N_TRD_tracklets; // total number of TRD tracklets
    Float_t   cent_class_ZNA; // ZDC neutral A
    Float_t   cent_class_ZNC; // ZDC neutral C
    Float_t   cent_class_V0A; // V0 A
    Float_t   cent_class_V0C; // V0 C
    Float_t   cent_class_V0M; // V0 average
    Float_t   cent_class_CL0; // clusters in layer 0
    Float_t   cent_class_CL1; // clusters in layer 1
    Float_t   cent_class_SPD; // SPD
    Float_t   cent_class_V0MEq; //
    Float_t   cent_class_V0AEq; //
    Float_t   cent_class_V0CEq; //
    Float_t   ADC_sum_det[540];


    Float_t BeamIntAA; // ZDC coincidence rate
    Float_t T0zVertex; // z-vertex position from VPD

    TString TriggerWord; // Trigger word

    UShort_t      fNumTRD_Photons; // number of photon conversions
    UShort_t      fNumTRD_Nuclear_interactions; // number of nuclear interactions
    
    TClonesArray* fPhotons;                 //->
    TClonesArray* fNuclear_interactions;    //->
    //TClonesArray* fKalman_tracks;                 //->
    //TClonesArray* fTPC_tracks;    //->


public:
    Ali_TRD_Self_Event() :
    eventNumber(-1),x(-1),y(-1),z(-1),id(-1),N_TPC_tracks(0),N_TRD_tracklets(0),
    cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
        cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),BeamIntAA(-1),T0zVertex(-1),TriggerWord(),fNumTRD_Photons(0),fNumTRD_Nuclear_interactions(0),
        ADC_sum_det()
    {
        fPhotons                = new TClonesArray( "Ali_TRD_Photon", 10 );
        fNuclear_interactions   = new TClonesArray( "Ali_TRD_Nuclear_interaction", 10 );
        //fKalman_tracks                = new TClonesArray( "Ali_Kalman_track", 10 );
        //fTPC_tracks   = new TClonesArray( "Ali_TPC_track", 10 );
    }
    ~Ali_TRD_Self_Event()
    {
        delete fPhotons;
        fPhotons = NULL;
        delete fNuclear_interactions;
        fNuclear_interactions = NULL;
    }

    void       setEventNumber(Int_t e)            { eventNumber = e;                         }
    Int_t      getEventNumber() const             { return eventNumber;                      }

    void       setx(Float_t r)                    { x = r;                         }
    Float_t    getx() const                       { return x;                      }

    void       sety(Float_t r)                    { y = r;                         }
    Float_t    gety() const                       { return y;                      }

    void       setz(Float_t r)                    { z = r;                         }
    Float_t    getz() const                       { return z;                      }

    void       setid(Int_t  r)                    { id = r;                        }
    Int_t      getid() const                      { return id;                     }

    void       setN_TPC_tracks(Int_t r)                 { N_TPC_tracks = r;                    }
    Int_t      getN_TPC_tracks() const                    { return N_TPC_tracks;                 }

    void       setN_TRD_tracklets(Int_t r)                 { N_TRD_tracklets = r;                    }
    Int_t      getN_TRD_tracklets() const                    { return N_TRD_tracklets;                 }

    void       setcent_class_ZNA(Float_t r)             { cent_class_ZNA = r;                }
    Float_t      getcent_class_ZNA() const              { return cent_class_ZNA;             }

    void       setcent_class_ZNC(Float_t r)             { cent_class_ZNC = r;                }
    Float_t      getcent_class_ZNC() const              { return cent_class_ZNC;             }

    void       setcent_class_V0A(Float_t r)             { cent_class_V0A = r;                }
    Float_t      getcent_class_V0A() const              { return cent_class_V0A;             }

    void       setcent_class_V0C(Float_t r)             { cent_class_V0C = r;                }
    Float_t      getcent_class_V0C() const              { return cent_class_V0C;             }

    void       setcent_class_V0M(Float_t r)             { cent_class_V0M = r;                }
    Float_t      getcent_class_V0M() const              { return cent_class_V0M;             }

    void       setcent_class_CL0(Float_t r)             { cent_class_CL0 = r;                }
    Float_t      getcent_class_CL0() const              { return cent_class_CL0;             }

    void       setcent_class_CL1(Float_t r)             { cent_class_CL1 = r;                }
    Float_t      getcent_class_CL1() const              { return cent_class_CL1;             }

    void       setcent_class_SPD(Float_t r)             { cent_class_SPD = r;                }
    Float_t      getcent_class_SPD() const              { return cent_class_SPD;             }

    void       setcent_class_V0MEq(Float_t r)             { cent_class_V0MEq = r;                }
    Float_t      getcent_class_V0MEq() const              { return cent_class_V0MEq;             }

    void       setcent_class_V0AEq(Float_t r)             { cent_class_V0AEq = r;                }
    Float_t      getcent_class_V0AEq() const              { return cent_class_V0AEq;             }

    void       setcent_class_V0CEq(Float_t r)             { cent_class_V0CEq = r;                }
    Float_t      getcent_class_V0CEq() const              { return cent_class_V0CEq;             }

    void       setBeamIntAA(Float_t r)                 { BeamIntAA = r;                      }
    Float_t    getBeamIntAA() const                    { return BeamIntAA;                   }

    void       setT0zVertex(Float_t r)            { T0zVertex = r;                     }
    Float_t    getT0zVertex() const               { return T0zVertex;                  }

    void       setTriggerWord(TString s)          { TriggerWord = s;}
        TString    getTriggerWord() const             { return TriggerWord; }

        void       setADC_sum_det(Int_t i_det, Float_t r) { ADC_sum_det[i_det] = r;}
        ULong64_t  getADC_sum_det(Int_t i_det) const   { return ADC_sum_det[i_det]; }


        //----------------------------
    
    //----------------------------

 
        //----------------------------
        Ali_TRD_Photon* createPhoton() // photon conversion
    {
        if (fNumTRD_Photons == fPhotons->GetSize())
        fPhotons->Expand( fNumTRD_Photons + 10 );
        if (fNumTRD_Photons >= 10000)
        {
        Fatal( "Ali_TRD_Self_Event::createPhoton()", "ERROR: Too many photons (>10000)!" );
        exit( 2 );
        }

        new((*fPhotons)[fNumTRD_Photons++]) Ali_TRD_Photon;
        return (Ali_TRD_Photon*)((*fPhotons)[fNumTRD_Photons - 1]);
    }
    void clearPhotonList()
    {
        fNumTRD_Photons   = 0;
        fPhotons      ->Clear();
    }
    UShort_t getNumPhotons() const
    {
        return fNumTRD_Photons;
    }
    Ali_TRD_Photon* getPhoton(UShort_t i) const
    {
        return i < fNumTRD_Photons ? (Ali_TRD_Photon*)((*fPhotons)[i]) : NULL;
        }
        //----------------------------


  //----------------------------
        Ali_TRD_Nuclear_interaction* createNucInteraction() // photon conversion
    {
        if (fNumTRD_Nuclear_interactions == fNuclear_interactions->GetSize())
        fNuclear_interactions->Expand( fNumTRD_Nuclear_interactions + 10 );
        if (fNumTRD_Nuclear_interactions >= 10000)
        {
        Fatal( "Ali_TRD_Self_Event::createNuclear_interaction()", "ERROR: Too many photons (>10000)!" );
        exit( 2 );
        }

        new((*fNuclear_interactions)[fNumTRD_Nuclear_interactions++]) Ali_TRD_Nuclear_interaction;
        return (Ali_TRD_Nuclear_interaction*)((*fNuclear_interactions)[fNumTRD_Nuclear_interactions - 1]);
    }
    void clearNucInteractionsList()
    {
        fNumTRD_Nuclear_interactions   = 0;
        fNuclear_interactions      ->Clear();
    }
    UShort_t getNumNucInteractions() const
    {
        return fNumTRD_Nuclear_interactions;
    }
    Ali_TRD_Nuclear_interaction* getNucInteraction(UShort_t i) const
    {
        return i < fNumTRD_Nuclear_interactions ? (Ali_TRD_Nuclear_interaction*)((*fNuclear_interactions)[i]) : NULL;
        }
        //----------------------------
      


ClassDef(Ali_TRD_Self_Event,1);  // A simple event compiled of tracks
};
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
class Ali_Helix_copy : public TObject
{
private:
    Float_t        aliHelix_params[6];


public:
    Ali_Helix_copy() :
        aliHelix_params()
    {}
        ~Ali_Helix_copy(){}

        void setHelix(Float_t a, Float_t b,Float_t c,Float_t d,Float_t e,Float_t f)
        {
            aliHelix_params[0] = a; // y0
            aliHelix_params[1] = b; // z0
            aliHelix_params[2] = c;
            aliHelix_params[3] = d;
            aliHelix_params[4] = e; // c curvature
            aliHelix_params[5] = f; // x0
        }

        Float_t getHelix_param(Int_t i_param) const              {return aliHelix_params[i_param]; }

        void Evaluate(Double_t t, // helix evaluation, taken from AliHelix
                      Double_t r[3]);  //radius vector

        ClassDef(Ali_Helix_copy,1);
};
//----------------------------------------------------------------------------------------



//________________________________________________________________________
void Ali_Helix_copy::Evaluate(Double_t t,Double_t r[3])  //radius vector
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives at given phase
  //--------------------------------------------------------------------
  float phase=aliHelix_params[4]*t+aliHelix_params[2];
  Double_t sn=sinf(phase), cs=cosf(phase);
  //  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = aliHelix_params[5] + sn/aliHelix_params[4];
  r[1] = aliHelix_params[0] - cs/aliHelix_params[4];
  r[2] = aliHelix_params[1] + aliHelix_params[3]*t;
}
//________________________________________________________________________


#endif // __ALI_TRD_SELF_EVENT_H__
