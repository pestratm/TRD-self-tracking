#ifndef __ALI_TRD_ST_H__
#define __ALI_TRD_ST_H__

#include "TObject.h"
#include "TClonesArray.h"

#include "TVector3.h"
#include "TLorentzVector.h"

//----------------------------------------------------------------------------------------
class Ali_MC_particle : public TObject
{
private:
    TVector3 TV3_particle_vertex; //
    TLorentzVector TLV_particle;
    Int_t PDGcode;
    Int_t index_particle;
    Int_t index_mother;
    Int_t N_daughters;
    Int_t arr_index_daughters[5];

public:
    Ali_MC_particle() :
        TV3_particle_vertex(),TLV_particle(),PDGcode(0),index_particle(0),index_mother(0),N_daughters(0),arr_index_daughters()
    {}
        ~Ali_MC_particle(){}

        void       set_TV3_particle_vertex(TVector3 TV3_particle_vertex_in)           { TV3_particle_vertex = TV3_particle_vertex_in;    }
        TVector3   get_TV3_particle_vertex() const                           { return TV3_particle_vertex;             }

        void set_TLV_particle(TLorentzVector tlv)     { TLV_particle = tlv; }
        TLorentzVector get_TLV_particle() const       { return TLV_particle;   }

        void       set_PDGcode(Int_t PDGcode_in)                    { PDGcode = PDGcode_in;          }
        Int_t      get_PDGcode() const                              { return PDGcode;                }

        void       set_index_mother(Int_t index_mother_in)                    { index_mother = index_mother_in;          }
        Int_t      get_index_mother() const                              { return index_mother;                }

        void       set_index_particle(Int_t index_particle_in)                    { index_particle = index_particle_in;          }
        Int_t      get_index_particle() const                              { return index_particle;                }

        void       set_N_daughters(Int_t N_daughters_in)                    { N_daughters = N_daughters_in;          }
        Int_t      get_N_daughters() const                              { return N_daughters;                }

        void       set_arr_index_daughters(Int_t index, Int_t index_daughter_in)  { arr_index_daughters[index] = index_daughter_in; }
        Int_t      get_arr_index_daughters(Int_t index) const                { return arr_index_daughters[index];      }


        ClassDef(Ali_MC_particle,1);
};
//----------------------------------------------------------------------------------------




//----------------------------------------------------------------------------------------
class Ali_TRD_ST_Tracklets : public TObject
{
private:
    TVector3 TV3_offset; //
    TVector3 TV3_dir;
    Int_t    TRD_det;
    Double_t ADC_val[30]; //
    UShort_t TPC_match;
    Int_t    n_tracklets_around;
    Double_t min_dist_to_next_trkl; //
    Int_t    index;


public:
    Ali_TRD_ST_Tracklets() :
        TV3_offset(),TV3_dir(),TRD_det(-4),ADC_val(),TPC_match(0),n_tracklets_around(0),min_dist_to_next_trkl(-5.0),index(0)
    {}
        ~Ali_TRD_ST_Tracklets(){}

        void       set_TV3_offset(TVector3 TV3_offset_in)           { TV3_offset = TV3_offset_in;    }
        TVector3   get_TV3_offset() const                           { return TV3_offset;             }

        void       set_TV3_dir(TVector3 TV3_dir_in)                 { TV3_dir = TV3_dir_in;          }
        TVector3   get_TV3_dir() const                              { return TV3_dir;                }

        void       set_TRD_det(Int_t TRD_det_in)                    { TRD_det = TRD_det_in;          }
        Int_t      get_TRD_det() const                              { return TRD_det;                }

        void       set_ADC_val(Int_t time_bin, Double_t ADC_value)  { ADC_val[time_bin] = ADC_value; }
        Double_t   get_ADC_val(Int_t time_bin) const                { return ADC_val[time_bin];      }


        void       set_TPC_match(UShort_t TPC_match_in)             { TPC_match = TPC_match_in;      }
        UShort_t   get_TPC_match() const      						{ return TPC_match;              }
                
        void       set_n_tracklets_around(Int_t n_tracklets_around_in)          { n_tracklets_around = n_tracklets_around_in;       }
        Int_t      get_n_tracklets_around() const                               { return n_tracklets_around;                        }

        void       set_min_dist_to_next_trkl(Double_t min_dist_to_next_trkl_in) { min_dist_to_next_trkl = min_dist_to_next_trkl_in; }
        Double_t   get_min_dist_to_next_trkl() const                            { return min_dist_to_next_trkl;                     }

        void       set_TRD_index(Int_t TRD_IND_in)                   { index = TRD_IND_in;          }
        Int_t      get_TRD_index() const                             { return index;                }


        ClassDef(Ali_TRD_ST_Tracklets,1);
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Ali_TRD_ST_TPC_Track : public TObject
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
    Float_t        TRD_signal; // TRD PID signal
    Float_t        TRDsumADC; // TRD electron PID probability
    Float_t        dca; // distance of closest approach of particle A
    TLorentzVector TLV_part; // Lorentz vector properties of this particle
    UShort_t       NTPCcls; // Number of TPC clusters
    UShort_t       NTRDcls; // Number of TRD clusters
    UShort_t       NITScls; // Number of TRD clusters
    UShort_t       status; // status of track: bit 0: ITS refit, bit1: TPC refit
    Float_t        TPCchi2; // TPC chi2
    Float_t        impact_angle_on_TRD; // Track impact angle on TRD
    Float_t        TPCdEdx; // Energy loss information of TPC
    Float_t        TOFsignal; // Time-of-flight
    Float_t        Track_length; // length of track

    Float_t        aliHelix_params[6];
    Float_t        aliHelix_TRD_params[6];

public:
    Ali_TRD_ST_TPC_Track() :
        nsigma_e_TPC(-3)
        ,nsigma_e_TOF(-3),nsigma_pi_TPC(-3),nsigma_pi_TOF(-3),nsigma_K_TPC(-3),nsigma_K_TOF(-3),nsigma_p_TPC(-3),nsigma_p_TOF(-3),TRD_signal(-3),
        TRDsumADC(-3),dca(-3),TLV_part(),NTPCcls(-3),NTRDcls(-3),NITScls(-3),status(-3),TPCchi2(-3),
        impact_angle_on_TRD(-3),TPCdEdx(-3),TOFsignal(-3),Track_length(-3),aliHelix_params(),aliHelix_TRD_params()
    {
    }
	~Ali_TRD_ST_TPC_Track()
	{
	}

	// setters
	void setnsigma_e_TPC(Float_t f)                     { nsigma_e_TPC = f;         }

        void setnsigma_e_TOF(Float_t f)                     { nsigma_e_TOF = f;         }
	void setnsigma_pi_TPC(Float_t f)                     { nsigma_pi_TPC = f;         }
	void setnsigma_pi_TOF(Float_t f)                     { nsigma_pi_TOF = f;         }
	void setnsigma_K_TPC(Float_t f)                     { nsigma_K_TPC = f;         }
	void setnsigma_K_TOF(Float_t f)                     { nsigma_K_TOF = f;         }
	void setnsigma_p_TPC(Float_t f)                     { nsigma_p_TPC = f;         }
	void setnsigma_p_TOF(Float_t f)                     { nsigma_p_TOF = f;         }
	void setTRDSignal(Float_t f)                     { TRD_signal = f;         }
	void setTRDsumADC(Float_t f)                     { TRDsumADC = f;         }
	void setdca(Float_t f)                    { dca = f;        }
	void set_TLV_part(TLorentzVector tlv)     { TLV_part = tlv; }
	void setNTPCcls(UShort_t s)               { NTPCcls = s;}
	void setNTRDcls(UShort_t s)               { NTRDcls = s;}
	void setNITScls(UShort_t s)               { NITScls = s;}
	void setStatus(UShort_t s)                { status = s;}
	void setTPCchi2(Float_t f)                { TPCchi2 = f;}
        void setimpact_angle_on_TRD(Float_t f)           {impact_angle_on_TRD = f;}
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
        void setHelix_TRD(Float_t a, Float_t b,Float_t c,Float_t d,Float_t e,Float_t f)
        {
            aliHelix_TRD_params[0] = a;
            aliHelix_TRD_params[1] = b;
            aliHelix_TRD_params[2] = c;
            aliHelix_TRD_params[3] = d;
            aliHelix_TRD_params[4] = e;
            aliHelix_TRD_params[5] = f;
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
	Float_t getTRDSignal() const                     { return TRD_signal;         }
	Float_t getTRDsumADC() const                     { return TRDsumADC;         }
	Float_t getdca() const                    { return dca;        }
	TLorentzVector get_TLV_part() const       { return TLV_part;   }
	UShort_t getNTPCcls() const               { return NTPCcls;    }
	UShort_t getNTRDcls() const               { return NTRDcls;    }
	UShort_t getNITScls() const               { return NITScls;    }
	UShort_t getStatus() const               { return status;    }
	Float_t  getTPCchi2() const              { return TPCchi2; }
        Float_t   getimpact_angle_on_TRD() const    { return impact_angle_on_TRD; }
	Float_t   getTPCdEdx() const                { return TPCdEdx; }
	Float_t   getTOFsignal() const              { return TOFsignal; }
        Float_t   getTrack_length() const           { return Track_length; }
        Float_t   getHelix_param(Int_t i_param) const              {return aliHelix_params[i_param]; }
        Float_t   getHelix_TRD_param(Int_t i_param) const              {return aliHelix_TRD_params[i_param]; }

        


        void Evaluate(Double_t t, // helix evaluation, taken from AliHelix
                     Double_t r[3]);  //radius vector
                     
        void Evaluate_TRD(Double_t t, Double_t r[3]);


        ClassDef(Ali_TRD_ST_TPC_Track,1);  // A simple track of a particle
};
//----------------------------------------------------------------------------------------



//________________________________________________________________________
void Ali_TRD_ST_TPC_Track::Evaluate(Double_t t,Double_t r[3])  //radius vector
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

//________________________________________________________________________
void Ali_TRD_ST_TPC_Track::Evaluate_TRD(Double_t t,Double_t r[3])  //radius vector
{
  //--------------------------------------------------------------------
  // Calculate position of a point on a track and some derivatives at given phase
  //--------------------------------------------------------------------
  float phase=aliHelix_TRD_params[4]*t+aliHelix_TRD_params[2];
  Double_t sn=sinf(phase), cs=cosf(phase);
  //  Double_t sn=TMath::Sin(phase), cs=TMath::Cos(phase);

  r[0] = aliHelix_TRD_params[5] + sn/aliHelix_TRD_params[4];
  r[1] = aliHelix_TRD_params[0] - cs/aliHelix_TRD_params[4];
  r[2] = aliHelix_TRD_params[1] + aliHelix_TRD_params[3]*t;
}
//________________________________________________________________________



//----------------------------------------------------------------------------------------
class Ali_TRD_ST_Event : public TObject
{
private:
    Float_t x; // Event vertex x
    Float_t y; // Event vertex y
    Float_t z; // Event vertex z
    Int_t   id; // Run id
    Int_t   N_tracks; // total number of tracks
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
    Short_t   flag_Data_MC; // 0 = data, 1 = MC


    Float_t BeamIntAA; // ZDC coincidence rate
    Float_t T0zVertex; // z-vertex position from VPD

    TString TriggerWord; // Trigger word
    UShort_t N_TRD_time_bins; // Number of TRD time bins

    UShort_t   fNumTracks; // number of tracks in event
    Int_t      fNumTracklets; // number of tracks in event
    Int_t      fNumMCparticles; // number of MC particles

    TClonesArray* fTracks;      //->
    TClonesArray* fTracklets;      //->
    TClonesArray* fMCparticles;  //->

public:
    Ali_TRD_ST_Event() :
        x(-3)
        ,y(-3),z(-3),id(-3),N_tracks(0),N_TRD_tracklets(0),
	cent_class_ZNA(0),cent_class_ZNC(0),cent_class_V0A(0),cent_class_V0C(0),cent_class_V0M(0),cent_class_CL0(0),cent_class_CL1(0),
        cent_class_SPD(0),cent_class_V0MEq(0),cent_class_V0AEq(0),cent_class_V0CEq(0),flag_Data_MC(0),BeamIntAA(-3),T0zVertex(-3),
        TriggerWord(),N_TRD_time_bins(0),fNumTracks(0),fNumTracklets(0),fNumMCparticles(0),fTracks(),fTracklets(),fMCparticles()
    {
        fTracks         = new TClonesArray( "Ali_TRD_ST_TPC_Track", 10 );
        fTracklets      = new TClonesArray( "Ali_TRD_ST_Tracklets", 10 );
        fMCparticles    = new TClonesArray( "Ali_MC_particle", 10 );
    }
	~Ali_TRD_ST_Event()
	{
	    delete fTracks;
            fTracks = NULL;
            delete fTracklets;
            fTracklets = NULL;
            delete fMCparticles;
            fMCparticles = NULL;
	}

	void       setx(Float_t r)                    { x = r;                         }
	Float_t    getx() const                       { return x;                      }


	void       sety(Float_t r)                    { y = r;                         }
	Float_t    gety() const                       { return y;                      }

	void       setz(Float_t r)                    { z = r;                         }
	Float_t    getz() const                       { return z;                      }

	void       setid(Int_t  r)                    { id = r;                        }
	Int_t      getid() const                      { return id;                     }

	void       setN_tracks(Int_t r)                 { N_tracks = r;                    }
	Int_t      getN_tracks() const                    { return N_tracks;                 }

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

        void       set_flag_Data_MC(UShort_t s)             { flag_Data_MC = s;                }
	UShort_t      get_flag_Data_MC() const              { return flag_Data_MC;             }

	void       setBeamIntAA(Float_t r)                 { BeamIntAA = r;                      }
	Float_t    getBeamIntAA() const                    { return BeamIntAA;                   }

	void       setT0zVertex(Float_t r)            { T0zVertex = r;                     }
	Float_t    getT0zVertex() const               { return T0zVertex;                  }

	void       setTriggerWord(TString s)          { TriggerWord = s;}
        TString    getTriggerWord() const             { return TriggerWord; }

        void       setN_TRD_time_bins(UShort_t s)             { N_TRD_time_bins = s;                }
	Float_t    getN_TRD_time_bins() const              { return  N_TRD_time_bins;            }


        //----------------------------
	Ali_TRD_ST_TPC_Track* createTrack()
	{
	    if (fNumTracks == fTracks->GetSize())
		fTracks->Expand( fNumTracks + 10 );
	    if (fNumTracks >= 65000)
	    {
		Fatal( "Ali_TRD_ST_Event::createTrack()", "ERROR: Too many tracks (>65000)!" );
		exit( 2 );
	    }

	    new((*fTracks)[fNumTracks++]) Ali_TRD_ST_TPC_Track;
	    return (Ali_TRD_ST_TPC_Track*)((*fTracks)[fNumTracks - 1]);
	}
	void clearTrackList()
	{
	    fNumTracks   = 0;
	    fTracks      ->Clear();
	}
	UShort_t getNumTracks() const
	{
	    return fNumTracks;
	}
	Ali_TRD_ST_TPC_Track* getTrack(UShort_t i) const
	{
	    return i < fNumTracks ? (Ali_TRD_ST_TPC_Track*)((*fTracks)[i]) : NULL;
        }
        //----------------------------


        //----------------------------
        Ali_TRD_ST_Tracklets* createTracklet() // online tracklet
	{
	    if (fNumTracklets == fTracklets->GetSize())
		fTracklets->Expand( fNumTracklets + 10 );
	    if (fNumTracklets >= 650000)
	    {
		Fatal( "Ali_TRD_ST_Event::createTracklet()", "ERROR: Too many tracklets (>650000)!" );
		exit( 2 );
	    }

	    new((*fTracklets)[fNumTracklets++]) Ali_TRD_ST_Tracklets;
	    return (Ali_TRD_ST_Tracklets*)((*fTracklets)[fNumTracklets - 1]);
	}
	void clearTrackletList()
	{
	    fNumTracklets   = 0;
	    fTracklets      ->Clear();
	}
	Int_t getNumTracklets() const
	{
	    return fNumTracklets;
	}
	Ali_TRD_ST_Tracklets* getTracklet(Int_t i) const
	{
	    return i < fNumTracklets ? (Ali_TRD_ST_Tracklets*)((*fTracklets)[i]) : NULL;
        }
        //----------------------------



        //----------------------------
        Ali_MC_particle* createMCparticle() // Monte Carlo particle
	{
            if (fNumMCparticles == fMCparticles->GetSize())
		fMCparticles->Expand( fNumMCparticles + 10 );
	    if (fNumMCparticles >= 650000)
	    {
		Fatal( "Ali_TRD_ST_Event::createMCparticle()", "ERROR: Too many MC particles (>650000)!" );
		exit( 2 );
	    }

	    new((*fMCparticles)[fNumMCparticles++]) Ali_MC_particle;
	    return (Ali_MC_particle*)((*fMCparticles)[fNumMCparticles - 1]);
	}
	void clearMCparticleList()
	{
	    fNumMCparticles   = 0;
	    fMCparticles      ->Clear();
	}
	Int_t getMCparticles() const
	{
	    return fNumMCparticles;
	}
	Ali_MC_particle* getMCparticle(Int_t i) const
	{
	    return i < fNumMCparticles ? (Ali_MC_particle*)((*fMCparticles)[i]) : NULL;
        }
        //----------------------------

ClassDef(Ali_TRD_ST_Event,1);  // A simple event compiled of tracks
};
//----------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------
class Ali_Helix : public TObject
{
private:
    Float_t        aliHelix_params[6];


public:
    Ali_Helix() :
        aliHelix_params()
    {}
        ~Ali_Helix(){}

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

        ClassDef(Ali_Helix,1);
};
//----------------------------------------------------------------------------------------



//________________________________________________________________________
void Ali_Helix::Evaluate(Double_t t,Double_t r[3])  //radius vector
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




#endif // __ALI_TRD_ST_H__
