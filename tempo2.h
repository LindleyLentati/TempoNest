//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

#ifndef __Tempo2_h
#include <stdio.h>
#include <time.h>
#define __Tempo2_h
#define TEMPO2_h_VER "$Revision: 1.65 $"
#define TSUN (4.925490947e-6L) // Solar constant for mass calculations.
#define MAX_FREQ_DERIVATIVES 13    /* F0 -> Fn   where n=10                            */
#define MAX_DM_DERIVATIVES   10    /* DM0 -> DMn where n=10                            */
#define MAX_PSR_VAL          40    /* Maximum number of pulsars                        */             
#define MAX_COMPANIONS       4     /* Maximum number of binary companions              */
#define NE_SW_DEFAULT        4     /* Default value for electron density (cm-3) at 1AU due to solar wind */
#define ECLIPTIC_OBLIQUITY_VAL 84381.4059 /* mean obliquity of ecliptic in arcsec        */
#define MAX_COEFF            5000  /* Maximum number of coefficients in polyco         */
#define MAX_CLKCORR          5000  /* Maximum number of lines in time.dat file         */
#define MAX_LEAPSEC          100   /* Maximum number of line in the leap second file   */
#define MAX_STRLEN           1000  /* Maximum length for strings                       */
#define MAX_FILELEN          500   /* Maximum filename length                          */
#define MAX_STOREPRECISION   50    /* How many routines in TEMPO2 store precision information */
#define MAX_OBSN_VAL         10000 /* Maximum number of TOAs                           */ 
#define MAX_SITE             100   /* Maximum number of observatory sites              */
#define MAX_PARAMS           500   /* Maximum number of parameters                     */
#define MAX_JUMPS            2000  /* Maximum number of phase jumps                    */
#define MAX_WHITE            100   /* Maximum number of parameters for whitening       */
#define MAX_IFUNC            500   /* Maximum number of parameters for interpolation function  */
#define MAX_TEL_CLK_OFFS     500   /* Maximum number of parameters for telescope clock offset */
#define MAX_TEL_DX           500   /* Maximum number of parameters for interpolation function  */
#define MAX_TEL_DY           500   /* Maximum number of parameters for interpolation function  */
#define MAX_TEL_DZ           500   /* Maximum number of parameters for interpolation function  */
#define MAX_FIT              10000  /* Maximum number of parameters to fit for */
#define MAX_T2EFAC           50    /* Maximum number of T2EFACs allowed                */
#define MAX_T2EQUAD          50    /* Maximum number of T2EQUADs allowed               */
#define MAX_BPJ_JUMPS        5     /* Maximum number of jumps in binary params - for BPJ model */
#define MAX_TOFFSET          10    /* Number of time jumps allowed in .par file        */
#define MAX_QUAD             150   /* Maximum number of frequency channels in quadrupolar function */
#define MAX_DMX              64    /* Max number of DM steps allowed */
#define MAX_FLAGS            20    /* Maximum number of flags in .tim file/observation */
#define MAX_FLAG_LEN         32    /* Maximum number of characters in each flag */
#define MAX_CLK_CORR         30    /* Maximum number of steps in the correction to TT  */ 
#define SECDAY               86400.0       /* Number of seconds in 1 day                 */
#define SPEED_LIGHT          299792458.0 /* Speed of light (m/s)                       */
#define SOLAR_MASS  1.98892e30           /* Mass of Sun (kg)                           */
#define SOLAR_RADIUS 6.96e8              /* Radius of the Sun (in meters)              */
#define BIG_G       6.673e-11            /* Gravitational constant                     */
#define GM          1.3271243999e20      /* Gravitational constant * mass sun          */
#define GM_C3       4.925490947e-6       /* GM_odot/c^3 (in seconds)                   */
#define GMJ_C3      4.70255e-9           /* GM_jupiter/c^3 (in seconds)                */
#define GMS_C3      1.40797e-9           /* GM_saturn/c^3 (in seconds)                 */
#define GMV_C3      1.2061e-11           /* GM_venus/c^3 (in seconds)                  */
#define GMU_C3      2.14539e-10          /* GM_uranus/c^3 (in seconds)                 */
#define GMN_C3      2.54488e-10          /* GM_neptune/c^3 (in seconds)                */
#define OBLQ        23.4458333333333333  /* Obliquity of the ecliptic                  */
#define AULTSC      499.00478364         /* Number of light seconds in 1 AU            */
#define AU_DIST     1.49598e11           /* 1 AU in m                                  */
#define DM_CONST    2.41e-4
#define DM_CONST_SI 7.436e6              /* Dispersion constant in SI units            */
#define PCM         3.08568025e16        /* one parsec in meters                       */
#define MASYR2RADS  1.53628185e-16       /* Converts from mas/yr to rad/s              */
#define MAX_MSG     50                   /* Maximum number of different warnings       */

/* Path for the file containing dates when leap seconds should be added */
#define LEAPSECOND_FILE "/clock/leap.sec"

/* Path for the file containing TAI-UT1 */
#define UT1_FILE "/clock/ut1.dat" 

/* Path for file containing IERS EOP C04 series */
#define EOPC04_FILE "earth/eopc04_IAU2000.62-now"

/* Path for file containing TDB-TDT ephemeris */
#define TDBTDT_FILE "/ephemeris/TDB.1950.2050"
#define IFTEPH_FILE "/ephemeris/TIMEEPH_short.te405"

/* Path for file containing Observatory data (obsys.dat) */
#define OBSSYS_FILE "/observatory/newobsys.dat"

/* TEMPO emulation modes */
#define SI_UNITS 1  /* New tempo2 mode */
#define TDB_UNITS 2 /* original tempo mode */

#define IF99_TIMEEPH 1  /* Irwin & Fukushima time ephemeris */
#define FB90_TIMEEPH 2  /* Fairhead & Bretagnon time ephemeris */

#define T2C_IAU2000B 1
#define T2C_TEMPO   2


/* define some functions for log message 
 * M.Keith 2012 - let me know if this fails to compile anywhere.
 * mkeith@pulsarastronomy.net
 **/
#define LOG_OUTFILE stdout
#define WHERESTR  "[%s:%d] "
#define WHEREARG  __FILE__, __LINE__
#define ENDL "\n"
#define WHEREERR "******\nERROR [%s:%d] "
#define WHERETCHK "[%s:%d] T=%.2f s: "
#define _LOG(...) fprintf(LOG_OUTFILE,__VA_ARGS__)
#define logmsg(_fmt, ...) _LOG(WHERESTR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logdbg(_fmt, ...)  if(debugFlag)logmsg(_fmt,##__VA_ARGS__)
#define logerr(_fmt, ...) _LOG(WHEREERR _fmt ENDL, WHEREARG,##__VA_ARGS__)
#define logtchk(_fmt, ...) if(tcheck)_LOG(WHERETCHK _fmt ENDL, WHEREARG,(clock()-timer_clk)/(float)CLOCKS_PER_SEC,##__VA_ARGS__)



/* Type for doing extra precision computations: longdouble */

/* OSes/architectures that have built-in long double */
//#if defined sun || defined linux || (defined __APPLE__ && defined __GNUC__ && __GNUC__ >= 4)
#ifndef TEMPO2_USE_QD
#define USE_BUILTIN_LONGDOUBLE
#endif

#ifdef sun
#include <sunmath.h> 
//* there is no such file ! J. Wang */
// Note: you sometimes need a compile line e.g. -I/opt/SUNWspro/WS6U2/include/cc
#endif

#ifdef USE_BUILTIN_LONGDOUBLE
typedef long double longdouble;
#define LONGDOUBLE_ONE 1.0L
/* OSes/architectures lacking built-in double; use "qd" library */
#else
#include "dd.h"
typedef dd_real longdouble;
#define LONGDOUBLE_ONE "1.0"
dd_real pow(const dd_real &a, const dd_real &b);
//operator float(const dd_real &a);
#endif

/* function to get longdouble as string (%g style). This returns
  std::string, from which you can get a normal C char * like this:
  print_longdouble(x).c_str(). Unfortunately we can't just return char *
  directly as it would end up pointing to de-allocated memory, and we
  can't do it statically in case you want to say call this twice for
  one printf */
#ifdef __cplusplus
#include <string>
std::string print_longdouble(const longdouble &ld);
#endif

/* function to parse a string as longdouble */
longdouble parse_longdouble(const char *str);

/* TEMPO2 environment variable */
extern char TEMPO2_ENVIRON[];

/* TEMPO2 error messages */
extern char TEMPO2_ERROR[];

enum label {param_raj,param_decj,param_f,param_pepoch,param_posepoch,
	    param_dmepoch,param_dm,param_pmra,param_pmdec,param_px,
	    param_sini,param_pb,param_fb,param_t0,param_a1,param_om,param_pmrv,
	    param_ecc,param_edot,param_xpbdot,param_pbdot,param_a1dot,
	    param_omdot,param_tasc,param_eps1,param_eps2,param_m2,param_gamma,
	    param_mtot,param_glep,param_glph,param_glf0,param_glf1,param_glf2,
	    param_glf0d,param_gltd,param_start,param_finish,param_track,param_bp,param_bpp,
	    param_tzrmjd,param_tzrfrq,param_fddc,param_fddi,param_dr,param_dtheta,param_tspan,
            param_bpjep,param_bpjph,param_bpja1,param_bpjec,param_bpjom,param_bpjpb,
            param_wave_om,param_kom,param_kin,param_shapmax,param_dth,param_a0,
	    param_b0,param_xomdot,param_afac,param_eps1dot,param_eps2dot,param_tres,
            param_dshk,param_ephver,param_daop,param_iperharm,param_dmassplanet,param_waveepoch,param_ifunc,param_clk_offs,
            param_dmx,param_dmxr1,param_dmxr2,param_dmmodel,param_gwsingle,param_cgw,param_quad_om,
	    param_h3,param_h4,param_nharm,param_stig,
            param_telx,param_tely,param_telz,param_telEpoch,param_quad_ifunc_p,
	    param_quad_ifunc_c,param_tel_dx,param_tel_dy,param_tel_dz,
	    param_tel_vx,param_tel_vy,param_tel_vz,param_tel_x0,param_tel_y0,param_tel_z0,param_gwm_amp};

/*
 * These represent the possible constraints to the fit that have been implemented.
 */
enum constraint {
	constraint_dmmodel_mean,
	constraint_dmmodel_cw_0,
	constraint_dmmodel_cw_1,
	constraint_dmmodel_cw_2,
	constraint_dmmodel_cw_3,
	constraint_ifunc_0,
	constraint_ifunc_1,
	constraint_ifunc_2,
	constraint_tel_dx_0,
	constraint_tel_dx_1,
	constraint_tel_dx_2,
	constraint_tel_dy_0,
	constraint_tel_dy_1,
	constraint_tel_dy_2,
	constraint_tel_dz_0,
	constraint_tel_dz_1,
	constraint_tel_dz_2,
	constraint_quad_ifunc_p_0,
	constraint_quad_ifunc_p_1,
	constraint_quad_ifunc_p_2,
	constraint_quad_ifunc_c_0,
	constraint_quad_ifunc_c_1,
	constraint_quad_ifunc_c_2,
	constraint_dmmodel_cw_year_sin,
	constraint_dmmodel_cw_year_cos,
	constraint_dmmodel_cw_year_xsin,
	constraint_dmmodel_cw_year_xcos,
	constraint_dmmodel_cw_year_sin2,
	constraint_dmmodel_cw_year_cos2,
	constraint_dmmodel_cw_px,
	constraint_ifunc_year_sin,
	constraint_ifunc_year_cos,
	constraint_ifunc_year_xsin,
	constraint_ifunc_year_xcos,
	constraint_ifunc_year_sin2,
	constraint_ifunc_year_cos2,
	constraint_qifunc_p_year_sin,
	constraint_qifunc_p_year_cos,
	constraint_qifunc_p_year_xsin,
	constraint_qifunc_p_year_xcos,
	constraint_qifunc_p_year_sin2,
	constraint_qifunc_p_year_cos2,
	constraint_qifunc_c_year_sin,
	constraint_qifunc_c_year_cos,
	constraint_qifunc_c_year_xsin,
	constraint_qifunc_c_year_xcos,
	constraint_qifunc_c_year_sin2,
	constraint_qifunc_c_year_cos2,
};


extern int MAX_PSR;
extern int MAX_OBSN;
extern double ECLIPTIC_OBLIQUITY;

extern int debugFlag;   /* Global = 1 if debug mode is running */
extern int forceGlobalFit;   /* Global = 1 if we are forcing a global fit */
extern int writeResiduals;   /* Global = 1 if we are writing out post-fit residuals */
extern int tcheck;   /* Global = 1 if time check message should be printed is running */
extern clock_t timer_clk;
extern int veryFast;    /* Global to run the code fast */
extern char tempo2MachineType[MAX_FILELEN];
extern int displayCVSversion; /* Display CVS version */

extern char dcmFile[MAX_FILELEN];
extern char covarFuncFile[MAX_FILELEN];;
extern char tempo2_plug_path[32][MAX_STRLEN];
extern int tempo2_plug_path_len;

typedef struct storePrecision {
  longdouble minPrec;
  char routine[100];
  char comment[MAX_STRLEN];
} storePrecision;

/* If this structure is modified - must update copyParam in tempo2Util.C */
typedef struct parameter {
  char **label;              /* Label about this parameter                         */
  char **shortlabel;         /* Label about this parameter without units           */
  longdouble *val;           /* Value of parameter                                 */
  longdouble *err;           /* Uncertainty on parameter value                     */
  int  *fitFlag;             /* = 1 if fitting required, = 2 for global fit        */
  int  *paramSet;            /* = 1 if parameter has been set                      */
  longdouble *prefit;        /* Pre-fit value of the parameter                     */
  longdouble *prefitErr;     /* Pre-fit value of the uncertainty                   */
  int aSize;                 /* Number of elements in the array for this parameter */
  int linkFrom[5];
  int linkTo[5];
  int nLinkTo;
  int nLinkFrom;
} parameter;

/* struct clock_correction : struct observation contains an array
   of these, which getClockCorrections() fills in */
typedef struct 
{
  double correction;
  char corrects_to[32];
} clock_correction;

typedef struct observation {
  longdouble sat;                 /* Site arrival time                                          */
  longdouble bat;                 /* Infinite frequency barycentric arrival time                */
  longdouble bbat;                /* Arrival time at binary barycentre                          */
  longdouble pet;                 /* Pulsar emission time                                       */
  int clockCorr;                  /* = 1 for clock corrections to be applied, = 0 for BAT       */
  int delayCorr;                  /* = 1 for time delay corrections to be applied, = 0 for BAT  */
  int deleted;                    /* = 1 if observation has been deleted                        */
                                  /* = -1 if not included in fit                                */
  longdouble prefitResidual;      /* Pre-fit residual                                           */
  longdouble residual;            /* residual                                                   */
  double      freq;               /* Frequency of observation (in MHz)                          */
  double      freqSSB;            /* Frequency of observation in barycentric frame (in Hz)      */
  double      toaErr;             /* Error on TOA (in us)                                       */
  double      toaDMErr;           /* Error on TOA due to DM (in us)                             */
  double      origErr;            /* Original error on TOA after reading tim file (in us)       */
  double      phaseOffset;        /* Phase offset                                               */
  char        fname[MAX_FILELEN]; /* Name of data file giving TOA                               */
  char        telID[100];         /* Telescope ID                                               */
  clock_correction correctionsTT[MAX_CLK_CORR]; /* chain of corrections from site TOA to chosen realisation of TT */
  int nclock_correction;

  longdouble correctionTT_TB;     /* Correction to TDB/TCB           */
  double einsteinRate;            /* Derivative of correctionTT_TB   */
  longdouble correctionTT_Teph;   /* Correction to Teph              */
  longdouble correctionUT1;       /* Correction from site TOA to UT1 */
  /*  long double a1utcf; */

  double sun_ssb[6];              /* Ephemeris values for Sun w.r.t SSB (sec)             (RCS) */
  double sun_earth[6];            /* Ephemeris values for Sun w.r.t Earth (sec)                 */
  double planet_ssb[9][6];        /* Ephemeris values for all planets w.r.t. SSB (sec)   */
  double jupiter_earth[6];        /* Ephemeris values for Jupiter w.r.t. Earth centre (sec)     */
  double saturn_earth[6];         /* Ephemeris values for Saturn w.r.t. Earth centre (sec)      */
  double venus_earth[6];          /* Ephemeris values for Venus w.r.t. Earth centre (sec)      */
  double uranus_earth[6];         /* Ephemeris values for Uranus w.r.t. Earth centre (sec)      */
  double neptune_earth[6];        /* Ephemeris values for Neptune w.r.t. Earth centre (sec)     */
  double earthMoonBary_ssb[6];    /* Ephem values for Earth-Moon barycentre wrt SSB (sec) (RCB) */
  double earthMoonBary_earth[6];  /* Position of Earth-Moon barycentre with respect to Earth (sec) (RBE) */
  double earth_ssb[6];            /* Centre of Earth w.r.t. SSB                                 */
  double observatory_earth[6];    /* Observatory site with respect to Earth centre (sec)  (REA) */  
  double psrPos[3];               /* Unit vector giving position of the pulsar at observation time from Earth */
  double zenith[3];               /* Zenith vector, in BC frame. Length=geodetic height */
  double nutations[6];
  double siteVel[3];              /* Observatory velocity w.r.t. geocentre                      */

  long double shklovskii;         /* Shklovskii delay term                                      */
  double shapiroDelaySun;         /* Shapiro Delay due to the Sun                               */
  double shapiroDelayJupiter;     /* Shapiro Delay due to Jupiter                               */
  double shapiroDelaySaturn;      /* Shapiro Delay due to Saturn                                */
  double shapiroDelayVenus;       /* Shapiro Delay due to Venus                                 */
  double shapiroDelayUranus;      /* Shapiro Delay due to Uranus                                */
  double shapiroDelayNeptune;     /* Shapiro Delay due to Neptune                               */
  double troposphericDelay;       /* Delay due to neutral refraction in atmosphere              */
  double tdis1;                   /* Interstellar dispersion measure delay                      */
  double tdis2;                   /* Dispersion measure delay due to solar system               */
  longdouble roemer;              /* Roemer delay                                               */
  longdouble torb;                /* Combined binary delays */
  longdouble nphase;              /* allows the pulse number to be determined                   */
  longdouble phase;               
  unsigned long long pulseN;                    /* Pulse number */

  char flagID[MAX_FLAGS][MAX_FLAG_LEN];     /* Flags in .tim file                                         */
  char flagVal[MAX_FLAGS][MAX_FLAG_LEN];
  int  nFlags;                   
  int  jump[MAX_FLAGS];           /* Jump region */
  int  obsNjump;                  /* Number of jumps for this observation */
  double efac;                    /* Error multiplication factor                                */
  double equad;                   /* Value to add in quadrature                                 */
} observation;

typedef struct pulsar {
  char  name[100];
  /*                                                                 */
  /* Note: when adding a new parameter, initialise it in intialise.c */
  /*                                                                 */
  int  fixedFormat;              /* = 0 for separate .par and .tim files, > 0 indicates number of lines to skip */
  parameter param[MAX_PARAMS];
  char rajStrPre[100],decjStrPre[100];   /* String containing RAJ and DECJ  (prefit)              */
  char rajStrPost[100],decjStrPost[100]; /* String containing RAJ and DECJ  (postfit)           */
  char binaryModel[100];                 /* Binary model e.g. BT/ELL1/BT2P etc.                        */


  double **ToAextraCovar;

  int    dmoffsDMnum;
  int    dmoffsCMnum;
  double dmoffsDM_mjd[MAX_IFUNC]; 
  double dmoffsDM[MAX_IFUNC];
  double dmoffsDM_error[MAX_IFUNC];
  double dmoffsDM_weight[MAX_IFUNC];
  double dmoffsCM_mjd[MAX_IFUNC]; 
  double dmoffsCM[MAX_IFUNC];
  double dmoffsCM_error[MAX_IFUNC];
  double dmoffsCM_weight[MAX_IFUNC];

  // Single source gravitational wave information
  double gwsrc_ra;
  double gwsrc_dec;
  double gwsrc_aplus_r,gwsrc_aplus_i;
  double gwsrc_across_r,gwsrc_across_i;
  double gwsrc_aplus_r_e,gwsrc_aplus_i_e;
  double gwsrc_across_r_e,gwsrc_across_i_e;
  double gwsrc_epoch;
  double gwsrc_psrdist;

  // Extra parameters for Xingjiang's updates to the continuous wave sources
  double cgw_h0;
  double cgw_cosinc;
  double cgw_angpol;
  double cgw_mc;

  // Gravitational wave memory information
  double gwm_raj;
  double gwm_decj;
  double gwm_epoch;
  double gwm_phi; // Polarisation angle

  // General pulsar information
  double posPulsar[3];            /* 3-vector pointing at pulsar                                */
  double velPulsar[3];            /* 3-vector giving pulsar's velocity                          */  
  long double phaseJump[MAX_JUMPS];    /* Time of phase jump                                         */
  int    phaseJumpDir[MAX_JUMPS]; /* Size and direction of phase jump                           */
  int    phaseJumpID[MAX_JUMPS];  /* ID of closest point to the phase jump */
  int    nPhaseJump;              /* Number of phase jumps                                      */
  double dmOffset;                /* Value to add to DM flags */
  double ne_sw;                   /* Electron density at 1AU due to the solar wind              */
  int    nCompanion;              /* Number of binary companions                                */
  int    eclCoord;                /* = 1 for ecliptic coords otherwise celestial coords         */

  int    nJumps;                  /* Number of jumps                                            */
  char fjumpID[16];
  double jumpVal[MAX_JUMPS];      /* Value of jump                                              */
  int    fitJump[MAX_JUMPS];      /* = 1 if fit for jump                                        */
  double jumpValErr[MAX_JUMPS];   /* Error on jump                                              */
  char   jumpStr[MAX_JUMPS][MAX_STRLEN]; /* String describing jump                              */
  char   filterStr[MAX_STRLEN];   /* String describing filters */
  char   passStr[MAX_STRLEN];   /* String describing filters */
  double tOffset[MAX_TOFFSET];    /* Offsets in TOAs in seconds                                 */ 
  double tOffset_f1[MAX_TOFFSET],tOffset_f2[MAX_TOFFSET];  /* Range for offset to be applied    */
  double tOffset_t1[MAX_TOFFSET],tOffset_t2[MAX_TOFFSET];
  char   tOffsetSite[MAX_TOFFSET][100],tOffsetFlags[MAX_TOFFSET][1000];
  int    nToffset;
  int    ndmx;                    /* Number of DM steps */
  double fitChisq;                /* Chisq value from the fit */
  int    fitNfree;                /* Number of degrees of freedom in fit */
  int    globalNfit;              /* Total number of parameters in the fit */
  int    globalNoConstrain;       /* Total number of points without constraints */
  int    nFit;                    /* Number of points in the fit */
  int    nParam;                  /* Number of parameters in the fit */
  int    nGlobal;                 /* Number of global parameters in the fit */
  int    fitParamI[MAX_FIT];
  int    fitParamK[MAX_FIT];
  int    fitMode;                 /* = 0 not fitting with errors, = 1 fitting with errors (MODE 1) */
  int    rescaleErrChisq;         /* = 1 to rescale errors based on the reduced chisq, = 0 not to do this */
  double offset;                  /* Offset, always fitted for */
  double offset_e;                /* Error in the offset */
  double **covar; //[MAX_PARAMS][MAX_PARAMS];

  int    calcShapiro;              /* = 1 Calculate Solar system Shapiro delay (otherwise -1)*/
  int    planetShapiro;            /* = 1 if included otherwise 0 */
  int    jboFormat;                /* = 1 => JBO arrival time format and file structure (not byte swapping) = 2 => JBO format with byte swapping */

  observation *obsn; /* [MAX_OBSN_VAL]; */
  int nobs;                       /* Number of observations in .tim file                        */
  int units;  /* TDB or SI units (tempo emulation mode uses TDB) 
                                     see #define definition above for possible units            */
  int tempo1; /* = 1 if tempo1 is emulated */
  int dilateFreq;  /* whether or not to apply SS time dilation to RFs */
  int timeEphemeris;              /* Which code to use for Einstein delay */
  int t2cMethod;  /* How to transform from terrestrial to celestial coords */
  int correctTroposphere;     /* whether or not do correct for tropospheric delay */
  int noWarnings;                 /* = 1, do not display warning messages                       */
  /* Path for the file containing the corrections between observatory clocks and UTC(NIST)      */
  /* - set in readParfile.C                                                                     */
/*   char OBSERVATORY_CLOCK_2_UTC_NIST[MAX_FILELEN];  */
  char clock[16]; /* Clock standard to use as "UTC" */
  char clockFromOverride[64];    /* Clock code to assume TOAs are measured against (e.g. UTC to turn off clock corrections, or TDB/TCG to turn off those + Einstein delay */
  char JPL_EPHEMERIS[MAX_FILELEN];
  char ephemeris[MAX_FILELEN];
  storePrecision storePrec[MAX_STOREPRECISION];
  int  nStorePrecision;
  int  bootStrap;           /* > 0 if calculating errors using bootstrap Monte-Carlo method */
  char tzrsite[100];        /* Site-code for polyco                                         */
  double rmsPre,rmsPost;
  char deleteFileName[100]; /* File name containing deleted points                          */
  int  nits;                /* Number of iterations for the fit                             */
  int  ipm;                 /* = 1 if use interplanetary medium DM correction, = 0 otherwise*/
  int  swm;                 /* = 0 for basic tempo2 solar wind model, = 1 for XPY Solar wind model */

  /* For whitening */
  double wave_sine[MAX_WHITE], wave_sine_err[MAX_WHITE];
  double wave_cos[MAX_WHITE],  wave_cos_err[MAX_WHITE];
  int    nWhite;
  double waveScale;

  // Quadrapolar function
  double quad_aplus_r[MAX_QUAD],quad_aplus_r_e[MAX_QUAD];
  double quad_aplus_i[MAX_QUAD],quad_aplus_i_e[MAX_QUAD];
  double quad_across_r[MAX_QUAD],quad_across_r_e[MAX_QUAD];
  double quad_across_i[MAX_QUAD],quad_across_i_e[MAX_QUAD];
  double quadEpoch;
  double quadRA,quadDEC;
  int    nQuad;

  // IFUNC

  double ifuncT[MAX_IFUNC], ifuncV[MAX_IFUNC],ifuncE[MAX_IFUNC],ifunc_weights[MAX_IFUNC];
  int    ifuncN;

  double clk_offsT[MAX_TEL_CLK_OFFS], clk_offsV[MAX_TEL_CLK_OFFS];
  double clk_offsE[MAX_TEL_CLK_OFFS];
  int    clkOffsN;

  // Quadrupolar IFUNC
  
  double quad_ifuncT_p[MAX_IFUNC], quad_ifuncV_p[MAX_IFUNC],quad_ifuncE_p[MAX_IFUNC];
  int    quad_ifuncN_p;
  double quad_ifuncT_c[MAX_IFUNC], quad_ifuncV_c[MAX_IFUNC],quad_ifuncE_c[MAX_IFUNC];
  int    quad_ifuncN_c;
  double quad_ifunc_p_RA,quad_ifunc_p_DEC,quad_ifunc_c_RA,quad_ifunc_c_DEC;
  double quad_ifunc_geom_p,quad_ifunc_geom_c;

  // Telescope position error
  int    nTelDX,setTelVelX;
  double telDX_t[MAX_TEL_DX];
  double telDX_v[MAX_TEL_DX];
  double telDX_e[MAX_TEL_DX];
  double telDX_vel[MAX_TEL_DX];
  double telDX_vel_e[MAX_TEL_DX];

  int    nTelDY,setTelVelY;
  double telDY_t[MAX_TEL_DY];
  double telDY_v[MAX_TEL_DY];
  double telDY_e[MAX_TEL_DY];
  double telDY_vel[MAX_TEL_DY];
  double telDY_vel_e[MAX_TEL_DY];

  int    nTelDZ,setTelVelZ;
  double telDZ_v[MAX_TEL_DZ];
  double telDZ_t[MAX_TEL_DZ];
  double telDZ_e[MAX_TEL_DZ];
  double telDZ_vel[MAX_TEL_DZ];
  double telDZ_vel_e[MAX_TEL_DZ];

  // T2EFAC/T2EQUAD
  int    nT2efac,nT2equad;
  char   T2efacFlagID[MAX_T2EFAC][MAX_FLAG_LEN],T2efacFlagVal[MAX_T2EFAC][MAX_FLAG_LEN];
  double T2efacVal[MAX_T2EFAC];
  char   T2equadFlagID[MAX_T2EQUAD][MAX_FLAG_LEN],T2equadFlagVal[MAX_T2EQUAD][MAX_FLAG_LEN];
  double T2equadVal[MAX_T2EQUAD];
  double T2globalEfac;
  

  /* Which fit function are we using */
  char fitFunc[MAX_FILELEN];


  int nconstraints;                       /* Number of fit constraints specified                      */
  enum constraint constraints[MAX_PARAMS];/* Which constraints are specified */
  char auto_constraints;
} pulsar;

// PLUGIN FUNCTIONS


/* FUNCTION DEFINITIONS */
#ifdef __cplusplus
extern "C" {
#endif
int id_residual(float xcurs,float ycurs);
float setStart(float xcurs,float ycurs,int flag);
int zoom_graphics(float xcurs2,float ycurs2,int flag);
void getInputs(pulsar *psr,int argc, char *argv[],char timFile[][MAX_FILELEN],
	       char parFile[][MAX_FILELEN],int *displayParams,int *npsr,
	       int *nGlobal,int *outRes,int *writeModel,char *outputSO,int *polyco,
	       char *polyco_args, int *newpar,int *onlypre,char *dcmFile,char *covarFuncFile,char* newparname);
void polyco(pulsar *psr,int npsr,longdouble polyco_MJD1,longdouble polyco_MJD2,int nspan,int ncoeff,
	    longdouble maxha,char *sitename,longdouble freq,longdouble coeff[MAX_COEFF],int trueDM);
void readParfile(pulsar *psr,char parFile[][MAX_FILELEN],char timFile[][MAX_FILELEN],int npsr);
void readParfileGlobal(pulsar *psr,int npsr,char tpar[MAX_STRLEN][MAX_FILELEN],
		       char ttim[MAX_STRLEN][MAX_FILELEN]);
int readSimpleParfile(FILE *fin,pulsar *p);
int setupParameterFileDefaults(pulsar *p);
void displayParameters(int pos,char timeFile[][MAX_FILELEN],char parFile[][MAX_FILELEN],pulsar *psr,int npsr);

void initialise (pulsar *psr, int noWarnings);

void initialiseOne (pulsar *psr, int noWarnings, int fullSetup);
void destroyOne (pulsar *psr);

void recordPrecision(pulsar *psr,longdouble prec,char *routine,char *comment);
void readTimfile(pulsar *psr,char timFile[][MAX_FILELEN],int npsr);
void formBats(pulsar *psr,int npsr);
void formBatsAll(pulsar *psr,int npsr);
void updateBatsAll(pulsar *psr,int npsr);
void formResiduals(pulsar *psr,int npsr,int removeMean);
int  bootstrap(pulsar *psr,int p,int npsr);
void doFitAll(pulsar *psr,int npsr,char *covarFuncFile);
void doFit(pulsar *psr,int npsr,int writeModel);
void doFitDCM(pulsar *psr,char *dcmFile,char *covarFuncFile,int npsr,int writeModel);
void doFitGlobal(pulsar *psr,int npsr,double *globalParameter,int nGlobal,int writeModel); 
void getCholeskyMatrix(double **uinv, char* fname, pulsar *psr, double *resx,double *resy,double *rese, int np, int nc, int* ip);
double getParamDeriv(pulsar *psr,int ipos,double x,int i,int k);
void textOutput(pulsar *psr,int npsr,double globalParameter,int nGlobal,int outRes,int newpar,char *fname);
void shapiro_delay(pulsar *psr,int npsr,int p,int i,double delt,double dt_SSB);
void dm_delays(pulsar *psr,int npsr,int p,int i,double delt,double dt_SSB);
void calculate_bclt(pulsar *psr,int npsr);
void secularMotion(pulsar *psr,int npsr);

void autoConstraints(pulsar* psr, int ipsr,int npsr);
void setPlugPath();

void preProcess(pulsar *psr,int npsr,int argc,char *argv[]);

/* used by psrchive to initialize a single psr struct */
void preProcessSimple (pulsar *psr);

/* three stages of preProcess that interleave with preProcess (complicated) */

void preProcessSimple1 (pulsar *psr, int tempo1, double thelast);

void preProcessSimple2 (pulsar *psr,
			float startdmmjd,
			int ndm,
			float* dmvals,
			int trimonly);

void preProcessSimple3 (pulsar *psr);



void useSelectFile(char *fname,pulsar *psr,int npsr);
void processSimultaneous(char *line,pulsar *psr, int npsr);
void processFlag(char *line,pulsar *psr,int npsr);
void toa2utc(pulsar *psr,int npsr);
void utc2tai(pulsar *psr,int npsr);
void tt2tb(pulsar *psr,int npsr);
void tai2tt(pulsar *psr,int npsr);
void tai2ut1(pulsar *psr,int npsr);
void vectorPulsar(pulsar *psr,int npsr);
void readEphemeris(pulsar *psr,int npsr,int addEphemNoise);
void get_obsCoord(pulsar *psr,int npsr);
double calcRMS(pulsar *psr,int p);

void allocateMemory(pulsar *psr,int realloc);
void destroyMemory(pulsar *psr);

void readJBO_bat(char *fname,pulsar *psr,int p);
void readObsFile(double alat[MAX_SITE],double along[MAX_SITE],
		 double elev[MAX_SITE],int icoord[MAX_SITE],
		 char obsnam[MAX_SITE][100],char obscode[MAX_SITE][100],
		 int *nobservatory,int obsnum[MAX_SITE]);
double dotproduct(double *v1,double *v2);
void vectorsum(double *res, double *v1,double *v2);
void vectorscale(double *v, double k);
void writeTim(char *timname,pulsar *psr,char *fileFormat);
int turn_hms(double turn, char *hms);
int turn_dms(double turn, char *dms);
double dms_turn(char *line);
double hms_turn(char *line);
double turn_deg(double turn);
longdouble fortran_mod(longdouble a,longdouble p);
int fortran_nint(double x);
void equ2ecl(double *x);
void copyParam(parameter p1,parameter *p2);
void copyPSR(pulsar *p,int p1,int p2);
longdouble getParameterValue(pulsar *psr,int param,int arr);
void simplePlot(pulsar *psr, double unitFlag);
float solarWindModel(pulsar psr,int iobs);

/* BINARY MODELS */
double MSSmodel(pulsar *psr,int p,int obs,int param);
void updateMSS(pulsar *psr,double val,double err,int pos);
double BTmodel(pulsar *psr,int p,int obs,int param);
void updateBT(pulsar *psr,double val,double err,int pos);
double BTJmodel(pulsar *psr,int p,int obs,int param,int arr);
void updateBTJ(pulsar *psr,double val,double err,int pos,int arr);
double BTXmodel(pulsar *psr,int p,int obs,int param,int arr);
void updateBTX(pulsar *psr,double val,double err,int pos,int arr);
double ELL1model(pulsar *psr,int p,int obs,int param);
void updateELL1(pulsar *psr,double val,double err,int pos);
long double DDmodel(pulsar *psr,int p,int obs,int param);
void updateDD(pulsar *psr,double val,double err,int pos);
double T2model(pulsar *psr,int p,int obs,int param,int arr);
void updateT2(pulsar *psr,double val,double err,int pos,int arr);
double JVmodel(pulsar *psr,int p,int obs,int param,int arr);
void updateJV(pulsar *psr,double val,double err,int pos,int arr);
double DDKmodel(pulsar *psr,int p,int obs,int param);
void updateDDK(pulsar *psr,double val,double err,int pos);
double DDSmodel(pulsar *psr,int p,int obs,int param);
void updateDDS(pulsar *psr,double val,double err,int pos);
double DDGRmodel(pulsar *psr,int p,int obs,int param);
void updateDDGR(pulsar *psr,double val,double err,int pos);
// Freire & Wex (2010; FW10) models:
double DDHmodel( pulsar *psr, int p, int obs, int param );
void updateDDH( pulsar *psr, double val, double err, int pos );
double ELL1Hmodel( pulsar *psr, int p, int obs, int param );
void updateELL1H( pulsar *psr, double val, double err, int pos );

void displayMsg(int type,char *key,char *searchStr,char *variableStr,int noWarnings);
void CVSdisplayVersion(char *file,char *func,const char *verNum);

/* stuff for SI/TDB units */
void transform_units(struct pulsar *psr, int from, int to);

/* This function uses the numerical recipes svdfit for the fitting */
void FITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos,int ipsr);
void updateParameters(pulsar *psr,int p,double *val,double *error);

/* defineClockCorrectionSequence: call to provide the clock correction
   module with a sequence of files to use for corrections. May be called
   multiple times for sequences with different start/end clocks (e.g. for
   multi-observatory fitting). */
void defineClockCorrectionSequence(char *fileList,int dispWarnings);

/* getClockCorrections : gets the sequence of corrections for a particular
   observation and stores them in obs->clock_corrections.  Uses one
   of the pre-defined sequences (from defineClockCorrectionSequence) if
   available, otherwise makes one automatically.  */
void getClockCorrections(observation *obs, char *clockFrom, char *clockTo,int warnings);

/* getCorrectionTT : convenience function to return the sum of all
   correctionsTT terms in an observation */
double getCorrectionTT(observation *obs);
/* convenience function to obtain correction to a named clock 
(for intermediate use e.g. in obtaining Earth orientation parameters;
  does not store steps used in obs->correctionsTT */
double getCorrection(observation *obs, char *clockFrom, char *clockTo, int warnings);

/* redwards stuff for tempo2 to look after the database of observatories */

typedef struct
{
  double x, y, z; /* Geocentric coordinates */
  double longitude_grs80, latitude_grs80, height_grs80; /* GRS80 geodetic coords */
  char name[32];
  char code[16];
  char clock_name[16];
} observatory;

observatory *getObservatory(char *code);
void lookup_observatory_alias(char *incode, char *outcode);
void 
get_obsCoord_IAU2000B(double observatory_trs[3],
		      double zenith_trs[3],
		      longdouble tt_mjd, longdouble utc_mjd,
		      double observatory_crs[3],
		      double zenith_crs[3],
		      double observatory_velocity_crs[3]);


/* redwards stuff to get earth orientation parameters */
void get_EOP(double mjd, double *xp, double *yp, double *dut1, 
	     double *dut1dot, int dispWarnings);
/* ... and tropospheric delays ... */
void compute_tropospheric_delays(pulsar *psr,int npsr);


#ifdef __cplusplus
}
#endif

#endif /* Defined __Tempo2_h */




