#ifndef HOKO_H
#define HOKO_H
#include <stdio.h>
//#include <conio.h>
//#include <alloc.h>
#include <math.h>
#include <time.h>
//#include <io.h>
#include <string.h>
//#include <dos.h>
#include <stdlib.h>
#include <time.h>
#define far                      // It is used for Windows compatibility
///#define HOKODLL 1                // for DLL creation
#define DEMO false               // It is used for creation Demo version
//#define EPOCHA_IS_QUASINERTIAL   // It is used for coordinate system determination
#define DUT1  0.6/86400.0        //UT1-UTC    in days
#define DAT  32.0/86400.0        //TAI-UTC    in days
#define TDT_TAI  32.184/86400.0  //TDT-TAI    in days
#define BOKO  1                  //Translation for high altitude orbits
#define HOKO  1                  //Translation for low altitude orbits
#define GEM_10                   //Choose Earth gravity model
#if defined(GEM_10)              //
 #define LMAX 8                  //Maximum allowable value of gravity field order
 #define LGEM 300
#else
 #define LMAX 16                  //Maximum allowable value of gravity field order
 #define LGEM (LMAX+1)*(LMAX+2)-6 //Dimension of array for store gravity field coefficients
#endif
#define sqr(x) ((x)*(x))
#define onedeg(x) (((x)%2)? -1:1)
#define MIN(x,y) ((x)<(y)? (x):(y))
#define MAX(x,y) ((x)<(y)? (y):(x))
#define QMAX 6                   //Maximum allowable order of expansion on eccentricity
#define LLUN 5                   //Maximum order of terms considered in expansions of the Moon disturbing function
#define LSUN 3                   //Maximum order of  terms  in expansions of the Sun disturbing function
#define PI M_PI
#define MAXREZ LMAX              //Maximum value of considered resonance order (<=LMAX)
#define	ALUN 384400.0            //Semimajor axis of the Moon
#define	AC 149.6E6               //Semimajor axis of the Sun
#define	ML 0.012300123           //Ratio of the Moon mass to the Earth mass
#define	MC 332946.0              //Ratio of the Sun mass to the Earth mass
#if defined(EPOCHA_IS_QUASINERTIAL)
 #define CBZ 6.300387486753      //Rotational velocity of the Earth in the quasi-inertial coordinate system
#else
 #define CBZ 6.30038738          //Rotational velocity of the Earth in the Mean Equinox coordinate system
#endif
#define	 RADIAN 57.29577951
#define	 PI2 6.28318530719
#define	 ECJ 0.003352813178      //Flattening of the Earth
#define	 CTPOBFORCE 0.005        //Parameter, which determines how often will be recalculated slow changing functions of averaged equations
#define	 CTPOBKP 0.005           //Parameter, which determines how often will be recalculated slow changing functions of short periodics
#define	 EPSKP 1.0E-8            //Parameter, which determines the value of short periodical amplitudes for taken into account

extern double HR,skor,xp[3],nx;
extern long Nk,NO;
extern double HK1[LGEM],GME,RZ,C20;
//Structure of input orbital data
 struct bxprog{/*Name of structure*/
 long NO,NH;   /*Satellite number and revolution number*/
 double  DT,   /*Initial time counted off from 0hUTC Dec 31 1957 in days*/
	 A,        /*Semimajor axis, in km*/
	 I,        /*Inclination, in radians*/
	 DBY,      /*Longitude of ascending nodes, in radians*/
	 L,        /*Parameter of Laplace vector = e*cosw*/
	 H,        /*Parameter of Laplace vector = e*sinw*/
	 U,        /*Argument of latitude*/
	 KB,       /*Ballistic coefficient*/
	 F135,     /*The weighting-averaged value of solar activity index F10,7 for preceding 81 days F10,7*/
	 FT,       /*The daily averaged value of solar activity index F10,7*/
	 KP;       /*The daily averaged planetary index of geomagnetic activity*/
 int PRMODEL,  /*Parameter for determining motion model*/
/*
 if PRMODEL=0 then it is chosen full model: (l,m)<=(LMAX,LMAX), atmospheric drag, attraction of the Moon and Sun, solar radiation pressure
 if PRMODEL=1 then it is taken into account perturbations due second zonal harmonic
 if PRMODEL=n then it is chosen full model: (l,m)<=(n,n), atmospheric drag, attraction of the Moon and Sun, solar radiation pressure
*/
	 POCK,     /*Parameter for determining form of input and ouput data*/
/*
 if POCK=0 then input and output elements are osculating ones
 if POCK=-1 then input elements are mean elements, output elements are osculating ones
 if POCK=-2 then input and output elements are mean ones
*/

	 PZADACHI; /*Parameter for determining the solvable prediction task*/
/*
 if PZADACHI=0 then "prediction to requestion time";
 if PZADACHI=1 then "prediction to ascending node of preset revolution number";
 if PZADACHI=2 then "prediction to ascending node of current revolution";
 if PZADACHI=3 then satellite track calculation;
 if PZADACHI=4 then transformation of RSSS input elements to Keplerian format;
*/

	 };


struct bixprog
	{
struct bxprog bix;/*structure of output data*/
double  PER,     /*Nodal period, in min*/
	DT,          /*Change in nodal period due to the atmospheric drag (min/revolution)*/
	LIFET,       /*Reentry time estimation*/
	E,           /*eccentricity*/
	W,           /*argument of perigee, in radians*/
	HA,          /*Apogee altitude, in km*/
	HP,          /*Perigee altitude, in km*/
	RAD,         /*Geocentric range, in km*/
	ALTITUDE,    /*Latitude*/
	LONGITUDE,   /*Longitude*/
	X[6],        /*Vector of position and velocity in the ECI*/
	Y[6],        /*Keplerian mean elements*/
	Z[6];        /*Nonsingular mean elements*/
	};


#ifdef HOKODLL
 #define JVC_DLL  _export _stdcall
extern "C"
 {
#else
 #define JVC_DLL
#endif
      ; bool   JVC_DLL demo(int i)
      ; int    JVC_DLL sign(double X)
	  ; double JVC_DLL KEPLER(double M,double EKC)
	  ; void   JVC_DLL ALFDEL(double T,double far *ALF,double far *DEL)
	  ; double JVC_DLL GHCK(double far E0[],double far EK[])
	  ; double JVC_DLL XBL(double far X[],double far E[])
	  ; double JVC_DLL LBK(double far E[],double far E1[])
	  ; double JVC_DLL ST(int D,double T0)
	  ; void   JVC_DLL KOOPLS(double TM,double far *XL,double far *YL,double far *ZL,double far *XC,double far *YC,double far *ZC)
	  ; double JVC_DLL ATM(double H1,double F135,double F,int F0,double KP,double T0)
	  ; double JVC_DLL KBL(double far E[],double far E1[])
	  ; void   JVC_DLL UV(double X,double Y,double Z,int LK,double far HLC[])
	  ; void   JVC_DLL FINC(int L,int M,int P,double S,double C,double far*FF1,double  far*FF2)
      ; double JVC_DLL BIN(double B,int A,int K)
	  ; double JVC_DLL XNEWC1(int n,int m,int k,double E)
	  ; void   JVC_DLL HANSEN(int N,int P,int Q,double E,double  far*A,double  far*A1)
	  ; void   JVC_DLL C202(void)
	  ; void   JVC_DLL QPRT(void)
	  ; void   JVC_DLL LS(void)
	  ; void   JVC_DLL REZ(void)
      ; double JVC_DLL I0(double X)
      ; double JVC_DLL I1(double X)
	  ; void   JVC_DLL KATM(double H,double B1,double B2)
	  ; void   JVC_DLL APPROK(double DT,double T,double TK,double  far*THAX,double  far*TKOH)
	  ; int    JVC_DLL FORCE(double T,double  far E1[],double  far E2[])
	  ; int    JVC_DLL SKP(int PKP,double  far EBX[],double T,double  far E1[])
	  ; int    JVC_DLL PROGNOZ(int NTOCH,struct bxprog far*BX,struct bixprog far*BIX,double  far TKNK[])
	  ; int    JVC_DLL PROGNOZCHICL(int NTOCH,struct bxprog far*BX,struct bixprog far*BIX,double TKNK[])
//	  ; void   error_print(char *s)
	  ; void   JVC_DLL x_k(double h[],double e[], int pech)
	  ; void   JVC_DLL km_x(double E0[],int prz,double EK[],int pech)
	  ; int    JVC_DLL FORCE_CHIS(double *x,double *V,double tm,double *F)
	  ; int    JVC_DLL chisint(double *x,double *V,double tI,double tF,double xL,int LL,int NV,int NI,int *NF,int *NS,int klass,int NoR,int pr)
	  ; void   JVC_DLL track(int NT, double T)
      ; int    JVC_DLL F107KP(double t,double *f107,double *f81,double *kp)
	  ; void   JVC_DLL close_prognoz(void)
      ; void   JVC_DLL SHADOW(double al,double del,double AP,double DBY,double *E1/*âõîä*/,double *E2/*âûõîä*/,double *S,double *T,double *W)
// -¨þ¡õôº¨_ ¦--
      ; void   JVC_DLL error_print(char *s)
      ; void   JVC_DLL task_exit(int pech)
      ; int    JVC_DLL obmatr(double *az,int jaz,int n,double *a ,int ja )
      ; void   JVC_DLL vpm(double el[],double pme[7][7],double delvre)
      ; void   JVC_DLL D_rad_xyz(double radius,double alfa,double delta,double A[2][3])
      ;	void   JVC_DLL KOPCT(double f,double l,double H,double *L)
      ; void   JVC_DLL D_xyz_alambda(double e[6],double D[6][6])
      ;	void   JVC_DLL REDUCTION(double alfa,double delta,double JD1957,int Epocha, int Psr,
			                     double *d_alfa,double *d_delta)
      ; void   JVC_DLL KVAZI_TO_T(double X[],double T,int p)
      ; void   JVC_DLL T_TO_KVAZI(double X[],double T,int p)
      ; void   JVC_DLL KEP_TO_LAP(double A,double I,double w,double DBY,double T,
				double *IL,double *wL,double *DBYL,double *LL)
      ; void   JVC_DLL PROGNOZ_KMO(double EL[],double KMO[7][7],double DELT,double EK[],
			    double radius,double alfa,double delta,int P,double C[7][7],double C2[2][2])
      ; void   JVC_DLL GRAD_RAD1(char *gmsds,double *rdr,int k);
      ; void   JVC_DLL STR_TIME(double t,char *dmg,char *ch_m_s,int k)
      ;	void   JVC_DLL GNSK_V_GPSK(double E0[],double S,int prz,double EK[])
      ; double JVC_DLL CH_M_SToTime(char ch_m_s[])
      ; double JVC_DLL CH_M_SToDelta(char s_delta[])
      ; bool JVC_DLL GEM_input (void)
     ;
#ifdef HOKODLL
 }
#endif

#endif
