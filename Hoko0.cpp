#include "hoko.h"
#if BOKO
#define		ALP 0.228027134951881    /*=13.0649924465/RADIAN*/
#define		A1P 0.0172019697686106   /*=0.985600267/RADIAN*/
#define		F1P 0.230895723247662    /*=13.229350449/RADIAN*/
#define		DP  0.212768711686213    /*=12.1907491914/RADIAN*/
#define		B3P 0.229971502953377    /*=13.1763965268/RADIAN*/
#define		ALSP 0.0172027912671643  /*=0.9856473354/RADIAN*/
#define		r5 4.84813681109762e-6    /*=1/206264.806247*/


extern int	   LM
		,  LT
		,  LKL
		,  LKS
		,  JREZ[MAXREZ]
		 ;
extern double

		   A
		,  E
		,  I
		,  SI
		,  CI
		,  SI2
		,  CI2
		,  CE
		,  SE
		,  EA
		,  IA
		,  CIA
		,  SIA
		,  DT0
		;


       double    far*HL
	      , far*DL
	      , far*D2L
	      , far*HL0
	      , far*HL1
	      , far*HL2
	      , far*HL3
	      , far*HL4
	      , far*DC
	      , far*D2C
	      , far*HC
	      , far*HC0
	      , far*HC1
	      , far*HC2
	      , far*HC3
	      , far*HC4
	      , far*LE   /* LE[LLUN-1][LLUN+1][LLUN+1]*/
	      , far*LI   /* LI[LLUN-1][LLUN+1][LLUN+1]*/
	      , far*LV   /* LV[LLUN-1][LLUN+1][LLUN+1]*/
	      , far*LPI  /* LPI[LLUN-1][LLUN+1][LLUN+1]*/
	      , far*LLAM /*LLAM[LLUN-1][LLUN+1][LLUN+1]*/
	      ,	far* QC
	      , far* QD
	      , far* PC
	      , far* PD
	      , far* RC
	      , far* RD
	      , far* TC
	      , far* TD

	      ;

 double

			 f1
		      ,  B3
		      ,  B3F1
		      ,  PS1
		      ,  B2
		      ,  SB2
		      ,  CB2
		      ,  AL
		      ,  A1
		      ,  D
		      ,  SL1
		      ,  CL1
		      ,  SL
		      ,  CL
		      ,  S2L
		      ,  C2L
		      ,  S3L
		      ,  C3L
		      ,  SF
		      ,  CF
		      ,  S2F
		      ,  C2F
		      ,  SD
		      ,  CD
		      ,  S2D
		      ,  C2D
		      ,  S4D
		      ,  C4D
		      ,  SLM2D
		      ,  CLM2D
		      ,  SLP2D
		      ,  CLP2D
		      ,  SLMD
		      ,  CLMD
		      ,  S2LMD
		      ,  C2LMD
		      ,  SPLM2D
		      ,  CPLM2D
		      ,  S1M2D
		      ,  C1M2D
		      ,  S1P2D
		      ,  C1P2D
		      ,  SLM1
		      ,  CLM1
		      ,  SLP1
		      ,  CLP1
		      ,  SLM4D
		      ,  CLM4D
		      ,  S2LM4D
		      ,  C2LM4D
		      ,  S2LPD
		      ,  C2LPD
		      ,  S1PD
		      ,  S3LM2D
		      ,  S2LM1
		      ,  SLM1M2
		      ,  CLM1M2
		      ,  SLM1P2
		      ,  CLM1P2
		      ,  SLM2F
		      ,  CLM2F
		      ,  SLP2F
		      ,  S2FMD
		      ,  SLM2FD
		      ,  SLPF
		      ,  SLMF
		      ,  SFM2D
		      ,  CFM2D
		      ,  SFP2D
		      ,  SF2DML
		      ,  SF2DPL
		      ,  S2LPF
		      ,  SFM2L
		      ,  SFML2D
		      ,  SFPL2D
		      ,  SFP12D
		      ,  S2LDF
		      ,  S2DL1F
		      ,  SFM2D1
		      ,  SFM12D
		      ,  dl
		      ,  S1
		      ,  O
		      ,  B4
		      ,  SB4
		      ,  CB4
		      ,  SS1
		      ,  CS1
		      ,  CSSB
		      ,  B
		      ,  S2
		      ,  C2
		      ,  S3
		      ,  C3
		      ,  S4
		      ,  C4
		      ,  DLS
		      ,  DR4
		      ,  R4
		      ,  R4D
		      ,  SB
		      ,  CB
		      ,  BET

		      ,  Xl
		      ,  Yl
		      ,  Zl
		      ,  Xc
		      ,  Yc
		      ,  Zc
		      , kappa
		      , omega
		      , nu
		      ;

/*########################## Function KOOPLS   ###########################*/
 /* ======================================================================
	 Назначение: Определение координат Луны и Солнца на заданный момент
							 времени.
	 Форма обращения: KOOPLS(TM,XL,YL,ZL,XC,YC,ZC)
	 Входная информация:
        ТM - время, отсчитываемое от 0h 31 декабря 1957 г. (сутки),
	 Выходная информация:
	 Координаты Луны и Солнца: XL,YL,ZL,XC,YC,ZC в ГНСК (км.)
     ------------------------------------------------------------------- */
/* ======================================================================
	 Destination: Calculation of the Moon and Sun coordinates at
							 requestion time.
	 Call: KOOPLS(TM,XL,YL,ZL,XC,YC,ZC).
	 Input data: ТM is the time in days from 0h UTC Dec 31, 1957.
	 Output data:
	 Coordinates of the Moon (L) and Sun (C) in km: XL,YL,ZL,XC,YC,ZC.
     ------------------------------------------------------------------- */
	     void JVC_DLL KOOPLS(double TM,double far *XL,double far *YL,double far *ZL,double far *XC,double far *YC,double far *ZC)
	       {

	       static double
		   AL0
		,  A10
		,  F10
		,  D0
		,  DT00
		,  B3O
		,  ALSO
		,  P1
		,  PS
		,  AM1
		,  AM2
		,  AM3
		,  AM4
		,  ANU
		,  AE1
		,  AE2
		,  AE3
		,  AE4
		,  sinB2
		,  cosB2
		,  DI;
		double  P[8],d_mu,d_nu,d_epsilon;

			if  (fabs(DT00-DT0)>1.0e-9)
		   {
			    double   B0
				  ,  T1
				  ,  T2
				  ,  T3
				  ,  EP
				  ,  ES
				  ,  AE
				  ,  EZ
				  ,  ED
				  ,  EV
			    ;  DT00=DT0
			    ;  B0=DT0+21183.5+TDT_TAI+DAT
			    ;  T1=B0/36525.0
			    ;  T2=T1*T1
			    ;  T3=T1*T2
			    ;  EP=(23.452294-0.0130125*T1-0.164*T2*1.0E-5+0.503*T3*1.0E-6)/RADIAN
			    ;  ES=0.01675104-0.0000418*T1-0.000000126*T2
			    ;  AL0=(296.104608+0.009192*T2+0.0000147*T3)/RADIAN
			    ;  A10=(358.475845-0.000150*T2-0.000003*T3)/RADIAN
			    ;  F10=(11.250889-0.003214*T2-0.3E-6*T3)/RADIAN
			    ;  D0=(350.737486-0.001436*T2+0.19E-5*T3)/RADIAN
			    ;  B3O=(270.434164-0.001133*T2+0.19E-5*T3)/RADIAN
			    ;  ALSO=(279.696678+0.302E-3*T2)/RADIAN
			    ;  AL0+=ALP*B0
			    ;  F10+=F1P*B0
			    ;  A10+=A1P*B0
			    ;  D0+=DP*B0
			    ;  B3O+=B3P*B0
			    ;  ALSO+=ALSP*B0
			    ;  PS=(17.2327+0.01737*T1)*r5
			    ;  DI=(9.21+0.00091*T1)*r5
			    ;  P1=20.46*AC/RZ*r5
			    ;  sinB2=sin(EP)
			    ;  cosB2=cos(EP)
			    ;  EZ=ES*ES
			    ;  ED=EZ*ES
			    ;  EV=ED*ES
			    ;  AM1=2*ES-0.25*ED
			    ;  AM2=1.25*EZ-11.0/24.0*EV
			    ;  AM3=13.0/12.0*ED
			    ;  AM4=103.0/96.0*EV
			    ;  AE=1.00000023*AC/RZ
			    ;  ANU=AE*(1+EZ*0.5)
			    ;  AE1=AE*(ES-0.375*ED)
			    ;  AE2=AE*(EZ*0.5-EV/3)
			    ;  AE3=AE*ED*0.375
			    ;  AE4=AE*EV/3  ;
		  }
		;  f1=F10+F1P*TM
		;  B3=B3O+B3P*TM
		;  B3F1=B3-f1
		;  PS1=PS*sin(B3F1)
		;  B2=DI*cos(B3F1)
		;  SB2=sinB2+B2*cosB2
		;  CB2=cosB2-B2*sinB2
		;  AL=AL0+ALP*TM
		;  A1=A10+A1P*TM
		;  D=D0+DP*TM
		;  SL1=sin(A1)
		;  CL1=cos(A1)
		;  SL=sin(AL)
		;  CL=cos(AL)
		;  S2L=2*CL*SL
		;  C2L=CL*CL-SL*SL
		;  S3L=S2L*CL+C2L*SL
		;  C3L=C2L*CL-S2L*SL
		;  SF=sin(f1)
		;  CF=cos(f1)
		;  S2F=2*CF*SF
		;  C2F=CF*CF-SF*SF
		;  SD=sin(D)
		;  CD=cos(D)
		;  S2D=2*SD*CD
		;  C2D=CD*CD-SD*SD
		;  S4D=2*S2D*C2D
		;  C4D=C2D*C2D-S2D*S2D
		;  SLM2D=SL*C2D-S2D*CL
		;  CLM2D=CL*C2D+SL*S2D
		;  SLP2D=SL*C2D+S2D*CL
		;  CLP2D=CL*C2D-SL*S2D
		;  SLMD=SL*CD-SD*CL
		;  CLMD=CL*CD+SL*SD
		;  S2LMD=2*SLMD*CLMD
		;  C2LMD=CLMD*CLMD-SLMD*SLMD
		;  SPLM2D=SL1*CLM2D+CL1*SLM2D
		;  CPLM2D=CL1*CLM2D-SL1*SLM2D
		;  S1M2D=SL1*C2D-S2D*CL1
		;  C1M2D=CL1*C2D+SL1*S2D
		;  S1P2D=SL1*C2D+CL1*S2D
		;  C1P2D=CL1*C2D-SL1*S2D
		;  SLM1=SL*CL1-CL*SL1
		;  CLM1=CL*CL1+SL*SL1
		;  SLP1=SL*CL1+CL*SL1
		;  CLP1=CL*CL1-SL*SL1
		;  SLM4D=SL*C4D-CL*S4D
		;  CLM4D=CL*C4D+SL*S4D
		;  S2LM4D=S2L*C4D-C2L*S4D
		;  C2LM4D=C2L*C4D+S2L*S4D
		;  S2LPD=S2L*C2D+C2L*S2D
		;  C2LPD=C2L*C2D-S2L*S2D
		;  S1PD=SL1*CD+CL1*SD
		;  S3LM2D=S3L*C2D-C3L*S2D
		;  S2LM1=S2L*CL1-C2L*SL1
		;  SLM1M2=SLM2D*CL1-CLM2D*SL1
		;  CLM1M2=CLM2D*CL1+SLM2D*SL1
		;  SLM1P2=SLP2D*CL1-CLP2D*SL1
		;  CLM1P2=CLP2D*CL1+SLP2D*SL1
		;  SLM2F=SL*C2F-CL*S2F
		;  CLM2F=CL*C2F+SL*S2F
		;  SLP2F=SL*C2F+CL*S2F
		;  S2FMD=S2F*C2D-C2F*C2D
		;  SLM2FD=SLM2D*C2F-CLM2D*S2F
		;  SLPF=SL*CF+CL*SF
		;  SLMF=SL*CF-CL*SF
		;  SFM2D=SF*C2D-CF*S2D
		;  CFM2D=CF*C2D+SF*S2D
		;  SFP2D=SF*C2D+CF*S2D
		;  SF2DML=SF*CLM2D-CF*SLM2D
		;  SF2DPL=SF*CLM2D+CF*SLM2D
		;  S2LPF=S2L*CF+C2L*SF
		;  SFM2L=SF*C2L-CF*S2L
		;  SFML2D=SF*CLP2D-CF*SLP2D
		;  SFPL2D=SF*CLP2D+CF*SLP2D
		;  SFP12D=SL1*CFM2D+CL1*SFM2D
		;  S2LDF=S2LMD*CF+C2LMD*SF
		;  S2DL1F=SF*CPLM2D-CF*SPLM2D
		;  SFM2D1=SF*C1P2D-CF*S1P2D
		;  SFM12D=SF*C1M2D-CF*S1M2D
		;  dl=(22639.58*SL-4586.438*SLM2D+2369.899*S2D+769.021*S2L-668.944*SL1-411.614*S2F-211.658*S2LMD-206.219*
       SPLM2D+191.954*SLP2D-165.351*S1M2D+147.878*SLM1-124.785*SD-109.804*SLP1-55.174*S2FMD-45.1*SLP2F-   /*PROVERIT ZNAK - */
       39.532*SLM2F-38.428*SLM4D+36.134*S3L-30.773*S2LM4D+28.511*SLM1M2-24.451*S1P2D+18.609*SLMD+18.023*
       S1PD+14.577*SLM1P2+14.387*S2LPD+13.902*S4D-13.193*S3LM2D+9.703*S2LM1+9.366*SLM2FD)*r5
		;  S1=(18461.48*SF+1010.18*SLPF+999.695*SLMF-623.658*SFM2D+199.485*SF2DML-166.577*SF2DPL+117.262*SFP2D+
       61.913                  *S2LPF-33.359*SFML2D-31.763*SFM2L-29.689*SFP12D-15.565*S2LDF+15.122*SFPL2D+8.902*S2DL1F+12.14*
       SFM2D1                  +8.001*SFM12D)*r5
		;  O=(3422.7+186.5398*CL+34.3117*CLM2D+28.2333*C2D+10.1657*C2L+3.0861*CLP2D+1.9202*C1M2D+1.4455*CPLM2D+
       1.1542                  *CLM1-0.9752*CD-0.9502*CLP1-0.7136*CLM2F+0.6215*C3L+0.6008*CLM4D-0.3997*CL1+0.3722*C2LM4D-
       0.3039                  *C2LMD-0.3*C1P2D+0.2833*C2LPD+0.2607*C4D-0.2257*CLM1M2+0.2302*CLM1P2)*r5
		;  B4=B3+dl-PS1
		;  SB4=sin(B4)
		;  CB4=cos(B4)
		;  SS1=sin(S1)
		;  CS1=cos(S1)
		;  A1=1/O
		;  CSSB=CS1*SB4
		;  Xl=CS1*CB4*A1*RZ
		;  Yl=(CSSB*CB2-SS1*SB2)*A1*RZ
		;  Zl=(CSSB*SB2+SS1*CB2)*A1*RZ
		;  B=ALSO+ALSP*TM-PS1
		;  S2=2*SL1*CL1
		;  C2=CL1*CL1-SL1*SL1
		;  S3=S2*CL1+C2*SL1
		;  C3=C2*CL1-S2*SL1
		;  S4=S3*CL1+C3*SL1
		;  C4=C3*CL1-S3*SL1
		;  DLS=AM1*SL1+AM2*S2+AM3*S3+AM4*S4
		;  DR4=AE1*CL1+AE2*C2+AE3*C3+AE4*C4
		;  R4=ANU-DR4
		;  R4D=1/R4
		;  B=B+DLS-P1*R4D
		;  SB=sin(B)
		;  CB=cos(B)
		;  BET=A1*SS1*R4D*0.0121506682868;

#if defined (EPOCHA_IS_QUASINERTIAL)
	     D=DT0+TM+2921;
	     B=(12.1128-0.0529539*D)/RADIAN;
     d_mu =-76.7e-6*sin(B)
	    +0.9e-6*sin(2*B)
	    -5.7e-6*sin(2*(280.00812+0.9856473*D)/RADIAN)
	    -0.9e-6*sin(2*(64.3824+13.176396*D)/RADIAN);
    d_nu  =-33.3e-6*cos(B)
	    +0.4e-6*cos(2*B)
	    -2.5e-6*cos(2*(280.00812+0.9856473*D)/RADIAN)
	    -0.4e-6*cos(2*(64.3824+13.176396*D)/RADIAN);

d_epsilon  = 44.7e-6*cos(B)
	    -0.4e-6*cos(2*B)
	    +2.7e-6*cos(2*(280.00812+0.9856473*D)/RADIAN)
	    +0.4e-6*cos(2*(64.3824+13.176396*D)/RADIAN);
	       Xc=Xl+Yl*d_mu+Zl*d_nu;
	       Yc=-Xl*d_mu+Yl+Zl*d_epsilon;
	       Zc=-Xl*d_nu-Yl*d_epsilon+Zl;
		kappa=0.063107*D/(3600*RADIAN);
		omega=kappa;
		nu=0.054875*D/(3600*RADIAN);
		P[0]=-sin(kappa)*sin(omega)+cos(kappa)*cos(omega)*cos(nu);
		P[1]=-cos(kappa)*sin(omega)-sin(kappa)*cos(omega)*cos(nu);
		P[2]=-cos(omega)*sin(nu);
		P[3]=sin(kappa)*cos(omega)+cos(kappa)*sin(omega)*cos(nu);
		P[4]=cos(kappa)*cos(omega)-sin(kappa)*sin(omega)*cos(nu);
		P[5]=-sin(omega)*sin(nu);
		P[6]=cos(kappa)*sin(nu);
		P[7]=-sin(kappa)*sin(nu);
		P[8]=cos(nu);
		/*
		P-матрица прецессии осуществляет переход от эпохи 1950г к
		среднему равноденствию эпохи DT0+T
        здесь нужно разобраться (почему не наоборот)
		*/
		/*Вычисление координат Луны и Солнца в динамической системе*/

		*XL=P[0]*Xc+P[3]*Yc+P[6]*Zc;
		*YL=P[1]*Xc+P[4]*Yc+P[7]*Zc;
//было	*ZL=P[2]*Xl+P[5]*Yc+P[8]*Zc;
		*ZL=P[2]*Xc+P[5]*Yc+P[8]*Zc;
		 Xl=CB*R4*RZ;
		 Yl=(SB*CB2-BET*SB2)*R4*RZ;
		 Zl=(SB*SB2+BET*CB2)*R4*RZ;
	     Xc=Xl+Yl*d_mu+Zl*d_nu;
	     Yc=-Xl*d_mu+Yl+Zl*d_epsilon;
	     Zc=-Xl*d_nu-Yl*d_epsilon+Zl;
		*XC=P[0]*Xc+P[3]*Yc+P[6]*Zc;
		*YC=P[1]*Xc+P[4]*Yc+P[7]*Zc;
		*ZC=P[2]*Xc+P[5]*Yc+P[8]*Zc;
#else
		;*XL=Xl;
		*YL=Yl;
		*ZL=Zl;
		*XC=CB*R4*RZ;
		*YC=(SB*CB2-BET*SB2)*R4*RZ;
		*ZC=(SB*SB2+BET*CB2)*R4*RZ;
#endif

	    }


/*########################## Function KOOPLS ###########################*/

/*########################## Function UV ###########################*/
/* ======================================================================
	 Назначение: Вычисление по координатам возмущающего тела значений
							 шаровых функций.
	 Форма обращения: UV(X,Y,Z,LK,HLC)
	 Входная информация:
		Координаты возмущающего тела: X,Y,Z,
		число членов в разложении по параллаксу - Lк.
	 Выходная информация:
		Массив значений шаровых функций HLC[1:M, 1:2],
		где M=Lк*(Lк+1)/2-2+Lк.
	 Контрольный вариант:
	 Входная информация:
		X=-1.03811685765, Y=-.101398483053, Z=.038430850163, Lк=5.
	 Выходная информация:
		HLC[1,1]=-.97921271, HLC[2,1]=-.124726493, HLC[3,1]=1.6685198677,
		HLC[4,1]=-.122835649, HLC[5,1]=1.348477187, HLC[6,1]=.155723805,
		HLC[7,1]=-1.6842145413, HLC[8,1]=8.9581971e-1, HLC[9,1]=2.0965921e-1,
		HLC[10,1]=-1.31375875 и т.д.
		HLC[1,2]=0., HLC[2,2]=-1.21827e-2, HLC[3,2]=3.29086e-1,
		HLC[4,2]=0., HLC[5,2]=1.31713e-1, HLC[6,2]=3.071379e-2,
		HLC[7,2]=-5.0644447e-1, HLC[8,2]=0., HLC[9,2]=2.04785e-2,
		HLC[10,2]=-2.591159e-1 и т.д.
	 ------------------------------------------------------------------- */
/* ======================================================================
	 Destination: Calculation of spherical functions value.
     (see formula (1.3.38), Topic 5.2, Report 1).
	 Call: UV(X,Y,Z,LK,HLC)
	 Input data:
        Coordinates: X,Y,Z.
		Maximum order of spherical functions to be calculated - Lк.
	 Output data:
		Array of values of spherical functions HLC[1:M, 1:2],
		where M=Lк*(Lк+1)/2-2+Lк.
 ------------------------------------------------------------------- */

              ; void JVC_DLL UV(double X,double Y,double Z,int LK,double far HLC[])
		      {
		       int
			       L
			    ,  N
			    ,  I
			    ,  T
			    ,  T1
			    ,  A
			    ,  B
			    ,  C
			    ,  K
			    ,  T2;
		    long       SA
			    ,  SB;
		  double
			       R
			    ,  UD
			    ,  VD
			    ,  F
			    ,  R1
			    ,  RK
			       ;
		    double      far*U
		      ;  if((U=(double far*) calloc(48,sizeof(double)))==NULL)
     {printf("\n ");
     exit(1);
     }
		      ;  R=sqrt(X*X+Y*Y+Z*Z)
		      ;  R1=1.0/(R*R)
		      ;  U[0]=1.0/R
		      ;  U[1]=0
		      ; for (L= 1;L<= LK;L++)
			   {  RK=(2*(L-1)+1)*R1
			    ;  T=L*(L+1)-2
			    ;  T1=T+2*(L+1)
			    ;  U[T1]=RK*(X*U[T]-Y*U[T+1])
			    ;  U[T1+1]=RK*(X*U[T+1]+Y*U[T])
			;}
		      ; for  (N=0;N<LK;N++)
			 {
			 for (L= N;L<LK;L++)
			       {  T=(L+1)*(L+2)+2*N
				  ;  T1=T-2*(L+1)
				  ;  T2=T1-2*L
				  ; if  (N>L-1)
					{    UD=0.0;
					     VD=0.0
					 ;}
				    else {
					   UD=U[T2]
					;  VD=U[T2+1]
					 ;}
				  ;  K=2*L+1
				  ;  U[T]=(K*Z*U[T1]*R1-(L+N)*UD*R1)/(L-N+1)
				  ;  U[T+1]=(K*Z*U[T1+1]*R1-(L+N)*VD*R1)/(L-N+1)
				 ;}
			    ;}
		      ; for (N=0;N<=LK;N++)
			   {
			      C=(N<=2)? 2:N
			    ; for (L= C;L<=LK;L++)
				 {  A=L-N
				     ;  B=L+N
				     ;  SA=1
				     ;  SB=1
				     ;  K=2*L+1
				     ; for (I= 2;I<=B;I++)
				      {if (I<=A) SA*=I ; SB*=I;}
				  ;  R= (N!=0) ? 2.0 : 1.0
				  ;  F=sqrt((R*K*SA)/SB);
				     T=L*(L+1)-6+2*N;
				     HLC[T]=U[T+6]*F;
				    HLC[T+1]=U[T+7]*F;
			      ;}
			;}
			free(U)
		  ;}
/*########################## Function UV ###########################*/

/*########################## Function REZ ###########################*/
/* ======================================================================
	 Назначение: Вычисление по значениям большой полуоси А, эксцентриситета E
							 и наклонения I значений функций, определяющих резонансные
							 возмущения от тессеральных и секториальных гармоник разло-
							 жения гравитационного поля Земли.
	 Форма обращения: REZ().
	 Входная информация:
	 Большая полуось А(км.), эксцентриситет E, наклонение I(рад.), порядок
	 учитываемых гармоник Lт, геоцентрическая гравитационная постоянная GME,
	 скорость  вращения Земли CBZ, экваториальный радиус Земли RZ, массив
	 коэффициентов разложения геопотенциала HK1[1:300,1:2].
	 Выходная информация: JREZ[1:7],
		QC[1:7,-20:20], QD[1:7,-20:20], PC[1:7,-20:20], PD[1:7,-20:20],
		RC[1:7,-20:20], RD[1:7,-20:20], TC[1:7,-20:20], TD[1:7,-20:20].
	 Используется обращение к процедурам FINC,HANSEN.
	 Контрольный вариант:
	 Входная информация:
	 A=26569.13206, E=.01, I=PI/6, LM=LT=4.
	 Выходная информация: JREZ[1]=2, JREZ[2]=4,
		массив QC:
		QC[1:-3]=2.90795e-7, QC[1:-2]=1.76828e-13, QC[1:-1]=6.62596e-10,
		QC[1:0]=-1 2288e-8,  QC[1:1]=-3.238e-10,   QC[1:2]=-2.6922e-13,
		QC[1:3]=-2.2085e-16,                       QC[2:-2]=8.36358e-14,
		QC[2:-1]=0.,         QC[2:0]=4.87e-10,     QC[2:1]=0.,
		QC[2:2]=2.587e-14,   QC[2:3]=0.,           QC[2:4]=3.409e-20,
		массив QD:
    QD[1:-3]=-1.514e-17, QD[1:-2]=-1.1655e-13, QD[1:-1]=1.1727e-9,
    QD[1:0]=8 0997e-9,   QD[1:1]=-4.852e-10,   QD[1:2]=1.775e-13,
    QD[1:3]=8.288e-17,                         QD[2:-2]=4.404e-14,
                         QD[2:0]=2.56579e-10,
    QD[2:2]=1.362e-14,
    массив PC:
    PC[1:-3]=4.575e-15,  PC[1:-2]=2.3475e-11, PC[1:-1]=-1.18e-7,
    PC[1:0]=-3 2178e-10, PC[1:1]=4.88e-8,     PC[1:2]=-3.576e-11,
    PC[1:3]=-2.505e-14,  PC[1:4]=-4.97e-19 ,  PC[2:-2]=-8.87e-12,
                         PC[2:0]=-5.097e-12,
    PC[2:2]=-2.745e-12,                       PC[2:4]=-7.23e-18
    и т.д.
   ------------------------------------------------------------------- */
/* ======================================================================
	 Destination: Calculation of slow-changing resonance functions
     (see formulas (1.3.27), (1.3.28) Topic 5.2, Report 1).
	 Call: REZ().
	 It is used variables, declared as global:
	 А(semimajor axis, km), E (eccentricity),
     SI2(sin(I/2)),CI2 (cos(I/2)), where I is inclination.
	 Otput data: static array JREZ[0:MAXREZ-1],
        dynamic arrays
		QC[0:MAXREZ*(2*LMAX+1)-1], QD[0:MAXREZ*(2*LMAX+1)-1],
        PC[0:MAXREZ*(2*LMAX+1)-1], PD[0:MAXREZ*(2*LMAX+1)-1],
		RC[0:MAXREZ*(2*LMAX+1)-1], RD[0:MAXREZ*(2*LMAX+1)-1],
        TC[0:MAXREZ*(2*LMAX+1)-1], TD[0:MAXREZ*(2*LMAX+1)-1].
   ------------------------------------------------------------------- */


	      ;  void JVC_DLL REZ(void)
/* выход признак резонанса */ {
		      double   N
			    ,  Z
			    ,  A1
			    ,  BD
			    ,  AD
			    ,  K1
			    ,  H2
			    ,  F1
			    ,  F2
			    ,  G1
			    ,  G2;

		      int      M
			    ,  S
			    ,  P
			    ,  K
			    ,  J
			    ,  J1=0
			    ,  Q
			    ,  L
			       ;

			 N=sqrt(GME/A)/A*86400.0
		      ;  H2=sqrt(1-E*E)
		      ; for (J=1;J<=MAXREZ;J++)
			for (M=1;M<=LT;M++)          
			    {    Z=J-M*CBZ/N
			       ; if  (fabs(Z)<5.0E-2)
				  {
				  JREZ[J-1]=M
			     ; for (Q=J-LT;Q<=J+LT;Q++)
				  {

				 ;  QC[J1]=QD[J1]=RC[J1]=RD[J1]=PC[J1]=PD[J1]=TC[J1]=TD[J1]=0.0
			       ;  S=abs(J-Q)
			       ;  S=MAX(M,S)
			       ; for  (L= S;L<=LM;L++)
				    {  A1=pow(RZ/A,L)
			       ;  P=L+Q-J
			       ; if (!(P%2) && (L!=1))
				      {    K=L*(L+1)-6+2*M
			       ;  P/=2
			       ; if  (!((L-M)%2))
				    {    BD=HK1[K+1]
				       ; AD=-HK1[K]
				    ; }  else
				 {
				     BD=HK1[K]
			      ;  AD=HK1[K+1]
				 ;}

			      ;  FINC(L,M,P,SI2,CI2,&F1,&F2)
			      ;   HANSEN(-L-1,J-Q,J,E,&G1,&G2);
			      ;   K1=A1*F1*G1
			      ;  QC[J1]+=K1*BD
			      ;  QD[J1]+=K1*AD
			      ;  K1=A1*F1*G2
			      ;  PC[J1]-=K1*AD
			      ;  PD[J1]+=K1*BD
			      ;  K1=0.5*A1*F2*G1/CI2/H2
			      ;  RC[J1]-=K1*AD
			      ;  RD[J1]+=K1*BD
			      ;  K1=2*A1*(L+1)*F1*G1
			      ;  TC[J1]-=K1*AD
			      ;  TD[J1]+=K1*BD
				      ;}
				    ;}
				  J1++
				  ;}


			      ;}
			;}
		  ;}

/*########################## Function REZ ###########################*/
/*########################## Function LS ###########################*/
/* ======================================================================
   Назначение: Вычисление значений функций,определяющих вековые и долго-
               периодические возмущения от притяжения внешних тел.
   Форма обращения:LS().
   Входная информация:
    Эксцентриситет Е, наклонение I,число параллактических членов - Lк.
   Выходная информация:
    Массивы Lе[2:Lк, 0:Lк,0:Lк,], Li[2:Lк, 0:Lк,0:Lк,], Lv[2:Lк, 0:Lк,0:Lк,],
    Lп[2:Lк, 0:Lк,0:Lк,],Lл[2:Lк, 0:Lк,0:Lк,].
   Используется обращение к процедурам FINC,HANSEN.
   Контрольный вариант:
   Входная информация:
    E=.0099371046, I=.523587802, Lк=5.
   Выходная информация:
    LE[2,0,0]=1.0414734e-2, LE[3,0,0]=-1.3393e-4, LE[4,0,0]=-1.5845e-6,
    LE[5,0,0]=1.795e-8, LE[4,1,0]=7.48018e-6, LE[5,1,0]=8.649e-8,
    LE[2,2,0]=-8.3753e-2, LE[3,2,0]=1.445e-3, LE[4,2,0]=1.974e-5,
		LE[5,2,0]=-2.44e-7, и т.д.
   ------------------------------------------------------------------- */
   /* ======================================================================
   Destination: Calculation of slow-changing functions for taken into
               account perturbations due to attraction of the Moon and Sun.
   Call:LS().
   It is used variables, declared as global:
   E (eccentricity), SI2(sin(I/2)), CI2 (cos(I/2)), where I is inclination,
   maximum order of parallactic terms LLUN.
   Otput data:
    Dynamic arrays: LE[0:(LLUN-1)*(LLUN+1)*(LLUN+1)-1], LI[0:(LLUN-1)*(LLUN+1)*(LLUN+1)-1],
    LV[0:(LLUN-1)*(LLUN+1)*(LLUN+1)-1],LPI[0:(LLUN-1)*(LLUN+1)*(LLUN+1)-1],
    LLAM[0:(LLUN-1)*(LLUN+1)*(LLUN+1)-1].
   ------------------------------------------------------------------- */

           ; void JVC_DLL LS(void)
		    {
		     double
			       H1
			    ,  H2
			    ,  S
			    ,  F1
			    ,  F2
			    ,  T1
                ,  T11
			    ,  T2
			    ,  T3
			    ,  T4
			    ,  T5
			    ,  T6
			    ,  T7;
	      int	   L
			    ,  T
			    ,  M
			    ,  P
			    ,  J=0  ;
			     S=sqrt(1-E*E)
		      ;  T1=-S/E
		      ;  T3=1.0/(4*S*SI2)
		      ;  T4=1.0/(2*S*CI2)
              ;  T7=2*SI2
		      ;  T5=E*T7
		      ;  T6=(S-1+E*E)/E
		      ; for  (L=2;L<=LLUN;L++)
			   { for (P=0;P<=L;P++)
				 {   T=L-2*P
				  ;  T2=T*(1-T7*SI2)
				  ;  HANSEN(L,T,0,E,&H1,&H2)
				  ; for (M=0;M<=L;M++)
				    {  FINC(L,M,P,SI2,CI2,&F1,&F2)
					;  LE[J]=T*F1*H1*T1
					;  LI[J]=(T2-M)*F1*H1*T3
					;  LV[J]=T11=F2*H1*T4
					;  LPI[J]=S*F1*H2+T5*T11
					;  LLAM[J++]=T6*F1*H2-(2*L)*F1*H1+T7*T11
				   ;}
			     }
			  ;}

		  ;}
/*########################## Function LS ###########################*/
/*########################## Function APPROK ###########################*/
/* ======================================================================
   Назначение: Аппроксимация шаровых функций*/
/* ======================================================================
   Destination: Approximation of spherical functions.
   Call: APPROK(DT,T,TK,THAX,TKOH).
   Input data:
         DT is step of approximation,
         T is initial time,
         TK is final time.
   Output data:
   Dynamical arrays HL,DL,D2L of dimension ((LLUN+1)*(LLUN+2)-6),
                    HC,DC,D2C of dimension ((LSUN+1)*(LSUN+2)-6).
   ------------------------------------------------------------------- */
 void JVC_DLL APPROK(double DT,double T,double TK,double far *THAX,double far *TKOH)
 {
      double far*ML1,
	     far*ML2,
	     far*ML3,
	     far*ML4,
	     far*MC1,
	     far*MC2,
	     far*MC3,
	     far*MC4;
      double
	     XL,
	     YL,
	     ZL,
	     XC,
	     YC,
	     ZC,
	     X ;
       int K;
			if  ((T<*THAX) || (T>*TKOH))
			    {   if  (TK>=0)
				 {   *THAX=T
				  ;  *TKOH=T+4*DT
				 ;} else {*THAX=T-4*DT
					   ;*TKOH=T
					 ;}

	       if(
		(ML1=(double far*) calloc(LKL,sizeof(double)))==NULL||
		(ML2=(double far*) calloc(LKL,sizeof(double)))==NULL||
		(ML3=(double far*) calloc(LKL,sizeof(double)))==NULL||
		(ML4=(double far*) calloc(LKL,sizeof(double)))==NULL||
		(MC1=(double far*) calloc(LKS,sizeof(double)))==NULL||
		(MC2=(double far*) calloc(LKS,sizeof(double)))==NULL||
		(MC3=(double far*) calloc(LKS,sizeof(double)))==NULL||
		(MC4=(double far*) calloc(LKS,sizeof(double)))==NULL)
		{printf("\nAPPROK");
		 exit(1);
		 }

			    ; for  (K=0;K<=4;K++)
				 {  KOOPLS(*THAX+K*DT,&XL,&YL,&ZL,&XC,&YC,&ZC)
				  ;  XL/=ALUN
				  ;  YL/=ALUN
				  ;  ZL/=ALUN
				  ;  XC/=AC
				  ;  YC/=AC
				  ;  ZC/=AC
				  ; switch(K)
				     {
				       case 0:{  UV(XL,YL,ZL,LLUN,HL0)
				     ;     UV(XC,YC,ZC,LSUN,HC0) ;break;}
				     ; case 1:{  UV(XL,YL,ZL,LLUN,ML1)
				     ;     UV(XC,YC,ZC,LSUN,MC1) ;break;}
				     ; case 2:{  UV(XL,YL,ZL,LLUN,ML2)
				     ;     UV(XC,YC,ZC,LSUN,MC2) ;break;}
				     ; case 3:{  UV(XL,YL,ZL,LLUN,ML3)
				     ;     UV(XC,YC,ZC,LSUN,MC3) ;break;}
				     ; case 4:{  UV(XL,YL,ZL,LLUN,ML4)
				     ;     UV(XC,YC,ZC,LSUN,MC4) ;break;}
				    ;}
			      ;}
			    ; for  (K=0;K<LKL;K++)
				  {  HL1[K]=(-25.0/12.0*HL0[K]+4*ML1[K]-3*ML2[K]+4.0/3.0*ML3[K]-0.25*ML4[K])/DT
				  ;  HL2[K]=0.5*(35.0/12.0*HL0[K]-26.0/3.0*ML1[K]+9.5*ML2[K]-14.0/3.0*ML3[K]+11.0/12.0*ML4[K])/(DT*DT)
				  ;  HL3[K]=(-2.5*HL0[K]+9*ML1[K]-12*ML2[K]+7*ML3[K]-1.5*ML4[K])/(6*DT*DT*DT)
				  ;  HL4[K]=(HL0[K]-4*ML1[K]+6*ML2[K]-4*ML3[K]+ML4[K])/(24*DT*DT*DT*DT)
				  ; if  (K<LKS)
				      {    HC1[K]=(-25.0/12.0*HC0[K]+4*MC1[K]-3*MC2[K]+4.0/3.0*MC3[K]-0.25*MC4[K])/DT
					;  HC2[K]=0.5*(35.0/12.0*HC0[K]-26.0/3.0*MC1[K]+9.5*MC2[K]-14.0/3.0*MC3[K]+11.0/12.0*MC4[K])/(DT*DT)
					;  HC3[K]=(-2.5*HC0[K]+9*MC1[K]-12*MC2[K]+7*MC3[K]-1.5*MC4[K])/(6*DT*DT*DT)
					;  HC4[K]=(HC0[K]-4*MC1[K]+6*MC2[K]-4*MC3[K]+MC4[K])/(24*DT*DT*DT*DT)
				      ;}
			      ;}
		free(ML1);
		free(ML2);
		free(ML3);
		free(ML4);
		free(MC1);
		free(MC2);
		free(MC3);
		free(MC4);

			;}
		      ;  X=T-*THAX
		      ; for   (K=0;K<LKL;K++)
			    {  HL[K]=HL0[K]+(((HL4[K]*X+HL3[K])*X+HL2[K])*X+HL1[K])*X
			    ;  DL[K]=HL1[K]+((4*HL4[K]*X+3*HL3[K])*X+2*HL2[K])*X
			    ;  D2L[K]=2*HL2[K]+(12*HL4[K]*X+6*HL3[K])*X
			    ; if  (K<LKS)
				{    HC[K]=HC0[K]+(((HC4[K]*X+HC3[K])*X+HC2[K])*X+HC1[K])*X
				  ;  DC[K]=HC1[K]+((4*HC4[K]*X+3*HC3[K])*X+2*HC2[K])*X
				  ;  D2C[K]=2*HC2[K]+(12*HC4[K]*X+6*HC3[K])*X
			       ;}
			    ;}
		  ;}

/*########################## Function APPROK ###########################*/
/*########################## Function SHADOW ###########################*/
/* ======================================================================
   Назначение: Вычисление углов входа и выхода КО из тени*/
/* ======================================================================
   Destination: Calculation of shadow entry and exit parameters
   -----------------------------------------------------------------------*/

              ; void JVC_DLL SHADOW(double al,double del,double AP,double DBY,double *E1/*вход*/,double *E2/*выход*/,double *S,double *T,double *W)

		      { double c1,c2,cc,u1,u2,r1,r2,be;
			c1=cos(al-DBY)*cos(del);
			c2=sin(al-DBY)*cos(del)*cos(I)+sin(del)*sin(I);
			*S=-c1*cos(AP)-c2*sin(AP);
			*T=c1*sin(AP)-c2*cos(AP);
			*W=-sin(del)*cos(I)+cos(del)*sin(I)*sin(al-DBY);
			cc=sqrt(c1*c1+c2*c2);
			be=atan2(c1,c2);
			u1=-be;
			u2=u1+0.5*PI2;
                        do
                        {
                        r1=A*(1-E*E)/(1+E*cos(u1-AP));
                        r2=A*(1-E*E)/(1+E*cos(u2-AP));
                        c1=-sqrt(1-sqr(RZ/r1));
                        c2=-sqrt(1-sqr(RZ/r2));
                        if(fabs(c1/cc)>=1.0||fabs(c2/cc)>=1.0)
                        { *E2=0;*E1=PI2;break;}
                        r1=u1;
                        r2=u2;
                        u1=-be+asin(c1/cc);
                        u2=0.5*PI2-be-asin(c2/cc);
                        r1-=u1;
                        r1=fmod(fabs(r1),PI2);
                        if(r1>0.5*PI2) r1-=PI2;
                        r2=fmod(fabs(r2),PI2);
                        r2-=u2;
                        if(r2>0.5*PI2) r2-=PI2;
                        if(fabs(u1-r1)>0.0001||fabs(u2-r2)>0.0001)
//			if(r1<0.01&&r2<0.01)
			{
			 u1-=AP;
			 u2-=AP;
			 r1=2*atan2(sqrt((1-E)/(1+E))*sin(0.5*u1),cos(0.5*u1));
			 r2=2*atan2(sqrt((1-E)/(1+E))*sin(0.5*u2),cos(0.5*u2));
			 r1=fmod(r1,PI2);
			 r1=fmod(r1,PI2);
			 *E1=(r1<r2) ? r1+PI2 : r2+PI2;
			 *E2=(r1<r2) ? r2 : r1;
			 if(*E2-*E1<0.5*PI2)
			 {u1=*E1-PI2;*E1=*E2;*E2=u1;}
			break;
			}


			}while (1>0) ;


		      }
/*########################## Function SHADOW ###########################*/


#endif
