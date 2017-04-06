#include "hoko.h"

FILE  far *MAEIKP, far *MULTIINDEX;


extern int	   LM
		,  LT
		,  LKL
		,  LKS
		,  JREZ[MAXREZ]
		,  L
		,  M
		,  PLUNASUN

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
		,  AAL
		,  AAC
		,  TK
		,  S0
		,  THAX
		,  TKOH
		   ;

 extern double
		 AX
	      ,  EX
	      ,  IX
#if BOKO
	      ,	far*HL
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
#endif

	      ;
double
			       G[LMAX+1][2*LMAX+1][2]
			     , far*FF
#if BOKO
			    ,  GK[LLUN+1][2*QMAX+1][2]
#endif
			    ,  DKP[12]
			    ,  MCM[LMAX+1]
			    ;
double
			       Y[6]      /*     [6]      */
			    ,  EK[6]      /*     [6]      */
			    ;
		 double        W
			    ,  U
			    ,  V
			    ,  X
			    ,  K13
			    ,  MK
			    ,  AK
			    ,  AP
			    ,  DBY
			    ,  N13
			    ,  CA
			    ,  H2
			    ,  N
			    ,  Z
			    ,  CM
			    ,  SK
			    ,  CK
			    ,  FF1
			    ,  H1
			    ,  H3
			    ,  N1
			    ,  N2
			    ,  N3
			    ,  N4
			    ,  N5
			    ,  N6
			    ,  N7
			    ,  N8
			    ,  N9
			    ,  N10
			    ,  FF2
			    ,  CLM
			    ,  SLM
			    ,  K0
			    ,  K1
			    ,  K2
			    ,  K3
			    ,  K4
			    ,  K5
			    ,  K6
			    ,  K7
			    ,  K8
			    ,  K9
			    ,  F1
			    ,  F2
			    ,  N11
			    ,  N12
			    ,  SN
			    ,  TI2
			    ,  DOW
			    ,  DW
			    ,  DM
			    ,  G1
			    ,  G2
			    ,  RADIUS
                ,  dS
			       ;
		    int
			       CQ
			    ,  P
			    ,  PF
			    ,  QMM
			    ,  J1
			    ,  K
			    ,  LK
			    ,  J
			    ,  L0
			    ,  L1
			    ,  Q
			    ,  S
			    ,  N_korper
			       ;

/*########################## НАЧАЛО FINC2 ###########################*/
/* ======================================================================
   Destination: Calculation of all inclination and its derivatives up to order L2.
   Call: FINC2(L2).
   It is used the global variables
    SI2=sin(.5*i); CI2=cos(.5*i);
   Output data: dynamic array FF with value of inclination function
                and its derivatives.
  ------------------------------------------------------------------- */

		     ;void  FINC2(int L2)
		     {
		     double
			    A1,A11,A2,A22,A3,A33,R1,R2,R3,R4,R5,R6;
		     int
			    J,N,K,I,M,P;

			R1=SI2*CI2;
			R2=SI2*SI2;
			R3=CI2*CI2;
			FF[0]=1.0;
			FF[1]=0.0;
			FF[9]=FF[2]=R1;
			FF[3]=(R3-R2)*0.5;
			FF[7]=FF[4]=-R1;
			FF[5]=-FF[3];
			FF[6]=R3;
			FF[8]=R2;
			R1=2*R1;R2=R3-R2;
			for (N=2;N<=L2;N++)
			for (M=0;M<=N;M++)
		 {
		  if (N!=M)
		       { K=N-M+1;R4=1;
			if (K%2)  R4=-R4;
			R4*=(2.0*N-1.0)/(2.0*(N-M));
			R6=(N+M-1.0)/(N-M);
			for (P=0;P<=N;P++)
			 {
			  if (P==0 || P>(N-1) || M==(N-1))
			  { A2=A22=0.0 ;}
			   else
			{
			  J=((N-2)*(N-1)*(2*N-3))/3+2*((N-1)*M+P-1)
			;
			 A2=FF[J];A22=FF[J+1]
			;}
			if (P>(N-1))
			 { A1=A11=0.0 ;}
			  else
			 {
			  J=((N-1)*N*(2*N-1))/3+2*(N*M+P)
			 ;
			  A1=FF[J];A11=FF[J+1] ;};
			if (P==0)
			 { A3=A33=0.0;}
			  else
			 {
			  J=((N-1)*N*(2*N-1))/3+2*(N*M+P-1)
			; A3=FF[J];A33=FF[J+1]
			 ;}
			 R5=R4*(A1-A3);
			;J=(N*(N+1)*(2*N+1))/3+2*((N+1)*M+P)
			;
			FF[J]=R5*R1-R6*A2;
			FF[J+1]=R5*R2+R4*R1*(A11-A33)-R6*A22
			;}
		  ;} else
			for (P=0;P<=N;P++)
			 {
			if (P>(N-1))
			 { A1=A11=0.0 ;}
			else
			 {
			 J=((N-1)*N*(2*N-1))/3+2*(N*(N-1)+P)
			;
			 A1=FF[J];A11=FF[J+1] ;};
			if (P==0)
			 { A3=A33=0.0 ;} else
			{
			 J=((N-1)*N*(2*N-1))/3+2*(N*(N-1)+P-1)
			; A3=FF[J];A33=FF[J+1] ;};
			;J=(N*(N+1)*(2*N+1))/3+2*((N+1)*M+P)
			;
			FF[J]=0.5*(2*N-1)*((1-R2)*A3+(1+R2)*A1);
			FF[J+1]=0.5*(2*N-1)*(R1*(A3-A1)+(1-R2)*A33+(1+R2)*A11)
			;}
		 ;};

			for (N=2;N<=L2;N++)
			for (M=0;M<=N;M++)
		   { R3= (M==0) ? 1.0 : 2.0;
			R4=R5=1.0;
			for (I=1;I<=(N-M);I++)
			    R4*=I;
			for (I=1;I<=(N+M);I++)
			    R5*=I;
			R3=sqrt(R3*R4*(2*N+1)/R5);
			for (P=0;P<=N;P++)
			{
			J=(N*(N+1)*(2*N+1))/3+2*((N+1)*M+P)
			;FF[J]*=R3;
			FF[J+1]*=R3
			;}
		   ;}

	  ;}
/*########################## КОНЕЦ FINC2 ###########################*/
/* ======================================================================
   Destination: Calculation of short periodical perturbations.
   Call: SKP(PKP,EBX,T,E1).
   Input data: PKP is option for calculation of short periodical perturbations
              at initial or final points,
              EBX[0:5] is vector of nonsingular elements,
              T is the time counted off from 0h UTC Dec 31 1957.
   Output data: E1[0:5] is vector of short periodical perturbations
  ------------------------------------------------------------------- */

		  int JVC_DLL SKP(int PKP,double far EBX[],double T,double far E1[])
		{
			 SN=S0+CBZ*T
		      ;  A=EBX[0]
		      ;  E=sqrt(sqr(EBX[1])+sqr(EBX[2]))
		      ;  SI2=sqr(EBX[3])+sqr(EBX[4])
		      ;  CI2=sqrt(1-SI2)
		      ;  CI=1-2*SI2
		      ;  SI2=sqrt(SI2)
		      ;  SI=2*SI2*CI2
		      ;  Z=SI*SI
		      ;  TI2=SI2/CI2
		      ;  DBY=atan(EBX[4]/EBX[3])
		      ;  if (EBX[3]<0)  DBY+=PI
		      ;  AP=atan(EBX[2]/EBX[1])
		      ;  if (EBX[1]<0)  AP+=PI
		      ;  CA=EBX[5]-AP
		      ;  AP-=DBY
		      ;  EA=KEPLER(CA,E)
		      ;  CE=cos(EA)
		      ;  SE=sin(EA)
		      ;  RADIUS=1-E*CE
		      ;  H2=sqrt(1-E*E)
		      ;  H1=1.0/H2
		      ;  H3=1.0/(1-E*E)
		      ;  SIA=H2*SE/RADIUS
		      ;  CIA=(CE-E)/RADIUS
		      ;  RADIUS*=A
		      ;  IA=atan(SIA/CIA)
		      ;  if (CIA<0)  IA+=PI
		      ;  U=AP+IA
		      ;  N=sqrt(GME/A)/A*86400.0
		      ;  N1=cos(AP)
		      ;  N2=sin(AP)
		      ;  G1=H3*H3
		      ;  K6=EBX[1]/E
		      ;  K7=EBX[2]/E
		      ;  K8=EBX[3]/SI2
		      ;  K9=EBX[4]/SI2
		      ;  CM=sqr(RZ/A)
		      ;  IA=fmod(IA,PI2)
		      ;  if (IA<0)  IA+=PI2
		      ;  CA=fmod(CA,PI2)
		      ;  if (CA<0)  CA+=PI2
		      ;  for (L=0;L<=5;E1[L]=0,L++)
			    ;  N5=N1*CIA-N2*SIA
			    ;  N6=CIA*N2+N1*SIA
			    ;  K1=C20*(sqr((1+E*CIA)*H3)*(1+E*CIA)*H3*(3*Z*N6*N6-1)-0.5*(3*Z-2)*H1*H3)
			    ;  K2=C20*(3*(N5*N5-N6*N6)+E*((4*N5*N5-3)*N5*N1+(3-4*N6*N6)*N6*N2)+3*E*(N5*N1-N2*N6))*H3*0.25
			    ;  N7=N1*N1-N2*N2
			    ;  N8=2*N1*N2
			    ;  G2=CM*C20
			    ;  FF1=Z*Z
			    ;  E1[0]=A*K1*CM*(1-K1*CM-0.75*G2*(3*Z-2)*H3*H1)-1.5*CM*G2*A*(4-5*Z)*Z*K2*G1*H1-3.0/64.0*G2*G2*A*G1*
       H3*H1*((-8+8*Z+5*FF1)*H2*H2+4*(4-12*Z+9*FF1)*H2-5*(-8+16*Z-7*FF1)+2*(-14*Z+15*FF1)*E*E*N7)
			    ;  F1=CM*(0.5*K1*H2*H2+Z*K2)/E
			    ;  FF1=Z/E
			    ;  FF2=Z*E
			    ;  N1=CIA*CIA-SIA*SIA
			    ;  N2=2*CIA*SIA
			    ;  N5=CIA*N1-SIA*N2
			    ;  N6=N1*SIA+N2*CIA
			    ;  N11=N8*N5+N6*N7
			    ;  N9=N8*CIA+SIA*N7
			    ;  N10=N8*CIA-SIA*N7
			    ;  K1=1.5*G2*CI*G1*(IA-CA+E*SIA-0.5*(N8*N1+N7*N2)-E/6*N11-E/2*N9)
			    ;  N12=N8*(N1*N5-N2*N6)+(N2*N5+N1*N6)*N7
			    ;  N13=N8*(N1*N1-N2*N2)+2*N2*N1*N7
    ;  K3=0.75*G2*G1*((5*Z-4)*(IA-CA)+(3*FF1+17.0/4.0*FF2-2/E-3.5*E)*SIA+0.125*FF2*N10+(0.5*FF1-15.0/8.0*FF2+E)*N9+
       (1.5*Z-1)*N2+(1-2.5*Z)*(N8*N1+N7*N2)+(0.25*FF2-E/6)*N6+(E/3-7*FF1/6-19.0/24.0*FF2)*N11-
       0.75*Z*(N13-N8)-0.125*FF2*N12)
			    ;  K13=K1+K3
			    ;  E1[1]=F1*K6-E*K13*K7
			    ;  E1[2]=F1*K7+E*K13*K6
			    ;  F1=-CM*SI*CI*K2*H3
			    ;  E1[3]=CI2*K8*F1/2-K1*SI2*K9
			    ;  E1[4]=CI2*K9*F1/2+K1*SI2*K8
; F1=0.75*G2*H3*H1*((0.75*FF2-3*FF1+2/E-0.5*E)*SIA+(1-1.5*Z)*N2+(E/6-FF2/4)*N6-0.125*FF2*N10-(0.5*FF1+
 5.0*FF2/8.0)*N9+(7.0*FF1/6.0-FF2/24)*N11+0.75*Z*(N13-N8)+0.125*FF2*N12)
			    ;  E1[5]=K13+F1
			    ;  CM=G2*G1
			    ;  DOW=1.5*CM*CI
			    ;  DW=0.75*CM*(5*Z-4)
			    ;  DM=0.75*CM*H2*(3*Z-2)

			    ; if  (PKP==0)
                {
              for (L=0;L<=5;L++) EK[L]=EBX[L]-E1[L]
		      ;  A=EK[0]
		      ;  E=sqrt(sqr(EK[1])+sqr(EK[2]))
		      ;  SI2=sqr(EK[3])+sqr(EK[4])
		      ;  CI2=sqrt(1-SI2)
		      ;  CI=1-2*SI2
		      ;  SI2=sqrt(SI2)
		      ;  SI=2*SI2*CI2
		      ;  Z=SI*SI
		      ;  TI2=SI2/CI2
		      ;  DBY=atan(EK[4]/EK[3])
		      ;  if (EK[3]<0)  DBY+=PI
		      ;  AP=atan(EK[2]/EK[1])
		      ;  if (EK[1]<0)  AP+=PI
		      ;  CA=EK[5]-AP
		      ;  AP-=DBY
		      ;  EA=KEPLER(CA,E)
		      ;  CE=cos(EA)
		      ;  SE=sin(EA)
		      ;  RADIUS=1-E*CE
		      ;  H2=sqrt(1-E*E)
		      ;  H1=1.0/H2
		      ;  H3=1.0/(1-E*E)
		      ;  SIA=H2*SE/RADIUS
		      ;  CIA=(CE-E)/RADIUS
		      ;  RADIUS*=A
		      ;  IA=atan(SIA/CIA)
		      ;  if (CIA<0)  IA+=PI
		      ;  U=AP+IA
		      ;  N=sqrt(GME/A)/A*86400.0
		      ;  N1=cos(AP)
		      ;  N2=sin(AP)
		      ;  G1=H3*H3
		      ;  K6=EK[1]/E
		      ;  K7=EK[2]/E
		      ;  K8=EK[3]/SI2
		      ;  K9=EK[4]/SI2
		      ;  CM=sqr(RZ/A)
		      ;  IA=fmod(IA,PI2)
		      ;  if (IA<0)  IA+=PI2
		      ;  CA=fmod(CA,PI2)
		      ;  if (CA<0)  CA+=PI2
		      ;  for (L=1;L<=5;E1[L]=0,L++)
			    ;  N5=N1*CIA-N2*SIA
			    ;  N6=CIA*N2+N1*SIA
			    ;  K1=C20*(sqr((1+E*CIA)*H3)*(1+E*CIA)*H3*(3*Z*N6*N6-1)-0.5*(3*Z-2)*H1*H3)
			    ;  K2=C20*(3*(N5*N5-N6*N6)+E*((4*N5*N5-3)*N5*N1+(3-4*N6*N6)*N6*N2)+3*E*(N5*N1-N2*N6))*H3*0.25
			    ;  N7=N1*N1-N2*N2
			    ;  N8=2*N1*N2
			    ;  G2=CM*C20
			    ;  FF1=Z*Z
//			    ;  E1[0]=A*K1*CM*(1-K1*CM-0.75*G2*(3*Z-2)*H3*H1)-1.5*CM*G2*A*(4-5*Z)*Z*K2*G1*H1-3.0/64.0*G2*G2*A*G1*
//       H3*H1*((-8+8*Z+5*FF1)*H2*H2+4*(4-12*Z+9*FF1)*H2-5*(-8+16*Z-7*FF1)+2*(-14*Z+15*FF1)*E*E*N7)
			    ;  F1=CM*(0.5*K1*H2*H2+Z*K2)/E
			    ;  FF1=Z/E
			    ;  FF2=Z*E
			    ;  N1=CIA*CIA-SIA*SIA
			    ;  N2=2*CIA*SIA
			    ;  N5=CIA*N1-SIA*N2
			    ;  N6=N1*SIA+N2*CIA
			    ;  N11=N8*N5+N6*N7
			    ;  N9=N8*CIA+SIA*N7
			    ;  N10=N8*CIA-SIA*N7
			    ;  K1=1.5*G2*CI*G1*(IA-CA+E*SIA-0.5*(N8*N1+N7*N2)-E/6*N11-E/2*N9)
			    ;  N12=N8*(N1*N5-N2*N6)+(N2*N5+N1*N6)*N7
			    ;  N13=N8*(N1*N1-N2*N2)+2*N2*N1*N7
    ;  K3=0.75*G2*G1*((5*Z-4)*(IA-CA)+(3*FF1+17.0/4.0*FF2-2/E-3.5*E)*SIA+0.125*FF2*N10+(0.5*FF1-15.0/8.0*FF2+E)*N9+
       (1.5*Z-1)*N2+(1-2.5*Z)*(N8*N1+N7*N2)+(0.25*FF2-E/6)*N6+(E/3-7*FF1/6-19.0/24.0*FF2)*N11-
       0.75*Z*(N13-N8)-0.125*FF2*N12)
			    ;  K13=K1+K3
			    ;  E1[1]=F1*K6-E*K13*K7
			    ;  E1[2]=F1*K7+E*K13*K6
			    ;  F1=-CM*SI*CI*K2*H3
			    ;  E1[3]=CI2*K8*F1/2-K1*SI2*K9
			    ;  E1[4]=CI2*K9*F1/2+K1*SI2*K8
; F1=0.75*G2*H3*H1*((0.75*FF2-3*FF1+2/E-0.5*E)*SIA+(1-1.5*Z)*N2+(E/6-FF2/4)*N6-0.125*FF2*N10-(0.5*FF1+
 5.0*FF2/8.0)*N9+(7.0*FF1/6.0-FF2/24)*N11+0.75*Z*(N13-N8)+0.125*FF2*N12)
			    ;  E1[5]=K13+F1
			    ;  CM=G2*G1
			    ;  DOW=1.5*CM*CI
			    ;  DW=0.75*CM*(5*Z-4)
			    ;  DM=0.75*CM*H2*(3*Z-2);
                }

			    ; if  (PKP!=0)
			       { for (L=0;L<=5;L++)
				     EK[L]=EBX[L]+E1[L]
				  ;  LBK(EK,Y)
				  ;  K0=Y[1]
				  ;  IA=Y[5]
				  ;  EA=KEPLER(IA,K0)
				  ;  CE=cos(EA)
				  ;  SE=sin(EA)
				  ;  K4=sqrt(1-K0*K0)
				  ;  K3=1/K4
				  ;  K5=1/(1-K0*K0)
				  ;  CIA=1-K0*CE
				  ;  SIA=K4*SE/CIA
				  ;  CIA=(CE-K0)/CIA
				  ;  IA=atan(SIA/CIA)
				  ;  if (CIA<0)  IA+=PI
				  ;  Z=sqr(sin(Y[2]))
				  ;  N1=cos(Y[3])
				  ;  N2=sin(Y[3])
				  ;  G1=K5*K5
				  ;  CM=sqr(RZ/Y[0])
				  ;  N5=N1*CIA-N2*SIA
				  ;  N6=CIA*N2+N1*SIA
				  ;  K1=C20*(sqr((1+K0*CIA)*K5)*((1+K0*CIA)*K5)*(3*Z*N6*N6-1)-0.5*(3*Z-2)*K3*K5)
				  ;  K2=C20*(3*(N5*N5-N6*N6)+K0*((4*N5*N5-3)*N5*N1+(3-4*N6*N6)*N6*N2)+3*K0*(N5*N1-N2*N6))*K5*0.25
				  ;  G2=CM*C20
				  ;  N8=2*N1*N2
				  ;  G2=CM*C20
				  ;  FF1=Z*Z
		     ;  E1[0]=Y[0]*K1*CM*(1-K1*CM-0.75*G2*(3*Z-2)*K5*K3)-1.5*CM*G2*Y[0]*(4-5*Z)*Z*K2*G1*K3-3.0/
		      64.0*G2*G2*Y[0]*G1*K5*K3*((-8+8*Z+5*FF1)*K4*K4+4*(4-12*Z+9*FF1)*K4-5*(-8+16*Z-7*FF1)+2*(-14*
		      Z+15*FF1)*K0*K0*N7)
			      ;}
 /*Обход расчета короткой периодики -24.02.93*/

				   ; if ((LM>2) || (LT!=0) || PLUNASUN)
				   {
				    L=((LM>LLUN)? LM:LLUN)
				  ;  J1=((L+1)*(L+2)*(2*L+3)) / 3
				  ;  FF=(double far*) calloc(J1,sizeof(double));
				  if (FF==NULL) {printf("\n Нет памяти для FF  в процедуре SKP\n");return(3);}

				  ; if (PLUNASUN)  FINC2(L)

			 ; if ((LM!=2) || (LT!=0))

			   {


			      if ((fabs(1-A/AX)>CTPOBKP)
				    || (fabs(1-E/(EX+10e-10))>CTPOBKP) || (fabs(1-SI2/(IX+10e-10))>CTPOBKP))
			       {
				  ; if (!PLUNASUN)  FINC2(L)
				  ;  rewind(MULTIINDEX)
				  ;  rewind(MAEIKP)
				  ;  N_korper=0
				  ;  AX=A
				  ;  EX=E
				  ;  IX=SI2
				  ; fwrite(&LM,sizeof(LM),1,MULTIINDEX)
				  ; fwrite(&LT,sizeof(LT),1,MULTIINDEX)
				  ; fseek(MULTIINDEX,3*sizeof(LM),SEEK_SET)
				  ; fwrite(&AX,sizeof(AX),1,MAEIKP)
				  ; fwrite(&EX,sizeof(EX),1,MAEIKP)
				  ; fwrite(&IX,sizeof(IX),1,MAEIKP)
				  ; G1=RZ/A
				  ; MCM[1]=G1*G1
				  ; for (L= 2;L<LM;L++) MCM[L]=G1*MCM[L-1]
				  ;  N1=H1/(2*CI2)
				  ;  N2=TI2*H1
				  ;  N3=E*N2
				  ;  N4=(H2-1+E*E)/E
				  ;  QMM=MIN((-15/log(E)+1),QMAX)
				  ; for (J=0;J<=QMM;J++)
				      { for (CQ=0;CQ<=MIN(1,J);CQ++)
					    {
						if  (CQ==0) Q=J;	else Q=-J
					      ; for  (S=-LM;S<=LM;S++)
						   {  N6=H2*(H2*(S+Q)-S)
						    ; for (M=0;M<=LT;M++)
							 {
							    if((CQ==0) && (M==0))
							     for  (L= 2;L<=LM;L++)
								  {  HANSEN(-L-1,S,S+Q,E,&G1,&G2)
								      ;  G[L-1][LM+S][0]=G1
								      ;  G[L-1][LM+S][1]=G2
								  ;}
							  ;  W=S*DW+(S+Q)*(1+DM)+M*(DOW-CBZ/N)
							  ; if  (fabs(S+Q-M*CBZ/N)>5.0E-2)
							     {    N7=(S*CI-M)/(SI*H2)
								;  N8=1/W
								;  N9=3*(S+Q)*N8
								;  for(P=0;P<=11;P++)   DKP[P]=0.0
								;  if  (M>=2) L=M; else L=3
								; if  (L<abs(S))  L=abs(S)
								; if ((L-S) % 2) L++
								; while (((M==0) && (L<=LM)) || ((M!=0) && (L<=LT)))
								   {  J1=L*(L+1) -6+2*M
								      ; if  ((L-M)%2)
									     { CLM=-HK1[J1+1]
									      ;SLM=HK1[J1]
									     ;}
									     else
									      { CLM=HK1[J1]
									     ;  SLM=HK1[J1+1]
									     ;}

								      ;  P=(L-S) / 2
								      ;  J1=(L*(L+1)*(2*L+1)) / 3+2*((L+1)*M+P)
								      ;  FF1=FF[J1]
								      ;  FF2=FF[J1+1]
								      ; if  (CQ==0)
									 { G1=G[L-1][LM+S][0]
									    ;  G2=G[L-1][LM+S][1]
									     ;}
									else { G1=G[L-1][LM-S][0]
									    ;  G2=G[L-1][LM-S][1]
									;}
								      ;  N10=MCM[L-1]*FF1*G1
								      ;  N11=MCM[L-1]*FF1*G2
								      ;  N12=MCM[L-1]*FF2*G1
								      ;  DKP[0]+=N10*CLM
								      ;  DKP[1]+=N10*SLM
								      ;  DKP[7]+=N12*CLM
								      ;  DKP[6]-=N12*SLM
								      ;  DKP[9]+=N11*CLM
								      ;  DKP[8]-=N11*SLM
								      ;  DKP[11]+=N10*(2*(L+1)-N9)*CLM
								      ;  DKP[10]-=N10*(2*(L+1)-N9)*SLM
								      ;  L+=2
								   ;}

								;  DKP[0]*=N8
								;  DKP[1]*=N8
								;  DKP[2]=DKP[0]*N6/E
								;  DKP[3]=DKP[1]*N6/E
								;  DKP[4]=DKP[0]*N7
								;  DKP[5]=DKP[1]*N7
								;  DKP[0]*=(S+Q)
								;  DKP[1]*=(S+Q)
								;  DKP[10]=(DKP[10]+DKP[8]*N4+DKP[6]*N2)*N8
								;  DKP[11]=(DKP[11]+DKP[9]*N4+DKP[7]*N2)*N8
								;  DKP[8]=(DKP[8]*H2+DKP[6]*N3)*N8
								;  DKP[9]=(DKP[9]*H2+DKP[7]*N3)*N8
								;  DKP[6]*=N1*N8
								;  DKP[7]*=N1*N8
								; if    (   ((fabs(DKP[0])<EPSKP)
									&& (fabs(DKP[1])<EPSKP)
									&& (fabs(DKP[2])<EPSKP)
									&& (fabs(DKP[3])<EPSKP)
									&& (fabs(DKP[4])<EPSKP)
									&& (fabs(DKP[5])<EPSKP)
									&& (fabs(DKP[6])<EPSKP)
									&& (fabs(DKP[7])<EPSKP)
									&& (fabs(DKP[8])<EPSKP)
									&& (fabs(DKP[9])<EPSKP)
									&& (fabs(DKP[10])<EPSKP)
									&& (fabs(DKP[11])<EPSKP)) )
								   { K0=(S+Q)*CA+S*AP+M*(DBY-SN)
								      ;  SK=sin(K0)
								      ;  CK=cos(K0)
								      ;  E1[0]+=2*A*(DKP[0]*CK+DKP[1]*SK)
									 ;}
								  else {
								; fwrite(&S,sizeof(S),1,MULTIINDEX)
								; fwrite(&Q,sizeof(Q),1,MULTIINDEX)
								; fwrite(&M,sizeof(M),1,MULTIINDEX)
								; fwrite(DKP,sizeof(DKP),1,MAEIKP)
								; N_korper++
								       ;}
							    ;}                              /*   { % W<1.0E-2}  */
						      ;}                                    /*  { %  M}         */
						;}                                          /*   {%  S  }       */
					  ;}                                                /*  { %  ЦQ}        */
				    ;}
#if PECHAT
;printf("\n N_korper=%d \n ", N_korper);
#endif
				fseek(MULTIINDEX,2*sizeof(S),SEEK_SET)
			      ; fwrite(& N_korper,sizeof( N_korper),1,MULTIINDEX)
			      ;}      /*{ %K ПPOB A,E,I}*/
				  ;  for(J=0;J<=5;J++)   EK[J]=0.0
				  ;  fseek(MAEIKP,3*sizeof(AX),SEEK_SET)
				  ;  fseek(MULTIINDEX,2*sizeof(S),SEEK_SET)
				  ;  fread(&N_korper,sizeof(N_korper),1,MULTIINDEX)
				  ;  for(L=0;L<N_korper;L++)
				     {
				      if(fread(&S,sizeof(S),1,MULTIINDEX)&&
				      fread(&Q,sizeof(Q),1,MULTIINDEX)&&
				      fread(&M,sizeof(M),1,MULTIINDEX)&&
				      fread(DKP,sizeof(DKP),1,MAEIKP))
				      {
				      ;  K0=(S+Q)*CA+S*AP+M*(DBY-SN)
				      ;  SK=sin(K0)
				      ;  CK=cos(K0)
				      ;  for(J=0;J<=5;J++) EK[J]+=DKP[2*J]*CK+DKP[2*J+1]*SK

/*
#if PECHAT
;printf("\nS=%d Q=%d M=%d SK=%11.8f CK=%11.8f \n DKP(0-05)=",S,Q,M,SK,CK);
;  for(J=0;J<=5;J++)  printf("%10.6e ",DKP[J]);
printf("\n DKP(6-11)=");
;  for(J=0;J<=5;J++)  printf("%10.6e ",DKP[J+6]);
N_korper++
#endif
*/
				      ;}
				      else return(4);
				    }
#if PECHAT
;printf("\nN_korper=%d \n ",N_korper)
#endif

			    ;  E1[0]+=2.0*A*EK[0]
			    ;  E1[1]+=(K6*EK[1]-K7*EK[4])
			    ;  E1[2]+=(K7*EK[1]+K6*EK[4])
			    ;  E1[3]+=(0.5*CI2*K8*EK[2]-K9*EK[3])
			    ;  E1[4]+=(0.5*CI2*K9*EK[2]+K8*EK[3])
			    ;  E1[5]+=EK[5]
			      ;}/*{ПРОВЕРКА LM=2 && LT=0 }*/

#if BOKO
			;  if (PLUNASUN)
			 {

	      /*	       ; if ((T<THAX) || (T>TKOH))
				       APPROK(0.25,T,TK,&THAX,&TKOH);
				 else
	       */
				 {
				    X=T-THAX
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
				 }
			    ;  LK=LLUN
			    ;  AAL=sqr(A/ALUN)
			    ;  AAC=sqr(A/AC)

			    ; for  (L= 2;L<=LK;L++)
			       {
				  ;  AAL*=A/ALUN
				  ;  AAC*=A/AC
				  ;  CLM=ML*AAL/(2*L+1)
				  ;  SLM=MC*AAC/(2*L+1)
				  ;  QMM=QMAX
				  ; for (M=0;M<=L;M++)
				       {  J1=(L*(L+1))/2-3+M
					;  L0=2*J1
					;  L1=L0+1
					;  J=(((L-M)%2) ? 0:1)
					; for (P=0;P<=L;P++)
					     { for (Q=-QMM;Q<=QMM;Q++)
						{ if (Q==-QMM)
						    {  J1=(L*(L+1)*(2*L+1)) / 3+2*(L+1)*M+2*P
						    ;  FF1=FF[J1]
						    ;  FF2=FF[J1+1];
						    }
						       S=L-2*P+Q
						    ; if  (L-2*P<0)
							   { G1=GK[L-P][-Q+QMAX][0]
							    ;G2=GK[L-P][-Q+QMAX][1]
							  ;}
                              else
                               { if (M==0)
							    {  HANSEN(L,S-Q,S,E,&G1,&G2)
								;  GK[P][Q+QMAX][0]=G1
								;  GK[P][Q+QMAX][1]=G2
								;}
						        else
                                {
							     G1=GK[P][Q+QMAX][0]
							    ;G2=GK[P][Q+QMAX][1]
						        ;}
                               }
						    ; if  (S!=0)
						       { K1=FF1*G1*CLM
							  ;  K2=FF1*G2*CLM
							  ;  K3=FF2*G1*CLM
							  ;  N1=FF1*G1*SLM
							  ;  N2=FF1*G2*SLM
							  ;  N3=FF2*G1*SLM ;
 	      ; if ((fabs(K1/E)>=1E-12) || (fabs(K1/SI)>=1E-12) || (fabs(K2)>=1E-8) || (fabs(K3)>=1E-8)
                ||(fabs(N1/E)>=1E-12) || (fabs(N1/SI)>=1E-12) || (fabs(N2)>=1E-8) || (fabs(N3)>=1E-8 ))
							     { K0=(S-Q)*AP+S*CA+M*DBY
                                ; dS=(S-Q)*DW+S*(1+DM)+M*DOW
//                                ; dS=S;
								; SK=sin(K0)
								; CK=cos(K0)
								; U=HL[L0]+DL[L1]/(N*dS)-D2L[L0]/(dS*dS*N*N)
								; V=HL[L1]-DL[L0]/(N*dS)-D2L[L1]/(dS*dS*N*N)
								; if  (L<=LSUN)
								   { W=HC[L0]+DC[L1]/(N*dS)-D2C[L0]/(dS*dS*N*N)
								    ;Z=HC[L1]-DC[L0]/(N*dS)-D2C[L1]/(dS*dS*N*N)
								  ;}
								  else { W=0; Z=0;}
								;  K4=U*(J*CK+(1-J)*SK)+V*(J*SK+(J-1)*CK)
								;  K5=U*(-J*SK+(1-J)*CK)+V*(J*CK+(1-J)*SK)
								;  N4=W*(J*CK+(1-J)*SK)+Z*(J*SK+(J-1)*CK)
								;  N5=W*(-J*SK+(1-J)*CK)+Z*(J*CK+(1-J)*SK)
								;  K0=K1*K5+N1*N5
								;  K1=K1*K4+N1*N4
								;  K2=K2*K5+N2*N5
								;  K3=K3*K5+N3*N5
//								;  E1[0]+=2*A*K1
								;  E1[0]+=2*A*K1*S/dS
								;  if (abs(Q)<=8)
//								   {     F1=H2*(H2+(2.0*P-L)/dS)/E*K1
								   {     F1=H2*(H2*S/dS+(2.0*P-L)/dS)/E*K1
								      ;  F2=-(H2*K2+E*TI2*K3/H2)/dS
								      ;  E1[1]+=(K6*F1-K7*F2)
								      ;  E1[2]+=(K7*F1+K6*F2)
								      ;  F1=K1*((L-2*P)*CI-M)/(SI*H2*dS*2)
								      ;  F2=-K3/(2*H2*CI2*dS)
								      ;  E1[3]+=(CI2*K8*F1-K9*F2)
								      ;  E1[4]+=(CI2*K9*F1+K8*F2)
//								      ;  E1[5]-=((H2+E*E-1)/E*K2+TI2/H2*K3-(2*L+3)*K0)/dS                                      
								      ;  E1[5]-=((H2+E*E-1)/E*K2+TI2/H2*K3-(2*L+3*S/dS)*K0)/dS
									;}
							      ;}
						      ;}
						;}
					  ;}
				    ;}
			      ;}
			;}

#endif
			free(FF)
			;      }
		  return(0)
		  ;}


