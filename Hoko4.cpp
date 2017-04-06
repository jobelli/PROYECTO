#include "hoko.h"
#include <time.h>
#include <string.h>
#define sqr(x) ((x)*(x))
#define PECHAT 0
#define TRACK 0
/*********************** Только для графики***************/
#if TRACK
int bd_prox,i_bd_prox;
void track(int NT, double T);
#endif
double xp[3],nx,skor,HR;
long NO,Nk;
/*********************** Конец только для графики***************/

extern double F135
       , FT
       , KP
       , far*AT1
       ;

extern double KB0,KOTP,DT0,S0,HP,ALTITUDE,LONGITUDE,DTREENTRY;
#if HOKO
extern int const MF0[6];
extern double   AT[450];
#endif
extern int PLUNASUN,PATMOSFERA,PSUNDAB;
void   JVC_DLL x_k(double h[],double e[], int pech);
//void x_k(double h[],double e[], int pech);
void   JVC_DLL km_x(double E0[],int prz,double EK[],int pech);
//void km_x(double E0[],int prz,double EK[],int pech);   /* у Васи GNSK */
int    JVC_DLL FORCE_CHIS(double *x,double *V,double tm,double *F);
int    JVC_DLL chisint(double *x,double *V,double tI,double tF,double xL,int LL,int NV,int NI,int *NF,int *NS,int klass,int NoR,int pr);
//int  FORCE_CHIS(double *x,double *V,double tm,double *F);
//int chisint(double *x,double *V,double tI,double tF,double xL,int LL,int NV,int NI,int *NF,int *NS,int klass,int NoR,int pr);
int XAXA=0,NZ,MTS,MTS1,NTS,Nya,MI1,Nm=LMAX+2;
double              ym,ReFL,SS,SWW,
                         far *C,
			 far *SS2,
			 far *Wm,
			 far *Vm,
 		         far *H_rada,
			 far *W_rada,
			 far *U_rada,
			 far *C_rada,
			 far *D_rada,
			 far *R_rada,
			 far *XI_rada,
			 far*F1_rada,
			 far*FJ_rada,
			 far*Y_rada,
			 far*Z_rada
			 ;
#ifdef HOKODLL
 extern "C" int _export _stdcall PROGNOZCHICL(int NTOCH,struct bxprog *BX,struct bixprog *BIX,double TKNK[])
#else
 int  PROGNOZCHICL(int NTOCH,struct bxprog far*BX,struct bixprog far*BIX,double TKNK[])
#endif
  {
/*   выход PROGNOZ=
      0 -все нормально
      1 -неправильно заданы исходные данные
      2 -ошибка открытия файлов
      3 -недостаточно памяти для решения задачи
      4 -ошибка чтения данных из файлов
      777 -сгорел
*/




int
        PRMODEL=0
      , PZADACHI=0
      , NT=0              /*номер точки прогноза*/
      , m,k,L,N=0
      , LL=7
      ,	NV=3
      , NI=2
      , NF=0
      , NS=0
      , N_iterac=3
      , klass=2
      , NoR=15
      , pr_prod_int=0
      , J
      , JD
      , L13
      , L5
      , J1
      , k1
      , k2
      , xL2
      , xL3
	;
double  xL=0.003
     ,	A0
     ,  BETTA
     ,  P
     ,  I0
     ,  EKC
     ,  L0
     ,  H0
     ,  U0
     ,  arg_per
     ,  EA
     ,  MA
     ,  DBY
     ,  arg_shiroti
     ,  arg_perigee
     ,  H
     ,  IA
     ,  SE
     ,  CE
     ,  CIA
     ,  DM
     ,  DA=0
     ,  DKB=0
     ,  DDT=0
     ,  DT=0
     ,  Dt=0
     ,  T=0
     ,  POCK
     ,  RAD
     ,  PER
     ,  TK
     ,  DTK
     ,  PADEHIE=0
     ,  LIFETIME=0
     ,  DTBAKB=0
     ,  DTGRAV=0
     ,  E0[7]
     ,  ER[6]
     ,  E1[7]
     ,  EK[6]
     ,  X[3]
     ,  V[3]
     ,  R1
     ,  Se
     ,  R2

	;

 time_t time1, time2;

char   *shablon,*shabout;
 /************************Ввод коэффициентов гравитационного поля земли*****************************/
#if !defined(GEM_10)
	 if (! GEM_input())

	 {perror("\nError of GEM coefficients input\n");close_prognoz();return(2);}
#endif
 /***************************************************************************/

     /*********BX-СТРУКТУРА bxprog,ГДЕ СФОРМИРОВАНЫ ВХОДНЫЕ ДАННЫЕ**********/
     /**********************************************************/
			if((C=(double far*) calloc(((LMAX+2)*(LMAX+3))/2,sizeof(double)))==NULL||
			  (SS2=(double far*) calloc(((LMAX+2)*(LMAX+3))/2,sizeof(double)))==NULL||
			  (Wm=(double far*) calloc(((LMAX+2)*(LMAX+3))/2+1,sizeof(double)))==NULL||
                          (Vm=(double far*) calloc(((LMAX+2)*(LMAX+3))/2+1,sizeof(double)))==NULL||
                          (H_rada=(double far*) calloc(13,sizeof(double)))==NULL||
			  (W_rada=(double far*) calloc(13,sizeof(double)))==NULL||
			  (U_rada=(double far*) calloc(13,sizeof(double)))==NULL||
                          (C_rada=(double far*) calloc(78,sizeof(double)))==NULL||
			  (D_rada=(double far*) calloc(78,sizeof(double)))==NULL||
			  (R_rada=(double far*) calloc(78,sizeof(double)))==NULL||
                          (XI_rada=(double far*) calloc(78,sizeof(double)))==NULL||
                          (F1_rada=(double far*) calloc(NV,sizeof(double)))==NULL||
			  (FJ_rada=(double far*) calloc(NV,sizeof(double)))==NULL||
			  (Y_rada=(double far*) calloc(NV,sizeof(double)))==NULL||
			  (Z_rada=(double far*) calloc(NV,sizeof(double)))==NULL)
     {printf("\n НЕТ ПАМЯТИ ДЛЯ ДИНАМИЧЕСКOГО МАССИВА В ПРОЦЕДУРЕ chisint");
     exit(1);
     }


     /*********BX-СТРУКТУРА bxprog,ГДЕ СФОРМИРОВАНЫ ВХОДНЫЕ ДАННЫЕ**********/
     /**********************************************************/
	NO=BX->NO;
	Nk=BX->NH;
	DT0=BX->DT;
	A0=E0[0]=BX->A;
	I0=E0[2]=BX->I;
	DBY=E0[4]=BX->DBY;
	L0=BX->L;
	H0=BX->H;
	U0=E0[5]=BX->U;
	KB0=BX->KB;
	F135=BX->F135;
	FT=BX->FT;
	KP=BX->KP;
	PRMODEL=BX->PRMODEL;
	POCK=BX->POCK;
	PZADACHI=BX->PZADACHI;
shablon="\n\nВХОД\n%d %d %lf %lf %lf %lf %lf %lf %lf %f %f %f %f %d %d %d";
#if PECHAT
printf( shablon,BX->NO,BX->NH,DT0,A0,I0,DBY,L0,H0,U0,KB0,F135,FT,KP,PRMODEL,POCK,PZADACHI);
#endif
    /*******************КОНЕЦ СЧИТЫВАНИЯ ВХОДНЫХ ДАННЫХ*************/
		         JD=DT0  /*12053-ЧИСЛО СУТОК ОТ 31.12.57 ДО 31.12.90*/
 /**************************************************************************/
		      ;  BETTA=L0*L0+H0*H0
		      ;  E0[1]=EKC=sqrt(BETTA)
		      ;  E0[3]=(L0==0&&H0==0)? 0 : atan2(H0,L0)
                      ;  arg_perigee=E0[3]
           	      ;  BETTA=1-BETTA

	;if (A0<1880 && A0>85)
	{ DTBAKB=1;
	  N_iterac=1;
	  DT=KB0/1440;
	  PER=A0/1440;
          KB0=0;
	  A0=pow(GME*sqr(PER*86400.0/PI2),1.0/3.0);
	  P=A0*BETTA;
	  E0[0]=A0*=1-C20*sqr(RZ/P)*((2-2.5*sqr(sin(I0)))*sqrt(pow(BETTA,3))/sqr(1+L0)+
		pow((1+L0),3)/BETTA);
	}
 /***************************************************************************/
 /*************************ПРОВЕРКА НА ДОСТОВЕРНОСТЬ ВХОДА*****************/
      ;	 BIX->HP=A0*BETTA/(1+EKC*cos(U0-arg_perigee))-RZ*(1-ECJ*sqr(sin(I0)*sin(U0)));
      ;  if (BIX->HP<0.0)
	 {
	 printf("\nВХОДНЫЕ ДАННЫЕ ВЫХОДЯТ ЗА ПРЕДЕЛЫ ДОПУСТИМОГО\n")
          ;  return(776);
	 }
 /**************************************************************************/
	BETTA=sqrt(BETTA)
       ; if  (I0<1.0E-8)  I0=1.0E-8


 /********************ОПРЕДЕЛЕНИЕ МОДЕЛИ ДВИЖЕНИЯ**********************************/
		       ; switch(PRMODEL)
			 {
			case 0:   {    PATMOSFERA=1
			      ;  NZ=NTS=MTS=LMAX
			      ;  if  (A0<10000)
				   { PSUNDAB=0;PLUNASUN=0  ;}
                 else { PSUNDAB=0;PLUNASUN=1;}
			     break
			     ;};
            case 1:   {    PATMOSFERA=0
			      ;  NZ=2;NTS=MTS=0
			      ;  PSUNDAB=0;PLUNASUN=0;
			     break
			     ;};
			 default:
			     {
                        PLUNASUN=0;PATMOSFERA=1;PSUNDAB=0;
                     ;  NZ=MTS=NTS=MIN(PRMODEL,LMAX)

			     ;}
			 }
 /****************КОРРЕКТИРОВКА МОДЕЛИ ДВИЖЕНИЯ************/
		    ; if  (KOTP==0)  PSUNDAB=0
            ; BIX->HP=A0*(1-EKC)-RZ;
		    ; if  ((KB0==0 && DT==0)||BIX->HP>1550)   PATMOSFERA=0 ;
 /**********************************************************/
 /****************************************************************************/
		      ;  for(m=0 ;m<=LMAX; m++)
			 {     for(L=m ;L<=LMAX; L++)
			       {     if    (L>=2)
				     {     k=L*(L+1)-6+2*m;
					   C[N]=HK1[k]*1.0e6;
					   SS2[N]=HK1[k+1]*1.0e6;
				     }
                                     else C[N]=SS2[N]=0;
				     N++;
			       }
                             N++;
			 }
			 C[0]=1.0e6 ;
			 NTS++;
			 NZ++;
			 MTS++;
			 ym=28899.4592;
			 ReFL=1.44*4.65e-4*8.64*8.64*KOTP/RZ;
			 MTS1=MTS+1;
			 Nya=NTS;
			 MI1=MTS1;
			 if (NTS<=NZ)  Nya=NZ;
			 Se=GME*86.4*86.4/(RZ*RZ*RZ);
			 SWW=Se*ML*1.0e6;
			 SS=Se*MC*1.0e6;
			 for(L=1 ;L<Nm; L++)
			 {     xL2=2*L-1;
			       R1=sqrt(xL2);
			       C[L-1]*=R1*Se*ym;
			 }
			 for(m=2 ;m<Nm; m++)
			 {     L5=Nm-m;
			       L13=(m-1)*Nm-((m-1)*m)/2+m-1;
			       for(k=1 ;k<=L5; k++)
			       {     xL3=4*(m+k)-6;
				     R1=sqrt(xL3);
				     J1=L13+k-1;
				     C[J1]*=R1*ym;
				     SS2[J1]*=R1*ym;
				     k1=k+2*m-3;
				     for(k2=k ;k2<= k1; k2++)
				     {     R2=sqrt(k2);
					   C[J1]/=R2;
					   SS2[J1]/=R2;
				     }
				     C[J1]*=Se;
				     SS2[J1]*=Se;
			       }
			 }
			 ym=1/ym;
 /****************************************************************************/

  /*****************РАСЧЕТ СРЕДНЕЙ АНОМАЛИИ********************/
		      ;  IA=U0-arg_perigee
		      ;  CIA=cos(IA)
		      ;  SE=BETTA*sin(IA)
		      ;  CE=EKC+CIA
		      ;  EA=atan2(SE,CE)
		      ;  MA=EA-EKC*SE/(1+EKC*CIA)
		      ;  MA=fmod(MA,PI2)
 /*************************************************************/

 /****************ВЫЧИСЛЕНИЕ ЗВЕЗДНОГО ВРЕМЕНИ-S0 НА НАЧАЛЬНЫЙ МОМЕНТ*********/
		      ;  S0=ST(JD,DT0-JD)
 /****************************************************************************/

		       ;  DM=0


  ; if(!DTBAKB)
  {
  ; PER=PI2*E0[0]/86400.0*sqrt(E0[0]/GME)
  ; P=(1-EKC*EKC)
  ; PER*=1+1.5*C20*sqr(RZ/(P*E0[0]))*sqrt(P)*((2-2.5*sqr(sin(E0[2])))*P/sqr(1+E0[1]*cos(E0[4]))+
	(1-1.5*sqr(sin(E0[2]))));

   }

                 km_x(E0,1,EK,0);   /* у Васи GNSK */
                 for(J=0;J<3;J++){X[J]=EK[J]/RZ;V[J]=EK[J+3]/RZ*86400.0;}
                 for(J=0;J<6;J++) ER[J]=E0[J];

 /************НАЧАЛО ЦИКЛА ПО MОМЕНТАМ ПРОГНОЗА****************************************************************/




    do
   {
   /*начало цикла по моментам прогноза при прогнозе на
		     время и номер витка */
/*% OПPEДEЛEHИE ШAГA H,BPEMEH HAЧAЛA И KOHЦA(T И TK)*/
	   if (PZADACHI==3)   NTOCH=1000;
		     if (DTBAKB==0)
		    {
		      if (PZADACHI==0)
			 {

			 DTK=TKNK[NT]
		      ;  TK=DTK-DT0-T
              ;  if(fabs(TK)>300*PER)
			  {
			 printf("\nНЕПРАВИЛЬНО ОПРЕДЕЛЕНЫ МОМЕНТЫ ПРОГНОЗИРОВАНИЯ(>10 сут)")
		      ;  return(1)
			 ;}
			  }
			else if (PZADACHI==1)
			 {
			   BIX->bix.NH=(long)TKNK[NT]
			 ; TK=PER*(BIX->bix.NH-BX->NH)-T
			 ; DTK=DT0+TK+T
             ;  if(fabs(TK)>1000*PER)
			  {
			 printf("\nНЕПРАВИЛЬНО ОПРЕДЕЛЕНЫ МОМЕНТЫ ПРОГНОЗИРОВАНИЯ(>10 сут)")
		      ;  return(1)
			  ;}
			  }
			else if (PZADACHI==2)
			      {
			       CIA=cos(arg_perigee)
			    ;  SE=-BETTA*sin(arg_perigee)
			    ;  CE=EKC+CIA
			    ;  EA=atan(SE/CE)
			    ;  if (CE<0)   EA+=PI
			    ;  DM=EA-EKC*SE/(1+EKC*CIA)
			    ;  DM=MA-DM
			    ;  DM=fmod(DM,PI2)
			    ;  if (DM<0) DM+=PI2
			    ;  TK=-PER/PI2*DM
			    ;  DTK=DT0+TK
			     ;}
		       else if (PZADACHI==3)
			 {
                 NTOCH=100*(TKNK[1]-TKNK[0])/PER
              ;  DTK=DT0+0.01*PER*NT+0.00000001
		      ;  TK=DTK-DT0-T
              ;}
		      }
		      else
		      {
			   TK=-PER
			 ; DTK=DT0+TK+T
			 ; NT=-1
			 ;  if(fabs(TK)>100)
			  {
			 printf("\nНЕПРАВИЛЬНО ОПРЕДЕЛЕНЫ МОМЕНТЫ ПРОГНОЗИРОВАНИЯ")
		      ;  return(1)
			 ;}
		      }

 /************НАЧАЛО ЦИКЛА ПО ВРЕМЕНИ****************************************************************/
	 ;  do     /*цикл для пересчета*/
	     {
            do    /*цикл для  выхода на два витка назад*/
	     {
            do    /*цикл для  выхода на U=0*/
	     {
#if PECHAT
		       ;  printf("\nВРЕМЯ СЧЕТА KYTTA : %f  sec\n",difftime(time2,time1));
		       ;  printf("\nPOSLE kytta\n")
		       ;  for(J=0;J<=5;J++) printf("%  lf",Q1[J]);
#endif

                      if(fabs(TK)>1.0e-9)
                      if((XAXA=chisint(X,V,T,T+TK,xL,LL,NV,NI,&NF,&NS,klass,NoR,pr_prod_int))) goto MREENTRY;
                      pr_prod_int=1;
                      T+=TK;
		      for(J=0;J<3;J++){BIX->X[J]=X[J]*RZ;BIX->X[J+3]=V[J]*RZ/86400.0;}
   		      x_k(BIX->X,E1,0)
			    ; if  ((PZADACHI!=0&&PZADACHI!=3) || DTBAKB)
				{
			       TK=-X[2]/V[2]
			    ;  DTK=DT0+T+TK
				;} else
				    TK=0;
	   }
	 while( fabs(TK)>1E-9);
            if (DTBAKB && N_iterac==1) {Dt=T ;TK=T;}
            if (DTBAKB && N_iterac==2) PADEHIE=-2*Dt+T ;
            N_iterac++;
          } /* конец цикла по прогнозам на два витка назад*/
         while( N_iterac<3);
         if(DTBAKB==1)

	 {
	   Dt+=PER;
     	   TK=-PER ;
	   T=0;
	 if(PATMOSFERA && DT<0 )
	 {
           if(DTGRAV==0)
            {DTGRAV=PADEHIE;DKB=0.001;}
           else
            {PADEHIE-=DTGRAV;DDT=DT-PADEHIE;DKB=DDT*KB0/PADEHIE;}
         } else {DDT=DKB=DT=0;};
	  KB0+=DKB;
	  if(KB0<0) KB0=0.0001;
	  if(KB0>200) KB0=200;
	  DA=2*ER[0]*Dt/(3*PER);
	  ER[0]+=DA;
          N_iterac=1;
          pr_prod_int=0;
          km_x(ER,1,EK,0);   /* у Васи GNSK */
          for(J=0;J<3;J++){X[J]=EK[J]/RZ;V[J]=EK[J+3]/RZ*86400.0;}
	  for (J=0;J<=5;J++)  E0[J]=ER[J];
	  }
	  }
	 while( fabs(Dt)>1E-8 || (fabs(DDT)>0.01*fabs(DT) && fabs(DDT)>1E-8));
	 if(DTBAKB==1 && PZADACHI<4) {DTBAKB=0;Nk=BX->NH;NT++;continue;}
         DTK=DT0+T
  ; if (((PZADACHI==0) && (fabs(DTK-DT0)<1.0E-12)) ||((PZADACHI==1) && (BIX->bix.NH==BX->NH) && U0==0) || ((PZADACHI==2) && (U0==0))||(PZADACHI>3))
	       {
	      E1[0]=A0
	    ; E1[1]=EKC
	    ; E1[2]=I0
	    ; E1[3]=arg_perigee
	    ; E1[4]=DBY
	    ; E1[5]=MA
		  ;}
		      ;  BIX->bix.A=E1[0]
		      ;  BIX->bix.I=E1[2]
		      ;  BIX->E=E1[1]
		      ;  BIX->W=E1[3]
		      ;  if (BIX->W<0) BIX->W+=PI2
		      ;  BIX->bix.L=BIX->E*cos(BIX->W)
		      ;  BIX->bix.H=BIX->E*sin(BIX->W)
		      ;  BIX->bix.DBY=E1[4]
		      ;  BIX->bix.DBY=fmod(BIX->bix.DBY,PI2)
		      ;  if (BIX->bix.DBY<0) BIX->bix.DBY+=PI2
		      ;  BIX->bix.U=E1[6]
              ;  if (BIX->bix.U<0) BIX->bix.U+=PI2
		      ;  KBL(E1,BIX->Z)
              ;  BIX->Y[0]=E1[0]
		      ;  BIX->Y[1]=E1[2]
		      ;  BIX->Y[2]=E1[4]
		      ;  BIX->Y[3]=E1[1]*cos(E1[3])
		      ;  BIX->Y[4]=E1[1]*sin(E1[3])
		      ;  BIX->Y[5]=E1[6]
		      ; if ((PZADACHI==1) || (PZADACHI==2))  BIX->bix.U=0
;
    if (!DTBAKB)
       PER=PI2/86400.0*E1[0]*sqrt(E1[0]/GME)*(1+1.5*C20*sqr(RZ/E1[0]/(1-EKC*EKC))*((2-2.5*sqr(sin(I0)))*
       sqrt(pow((1-EKC*EKC),3))/sqr(1+L0)+pow((1+L0),3)/(1-EKC*EKC)));
			if  (PZADACHI==0)
			  { BIX->bix.NH=BX->NH+(T/PER)
			    ;  BIX->bix.NH=BX->NH+(T/(PER+0.5*PADEHIE*(BIX->bix.NH-BX->NH)))
			      ;}
		       /*   ; if  ((TK>0) && (U0<BX[8]))
				BIX->bix.NH=BIX->bix.NH+1
			    ; if  (TK<0) && (U0>BX[8])
				BIX->bix.NH=BIX->bix.NH-1
			      else DTK=DT0+TK   */
			    ; if  (PZADACHI==2||PZADACHI>3)
				BIX->bix.NH=BX->NH
/***************ЗАПИСЬ РЕЗУЛЬТАТА***************************************/
          ;RAD=sqrt(BIX->X[0]*BIX->X[0]+BIX->X[1]*BIX->X[1]+BIX->X[2]*BIX->X[2])
	  ;HR=RAD-RZ+ECJ*BIX->X[2]/RAD;
	  ;HP=BIX->bix.A*(1+BIX->E)-RZ*(1.0-ECJ*sin(BIX->bix.I)*sin(BIX->W))
          ;BIX->HA=BIX->bix.A*(1-BIX->E)-RZ*(1.0-ECJ*sin(BIX->bix.I)*sin(BIX->W+0.5*PI2))
	  ;ALTITUDE=asin(BIX->X[2]/RAD)*RADIAN
	  ;LONGITUDE=(BIX->bix.DBY-S0-CBZ*T+atan2(sin(BIX->bix.U)*cos(BIX->bix.I),cos(BIX->bix.U)))*RADIAN
	  ;LONGITUDE=fmod(LONGITUDE,360)
	  ;if (LONGITUDE<-180) LONGITUDE+=360
     ;if (LONGITUDE>180) LONGITUDE-=360
	  ;
;MREENTRY:
	BIX->bix.NO=(XAXA) ? 777 : BX->NO;
	BIX->bix.NH=Nk;
	BIX->bix.DT=(XAXA) ? DTREENTRY : DTK;
	BIX->bix.KB=KB0;
	BIX->bix.F135=F135;
	BIX->bix.FT=FT;
	BIX->bix.KP=KP;
	BIX->bix.PRMODEL=PRMODEL;
	BIX->bix.POCK=POCK;
	BIX->bix.PZADACHI=PZADACHI;
	BIX->PER=PER*1440;
	BIX->DT=PADEHIE*1440;
	BIX->LIFET=DTREENTRY;
	BIX->RAD=RAD;
	BIX->ALTITUDE= ALTITUDE;
	BIX->LONGITUDE=LONGITUDE;
if(LONGITUDE>180) LONGITUDE-=360;

#if TRACK
/************************ТРАССА************************/
if(PZADACHI==3)
{   xp[0]=RAD*cos(LONGITUDE/RADIAN)*cos(ALTITUDE/RADIAN);
    xp[1]=RAD*sin(LONGITUDE/RADIAN)*cos(ALTITUDE/RADIAN);
    xp[2]=RAD*sin(ALTITUDE/RADIAN);
   track(NT,DT0+T);
}
/************************КОНЕЦ ТРАССЫ************************/
#endif

	if (XAXA)   return(XAXA);
	NT++;
	if ((PZADACHI<2) && NT<NTOCH) BIX++;
	  }
    /*конец цикла по моментам прогноза при прогнозе на
		     время и номер витка */
while  (PZADACHI<4 && NT<NTOCH);
			 free(C);
			 free(SS2);
			 free(Wm);
			 free(Vm);
			 free(H_rada);
			 free(W_rada);
			 free(U_rada);
			 free(C_rada);
			 free(D_rada);
			 free(R_rada);
			 free(XI_rada);
			 free(F1_rada);
			 free(FJ_rada);
			 free(Y_rada);
			 free(Z_rada)
			 ;

#if TRACK
/************************ТРАССА************************/
if(PZADACHI==3)  {d4close();while(!kbhit());closegraph();}
/************************КОНЕЦ ТРАССЫ************************/
#endif

  return(XAXA);

  ;}

/* ***************************************      	 */
		    double intst(int x,double a)   	        	        	        	     /*< */
/* ЦEЛCT: A B CTEПEHИ X (X-ЦEЛOE (ЛЮБOГO ЗHAKA)   (A**X) */
		   {
		   int I;
		   double S=1.0,b=(x>=0) ? a :(a==0)? 0:1.0/a ;
		   for(I=1 ;I<=abs(x); I++) S*=b;
		   return(S);
		   }
/* ******************** 				 */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      */


//  int chisint(double *x,double *V,double tI,double tF,double xL,int LL,int NV,int NI,int *NF,int *NS,int klass,int NoR,int pr)
 int    JVC_DLL chisint(double *x,double *V,double tI,double tF,double xL,int LL,int NV,int NI,int *NF,int *NS,int klass,int NoR,int pr)
/*  ,pr,F1_rada,FJ_rada,Y,Z,veii,vtii,vii,a,n,U,Wm,s,R,D,xI,FORCE_CHIS);*/
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       */
		   {
	     double      oNe=1.0,
			 ZeRo=0;
	     static
	     int
			 kD,
			 kD2,
			 ke,
			 kF,
			 Nper,
			 NSF,
			 NsL,
			 NPQ,
			 NeS;

	     static
	     double
			 PW,
			 SS_rada,
			 SR,
			 tDIF,
			 DIR=0,
			 t,
			 tm,
			 hsum,
			 W1,
                         B_rada[13][3],
                         BE_rada[13][3],
                         BT_rada[13][3],
                         Z0
			 ;
	     int         I,
			 J,
			 k,
			 L,
			 La,
			 Lv,
			 LD,
			 m,
			 JD,
			 Ls,
			 JDm,
			 NL,
			 J2,
			 Le,
			 N,
			 Ncount,
			 JR;
	     double      tR,
			 Sm,
			 tval,
			 t2,
			 ReS,
			 temr,
			 val,
			 bdubl,
			 S,
                         Q,
			 tr
                         ;
      const int NW[]={0,0,1,3,6,10,15,21,28,36,45,55,66,78},
		MCONST[]={1,13,24,34,43,51,58,64,69,73,76,78},
 	       NxI[]={2,3,4,5,6,7,8,9,10,11,12,13,3,6,10,15,21,28,36,45,55,66,78,4,10,20,35,56,84,
		      120,165,220,286,5,15,35,70,126,210,330,495,715,6,21,56,126,252,462,792,1287,7,28,84,210,462,
		      924,1716,8,36,120,330,792,1716,9,45,165,495,1287,10,55,220,715,11,66,286,12,78,13};
   const double HH[]={0.212340538239152,0.590533135559265,0.911412040487296,0.098535085798826426,
		      0.304535726646363905,0.562025189752613855,0.801986582126391827,0.960190142948531257,
		      0.056262560536922146,0.180240691736892364,0.352624717113169637,0.547153626330555383,
		      0.734210177215410531,0.885320946839095768,0.977520613561287501,0.036257812883209460,
		      0.118078978789998700,0.237176984814960385,0.381882765304705975,0.538029598918989065,
		      0.690332420072362182,0.823883343837004718,0.925612610290803955,0.985587590351123451,
		      0.025273620397520349,0.083041613447405146,0.169175100377181425,0.277796715109032074,
		      0.401502720232860816,0.531862386910415957,0.659991842085334811,0.777159392956162144,
		      0.875380774855556926,0.947964548872819447,0.989981719538319594,0.018610365010987851,
		      0.061475540899268987,0.126305178693310580,0.209842971726562514,0.307898998280398343,
		      0.415556035978659544,0.527415613995882274,0.637868602717761199,0.741376459294237483,
		      0.832748988608442268,0.907404775300997364,0.961601861260321649,0.992635348973910678};


		      if (pr==0)
	       {
			 Nper=0;
			 NSF=0;
                         Z0=x[2];
			 kD=(NoR-3)/2;
			 kD2=kD/2;
			 ke=kD+1;
			 kF=kD+2;
			 PW=oNe/(kD+3.0);
			 NsL=(klass==1) ? 1 : 0;
			 NPQ=(klass<2) ? 1:0;
			 NeS=(LL<0) ? 1 :0;
			 SR=(NV==1) ? 1.2 : 1.5;
			 tDIF=tF-tI;
		         if(tDIF==0) return(0); /*????????????*/
			 DIR=sign(tDIF);
		      if (NeS) xL=fabs(xL)*DIR;
			 klass=abs(klass);
			 La=kD2*kD2-2;
			 for(I=2; I<= kF; I++)
			{ La++;
			  H_rada[I-1]=HH[La];
			  W_rada[I-2]=oNe/(I+I*I*(klass-1.0));
			  U_rada[I-2]=I+1;
			}
			  for(k=0; k<NV; k++)
			{ if (NsL) V[k]=ZeRo;
			  for(L=0; L<=kD; L++)
                           BT_rada[L][k]=B_rada[L][k]=ZeRo;
			}
			  W1=oNe/klass;
			  for(J=0; J<=kD-1; J++)
			 { m=MCONST[J]-1;
                           JD=J+1;
			   for(L=JD; L<=kD; L++)
			  { XI_rada[m]=NxI[m]*W_rada[J]/W_rada[L];
				     m++;
			  }
			 }
			  C_rada[0]=-H_rada[1]*W_rada[0];
			  D_rada[0]=H_rada[1]/W_rada[1];
			  R_rada[0]=oNe/(H_rada[2]-H_rada[1]);
			  La=0;
			  Ls=0;
			  for(k=3; k<=ke; k++)
			{ Lv=La;
			  La=Ls+1;
			  Ls=NW[k]-1;
			  JD=Ls-La;
			  C_rada[La]=-H_rada[k-1]*C_rada[Lv];
			  C_rada[Ls]=(C_rada[La-1]/W_rada[JD-1]-H_rada[k-1])*W_rada[JD];
			  D_rada[La]=H_rada[1]*D_rada[Lv]*W_rada[k-2]/W_rada[k-1];
			  D_rada[Ls]=(D_rada[La-1]*W_rada[k-2]+H_rada[k-1])/W_rada[k-1];
			  R_rada[La]=oNe/(H_rada[k]-H_rada[1]);
			  R_rada[Ls]=oNe/(H_rada[k]-H_rada[k-1]);
			if (k!=3)
			 { for(L=4; L<=k; L++)
			  {LD=La+L-3;
			   Le=Lv+L-4;
			   JDm=LD-La;
			   C_rada[LD]=W_rada[JDm]*C_rada[Le]/W_rada[JDm-1]-H_rada[k-1]*C_rada[Le+1];
			   D_rada[LD]=(D_rada[Le]+H_rada[L-2]*D_rada[Le+1])*W_rada[k-2]/W_rada[k-1];
			   R_rada[LD]=oNe/(H_rada[k]-H_rada[L-2]);
			  }
			 }
			}
			   SS_rada=pow(10.0,-LL);
			   NL=NI+30;
			   tr=((NoR/11.0)*pow(0.5,0.4*LL))*DIR/2;
			  if (NeS)  tr=xL;
			  if (tr/tDIF>0.5)  tr=0.5*tDIF;
			   *NF=0;
			   Ncount=0;
		   m1:     *NS=0;
			   tm=tI;
			   Sm=1e4;
             if((XAXA=FORCE_CHIS(x,V,tm,F1_rada))) return(XAXA);
			   *NF=*NF+1;
		   m2:
			   for(k=0; k<NV; k++)
		       {   BE_rada[kD][k]=B_rada[kD][k]/W_rada[kD];
			   for(J=0 ;J<kD; J++)
		        {  JD=J+1;
                           BE_rada[J][k]=B_rada[J][k]/W_rada[J];
        		   for(L=JD; L<=kD; L++)
			 { N=NW[L]+J;
			   BE_rada[J][k]+=D_rada[N]*B_rada[L][k];
			 }
			}
		       }
			   t=tr;
			   tval=fabs(t);
			   t2=intst(klass,t);
			   for(m=1; m<=NL; m++)
		     {     J2=1;
			   for(J=1; J<=ke; J++)
		      {     JD=J-1;
			    La=NW[JD]-1;
			    JDm=J-2;
			    S=H_rada[J];
			    Q=intst((klass-1),S);
			  if (NPQ==0)
			{  for (k=0; k<NV; k++)
			 {  ReS=B_rada[kD][k];
			    temr=ReS*U_rada[kD];
			   for(L=1; L<=kD; L++)
			  { JR=kD-L;
			    ReS=B_rada[JR][k]+S*ReS;
			    temr=B_rada[JR][k]*U_rada[JR]+S*temr;
			  }
			    Y_rada[k]=x[k]+Q*(t*V[k]+t2*S*(F1_rada[k]*W1+S*ReS));
			    Z_rada[k]=V[k]+S*t*(F1_rada[k]+S*temr);
			}
		       }
			  else
		       {  for(k=0; k<NV; k++)
			{  ReS=B_rada[kD][k];
			   for(L=1; L<=kD; L++)
			 { JR=kD-L;
			   ReS=B_rada[JR][k]+S*ReS;
			 }
			   Y_rada[k]=x[k]+Q*(t*V[k]+t2*S*(F1_rada[k]*W1+S*ReS));
			}
		       }
             if((XAXA=FORCE_CHIS(Y_rada,Z_rada,tm+S*t,FJ_rada))) return(XAXA);
			   *NF=*NF+1;
			 if (J2==0)
		      {    for(k=0; k<NV; k++)
		       {   temr=BE_rada[JD][k];
			   ReS=(FJ_rada[k]-F1_rada[k])/S;
			   N=La;
			   for(L=0; L<=JDm; L++)
			{  N++;
			   ReS=(ReS-BE_rada[L][k])*R_rada[N];
			}
			   BE_rada[JD][k]=ReS;
			   temr=ReS-temr;
			   B_rada[JD][k]+=temr*W_rada[JD];
			   N=La;
			   for(L=0; L<=JDm; L++)
			 { N++;
			   B_rada[L][k]+=C_rada[N]*temr;
			 }
		       }
		      }
			 else
		       {   J2=0;
			   for(k=0; k<NV; k++)
			 { temr=BE_rada[0][k];
			   ReS=(FJ_rada[k]-F1_rada[k])/S;
			   BE_rada[0][k]=ReS;
			   B_rada[0][k]+=(ReS-temr)*W_rada[0];
			 }
		       }
		      }
			 if (m>=NI)
			{  hsum=ZeRo;
			   val=intst(-ke,tval);
			  for(k=0; k<=NV-1; k++)
			{ bdubl=B_rada[kD][k];
			  hsum+=bdubl*bdubl;
			}
			  hsum=val*sqrt(hsum);
			if (NSF==0)
		       {if (fabs((hsum-Sm)/hsum)<0.01)
			   goto m3;
			   Sm=hsum;
		       }
		      }
		     }
		   m3:
			if (NSF==0)
		     {  tr=pow(SS_rada/hsum,PW)*DIR;
			if (NeS)  tr=xL;
			if ((NeS==0) && (tr/t<=oNe))
		      {   tr*=0.8;
			  Ncount++;
			if (Ncount>10)  goto konets;
			  goto m1;
		      }
			  NSF=1;
		     }
			  for(k=0; k<NV; k++)
		       {  ReS=B_rada[kD][k];
			  for(L=0; L<kD; L++)
			  ReS+=B_rada[L][k];
			  x[k]+=V[k]*t+t2*(F1_rada[k]*W1+ReS);
			if (NsL==0)
			{ ReS=B_rada[kD][k]*U_rada[kD];
			  for(L=0; L<kD; L++)
			  ReS+=B_rada[L][k]*U_rada[L];
			  V[k]+=t*(F1_rada[k]+ReS);
			}
		       }
			  tm+=t;
			  *NS=*NS+1;
		  	if ((x[2]*DIR>0)  && (Z0*DIR<0)) /* && (pz!=1)  && (pz!=2))*/
			  Nk+=DIR;
			  Z0=x[2];
			if (Nper)  goto konets;
			 goto m4;
	       }
			  Nper=0;
            m4:  if((XAXA=FORCE_CHIS(x,V,tm,F1_rada))) return(XAXA);
			  *NF=*NF+1;
			if (NeS==0)
		      {	  tr=pow(SS_rada/hsum,PW)*DIR;
			if (tr/t>SR) tr=SR*t;
		      }
			if (NeS)  tr=xL;
			if ((DIR*(tm+tr))>=(DIR*tF-1e-10))
		      {	  tr=tF-tm;
			  Nper=1;
		      }
			  Q=tr/t;
			  for(k=0; k<NV; k++)
		      {   ReS=oNe;
			  for(J=0; J<=kD; J++)
		       { if (*NS>1)
			  BT_rada[J][k]=B_rada[J][k]-BT_rada[J][k];
			 if (J!=kD)
			{ m=MCONST[J]-1;
			  JD=J+1;
			  for(L=JD; L<=kD; L++)
			 { B_rada[J][k]+=XI_rada[m]*B_rada[L][k];
			  m++;
			 }
			}
			  ReS*=Q;
			  temr=ReS*B_rada[J][k];
			  B_rada[J][k]=temr+BT_rada [J][k];
			  BT_rada[J][k]=temr;
/* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& */
		       }
		      }
			   NL=NI;
			   goto m2;
	     konets:     /*  pr=5;*/
                           return(0);
		   }





/* *********************************************	 */
/*  ------------------------------------------------------  */
/*   BЫЧИCЛEHИE BЛИЯHИЯ ГPABИTAЦИOHHOГO ПOЛЯ ЗEMЛИ,         */
/*   ПPИTЯЖEHИЯ ЛУHЫ И COЛHЦA,  			    */
/*   A TAKЖE CBETOBOГO ДABЛEHИЯ,ATMOCФEPЫ,ПPИЛИBHЫX ГOPБOB  */
/*  ------------------------------------------------------- */
//	  int   FORCE_CHIS(double *x,double *V,double tm,double *F)
 int    JVC_DLL FORCE_CHIS(double *x,double *V,double tm,double *F)
 {
		  double       RR
			    ,  RQ
			    ,  Rx
			    ,  Ry
			    ,  RZZ
			    ,  a
			    ,  v
			    ,  v5
			    ,  vR
			    ,  v1
			    ,  ax
			    ,  ay
			    ,  aZ
			    ,  R4
			    ,  R4D
			    ,  SFI2
			    ,  Rsol
			    ,  alfa
			    ,  G
			    ,  SG
			    ,  sF
			    ,  sr
			    ,  Rs
			    ,  Dm
			    ,  sG
			    ,  delta
			    ,  V1
			    ,  VR
			    ,  xS
			    ,  yS
			    ,  ZS
			    ,  xm
			    ,  Zm
			    ,  yN
			    ,  Vx
			    ,  Vy
			    ,  VZ
			    ,  xn
			    ,  a1
			    ,  s1
			    ,  xS1
			    ,  xL2
			    ,  SZ
			    ,  v2
			    ,  v3
			    ,  xp
			    ,  yp
			    ,  Zp
			    ,  R1
			    ,  U
			    ,  o
			    ,  QL=0
			    ,  Ro
			    ,  F0
                            ;
		  int          L10
			    ,  I1I2
			    ,  I3
			    ,  L1
			    ,  L6
			    ,  L7
			    ,  L9
			    ,  L4
			    ,  L8
			    ,  L13
			    ,  L5
			    ,  I1
			    ,  I2
			    ,  J1
			    ,  J2
			    ,  J0
			    ,  k1
			    ,  m
			    ,  L
			    ,  k
			    ;

			SZ=S0+tm*CBZ;
			v2=sin(SZ);
			v3=cos(SZ);
/*    BЫЧИCЛEHИE KOOPДИHAT B ГPИHBИЧCKOЙ CИCTEME	 */
			xp=x[0]*v3+x[1]*v2;
			yp=x[1]*v3-x[0]*v2;
			Zp=x[2];
			R1=sqrt(xp*xp+yp*yp+Zp*Zp);
			RR=1/R1;
			RQ=RR*RR;
			Rx=xp*RQ;
			Ry=yp*RQ;
			RZZ=Zp*RQ;
			a=0;
			vR=0;
			v1=0;
/*    BЫЧИCЛEHИE CФEPИЧECKИX ПOЛИHOMOB V,Wm      	 */
			Vm[0]=RR;
			Wm[0]=0;
			for(m=1; m<=MI1; m++)
		      { xL2=2*m-1;
			ax=xL2*Rx;
			ay=xL2*Ry;
			aZ=xL2*RZZ;
			J1=(m-1)*Nm-((m-1)*m)/2+m-1;
			L10=(m-1)*Nm-((m-1)*m)/2+m;
			Vm[L10]=aZ*Vm[J1];
			Wm[L10]=aZ*Wm[J1];
			I1=m*Nm-((m+1)*(m-2))/2-1;
			Vm[I1]=ax*Vm[J1]-ay*Wm[J1];
			Wm[I1]=ay*Vm[J1]+ax*Wm[J1];
		      }
		      if (Nya>=2)
		    {  for(L=2; L<=Nya; L++)
		     { I1=L;
		       J1=I1-1;
		       k1=I1-2;
		       xn=RQ/L;
		       a1=(2*L-1)*Zp*xn;
		       s1=(L-1)*xn;
		       Vm[I1]=a1*Vm[J1]-s1*Vm[k1];
		       Wm[I1]=a1*Wm[J1]-s1*Wm[k1];
		     }
		    }
		     if (Nya>=3)
		    {  for(L=3; L<=Nya; L++)
		     { I1=Nm+L-1;
		       J1=I1-1;
		       k1=I1-2;
		       xn=RQ/(L-1);
		       a1=(2*L-1)*Zp*xn;
		       s1=L*xn;
		       Vm[I1]=a1*Vm[J1]-s1*Vm[k1];
		       Wm[I1]=a1*Wm[J1]-s1*Wm[k1];
		      }
		     }
		     if (MTS1>=3)
		  {   for(m=3; m<=MTS1; m++)
		   {  L1=m+1;
		      L6=m-1;
		      L7=m-2;
		      L9=(m-1)*Nm-((m-1)*m)/2+1;
		      for(L=L1; L<=NTS; L++)
		    { I1=L9+L-1;
		      k1=I1-2;
		      J1=I1-1;
		      xn=RQ/(L-L6);
		      a1=(2*L-1)*Zp*xn;
		      s1=(L7+L)*xn;
		      Vm[I1]=a1*Vm[J1]-s1*Vm[k1];
		      Wm[I1]=a1*Wm[J1]-s1*Wm[k1];
		    }
		   }
		  }
/*    BЫЧИCЛEHИE ПPABЫX ЧACTEЙ OT HEЦEHTPAЛЬHOCTИ ПPИ M:=1 */
		      for(L=0; L<NZ; L++)
		    { I1=Nm+L;
		      a-=C[L]*Vm[I1];
		      vR-=C[L]*Wm[I1];
		      v1-=(L+1)*C[L]*Vm[L+1];
		    }
/*    BЫЧИCЛEHИE ПPABЫX ЧACTEЙ  			 */
		     if (MTS>=2)
		  {   for(m=2; m<=MTS; m++)
		   {  L5=NTS-m+1;
		      L13=(m-1)*Nm-((m-1)*m)/2+m-2;
		      L4=m*Nm-((m+1)*(m-2))/2-2;
		      L8=(m-2)*Nm-((m-2)*(m-1))/2+m-1;
		      for(k=1; k<=L5; k++)
		    { J1=L13+k;
		      I1=L4+k;
		      I2=J1+1;
		      I3=L8+k;
		      xS1=k*(k+1);
		      a+=(C[J1]*(xS1*Vm[I3]-  	        	      /* */
			 Vm[I1])+SS2[J1]*(xS1*Wm[I3]-Wm[I1]))*0.5;
		      vR+=(SS2[J1]*(xS1*Vm[I3]+       	              /* */
			 Vm[I1])-C[J1]*(xS1*Wm[I3]+
			 Wm[I1]))*0.5;
		      v1-=k*(C[J1]*Vm[I2]+   	        	      /* */
			 SS2[J1]*Wm[I2]);
		    }
		   }
		  }
/*    ПPABЫE ЧACTИ OT HEЦEHTPAЛЬHOCTИ ГPAB ПOЛЯ ЗEMЛИ    */
		      F[0]=(v3*a-v2*vR)*ym;
		      F[1]=(v2*a+v3*vR)*ym;
		      F[2]=v1*ym;
#if BOKO
		     if (PATMOSFERA+PSUNDAB+PLUNASUN!=0)
		    {
		      KOOPLS(tm, &xm, &yN, &Zm, &xS, &yS, &ZS);
		      xm/=RZ;
		      yN/=RZ;
		      Zm/=RZ;
		      xS/=RZ;
		      yS/=RZ;
		      ZS/=RZ;
		      o=1/sqrt(xm*xm+yN*yN+Zm*Zm);
		      R4=sqrt(xS*xS+yS*yS+ZS*ZS);
		      R4D=1/R4;
		    }
		     if (PLUNASUN!=0)
		    { v=(x[0]-xm)*(x[0]-xm)+(x[1]-yN)*(x[1]-yN)+(x[2]-Zm)*(x[2]-Zm);
		      v=SWW/(v*sqrt(v));
		      v5=SWW*o*o*o-v;
		      F[0]-=xm*v5+x[0]*v;
		      F[1]-=yN*v5+x[1]*v;
		      F[2]-=Zm*v5+x[2]*v;
		      v=(x[0]-xS)*(x[0]-xS)+(x[1]-yS)*(x[1]-yS)+(x[2]-ZS)*(x[2]-ZS);
		      U=SS/(v*sqrt(v));
		      v3=SS*R4D*R4D*R4D-U;
		      F[0]-=xS*v3+x[0]*U;
		      F[1]-=yS*v3+x[1]*U;
		      F[2]-=ZS*v3+x[2]*U;
		    }
#endif
		    if (PLUNASUN+PSUNDAB!=0)
		   {  delta=asin(ZS*R4D);
		      alfa=atan2(yS/xS,xS);
		   }
#if HOKO
		     if (PATMOSFERA!=0)
	       {      SFI2=RZZ*x[2];
		      HR=(R1-1+ECJ*SFI2)*RZ;
                      QL=KB0;
 /*
	      viichkb(pkb,kb0,W1,W2,W3,W4,DT0+tm-dzv,Dt1,Dt2,Dt3,&QL);
 */
		     if ((HR<=1500)&&(QL>=1.0e-16))
		{    if (HR<0) HR=0.001;
		      Vx=V[0]+CBZ*x[1];
		      Vy=V[1]-CBZ*x[0];
		      VZ=V[2];
		      VR=sqrt(Vx*Vx+Vy*Vy+VZ*VZ);
		      skor=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])*RZ/86400.0;

 		      J0=(HR<120)?-1:((HR<180)?0:((HR<600)?1:2));
		     sr=1.0;
		     if (J0!=-1)
		  {  	        	        	        	       /* PACЧET */
  /*
		      oprFkr(tm);
  */
						 /* PACЧET ПЛOTH. */
                    F107KP(DT0+tm,&FT,&F135,&KP);
		    for(k=0;(fabs(F135-MF0[k])>fabs(F135-MF0[k+1]) && (k<=5));k++);
		    F0=MF0[k];
		    AT1=(AT+k*75+J0*25);
#if BOKO

//		     else
#endif
		     ALFDEL(DT0+tm, &alfa, &delta);
		     G=alfa+AT1[14]-SZ;      	        	      /* COЛHEЧHOГO */
		     SG=sin(G);    	        	        	      /* BЗДУTИЯ */
		     sG=cos(G);    	        	        	      /* [ J0=-1] */
		     alfa=sin(delta);       	        	              /* ПЛOTHOCTИ */
		     delta=cos(delta);       	        	              /* C УЧETOM */
		     sF=(Zp*alfa+delta*(xp*sG+yp*SG))*RR;
		     sr=pow((1+sF)*0.5,(AT1[13]+0.006*HR)*0.5);
		     sr*=AT1[9]+AT1[10]*HR+AT1[11]*HR*HR+AT1[12]*HR*HR*HR;
                     sr+=1.0;
/* BЫЧЛ ПPABЫX ЧACTEЙ ДЛЯ OПP.  			 */
/* BOЗMУЩEHИЙ OT ATMOCФEPЫ      			 */
		  }
		     Ro=ATM(HR,F135,FT,F0,KP,DT0+tm);   	        	/* БEЗ BЗДУTИЯ */
		     Ro*=sr;
		     Ro=Ro*QL*RZ*VR;
		     F[0]-=Ro*Vx;
		     F[1]-=Ro*Vy;
		     F[2]-=Ro*VZ;
              nx=Ro*VR*RZ/(86.400*86400*9.81);
		     if((HR<30 && HR>0.01 && skor<0.98*sqrt(GME/(R1*RZ)))|| skor<0.6 )
		   {
		      DTREENTRY=DT0+tm;
		      ALTITUDE=asin(Zp/R1)*RADIAN;
		      LONGITUDE=atan2(yp,xp)*RADIAN;
            	      return(777);
		   }
		}
               }
#endif
#if BOKO
		    if (PSUNDAB!=0)
		  {  v=(x[0]-xS)*(x[0]-xS)+(x[1]-yS)*(x[1]-yS)+(x[2]-ZS)*(x[2]-ZS);
		     v=1/(v*sqrt(v));
		     Rs=RR*R4D;
		     v1=(x[0]*xS+x[1]*yS+x[2]*ZS)*Rs;
		     v5=sqrt(1-v1*v1);
		    if ((v1>=0)||(v5>=RR))
		   { Dm=ReFL*R4*R4*v;
		     F[0]-=(xS-x[0])*Dm;
		     F[1]-=(yS-x[1])*Dm;
		     F[2]-=(ZS-x[2])*Dm;
		   }
		  }
#endif

	      return(0);
	      }

/* *********************************************	 */
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      */
/* ================================================================ */
//   void km_x(double E0[],int prz,double EK[],int pech)   /* у Васи GNSK */
   void   JVC_DLL km_x(double E0[],int prz,double EK[],int pech)
/* НАЗНАЧЕНИЕ: Пересчет кеплеровых(k) элементов (E0[])
	       в ГНСК(x) (EK[])
   Кеплеровы : a , e , i , w , дву , М или u в зависимости от prz.
	       Если prz=0 то М ,
		    prz=1 то u.
   ГНСК      : X , Y , Z , Vx , Vy , Vz
*/
	    {
	     double
			 E
		      ,  W
		      ,  U
		      ,  K1
		      ,  K2
		      ,  K3
		      ,  K4
		      ,  F1
		      ,  F2
		      ,  A
		      ,  B
		;  if   (prz==0)
		   {    E=KEPLER(E0[5],E0[1]);
			W=2*atan2(sqrt((1+E0[1])/(1-E0[1]))*sin(0.5*E),cos(0.5*E));
			U=E0[3]+W;
		   }
		   else
		   {    U=E0[5];
			W=U-E0[3];
			E=2*atan2(sqrt((1-E0[1])/(1+E0[1]))*sin(0.5*W),cos(0.5*W));
		   }
		;  K3=cos(E0[2])
		;  K4=sin(E0[2])
		;  K1=cos(E0[4])
		;  K2=sin(E0[4])
		;  F1=cos(U)
		;  F2=sin(U)
		;  A=E0[0]*(1-E0[1]*cos(E))
		;  EK[0]=A*(K1*F1-K2*F2*K3)
		;  EK[1]=A*(K2*F1+K1*F2*K3)
		;  EK[2]=A*F2*K4
		;  A=sqrt(GME/(E0[0]*(1-E0[1]*E0[1])))
		;  B=A*(1+E0[1]*cos(W))
		;  A=A*E0[1]*sin(W)
		;  EK[3]=A*(K1*F1-K2*F2*K3)-B*(K1*F2+K2*F1*K3)
		;  EK[4]=A*(K2*F1+K1*F2*K3)-B*(K2*F2-K1*F1*K3)
		;  EK[5]=A*F2*K4+B*F1*K4
		;  if   (pech==1)
		   {
		   /*
		       if   (prz==0)
			     pechvekt("\nKm_X.   Кеплер Э <=> (a,e,i,w,дву,M):\n_rada",E0,0,5);
			else
			     pechvekt("\nKm_X.   Кеплер Э <=> ( a e i w дву u ) :\n_rada",E0,0,5);
			pechvekt("\n_rada               X <=> ( x y z Vx Vy Vz ):\n_rada",EK,0,5);
		   */
		   }
	    }
/* ================================================================ */
/* ================================================================ */
//   void x_k(double h[],double e[], int pech)
 void   JVC_DLL x_k(double h[],double e[], int pech)
/* НАЗНАЧЕНИЕ: Пересчет элементов (h) в ГНСК (h[6])
	       в кеплеровы (k) элементы и u (e[7])
   ГНСК      : ( X , Y , Z , Vx , Vy , Vz )
   Кеплеровы : ( a , e , i , w , дву ,  М ) + u
	       ( u здесь - дополнительная информация )
*/
       {
	     double   c1
	       ,   c2
	       ,   c3
	       ,   c
	       ,   r[3]
	       ,   a[3]
	       ,   rv
	       ,   dr
	       ,   u
	       ,   E
	       ,   k
	       ,   costeta
	       ,   cos2teta;
	     c1=h[1]*h[5]-h[2]*h[4];
	     c2=h[2]*h[3]-h[0]*h[5];
	     c3=h[0]*h[4]-h[1]*h[3];
	     c=sqrt(c1*c1+c2*c2+c3*c3);
	     e[4]=atan2( c1 , -c2 );   /* дву */
	     e[2]=acos( c3/c );         /*  i  */
	     if   (e[2]==0)
		  e[2]=1e-10;
	     r[0]=sqrt(h[0]*h[0]+h[1]*h[1]+h[2]*h[2]); /*  r  */
	     r[1]=sqrt(h[3]*h[3]+h[4]*h[4]+h[5]*h[5]);
	     dr=h[0]*h[3]+h[1]*h[4]+h[2]*h[5];
	     rv=r[0]*r[1];
	     r[2]=atan2( dr , sqrt(rv*rv-dr*dr) );
	     k=r[0]*r[1]*r[1]/GME;
	     costeta=cos(r[2]);
	     cos2teta=costeta*costeta;
	     a[0]=r[0]/(2-k);
	     a[1]=sqrt(1-k*(2-k)*cos2teta);
	     a[2]=atan2(k*sin(r[2])*costeta , k*cos2teta-1);
	     e[0]=a[0];  /*  a  */
	     e[1]=a[1];  /*  e  */
	     e[6]=u=atan2(h[2]/sin(e[2]),h[0]*cos(e[4])+h[1]*sin(e[4]));/*  u  */
	     e[3]=u-a[2];/*  w  */        /* a[2] - ист. аномал */
	     fmod(e[3],PI2);  /*  w  */
	     E=2*atan2(sqrt((1-a[1])/(1+a[1]))*sin(a[2]/2),cos(a[2]/2));
	     e[5]=E-a[1]*sin(E);
	     fmod(e[5],PI2);  /*  M  */
	     if    (pech==1)
	     {
		/*
		   pechvekt("\nXVK.   X<=> (xyzVxVyVz    )    :\n",h ,0,5);
		   pechvekt("\nКеплер Э<=> (a,e,i,w,дву,M) + u:\n",e  ,0,6);
		*/
	     }
       }
/* ================================================================ */
/*  ------------------------------------------------------  */
