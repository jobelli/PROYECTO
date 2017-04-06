#define PECHAT 0
#include "hoko.h"
//#include <graphics.h>                   0
#include <time.h>
#define TRACK 0
/*********************** Только для графики***************/
#if TRACK
extern double HR,skor,xp[3],nx;
extern long Nk,NO;
void track(int NT, double T);
extern bd_prox,i_bd_prox;
#endif
/*********************** Конец только для графики***************/
extern  FILE  far *MAEIKP,far *MULTIINDEX;
	FILE  far *MAEIFORCE, far *MREZ,far *MEMORY;


extern double far *AT1;
#if HOKO
extern double AT[450];
extern int MF0[6];
#endif


extern int LM
		,  LT
		,  LKL
		,  LKS
		,  JREZ[MAXREZ]
		,  L
		,  M
		,  PSUNDAB
		,  PATMOSFERA
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
		,  DT0
		,  RO
		,  EK[6]
		,  KA[6][5]
		,  KC[5][2]
        ,  S0
		,  HP
		,  ALTITUDE
		,  LONGITUDE
		,  DTREENTRY
		;

extern double
		F135
	      ,  FT
	      ,  KP
	      ,  KB0
	      ,  AX
	      ,  EX
	      ,  IX
	      ,  AG
	      ,  EG
	      ,  IG

	      ;


#if BOKO
extern double  	far*HL
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
	      , far* QC
	      , far* QD
	      , far* PC
	      , far* PD
	      , far* RC
	      , far* RD
	      , far* TC
	      , far* TD
	      ;
#endif

extern double
		 QQ[LMAX-1]
	      ,  PP[LMAX-1]
	      ,  RR[LMAX-1]
	      ,  TT[LMAX-1]

	      ;


 int      F0=150
		, PZADACHI
        , NPRABCHACTEI;

double
		   RAD
		,  DTBAKB=0
		,  PADEHIE=0
		,  LIFETIME
		,  DTSTART
		,  CE
		,  SE
		,  EA
		,  IA
		,  CIA
		,  SIA
		,  PER
		,  DTK
		,  AAL
		,  AAC
		,  TK
		,  THAX
		,  TKOH
        ,  K1_SUN,K2_SUN,K3_SUN,K4_SUN,K5_SUN,T_SUN=1000000
        ,  SHAG
#if HOKO
//		,  SHAG=10.0 /*ДЛЯ НИЗКИХ 10*/
#else
//		,  SHAG=1.0 /*ДЛЯ ВЫСОКИХ 0.5*/
#endif
		,  CH=0.5
		,  ER[6]
		,  PQ0[6]
		,  PQ1[6]
		,  PQ2[6]
		,  PQ3[6]
		,  PQ4[6]
		,  Q1[6]
		,  E1[6]
		,  ER0[6]
		,  EKP[6]
		,  CW[MAX(2*LMAX+1,7)]
		,  SW[MAX(2*LMAX+1,7)]
		;

/* ########################## НАЧАЛО PROGNOZ ###################*/
/* ########################## НАЧАЛО PROGNOZ ###################*/
/* ########################## НАЧАЛО PROGNOZ ###################*/
/* ======================================================================
   Destination: Prediction of satellite motion.
   Call: PROGNOZ(NTOCH,BX,BIX,TKNK).
   Input data: NTOCH is number of prediction points,
              BX is initial vector of orbital elements,
              TKNK[0:NTOCH-1] is prediction mooments.
   Output data: BIX[0:NTOCH-1] is array of predicted orbital elements.
  ------------------------------------------------------------------- */

#ifdef HOKODLL
  extern "C" int _export _stdcall PROGNOZ(int NTOCH,struct bxprog *BX,struct bixprog *BIX,double TKNK[])
#else
  int  PROGNOZ(int NTOCH,struct bxprog *BX,struct bixprog *BIX,double TKNK[])
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



/*	LABEL NZIKLA,MREENTRY,KONEZ;  */

int
        XAXA=0
      , POCK=0
      , PRMODEL=0
      , J
      , JD
      , DIR
      , DIRSHAGA=0
      , NT=0              /*номер точки прогноза*/
      , LMOLD=0
      , LTOLD=0
      , OLDMODEL
      , JQC=MAXREZ*(2*LMAX+1)
      , JQC_double=JQC*sizeof(double)
      , JLE=(LLUN-1)*(LLUN+1)*(LLUN+1)
      , JLE_double=JLE*sizeof(double)
	;
double   FOLDMODEL;
double	A0
     ,  BETTA
     ,  P
     ,  I0
     ,  EKC
     ,  L0
     ,  H0
     ,  U0
     ,  W
     ,  MA
     ,  DBY
     ,  U
     ,  H
     ,  DA=0
     ,  DKB=0
     ,  DDT=0
     ,  DT=0
     ,  Dt=0
     ,  T=0
     ,  CH2
     ,  DTFINISH
     ,  DM
     ,  DH
     ,  C1
     ,  C2
     ,  C3
     ,  TKP=10E10
     ,  E0[6]
	;

 time_t time1, time2;

int	PKYTTA=0
       ,POCKYL=0
       ,POCKYLOUT=0
       ,PPC=1;

     /*********BX-СТРУКТУРА bxprog,ГДЕ СФОРМИРОВАНЫ ВХОДНЫЕ ДАННЫЕ**********/
     /**********************************************************/
#if TRACK
NO=BX->NO;
#endif
DT0=BX->DT;
A0=BX->A;
I0=BX->I;
DBY=BX->DBY;
L0=BX->L;
H0=BX->H;
U0=BX->U;
KB0=BX->KB;
F135=BX->F135;
FT=BX->FT;
KP=BX->KP;
PRMODEL=BX->PRMODEL;
POCK=BX->POCK;
PZADACHI=BX->PZADACHI;
NPRABCHACTEI=0;


                 JD=DT0  /*12053-ЧИСЛО СУТОК ОТ 31.12.57 ДО 31.12.90*/
    /*******************КОНЕЦ СЧИТЫВАНИЯ ВХОДНЫХ ДАННЫХ*************/
		      ;  if (POCK>=0)   POCKYL=1
		      ;  if (POCK>-2)   POCKYLOUT=1
		      ;  BETTA=L0*L0+H0*H0
		      ;  EKC=sqrt(BETTA)
		      ;  BETTA=1-BETTA
;if (A0<1880 && A0>85)
{ DTBAKB=1;
  DT=KB0/1440;
  if(DT>-1E-9) {KB0=0;DT=0;} else KB0=0.001;
  PER=A0/1440;
  A0=pow(GME*sqr(PER*86400.0/PI2),1.0/3.0);
  P=A0*BETTA;
  A0*=1-C20*sqr(RZ/P)*((2-2.5*sqr(sin(I0)))*sqrt(pow(BETTA,3))/sqr(1+L0)+
	pow((1+L0),3)/BETTA);
}

 /*****************РАСЧЕТ АРГУМЕНТА ПЕРИГЕЯ********************/
		      ; if(fabs(L0)<1E-8) L0=1E-8
		      ; W=atan2(H0,L0)
  /*************************************************************/
		      ;	 BIX->HP=A0*(1-EKC)-RZ*(1-ECJ*sqr(sin(I0)*sin(W)))
 /*************************ПРОВЕРКА НА ДОСТОВЕРНОСТЬ ВХОДА*****************/
		      ;  if (BIX->HP<90.0)
			 {
			 printf("\nВХОДНЫЕ ДАННЫЕ ВЫХОДЯТ ЗА ПРЕДЕЛЫ ДОПУСТИМОГО\n");
                          close_prognoz();
		      return(1);
			 }
 /**************************************************************************/
#if !defined(GEM_10)
 /************************Ввод коэффициентов гравитационного поля земли*****************************/
	 if (! GEM_input())

	 {perror("\nError of GEM coefficients input\n");close_prognoz();return(2);}

 /***************************************************************************/
#endif 

 /****************ОПРЕДЕЛЕНИЕ ФАЙЛОВ ХРАНЕНИЯ РЕЗУЛЬТАТОВ ПРОГНОЗА***********/
 /************************УСТАНОВКА ПАМЯТИ*****************************/
	 if(MEMORY==NULL)
	 {
     if((MEMORY=fopen("CMEMORY.jvc","r+b"))!=NULL) ;
	 else if((MEMORY=fopen("CMEMORY.jvc","w+b"))!=NULL) ;
	 else {perror("\nОШИБКА ОТКРЫТИЯ ФАЙЛА MEMORY\n");close_prognoz();return(2);}
	 }
	 else  fseek(MEMORY,0L,SEEK_SET);
 /***************************************************************************/
			BETTA=sqrt(BETTA)
		      ; if  (I0<1.0E-8)  I0=1.0E-8
 /********************ОПРЕДЕЛЕНИЕ МОДЕЛИ ДВИЖЕНИЯ**********************************/
		       ; switch(PRMODEL)
			 {
            case 0:   {    PATMOSFERA=1
			      ;  LM=LMAX;LT=LMAX
			      ;  if  (A0<10000)
				   { PSUNDAB=0;PLUNASUN=0  ;}
                 else { PSUNDAB=0;PLUNASUN=1;}
			     break
			     ;};
            case 1:   {    PATMOSFERA=0
			      ;  LM=2;LT=0
			      ;  PSUNDAB=0;PLUNASUN=0;
			     break
			     ;};
			 default:
			     {
                PLUNASUN=0;PATMOSFERA=1;PSUNDAB=0;
                       LM=MIN(PRMODEL,LMAX);LT=MIN(PRMODEL,LMAX)
			     ;}
			 }
 /****************КОРРЕКТИРОВКА МОДЕЛИ ДВИЖЕНИЯ************/
		    ; if  (KB0==0)  PSUNDAB=0
		    ; if  (KB0==0 || BIX->HP>1500  )   PATMOSFERA=0 ;
 /**********************************************************/
 /****************КОНЕЦ ОПРЕДЕЛЕНИЯ МОДЕЛИ ДВИЖЕНИЯ**********************************/

#if BOKO
	       if(
		(QC=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(QD=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(PC=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(PD=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(RC=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(RD=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(TC=(double far*) calloc(JQC,sizeof(double)))==NULL||
		(TD=(double far*) calloc(JQC,sizeof(double)))==NULL)
		{printf("\n НЕТ ПАМЯТИ ДЛЯ ДИНАМИЧЕСКИХ МАССИВОВ QC,QD...В ПРОЦЕДУРЕ PROGNOZ");
		 close_prognoz();return(2);
		 }
				      if  (PLUNASUN)
				    {
	  THAX=1000;
	  TKOH=-1000;
     if  ((HL=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(DL=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(D2L=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(HL0=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(HL1=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(HL2=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(HL3=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(HL4=(double far*) calloc(LKL,sizeof(double)))==NULL
	||(HC=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(DC=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(D2C=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(HC0=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(HC1=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(HC2=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(HC3=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(HC4=(double far*) calloc(LKS,sizeof(double)))==NULL
	||(LE=(double far*) calloc(JLE,sizeof(double)))==NULL
	||(LI=(double far*) calloc(JLE,sizeof(double)))==NULL
	||(LV=(double far*) calloc(JLE,sizeof(double)))==NULL
	||(LPI=(double far*) calloc(JLE,sizeof(double)))==NULL
	||(LLAM=(double far*) calloc(JLE,sizeof(double)))==NULL)
	{printf("\n НЕТ ПАМЯТИ ДЛЯ ДИНАМИЧЕСКИХ МАССИВОВ HL,HC\n");
	 close_prognoz();return(2);
	 }
				    }
#endif



  /*****************РАСЧЕТ СРЕДНЕЙ АНОМАЛИИ********************/
		      ;  IA=U0-W
		      ;  CIA=cos(IA)
		      ;  SE=BETTA*sin(IA)
		      ;  CE=EKC+CIA
		      ;  EA=atan2(SE,CE)
		      ;  MA=EA-EKC*SE/(1+EKC*CIA)
		      ;  MA=fmod(MA,PI2)
 /*************************************************************/
 /***************ФОРМИРОВАНИЕ МАССИВА ЛАГРАНЖЕВЫХ ЭЛЕМЕНТОВ***************/
		      ;  E0[0]=A0
		      ;  E0[1]=L0*cos(DBY)-H0*sin(DBY)
		      ;  E0[2]=H0*cos(DBY)+L0*sin(DBY)
		      ;  E0[3]=sin(0.5*I0)*cos(DBY)
		      ;  E0[4]=sin(0.5*I0)*sin(DBY)
		      ;  E0[5]=MA+W+DBY
		      ;  E0[5]=fmod(E0[5],PI2)
 /************************************************************************/

 /****************ВЫЧИСЛЕНИЕ ЗВЕЗДНОГО ВРЕМЕНИ-S0 НА НАЧАЛЬНЫЙ МОМЕНТ*********/
		      ;  S0=ST(JD,DT0-JD)
 /****************************************************************************/


 /*******************ПРОВЕРКА ВОЗМОЖНОСТИ ИСПОЛЬЗОВАНИЯ ДАННЫХ ПАМЯТИ*************/
		      ; if (fread(PQ0,sizeof(PQ0),1,MEMORY))

			{

			fread(PQ1,sizeof(PQ1),1,MEMORY)
		      ; if (ferror(MEMORY))
			printf("\nОШИБКА ЧТЕНИЯ ФАЙЛА MEMORY\n");
		      ; OLDMODEL=PQ1[0]
	 /********ПРОВЕРКА СОВПАДЕНИЯ НАЧАЛЬНЫХ ДАННЫХ***********/
		      ; for (J=0;J<=5;J++)
			if (fabs(E0[J]-PQ0[J])>1.0E-9)
			 {
			 PPC=0
			 ;break
			 ;}
			if (fabs(DT0-PQ1[2])>1.0E-8)   PPC=0
		      ;	if (fabs(KB0-PQ1[3])>1.0E-8)   PPC=0
                      ;	if (fabs(KB0-PQ0[5])>1.0E-8)   PPC=0
	 /*********ПРОВЕРКА СОВПАДЕНИЯ МОДЕЛЕЙ ДВИЖЕНИЯ***********/
		      ; if (PRMODEL!=OLDMODEL)   PPC=0
	 /*********ПРОВЕРКА СОВПАДЕНИЯ ТИПА ВХОДНЫХ ЭЛЕМЕНТОВ*****/
		      ; if (fabs(POCK-PQ1[1])>1.0E-3)   PPC=0

			    ;} else PPC=0
 /**********************************???????????????**********************************************/
		      ;  if ( PPC && fread(EKP,sizeof(EKP),1,MEMORY))
			 {
			 fread(PQ0,sizeof(PQ0),1,MEMORY)
		       ;  DTSTART=PQ0[0]
		       ;  DTFINISH=PQ0[1]
		       ;  LIFETIME=PQ0[2]
		       ;  PADEHIE=PQ0[3]
		       ;  CH=PQ0[4]           /*ШАГ*/
		       ;  KB0=PQ0[5]           /*kb*/
			 ;}  else
			  {
			   PQ1[0]=PRMODEL
		       ;   PQ1[1]=POCK
		       ;   PQ1[2]=DT0
		       ;   PQ1[3]=KB0
		       ;   PQ1[4]=DT0
		       ;   PQ1[5]=POCK
		       ;  fseek(MEMORY,0L,SEEK_SET)
		       ;  fwrite(E0,sizeof(E0),1,MEMORY)
		       ;  fwrite(PQ1,sizeof(PQ1),1,MEMORY)
		       ;  PKYTTA=1
			   ;}
			 if ((POCKYL || POCKYLOUT) && (LM>2 || LT>0))
			  {
			  if(MAEIKP==NULL)
              if((MAEIKP=fopen("CMAEIKP.jvc","r+b"))!=NULL)
			  ;else if((MAEIKP=fopen("CMAEIKP.jvc","w+b"))!=NULL)
			  ;else {perror("\nОШИБКА ОТКРЫТИЯ ФАЙЛА МАЕIKP\n");close_prognoz();return(2);}
			   else fseek(MAEIKP,0L,SEEK_SET);
			  if(MULTIINDEX==NULL)
              if((MULTIINDEX=fopen("CMULTI.jvc","r+b"))!=NULL)
			  ;else if((MULTIINDEX=fopen("CMULTI.jvc","w+b"))!=NULL)
			  ;else {perror("\nОШИБКА ОТКРЫТИЯ ФАЙЛА MULTIINDEX\n");close_prognoz();return(2);}
			   else fseek(MULTIINDEX,0L,SEEK_SET)
		      ;  if ( fread(&LMOLD,sizeof(LMOLD),1,MULTIINDEX))
			 fread(&LTOLD,sizeof(LTOLD),1,MULTIINDEX)
		      ;  if  ((LM==LMOLD) || (LT==LTOLD))
			     {
			 ;  fread(&AX,sizeof(AX),1,MAEIKP)
			 ;  fread(&EX,sizeof(EX),1,MAEIKP)
			 ;  fread(&IX,sizeof(IX),1,MAEIKP)
			    ; }
			  }

			if (POCKYL)   {
			if (PPC)    TKP=10e10;
			     else {
			 time1=time(NULL)
#if BOKO
		    ;	if  (PLUNASUN && ((T>TKOH) || (T<THAX)))
                APPROK(1.00,T,T+1.0,&THAX,&TKOH);
#endif
		    ;   XAXA=SKP(0,E0,T,EKP)
		    ;   if(XAXA!=0) {close_prognoz();return(XAXA);}
		    ;   time2=time(NULL)
#if PECHAT
		    ;	printf("\nВРЕМЯ СЧЕТАT КР ПРИ Т=0 : %f  sec\n",difftime(time2,time1));
#endif
		    ;   TKP=T
				    ;}
#if PECHAT
		       printf("\nKP ПРИ T=0\n")
		    ; for(J=0;J<=5;J++) printf("%lf",EKP[J])
#endif
		    ; for (J=0;J<=5;J++)  ER0[J]=ER[J]=E0[J]-EKP[J]
			;}
			else for   (J=0;J<=5;J++) {ER0[J]=ER[J]=E0[J];EKP[J]=0.0;}
		       ;   fseek(MEMORY,(long)12*sizeof(double),0)
		       ;   fwrite(EKP,sizeof(EKP),1,MEMORY)
		       ;  DM=0
#if PECHAT
		       ;  printf("\nдо kytta\n")
		       ; for(J=0;J<=5;J++) printf("%  lf",ER[J])
#endif
		       ; for (J=0;J<MAXREZ;J++) JREZ[J]=0

/*% OПPEДEЛEHИE ШAГA H,BPEMEH HAЧAЛA И KOHЦA(T И TK)*/
  ; if(!DTBAKB)
  {
    LBK(ER,Q1)
  ; PER=PI2*ER[0]/86400.0*sqrt(ER[0]/GME)
  ; P=(1-Q1[1]*Q1[1])
  ; PER*=1+1.5*C20*sqr(RZ/(P*ER[0]))*sqrt(P)*((2-2.5*sqr(sin(Q1[2])))*P/sqr(1+Q1[1]*cos(Q1[4]))+
	(1-1.5*sqr(sin(Q1[2]))));

   }
   SHAG=10;
   if (SHAG*PER>1.0)  SHAG=1;
 /************НАЧАЛО ЦИКЛА ПО MОМЕНТАМ ПРОГНОЗА****************************************************************/




    do
   {
   /*начало цикла по моментам прогноза при прогнозе на
		     время и номер витка */
		     if (DTBAKB==0)
		    {
		      if (PZADACHI==0)
			 {

			 DTK=TKNK[NT]
		      ;  TK=DTK-DT0-T
		      ;  if(fabs(TK)>100000)
			  {
			 printf("\nНЕПРАВИЛЬНО ОПРЕДЕЛЕНЫ МОМЕНТЫ ПРОГНОЗИРОВАНИЯ")
		      ;  close_prognoz();return(1)
			 ;}
			  }
			else if (PZADACHI==1)
			 {
			   BIX->bix.NH=(long)TKNK[NT]
			 ; TK=PER*(BIX->bix.NH-BX->NH)-T
			 ; DTK=DT0+TK+T
			 ;  if(fabs(TK)>100000)
			  {
			 printf("\nНЕПРАВИЛЬНО ОПРЕДЕЛЕНЫ МОМЕНТЫ ПРОГНОЗИРОВАНИЯ")
		      ; close_prognoz(); return(1)
			  ;}
			  }
			else if (PZADACHI==2)
			      {
			       CIA=cos(W)
			    ;  SE=-BETTA*sin(W)
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
		      ;	 DTK=TKNK[0]+0.01*PER*NT+0.00000001
		      ;  TK=DTK-DT0-T

			 ;}
		      }
		      else
		      {
			   TK=-PER
			 ; DTK=DT0+TK+T
			 ; NT=-1
			 ;  if(fabs(TK)>100000)
			  {
			 printf("\nНЕПРАВИЛЬНО ОПРЕДЕЛЕНЫ МОМЕНТЫ ПРОГНОЗИРОВАНИЯ")
		      ;  close_prognoz();return(1)
			 ;}
		      }

 /************НАЧАЛО ЦИКЛА ПО ВРЕМЕНИ****************************************************************/
	 ;  do     /*цикл для пересчета*/
	     {
	    do    /*цикл для  выхода на U=0*/
	     {
				time1=time(NULL)
		       ; DIR=(TK<0) ? -1:1
		       ; H=SHAG*DIR*PER
		       ; if (DIRSHAGA==0) DIRSHAGA=DIR;
			 else DIRSHAGA=sign(DTFINISH-DTSTART);
		       ;  if ((PKYTTA) || ((DT0+T+TK-DTSTART)*DIRSHAGA<0) || ((DTFINISH-DT0-T-TK)*DIRSHAGA<0)
			  || (PRMODEL!=OLDMODEL))
			     {

				    if(NPRABCHACTEI!=0 && DTBAKB==0)
				  { if(DIRSHAGA==DIR)
				    {T=DTSTART-DT0+CH;TK=DTK-DT0-T;
				     for (J=0;J<=5;J++)  ER[J]=ER0[J];
				    } else
				    {T=DTSTART-DT0;TK=DTK-DT0-T;}
				  }
				 ; if(MAEIFORCE==NULL )
				 {
                 if((MAEIFORCE=fopen("CMAEIFORCE.jvc","r+b"))!=NULL)
				 ;else if((MAEIFORCE=fopen("CMAEIFORCE.jvc","w+b"))!=NULL)
				 ;else {perror("\nОШИБКА ОТКРЫТИЯ ФАЙЛА MAEIFORCE\n");close_prognoz();return(2);}

                 if((MREZ=fopen("CMREZ.jvc","r+b"))!=NULL)
				 ;else if((MREZ=fopen("CMREZ.jvc","w+b"))!=NULL)
				 ;else {perror("\nОШИБКА ОТКРЫТИЯ ФАЙЛА MREZ\n");close_prognoz();return(2);}
				 }
				  else {fseek(MAEIFORCE,0L,SEEK_SET);
					fseek(MREZ,0L,SEEK_SET);
				       }
				 ; if (fread(&FOLDMODEL,sizeof(FOLDMODEL),1,MAEIFORCE))

				   OLDMODEL=FOLDMODEL;
				     else OLDMODEL=-10

				 ; if  (PRMODEL==OLDMODEL)
				  {
				    fread(&AG,sizeof(AG),1,MAEIFORCE)
				 ;  fread(&EG,sizeof(EG),1,MAEIFORCE)
				 ;  fread(&IG,sizeof(IG),1,MAEIFORCE)
				 ;  fread(KC,sizeof(KC),1,MAEIFORCE)
				 ;  if  (LM>1)
				   {
				    fread(QQ,sizeof(QQ),1,MAEIFORCE)
				 ;  fread(PP,sizeof(PP),1,MAEIFORCE)
				 ;  fread(RR,sizeof(RR),1,MAEIFORCE)
				 ;  fread(TT,sizeof(TT),1,MAEIFORCE)
				 ;
#if BOKO
				    if ((LT>0)&&fread(JREZ,sizeof(JREZ),1,MREZ))
				    {


				     if ((JREZ[0]!=0) && fread(QC,JQC_double,1,MAEIFORCE))
				       {
				    fread(QD,JQC_double,1,MAEIFORCE)
				 ;  fread(PC,JQC_double,1,MAEIFORCE)
				 ;  fread(PD,JQC_double,1,MAEIFORCE)
				 ;  fread(RC,JQC_double,1,MAEIFORCE)
				 ;  fread(RD,JQC_double,1,MAEIFORCE)
				 ;  fread(TC,JQC_double,1,MAEIFORCE)
				 ;  fread(TD,JQC_double,1,MAEIFORCE)
				       ;}
				    }
#endif
				   }
#if BOKO
				    if  (PLUNASUN)
				    {
				   fread(LE,JLE_double,1,MAEIFORCE)
				;  fread(LI,JLE_double,1,MAEIFORCE)
				;  fread(LV,JLE_double,1,MAEIFORCE)
				;  fread(LPI,JLE_double,1,MAEIFORCE)
				;  fread(LLAM,JLE_double,1,MAEIFORCE)

				    ;}
#endif
#if HOKO

				   if  (PATMOSFERA && (BIX->HP<1500) && (KB0>0))

				   fread(KA,sizeof(KA),1,MAEIFORCE);

#endif
				  }
				     do
				  {
				     CH=H
#if BOKO
				  ;   if  (PLUNASUN && ((T>TKOH) || (T<THAX)||(T+H>TKOH) || (T+H<THAX)))
                      APPROK(1.0,T,TK,&THAX,&TKOH);
#endif
#if PECHAT1
		       ;  printf("\nдо force1\n")
		       ; for(J=0;J<=5;J++) printf("%  lf",ER[J])
#endif
				  ;  if(XAXA=FORCE(T,ER,PQ1))  goto MREENTRY
				  ;  C1=LIFETIME-DT0
//********************************************8.05.97********************
              ;  C1=LIFETIME-DT0-T
   			  ; if (fabs(C1)<10 && PATMOSFERA)          /*Р А Б О Т А  П О   С Г О Р А Ю Щ И М*/
				    H=0.5*DIR*fabs(C1)*PER
				  ;
//********************************************8.05.97********************
				  ; if ((T+H)>C1 && (T+H)>TK && PATMOSFERA)          /*Р А Б О Т А  П О   С Г О Р А Ю Щ И М*/
				      CH=0.1*PER
				    ;else CH=H
				  ; CH2=0.5*CH
				  ; for (J=0;J<=5;J++)
				      Q1[J]=ER[J]+CH2*PQ1[J]
#if PECHAT1
		       ;  printf("\nдо force2\n")
		       ; for(J=0;J<=5;J++) printf("%  lf",Q1[J])
#endif
				  ; if(XAXA=FORCE(T+CH2,Q1,PQ2))  goto MREENTRY
				  ; for (J=0;J<=5;J++)
				      Q1[J]=ER[J]+CH2*PQ2[J]
#if PECHAT1
		       ;  printf("\nдо force3\n")
		       ; for(J=0;J<=5;J++) printf("%  lf",Q1[J])
#endif
				  ; if(XAXA=FORCE(T+CH2,Q1,PQ3))  goto MREENTRY
				  ; for (J=0;J<=5;J++)
				      Q1[J]=ER[J]+CH*PQ3[J]
#if PECHAT1
		       ;  printf("\nдо force4\n")
		       ; for(J=0;J<=5;J++) printf("%  lf",Q1[J])
#endif
				  ; if(XAXA=FORCE(T+CH,Q1,PQ4))   goto MREENTRY
				  ; DTSTART=DT0+T
				  ; DTFINISH=DTSTART+CH
				  ; if ((CH-TK)*DIR<-1.0e-9)
				     {
				  ; for (J=0;J<=5;J++)
				       ER[J]+=CH*(PQ1[J]+2.0*PQ2[J]+2.0*PQ3[J]+PQ4[J])/6
				  ; ER[5]=fmod(ER[5],PI2)
				  ; TK-=CH
				  ; T+=CH
				      ;}
				  }
				  while ((DTFINISH-DTK)*DIR<-1.0e-9)
				  ;   PQ0[0]=DTSTART
				  ;   PQ0[1]=DTFINISH
				  ;   PQ0[2]=LIFETIME
				  ;   PQ0[3]=PADEHIE
				  ;   PQ0[4]=CH
				  ;   PQ0[5]=KB0
				  ;   OLDMODEL=PRMODEL
				  ;   fseek(MEMORY,(long)18*sizeof(double),SEEK_SET)
				  ;   fwrite(PQ0,sizeof(PQ0),1,MEMORY)
				  ;   fwrite(ER,sizeof(ER),1,MEMORY)
				  ;   fwrite(PQ1,sizeof(PQ1),1,MEMORY)
				  ;   fwrite(PQ2,sizeof(PQ2),1,MEMORY)
				  ;   fwrite(PQ3,sizeof(PQ3),1,MEMORY)
				  ;   fwrite(PQ4,sizeof(PQ4),1,MEMORY)
				  ;  rewind(MAEIFORCE)
				  ;  FOLDMODEL=PRMODEL
				  ;  fwrite(&FOLDMODEL,sizeof(FOLDMODEL),1,MAEIFORCE)
				  ;  fwrite(&AG,sizeof(AG),1,MAEIFORCE)
				  ;  fwrite(&EG,sizeof(AG),1,MAEIFORCE)
				  ;  fwrite(&IG,sizeof(AG),1,MAEIFORCE)
				  ;  if  (LM>1)    fwrite(KC,sizeof(KC),1,MAEIFORCE);
				  ;  if  (LM>1)
				   {

				   fwrite(QQ,sizeof(QQ),1,MAEIFORCE)
				  ;  fwrite(PP,sizeof(PP),1,MAEIFORCE)
				  ;  fwrite(RR,sizeof(RR),1,MAEIFORCE)
				  ;  fwrite(TT,sizeof(TT),1,MAEIFORCE)
				   ;}
#if BOKO
				   ;  if  (LT>0)
				   {
				     rewind(MREZ)
				  ;  fwrite(JREZ,sizeof(JREZ),1,MREZ)
				  ;  if   (JREZ[0]!=0)
				   {
				     fwrite(QC,JQC_double,1,MAEIFORCE)
				  ;  fwrite(QD,JQC_double,1,MAEIFORCE)
				  ;  fwrite(PC,JQC_double,1,MAEIFORCE)
				  ;  fwrite(PD,JQC_double,1,MAEIFORCE)
				  ;  fwrite(RC,JQC_double,1,MAEIFORCE)
				  ;  fwrite(RD,JQC_double,1,MAEIFORCE)
				  ;  fwrite(TC,JQC_double,1,MAEIFORCE)
				  ;  fwrite(TD,JQC_double,1,MAEIFORCE)
				  ;}
				   }
				     if  (PLUNASUN)
				    {
				     fwrite(LE,JLE_double,1,MAEIFORCE)
				  ;  fwrite(LI,JLE_double,1,MAEIFORCE)
				  ;  fwrite(LV,JLE_double,1,MAEIFORCE)
				  ;  fwrite(LPI,JLE_double,1,MAEIFORCE)
				  ;  fwrite(LLAM,JLE_double,1,MAEIFORCE)

				    ;}
#endif
#if HOKO
				    if  (PATMOSFERA && (BIX->HP<1500) && (KB0>0))

				   fwrite(KA,sizeof(KA),1,MAEIFORCE);

#endif
				 } else {
				      if(NPRABCHACTEI==0)
				      {
				      fseek(MEMORY,(long)24*sizeof(double),SEEK_SET)
				  ;   fread(ER,sizeof(ER),1,MEMORY)
				  ;   fread(PQ1,sizeof(PQ1),1,MEMORY)
				  ;   fread(PQ2,sizeof(PQ2),1,MEMORY)
				  ;   fread(PQ3,sizeof(PQ3),1,MEMORY)
				  ;   fread(PQ4,sizeof(PQ4),1,MEMORY)
				  ;   NPRABCHACTEI++
				       ;}
				       }
				  ;   PKYTTA=0
				  ;   T+=TK
				  ;   DH=(DT0+T-DTSTART)/CH
				  ;   C1=DH*(2.0*DH*DH/3.0-1.5*DH+1.0)
				  ;   C2=DH*(DH-2.0*DH*DH/3.0)
				  ;   C3=DH*(2.0*DH*DH/3.0-0.5*DH);
				  ; for (J=0;J<=5;J++)
				    Q1[J]=ER[J]+CH*(C1*PQ1[J]+C2*(PQ2[J]+PQ3[J])+C3*PQ4[J])
				  ; Q1[5]=fmod(Q1[5],PI2)
		       ;        time2=time(NULL)
#if PECHAT
		       ;	printf("\nВРЕМЯ СЧЕТА KYTTA : %f  sec\n",difftime(time2,time1));
		       ;  printf("\nPOSLE kytta T=%lf \n",T)
		       ;  for(J=0;J<=5;J++) printf("%  lf",Q1[J])
#endif
			    ; if (POCKYLOUT)   {
				if (fabs(TKP-T)>0.000001*PER)
				     {
					time1=time(NULL)
#if BOKO
				     ;  if  (PLUNASUN && ((T>TKOH) || (T<THAX)))
					APPROK(1.00,T,T+1.0,&THAX,&TKOH);
#endif
				     ;  XAXA=SKP(1,Q1,T,EKP)
				     ;  if(XAXA!=0) {close_prognoz();return(XAXA);}
				     ;  time2=time(NULL)
#if PECHAT
				     ;  printf("\nВРЕМЯ СЧЕТА KP при Т=ТК : %f  sec\n",difftime(time2,time1));
				     ;  printf("\nKP при Т=ТК=%lf\n",T)
				     ;  for(J=0;J<=5;J++) printf("%lf",EKP[J])
#endif
				     ;  TKP=T
				     ;  }
				       for (J=0;J<=5;J++) E0[J]=Q1[J]+EKP[J]
				   ;}
				    else
				      for (J=0;J<=5;J++) E0[J]=Q1[J]
				   ;  LBK(E0,E1)
			    ; if  (PZADACHI!=0&&PZADACHI!=3 || DTBAKB)
				{
			       EA=KEPLER(E1[5],E1[1])
			    ;  CE=cos(EA)
			    ;  SE=sin(EA)
			    ;  RAD=E1[0]*(1-E1[1]*CE)
			    ;  BETTA=sqrt(1-E1[1]*E1[1])
			    ;  SIA=BETTA*SE
			    ;  CIA=CE-E1[1]
			    ;  IA=(SIA==0&&CIA==0) ? 0 : atan2(SIA,CIA)
			    ;  U=E1[3]+IA
			    ;  U=fmod(U,PI2)
			    ;  if(fabs(U)>3.0) U-=sign(U)*PI2
			    ;  TK=-U*(RAD/86400.0)*(RAD/sqrt(GME*E1[0]*BETTA))
			    ;
/*
			       }
*/
			      DTK=DT0+T+TK
				    ;} else
				    TK=0;
	   }
	 while( fabs(TK)>1E-7);
	 if(DTBAKB==1)
	 {
	   Dt=PER+T;
	   TK=-PER ;
	   T=0;
	 if(PATMOSFERA && DT<0)
	 {DDT=DT-PADEHIE;DKB=DDT*KB0/PADEHIE/*DT*/;}else{DDT=DKB=DT=0;};
	  KB0+=DKB;
	  if(KB0<0) KB0=0.0001;
      /*  if(KB0>200) KB0=200;*/
	  DA=2*ER[0]*Dt/(3*PER);
	  ER0[0]+=DA;
	  for (J=0;J<=5;J++)  ER[J]=ER0[J];
	  PKYTTA=1;
	  }
	  }
	 while( fabs(Dt)>1E-7 || (fabs(DDT)>0.01*fabs(DT) && fabs(DDT)>1E-8));
	 ; for (J=0;J<=5;J++) ER0[J]=ER[J]+CH*(PQ1[J]+2.0*PQ2[J]+2.0*PQ3[J]+PQ4[J])/6
	 ; ER0[5]=fmod(ER0[5],PI2);
	 if(DTBAKB==1 && PZADACHI<4) {DTBAKB=0;PKYTTA=0;NT++;continue;}
	 if (PZADACHI==4)
	      {
		 for (J=0;J<=5;J++) Q1[J]=ER[J]
	      ;  XAXA=SKP(1,Q1,0,EKP)
	      ;  if(XAXA!=0) {close_prognoz();return(XAXA);}
	      ;  TKP=0
	      ;  A0=Q1[0]+EKP[0];
	       }
	  DTK=DT0+T
  ; if ((((PZADACHI==0) && (fabs(DTK-DT0)<0.5E-9)) ||((PZADACHI==1) && (BIX->bix.NH==BX->NH) && U0==0) || ((PZADACHI==2) && (U0==0))||(PZADACHI==4))&&POCKYL==POCKYLOUT)
	       {
	      E1[0]=A0
	    ; E1[1]=EKC
	    ; E1[2]=I0
	    ; E1[3]=W
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
		      ;  BIX->bix.U=GHCK(E1,BIX->X)
		      ;  BIX->bix.U=fmod(BIX->bix.U,PI2)
		      ;  if (BIX->bix.U<0) BIX->bix.U+=PI2
		      ; if ((PZADACHI==1) || (PZADACHI==2))  BIX->bix.U=0
		      ; for (J=0;J<=5;J++) BIX->Z[J]=Q1[J]
			    ;  LBK(Q1,EK)
			    ;  BIX->Y[0]=EK[0]
			    ;  BIX->Y[1]=EK[2]
			    ;  BIX->Y[2]=EK[4]
			    ;  BIX->Y[3]=EK[1]*cos(EK[3])
			    ;  BIX->Y[4]=EK[1]*sin(EK[3])
			    ;  EA=KEPLER(EK[5],EK[1])
			    ;  BIX->Y[5]=EK[3]+2*atan2(sqrt((1+EK[1])/(1-EK[1]))*sin(EA/2),cos(EA/2))
			    ;  BIX->Y[5]=fmod(BIX->Y[5],PI2)
;
if (!DTBAKB)
/*PER=PI2/86400*E1[0]*sqrt(E1[0]/GME)*(1+1.5*C20*sqr(RZ/E1[0]/(1-EKC*EKC))*((2-2.5*sqr(sin(I0)))*
sqrt(pow((1-EKC*EKC),3))/sqr(1+L0)+pow((1+L0),3)/(1-EKC*EKC)));
*/
  {  PER=PI2*EK[0]/86400.0*sqrt(EK[0]/GME)
  ; P=(1-EK[1]*EK[1])
  ; PER*=1+1.5*C20*sqr(RZ/(P*EK[0]))*sqrt(P)*((2-2.5*sqr(sin(EK[2])))*P/sqr(1+BIX->Y[3])+
	(1-1.5*sqr(sin(EK[2]))));
   }
			if  (PZADACHI==0||PZADACHI==3)
			  { BIX->bix.NH=BX->NH+(T/PER)
			    ;  BIX->bix.NH=BX->NH+(T/(PER+0.5*PADEHIE*(BIX->bix.NH-BX->NH)))
			      ;}
			/*    ; if  ((TK>0) && (U0<BX[8]))
				BIX->bix.NH=BIX->bix.NH+1
			    ; if  (TK<0) && (U0>BX[8])
				BIX->bix.NH=BIX->bix.NH-1
			      else DTK=DT0+TK   */
			    ; if  (PZADACHI==2||PZADACHI>3)
				BIX->bix.NH=BX->NH


/***************ЗАПИСЬ РЕЗУЛЬТАТА***************************************/
  ;BIX->HA=BIX->bix.A*(1+BIX->E)-RZ*(1-ECJ*sqr(sin(I0)*sin(BIX->W)));HP=BIX->bix.A*(1-BIX->E)-RZ*(1-ECJ*sqr(sin(I0)*sin(BIX->W)))
  ;RAD=sqrt(BIX->X[0]*BIX->X[0]+BIX->X[1]*BIX->X[1]+BIX->X[2]*BIX->X[2])
  ;ALTITUDE=asin(BIX->X[2]/RAD)*RADIAN
/*
  ;LONGITUDE=(atan2(BIX->X[1],BIX->X[0])-S0-CBZ*T)*RADIAN;
  ;LONGITUDE=fmod(LONGITUDE,360)
  ;if (LONGITUDE<0) LONGITUDE+=360
  ;
*/
/*  ;LONGITUDE=(BIX->bix.DBY-S0-CBZ*T+atan2(sin(BIX->bix.U)*cos(BIX->bix.I),cos(BIX->bix.U)))*RADIAN*/
  ;JD=DTK
/*  ;LONGITUDE=(BIX->bix.DBY-ST(JD,DTK-JD)+atan2(sin(BIX->bix.U)*cos(BIX->bix.I),cos(BIX->bix.U)))*RADIAN*/
  ; LONGITUDE=(atan2(BIX->X[1],BIX->X[0])-ST(JD,DTK-JD))*RADIAN;
  ; LONGITUDE=fmod(LONGITUDE,360)
  ; if (LONGITUDE>180) LONGITUDE-=360
  ; if (LONGITUDE<-180) LONGITUDE+=360
  ;

 ;MREENTRY:  ;
BIX->bix.NO=(XAXA) ? 777 : BX->NO;
BIX->bix.DT=DTK;
BIX->bix.KB=KB0;
BIX->bix.F135=F135;
BIX->bix.FT=FT;
BIX->bix.KP=KP;
BIX->bix.PRMODEL=PRMODEL;
BIX->bix.POCK=POCK;
BIX->bix.PZADACHI=PZADACHI;
BIX->PER=PER*1440;
BIX->DT=PADEHIE*1440;
BIX->LIFET=LIFETIME;
BIX->HP= HP;
BIX->RAD= RAD;
BIX->ALTITUDE= ALTITUDE;
BIX->LONGITUDE=LONGITUDE;
#if TRACK
/************************ТРАССА************************/
if(PZADACHI==3)
{   xp[0]=RAD*cos(LONGITUDE/RADIAN)*cos(ALTITUDE/RADIAN);
    xp[1]=RAD*sin(LONGITUDE/RADIAN)*cos(ALTITUDE/RADIAN);
    xp[2]=RAD*sin(ALTITUDE/RADIAN);
 HR=RAD-RZ*(1-ECJ*sqr(sin(I0)*sin(BIX->bix.U)));
 skor=sqrt(sqr(BIX->X[3])+sqr(BIX->X[4])+sqr(BIX->X[5]));
 nx=1000.0*KB0*RO*skor*skor/9.81;
 if(NT==0)
 Nk=BX->NH+((DTK-DT0)/(PER+0.5*PADEHIE*(DTK-DT0)/PER));
else
 if ((xp[2]*sign(DTK-DT0)>0)  && (Z0*sign(DTK-DT0)<0)) /* && (pz!=1)  && (pz!=2))*/
   Nk+=DIR;
   Z0=xp[2];
 if (LONGITUDE>180) LONGITUDE-=360;
 if (LONGITUDE<-180) LONGITUDE+=360;
 track(NT,DT0+T);
}
/************************ТРАССА************************/
if(PZADACHI==3&&XAXA!=0)  {d4close_all();fcloseall();
while(!kbhit());closegraph();}
/************************КОНЕЦ ТРАССЫ************************/

/************************КОНЕЦ ТРАССЫ************************/
#endif
if (XAXA)
{
  fseek(MEMORY,0L,SEEK_SET);
  for (J=0;J<6;J++) E0[0]=0;
  fwrite(E0,sizeof(E0),1,MEMORY);
  close_prognoz();return(XAXA);
}
NT++;

if (PZADACHI<2&& NT<NTOCH) BIX++;
  }
    /*конец цикла по моментам прогноза при прогнозе на
		     время и номер витка */
  while  (PZADACHI<4 && NT<NTOCH);
#if BOKO
	      ;	free(QC)
	      ; free(QD)
	      ; free(PC)
	      ; free(PD)
	      ; free(RC)
	      ; free(RD)
	      ; free(TC)
	      ; free(TD)
	      ; if  (PLUNASUN)
		  {
		free(HL)
	      ; free(DL)
	      ; free(D2L)
	      ; free(HL0)
	      ; free(HL1)
	      ; free(HL2)
	      ; free(HL3)
	      ; free(HL4)
	      ; free(DC)
	      ; free(D2C)
	      ; free(HC)
	      ; free(HC0)
	      ; free(HC1)
	      ; free(HC2)
	      ; free(HC3)
	      ; free(HC4)
	      ; free(LE)
	      ; free(LI)
	      ; free(LV)
	      ; free(LPI)
	      ; free(LLAM)
	      ;}
#endif

#if TRACK
/************************ТРАССА************************/
if(PZADACHI==3)  {d4close();while(!kbhit());closegraph();}
/************************КОНЕЦ ТРАССЫ************************/
#endif

   close_prognoz();return(XAXA)
  ;}

/*******************KOHEЦ ЧИCЛEHHO-AHAЛИT. MODYLЯ ********************/
/* ======================================================================
   Destination: Calculation of right-hand side of averaged equations.
   Call: FORCE(T1,E1,E2).
   Input data: T1 is the time count off from initial,
              E1[0:5] is vector of nonsingular elements.
   Output data: E2[0:5] is vector of right-hand side of averaged equations value.
  ------------------------------------------------------------------- */



		  int JVC_DLL FORCE(double T,double E1[],double  E2[])
		  {
		   double      Z
			    ,  X
			    ,  H1
			    ,  H2
			    ,  A1
			    ,  A2
			    ,  GR
			    ,  L1
			    ,  L2
			    ,  L3
			    ,  L4
			    ,  L5
			    ,  L6
			    ,  N
			    ,  MK
			    ,  AK
			    ,  K1
			    ,  SK
			    ,  CK
			    ,  FR
			    ,  SN
			    ,  K2
			    ,  K3
			    ,  K4
			    ,  K5
			    ,  AP
			    ,  DBY
			    ,  CA
			    ,  DOM
			    ,  DOW
			    ,  DM
			    ,  C1
			    ,  S1
			    ,  JK
			    ,  KJ
			    ,  UG
			    ,  EATM
			    ,  HODH
			    ,  B
			    ,  C2
			    ,  C3
			    ,  B1
			    ,  B2
			    ,  FF
			    ,  N2
			    ;
	      int	       K
			    ,  J0
			    ,  J1
			    ,  J
			    ,  L
			    ,  M
			    ,  P
			    ,  WQ
			    ,  J10
			    ,  J11
			    ,  Q
			    ;

			 A=E1[0]
		      ;  H1=sqr(E1[1])+sqr(E1[2])
		      ;  E=sqrt(H1)
		      ;  H1=1-H1
		      ;  H2=sqrt(H1)
		      ;  SI2=sqr(E1[3])+sqr(E1[4])
		      ;  CI2=1-SI2
		      ;  CI=CI2-SI2
		      ;  SI2=sqrt(SI2)
		      ;  CI2=sqrt(CI2)
		      ;  SI=2*SI2*CI2
		      ;  I=atan(SI/CI)
		      ;  if (I<0)   I+=PI
		      ;  Z=SI*SI
		      ;  AP=atan(E1[2]/E1[1])
		      ;  if (E1[1]<0)   AP+=PI
		      ;  DBY=atan(E1[4]/E1[3])
		      ;  if (E1[3]<0)   DBY+=PI
		      ;  CA=E1[5]-AP
		      ;  AP-=DBY
		      ;  N=sqrt(GME/A)/A*86400.0
		      ;  L1=0.0
		      ;  L2=0.0
		      ;  L3=0.0
		      ;  L4=0.0
		      ;  L5=0.0
		      ;  L6=1.0
		      ;  A1=RZ/A
		      ;  GR=C20*A1*A1/(H1*H1)
		      ;  CW[0]=1.0
		      ;  SW[0]=0.0
		      ;  C1=cos(AP)
		      ;  S1=sin(AP)
		      ; for  (J= 1;J<MAX((2*LMAX),6);J++)
			  {  CW[J]=CW[J-1]*C1-SW[J-1]*S1
			    ;  SW[J]=SW[J-1]*C1+CW[J-1]*S1
			;}
#if HOKO
			if  (PATMOSFERA)
			   {
			       HP=A*(1-E)-RZ*(1-ECJ*Z*S1*S1)
			    ; if  (LM>1)
				HP=HP-0.25*GR*A*H1*((3*Z-2)*(1+2*(1-E)/H2+E/(1+H2))+Z*CW[2])
			    ; if (HP <120)
				J0=-1;
			      else if  (HP<180)
				J0=0;
			      else if  (HP<600)
				J0=1;
			      else (J0=2)
			    ; if  (J0==-1)
				{    B1=0
				  ; if  (HP<80)
				      {
					      B=KEPLER(CA,E)
					      ;  RAD=A*(1-E*cos(B))
					      ;  UG=AP+2*atan2(sqrt((1+E)/(1-E))*sin(B/2),cos(B/2))
					      ;  DTREENTRY=DT0+T
					      ;  ALTITUDE=asin(SI*sin(UG))
					      ;  LONGITUDE=DBY-S0-CBZ*T+atan2(sin(UG)*CI,cos(UG)) /* ЛGROTA*/
                     ;  if (LONGITUDE<-180) LONGITUDE+=360
                     ;  if (LONGITUDE>180) LONGITUDE-=360
					      ;  return(777)
					 ;}
				   if  (HP<100)
				      {    HODH=1/(0.12378+0.00087*(HP-60))
					;  B2=0.0 ;}
				    else { HODH=1/(0.17527-0.001287*(HP-100))
					;  B2=0.001287
					 ;}
			   }    else {
                             F107KP(DT0+T,&FT,&F135,&KP);
				    for(K=0;(fabs(F135-MF0[K])>fabs(F135-MF0[K+1]) && (K<=6));K++);
					  F0=MF0[K]
				  ;  AT1=(AT+K*75+J0*25)
				  ;  B2=AT1[1]
				  ;  HODH=2/B2*sqrt(HP-AT1[2])
				  ;  B2=1/(HODH*HODH*HODH*B2*B2)
				  ;  B1=AT1[10]+2*AT1[11]*HP+3*AT1[12]*HP*HP
				  ;  FF=AT1[9]+AT1[10]*HP+AT1[11]*HP*HP+AT1[12]*HP*HP*HP
//				  ;  N2=0.5*AT1[13]     //было
				  ;  N2=0.5*(AT1[13]+0.006*(HP+0.5*A*E))     //здесь коэффициент зависимости n от высоты n1 задан константой 0.006
				  ;  SK=pow(0.5,N2)
				  ;  FF*=SK
				  ;  B1*=SK/(1+FF)
			      ;}
			 }
#endif
			if ((fabs(1-A/AG)>CTPOBFORCE) || (fabs(1-E/(EG+10e-10))>CTPOBFORCE) || (fabs(1-SI2/(IG+10e-10))>CTPOBFORCE))
			  {
			       AG=A
			    ;  EG=E
			    ;  IG=SI2
			    ;  JREZ[0]=0;
			       if  (LM>1)    C202();
			       if  (LM>2)    QPRT();
#if BOKO
			       if ((LT>0) && (N<=(LT+0.2)*CBZ))   REZ() ;
                   if  (PLUNASUN) LS();
#endif
#if HOKO
			       if  (PATMOSFERA && (HP<1500) && (KB0>0))  KATM(HODH,B1,B2)
#endif
			  ;}
			if  (LM>=2)
			   {
			       DOM=1.5*GR*CI+KC[3][0]+KC[3][1]*CW[2]
			    ;  DOW=0.75*GR*(5*Z-4)+KC[2][0]+KC[2][1]*CW[2]
			    ;  DM=0.75*GR*H2*(3*Z-2)+KC[4][0]+KC[4][1]*CW[2]
			    ;  L4+=E*(DOM+DOW)
			    ;  L5+=SI2*DOM
			    ;  L6+=DOM+DOW+DM
			    ;  L2+=KC[0][1]*SW[2]
			    ;  L3+=KC[1][1]*SW[2]

			    ; if  (LM>=3)
				{
					  for   (K=0;K<=(LM-2);K++)
					     { if (K%2)
					    {    JK=CW[K]
					      ;  KJ=SW[K]
					       ;}
					  else { JK=SW[K]
					      ;  KJ=CW[K]
					       ;}

					   L2+=H2*QQ[K]*JK
					;  L3-=E*CI/(4*SI2*H2)*JK*QQ[K]
					;  DOM=RR[K]*KJ
					;  DOW=PP[K]*KJ
					;  L5+=DOM
					;  L4+=H2*DOW+2*E*SI2*DOM
					;  L6+=(H2-H1)/E*DOW+TT[K]*KJ+2*SI2*DOM
				       ;}
				   }
#if BOKO
			      if ((LT>0) && JREZ[0]!=0 )
				      {
					   SN=S0+CBZ*T
					;  J=1
					;  J1=0
					;  while ((J<=MAXREZ) && (JREZ[J-1]>0))
					    {
						 M=JREZ[J-1]
					      ;  FR=J*(AP+CA)+M*(DBY-SN)
					      ;  C1=cos(FR)
					      ;  S1=sin(FR)
					      ; for  (Q=J-LT;Q<=(J+LT);Q++)
						  {
						       WQ=1
						    ;  if (Q<0)   WQ=-1; else if (Q==0)   WQ=0
						    ;  K=abs(Q)
						    ;  SK=S1*CW[K]-C1*WQ*SW[K]
						    ;  CK=C1*CW[K]+S1*WQ*SW[K]
						    ;  K1=QC[J1]*CK+QD[J1]*SK
						    ;  L1+=2*A*J*K1
						    ;  L2+=H2*(J*(H2-1)+Q)*K1/E
						    ;  L3+=((J-Q)*CI-M)/(4*SI2*H2)*K1
						    ;  DOM=RC[J1]*CK+RD[J1]*SK
						    ;  DOW=PC[J1]*CK+PD[J1]*SK
						    ;  DM=TC[J1]*CK+TD[J1]*SK
						    ;  L4+=H2*DOW+E*(CK=2*SI2*DOM)
						    ;  L5+=DOM
						    ;  L6+=(H2-H1)/E*DOW+CK+DM
						    ;  J1++
						   ;}        /* ЦИКЛ ПО Q*/

						      J++

					     ;}              /* FR=0*/
					  }                  /* LT>0 && JREZ=1*/
#endif
			       }                             /*LM>2*/
#if BOKO
			if  (PLUNASUN)
			  {
 /*
			       if ((T>TKOH) || (T<THAX))
			       APPROK(0.25,T,TK,&THAX,&TKOH);
			       else
 */
				    {
			       X=T-THAX
			    ; for   (K=0;K<LKL;K++)
				{  HL[K]=HL0[K]+(((HL4[K]*X+HL3[K])*X+HL2[K])*X+HL1[K])*X
			    ; if  (K<LKS)
				    HC[K]=HC0[K]+(((HC4[K]*X+HC3[K])*X+HC2[K])*X+HC1[K])*X
				 ;}
				     }
			      AAL=sqr(A/ALUN)
			    ; AAC=sqr(A/AC)
			    ; Q=0
			    ; for   (L=2;L<=LLUN;L++)
				{
				     AAL*=A/ALUN
				  ;  AAC*=A/AC
				  ;  A1=ML*AAL/(2*L+1)
				  ;  A2=MC*AAC/(2*L+1)
				  ; for (P=0;P<=L;P++)
				      {
						 FR=(L-2*P)*AP-DBY
					      ;   for (M=0;M<=L;M++)
						{
						 J10=L*(L+1)-6+2*M
					      ;  J11=J10+1
					      ;  J=(L-M)%2 ? 0 : 1
					      ;  FR+=DBY
					      ;  SK=sin(FR)
					      ;  CK=cos(FR)
					      ;  S1=J*CK+(1-J)*SK
					      ;  C1=J*SK+(J-1)*CK
					      ;  K4=A1*(HL[J10]*S1+HL[J11]*C1)
					      ;  K5=A1*(-HL[J10]*C1+HL[J11]*S1)
					      ; if  (L<=LSUN)
						  {    K4+=A2*(HC[J10]*S1+HC[J11]*C1)
						    ;  K5+=A2*(-HC[J10]*C1+HC[J11]*S1)
						  ;}
 						     L2+=  LE[Q]*K5
					      ;  L3+=  LI[Q]*K5
					      ;  L4+= LPI[Q]*K4
					      ;  L5+=  LV[Q]*K4
					      ;  L6+=LLAM[Q++]*K4
					  ;}
				     }
			       }
			 }
#endif

			if ( (PATMOSFERA && HP<1500) || PSUNDAB )
			  {
				 ALFDEL(DT0+T,&AK,&MK)
#if BOKO
				; if  (PSUNDAB>0)
				{
                    if (fabs(T-T_SUN)>10.0)
                    {
				     SHADOW(AK,MK,AP,DBY,&K1_SUN,&K2_SUN,&K3_SUN,&K4_SUN,&K5_SUN);
                     T_SUN=T;
                    }
				  ;SK=sin(K1_SUN)-sin(K2_SUN)
				  ;CK=cos(K1_SUN)-cos(K2_SUN)
				  ;S1=sin(2*K1_SUN)-sin(2*K2_SUN)
				  ;C1=cos(2*K1_SUN)-cos(2*K2_SUN)
				  ;K1_SUN-=K2_SUN
				  ;GR=1.44*KB0*4.65E-9*(A/GME)*A/PI2
				  ;L1+=2*GR*A*(K3_SUN*CK+K4_SUN*H2*SK)
				  ;L3+=GR*K5_SUN*0.5*CI2/H2*(((1+E*E)*SK-0.25*E*S1)*CW[1]+H2*(CK-0.25*E*C1)*SW[1]-1.5*E*CW[1]*K1_SUN)
				  ;DOM=GR*K5_SUN/(H2*2*CI2)*(((1+E*E)*SK-0.25*E*S1)*SW[1]-H2*(CK-0.25*E*C1)*CW[1]-1.5*E*K1_SUN*SW[1])
				  ;L5+=DOM
				  ;L2+=GR*H2*(0.25*K3_SUN*H2*C1+K4_SUN*(0.25*S1-E*SK)+1.5*K4_SUN*K1_SUN)
				  ;JK=2*SI2*DOM
				  ;K2=GR*H2*(K3_SUN*(E*SK+0.25*S1)+K4_SUN/H2*(E*CK-0.25*C1)-1.5*K3_SUN*K1_SUN)
				  ;KJ=-2*GR*(K3_SUN*((1+E*E)*SK-0.25*E*S1)-K4_SUN*H2*(CK-0.25*E*C1)-1.5*E*K3_SUN*K1_SUN)
				  ;L6+=((1-H2)/E*K2+JK+KJ)
				  ;L4+=(E*JK+K2)
			      ;}
#endif
#if HOKO
			     ; if  (PATMOSFERA && (HP<=1500))
				      {
			  K4=1;K1=K2=K3_SUN=K5=C1=S1=CK=SK=0
					; if  (J0!=-1)
					    {    AK+=AT1[14]
					      ;  K1=cos(DBY-AK)*cos(MK)
					      ;  K2=SI*sin(MK)-sin(DBY-AK)*CI*cos(MK)
					      ;  K3_SUN=sqrt(K1*K1+K2*K2)
					      ;  C1=(CW[1]*K1+SW[1]*K2)/K3_SUN
					      ;  S1=(SW[1]*K1-CW[1]*K2)/K3_SUN
					      ;  C2=0.25*N2*(N2-1)*K3_SUN*K3_SUN
                          ;  C3=N2*(1+0.125*(N2-1)*(N2-2)*K3_SUN*K3_SUN)*K3_SUN
                          ;  K4=1+(1+C2)*FF
					      ;  K1=FF*C3/K4
					      ;  K2=FF*C2/K4
					      ;  SK=2*S1*C1
					      ;  CK=C1*C1-S1*S1
					      ;  K5=CW[2]*C1+SW[2]*S1
					      ;  K3_SUN=SW[2]*C1-CW[2]*S1
					  ;}
					   RO=ATM(HP,F135,FT,F0,KP,DT0+T)
					;  B=-0.25*(RZ/(A*H1*H1)*C20+2*ECJ)*RZ*Z/HODH
					;  EATM=2*KB0*RO*A*K4*exp(B*CW[2])
					;  C3=0.25*B*B
					;  FR=1+C3
					;  C2=B*(1+0.125*B*B)
					;  EK[0]=FR-0.5*K2*C2*(CW[2]*CK+CW[2]*SK)
					;  EK[1]=K1*(FR*C1-0.5*C2*K5)
					;  EK[2]=FR*K2*CK-C2*CW[2]+0.5*C3*K2*(1-2*K5*K5)
					;  EK[3]=K1*(0.5*C3*(CW[4]*C1-SW[4]*S1)-0.5*C2*(CW[2]*C1-SW[2]*S1))
					;  EK[4]=-0.5*K2*C2*(CW[2]*CK-SW[2]*SK)+C3*CW[4]
					;  for  (J=0;J<=4;J++) { L1-=A*EATM*KA[0][J]*EK[J]
								 ;L2-=EATM*KA[1][J]*EK[J]
								;}
					   EK[0]=SW[2]
					;  EK[1]=K1*(FR*S1-0.5*C2*K3_SUN)
					;  EK[2]=FR*K2*SK-C2*SW[2]+C3*K2*K3_SUN*K5
					;  EK[3]=K1*(0.5*C3*(SW[4]*C1+CW[4]*S1)-0.5*C2*(SW[2]*C1+CW[2]*S1))
					;  EK[4]=-0.5*K2*C2*(SW[2]*CK+CW[2]*SK)+C3*SW[4]

					;  for  (J=0;J<=4;J++)
						{
						L4-=EATM*KA[3][J]*EK[J] ;
						L6-=EATM*KA[5][J]*EK[J] ;
						 }
		  /*Учет скорости сверхвращения атмосферы */
					;  C1=(HP<200||HP>400) ? 1.0 : 1.0  /* (1.0+0.0025*(HP-200))*/
					;  L3-=C1*EATM*FR*(KA[2][0]+KA[2][1]*CW[2])
					;  L5-=EATM*FR*KA[4][0]*SW[2]


				    ;}

#endif
			}
			 L1*=N
		      ;  L2*=N
		      ;  L3*=N
		      ;  L4*=N
		      ;  L5*=N
		      ;  L6*=N;

/* Возмущения обусловленные неинерциальностью систем координат
//                  if (P_MOVE_CK)
                  {
                  double t,l_moon,l_sun,N_moon,DN_moon,n_moon,n_sun,L_moon,
                  L_sun,k_sec,t_since_1950,de1_e0,e1_e0,dpsi,psi,TE,
                  d_teta_sina,d_teta_cosa,teta_sina,teta_cosa;
// l_moon - средняя аномалия Луны;  /
// 1_sun- средняя аномалия Солнца; /
// L_moon - средняя долгота Луны;  /
// L_sun- средняя долгота Солнца; /
// n_moon-среднее движение Луны;
// n_sun-среднее движение Солнца ;
// F - средний аргумент широты Луны;
// D- разность средних долгот Луны и Солнца;
// N_moon - средняя долгота восходящего узла орбиты Луны на эклиптике.
// Разложения фундаментальных аргументов имеют вид:
// 15341 -это число суток от 31.12.1957 до 01.01.2000
                    t=(DT0+T-15341)/36525;
					l_moon =2.355548393 +t*(8328.69142288+t*(1.517952e-4 +t*3.103e-7));
					l_sun = 6.24003594 + t*(628.30195602*t - t*(2.7974e-6+ t*5.82e-8));
					N_moon = 2.182438624 + t*(- 33.757045936 +t*(3.61429e-5+ 3.88e-8*t));
                    DN_moon=  (- 33.757045936 +t*(2*3.61429e-5 +3*3.88e-8*t))/36525.0;
                    n_moon = 13.1763965/RADIAN;
                    n_sun =  0.9856262833676/RADIAN; //??
             		TE=(21183.5+DT0+T)/36525.0;
                    L_moon= 4.719966501865+(8399.70915755+ 3.488127776573e-5*TE)*TE;      //???
  					L_sun=6.2565835+628.30193*TE;
                    k_sec=1/(3600*RADIAN);//Число радиан в 1 секунде
                    t_since_1950=DT0+T+2921;
                    de1_e0=k_sec*(-9.21*DN_moon*sin(N_moon)+0.18*DN_moon*sin(2*N_moon)-
                           1.1*n_sun*sin(2*L_sun)-0.18*n_moon*sin(2*L_moon)-0.001281);
                    e1_e0 =k_sec*(9.21*cos(N_moon)-0.09*cos(2*N_moon)+
                           0.55*cos(2*L_sun)+0.09*cos(2*L_moon)-0.001281*t_since_1950);
                    dpsi  =k_sec*(-17.24*DN_moon*cos(N_moon)+0.42*DN_moon*cos(2*N_moon)-
                           2.54*n_sun*cos(2*L_sun)+0.13*n_sun*cos(l_sun)-
                           0.4*n_moon*cos(2*L_moon)+0.07*n_moon*cos(l_moon)+0.1379146);
                    psi   =k_sec*(-17.24*sin(N_moon)+0.21*sin(2*N_moon)-
                           1.27*cos(2*L_sun)+0.13*sin(l_sun)-
                           0.2*sin(2*L_moon)+0.07*cos(l_moon)+0.1379146*t_since_1950);
                    d_teta_sina=0.9175*sin(psi)*de1_e0+0.3979*cos(psi)*dpsi;
                    d_teta_cosa=-(0.1583+0.8418*cos(psi))*de1_e0+0.3651*sin(psi)*dpsi;
                    teta_sina=(0.3979+e1_e0)*sin(psi);
                    teta_cosa=0.3651*(1-cos(psi))-e1_e0;
                    L3+=0.5*cos(0.5*I)*(-d_teta_cosa*cos(DBY)-d_teta_sina*sin(DBY));
                    L5+=-cos(I)/(0.5*cos(0.5*I))*(d_teta_sina*cos(DBY)-d_teta_cosa*sin(DBY))+
                        0.5*sin(0.5*I)*(d_teta_sina*teta_cosa-d_teta_cosa*teta_sina);
                    DM=tan(0.5*I)*(d_teta_sina*cos(DBY)-d_teta_cosa*sin(DBY))+
                       0.5*(d_teta_sina*teta_cosa-d_teta_cosa*teta_sina);
                    L6+=DM;
                    L4+=E*DM;

                  }
*/




		      ;  if  (PATMOSFERA && (L1+1.0E-40<0))
             {
			  if  (A*E<HODH*20)
				H1=-HODH/L1;
			      else H1=-A*E/L1;
              }
			      else H1=80000.0;
		      ;  if(DTBAKB==1)
			 {
			 if(T==0) PADEHIE=1.5*L1*sqr(PI2/N)/A;
			 }
			 else

			 PADEHIE=1.5*L1*sqr(PI2/N)/A
                           ;  if (PZADACHI==1&&NPRABCHACTEI==0&&PADEHIE!=0)
                               {
                               C1=0.5*PADEHIE*TK*TK/(PER*PER);
                               TK+=C1;
                               DTK+=C1;
                               }
		      ;  LIFETIME=DT0+T+H1
		      ;  E2[0]=L1
		      ;  E2[1]=(L2*E1[1]-L4*E1[2])/E
		      ;  E2[2]=(L2*E1[2]+L4*E1[1])/E
		      ;  E2[3]=(L3*E1[3]-L5*E1[4])/SI2
		      ;  E2[4]=(L3*E1[4]+L5*E1[3])/SI2
		      ;  E2[5]=L6
		      ;  NPRABCHACTEI++
		      ;  return(0)
		  ;}
/* ***************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  void JVC_DLL close_prognoz(void)
  {
   if (MAEIKP!=NULL) fclose(MAEIKP);
   if (MULTIINDEX!=NULL) fclose(MULTIINDEX);
   if (MAEIFORCE!=NULL) fclose(MAEIFORCE);
   if (MREZ!=NULL) fclose(MREZ);
   if (MEMORY!=NULL) fclose(MEMORY);
   MAEIKP=NULL;
   MULTIINDEX=NULL;
   MAEIFORCE=NULL;
   MREZ=NULL;
   MEMORY=NULL;
  }
