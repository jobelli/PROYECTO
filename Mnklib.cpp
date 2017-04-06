#include "hoko.h"
//#include "GeoForb.h"
//#include <dos.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

FILE *F_doc_vtor;
//extern DT0;
/*
	void JVC_DLL error_print(char *s)
	{
         Application->MessageBox(s,"+°шсър",MB_OK);
	     fclose(F_doc_vtor);
         exit(1);
    }
	   void JVC_DLL task_exit(int pech)
	   { struct tm *tprognoza;
	     time_t delta_t;
	     char str[90];
	   if(getch()=='\x1B')
	   {
              if (Application->MessageBox("¦ртхЁ°шЄ№ чрфрўє ?","-эшьрэшх!",MB_YESNO)==IDYES)
               {
     	         fclose(F_doc_vtor);
	      if(pech)
		{
		 time(&delta_t);
		 tprognoza=gmtime(&delta_t);
		 strftime(str, 90,  "\n %c Принудительное завершение задачи ",tprognoza);
		 fprintf(F_doc_vtor,"\n__________________________________________________________________________________");
		 fprintf(F_doc_vtor,"\%s",str);
		 fprintf(F_doc_vtor,"\n__________________________________________________________________________________");
		}
			 exit(1);
               }
	   }
	   }
    TDateTime YYYYMMDDToTDateTime(int NDMG)
    {
     int     chdnei[12]={0,31,59,90,120,151,181,212,243,273,304,334};
     int I,L1,J,M,K,dday;
      I= (int) (NDMG/10000);
      L1=(int) (NDMG-I*10000);
      I-=1958;
      J=I*365+(I+2)/4;
      M=L1/100-(L1/10000)*100;
      K=L1-(L1/100)*100;
      dday=J+chdnei[M-1]+K;
      if(((I+2)/4==(I+2)/4.0)&&((chdnei[M-1]+K)<=60)) dday-=1;
      if(dday>15400) dday-=1;
      return (dday+21185);
 }
*/

/* ======================================================================
   Destination: matrix inversion.
   Call: P=obmatr(az,jaz,n,a ,ja )
   Input data:
   az is square matrix of dimension jaz x jaz ;
   n is demension of left-top submatrix of matrix az to be inversed;
   Output data: a is square matrix of demension ja x ja, where are placed the inversed matrix;
   P=1 if matrix is singular, else P=0.
  ------------------------------------------------------------------- */
int JVC_DLL obmatr(double *az,int jaz,int n,double *a ,int ja )
 /* обращение произв. матрицы */
       {     double   o,s,*p1,*p2; int      k
			 ,      v
               ,      c
	       ,      i
               ,      j
                 , r[20]    ;
//               ,      *r
//           ; r=(int*)calloc(n,2);
		  for(i=0 ;i<= n-1; i++)
             {     p1=a+i*ja;
		   p2=az+i*jaz;
                   for(j=0 ;j<= n-1; j++)
			 *(p1++)=*(p2++);
             }
	     for(i=0 ;i<= n-1; i++)
                   r[i]=i;
             for(k=0 ;k<= n-1; k++)
             {     o=0;
                   for(i=k ;i<= n-1; i++)
		   {     if    ( (s=fabs(*(a+i*ja+k))) >o )
			 {     o=s;
			       v=i;
			 }
                   }
			if    (o==0)
                   {   /*  puts("\n ******* MATPИЦA BЫPOЖДEHA ******");*/
//                         free(r);
                         return(1);
		   }
                   if    (k<v)
		   {     for(j=0 ;j<= n-1; j++)
                         {     s=*(a+k*ja+j);
					 *(a+k*ja+j)=*(a+v*ja+j);
                               *(a+v*ja+j)=s;
			 }
                         s=r[k];
                         r[k]=r[v];
                         r[v]=s;
                   }
		   *(a+k*ja+k)=1/(*(a+k*ja+k));
		   for(j=0 ;j<= n-1; j++)
                   {     if    (j!=k)
			 {     s=*(a+k*ja+j)**(a+k*ja+k);
                               for(i=0 ;i<= n-1; i++)
			       {     if    (i!=k)
                                           *(a+i*ja+j)-=*(a+i*ja+k)*s;
										 }
                               *(a+k*ja+j)=s;
			 }
                   }
		   for(i=0 ;i<= n-1; i++)
                   {     if    (i!=k)
			       *(a+i*ja+k)*=-*(a+k*ja+k);
                   }
	     }
             for(k=0 ;k<= n-1; k++)
				 {
               while (k<r[k])
		   {     i=r[k];
			 for(j=0 ;j<= n-1; j++)
                         {     s=*(a+j*ja+k);
                               *(a+j*ja+k)=*(a+j*ja+i);
			       *(a+j*ja+i)=s;
                         }
			 c=r[i];
                         r[i]=r[k];
                         r[k]=c;
                   }
		  }
//             free(r);
	     return(0);
       }

/* ================================================================ */
/* ======================================================================
   Destination: transition matrix calculation.
   Call: vpm(E[],P[7][7],DT).
   Input data:E[0:5] is vector of nonsingular elements;
   DT is prediction interval in days;
   Output data: transition matrix of dimension 7x7;
  ------------------------------------------------------------------- */

   void JVC_DLL vpm(double el[],double pme[7][7],double delvre)
/*  ПPOЦEДУPA PACЧETA MATPИЦЫ ЧACTHЫX                    */
/*  ПPOИЗBOДHЫX (ПEPEXOДHOЙ MATPИЦЫ)                     */
/*         << BПM (ЭЛ,ПMЭ,ДEЛBPE)>>                      */
/*    BXOД: ЭЛ[0:6]- BEKTOP ЛAГPAHЖEBЫX ЭO               */
/*          ДEЛBPE - ИHTEPBAЛ ПPOГHOЗ-Я                  */
/*    BЫXOД: ПMЭ[6,6]-MATPИЦA ЧACTHЫX                    */
/*                    ПPOИЗBOДHЫX                        */
       {     int   i,j;
	     double   el0[6]
	       ,   am
	       ,   ekc
	       ,   nakl
	       ,   ash
	       ,   ash2
		   ,   pe
	       ,   n
           ,   n1
	       ,   k0
	       ,   k1
	       ,   k2
	       ,   k3
	       ,   k4
	       ,   k5
	       ,   k6
	       ,   k7
	       ,   k8
	       ,   k9
	       ,   k10
		   ,   z,
           J,
           kk,
           t,
           eps,
           N_moon,
           sin2I,
           n_moon_sqr_to_n,
           n_sun_sqr_to_n,
           sin2eps
		; LBK(el,el0);
	     am=el0[0];
	     ekc=el0[1];
	     nakl=el0[2];
	     ash=1-ekc*ekc;
	     ash2=sqrt(ash);
	     pe=am*(1-ekc*ekc);
	     n1=sqrt(GME/(am*am*am))*86400.0;
   	     n=n1*delvre;
		  k0=C20*n*(sqr(RZ/pe));
	     z=sqr(sin(nakl));
	     k1=1.5*k0;
	     k2=k1*cos(nakl);
	     k3=0.75*k0*(5*z-4);
	     k4=0.75*k0*ash2*(3*z-2);
	     k5=k2+k3;
	     k6=k5+k4;
/*
         if (am>20000.0)         //Учет влияния Луны и Солнца
         {
         	J=5.145396/RADIAN;
	        t=(DT0+delvre+21183.5)/36525;
    	    eps=(23.45229444-0.0130125*t)/RADIAN;
        	N_moon = (259.182328+ t*(- 1934.14201 +t*0.00208))/RADIAN;
	        sin2I=0.5*(sqr(sin(J))*(1+sqr(cos(eps)))+2*sqr(sin(eps)*cos(J))+
       		      sin(2*eps)*sin(2*J)*cos(N_moon)-sqr(sin(J)*sin(eps))*cos(2*N_moon));
        	n_moon_sqr_to_n=0.75*sqr(13.1763965/RADIAN)/(n1*ash2)*ML;
         	n_sun_sqr_to_n=0.75*sqr(0.9856262833676/RADIAN)/(n1*ash2);
	        sin2I=1-1.5*sin2I;
            sin2eps=1-1.5*sqr(sin(eps));
    	    kk=(n_moon_sqr_to_n*sin2I+n_sun_sqr_to_n*sin2eps)*delvre;
        	k2-=kk*cos(nakl)*(1+1.5*ekc*ekc);
	        k5+=kk*(2-2.5*z+0.5*ekc*ekc-(1+1.5*ekc*ekc)*cos(nakl));
    	    k6+=kk*(2-2.5*z+0.5*ekc*ekc-(1+1.5*ekc*ekc)*cos(nakl)-
                 1/3*ash2*(1-1.5*z)*(7+3*ekc*ekc));
         }
*/

	     k7=el[1]*cos(k5)-el[2]*sin(k5);
	     k8=el[2]*cos(k5)+el[1]*sin(k5);
	     k9=el[3]*cos(k2)-el[4]*sin(k2);
	     k10=el[4]*cos(k2)+el[3]*sin(k2);

	     for(j=0 ;j<= 5; j++)
	     for(i=0 ;i<= 6; i++)  pme[j][i]=0;
	     pme[0][0]=1;
	     pme[1][0]=3.5*k8*k5/am;
		  pme[1][1]=cos(k5)-4*el[1]*k8*k5/ash;
	     pme[1][2]=-sin(k5)-4*el[2]*k8*k5/ash;
	     pme[1][3]=-el[3]*k8*(20*k2-4*k1);
	     pme[1][4]=-el[4]*k8*(20*k2-4*k1);
	     pme[2][0]=-3.5*k7*k5/am;
	     pme[2][1]=sin(k5)+4*el[1]*k7*k5/ash;
		  pme[2][2]=cos(k5)+4*el[2]*k7*k5/ash;
	     pme[2][3]=el[3]*k7*(20*k2-4*k1);
	     pme[2][4]=el[4]*k7*(20*k2-4*k1);
	     pme[3][0]=3.5*k10*k2/am;
	     pme[3][1]=-4*el[1]*k10*k2/ash;
	     pme[3][2]=-4*el[2]*k10*k2/ash;
	     pme[3][3]=cos(k2)+4*el[3]*k10*k1;
	     pme[3][4]=-sin(k2)+4*el[4]*k10*k1;
	     pme[4][0]=-3.5*k9*k2/am;
	     pme[4][1]=4*el[1]*k9*k2/ash;
	     pme[4][2]=4*el[2]*k9*k2/ash;
	     pme[4][3]=sin(k2)-4*el[3]*k9*k1;
	     pme[4][4]=cos(k2)-4*el[4]*k9*k1;
	     pme[5][0]=-1.5*n/am-3.5*k6/am;
	     pme[5][1]=el[1]*(3*k6+k5)/ash;
	     pme[5][2]=el[2]*(3*k6+k5)/ash;
	     pme[5][3]=el[3]*((20+12*ash2)*k2-4*k1);
	     pme[5][4]=el[4]*((20+12*ash2)*k2-4*k1);
		 pme[5][5]=1;
       }
/* ================================================================ */

	void JVC_DLL D_rad_xyz(double radius,double alfa,double delta,double A[2][3])
		 {
        double r,ca,sa,cd,sd;
        r=1/radius;
	ca=cos(alfa)*r;
	sa=sin(alfa)*r;
	cd=cos(delta);
	sd=sin(delta);
	A[0][0]=-sa;
	A[0][1]=ca;
	A[0][2]=0.0;
	A[1][0]=-ca*sd;
	A[1][1]=-sa*sd;
	A[1][2]=cd*r;
       }
/* ================================================================ */
	void JVC_DLL KOPCT(double f,double l,double H,double *L)
       {
	double N,sf,cf,E=1.-ECJ;
	sf=sin(f);
	cf=cos(f);
	N=RZ/sqrt(cf*cf+E*sf*E*sf);
	L[0]=(N+H)*cf*cos(l);
	L[1]=(N+H)*cf*sin(l);
	L[2]=(N*E*E+H)*sf;
       }
/* ================================================================ */
/* ======================================================================
   Destination: the partial derivative matrix calculation.
   Call: D_xyz_alambda(E[6],D[6][6]).
   Input data:E[0:5] is vector of nonsingular elements;
   Output data: the partial derivative matrix D[6][6].
  ------------------------------------------------------------------- */

	void JVC_DLL D_xyz_alambda(double e[6],double D[6][6])
       {
	double
	      CC[3],SS[3],DCP[3],DCQ[3],DSP[3],DSQ[3],
	      A,L,H,Q,P,l,K,E,KCI,ETA,N,NU,R1,T,E0,H1,S,C,C1,
	      S1,C2,S2,CH,CL,SH,SL,DCH,DSH,DCL,DSL,CE,SE;
	int I;
	A=e[0];
	L=e[1];
	H=e[2];
	Q=e[3];
	P=e[4];
	l=e[5];
	N=sqrt(GME/A);
	KCI=sqrt(1-L*L-H*H);
	ETA=sqrt(1-P*P-Q*Q);
	NU=1/(1+KCI);
	T=NU*NU/KCI;
		 for(E0=l,E=l+0.01;fabs(E-E0)>1.E-8; E=l+L*sin(E0)-H*cos(E0)) E0=E;
       CE=cos(E);
       SE=sin(E);
       H1=L*SE-H*CE;
       K=L*CE+H*SE;
       R1=1.0/(1-K);
       C=CE-L+H*H1*NU;
       S=SE-H-L*H1*NU;
       C1=R1*(-SE+H*K*NU);
       S1=R1*(CE-L*K*NU);
       C2=-R1*R1*R1*C;
       S2=-R1*R1*R1*S;
      CH=A*((NU+T*H*H)*H1-(NU*H+C1)*CE);
      SH=A*(-1-T*L*H*H1+(NU*L-S1)*CE);
      CL=A*(-1+T*H*L*H1+(NU*H+C1)*SE);
      SL=A*(-(NU+T*L*L)*H1-(NU*L-S1)*SE);
      DCH=R1*((NU+T*H*H)*K+(NU*H+C1)*SE)-C2*CE;
      DSH=R1*(-T*L*H*K-(NU*L-S1)*SE)-S2*CE;
      DCL=R1*(T*H*L*K+(NU*H+C1)*CE)+C2*SE;
      DSL=R1*(-(NU+T*L*L)*K-(NU*L-S1)*CE)+S2*SE;
      CC[0]=1-2*P*P;
      CC[1]=2*P*Q;
      CC[2]=-2*P*ETA;
      SS[0]=2*P*Q;
		SS[1]=1-2*Q*Q;
      SS[2]=2*Q*ETA;
      DCP[0]=-4*P;
      DCP[1]=2*Q;
      DCP[2]=-2*(ETA-P*P/ETA);
      DCQ[0]=0.0;
      DCQ[1]=2*P;
      DCQ[2]=2*P*Q/ETA;
      DSP[0]=2*Q;
      DSP[1]=0.0;
      DSP[2]=-2*P*Q/ETA;
      DSQ[0]=2*P;
      DSQ[1]=-4*Q;
      DSQ[2]=2*(ETA-Q*Q/ETA);
      for(I=0; I<=2; I++)
      {
      D[I][0]=C*CC[I]+S*SS[I];
      D[I][1]=CL*CC[I]+SL*SS[I];
      D[I][2]=CH*CC[I]+SH*SS[I];
      D[I][3]=A*(C*DCQ[I]+S*DSQ[I]);
      D[I][4]=A*(C*DCP[I]+S*DSP[I]);
      D[I][5]=A*(C1*CC[I]+S1*SS[I]);
      D[I+3][0]=-0.5*N*(C1*CC[I]+S1*SS[I]);
      D[I+3][1]=N*(DCL*CC[I]+DSL*SS[I]);
		D[I+3][2]=N*(DCH*CC[I]+DSH*SS[I]);
      D[I+3][3]=N*(C1*DCQ[I]+S1*DSQ[I]);
      D[I+3][4]=N*(C1*DCP[I]+S1*DSP[I]);
		D[I+3][5]=N*(C2*CC[I]+S2*SS[I]);
		}
	  }
/* ================================================================ */
	void JVC_DLL REDUCTION(double alfa,double delta,double JD1957,int Epocha, int Psr,
			 double *d_alfa,double *d_delta)
		 {
	 /*
		 Процедура пересчета угловых координат из эпохи 1950г
		 в квазиинерциальную систему координат
		 Вход:
		 alfa-прямое восхождение на эпоху 1950 г (2000,JD1957)
		 delta-склонение на эпоху 1950 г (2000,JD1957)
		 Psr-признак
		 Epocha =
		 0-1950
		 1-2000
		 2-tekyschaya
		  Выход:
		 d_alfa-прямое восхождение в квазиинерциальной системе координат
		 d_delta-склонение в квазиинерциальной системе координат

	 */
	double d,T,L,sin_L,cos_L,sin_alfa,cos_alfa,sin_del,cos_del,
					kappa,d_kappa,d_epsilon,k,d_k;

/*Редукция измерений на эпоху 1950 года */
		T=JD1957+2921;
	 L=(12.1128-0.0529539*T)/RADIAN;

	if (Epocha==0)        // Epocha 1950
	{
		kappa=0.054875*T/(3600*RADIAN);
	 d_kappa=-33.3e-6*sin(L)
		 +0.4e-6*sin(2*L)
		 -2.5e-6*sin(2*(280.0812+0.9856473*T)/RADIAN)
		 -0.4e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
	 d_epsilon=44.7e-6*cos(L)
		 -0.4e-6*cos(2*L)
		 +2.7e-6*cos(2*(280.0812+0.9856473*T)/RADIAN)
		 +0.4e-6*cos(2*(64.3824+13.176396*T)/RADIAN);
	 *d_alfa =((kappa+d_kappa)*sin(alfa)-d_epsilon*cos(alfa))*tan(delta);
	 *d_delta=(kappa+d_kappa)*cos(alfa)+d_epsilon*sin(alfa);
	}

			if (Epocha==1) // Epocha 2000
	{
    		T=18263.0;
	 L=(12.1128-0.0529539*T)/RADIAN;
//  	 kappa=0.054875*T/(3600*RADIAN);
     kappa=0;
	 d_kappa=-33.3e-6*sin(L)
		 +0.4e-6*sin(2*L)
		 -2.5e-6*sin(2*(280.0812+0.9856473*T)/RADIAN)
		 -0.4e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
	 d_epsilon=44.7e-6*cos(L)
		 -0.4e-6*cos(2*L)
		 +2.7e-6*cos(2*(280.0812+0.9856473*T)/RADIAN)
		 +0.4e-6*cos(2*(64.3824+13.176396*T)/RADIAN);
	 *d_alfa =((kappa+d_kappa)*sin(alfa)-d_epsilon*cos(alfa))*tan(delta);
	 *d_delta=(kappa+d_kappa)*cos(alfa)+d_epsilon*sin(alfa);
  		 k    =(3.505971E-5+0.290975E-12*T+0.2071E-18*T*T)*T/RADIAN;/*прецессия по прямому восхождению*/
 	  d_k    =-76.7e-6*sin(L)   /*нутация в прямом восхождении*/
		  +0.9e-6*sin(2*L)
		  -5.7e-6*sin(2*(280.00812+0.9856473*T)/RADIAN)
		  -0.9e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
	 *d_alfa+=-k-d_k;

//	 *d_delta=0;

/*
		kappa=0.054875*T/(3600*RADIAN);
	 d_kappa=-33.6e-6*sin(L)
		 +0.4e-6*sin(2*L)
		 -2.5e-6*sin(2*(280.0812+0.9856473*T)/RADIAN)
		 -0.4e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
	 d_epsilon=44.7e-6*cos(L)
		 -0.4e-6*cos(2*L)
		 +2.7e-6*cos(2*(280.0812+0.9856473*T)/RADIAN)
		 +0.4e-6*cos(2*(64.3824+13.176396*T)/RADIAN);
	 *d_alfa =((kappa+d_kappa)*sin(alfa)-d_epsilon*cos(alfa))*tan(delta);
	 *d_delta=(kappa+d_kappa)*cos(alfa)+d_epsilon*sin(alfa);
*/
	}

		if (Epocha==2)    //tekyschaya Epocha
	{
		 k    =(3.505971E-5+0.290975E-12*T+0.2071E-18*T*T)*T/RADIAN;/*прецессия по прямому восхождению*/
	  d_k    =-76.7e-6*sin(L)   /*нутация в прямом восхождении*/
		  +0.9e-6*sin(2*L)
		  -5.7e-6*sin(2*(280.00812+0.9856473*T)/RADIAN)
		  -0.9e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
	 *d_alfa =-k-d_k;
	 *d_delta=0;
	}
/*Учет годичной аберрации */

	if(Psr)
		{
	d=JD1957+21183.5;
	T=d/36525.0;
	L=(297.696678+0.9856473354*d+0.000303*T*T+1.3966288*T)/360;
	sin_L=sin(L);
	cos_L=cos(L);
	sin_alfa=sin(alfa);
	cos_alfa=cos(alfa);
	sin_del=sin(delta);
	cos_del=cos(delta);
	*d_alfa-=0.48481368e-5*(20.47*sin_alfa*sin_L+18.87*cos_L*cos_alfa)/cos_del;
		 *d_delta-=0.48481368e-5*(20.47*sin_del*sin_L*cos_alfa+18.87*cos_L*(0.4336661*cos_del-sin_del*sin_alfa ));
		}

		 }
/* ================================================================ */

	void JVC_DLL KVAZI_TO_T(double X[],double T,int p)
		 {
/*
		 Процедура пересчета прямоугольных координат из квазиинерциальной
	системы координат к систему координат, связанную со средним
	экватором и средней ТВР текущей эпохи Т
       Вход:
       X[]-прямоугольные координаты в квазиинерциальной системе координат
		 T-момент тукущей эпохи
		 Выход:
		 если p=0, тo
		 X[]-прямоугольные координаты в системе координат, связанной с истинным
		  экватором и истинной ТВР текущей эпохи Т
		 если p=1, тo
		 X[]-прямоугольные координаты в системе координат, связанной со средним
		  экватором и средней ТВР текущей эпохи Т
*/
/* ======================================================================
   Destination: transformation of the position and velocity in the
   quasinertial coordinate system to the coordinate system  related to the epoch of date
   Call: KVAZI_TO_T(X,T,P);
   Input data: X[0:5] is the position and velocity vector in the
   quasinertial system coordinate;
   T is the time counted off from 0h UTC Dec 31 1957.
   if P=0 then the output data should be produced in the True Eqinox and True Equator of Date frame,
   else the output data should be produced in the Mean Eqinox and Mean Equator of Date frame.
   Output data:X[0:5] is the position and velocity vector in the
   coordinate system  related to the Epoch of Date.
  ------------------------------------------------------------------- */

	double L,x,y,z,dx,dy,dz,k,dk_dt,d_k,d_kappa,d_epsilon;

/*Редукция измерений на эпоху 1950 года */
      T+=2921;
		k    =(3.505971E-5+0.290975E-12*T+0.2071E-18*T*T)*T/RADIAN;/*прецессия по прямому восхождению*/
  dk_dt    =(3.505971E-5+2*0.290975E-12*T+3*0.2071E-18*T*T)/(RADIAN*86400.0);/*скорость прецессии по прямому восхождению*/
		L    =(12.1128-0.0529539*T)/RADIAN;
	 d_k    =-76.7e-6*sin(L)   /*нутация в прямом восхождении*/
		  +0.9e-6*sin(2*L)
		  -5.7e-6*sin(2*(280.00812+0.9856473*T)/RADIAN)
		  -0.9e-6*sin(2*(64.3824+13.176396*T)/RADIAN);

    x=cos(k+d_k)*X[0]-sin(k+d_k)*X[1];
    y=sin(k+d_k)*X[0]+cos(k+d_k)*X[1];
   dx=cos(k+d_k)*X[3]-sin(k+d_k)*X[4]-dk_dt*y;
	dy=sin(k+d_k)*X[3]+cos(k+d_k)*X[4]+dk_dt*x;
    X[0]=x;
    X[1]=y;
    X[3]=dx;
    X[4]=dy;
    z=X[2];
   dz=X[5];
	 if(p)
    {
    d_kappa=-33.6e-6*sin(L)  /*нутация в склонении*/
	     +0.4e-6*sin(2*L)
	     -2.5e-6*sin(2*(280.00812+0.9856473*T)/RADIAN)
	     -0.4e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
  d_epsilon= 44.7e-6*cos(L) /*нутация в наклоне экватора*/
	     -0.4e-6*cos(2*L)
	     +2.7e-6*cos(2*(280.00812+0.9856473*T)/RADIAN)
	     +0.4e-6*cos(2*(64.3824+13.176396*T)/RADIAN);
    X[0]=         x+      d_k*y+  d_kappa*z;
    X[1]=    -d_k*x+          y+d_epsilon*z;
    X[2]=-d_kappa*x-d_epsilon*y+          z;
    X[3]=         dx+      d_k*dy+  d_kappa*dz;
    X[4]=    -d_k*dx+          dy+d_epsilon*dz;
    X[5]=-d_kappa*dx-d_epsilon*dy+          dz;
    }
 }
/* ================================================================ */
/* ======================================================================
   Destination: transformation of the position and velocity in the
   quasinertial coordinate system to the coordinate system  related to the epoch of date
   Call: T_TO_KVAZI(X,T,P)
   Input data: X[0:5] is the position and velocity vector in the
   coordinate system  related to the Epoch of Date.
   T is the time counted off from 0h UTC Dec 31 1957.
   if P=0 then the input data are presented in the True Eqinox and Equator of Date,
   else the input data are presented in the Mean Equator of Date.
   Output data:X[0:5] is the position and velocity vector in the
   quasinertial system coordinate.
  ------------------------------------------------------------------- */
	void JVC_DLL T_TO_KVAZI(double X[],double T,int p)
		 {
/*
		 Процедура пересчета прямоугольных координат из систем координат, связанную со средним
	экватором и средней ТВР текущей эпохи Т v квазиинерциальнuyu
	системы координат к
       Вход:
       		 если p=0, тo
		 X[]-прямоугольные координаты в системе координат, связанной со истинным
		  экватором и истинной ТВР текущей эпохи Т
		 если p=1, тo
		 X[]-прямоугольные координаты в системе координат, связанной со средним
		  экватором и средней ТВР текущей эпохи Т
      Выход:
       X[]-прямоугольные координаты в квазиинерциальной системе координат
		 T-момент тукущей эпохи
*/
	double L,x,y,z,dx,dy,dz,k,dk_dt,d_k,d_kappa,d_epsilon;

/*Редукция измерений на эпоху 1950 года */
      T+=2921;
		k    =(3.505971E-5+0.290975E-12*T+0.2071E-18*T*T)*T/RADIAN;/*прецессия по прямому восхождению*/
  dk_dt    =(3.505971E-5+2*0.290975E-12*T+3*0.2071E-18*T*T)/(RADIAN*86400.0);/*скорость прецессии по прямому восхождению*/
		L    =(12.1128-0.0529539*T)/RADIAN;
	 d_k    =-76.7e-6*sin(L)   /*нутация в прямом восхождении*/
		  +0.9e-6*sin(2*L)
		  -5.7e-6*sin(2*(280.00812+0.9856473*T)/RADIAN)
		  -0.9e-6*sin(2*(64.3824+13.176396*T)/RADIAN);

    x=X[0];
    y=X[1];
    z=X[2];
    dx=X[3];
    dy=X[4];
    dz=X[5];

	 if(p)
    {
    d_kappa=-33.6e-6*sin(L)  /*нутация в склонении*/
	     +0.4e-6*sin(2*L)
	     -2.5e-6*sin(2*(280.00812+0.9856473*T)/RADIAN)
	     -0.4e-6*sin(2*(64.3824+13.176396*T)/RADIAN);
  d_epsilon= 44.7e-6*cos(L) /*нутация в наклоне экватора*/
	     -0.4e-6*cos(2*L)
	     +2.7e-6*cos(2*(280.00812+0.9856473*T)/RADIAN)
	     +0.4e-6*cos(2*(64.3824+13.176396*T)/RADIAN);
    X[0]=         x-      d_k*y-  d_kappa*z;
    X[1]=     d_k*x+          y-d_epsilon*z;
    X[2]= d_kappa*x+d_epsilon*y+          z;
    X[3]=         dx-      d_k*dy-  d_kappa*dz;
    X[4]=     d_k*dx+          dy-d_epsilon*dz;
    X[5]= d_kappa*dx+ d_epsilon*dy+         dz;
    }

    x=cos(k+d_k)*X[0]+sin(k+d_k)*X[1];
    y=-sin(k+d_k)*X[0]+cos(k+d_k)*X[1];
   dx=cos(k+d_k)*X[3]+sin(k+d_k)*X[4]+dk_dt*y;
   dy=-sin(k+d_k)*X[3]+cos(k+d_k)*X[4]-dk_dt*x;
    X[0]=x;
    X[1]=y;
    X[3]=dx;
    X[4]=dy;

 }
/* ================================================================ */



/* ================================================================ */
	void JVC_DLL KEP_TO_LAP(double A,double I,double w,double DBY,double T,
			double *IL,double *wL,double *DBYL,double *LL)
	 /*
       Процедура пересчета ЭО экваториальной системы координат к систему координат,
       связанную с плоскостью Лапласа
       Вход:
       A,I,DBY,w-наклонение, долгота восходящего узла и аргумент перигея
       в экваториальной системе координат
       T-момент тукущей эпохи
       Выход:
       IL,DBYL,wL,LL-наклонение, долгота восходящего узла, аргумент перигея
        в системе координат, связанной с плоскостью Лапласа и наклон плоскости Лапласа к экватору
     */
     {  double k=9.353e-4,n,eekl,c1,c2;
        n=sqrt(GME/A)/A*86400.0;
   		eekl=0.40931975-0.00022711097*(21183.5+T)/36525.5;
        *LL=0.5*atan(k*sin(2*eekl)/(-2*C20*sqr(RZ/A*n)+k*cos(2*eekl)));
        c1=sin(I)*sin(DBY);
        c2=-cos(I)*sin(*LL)+sin(I)*cos(*LL)*cos(DBY);
		*DBYL=atan2(c1,c2);
        *IL=asin(sqrt(c1*c1+c2*c2));
        c1=sin(*LL)*sin(DBY);
        c2=sin(I)*cos(*LL)-cos(I)*sin(*LL)*cos(DBY);
        *wL=w-atan2(c1,c2);
     }
/* ================================================================ */
void JVC_DLL PROGNOZ_KMO(double EL[],double KMO[7][7],double DELT,double EK[],
    double radius,double alfa,double delta,int P,double C[7][7],double C2[2][2])
    /*
       Процедура прогнозирования КМО ЭО
       Вход:
       EL[6]-вектор Лагранжевых элементов на момент t0
       KMO[49]-КМО определения ЭО КО на момент t0
       EL[6]-вектор Лагранжевых элементов на момент tk
       DELT -интервал прогноза DELT=tk-t0
       radius-наклонная дальность
       alfa - прямое восхождение
       delta-склонение
       P-признак =0, если осуществляется прогноз только КМО ЭО
                 =1, если рассчитывается КМО параметров измерений
		 Выход:
		 C[49]-КМО определения ЭО КО на момент tk
		 C2[2][2]-КМО определения alfa,delta на момент tk
	 */
{
  double A[2][3],PM[7][7],D[6][6],D2[3][6];
  int I,J,L,M;
  vpm(EL,PM,DELT) ;
  /*
  PM[0]+=PQ1[0]/BX.A*DELT;
  PM[30]+=15/8*PQ1[5]*PQ1[0]*DELT*DELT/(EL[0]*EL[0]);
  */
  for(M=0;M<6;M++)
  for(I=M;I<6;I++)
  {
  C[M][I]=0;
  if(I==M)
   {
   if(I<5)
   C[I][I]+=(I>0) ?  1.0e-12*DELT*DELT  : 16.0e-8*DELT*DELT;
   else
   C[I][I]+=1.0e-10*DELT*DELT;
   }
  for(J=0;J<6;J++)
  for(L=0;L<6;L++)
  C[M][I]+=PM[M][J]*KMO[J][L]*PM[I][L];
  C[I][M]=C[M][I];
  }
   if(P)
  {
   D_rad_xyz(radius,alfa,delta,A);
   D_xyz_alambda(EK,D);
  for(M=0;M<2;M++)
  for(I=0;I<6;I++)
  {
  D2[M][I]=0;
  for(J=0;J<3;J++)
  D2[M][I]+=A[M][J]*D[J][I];
  }
  for(M=0;M<2;M++)
  for(I=0;I<2;I++)
  {
  C2[M][I]=0;
  for(J=0;J<6;J++)
  for(L=0;L<6;L++)
  C2[M][I]+=D2[M][J]*C[J][L]*D2[I][L];
  }
  }
}
/* ================================================================ */


   void JVC_DLL GRAD_RAD1(char *gmsds,double *rdr,int k)
/* НАЗНАЧЕНИЕ: При к=1 : Пересчет gmsds из размерности градусы (г),минуты (м),
			 секунды (с),доли секунд (д)[ггг.мм.сс.дд] в rdr [рад]
		   k=0         [рад] - > [ггг.мм.сс.дд]
   Пример    : gmsds = 135.54.39.98 - это 135 град., 54 мин ,39.98 сек .
   !!!! Угол от -90 до 90 градусов
*/
   {         int I,g,m ;
	     double sds;
	     char string[15];
	     if   (k==1)
	     {

	      if(gmsds[0]=='-')
	       I=-1;
	      else
	       I=1;
	      strncpy(string,gmsds+1,2);
	      string[2] = '\0';
	      *rdr=atof(string)/360.0;
	      strncpy(string, gmsds+4,2);
	      string[2] = '\0';
	      *rdr=(*rdr+atol(string)/21600.0);
			strncpy(string, gmsds+7,5);
	      string[6] = '\0';
	      *rdr=I*(*rdr+atof(string)/1296000.0)*PI2;
	     }
	     else
	     {    sds=fabs(*rdr)*RADIAN;
		  g=sds;
		  m=(sds-g)*60;
		  sds=(sds-g)*3600.0-m*60.0;
		  gmsds[0] = '\0';
		  if(*rdr<0) strcat(gmsds,"-"); else strcat(gmsds," ");
		  if(g<10) strcat(gmsds,"0");
//		  itoa(g,string,10);
		  strcat(gmsds,string);
		  strcat(gmsds,".");
//		  itoa(m,string,10);
		  if(m<10) strcat(gmsds,"0");
		  strcat(gmsds,string);
		  strcat(gmsds,".");
		  m=sds;
//		  itoa(m,string,10);
		  if(m<10) strcat(gmsds,"0");
		  strcat(gmsds,string);
		  strcat(gmsds,".");
		  m=(sds-m)*100;
//		  itoa(m,string,10);
		  if(m<10) strcat(gmsds,"0");
		  strcat(gmsds,string);
		  gmsds[12] = '\0';
	     }
   }
/* ================================================================ */
/* ================================================================ */
   void JVC_DLL STR_TIME(double t,char *dmg,char *ch_m_s,int k)
/* НАЗНАЧЕНИЕ: Пересчет t из размерности JD1957ОВ в дмг и часы минуты секунды
Если к=0		    [сутки] - >[ггггммдд] [чч.мм.сс.ссс]
Если к=1		    [сутки] - >[дд/мм/гг] [чч.мм.сс.ссс]
*/
   {    int J; char str2[3];long msec;
	struct tm *tp;
	time_t d_t;
   if(t>4384) t-=4384;
   msec=t*8640000;
   d_t=t*86400+0.0005;
   tp=gmtime(&d_t);
  if (k)
  strftime(dmg,9,"%d/%m/%y",tp);
  else
  strftime(dmg,9,"%Y%m%d",tp);
  strftime(ch_m_s,9,"%H:%M:%S",tp);
  J=msec-d_t*100;
//  itoa(J,str2,10);
  strcat(ch_m_s,".");
  if(J<10) strcat(ch_m_s,"0");
  strcat(ch_m_s,str2);
  ch_m_s[11]='\0';
	}
/* ================================================================ */

/* ================================================================ */
	void JVC_DLL GNSK_V_GPSK(double E0[],double S,int prz,double EK[])
/* НАЗНАЧЕНИЕ: Пересчет XYZVxVyVz из ГНСК  (E0[])
			 в ГПСК(EK[]) при prz=0
			 Пересчет XYZVxVyVz из ГПСК  (E0[])
			 в ГНСК(EK[]) при prz=1

	ГНСК: км, км/с.
	ГПСК: км, км/с.
	S-звездное время в радианах
*/
		 {
		  double CS,SS;
		  int G;
		  G=(prz) ? 1 : -1;
		  CS=cos(G*S);
		  SS=sin(G*S);
		  EK[2]=E0[2];
		  EK[5]=E0[5];
			EK[0]=CS*E0[0]-SS*E0[1];
			EK[1]=SS*E0[0]+CS*E0[1];
			EK[3]=CS*E0[3]-SS*E0[4]-G*CBZ*EK[1]/86400.0;
			EK[4]=SS*E0[3]+CS*E0[4]+G*CBZ*EK[0]/86400.0;
	    }
/* ================================================================ */
/* ================================================================ */
double JVC_DLL CH_M_SToTime(char ch_m_s[])
//-ЁюуЁрььр яхЁхёўхЄр ++.¦¦.--.ёёё т фюыш ёєЄюъ
{
    double t;
    char promchar[15];
	strncpy(promchar,ch_m_s,2);
	promchar[2] = '\0';
	t=atol(promchar)/24.0;
	strncpy(promchar,ch_m_s+3,2);
	promchar[2] = '\0';
	t+=atol(promchar)/1440.0;
	strncpy(promchar,ch_m_s+6,6);
	promchar[6] = '\0';
	t+=atof(promchar)/86400.0;
    return (t);
}
/* ================================================================ */
double JVC_DLL CH_M_SToDelta(char s_delta[])
//-ЁюуЁрььр яхЁхёўхЄр ёъыюэхэш_ -++.¦¦.--.ёёё т Ёрфшррэ_
{
    double delta;
    int I;
    char promchar[15];
    if(s_delta[0]=='-')
	I=-1;
	else
	I=1;
	strncpy(promchar,s_delta+1,2);
	promchar[2] = '\0';
	delta=atof(promchar)/360.0;
	strncpy(promchar,s_delta+4,2);
	promchar[2] = '\0';
	delta=(delta+atol(promchar)/21600.0);
	strncpy(promchar,s_delta+7,6);
	promchar[6] = '\0';
	delta=I*(delta+atof(promchar)/1296000.0)*PI2;
    return (delta);
}
/* ================================================================ */

