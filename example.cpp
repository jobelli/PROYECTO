#include "hoko.h"
#include <math.h>
#include "aux_fun.h"
extern int Yes_return;

#define EXAMPLE_PI 3.1415926535898
#define EARTH_RADIUS 6378.137 /* This is the WGS84 value in km. */
#define MJD_AT_DEC_31_1957 36203

#define ALL_OUTPUT 1
#define XYZ_OUTPUT 2
#define KEPLER_MEAN_OUTPUT 3
#define KEPLER_OSC_OUTPUT 4

int main (void)
{
  int L, NTOCH;
  double MIZ[6];
  struct bxprog BX;
  struct bixprog * BIX;
  double * TKNK;
  double epoch_in, inclination_deg, ra_asc_node_deg,
    eccentricity, arg_perigee_deg, mean_anomaly_deg, 
    semimajor_axis_earth_radii, deg2rad, rad2deg, nu, snu, cnu, ballistic_coef;
  int object_in, rev_num_in, I,J,i,j;  
  char input_file[256];
  int NUM_STEPS_TO_PROPAGATE;
  double STEP_SIZE_IN_DAYS;
  int output_type;

  rad2deg = 180 / EXAMPLE_PI;
  deg2rad = EXAMPLE_PI / 180;

  /* Read the data from the input.txt file. */
  snprintf(input_file, 256, "%s", "input.txt");

  if (get_elset_input(input_file,
		      &epoch_in,
		      &object_in,
		      &rev_num_in,
		      &STEP_SIZE_IN_DAYS,
		      &NUM_STEPS_TO_PROPAGATE,
		      &semimajor_axis_earth_radii,
		      &inclination_deg,
		      &ra_asc_node_deg,
		      &eccentricity,
		      &arg_perigee_deg,
		      &mean_anomaly_deg,
		      &ballistic_coef,
		      &output_type) == 1)
    {
      printf("Error reading input from %s.\n");
      return 1;
    }

  /* Allocate memory for the TKNK and BIX arrays. */
  if ((TKNK = (double *) calloc(NUM_STEPS_TO_PROPAGATE, sizeof(double))) == NULL)
    {
      printf("Cannot Allocate space for TKNK array.\n");
      return 1;
    }

  if ((BIX = (bixprog *) calloc(NUM_STEPS_TO_PROPAGATE, sizeof(bixprog))) == NULL)
    {
      printf("Cannot Allocate space for BIX array.\n");
      return 1;
    }

  /* Populate the MIZ array with test data. */
  MIZ[0] = semimajor_axis_earth_radii * EARTH_RADIUS;
  MIZ[1] = eccentricity;
  MIZ[2] = inclination_deg * deg2rad;
  MIZ[3] = arg_perigee_deg * deg2rad;
  MIZ[4] = ra_asc_node_deg * deg2rad;
  MIZ[5] = mean_anomaly_deg * deg2rad;

  if (output_type == ALL_OUTPUT)
    {

      printf("These are the original Keplerian Elements in radians.\n");
      for (i = 0; i < 6; i++)
	{
	  printf("%lf\n", MIZ[i]);
	}
      printf("\n");
    }

  /************************************************************
   Prepare for the first call to PROGNOZ.
  *************************************************************/

  /* Populate TKNK which stores epochs. */
  for (I = 0; I < 1; I++)
    {
      TKNK[I] =  epoch_in - MJD_AT_DEC_31_1957;
    }

  /* Initial orbit determination as first measurement. */
  BX.DT       =  TKNK[0];
  BX.NO       =  object_in;
  BX.PZADACHI =  0;
  BX.POCK     =  -2;
  BX.DT       =  TKNK[0];
  BX.U        =  arg_lat(MIZ[1], MIZ[5], MIZ[3]);
  BX.A        =  MIZ[0];
  BX.I        =  MIZ[2];
  BX.DBY      =  MIZ[4];
  BX.L        =  MIZ[1]*cos(MIZ[3]);
  BX.H        =  MIZ[1]*sin(MIZ[3]);
  BX.KB       =  ballistic_coef;

  L           =  PROGNOZ(1,&BX,BIX,TKNK);

  /* Check for errors. */
  if (L > 0) 
    {
      printf("Error return %d from PROGNOZ.\n", L);
    }
  else if (output_type == ALL_OUTPUT)
    {
      printf("Here are the contents of the output structure from PROGNOZ.\n");
      printf("\n");
      printf("Output Elements --------------------------\n");
      printf("Satellite %ld, Epoch: %lf\n", BIX[0].bix.NO, BIX[0].bix.DT);
      printf("\nECI Vector:\n");
      for (i = 0; i < 6; i++)
	printf("%lf\n",BIX[0].X[i]);

      printf("\nKeplerian Elements:\n");
      for (i = 0; i < 6; i++)
	printf("%lf\n", BIX[0].Y[i]);

      printf("\nNonsingular Mean Elements:\n");
      for ( i = 0; i < 6; i++)
	printf("%lf\n", BIX[0].Z[i]);

      printf("\nEccentricity: %lf\n", BIX[0].E);
      printf("Argument of Perigee, radians: %lf\n", BIX[0].W);
      printf("Apogee Altitude: %lf\n", BIX[0].HA);
      printf("Perigee Altitude: %lf\n", BIX[0].HP);
      printf("Geocentric range, in km: %lf\n", BIX[0].RAD);
      printf("Latitude: %lf\n", BIX[0].ALTITUDE);
      printf("Longitude: %lf\n", BIX[0].LONGITUDE);
    }

  /**********************************************************
   Begin subsequent calls to PROGNOZ.
  **********************************************************/

  /* Populate TKNK[] with additional epochs. */
  for (I = 0; I < NUM_STEPS_TO_PROPAGATE; I++)
    {
      TKNK[I] = epoch_in - MJD_AT_DEC_31_1957 + STEP_SIZE_IN_DAYS * I;
    }

  /* Populate the input structure. */
  BX.DT       =  TKNK[0];
  BX.NO       =  object_in;
  BX.PZADACHI =  0;
  BX.POCK     =  -2;
  BX.DT       =  TKNK[0];
  BX.U        =  arg_lat(MIZ[1], MIZ[5], MIZ[3]);
  BX.A        =  MIZ[0];
  BX.I        =  MIZ[2];
  BX.DBY      =  MIZ[4];
  BX.L        =  MIZ[1]*cos(MIZ[3]);
  BX.H        =  MIZ[1]*sin(MIZ[3]);
  BX.KB       =  ballistic_coef;

  L           =  PROGNOZ(NUM_STEPS_TO_PROPAGATE,&BX,BIX,TKNK);

  /* Check for errors. */
  if (L > 0) 
    {
      printf("Error return %d from PROGNOZ.\n", L);
    }
  else if (output_type == ALL_OUTPUT)
    {
      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("\nHere are the contents of the output structure from PROGNOZ.\n");
	  printf("\n");
	  printf("Output Elements --------------------------\n");
	  printf("Satellite %ld, Epoch: %lf\n", BIX[i].bix.NO, BIX[i].bix.DT);
	  printf("\nECI Vector:\n");
	  for (j = 0; j < 6; j++)
	    printf("%lf\n",BIX[i].X[j]);

	  printf("\nKeplerian Elements:\n");
	  for (j = 0; j < 6; j++)
	    printf("%lf\n", BIX[i].Y[j]);

	  printf("\nNonsingular Mean Elements:\n");
	  for (j = 0; j < 6; j++)
	    printf("%lf\n", BIX[i].Z[j]);

	  printf("\nEccentricity: %lf\n", BIX[i].E);
	  printf("Argument of Perigee, radians: %lf\n", BIX[i].W);
	  printf("Apogee Altitude: %lf\n", BIX[i].HA);
	  printf("Perigee Altitude: %lf\n", BIX[i].HP);
	  printf("Geocentric range, in km: %lf\n", BIX[i].RAD);
	  printf("Latitude: %lf\n", BIX[i].ALTITUDE);
	  printf("Longitude: %lf\n", BIX[i].LONGITUDE);
	}
    }
  else if (output_type == XYZ_OUTPUT)
    {
      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("%09.9lf %09.9lf %09.9lf %09.9lf\n", 
		 BIX[i].bix.DT,
		 BIX[i].X[0], 
		 BIX[i].X[1], 
		 BIX[i].X[2]);
	}
    }
  else if (output_type == KEPLER_MEAN_OUTPUT)
    {
      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("%09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf\n",
		 BIX[i].bix.DT,
		 BIX[i].Y[0],
		 BIX[i].Y[1],
		 BIX[i].Y[2],
		 BIX[i].Y[3],
		 BIX[i].Y[4],
		 BIX[i].Y[5]);
        } 
    }
  else if (output_type == KEPLER_OSC_OUTPUT)
    {
      /*       printf("%09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf\n",
	     TKNK[0] + MJD_AT_DEC_31_1957,
	     MIZ[0],
	     MIZ[2] * rad2deg,
	     MIZ[4] * rad2deg,
	     MIZ[1]*cos(MIZ[3]),
	     MIZ[1]*sin(MIZ[3]),
	     arg_lat(MIZ[1], MIZ[5], MIZ[3]) * rad2deg); */

      for (i = 0; i < NUM_STEPS_TO_PROPAGATE; i++)
	{
	  printf("%09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf %09.9lf\n",
		 BIX[i].bix.DT + MJD_AT_DEC_31_1957,
		 BIX[i].bix.A,
		 BIX[i].bix.I * rad2deg,
		 BIX[i].bix.DBY * rad2deg,
		 BIX[i].bix.L,
		 BIX[i].bix.H,
		 BIX[i].bix.U * rad2deg);
	}
    }

  free(BIX);
  free(TKNK);

  return 0;
}




