/**************************************************
 * This header file contains auxiliary routine
 * prototypes for functions needed in order to make
 * a functioning example of the USM propagator in action.
 *******************************************************/

/*******************************************************
 * function arg_lat
 * Arguments: double eccentricity - unitless,
 *            double mean anomaly - radians,
 *            double argument of perigee - radians.
 * output   : double argument of latitude - radians.
 *
 * This routine takes the Eccentricity, Mean Anomaly and
 * Argument of Perigee Keplerian elements to calculate
 * the Argument of Latitude which it returns.
 * This algorithm was taken from:
 * http://home-2.worldonline.nl/%7Esamsvl/satpos.htm
 * This algorithm has not yet been checked.
 *******************************************************/

double arg_lat(double eccentricity_in, 
	       double mean_anomaly_in, 
	       double arg_perigee_in);

/********************************************************
 * function get_elset_input
 * Arguments: char * filename_in - filename from which
 *                                 to get the parameters.
 *            int * epoch_in - Modified Julian Day.
 *            int * object_in - object number.
 *            int * rev_num_in - revolution number.
 *            double * step_size_in_days - fraction of a 
 *                                         day to make each
 *                                         step.
 *            int * number_of_steps - number of steps to 
 *                                    propagate for.
 *            double * semimajor_axis_earth_radii.
 *            double * inclination_deg.
 *            double * ra_asc_node_deg.
 *            double * eccentricity.
 *            double * arg_perigee_deg.
 *            double * mean_anomaly_deg.
 *            double * ballistic_coef.
 *            int    * prmodel. - This is the model switch
 *                                passed to the bxprog 
 *                                structure passed to PROGNOZ.
 *                                (see bxprog structure)
 *            int    * pock. - This is the in/out type switch 
 *                             passed to the bxprog structure
 *                             passed to PROGNOZ.  (see bxprog
 *                             structure)
 *            int    * pzadachi - This is another switch used in
 *                                the bxprog structure.  (see bxprog
 *                                structure)
 *            int    * output_type.
 * output:    1 if error, 0 if success.
 *********************************************************/

int get_elset_input(char * filename_in,
		    double * epoch_in,
		    int * object_in,
		    int * rev_num_in,
		    double * step_size_in_days,
		    int * number_of_steps,
		    double * semimajor_axis_earth_radii,
		    double * inclination_deg,
		    double * ra_asc_node_deg,
		    double * eccentricity,
		    double * arg_perigee_deg,
		    double * mean_anomaly_deg,
		    double * ballistic_coef,
		    int * prmodel,
		    int * pock,
		    int * pzadachi,
		    int * output_type);
