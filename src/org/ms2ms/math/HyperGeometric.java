package org.ms2ms.math;

public class HyperGeometric
{
	public final double DEFAULT_RET_VAL = 99999;
	public static double factorial(long n) {
		double prod = 1;
		for (long k = 1; k <= n; ++k) prod *= k;
		return prod;
	}

	public static double ln_factorial(long n) {
	  double prod;

	  if (n > 15) prod = 0.5 * Math.log((double )(2 * n * Math.PI)) + n * Math.log((double )n) - n;
	  else        prod = Math.log((double )factorial(n));

	  return prod;
	}

	public static double ln_combination(long n, long k) {

	  double comb = ln_factorial(n) - ln_factorial(k) - ln_factorial(n-k);
	  return comb;
	}
	/**  score  = ln_hypergeometric(matched_pairs, n_obs, n_expected_ions, n_bins);   where
	     n_bins =  mz_range / mz_tolerance;
	 *
	 */
	public static double ln_hypergeometric(long n_success,
	                         long n_trials,
	                         long n_sucess_population,
	                         long n_population) {

//		 check for the validity of the arguments per Excel's way
	  double t1 = (double )n_trials - (double )n_population + (double )n_sucess_population;
	  if ( n_success > Math.min(n_trials, n_sucess_population) ||
	      (double )n_success < t1                                      ||
	      n_trials            > n_population                          ||
	      n_sucess_population > n_population) 
		  return Double.MAX_VALUE; // a bogus numbers

	  double p =  ln_combination(n_sucess_population, n_success) +
	              ln_combination(n_population - n_sucess_population, n_trials - n_success) -
	              ln_combination(n_population, n_trials);
	  return p;
	}
	
	// probability of having n_success or more
	public static double hypergeometric_oneSidedPValue(long n_success,
            long n_trials,
            long n_sucess_population,
            long n_population) {

		//check for the validity of the arguments per Excel's way
		double t1 = (double )n_trials - (double )n_population + (double )n_sucess_population;
		if ( n_success > Math.min(n_trials, n_sucess_population) ||
			(double )n_success < t1                                      ||
			n_trials            > n_population                          ||
			n_sucess_population > n_population) 
		return Double.MAX_VALUE; // a bogus numbers
		
		double totalProb=0.0;
		for(long success=n_success; success<= Math.min(n_trials, n_sucess_population); success++){
			
			double p =  ln_combination(n_sucess_population, success) +
			 ln_combination(n_population - n_sucess_population, n_trials - success) -
			 ln_combination(n_population, n_trials);			
			
			totalProb+=Math.exp(p);
		}
		return totalProb;
	}
}