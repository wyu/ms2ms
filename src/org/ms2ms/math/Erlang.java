package org.ms2ms.math;

public class Erlang
{	
	public static final double DEFAULT_RET_VAL=99999;
	
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
	
	public static double erlangPValue(double lambda, int k, double x){
		// = IncompleteGamma(k, x*lambda) / Gamma(k)		
		
		// = e^-b . Summation (for i=0 to k-1) (b^i / i!)  ... where b = x*lambda
		
		double result=0;
		double b= x*lambda;
		if(b < 0.00001){return DEFAULT_RET_VAL;}
		for(int i=0; i<k; i++){			
			double log_b_raisedto_i = i*Math.log(b);
			double log_i_factorial = ln_factorial((long) i);
			result+= Math.exp(log_b_raisedto_i - log_i_factorial);
		}
		
		result = result* Math.exp(-1*b);
		return result;		
	}
	
}
