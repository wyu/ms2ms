package org.ms2ms.math;

public class Binomial
{
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
	
	// probability of k successes in n trials, given background probability p
	public static double binomialPDF(long n, double p, long k){
		double logPDF=0;
		if(k>n)return 0;
		if(k==0) return 0;
		
		logPDF = ln_combination(n,k) + k*Math.log(p) + (n-k)*Math.log(1-p);
				
		return Math.exp(logPDF);
	}
	
	// probability of k or less successes in n trials, given background probability p
	public static double binomialCDF(long n, double p, long k){
		double cdf=0;
		if(k>n)return 1.0;
		if(k==0) return 0;
		for(int i=1; i<=k; i++){
			cdf+=binomialPDF(n,p,k);
		}
		return cdf;
	}
}