package org.ms2ms.utils;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by wyu on 4/26/14.
 */
public class Stats
{
  static private Map<Long, Double> sLnFactorials = new HashMap<Long, Double>();

  static
  {
    for (long i=0l; i<18l; i++) sLnFactorials.put(i, Math.log(factorial(i)));
  }
  public static double mean(Collection<Double> s)
  {
    if (!Tools.isSet(s)) return 0;

    double avg = 0d;
    for (Double v : s) avg+=v;
    return avg/(double )s.size();
  }
  public static double stdev(Collection<Double> s) { return Math.sqrt(variance(s)); }
  public static double variance(Collection<Double> s)
  {
    if (!Tools.isSet(s) || s.size()==1) return 0;

    double avg=mean(s), var=0d;
    for (Double v : s) { var+= (v-avg)*(v-avg); }

    return var/((double )s.size()-1d);
  }
  public static long factorial(long n)
  {
    long prod = 1l;
    for (long k=1; k<=n; ++k) prod *= k;
    return prod;
  }
  public static double ln_combination(long n, long k) { return ln_factorial(n)-ln_factorial(k)-ln_factorial(n-k); }
  public static double ln_factorial(long n)
  {
    if (n>18) return 0.5d*Math.log(2d*(double )n*3.14) + (double )n*Math.log((double )n) - (double )n;
    return sLnFactorials.get(n);
  }
  public static double hypergeometricPval1(long success, long trials, long success_population, long population)
  {
    long min_trial_pop = trials<success_population?trials:success_population;
    double t1 = (double )trials-(double )population+(double )success_population;

    if (success>min_trial_pop || (double )success<t1 || trials>population || success_population>population) return 1;

    double ln_pop_trials   =ln_combination(population, trials),
           lnfac_suc_pop   =ln_factorial(success_population),
           lnfac_pop_sucpop=ln_factorial(population-success_population), prob=0d;
    for (long suc=success; suc<=min_trial_pop; suc++)
    {
      double p=lnfac_suc_pop - ln_factorial(success) - ln_factorial(success_population-success) +
               lnfac_pop_sucpop - ln_factorial(trials-success) - ln_factorial(population-success_population-trials+success) - ln_pop_trials;
      prob += Math.exp(p);
    }
    return Math.log(prob);
  }
}
