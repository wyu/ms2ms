package org.ms2ms.utils;

import com.google.common.collect.*;
import com.google.common.primitives.Chars;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.math.NumberUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.stats.Histogram;
import org.expasy.mzjava.stats.HistogramImpl;

import java.util.*;

/**
 * Created by wyu on 4/26/14.
 */
public class Stats
{
  public enum Aggregator { MEAN, MEDIAN, STDEV, COUNT }

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
  public static double mean(double[] s, int up_to)
  {
    if (!Tools.isSet(s)) return 0;

    double avg = 0d;
    for (int i=0; i<up_to; i++) avg+=s[i];
    return avg/(double )s.length;
  }
  public static double stdev(Collection<Double> s) { return Math.sqrt(variance(s)); }
  public static double stdev(double[] s, int up_to) { return Math.sqrt(variance(s, up_to)); }
  public static double variance(Collection<Double> s)
  {
    if (!Tools.isSet(s) || s.size()==1) return 0;

    double avg=mean(s), var=0d;
    for (Double v : s) { var+= (v-avg)*(v-avg); }

    return var/((double )s.size()-1d);
  }
  public static double variance(double[] s, int up_to)
  {
    if (!Tools.isSet(s) || s.length==1) return 0;

    double avg=mean(s, up_to), var=0d;
    for (int i=0; i<up_to; i++) { var+= (s[i]-avg)*(s[i]-avg); }

    return var/((double )s.length-1d);
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
  public static Number aggregate(Collection data, Aggregator func)
  {
    if (!Tools.isSet(data)) return 0;
    if (func.equals(Aggregator.COUNT)) return data.size();

    Collection<Double> ns = new ArrayList<Double>();
    for (Object d : data) ns.add(toDouble(d));

    if (func.equals(Aggregator.MEAN)) return mean(ns);
    //else if (func.equals(Dataframes.Func.MEAN)) return mean(ns);

    return 0;
  }
  public static boolean isNumeric(Object s)
  {
    if (s==null) return false;
    if (s instanceof Double || s instanceof Float || s instanceof Integer || s instanceof Long) return true;
    return NumberUtils.isNumber(s.toString());
  }
  // convert the Object to Number if possible
  public static Object toNumber(Object s)
  {
    try
    {
      // quotes? must remain a string if so
      if (s instanceof String)
      {
        String val = Strs.trim((String )s);
        if ((val.charAt(0)=='"'  && val.charAt(val.length()-1)=='"') ||
            (val.charAt(0)=='\'' && val.charAt(val.length()-1)=='\'')) return val.substring(1, val.length()-1);

        boolean isNum = (val.charAt(0)>='0' && val.charAt(0)<='9');
        return isNum && val.indexOf('.')>=0 ? NumberUtils.createDouble(val) : (isNum?NumberUtils.createLong(val):val);
      }
    }
    catch (Exception e) {}
    return s;
  }
  public static Double toDouble(Object s)
  {
    if (s==null) return null;
    if      (s instanceof String)  return NumberUtils.createDouble((String )s);
    else if (s instanceof Double)  return (Double  )s;
    else if (s instanceof Float )  return ((Float  )s).doubleValue();
    else if (s instanceof Long  )  return ((Long   )s).doubleValue();
    else if (s instanceof Integer) return ((Integer)s).doubleValue();

    return null;
  }
  public static Long toLong(Object s)
  {
    if (s==null) return null;
    if      (s instanceof String)  return NumberUtils.createLong((String) s);
    else if (s instanceof Long  )  return ((Long   )s);
    else if (s instanceof Integer) return ((Integer )s).longValue();

    return null;
  }
  public static Range<Double> closed(double[] s)
  {
    if (s==null) return null;
    double lower=Double.MAX_VALUE, upper = Double.MAX_VALUE*-1;
    for (double x : s)
    {
      if (x<lower) lower=x;
      if (x>upper) upper=x;
    }
    return Range.closed(lower, upper);
  }
  /**
   *
   * @param Xs is the variable
   * @param xs and ys represent the X-Y curve on which the interpolation will be based
   * @param bandwidth is the fraction of source points closest to the current point is taken into account for computing
   *                  a least-squares regression when computing the loess fit at a particular point. A sensible value is
   *                  usually 0.25 to 0.5, the default value is 0.3.
   * @return
   */
  public static double[] interpolate(double[] xs, double[] ys, double bandwidth, double... Xs)
  {
    if (!Tools.isSet(xs) || !Tools.isSet(ys) || xs.length!=ys.length) return null;

    try
    {
      // average the values with duplicated x
      Multimap<Double, Double> xy = TreeMultimap.create();
      for (int i=0; i<xs.length; i++) xy.put(xs[i], ys[i]);
      int i=0; xs=new double[xy.keySet().size()]; ys=new double[xy.keySet().size()];
      for (Double x : xy.keySet()) { xs[i]=x; ys[i]=Stats.mean(xy.get(x)); i++; }
      Tools.dispose(xy);

      double[]                   Ys = new double[Xs.length];
      Range<Double>           bound = Stats.closed(xs);
      PolynomialSplineFunction poly = new LoessInterpolator(bandwidth, 2).interpolate(xs, ys);
      // compute the interpolated value
      for (i=0; i<Xs.length; i++)
      {
        double x=Xs[i];
        if      (x<bound.lowerEndpoint()) x=bound.lowerEndpoint();
        else if (x>bound.upperEndpoint()) x=bound.upperEndpoint();
        Ys[i] = poly.value(x);
      }

      return Ys;
    }
    catch (MathException me)
    {
      throw new RuntimeException("Not able to interpolate: ", me);
    }
  }
  public static double[] sum(double[]... ys)
  {
    if (!Tools.isSet(ys)) return null;

    double[] out = new double[ys[0].length];
    for (int i=0; i<ys.length; i++)
    {
      for (int j=0; j<out.length; j++)
      {
        out[j] += ys[i][j];
      }
    }
    return out;
  }
  public static Histogram newHistogram(String title, int bins, Range<Double> bound)
  {
    HistogramImpl hist=new HistogramImpl(bins, bound!=null?bound.lowerEndpoint():0d, bound!=null?bound.upperEndpoint():1d);
    hist.setName(title);

    return hist;
  }
  public static Integer[] newIntArray(int start, int end)
  {
    Integer[] out = new Integer[end-start];
    for (int i=start; i<end; i++) out[i-start]=i;

    return out;
  }
  public static <T extends Peak> Peak upperOutliers(Collection<T> A, double stdev_multiple, int rounds)
  {
    int      ys_i = 0;
    double[] ys   = new double[A.size()];
    double    avg = 0, bound = Double.MAX_VALUE;

    for (T t : A) t.setMzAndCharge(Math.abs(t.getMz()));

    for (int i = 0; i < rounds; i++)
    {
      ys_i = 0;
      for (T t : A) if (t.getMz() >= 0d && Math.abs(t.getIntensity() - avg) <= bound) ys[ys_i++] = t.getIntensity();

      avg   = mean( ys, ys_i);
      bound = stdev(ys, ys_i) * stdev_multiple;
      // are you converge yet?

      // remove the outlier from the consideration
      for (T t : A)
        if (t.getIntensity() - avg > bound)
          t.setMzAndCharge(Math.abs(t.getMz()) * -1d);
    }

    for (T t : A) t.setMzAndCharge(Math.abs(t.getMz()));

    return new Peak(avg, bound);
  }
}
