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
  /** Estimates the mean and stdev of intensities after 'rounds' of outlier rejections.
   *  The outliers are defined as those outside of the 'stdevs' multiples of the calculated stdev
   *
   * @param intensities are the intensities of the points to be considered
   * @param stdevs is the multiple of stdev that define the boundary for the outliers
   * @param rounds is the number of times the outlier rejections will be attempted.
   * @return a double array where mean is followed by stdev of the intensities excluding the outliers
   */
  public static double[] outliers_rejected(Collection<Double> intensities, double stdevs, int rounds)
  {
    // deal with null or singleton first
    if (intensities==null)    return null;
    if (intensities.size()<2) return new double[] { Tools.front(intensities), 0d};

    double avg=0, bound=Double.MAX_VALUE;

    for (int i = 0; i < rounds; i++)
    {
      Iterator<Double> itr = intensities.iterator();
      while (itr.hasNext())
        if (Math.abs(itr.next()-avg)>bound) itr.remove();

      avg   = mean(intensities);
      bound = stdev(intensities) * stdevs;
    }

    return new double[] {avg, bound};
  }
  public static Peak accumulate(Peak A, Peak B)
  {
    if (A == null || B == null) return A;

    // assuming the check and balance already done entering the call
    double    sum = A.getIntensity()+B.getIntensity(),
      sum_product = A.getIntensity()*A.getMz()+B.getIntensity()*B.getMz();

    A.setMzAndCharge(sum_product / sum);
    A.setIntensity(sum);

    return A;
  }
  public static double median(double[] ys)
  {
    if (ys == null || ys.length == 0) return Double.NaN;
    if (ys.length == 1) return ys[0];
    Arrays.sort(ys);
    if (ys.length % 2 == 0) return (ys[(int )(ys.length * 0.5)    ] +
        ys[(int )(ys.length * 0.5) - 1]) * 0.5;
    return ys[(int )(ys.length * 0.5)];
  }
  public static double median(List<Double> ys)
  {
    if (ys.size() == 1) return Tools.front(ys);
    Collections.sort(ys);
    if (ys.size() % 2 == 0) return (ys.get((int )(ys.size() * 0.5)    ) +
        ys.get((int )(ys.size() * 0.5) - 1)) * 0.5;
    return ys.get((int )(ys.size() * 0.5));
  }
  public static <T extends Peak> double median(Collection<T> ys)
  {
    if (ys.size() == 1) return Tools.front(ys).getIntensity();

    double[] pts = new double[ys.size()];
    int    order = 0;
    for (T t : ys) pts[order++] = t.getIntensity();

    return median(pts);
  }
  public static <T extends Peak> double stdevY(Collection<T> ys) { return Math.sqrt(varianceY(ys)); }
  public static <T extends Peak> double varianceY(Collection<T> ys)
  {
    if (ys.size() == 1) return -1;
    if (ys.size() == 1) return 0;

    double s = 0, ss = 0;
    for (T y : ys) { s += y.getIntensity(); ss += y.getIntensity() * y.getIntensity(); }
    return (ys.size() * ss - s * s) / (ys.size() * (ys.size() - 1));
  }
  public static <T extends Peak> Double meanIntensity(Collection<T> ys)
  {
    if (ys        == null) return null;
    if (ys.size() == 0)    return 0d;

    double s = 0;
    for (T y : ys) { s += y.getIntensity(); }
    return s / ys.size();
  }
  public static <T extends Peak> Double meanMz(Collection<T> ys)
  {
    if (ys        == null) return null;
    if (ys.size() == 0)    return 0d;

    double s = 0;
    for (T y : ys) { s += y.getMz(); }
    return s / ys.size();
  }
  public static <T extends Peak> Double centroid(Collection<T> points)
  {
    return centroid(points, null, null);
  }
  public static <T extends Peak> Double centroid(Collection<T> points, Double x0, Double x1)
  {
    if (! Tools.isSet(points)) return null;

    double sumXY = 0, sumY = 0;
    for (Peak xy : points)
    {
      if ((x0 == null || xy.getMz() >= x0) &&
          (x1 == null || xy.getMz() <= x1))
      {
        sumXY += xy.getMz() * xy.getIntensity();
        sumY  += xy.getMz();
      }
    }
    return sumY != 0 ? sumXY / sumY : null;
  }
  public static double filter(List<Double> A, int index_begin, double[] filters)
  {
    double Y = 0d;
    for (int i = 0; i < filters.length; i++)
      Y += A.get(i + index_begin) * filters[i];

    return Y;
  }
  public static List<Double> smoothBySG5(List<Double> A)
  {
    // Do nothing if the set isn't big enough to smooth.
    if (null==A || A.size()<6) return A;

    // store the smoothed data separately
    List<Double> smoothed = new ArrayList<>(A.size());

    // special handling for the edge points
    smoothed.add(filter(A, 0, new double[] {0.886,0.257,-0.086,-0.143,0.086}));
    smoothed.add(filter(A, 0, new double[] {0.257,0.371,0.343,0.171,-0.143}));

    // the mid section
    for (int i=2; i < A.size()-2; i++)
      smoothed.add(filter(A, i-2, new double[] {-0.086,0.343,0.486,0.343,-0.086}));

    // special handling for the edge points
    smoothed.add(filter(A, A.size()-6, new double[] {-0.143,0.171,0.343,0.371,0.257}));
    smoothed.add(filter(A, A.size()-6, new double[] {0.086,-0.143,-0.086,0.257,0.886}));

    return smoothed;
  }
}
