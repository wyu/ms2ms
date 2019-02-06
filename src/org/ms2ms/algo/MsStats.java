package org.ms2ms.algo;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.stats.Histogram;
import org.expasy.mzjava.stats.HistogramImpl;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * Created by wyu on 4/26/14.
 */
public class MsStats extends Stats
{
  public static Histogram newHistogram(String title, int bins, Range<Double> bound)
  {
    HistogramImpl hist = (bound!=null?new HistogramImpl(bins,
        (bound.upperEndpoint()-bound.lowerEndpoint())/(double )bins, bound.lowerEndpoint()): new HistogramImpl());
    hist.setName(title);

    return hist;
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
  public static <T extends Peak> double median(Collection<T> ys)
  {
    if (ys.size() == 1) return Tools.front(ys).getIntensity();

    double[] pts = new double[ys.size()];
    int    order = 0;
    for (T t : ys) pts[order++] = t.getIntensity();

    return median(pts);
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
  public static String asDeviation(double delta, double m, double max_ppm)
  {
    double ppm = m!=0?1E6*delta/m:Double.MAX_VALUE;
    return (Math.abs(ppm)<max_ppm?(Tools.d2s(ppm,1)+"ppm"):Tools.d2s(delta,4));
  }
}
