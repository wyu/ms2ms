package org.ms2ms.alg;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import com.google.common.collect.TreeMultimap;
import org.apache.hadoop.hbase.util.Hash;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.RetentionTime;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeList;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/23/14
 * Time: 5:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class Spectra
{
  public static boolean before(RetentionTimeList rt, double limit)
  {
    if (rt==null) return false;

    Iterator<RetentionTime> itr = rt.listIterator();
    while (itr.hasNext())
    {
      RetentionTime r = itr.next();
      // quit if any lower bound equals or exceeds the limit
      if (r.getMinRetentionTime()/60d>=limit) return false;
    }
    return true;
  }
  public static boolean after(RetentionTimeList rt, double limit)
  {
    if (rt==null) return false;

    Iterator<RetentionTime> itr = rt.listIterator();
    while (itr.hasNext())
    {
      RetentionTime r = itr.next();
      // quit if any upper bound equals or less than the limit
      if (r.getMaxRetentionTime()/60d<=limit) return false;
    }
    return true;
  }

  public static boolean contains(RetentionTimeList rt, Range<Double> bound)
  {
    if (rt==null || bound==null) return false;

    Iterator<RetentionTime> itr = rt.listIterator();
    while (itr.hasNext())
    {
      RetentionTime r = itr.next();
      if (bound.contains(r.getMaxRetentionTime()/60d) ||
          bound.contains(r.getMinRetentionTime()/60d)) return true;
    }
    return false;
  }
  public static <S extends PeakList<PeakAnnotation>> Multimap<Integer, S> toChargePeakList(Collection<S> spectra)
  {
    if (!Tools.isSet(spectra)) return null;

    Multimap<Integer, S> z_spec = HashMultimap.create();
    for (S spec : spectra) z_spec.put(spec.getPrecursor().getCharge(), spec);

    return z_spec;
  }
  public static PeakList invalidate(PeakList A)
  {
    if (A!=null && A.size()>0)
      for (int i=0; i<A.size(); i++)
        A.setIntensityAt(Math.abs(A.getIntensity(i))*-1d, i);

    return A;
  }
  public static PeakList validate(PeakList A)
  {
    if (A!=null && A.size()>0)
      for (int i=0; i<A.size(); i++)
        A.setIntensityAt(Math.abs(A.getIntensity(i)), i);

    return A;
  }
  public static boolean isValid(PeakList A, int i)
  {
    return (A!=null && A.size()>0 && A.size()>i) ? A.getIntensity(i)>0 : false;
  }
  public static PeakList invalidate(PeakList A, int i)
  {
    if (A!=null && A.size()>0 && A.size()>i)
      A.setIntensityAt(Math.abs(A.getIntensity(i))*-1d, i);

    return A;
  }
  public static PeakList validate(PeakList A, int i)
  {
    if (A!=null && A.size()>0 && A.size()>i)
      A.setIntensityAt(Math.abs(A.getIntensity(i)), i);

    return A;
  }
  public static List<Double> intensities(PeakList spec, List<Integer> index)
  {
    if (spec==null || index==null || spec.size()<index.size()) return null;

    List<Double> ais = new ArrayList<>(index.size());
    for (int i=0; i<index.size(); i++)
      ais.add(spec.getIntensity(i));

    return ais;
  }
  public static PeakList denoise_local(PeakList A, double span, Integer cream_of_top, double sigma, boolean set_bias)
  {
    if (A==null || span == 0 || cream_of_top == null) return A;

    // the range of the locals in Z
    int          size_a = A.size();
    Range<Double> slice = null;
    Double         base = 0d; // set the lowest base at 0 to validate sparse peaks at the front
    double         xmin = A.getMz(0), xmax = A.getMz(A.size() - 1);
    Double    min_noise = Double.MAX_VALUE, min_med = Double.MAX_VALUE;

    // make sure the span is smaller than the total mass range
    if (xmax-xmin<=span) return A;

    invalidate(A);

    // make a queue for the locals
    LinkedList<Integer> locals      = new LinkedList<>();
    LinkedList<Integer> tops        = new LinkedList<>();
    Collection<Integer> locals_prev = new HashSet<>();

    // going through the points
    for (int i = 0; i < size_a; i++)
    {
      // setup the local window in X
      slice = Tools.window(slice, A.getMz(i), span, xmin, xmax);
      // strip the front if necessary
      while (locals.size() > 0 &&
          A.getMz(locals.get(0)) < slice.lowerEndpoint())
      {
        if (Tools.isSet(tops)) tops.remove(locals.getFirst());
        locals.removeFirst();
      }
      // refill the locals
      for (int j = i; j < size_a; j++)
      {
        if (Tools.isSet(locals) && A.getMz(j)<=A.getMz(locals.getLast())) continue;
        if (A.getMz(j) > slice.upperEndpoint()) break;
        locals.add(j);

        // populate the "tops" list, which contains the same neighboring peaks in the order of their intensities
        boolean added = false;
        if (!Tools.isSet(tops)) { tops.add(j); added = true; }
        else
        {
          for (int k = 0; k < tops.size(); k++)
            if (A.getIntensity(j) < A.getIntensity(tops.get(k)))
            {
              tops.add(k, j); added = true;
              break;
            }
        }
        if (!added) tops.add(j);
      }

      // tag any points in top-xx as valid
      if (locals.size() > cream_of_top)
      {
        double[] bound = Stats.outliers_rejected(intensities(A, tops), sigma, 2);
        base = bound[0]+bound[1];

        if (base     < min_noise) min_noise = base;
        if (bound[0] < min_med)   min_med   = bound[0];

        // clear out the prior peaks that couldn't made the cut
        if (Tools.isSet(locals_prev))
        {
          for (Integer t : locals_prev) if (A.getIntensity(t)>base) validate(A, t);
          locals_prev.clear();
        }
      }
      else
      {
        //for (XYPoint t : locals) t.validate();
      }
      if (base!=null && Tools.isSet(locals))
      {
        for (Integer t : locals) if (A.getIntensity(t)>base) validate(A, t);
      }
      else
      {
        locals_prev.addAll(locals);
      }
    }
    // any leftover should be kept around by default
    if (Tools.isSet(locals_prev))
      for (Integer t : locals_prev) validate(A, t);

    if (set_bias)
    {
      // setup the baseline
      List<Peak> baseline = new ArrayList<>(size_a);
      Integer   ion      = null;
      Range<Peak> bound    = null;
      double   min_base = Double.MAX_VALUE;

      for (int i = 0; i < size_a; i++)
      {
        if (isValid(A,i)) continue;
        if (A.getIntensity(i)<min_base) min_base = A.getIntensity(i);
        baseline.add(new Peak(A.getMz(i), A.getIntensity(i)));
      }
      Peaks.smoothBySG5(baseline);

      // use the min_noise instead
      min_base = min_med;

      if (min_base < Double.MAX_VALUE)
      {
        int starter = 0;
        for (int i = 0; i < size_a; i++)
        {
          try                 { ion = i; }
          catch (Exception e) { ion = null; }

          if (ion != null)
          {
            bound=null;
//            bound.setLower(null);
//            bound.setUpper(null);
            for (int k = starter; k < baseline.size() - 1; k++)
            {
              if (A.getMz(ion)>=baseline.get(k).getMz() && A.getMz(ion)<baseline.get(k+1).getMz())
              {
                starter = k;
                bound = Range.closed((Peak )baseline.get(k),(Peak )baseline.get(k+1));
//                bound.setLower((T )baseline.get(k));
//                bound.setUpper((T )baseline.get(k+1));
                break;
              }
            }
            if (bound==null) bound = Range.closed(new Peak(A.getMz(0), min_base),new Peak(A.getMz(size_a-1), min_base));
            base = Peaks.interpolateForY(bound, A.getMz(i));
            base = ((base != null && base > min_base) ? base : min_base);
            if (base != 0) A.setIntensityAt(A.getIntensity(ion)/base, ion);
          }
        }
      }
    }

    return A;
  }
  public static void notch(PeakList A, Range<Double> bound)
  {
    if (A!=null && Tools.isSet(bound))
      for (int i=0; i<A.size(); i++)
        if (bound.contains(A.getMz(i))) invalidate(A,i);
  }
  public static void unnotch(PeakList A, Range<Double> bound)
  {
    if (A!=null && Tools.isSet(bound))
      for (int i=0; i<A.size(); i++)
        if (!bound.contains(A.getMz(i))) validate(A,i);
  }
  public static void notchUpto(PeakList A, double x) { notch(A, Range.closed(0d, x)); }

  /** Find an abundance threshold given a max peak count */
  public static double threshold(PeakList data,
                                 double base, int max_count, double start_mz,
                                 double first_cut, double ratio, int repeat)
  {
    int     count = 0;
    double  cut   = first_cut * base / 100;

    for (int i = 0; i < repeat; i++)
    {
      count = 0;
      for (int p=0; p<data.size(); p++)
        if ((isValid(data,p) && data.getIntensity(p)>cut) && (data.getMz(p)>start_mz)) ++count;
      if (count >= max_count) break;
      cut /= ratio;
    }
    if (count >= max_count) return cut; else return -1;
  }

  // Adjust the base peak to avoid so-called "towering peak" issue.
  public static boolean capBasePeak(PeakList data, double rel_cut, int n_peaks)
  {
    if (data==null) return false;

    // find the base peak first
    int base = data.getMostIntenseIndex();

    if (base == -1) return false;

    double old_cut = 0.01 * rel_cut * data.getBasePeakIntensity(),
           result  = threshold(data, data.getBasePeakIntensity(), n_peaks, 0, rel_cut, 1.4, 7);

    if ((result != -1) && (result < old_cut - 0.1))
    {
      old_cut = result * data.getBasePeakIntensity() / old_cut;
      //for (unsigned int i = 0; i < getPeaks().size(); i++)
      for (int p=0; p<data.size(); p++)
        if (data.getIntensity(p)>old_cut) data.setIntensityAt(old_cut, p);
    }

    return result==-1 ? false : true;
  }
  public static int validPeaks(PeakList data)
  {
    if (data==null || data.size()==0) return 0;

    int counts=0;
    for (int i=0; i<data.size(); i++)
      if (isValid(data, i)) counts++;

    return counts;
  }
}
