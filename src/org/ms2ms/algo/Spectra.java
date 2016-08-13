package org.ms2ms.algo;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.*;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.core.ms.spectrum.RetentionTime;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeList;
import org.ms2ms.data.ms.IsoEnvelope;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.math.Histogram;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
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
        double[] bound = MsStats.outliers_rejected(intensities(A, tops), sigma, 2);
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
        if (!bound.contains(A.getMz(i))) validate(A, i);
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
  // estimate the noise floor of the spectrum
  public static double floor(PeakList A, Range<Double> bound, int poolsize, double pct)
  {
    if (A!=null && Tools.isSet(bound))
    {
      List<Double> AIs = new ArrayList<>();
      for (int i=0; i<A.size(); i++)
        if (bound.contains(A.getMz(i))) AIs.add(A.getIntensity(i));

      Collections.sort(AIs);
      int max = (int )Math.ceil(A.size()*pct*0.01);
      return Stats.mean(AIs.subList(0, poolsize>max?max:poolsize));
    }
    // signal for invalid result
    return -1;
  }
  public static Map<String, Double> survey(PeakList A, Range<Double> bound, int poolsize, double pct)
  {
    Map<String, Double> stats = new HashMap<>();
    List<Double> AIs = new ArrayList<>();

    if (A!=null && Tools.isSet(bound))
    {
      for (int i=0; i<A.size(); i++)
        if (bound.contains(A.getMz(i))) AIs.add(A.getIntensity(i));

      Collections.sort(AIs);
      int max = (int )Math.ceil(AIs.size()*pct*0.01);
      stats.put("LocalFloor",Stats.mean(AIs.subList(0, poolsize > max ? max : poolsize)));
      stats.put("LocalMin", Tools.front(AIs));
      stats.put("LocalMax", Tools.back(AIs));
    }
    if (A!=null && A.size()>1)
    {
      AIs.clear();
      for (int i=0; i<A.size(); i++) AIs.add(A.getIntensity(i));
      Collections.sort(AIs);

      int max = (int )Math.ceil(A.size()*pct*0.01);
      stats.put("GlobalFloor",Stats.mean(AIs.subList(0, poolsize > max ? max : poolsize)));
      stats.put("GlobalMin", Tools.front(AIs));
      stats.put("GlobalMax", Tools.back (AIs));
      stats.put("GlobalMedian", AIs.get((int )Math.round((AIs.size()-1)*0.5)));
      stats.put("Global95pct",  AIs.get((int )Math.round((AIs.size()-1)*0.95)));
      stats.put("Global5pct",   AIs.get((int )Math.round((AIs.size()-1)*0.05)));
    }
    // signal for invalid result
    return stats;
  }
  public static MsnSpectrum consolidate(MsnSpectrum ms, OffsetPpmTolerance tol)
  {
    if (ms==null || ms.size()<2) return ms;

    PeakList m = Peaks.consolidate(ms, tol);
    ms.clear(); ms.addPeaks(m);

    return ms;
  }
  public static PeakList deisotope(PeakList peaks, Tolerance tol, int maxcharge, double zzstart)
  {
    PeakList out = peaks.copy(new PurgingPeakProcessor()); out.clear();

    Set<Integer> searched = new HashSet<>(peaks.size());
    for (int i=0; i<peaks.size(); i++)
    {
      if (searched.contains(i)) continue;

      int charge=1; IsoEnvelope best=null;
      if (peaks.getMz(i)>zzstart)
        for (int z=2; z<=maxcharge; z++)
        {
          // predicted isotope envelop above 10% ri.
          IsoEnvelope ev = Isotopes.calcIsotopesByMz(peaks.getMz(i), z, 10d, peaks.getIntensity(i));
          // for the charge state to be real, we must detect the same numbers of the isotopes as the charge state
          boolean ok=true;
          for (int k=1; k<z; k++)
          {
            if (find(peaks, ev.getPredicted(k).getMz(), tol, i, ev.getPredicted(k).getIntensity(), 50d)<=0) { ok=false; break; }
          }
          // the charge envelop if good. Let's copy the peaks
          if (ok)
          {
            charge=z; best=ev; break;
          }
        }

      // transfer the peak(s)
      if (best==null) best=Isotopes.calcIsotopesByMz(peaks.getMz(i), 1, 10d, peaks.getIntensity(i));
      int order=0;
      for (Peak iso : best.getPredicted())
      {
        int found=find(peaks, iso.getMz(), tol, i, iso.getIntensity(), 25d);
        if (found>=i && !searched.contains(found))
        {
          out.add(peaks.getMz(found) * charge - (charge - 1) * 1.007825d, peaks.getIntensity(found),
                  new IsotopePeakAnnotation(iso.getCharge(), order++));
          searched.add(found);
        }
        else break;
      }
    }

    return Peaks.consolidate(out, tol);
  }
  public static PeakList toSNR(PeakList peaks)
  {
    double noise=Double.MAX_VALUE;
    for (int i=0; i<peaks.size(); i++)
      // without the actual reading from XCalibur, use the lowest intensity as the substitute
      if (peaks.getIntensity(i)>=0 && peaks.getIntensity(i)<noise) noise = peaks.getIntensity(i);

    for (int i=0; i<peaks.size(); i++) peaks.setIntensityAt(peaks.getIntensity(i)/noise, i);

    return peaks;
  }
  public static PeakList toLocalSNR(PeakList peaks, double window, double percentile)
  {
    SortedMap<Double, Double> background = new TreeMap<>();
    for (int i=0; i<peaks.size(); i++) background.put(peaks.getMz(i), peaks.getIntensity(i));

    double noise=Double.MAX_VALUE, min=peaks.getMz(0), max=peaks.getMz(peaks.size()-1), left=0, right=0;
    // window can be no larger than the spectrum itself
    window = Math.min(window, 0.5*(max-min));

    for (int i=0; i<peaks.size(); i++)
    {
      left  = Math.max(i-window, min); right=left+2*window;
      if (right>max) { right=max; left=right-2*window; }

      noise = Stats.percentile(background.subMap(left, right).values(), percentile);
      peaks.setIntensityAt(peaks.getIntensity(i)/noise, i);
    }

    return peaks;
  }
  public static PeakList toLocalNorm(PeakList peaks, double window)
  {
    SortedMap<Double, Double> background = new TreeMap<>();
    for (int i=0; i<peaks.size(); i++) background.put(peaks.getMz(i), peaks.getIntensity(i));

    double min=peaks.getMz(0), max=peaks.getMz(peaks.size()-1), left=0, right=0;
    // window can be no larger than the spectrum itself
    window = Math.min(window, 0.5*(max-min));

    for (int i=0; i<peaks.size(); i++)
    {
      left  = Math.max(i-window, min); right=left+2*window;
      if (right>max) { right=max; left=right-2*window; }

      peaks.setIntensityAt(100d*peaks.getIntensity(i)/Collections.max(background.subMap(left, right).values()), i);
    }

    return peaks;
  }
  public static double maxIntensity(PeakList peaks, int left, int right)
  {
    if (peaks==null || right<=left || left<0 || right>peaks.size()) return 0;

    double top=0;
    for (int i=left; i<right; i++)
      if (peaks.getIntensity(i)>top) top=peaks.getIntensity(i);

    return top;
  }
  public static PeakList toRegionalNorm(PeakList peaks, int regions)
  {
    int bundle = (int )Math.round(0.5*peaks.size()/regions), left, right, max=peaks.size()-1;

    // keep a copy of the intensities
    List<Double> ai = new ArrayList<>();
    for (int i=0; i<peaks.size(); i++) ai.add(peaks.getIntensity(i));

    for (int i=0; i<peaks.size(); i++)
    {
      left = Math.max(i-bundle, 0); right=left+2*bundle;
      if (right>max) { right=max; left=right-2*bundle; }

      double top = Collections.max(ai.subList(left, right));
      peaks.setIntensityAt(100d*peaks.getIntensity(i)/(top>0?top:1d), i);
    }

    return peaks;
  }
  // locate the std peak and re-calibrate the rest of the peaks
  public static PeakList self_calibrate(PeakList peaks, double std, Tolerance tol)
  {
    PeakList out = peaks.copy(new PurgingPeakProcessor()); out.clear();

    int pos = find(peaks, std, tol, 0, 0d, 0d);
    if (pos>=0)
    {
      double offset = (std-peaks.getMz(pos))/std;
      for (int i=0; i<peaks.size(); i++)
        out.add(peaks.getMz(i)+peaks.getMz(i)*offset, peaks.getIntensity(i), peaks.getAnnotations(i));
    }
    else return peaks;

    return out;
  }
  public static int find(PeakList peaks, double mz, Tolerance tol, int start, double targetAI, double ai_pct_tol)
  {
    for (int i=start; i<peaks.size(); i++)
      if ( tol.withinTolerance(peaks.getMz(i), mz) &&
          (targetAI<=0 || (100*(peaks.getIntensity(i)-targetAI)/targetAI)<=ai_pct_tol)) return i;
      else if (peaks.getMz(i)>tol.getMax(mz)) break;

    return -1;
  }
  public static int precursorByComplements(PeakList ms, Tolerance tol)
  {
    Histogram     precursor = new Histogram("");
    double               mh = Peaks.toMH(ms.getPrecursor());
    Range<Double> isolation = Range.closed(mh-1.5, mh+2.5);

//    List<Double> masses = new ArrayList<>();
    for (int i=0; i<ms.size(); i++)
      for (int j=0; j< ms.size(); j++)
        if (i!=j)
        {
          double m = ms.getMz(i)+ms.getMz(j) + 1.007825d;
          if (isolation.contains(m)) precursor.add(m);
        }

    precursor.generate(50);

    return 0;
  }
}
