package org.ms2ms.algo;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import org.expasy.mzjava.core.io.ms.spectrum.MgfWriter;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.*;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.core.ms.spectrum.RetentionTime;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeList;
import org.ms2ms.data.ms.IsoEnvelope;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.io.MsReaders;
import org.ms2ms.math.Histogram;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.IOException;
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
  public static boolean deisotope(PeakList peaks, Tolerance tol, int maxcharge, double zzstart)
  {
//    System.out.println();

    PeakList out = peaks.copy(new PurgingPeakProcessor()); out.clear();

    Set<Integer> searched = new HashSet<>(peaks.size()); double isos=0, isobad=0;
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
      if (best==null) best = new IsoEnvelope(peaks.getMz(i), 1, 10d, peaks.getIntensity(i));
      int order=0;
      for (Peak iso : best.getPredicted())
      {
        int found=find(peaks, iso.getMz(), tol, i, iso.getIntensity(), 25d);
        if (found>=i) best.addIsotope(new Peak(peaks.getMz(found), peaks.getIntensity(found), charge));
        if (found>=i && !searched.contains(found))
        {
          out.add(peaks.getMz(found)*charge - (charge-1)*1.007825d, peaks.getIntensity(found),
                  new IsotopePeakAnnotation(iso.getCharge(), order++));
          searched.add(found);
        }
        else break;
      }
      // check the accuracy of the isotope distribution
      if (order>1 && best.getPredicted(0).getMz()>peaks.getPrecursor().getMz())
      {
        isos++;
        double dp = Similarity.dp(best.getPredicted().subList(0,order), best.getIsotopes().subList(0,order));
//        System.out.print(Tools.d2s(dp, 3)+",");
        Isotopes.dp_prediction.add(Math.log10(dp)*10);
        if (dp<0.95) isobad++;
      }
    }

//    System.out.println();
    peaks.clear(); peaks.addPeaks(Peaks.consolidate(out, tol));
    // indicating rejection due to skewed isotope envelop
    return (isobad/isos>=0.5);
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
  public static PeakList toRegionalNorm(PeakList peaks, int regions, double power)
  {
    int bundle = (int )Math.round(0.5*peaks.size()/regions), left, right, max=peaks.size()-1;

    // keep a copy of the intensities
    List<Double> ai = new ArrayList<>();
    for (int i=0; i<peaks.size(); i++) ai.add(Math.pow(peaks.getIntensity(i), power));

    for (int i=0; i<peaks.size(); i++)
    {
      left = Math.max(i-bundle, 0); right=left+2*bundle;
      if (right>max) { right=max; left=right-2*bundle; }

      double top = Collections.max(ai.subList(left, right));
      peaks.setIntensityAt(100d*Math.pow(peaks.getIntensity(i), power)/(top>0?top:1d), i);
    }

    return peaks;
  }
  public static PeakList toRegionalRanks(PeakList peaks, int regions)
  {
    int bundle = (int )Math.round(0.5*peaks.size()/regions), left, right, max=peaks.size()-1, found=0;

    List<Double> ai = new ArrayList<>();
    for (int i=0; i<peaks.size(); i++) ai.add(peaks.getIntensity(i));

    for (int i=0; i<peaks.size(); i++)
    {
      left = Math.max(i-bundle, 0); right=left+2*bundle;
      if (right>max) { right=max; left=right-2*bundle; }

      List<Double> region = Tools.copyOf(ai, left, right+1); Collections.sort(region); found=-1;
      for (int k=0; k<region.size(); k++) if (Math.abs(region.get(k)-peaks.getIntensity(i))<0.1) { found=k; break; }
      if (found>0) peaks.setIntensityAt(100d*(1d/(double )found), i);
    }

    return peaks;
  }
  public static PeakList toRegionalPercentile(PeakList peaks, int regions)
  {
    int bundle = (int )Math.round(0.5*peaks.size()/regions), left, right, max=peaks.size()-1, found=0;

    List<Double> ai = new ArrayList<>();
    for (int i=0; i<peaks.size(); i++) ai.add(peaks.getIntensity(i));

    for (int i=0; i<peaks.size(); i++)
    {
      left = Math.max(i-bundle, 0); right=left+2*bundle;
      if (right>max) { right=max; left=right-2*bundle; }

      List<Double> region = Tools.copyOf(ai, left, right+1); Collections.sort(region); found=-1;
      for (int k=0; k<region.size(); k++) if (Math.abs(region.get(k)-peaks.getIntensity(i))<0.1) { found=k; break; }
      if (found>=0) peaks.setIntensityAt(100d*(found+1)/(double )region.size(), i);
//      else
//        System.out.println();
    }

    return peaks;
  }
  // locate the std peak and re-calibrate the rest of the peaks
  public static PeakList self_calibrate(PeakList peaks, Tolerance tol, double... stds)
  {
    if (peaks==null || peaks.size()<3 || !Tools.isSet(stds) || tol==null) return peaks;

    PeakList out = peaks.copy(new PurgingPeakProcessor()); out.clear();

    Collection<Double> offsets = new ArrayList<>();
    for (double std : stds)
    {
      int pos = find(peaks, std, tol, 0, 0d, 0d);
      if (pos>=0) offsets.add(1E6*(std-peaks.getMz(pos))/std);
    }
    if (Tools.isSet(offsets))
    {
      Double offset = 1E-6*Stats.mean(offsets);
      for (int i=0; i<peaks.size(); i++)
        out.add(peaks.getMz(i)+peaks.getMz(i)*offset, peaks.getIntensity(i), peaks.getAnnotations(i));
    }

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
  public static double intensity(PeakList peaks, double mz, Tolerance tol)
  {
    double abund=0d;
    for (int i=0; i<peaks.size(); i++)
      if (tol.withinTolerance(peaks.getMz(i), mz)) abund+=peaks.getIntensity(i);
      else if (peaks.getMz(i)>tol.getMax(mz)) break;

    return abund;
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
  public static SortedMap<Double, Peak> toPeaks(PeakList ms)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double,Peak> peaks = new TreeMap<>();

    for (int i=0; i<ms.size(); i++)
      peaks.put(ms.getMz(i), new Peak(ms.getMz(i), ms.getIntensity(i)));

    return peaks;
  }

  public static SortedMap<Double, AnnotatedPeak> toPeaksWithExclusion(PeakList ms, Range<Double>... exclusion)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double, AnnotatedPeak> peaks = new TreeMap<>();

    int i=0, left,right;
    while (i<ms.size())
    {
      // skipping the c13 isotopes
      if (Peaks.hasC13(ms.getAnnotations(i))) { i++; continue; }

      double mz=ms.getMz(i);

      // check for exclusion
      if (Tools.isSet(exclusion))
      {
        boolean found=false;
        for (Range<Double> excluded : exclusion)
          if (excluded.contains(mz)) { found=true; break; }
        if (found) { i++; continue; }
      }

      // estimates the local frequency
      left=i-10;right=i+10;
      if      (left <0)           { left =0;          right=Math.min(left+20,ms.size()); }
      else if (right>ms.size()-1) { right=ms.size()-1; left=Math.max(right-20, 0); }

      AnnotatedPeak pk = new AnnotatedPeak(ms.getMz(i), ms.getIntensity(i));
      pk.setFrequency(Peaks.countC12(ms, left, right)/(ms.getMz(right)-ms.getMz(left)));
      peaks.put(mz, pk);
      // advance the pointer
      i++;
      // need to remove the rest of the isotopes from future consideration
      while (i<ms.size())
      {
        if (!Peaks.hasC13(ms.getAnnotations(i))) break; else i++;
      }
    }
    return peaks;
  }

  // spectra quality calls if the precursor was taken from high-c13 peaks
//  public static boolean isC13Skew(PeakList ms, double ppm, double ratio, double pct)
//  {
//    // when an high-c13 peak was isolated as precursor, we will see more c13 peaks even from low-mass fragment
//    // require de-isotoping call first to ID likely c13 peaks
//    boolean isC13=false, wasC12=true;
//    for (int i=0; i<ms.size(); i++)
//    {
//      if (ms.getAnnotations())
//      isC13 = Peaks.hasC13(ms.getAnnotations(i));
//
//      if (isC13 && wasC12) goods++;
//
//      if (ms.getMz(i)>ms.getPrecursor().getMz()+2)
//      {
//        above_precursor++;
//        if (isC13 && wasC12) above_precursor_good++;
//      }
//      wasC12=(!isC13);
//
//      if (reporter_range!=null && reporter_range.contains(ms.getMz(i)))
//      {
//        reporters++;
//        if (ms.getIntensity(i)<reporter_low) reporter_low=ms.getIntensity(i);
//      }
//      if (ms.getIntensity(i)<global_low) global_low=ms.getIntensity(i);
//    }
//    return false;
//  }
  // spectra quality calls if there are >times pairs of peaks <split away while having AI-ratio<ratio
  public static boolean hasPeakSplitting(PeakList ms, double ppm, double ratio, double pct)
  {
    // an example of splitting due to MS calibration artifact
    //    230.1674194336  116523.0859375000
    //    230.1746826172  110899.7343750000
    int occurances=0, goods=0; double r = Math.abs(Math.log10(ratio)), ai0 = Peaks.getMinIntensity(ms)*1.5d;
    if (ms!=null && ms.size()>10)
    {
      for (int i=0; i<ms.size()-1; i++)
      {
        if (ms.getIntensity(i)>ai0) goods++;
        // is there a doublet?
        if (Math.abs(Peaks.toPPM(ms.getMz(i+1),ms.getMz(i)))<ppm &&
            Math.abs(Math.log10(ms.getIntensity(i+1)/ms.getIntensity(i)))<r)
          occurances++;
      }
      if (goods>0 && 200*occurances/goods>pct) return true;
//      System.out.println("   " + occurances + "/" + (int) (ms.size() * 0.5) + "/" + (Tools.d2s(200 * occurances/ms.size(), 1)));
    }
    return false;
  }
  // for a de-isotoped, charged determined peak list, return the counts for quality assessment
  public static Map<String, Object> countPeaks(PeakList ms, Range<Double> reporter_range)
  {
    Map<String, Object> stats = new HashMap<>();

    if (ms!=null && ms.size()>0)
    {
      boolean wasC12=false, isC13=false;
      int above_precursor_good=0, above_precursor=0, reporters=0, goods=0, maxz=7;
      double reporter_low=Double.MAX_VALUE, global_low=Double.MAX_VALUE, mz_low=Double.MAX_VALUE, mz_high=0;
      double[] zmz = new double[maxz]; int[] nmz = new int[maxz];

      // setup the region unique to each charge state
      for (int z=2; z<=maxz; z++) zmz[z-2] = ms.getPrecursor().getMz()*(z-1)+2;
      for (int i=0; i<ms.size(); i++)
      {
        if (ms.getMz(i)<mz_low)  mz_low =ms.getMz(i);
        if (ms.getMz(i)>mz_high) mz_high=ms.getMz(i);

        isC13 = Peaks.hasC13(ms.getAnnotations(i));

        if (isC13 && wasC12) goods++;

        if (ms.getMz(i)>ms.getPrecursor().getMz()+2)
        {
          above_precursor++;
          if (isC13 && wasC12) above_precursor_good++;
        }
        for (int z=maxz; z>1; z--)
        {
          // starting from the upper end and check for any evidence toward the charge state
          if (ms.getMz(i)>=zmz[z-2] && ms.getMz(i)<zmz[z-1]) nmz[z-2]++;
        }
        wasC12=(!isC13);

        if (reporter_range!=null && reporter_range.contains(ms.getMz(i)))
        {
          reporters++;
          if (ms.getIntensity(i)<reporter_low) reporter_low=ms.getIntensity(i);
        }
        if (ms.getIntensity(i)<global_low) global_low=ms.getIntensity(i);
      }
      // fill-in the final stats
      stats.put(Peaks.CNT_PRECURSOR_2,      above_precursor);
      stats.put(Peaks.CNT_PRECURSOR_2_GOOD, above_precursor_good);
      stats.put(Peaks.CNT_REPORTER,         reporters);
      stats.put(Peaks.CNT_GOOD,             goods);
      stats.put(Peaks.CNT_MULTIZ,           nmz);

      stats.put(Peaks.AI_MIN,          global_low);
      stats.put(Peaks.AI_MIN_REPORTER, reporter_low);

      stats.put(Peaks.CNT_GLOBAL, ms.size());
      if (mz_high>mz_low) stats.put(Peaks.CNT_PER_10AA, (int )(1150*ms.size()/(mz_high-mz_low)));
    }

    return stats;
  }
  public static MsnSpectrum prepare(MsnSpectrum ms, Tolerance precision, boolean self_calibration, double peak_transform, boolean verbose)
  {
//    PeakList deisotoped = Spectra.deisotope(ms, precision, 3, 350d);
    boolean skewed_iso = Spectra.deisotope(ms, precision, 3, 350d);

//    PeakList deisotoped = Spectra.toRegionalRanks(ms.copy(new PurgingPeakProcessor()), 7);
//    PeakList deisotoped = Spectra.toRegionalPercentile(ms.copy(new PurgingPeakProcessor()), 7);
    PeakList deisotoped = Spectra.toRegionalNorm(ms.copy(new PurgingPeakProcessor()), 7, peak_transform); // with sqrt transform
    // perform self-calibration using small fragments that are likely to show up in most of the TMT-labelled spectra
    if (self_calibration)
      deisotoped = Spectra.self_calibrate(deisotoped, new OffsetPpmTolerance(15d), 175.11949,230.170757d,376.27627);

    if (verbose)
    {
      System.out.println("#peaks: " + deisotoped.size());
      for (int i=0; i<deisotoped.size(); i++)
        System.out.println(i+"\t"+ Strs.toStrings(deisotoped.getAnnotations(i))+"\t\t"+
            Tools.d2s(deisotoped.getMz(i),8)+"\t"+Tools.d2s(deisotoped.getIntensity(i), 2));
    }

    ms.clear(); ms.addPeaks(deisotoped); if (skewed_iso) ms.setMsLevel(-1);

    return ms;
  }
  public static MsnSpectrum purgeC13(MsnSpectrum ms)
  {
    if (ms!=null && ms.size()>0)
    {
      for (int i=0; i<ms.size(); i++)
        if (Peaks.hasC13(ms.getAnnotations(i))) ms.setIntensityAt(-1, i);

      return ms.copy(new PurgingPeakProcessor<>());
    }
    return ms;
  }
  // -2.3,2.9
  public static MsnSpectrum setIsolationComment(MsnSpectrum ms, double left, double right, MsnSpectrum ms10, MsnSpectrum ms11)
  {
    // parse the isolated m/z. It maybe different from the precursor m/z
    String[] strs = ms.getComment().split("@")[0].split(" ");
    double center = (Tools.isSet(strs)? Stats.toDouble(strs[strs.length-1]):ms.getPrecursor().getMz());

    // save the precursor isolation region
    String isolation = "isolation/"+ms10.getScanNumbers().getFirst().getValue()+"/"+center+
        (ms11!=null?("/"+ms11.getScanNumbers().getFirst().getValue()):""),
        line = MsReaders.Peaks2Str(Peaks.isolate(ms10, center+left, center+right)) +
            (ms11!=null?(","+MsReaders.Peaks2Str(Peaks.isolate(ms11, center+left, center+right))):"");

    if (Strs.isSet(line)) ms.setComment(isolation+","+line);

    return ms;
  }
}
