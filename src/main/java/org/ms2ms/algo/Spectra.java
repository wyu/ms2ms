package org.ms2ms.algo;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.*;
import org.expasy.mzjava.core.ms.spectrum.LibPeakAnnotation;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.core.ms.spectrum.RetentionTime;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeList;
import org.ms2ms.data.ms.IsoEnvelope;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.data.ms.PeakMatch;
import org.ms2ms.io.MsReaders;
import org.ms2ms.math.Histogram;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;

import java.io.FileWriter;
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
  public enum PrepProcessor { LOCAL_NORM, LOCAL_SNR, NONE }

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
      ais.add(Math.abs(spec.getIntensity(i)));

    return ais;
  }
  public static PeakList denoise_local(PeakList A, double span, Integer cream_of_top, double sigma, boolean set_bias, boolean set_density)
  {
    if (A==null || A.size()==0 || span == 0 || cream_of_top == null) return A;

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
    Map<Integer, Double> density    = new TreeMap<>();

    // going through the points
    for (int i = 0; i < size_a; i++)
    {
      // setup the local window in X
      slice = Tools.window(A.getMz(i), span, xmin, xmax);
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
            if (Math.abs(A.getIntensity(j)) < Math.abs(A.getIntensity(tops.get(k))))
            {
              tops.add(k, j); added = true;
              break;
            }
        }
        if (!added) tops.add(j);
      }

      density.put(i, (double )locals.size()/(slice.upperEndpoint()-slice.lowerEndpoint()));

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
          for (Integer t : locals_prev) if (Math.abs(A.getIntensity(t))>base) validate(A, t);
          locals_prev.clear();
        }
      }
      else
      {
        // make everything in this local window 'good'
        for (Integer t : locals) validate(A,t);
      }
      if (base!=null && Tools.isSet(locals))
      {
        for (Integer t : locals) if (Math.abs(A.getIntensity(t))>base) validate(A, t);
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
        if (Math.abs(A.getIntensity(i))<min_base) min_base = A.getIntensity(i);
        baseline.add(new Peak(A.getMz(i), Math.abs(A.getIntensity(i))));
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
            for (int k = starter; k < baseline.size() - 1; k++)
            {
              if (A.getMz(ion)>=baseline.get(k).getMz() && A.getMz(ion)<baseline.get(k+1).getMz())
              {
                starter = k;
                bound = Range.closed(baseline.get(k),baseline.get(k+1));
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
    if (set_density)
      for (int i = 0; i < size_a; i++)
        A.setIntensityAt(A.getIntensity(i)/density.get(i), i);

    locals      = (LinkedList )Tools.dispose(locals);
    tops        = (LinkedList )Tools.dispose(tops);
    locals_prev = Tools.dispose(locals_prev);
    density     = Tools.dispose(density);

    return A;
  }
  public static void notch(PeakList A, Range<Double> bound)
  {
    if (A!=null && Tools.isSet(bound))
      for (int i=0; i<A.size(); i++)
        if (bound.contains(A.getMz(i))) invalidate(A,i);
  }
  public static double sum(PeakList A, Range<Double> bound)
  {
    double sum=0d;
    if (A!=null && Tools.isSet(bound))
      for (int i=0; i<A.size(); i++)
        if (bound.contains(A.getMz(i))) sum+=A.getIntensity(i);

    return sum;
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
  // rel_cut is the hightest rel% we'd tolerate
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

    PeakList m = Peaks.consolidate(ms, tol, 0);
    ms.clear(); ms.addPeaks(m);

    return ms;
  }
  public static PeakList accumulate(Tolerance tol, PeakList... As)
  {
    if (!Tools.isSet(As)) return null;
    if (As.length>1 && As[0]==null && As[1]!=null)
      return As[1].copy(new PurgingPeakProcessor());

    PeakList out = As[0].copy(new PurgingPeakProcessor());

    for (int i=1; i<As.length; i++) out.addPeaks(As[i]);

    PeakList merged = As[0].copy(new PurgingPeakProcessor());
    merged.clear(); merged.addPeaks(Peaks.consolidate(out, tol, 0));

    return merged;
  }
  public static double deviation(Tolerance tol, SortedMap<Double, PeakMatch> As, PeakList B)
  {
    Collection<Double> ppms = new ArrayList<>();
    for (int i=0; i<B.size(); i++)
    {
      SortedMap<Double, PeakMatch> slice = As.subMap(tol.getMin(B.getMz(i)), tol.getMax(B.getMz(i)));
      PeakMatch best=null;
      if (Tools.isSet(slice))
        for (PeakMatch pm : slice.values())
          if (best==null || pm.getIntensity()>best.getIntensity()) best=pm;

      if (best!=null) ppms.add(1E6*(best.getMz()-B.getMz(i))/B.getMz(i));
    }
    return Tools.isSet(ppms)?Stats.mean(ppms):0d;
  }
  public static MsnSpectrum accumulate(MsnSpectrum lead, OffsetPpmTolerance tol, float fr, Collection<MsnSpectrum> As)
  {
    if (!Tools.isSet(As)) return lead;

    SortedMap<Double, PeakMatch> pivots = PeakMatch.toPeaksWithExclusion(lead);
    // initiate the peaks
    for (PeakMatch pm : pivots.values()) pm.setSNR(pm.getMz()*pm.getIntensity());

//    MsnSpectrum out = lead.copy(new PurgingPeakProcessor());
//
    PpmTolerance tol2 = tol.scale(2); Collection<Peak> pcs = new ArrayList<>(); pcs.add(lead.getPrecursor());
    for (MsnSpectrum A : As)
    {
      if (A==lead) continue;

      pcs.add(A.getPrecursor());

      double delta = deviation(tol2, pivots, A), mz=0, err=0;
      // add to the pivots
      for (int i=0; i<A.size(); i++)
      {
        mz = A.getMz(i)+delta*A.getMz(i)*1E-6;
        SortedMap<Double, PeakMatch> slice = pivots.subMap(tol.getMin(mz), tol.getMax(mz));
        PeakMatch best=null;
        if (Tools.isSet(slice))
          for (PeakMatch pm : slice.values())
            if (best==null || Math.abs(pm.getMz()-mz)<err) { best=pm; err=Math.abs(pm.getMz()-mz); }

        if (best!=null)
        {
          PeakMatch pm = pivots.get(best.getMz());
          pm.setSNR(pm.getSNR()+A.getMz(i)*A.getIntensity(i));
          pm.setIntensity(pm.getIntensity()+A.getIntensity(i));
          pm.setCounts(pm.getCounts()+1);
        }
        else
        {
          pivots.put(A.getMz(i), new PeakMatch(A.getMz(i), A.getIntensity(i)).setSNR(A.getMz(i)*A.getIntensity(i)));
        }
      }
    }
    // deposit the merged peaks
    lead.clear(); int n=As.size()>1?Math.max(2, Math.round(As.size()*fr)):1;
    for (Double mz : pivots.keySet())
    {
      PeakMatch pm = pivots.get(mz);
      if (pm.getCounts()>=n)
      {
        lead.add(pm.getSNR()/pm.getIntensity(), pm.getIntensity(), new LibPeakAnnotation((int) pm.getCounts(), 0d, 0d));
      }
    }
//    if (out!=null && lead!=null && lead.getPrecursor()!=null && Tools.isSet(pcs))
    try
    {
      lead.setPrecursor(new Peak(Peaks.centroid(pcs), Peaks.AbsIntensitySum(pcs), lead.getPrecursor().getCharge()));
    }
    catch (Exception e)
    {
      System.out.println();
    }

    return lead;
  }
//  public static boolean deisotope(PeakList peaks, Tolerance tol, int maxcharge, double zzstart)
//  {
//    PeakList out = peaks.copy(new PurgingPeakProcessor()); out.clear();
//
//    Set<Integer> searched = new HashSet<>(peaks.size()); double isos=0, isobad=0;
//    for (int i=0; i<peaks.size(); i++)
//    {
//      if (searched.contains(i)) continue;
//
//      int charge=1; IsoEnvelope best=null;
//      if (peaks.getMz(i)>zzstart)
//        for (int z=2; z<=maxcharge; z++)
//        {
//          // predicted isotope envelop above 10% ri.
//          IsoEnvelope ev = Isotopes.calcIsotopesByMz(peaks.getMz(i), z, 10d, peaks.getIntensity(i));
//          // for the charge state to be real, we must detect the same numbers of the isotopes as the charge state
//          boolean ok=true;
//          for (int k=1; k<z; k++)
//          {
//            if (find(peaks, ev.getPredicted(k).getMz(), tol, i, ev.getPredicted(k).getIntensity(), 50d)<=0) { ok=false; break; }
//          }
//          // the charge envelop if good. Let's copy the peaks
//          if (ok)
//          {
//            charge=z; best=ev; break;
//          }
//        }
//
//      // transfer the peak(s)
//      if (best==null) best = new IsoEnvelope(peaks.getMz(i), 1, 10d, peaks.getIntensity(i));
//      int order=0;
//      for (Peak iso : best.getPredicted())
//      {
//        int found=find(peaks, iso.getMz(), tol, i, iso.getIntensity(), 25d);
//        if (found>=i) best.addIsotope(new Peak(peaks.getMz(found), peaks.getIntensity(found), charge));
//        if (found>=i && !searched.contains(found))
//        {
//          out.add(peaks.getMz(found)*charge - (charge-1)*1.007825d, peaks.getIntensity(found),
//                  new IsotopePeakAnnotation(iso.getCharge(), order++, peaks.getIntensity(found)));
//          searched.add(found);
//        }
//        else break;
//      }
//      // check the accuracy of the isotope distribution
//      if (order>1 && best.getPredicted(0).getMz()>peaks.getPrecursor().getMz())
//      {
//        isos++;
//        double dp = Similarity.dp(best.getPredicted().subList(0,order), best.getIsotopes().subList(0,order));
////        System.out.print(Tools.d2s(dp, 3)+",");
//        Isotopes.dp_prediction.add(Math.log10(dp)*10);
//        if (dp<0.95) isobad++;
//      }
//    }
//
////    System.out.println();
//    peaks.clear(); peaks.addPeaks(Peaks.consolidate(out, tol, 0));
//    // indicating rejection due to skewed isotope envelop
//    return (isobad/isos>=0.5);
//  }

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

      double top = (right>left && right<=max ?Collections.max(ai.subList(left, right)):0);
      peaks.setIntensityAt(100d*Math.pow(peaks.getIntensity(i), power)/(top>0?top:1d), i);
    }
    Tools.dispose(ai);

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
  public static Double self_calibrate_offset(PeakList peaks, Tolerance tol, double... stds)
  {
    if (peaks==null || peaks.size()<3 || !Tools.isSet(stds) || tol==null) return null;


    Collection<Double> offsets = new ArrayList<>();
    for (double std : stds)
    {
      int pos = find(peaks, std, tol, 0, 0d, 0d);
      if (pos>=0) offsets.add(1E6*(std-peaks.getMz(pos))/std);
    }
    return Tools.isSet(offsets)?1E-6*Stats.mean(offsets):null;
  }

  // locate the std peak and re-calibrate the rest of the peaks
  public static PeakList calibrate(PeakList peaks, Double offset)
  {
    if (offset==null || offset==0) return peaks;

    PeakList out = peaks.copy(new PurgingPeakProcessor()); out.clear();
    for (int i=0; i<peaks.size(); i++)
      out.add(peaks.getMz(i)+peaks.getMz(i)*offset, peaks.getIntensity(i), peaks.getAnnotations(i));

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
  public static double intensity(SortedMap<Double, Peak> peaks, double mz, Tolerance tol)
  {
    double abund=0d;
    SortedMap<Double, Peak> slice = peaks.subMap(tol.getMin(mz), tol.getMax(mz));
    if (Tools.isSet(slice))
      for (Peak p : slice.values())
        abund+=p.getIntensity();

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
          double m = ms.getMz(i)+ms.getMz(j) - 1.007825d;
          if (isolation.contains(m)) precursor.add(m);
        }

    precursor.generate(50);

    return 0;
  }
//  // to reproduce the QCorr algorithm from Yates's lab
//  public static int QCorr(PeakList ms, Tolerance tol, Range<Double> isolation, double step)
//  {
//    List<Peak> forward = toListOfPeaks(ms), reversed = Lists.reverse(forward);
//    for (double p = isolation.lowerEndpoint(); p<=isolation.upperEndpoint(); p+=step)
//    {
//
//    }
//    Histogram     precursor = new Histogram("");
//    double               mh = Peaks.toMH(ms.getPrecursor());
//    Range<Double> isolation = Range.closed(mh-1.5, mh+2.5);
//
////    List<Double> masses = new ArrayList<>();
//    for (int i=0; i<ms.size(); i++)
//      for (int j=0; j< ms.size(); j++)
//        if (i!=j)
//        {
//          double m = ms.getMz(i)+ms.getMz(j) + 1.007825d;
//          if (isolation.contains(m)) precursor.add(m);
//        }
//
//    precursor.generate(50);
//
//    return 0;
//  }
  public static SortedMap<Double, Peak> toPeaks(PeakList ms)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double,Peak> peaks = new TreeMap<>();

    for (int i=0; i<ms.size(); i++)
      peaks.put(ms.getMz(i), new Peak(ms.getMz(i), ms.getIntensity(i)));

    return peaks;
  }
  public static SortedMap<Double,Double> toIonMap(PeakList ms)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double,Double> peaks = new TreeMap<>();

    for (int i=0; i<ms.size(); i++)
      peaks.put(ms.getMz(i), ms.getIntensity(i));

    return peaks;
  }
  public static List<Peak> toListOfPeaks(PeakList ms) { return toListOfPeaks(ms,null); }
  public static List<Peak> toListOfPeaks(PeakList ms, Float min_mz)
  {
    // walking thro the peaks and recording the matching peaks
    List<Peak> peaks = new ArrayList<>();

    for (int i=0; i<ms.size(); i++)
      if (min_mz==null||ms.getMz(i)>=min_mz) peaks.add(new Peak(ms.getMz(i), ms.getIntensity(i)));

    return peaks;
  }
  public static List<Peak> toLMPeaks(PeakList ms, Float min_mz)
  {
    // walking thro the peaks and recording the matching peaks
    List<Peak> peaks = new ArrayList<>();

    for (int i=0; i<ms.size(); i++)
      if (min_mz==null||ms.getMz(i)<=min_mz) peaks.add(new Peak(ms.getMz(i), ms.getIntensity(i)));

    return peaks;
  }

  public static SortedMap<Double, AnnotatedPeak> toPeaksWithExclusion(PeakList ms, Range<Double>... exclusion)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double, AnnotatedPeak> peaks = new TreeMap<>();

    int i=0, left,right;
    while (i<ms.size())
    {
      // save the c13 with negative intensity, WYU 20170318
      if (Peaks.hasC13(ms.getAnnotations(i)))
      {
        peaks.put(ms.getMz(i), new AnnotatedPeak(ms.getMz(i), ((IsotopePeakAnnotation )Tools.front(ms.getAnnotations(i))).getIntensity()*-1d));
        i++; continue;
      }

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
      if      (left <0)           { left =0;          right=Math.min(left+20,ms.size()-1); }
      else if (right>ms.size()-1) { right=ms.size()-1; left=Math.max(right-20, 0); }

      AnnotatedPeak pk = new AnnotatedPeak(ms.getMz(i), ms.getIntensity(i));
      pk.setFrequency(Peaks.countC12(ms, left, right)/(ms.getMz(right)-ms.getMz(left)));
      peaks.put(mz, pk);
      // advance the pointer
      i++;
    }
    return peaks;
  }
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
  public static MsnSpectrum purgeC13(MsnSpectrum ms)
  {
    if (ms!=null && ms.size()>0)
    {
      for (int i=0; i<ms.size(); i++)
        if (Peaks.hasC13(ms.getAnnotations(i)))
        {
          ((IsotopePeakAnnotation )ms.getFirstAnnotation(i).get()).setIntensity(ms.getIntensity(i));
          ms.setIntensityAt(-1, i);
        }

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
  public static Peak base2next(MsnSpectrum ms, int bunch)
  {
    if (ms!=null && ms.size()>3)
    {
      List<Peak> peaks = new ArrayList<>();
      for (int i=0; i<ms.size(); i++)
        peaks.add(new Peak(ms.getMz(i),ms.getIntensity(i)));

      Collections.sort(peaks, new Peaks.IntensityDesendComparator());

      Collection<Double> next = new ArrayList<>();
      for (int i=1; i<peaks.size(); i++)
      {
        // look at the peaks that are outside of the base or precursor
        if (Math.abs(peaks.get(0).getMz()-peaks.get(i).getMz())>4 &&
            Math.abs(ms.getPrecursor().getMz()-peaks.get(i).getMz())>4) next.add(peaks.get(i).getIntensity());
        if (next.size()>bunch) break;
      }

      return new Peak(peaks.get(0).getMz(), peaks.get(0).getIntensity()/Stats.mean(next));
    }
    return null;
  }
  public static MsnSpectrum combineCharges(Map<String, MsnSpectrum> spectra, OffsetPpmTolerance precursor, Collection<Integer> scans)
  {
    return combineCharges(spectra, precursor, Strs.toStringArray(scans));
  }
  public static MsnSpectrum combineCharges(Map<String, MsnSpectrum> spectra, OffsetPpmTolerance precursor, String... scans)
  {
    PeakList[] specs = new PeakList[scans.length]; String line=null;
    for (int i=0; i<scans.length; i++)
    {
      specs[i] = spectra.get(scans[i]);
      line = Strs.extend(line, spectra.get(scans[i]).getScanNumbers().getFirst().getValue()+"@"+
          (specs[i].getPrecursor().getCharge()>0?"+":"")+specs[i].getPrecursor().getCharge(),";");
    }

    MsnSpectrum combined = (MsnSpectrum )Spectra.accumulate(precursor, specs);
    combined.setFragMethod(line);

    return combined;
  }
  public static Map<String, MsnSpectrum> combineCharges(Map<String, MsnSpectrum> spectra, OffsetPpmTolerance precursor, double rt_sec, boolean verbose)
  {
    // looking for charge combo
    int z3=0;
    TreeMap<               Double, MsnSpectrum>    ai_ms = new TreeMap<>(Ordering.natural().reverse());
    TreeBasedTable<Double, Double, MsnSpectrum> mh_rt_ms = TreeBasedTable.create();
    for (MsnSpectrum ms : spectra.values())
    {
      if (ms.getMsLevel()<2) ms.setMsLevel(2);
      ai_ms.put(ms.getTotalIonCurrent(), ms);
      mh_rt_ms.put(Peaks.toMH(ms.getPrecursor().getMz(), ms.getPrecursor().getCharge()), ms.getRetentionTimes().getFirst().getTime(), ms);
    }

    // go down the intensity ladder
    int cluster_counts=0;
    TreeMultimap<Double, Integer> clusters = TreeMultimap.create();
    for (Double ai : ai_ms.keySet())
    {
      MsnSpectrum M = ai_ms.get(ai);

      if (M==null || M.getMsLevel()==0) continue;

      // fetch the spectra in the vicinity
      Double        rt0 = M.getRetentionTimes().getFirst().getTime(), mh0 = Peaks.toMH(M.getPrecursor());
      Range<Double> mzr = precursor.getBoundary(mh0);
      SortedMap<Double, Map<Double, MsnSpectrum>> slice = mh_rt_ms.rowMap().subMap(mzr.lowerEndpoint(), mzr.upperEndpoint());
      if (Tools.isSet(slice))
      {
        for (Map<Double, MsnSpectrum> mm : slice.values())
          for (MsnSpectrum m : mm.values())
            if (Math.abs(m.getRetentionTimes().getFirst().getTime()-rt0)<=rt_sec) // secs
            {
              clusters.put(mh0, m.getScanNumbers().getFirst().getValue());
              m.setMsLevel(0);
            }

        if (Tools.isSet(clusters.get(mh0)) && clusters.get(mh0).size()>1) cluster_counts++;
      }
    }
    // preduce the output
    Map<String, MsnSpectrum> combo = new TreeMap<>();

    // check to see if there are any 3+ or higher MS/MS without siblings
    if (verbose) System.out.println("M+H\tScan\tz\tm/z");
    for (Double mh : clusters.keySet())
    {
      MsnSpectrum combined=combineCharges(spectra, precursor, clusters.get(mh));

      // deposit the combined for each scan
      for (Integer scan : clusters.get(mh))
        combo.put(scan+"", combined);

      if (verbose && clusters.get(mh).size()==1)
      {
        MsnSpectrum ms=spectra.get(Tools.front(clusters.get(mh))+"");
        if (ms!=null&&ms.getPrecursor().getCharge()>2)
          System.out.println(Tools.d2s(mh, 4)+"\t"+ms.getScanNumbers().getFirst().getValue()+"\t"+ms.getPrecursor().getCharge()+"\t"+Tools.d2s(ms.getPrecursor().getMz(), 4));
      }
    }
    if (verbose) System.out.println("Spectra counts (multi/cluster/total): " + cluster_counts+"/"+clusters.keySet().size()+"/"+spectra.values().size());

    return combo;
  }

  public static TreeMultimap<Double, Integer> locChargeClusters(MzMLObjectIterator<Spectrum> spectrumIterator, OffsetPpmTolerance tol, double rt_sec) throws IOException
  {
    // looping through the scans
    int counts=0;

    TreeMap<               Double, MsnSpectrum>    ai_ms = new TreeMap<>(Ordering.natural().reverse());
    TreeBasedTable<Double, Double, MsnSpectrum> mh_rt_ms = TreeBasedTable.create();
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next(), false);
      if (ms.getMsLevel()==2)
      {
        if (++counts % 1000 == 0) System.out.print(".");
        // tag the MS level to indicate unaffected spectrum
        if (ms.getMsLevel()<2) ms.setMsLevel(2);
        ai_ms.put(ms.getTotalIonCurrent(), ms);
        mh_rt_ms.put(Peaks.toMH(ms.getPrecursor().getMz(), ms.getPrecursor().getCharge()), ms.getRetentionTimes().getFirst().getTime(), ms);
      }
    }

    // go down the intensity ladder
    TreeMultimap<Double, Integer> clusters = cluster(mh_rt_ms, ai_ms, tol, rt_sec);

    ai_ms   =(TreeMap )Tools.dispose(ai_ms);
    mh_rt_ms=(TreeBasedTable )Tools.dispose(mh_rt_ms);
    return clusters;
  }
  public static TreeMultimap<Double, Integer> cluster(TreeBasedTable<Double, Double, MsnSpectrum> mh_rt_ms,
                                                      TreeMap<               Double, MsnSpectrum>    ai_ms,
                                                      OffsetPpmTolerance tol, double rt_sec)
  {
    // go down the intensity ladder
    TreeMultimap<Double, Integer> clusters = TreeMultimap.create();

    for (Double ai : ai_ms.keySet())
    {
      MsnSpectrum M = ai_ms.get(ai);

      if (M==null || M.getMsLevel()==0) continue;

      // fetch the spectra in the vicinity
      Double        rt0 = M.getRetentionTimes().getFirst().getTime(), mh0 = Peaks.toMH(M.getPrecursor());
      Range<Double> mzr = tol.getBoundary(mh0);
      SortedMap<Double, Map<Double, MsnSpectrum>> slice = mh_rt_ms.rowMap().subMap(mzr.lowerEndpoint(), mzr.upperEndpoint());
      if (Tools.isSet(slice))
      {
        for (Map<Double, MsnSpectrum> mm : slice.values())
          for (MsnSpectrum m : mm.values())
            if (Math.abs(m.getRetentionTimes().getFirst().getTime()-rt0)<=rt_sec) // secs
            {
              clusters.put(mh0, m.getScanNumbers().getFirst().getValue());
              m.setMsLevel(0);
            }
      }
    }

    return clusters;
  }
  public static boolean deisotope(PeakList peaks, Tolerance tol, int maxcharge, double zzstart, Isotopics isotopes)
  {
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
          IsoEnvelope ev = isotopes.calcIsotopesByMz(peaks.getMz(i), z, 0.05d, peaks.getIntensity(i));
          // for the charge state to be real, we must detect the same numbers of the isotopes as the charge state
          boolean ok=true;
          for (int k=1; k<z; k++)
          {
            if (Spectra.find(peaks, ev.getPredicted(k).getMz(), tol, i, ev.getPredicted(k).getIntensity(), 50d)<=0) { ok=false; break; }
          }
          // the charge envelop if good. Let's copy the peaks
          if (ok)
          {
            charge=z; best=ev; break;
          }
        }

      // transfer the peak(s)
//      if (best==null) best = new IsoEnvelope(peaks.getMz(i), 1, 10d, peaks.getIntensity(i), isotopes);
      if (best==null) best = isotopes.calcIsotopesByMz(peaks.getMz(i), 1, 0.05d, peaks.getIntensity(i)); // 20170511

      int order=0;
      for (Peak iso : best.getPredicted())
      {
        int found=Spectra.find(peaks, iso.getMz(), tol, i, iso.getIntensity(), 25d);
        if (found>=i) best.addIsotope(new Peak(peaks.getMz(found), peaks.getIntensity(found), charge));
        if (found>=i && !searched.contains(found))
        {
          out.add(peaks.getMz(found)*charge - (charge-1)*1.007825d, peaks.getIntensity(found),
              new IsotopePeakAnnotation(iso.getCharge(), order++, peaks.getIntensity(found)));
          searched.add(found);
        }
        else break;
      }
      // check the accuracy of the isotope distribution
      if (order>1 && best.getPredicted(0).getMz()>peaks.getPrecursor().getMz())
      {
        isos++;
        double dp = Similarity.dp(best.getPredicted().subList(0,order), best.getIsotopes().subList(0,order), false);
//        System.out.print(Tools.d2s(dp, 3)+",");
        Isotopes.dp_prediction.add(Math.log10(dp)*10);
        if (dp<0.95) isobad++;
      }
    }

//    System.out.println();
    peaks.clear(); peaks.addPeaks(Peaks.consolidate(out, tol, 0));
    // indicating rejection due to skewed isotope envelop
    return (isobad/isos>=0.5);
  }
  public static PeakList invalidateByPPM(PeakList A, double ppm)
  {
    if (A==null || A.size()<2) return A;

    // need at least 50ppm to avoid being labelled as "not centroided" by MSGF+
    for (int i=0; i<A.size(); i++)
      if (isValid(A, i))
        for (int j=i; j<A.size(); j++)
          if (isValid(A, j) && 1E6*Math.abs(A.getMz(j)-A.getMz(i))/A.getMz(i)>ppm)
          { invalidate(A, (A.getIntensity(i)>A.getIntensity(j))?j:i); break; }

    return A;
  }
  public static StringBuffer printProcessed(StringBuffer buf, MsnSpectrum ms, Map<Float,Float> clean, Set<Float> index)
  {
    if (buf==null) buf = new StringBuffer();

    for (int i=0; i<ms.size(); i++)
    {
      Float m = (float )ms.getMz(i);
      buf.append(i+"\t"+Tools.d2s(m,4)+"\t"+Tools.d2s(ms.getIntensity(i),2)+"\t");
      buf.append((clean.containsKey(m)?"T":"")+"\t");
      buf.append((index.contains(m)?"T":"")+"\t"+ms.getComment()+"\n");
    }

    return buf;
  }
  public static void write(FileWriter w, MsnSpectrum ms, String label, String tag) throws IOException
  {
    // scan mz  ai  z label tag
    for (int i=0; i<ms.size(); i++)
    {
      w.write(ms.getScanNumbers().getFirst().getValue()+"\t"+
          Tools.d2s(ms.getMz(i),4)+"\t"+Tools.d2s(ms.getIntensity(i),2)+"\t");
      int z=1, order=0;
      List<PeakAnnotation> annots = ms.getAnnotations(i);
      if (Tools.isSet(annots))
        for (PeakAnnotation pk : annots)
        {
          if (pk.getCharge()!=0) z=pk.getCharge();
          if (pk instanceof IsotopePeakAnnotation)
            order = ((IsotopePeakAnnotation )pk).getIsotopeOrder();
        }

      w.write(z+"\t"+(Strs.isSet(label)?label:(order>0?order+"":""))+"\t"+tag+"\n");
    }
  }
}
