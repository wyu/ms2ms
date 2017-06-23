package org.ms2ms.algo;

import com.google.common.collect.*;
import org.expasy.mzjava.core.mol.Mass;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.DoublePeakList;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.ms2ms.data.Point;
import org.ms2ms.data.collect.ImmutableNavigableMap;
import org.ms2ms.data.collect.NavigableMultimap;
import org.ms2ms.data.collect.TreeListMultimap;
import org.ms2ms.data.ms.FragmentEntry;
import org.ms2ms.data.ms.IsoEnvelope;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.data.ms.PeakMatch;
import org.ms2ms.math.Points;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentNavigableMap;
import java.util.concurrent.ConcurrentSkipListMap;

/** Collection of algorithms pertaining to the MS Peak
 *
 * Created by wyu on 4/23/14.
 */
public class Peaks {
  public static final String OBT_CID = "OrbitrapCID";
  public static final String LTQ_CID = "LTQCID";
  public static final String OBT_HR_CID = "OrbitrapCIDHr";
  public static final String OBT_HCD = "OrbitrapHCD";
  public static final String QTOF = "QTOF";

  public static final String CID = "cid";
  public static final String HCD = "hcd";
  public static final String ETD = "etd";

  public static final String CNT_MULTIZ = "peak counts toward each charge state";
  public static final String CNT_GLOBAL = "total peak counts";
  public static final String CNT_GOOD = "good peak counts";
  public static final String CNT_PRECURSOR_2_GOOD = "good peak counts above precursor";
  public static final String CNT_PRECURSOR_2 = "peak counts above precursor";
  public static final String CNT_REPORTER = "reporter peak counts";
  public static final String CNT_PER_10AA = "peak counts per 10 residues";
  public static final String AI_MIN_REPORTER = "min intensity in the reporter region";
  public static final String AI_MIN = "min intensity";
  public static final String PK_FACILE = "peak following facile neutral loss";
  public static final String PK_DOM = "dominance";
  public static final String MSMS = "MSMS";
  public static final String REJECTED = "Rejected";

  public static Map<String, String> sAbbr;

  static
  {
    sAbbr = new HashMap<>();

    sAbbr.put(CNT_GLOBAL,           "n.total");
    sAbbr.put(CNT_GOOD,             "n.gd");
    sAbbr.put(CNT_MULTIZ,           "n.z");
    sAbbr.put(CNT_PER_10AA,         "n.10AA");
    sAbbr.put(CNT_PRECURSOR_2,      "n>P");
    sAbbr.put(CNT_PRECURSOR_2_GOOD, "n.gd>P");
    sAbbr.put(CNT_REPORTER,         "n.rpt");
    sAbbr.put(AI_MIN,              "ai.min");
    sAbbr.put(AI_MIN_REPORTER,    "rpt.min");
    sAbbr.put(PK_DOM,             "domn");
    sAbbr.put(PK_FACILE,          "NL>");
  }
  public static String abbr(String s) { return sAbbr.get(s)!=null?sAbbr.get(s):s; }
  public static class IntensityDesendComparator implements Comparator<Peak> {
    public int compare(Peak o1, Peak o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getIntensity(), o1.getIntensity()) : 0;
    }
  }

  public static class IntensityAscendComparator implements Comparator<Peak> {
    public int compare(Peak o1, Peak o2) {
      return o1 != null && o2 != null ? Double.compare(o1.getIntensity(), o2.getIntensity()) : 0;
    }
  }

  public static class MzAscendComparator implements Comparator<Peak> {
    public int compare(Peak o1, Peak o2) {
      return o1 != null && o2 != null ? Double.compare(o1.getMz(), o2.getMz()) : 0;
    }
  }

  public static class MzDesendComparator implements Comparator<Peak> {
    public int compare(Peak o1, Peak o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getMz(), o1.getMz()) : 0;
    }
  }

  public static boolean isType(PepLibPeakAnnotation s, IonType... types) {
    try {
      IonType ion = s.getOptFragmentAnnotation().get().getIonType();
      if (s != null && Tools.isSet(types))
        for (IonType t : types)
          if (t.equals(ion)) return true;
    } catch (IllegalStateException e) {

    }

    return false;
  }

  public static boolean hasType(Collection<PepLibPeakAnnotation> s, IonType... types) {
    if (Tools.isSet(s))
      for (PepLibPeakAnnotation A : s)
        if (isType(A, types)) return true;

    return false;
  }

  /**
   * Count the number of valid peaks between x1 and x2
   *
   * @param ms  is a PeakList object. zero or negative intensity indicates invalid peak
   * @param x1, x2 define the m/z bound where the counting is to take place
   * @return the counts of peaks bound by (x1, x2) that's above 0
   */
  public static int countValid(PeakList ms, double x1, double x2) {
    if (ms == null || ms.size() == 0) return 0;

    int counts = 0;
    for (int i = 0; i < ms.size(); i++)
      if (ms.getMz(i) >= x1 && ms.getMz(i) <= x2 && isValidIntensity(ms.getIntensity(i))) counts++;

    return counts;
  }

  public static double getMinIntensity(PeakList msms) {
    return getMinIntensity(msms, Double.MIN_VALUE, Double.MAX_VALUE);
  }

  public static double getMinIntensity(PeakList msms, double x0, double x1) {
    if (msms == null || msms.size() == 0) return 0;

    double baseline = Double.MAX_VALUE;
    for (int i = 0; i < msms.size(); i++)
      if (msms.getMz(i) >= x0 && msms.getMz(i) <= x1 && msms.getIntensity(i) < baseline)
        baseline = msms.getIntensity(i);

    return baseline;
  }

  public static double getmeanIntensity(PeakList msms) {
    return getmeanIntensity(msms, Double.MIN_VALUE, Double.MAX_VALUE);
  }

  public static double getmeanIntensity(PeakList msms, double x0, double x1) {
    if (msms == null || msms.size() == 0) return 0;

    double sum = 0d;
    for (int i = 0; i < msms.size(); i++)
      if (msms.getMz(i) >= x0 && msms.getMz(i) <= x1) sum += msms.getIntensity(i);

    return sum / (double) msms.size();
  }

  public static double getMeanIntensity(Collection<Peak> msms) {
    if (msms == null || msms.size() == 0) return 0;

    double sum = 0d;
    for (Peak p : msms) sum += p.getIntensity();

    return sum / (double) msms.size();
  }

  public static double getBaseline(PeakList msms, double x0, double x1, int top, boolean exclude_precursor) {
    List<Peak> baselines = new ArrayList<Peak>();
    for (int i = 0; i < msms.size(); i++)
      if (isValidIntensity(msms.getIntensity(i)) &&
          (!exclude_precursor || !hasType(msms.getAnnotations(i), IonType.p)) &&
          msms.getMz(i) >= x0 && msms.getMz(i) <= x1)
        baselines.add(new Peak(msms.getMz(i), msms.getIntensity(i), 1));

    Collections.sort(baselines, new IntensityAscendComparator());
    if (baselines.size() > top) {
      Collection<Double> pts = new ArrayList<Double>();
      for (int i = 0; i < baselines.size() - top; i++) pts.add(baselines.get(i).getIntensity());
      return MsStats.mean(pts) + MsStats.stdev(pts);
    }
    return Tools.isSet(baselines) ? -1d * getMeanIntensity(Lists.partition(baselines, 3).get(0)) : 0;
  }

  public static boolean isValidIntensity(double s) {
    return s > 0;
  }

  public static void invalidate(Peak... s) {
    if (s != null) for (Peak a : s) a.setIntensity(Math.abs(a.getIntensity()) * -1d);
  }

  public static void validate(Peak... s) {
    if (s != null) for (Peak a : s) a.setIntensity(Math.abs(a.getIntensity()));
  }

  public static void validate(List<? extends Peak>... As) {
    if (As != null)
      for (List<? extends Peak> A : As)
        for (Peak a : A) validate(a);
  }

  public static void invalidate(List<? extends Peak>... As) {
    if (As != null)
      for (List<? extends Peak> A : As)
        for (Peak a : A) invalidate(a);
  }

  public static boolean isValid(Peak s) {
    return s.getIntensity() > 0;
  }

  public static StringBuffer printIon(StringBuffer buf, double mz, double ai, int z) {
    if (buf == null) buf = new StringBuffer();
    buf.append(Tools.d2s(mz, 4) + ", " + Tools.d2s(ai, 1) + ", +" + z);
    return buf;
  }

  public static StringBuffer print(StringBuffer buf, Peak p) {
    return printIon(buf, p.getMz(), p.getIntensity(), p.getCharge());
  }

  public static StringBuffer print(StringBuffer buf, AnnotatedPeak p) {
    buf = printIon(buf, p.getMz(), p.getIntensity(), p.getCharge());
    buf.append(", " + Tools.d2s(p.getSNR(), 2));
    return buf;
  }

  public static StringBuffer printAnnot(StringBuffer buf, Collection<PepLibPeakAnnotation> annotations) {
    if (annotations != null) {
      if (buf == null) buf = new StringBuffer();
      for (PepLibPeakAnnotation anno : annotations) {
        PepFragAnnotation f = anno.getOptFragmentAnnotation().get();
        double loss = f.getNeutralLoss().getMolecularMass();
        buf.append("+" + f.getCharge() +
            "," + f.getIonType() + (loss != 0 ? (loss > 0 ? "+" : "") + Math.round(loss) : "") + ", "
            + Tools.d2s(f.getTheoreticalMz(), 4) + ", " + f.getFragment() + ", " + f
            .getIsotopeComposition() + "; ");
      }
    }
    return buf;
  }

  public static StringBuffer print(StringBuffer buf, PeakList<PepLibPeakAnnotation> msms, boolean annot) {
    if (msms == null || msms.size() == 0) return null;

    if (buf == null) buf = new StringBuffer();
    buf.append("Precursor: ");
    buf = printIon(buf, msms.getPrecursor().getMz(), msms.getPrecursor().getIntensity(), msms.getPrecursor().getCharge());
    buf.append("\n");
    for (int i = 0; i < msms.size(); i++) {
      // skip the peak if annotation is required
      if (annot && msms.getAnnotations(i) == null) continue;

      buf = printIon(buf, msms.getMz(i), msms.getIntensity(i), 0);
      buf.append(", ");
      if (annot) buf = printAnnot(buf, msms.getAnnotations(i));
      buf.append("\n");
    }

    return buf;
  }

  public static boolean isCtermFrag(IonType s) {
    return s != null ? (s.equals(IonType.y) || s.equals(IonType.x) || s.equals(IonType.z) || s.equals(IonType.w)) : false;
  }

  public static boolean isNtermFrag(IonType s) {
    return s != null ? (s.equals(IonType.b) || s.equals(IonType.a) || s.equals(IonType.c) || s.equals(IonType.d)) : false;
  }

  public static double toMH(Peak p) {
    return toMH(p.getMz(), p.getCharge());
  }

  public static double toMH(double mz, int z) { return (mz * z - (z - 1) * 1.007825); }
  public static double MnH2MnH(double mh, int z1, int z2) { return (toMass(mh,z1)+z2*Peptides.H)/(double )z2; }
  public static Peak   MnH2MnH(Peak mh, int z1, int z2) { return new Peak(MnH2MnH(mh.getMz(), z1, z2), mh.getIntensity(), mh.getCharge()); }
  public static float  MnH2MnH(float mh, int z1, int z2) { return (float )(toMass((double )mh,z1)+z2*Peptides.H)/(float )z2; }

  public static double toMass(double mz, int z) {
    return (mz - 1.007825) * z;
  }

  public static double toMass(Peak p) {
    return toMass(p.getMz(), p.getCharge());
  }

  public static double toPPM(double m0, double m1) {
    return 1E6 * (m1 - m0) / m0;
  }

  public static double MH2Mz(double mh, int z) {
    return (mh + 1.007825 * (z - 1)) / (double) z;
  }

  /**
   * makeup a peaklist using string shorthand
   *
   * @param mz    and z are the m/z value and charge state of the precursor.
   * @param frags variable number of fragment ions. e.g. "334.5", "562/23", "mz/ai". Only the mz is required
   * @return
   */
  public static PeakList<PepLibPeakAnnotation> newPeakList(double mz, int z, String... frags)
  {
    PeakList<PepLibPeakAnnotation> spec = new DoublePeakList<>();
    spec.setPrecursor(new Peak(mz, 0d, z));
    // go thro the fragment ions if any
    for (String f : frags) {
      try {
        String[] strs = f.split("/");
        spec.add(Double.valueOf(strs[0]), strs.length > 2 ? Double.valueOf(strs[2]) : 0d);
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
    return spec;
  }

  /**
   * Parse the precursor specification from the web form
   *
   * @param line : "501.1/2 5051/3 505.6"
   * @param z    : the default charge if not specified with the m/z above
   * @return array of precursor peaks
   */
  public static Peak[] newPeaks(String line, int z) {
    if (!Strs.isSet(line)) return null;
    String[] strs = line.split(";|,|\\s");
    Peak[] peaks = new Peak[strs.length];
    for (int i = 0; i < strs.length; i++) {
      String[] items = strs[i].split("/");
      peaks[i] = new Peak(Double.valueOf(items[0]), 0d, items.length > 1 ? Integer.valueOf(items[1]) : z);
    }
    return peaks;
  }

  public static <T extends Peak> T getBasePeak(Collection<T> data) {
    if (!Tools.isSet(data)) return null;

    T base = null;
    for (T datum : data)
      if (base == null || (datum.getIntensity() > base.getIntensity())) base = datum;

    // send the base peak back
    return base;
  }

  public static int getBasePeak(double[] ais) {
    if (!Tools.isSet(ais)) return -1;

    int base = -1;
    for (int i = 0; i < ais.length; i++)
      if (base == -1 || (ais[i] > ais[base])) base = i;

    // send the base peak back
    return base;
  }

  public static <T extends Peak> T find(List<T> ms, T m, Tolerance tol) {
    // locate the point that's the cloest to 'm'
    T found = null;
    double best = Double.MAX_VALUE, fold;
    for (T t : ms)
      if (tol.withinTolerance(t.getMz(), m.getMz())) {
        fold = Math.abs(m.getIntensity() != 0d ? Math.log(t.getIntensity() / m.getIntensity()) : (t.getMz() - m.getMz()));
        if (found == null || fold < best) {
          found = t;
          best = fold;
        }
      }

    return found;
  }

  public static int find(double[] mzs, double m, Tolerance tol) {
    // locate the point that's the cloest to 'm'
    int found = -1;
    double best = Double.MAX_VALUE, delta;
    for (int t = 0; t < mzs.length; t++)
      if (tol.withinTolerance(mzs[t], m)) {
        if (found == -1 || Math.abs(mzs[t] - m) < best) {
          found = t;
          best = Math.abs(mzs[t] - m);
        }
      }

    return found;
  }

  /**
   * Remove the insignificant peaks from further consideration
   * The peak significance is local context dependent.
   *
   * @param A            = Profile
   * @param span         = half-width of the local window in dalton
   * @param cream_of_top = # of top peaks in each local window to be kept
   * @param sigma
   * @param set_bias     = whether to set the local noise level @return Numbers of the validated points from the profile
   */
  public static <T extends Peak> List<T> local_noise(List<T> A, double span, Integer cream_of_top, double sigma, boolean set_bias) {
    if (!Tools.isSet(A) || span == 0 || cream_of_top == null) return A;

    // make sure the points are in order
    Collections.sort(A, new MzAscendComparator());
    // make sure the span is smaller than the total mass range
    if (Tools.back(A).getMz() - Tools.front(A).getMz() <= span) return A;

    invalidate(A);

    // make a queue for the locals
    LinkedList<T> locals = new LinkedList<>();
    LinkedList<T> tops = new LinkedList<>();

    // the range of the locals in Z
    int size_a = A.size();
    Range<Double> slice = null;
    Double base = 0d; // set the lowest base at 0 to validate sparse peaks at the front
    double xmin = Tools.front(A).getMz(), xmax = Tools.back(A).getMz();
    Double min_noise = Double.MAX_VALUE, min_med = Double.MAX_VALUE;

    Collection<Peak> locals_prev = new HashSet<>();

    // going through the points
    for (int i = 0; i < size_a; i++) {
      // setup the local window in X
      slice = Tools.window(slice, A.get(i).getMz(), span, xmin, xmax);
      // strip the front if necessary
      while (locals.size() > 0 &&
          locals.get(0).getMz() < slice.lowerEndpoint()) {
        if (Tools.isSet(tops)) tops.remove(locals.getFirst());
        locals.removeFirst();
      }
      // refill the locals
      for (int j = i; j < size_a; j++) {
        if (Tools.isSet(locals) && A.get(j).getMz() <= locals.getLast().getMz()) continue;
        if (A.get(j).getMz() > slice.upperEndpoint()) break;
        locals.add(A.get(j));

        // populate the "tops" list, which contains the same neighboring peaks in the order of their intensities
        boolean added = false;
        if (!Tools.isSet(tops)) {
          tops.add(A.get(j));
          added = true;
        } else {
          for (int k = 0; k < tops.size(); k++)
            if (A.get(j).getIntensity() < tops.get(k).getIntensity()) {
              tops.add(k, A.get(j));
              added = true;
              break;
            }
        }
        if (!added) tops.add(A.get(j));
      }

      // tag any points in top-xx as valid
      if (locals.size() > cream_of_top) {
        Peak bound = outliers_rejected(tops, sigma, 2);
        base = bound.getMz() + bound.getIntensity();

        if (base < min_noise) min_noise = base;
        if (bound.getMz() < min_med) min_med = bound.getMz();

        // clear out the prior peaks that couldn't made the cut
        if (Tools.isSet(locals_prev)) {
          for (Peak t : locals_prev) if (t.getIntensity() > base) validate(t);
          locals_prev.clear();
        }
      } else {
        //for (XYPoint t : locals) t.validate();
      }
      if (base != null && Tools.isSet(locals)) {
        for (Peak t : locals) if (t.getIntensity() > base) validate(t);
      } else {
        locals_prev.addAll(locals);
      }
    }
    // any leftover should be kept around by default
    if (Tools.isSet(locals_prev))
      for (Peak t : locals_prev) validate(t);

    if (set_bias) {
      // setup the baseline
      List<Peak> baseline = new ArrayList<>(size_a);
      Peak ion = null;
      Range<T> bound = null;
      double min_base = Double.MAX_VALUE;

      for (int i = 0; i < size_a; i++) {
        if (isValid(A.get(i))) continue;
        if (A.get(i).getIntensity() < min_base) min_base = A.get(i).getIntensity();
        baseline.add(new Peak(A.get(i).getMz(), A.get(i).getIntensity()));
      }
      smoothBySG5(baseline);

      // use the min_noise instead
      min_base = min_med;

      if (min_base < Double.MAX_VALUE) {
        int starter = 0;
        for (int i = 0; i < size_a; i++) {
          try {
            ion = A.get(i);
          } catch (Exception e) {
            ion = null;
          }

          if (ion != null) {
            bound = null;
//            bound.setLower(null);
//            bound.setUpper(null);
            for (int k = starter; k < baseline.size() - 1; k++) {
              if (ion.getMz() >= baseline.get(k).getMz() && ion.getMz() < baseline.get(k + 1).getMz()) {
                starter = k;
                bound = Range.closed((T) baseline.get(k), (T) baseline.get(k + 1));
//                bound.setLower((T )baseline.get(k));
//                bound.setUpper((T )baseline.get(k+1));
                break;
              }
            }
            if (bound == null)
              bound = Range.closed((T) new Peak(A.get(0).getMz(), min_base), (T) new Peak(A.get(size_a - 1).getMz(), min_base));
            base = interpolateForY(bound, A.get(i).getMz());
            base = ((base != null && base > min_base) ? base : min_base);
            if (base != 0) ion.setIntensity(ion.getIntensity() / base);
          }
        }
      }
    }

    return A;
  }

  public static <T extends Peak> double filter(List<T> A, int index_begin, double[] filters) {
    double Y = 0d;
    for (int i = 0; i < filters.length; i++)
      Y += A.get(i + index_begin).getIntensity() * filters[i];

    return Y;
  }

  //--------------------------------------------------------------------------

  /**
   * Apply the pre-calculated coeffs to the data according to "http://www.ma.utexas.edu/documentation/nr/bookcpdf/c14-8.pdf"
   * <p>
   * M nL nR Sample Savitzky-Golay Coefficients
   * 2 2 2 ?0.086 0.343 0.486 0.343 ?0.086
   * 2 3 1 ?0.143 0.171 0.343 0.371 0.257
   * 2 4 0 0.086 ?0.143 ?0.086 0.257 0.886
   * 2 5 5 ?0.084 0.021 0.103 0.161 0.196 0.207 0.196 0.161 0.103 0.021 ?0.084
   * 4 4 4 0.035 ?0.128 0.070 0.315 0.417 0.315 0.070 ?0.128 0.035
   * 4 5 5 0.042 ?0.105 ?0.023 0.140 0.280 0.333 0.280 0.140 ?0.023 ?0.105 0.042
   *
   * @param A
   */
  public static <T extends Peak> List<T> smoothBySG5(List<T> A) {
    // Do nothing if the set isn't big enough to smooth.
    if (null == A || A.size() < 6) return A;

    // store the smoothed data separately
    List<T> smoothed = new ArrayList<T>(A.size());
    T cloned = null;

    // special handling for the edge points
    cloned = (T) new Peak(A.get(0));
    cloned.setIntensity(filter(A, 0, new double[]{0.886, 0.257, -0.086, -0.143, 0.086}));
    smoothed.add(cloned);

    cloned = (T) new Peak(A.get(1));
    cloned.setIntensity(filter(A, 0, new double[]{0.257, 0.371, 0.343, 0.171, -0.143}));
    smoothed.add(cloned);

    for (int i = 2; i < A.size() - 2; i++) {
      cloned = (T) new Peak(A.get(i));
      cloned.setIntensity(filter(A, i - 2, new double[]{-0.086, 0.343, 0.486, 0.343, -0.086}));
      smoothed.add(cloned);
    }

    // special handling for the edge points
    cloned = (T) new Peak(A.get(A.size() - 2));
    cloned.setIntensity(filter(A, A.size() - 6, new double[]{-0.143, 0.171, 0.343, 0.371, 0.257}));
    smoothed.add(cloned);

    cloned = (T) new Peak(A.get(A.size() - 1));
    cloned.setIntensity(filter(A, A.size() - 6, new double[]{0.086, -0.143, -0.086, 0.257, 0.886}));
    smoothed.add(cloned);

    for (int i = 0; i < A.size(); i++)
      A.get(i).setIntensity(smoothed.get(i).getIntensity());

    return A;
  }

  public static <T extends Peak> Double interpolateForY(Range<T> bound, Double x) {
    if (bound == null) return null;
    return interpolateForY(bound.lowerEndpoint(), bound.upperEndpoint(), x);
  }

  public static <T extends Peak> Double interpolateForY(T t0, T t1, Double x) {
    if (t0 == null ||
        t1 == null ||
        x == null) return null;

    if (x < t0.getMz() || x > t1.getMz()) return null;

    return (t0.getIntensity() + (t1.getIntensity() - t0.getIntensity()) * (x - t0.getMz()) / (t1.getMz() - t0.getMz()));
  }

  public static <T extends Peak> Double interpolateForX(T t0, T t1, Double y) {
    if (t0 == null ||
        t1 == null ||
        y == null || t1.getIntensity() == t0.getIntensity()) return null;

    return (t0.getMz() + (t1.getMz() - t0.getMz()) * (y - t0.getIntensity()) / (t1.getIntensity() - t0.getIntensity()));
  }

  public static Map<Peak, Peak> overlap(List<? extends Peak> A, List<? extends Peak> B, double delta,
                                        boolean matchHighest, boolean keep_unmatched, Map<Peak, Peak> outcomes) {
    boolean verbose = false;

    // false if one of the spectra is empty
    if (A == null || A.isEmpty() || B == null || B.isEmpty()) return null;

    if (outcomes == null) outcomes = new TreeMap<Peak, Peak>();
    else outcomes.clear();

    int ia = 0, ib = 0, start_b;
    double best_Y = 0, best_error = 1E23, abs_err;

    while (ia < A.size()) {
      // ignore any invalidated point
      if (!isValid(A.get(ia))) {
        ++ia;
        continue;
      }
//      if (verbose) System.out.print(">>" + A.get(ia).getX() + ", " + A.get(ia).getY());
      // reset the goalposts
      best_Y = 0;
      best_error = 1E23;
      start_b = -1;
      // recycle to avoid access violation
      if (ib >= B.size()) ib = 0;
      while (ib < B.size()) {
        // skip the minor points that's labled
        if (!isValid(B.get(ib))) {
          ++ib;
          continue;
        }
        if (B.get(ib).getMz() > A.get(ia).getMz() + delta + 0.1) break;
        if (!equivalent(A.get(ia), B.get(ib), (float) delta)) {
          ++ib;
          continue;
        }
        if (start_b == -1) start_b = ib;
        // find the highest peak in the matching window
        if (matchHighest && (B.get(ib).getIntensity() > best_Y)) {
          best_Y = B.get(ib).getIntensity();
          if (outcomes != null) outcomes.put(A.get(ia), B.get(ib));
        }
        // WYU 10-13-2000, option for matching the closest
        else {
          abs_err = Math.abs(B.get(ib).getMz() - A.get(ia).getMz());
          if (!matchHighest && (abs_err < best_error)) {
            best_Y = B.get(ib).getIntensity();
            best_error = abs_err;
            if (outcomes != null) outcomes.put(A.get(ia), B.get(ib));
          }
        }
        ++ib;
      }
      if (best_Y == 0 && keep_unmatched) {
        if (outcomes != null) outcomes.put(A.get(ia), new Peak(0d, 0d));
      }

      // increment the outer loop
      ++ia;
      // move the b iterator back a little
      if (start_b > 0) {
        ib = start_b - 1;
        start_b = -1;
      } else ib = 0;
    }

    return outcomes;
  }

  public static boolean equivalent(Peak A, Peak B, float tol) {
    return A != null && B != null && Math.abs(A.getMz() - B.getMz()) <= tol;
  }

  public static boolean isC12(Peak p) {
    return p.getCharge() >= 0;
  }

  public static <T extends Peak> Peak outliers_rejected(Collection<T> A, double stdev_multiple, int rounds) {
    int ys_i = 0;
    double[] ys = new double[A.size()];
    double avg = 0, bound = Double.MAX_VALUE;

    for (T t : A) t.setMzAndCharge(Math.abs(t.getMz()));

    for (int i = 0; i < rounds; i++) {
      ys_i = 0;
      for (T t : A) if (t.getMz() >= 0d && Math.abs(t.getIntensity() - avg) <= bound) ys[ys_i++] = t.getIntensity();

      avg = MsStats.mean(ys, ys_i);
      bound = MsStats.stdev(ys, ys_i) * stdev_multiple;
      // are you converge yet?

      // remove the outlier from the consideration
      for (T t : A)
        if (t.getIntensity() - avg > bound)
          t.setMzAndCharge(Math.abs(t.getMz()) * -1d);
    }

    for (T t : A) t.setMzAndCharge(Math.abs(t.getMz()));

    return new Peak(avg, bound);
  }

  /**
   * Tag the points around a location as invalid, in the manner of a notch filter
   * Typically used to remove precursor ion from further consideration
   *
   * @param A     is a list of peaks, not necessary sorted by m/z
   * @param bound define the lower and upper limit of the m/z
   */
  public static void notch(Collection<? extends Peak> A, Range<Double> bound) {
    if (Tools.isSet(A) && Tools.isSet(bound))
      for (Peak p : A)
        if (bound.contains(p.getMz())) invalidate(p);
  }

  public static void unnotch(Collection<? extends Peak> A, Range<Double> bound) {
    if (Tools.isSet(A) && Tools.isSet(bound))
      for (Peak p : A)
        if (!bound.contains(p.getMz())) validate(p);
  }

  public static void notchUpto(Collection<? extends Peak> A, double x) {
    notch(A, Range.closed(0d, x));
  }

  /**
   * Find an abundance threshold given a max peak count
   */
  public static double threshold(List<? extends Peak> data,
                                 double base, int max_count, double start_mz,
                                 double first_cut, double ratio, int repeat) {
    int count = 0;
    double cut = first_cut * base / 100;

    for (int i = 0; i < repeat; i++) {
      count = 0;
      for (Peak p : data)
        if ((isValid(p) && p.getIntensity() > cut) && (p.getMz() > start_mz)) ++count;
      if (count >= max_count) break;
      cut /= ratio;
    }
    if (count >= max_count) return cut;
    else return -1;
  }

  public static Range<Double> mz2window(double mz, Tolerance tol) {
    return Range.closed(tol.getMin(mz), tol.getMax(mz));
  }

  public static TreeMap<Double, Range<Double>> peaks2windows(PeakList peaks, Tolerance tol) {
    TreeMap<Double, Range<Double>> windows = new TreeMap<>();
    for (int i = 0; i < peaks.size(); i++) {
      windows.put(peaks.getMz(i), mz2window(peaks.getMz(i), tol));
    }

    return windows;
  }

  public static TreeMap<Double, Double> peaks2dic(PeakList peaks) {
    TreeMap<Double, Double> windows = new TreeMap<>();
    for (int i = 0; i < peaks.size(); i++) {
      windows.put(peaks.getMz(i), peaks.getIntensity(i));
    }

    return windows;
  }

  //** Query functions **//
  public static Collection<FragmentEntry> query(TreeMultimap<Float, FragmentEntry> indices, Range<Double> mz) {
    return Tools.isSet(mz) ? Tools.flatten(indices.asMap().subMap(mz.lowerEndpoint().floatValue(), mz.upperEndpoint().floatValue()).values()) : null;
  }

  public static boolean match(TreeMultimap<Float, FragmentEntry> indices, Range<Double> mz) {
    return Tools.isSet(mz) ? Tools.isSet(indices.asMap().subMap(mz.lowerEndpoint().floatValue(), mz.upperEndpoint().floatValue())) : false;
  }

  // TODO need to clear the list when done!
  public static List<FragmentEntry> query(NavigableMultimap<Float, FragmentEntry> indices, float[] range)
  {
    return range != null ? indices.subList(range[0], range[1]) : null;
  }
//  public static List<FragmentEntry> query(NavigableMultimap<Float, FragmentEntry> indices, float mh, OffsetPpmTolerance tol, int maxZ)
//  {
//    // setup the 1+ fragments
//    float[]               R = tol.toExpectedBoundary(mh);
//    List<FragmentEntry> out = indices.subList(R[0], R[1]);
//
//    if (maxZ>1 && out!=null)
//      for (int z=2; z<=maxZ; z++)
//      {
//        R = tol.toExpectedBoundary((float )Peaks.MnH2MnH(mh,1,z));
//        Tools.addAll(out, indices.subList(R[0], R[1]));
//      }
//
//    return out;
//  }

  public static boolean match(TreeListMultimap<Double, FragmentEntry> indices, Range<Double> mz) {
    return Tools.isSet(mz) ? indices.containsKey(mz.lowerEndpoint(), mz.upperEndpoint()) : false;
  }

  public static PeakList consolidate(PeakList peaks, Tolerance tol, int min_pks)
  {
    if (peaks == null || peaks.size() < 2) return peaks;

    Collection<Point> pts = new ArrayList<>();
    Collection<Peak> news = new ArrayList<>();
    Multimap<Integer, PeakAnnotation> pa = HashMultimap.create();
    for (int i = 0; i < peaks.size(); i++)
    {
      double max = tol.getMax(peaks.getMz(i));
      pts.clear();
      for (int j = i + 1; j < peaks.size(); j++)
      {
        if (j < peaks.size() && peaks.getMz(j) <= max)
        {
          pts.add(new Point(peaks.getMz(j), peaks.getIntensity(j)));
          peaks.setIntensityAt(-1d, j);
          if (peaks.getAnnotations(j) != null) pa.putAll(i, peaks.getAnnotations(j));
        } else break;
      }
      // require a min number of peaks
      if (pts.size() > min_pks)
      {
        pts.add(new Point(peaks.getMz(i), peaks.getIntensity(i)));
        if (peaks.getAnnotations(i) != null) pa.putAll(i, peaks.getAnnotations(i));

        peaks.setIntensityAt(-1d, i);
        news.add(new Peak(Points.centroid(pts), Points.sumY(pts), i));
      }
    }
    if (Tools.isSet(news))
      for (Peak xy : news)
        peaks.add(xy.getMz(), xy.getIntensity(), pa.get(xy.getCharge()));

    // dispose the intermediate objects
    Tools.dispose(pts, news); Tools.dispose(pa);

    // keeping just the peaks with positive intensities
    return peaks.copy(new PurgingPeakProcessor());
  }

  public static <T extends Peak> Double centroid(Collection<T> points) {
    return centroid(points, null, null);
  }

  public static <T extends Peak> Double centroid(Collection<T> points, Double x0, Double x1)
  {
    if (!Tools.isSet(points)) return null;

    double sumXY = 0, sumY = 0;
    for (Peak xy : points) {
      if ((x0 == null || xy.getMz() >= x0) &&
          (x1 == null || xy.getMz() <= x1)) {
        sumXY += xy.getMz() * xy.getIntensity();
        sumY += xy.getIntensity();
      }
    }
    return sumY != 0 ? sumXY / sumY : null;
  }
  public static <T extends Peak> Double centroid(T[] points, Double x0, Double x1) {
    if (!Tools.isSet(points)) return null;

    double sumXY = 0, sumY = 0;
    for (Peak xy : points) {
      if ((x0 == null || xy.getMz() >= x0) &&
          (x1 == null || xy.getMz() <= x1)) {
        sumXY += xy.getMz() * xy.getIntensity();
        sumY += xy.getIntensity();
      }
    }
    return sumY != 0 ? sumXY / sumY : null;
  }

  // any c13 isotope in the annotations?
  public static boolean hasC13(Collection<PeakAnnotation> pa) {
    if (Tools.isSet(pa))
      for (PeakAnnotation p : pa)
        if (p instanceof IsotopePeakAnnotation &&
            ((IsotopePeakAnnotation) p).getIsotopeOrder() > 0) return true;

    return false;
  }

  // any c13 isotope in the annotations?
  public static boolean hasCharge(Collection<PeakAnnotation> pa) {
    if (Tools.isSet(pa))
      for (PeakAnnotation p : pa)
        if (p.getCharge() > 0) return true;

    return false;
  }

  public static boolean isc12(Collection<IsotopePeakAnnotation> annotations) {
    if (Tools.isSet(annotations))
      for (IsotopePeakAnnotation iso : annotations) if (iso.getIsotopeOrder() == 0) return true;

    return false;
  }

  public static <T extends Peak> double[] asMz(List<T> peaks) {
    if (!Tools.isSet(peaks)) return null;
    double[] mz = new double[peaks.size()];

    for (int i = 0; i < peaks.size(); i++) mz[i] = peaks.get(i).getMz();

    return mz;
  }

  public static <T extends Peak> List<T> log10(List<T> peaks) {
    if (!Tools.isSet(peaks)) return peaks;

    List<T> out = new ArrayList<>(peaks.size());
    for (T p : peaks) {
      T p1 = (T) p.copy();
      p1.setIntensity(Math.log10(p1.getIntensity()));
      out.add(p1);
    }

    return out;
  }

  public static int countC12(PeakList peaks, int left, int right) {
    if (peaks == null || left < 0 || right < left) return 0;

    int c12 = 0;
    for (int i = left; i < right; i++)
      if (!hasC13(peaks.getAnnotations(i))) c12++;

    return c12;
  }

  public static <T extends Peak> double AbsIntensitySum(Collection<T> peaks) {
    if (!Tools.isSet(peaks)) return 0;

    double sum = 0d;
    for (T peak : peaks) sum += Math.abs(peak.getIntensity());

    return sum;
  }
  public static <T extends Peak> double AbsIntensitySum(T... peaks)
  {
    if (!Tools.isSet(peaks)) return 0;

    double sum = 0d;
    for (T peak : peaks) sum += Math.abs(peak.getIntensity());

    return sum;
  }
  public static <T extends Peak> double IntensitySum(T... peaks)
  {
    if (!Tools.isSet(peaks)) return 0;

    double sum = 0d;
    for (T peak : peaks) sum += peak.getIntensity();

    return sum;
  }
  public static <T extends Peak> boolean hasNegativeIntensity(T... peaks)
  {
    if (Tools.isSet(peaks))
      for (T peak : peaks) if (peak.getIntensity()<0) return true;

    return false;
  }
  public static int[] charges(PeakList peaks, int i) {
    List<Integer> zs = new ArrayList<>();
    if (peaks != null && i >= 0 && i < peaks.size()) {
      Collection<PeakAnnotation> annots = peaks.getAnnotations(i);
      if (Tools.isSet(annots))
        for (PeakAnnotation annot : annots)
          if (annot.getCharge() != 0) zs.add(annot.getCharge());
    }
    if (zs.size() > 0) {
      int[] zz = new int[zs.size()];
      for (int k = 0; k < zs.size(); k++) zz[k] = zs.get(k);
      return zz;
    }
    return new int[]{};
  }

  public static List<Peak> isolate(PeakList peaks, double left, double right) {
    if (peaks == null) return null;

    List<Peak> isolated = new ArrayList<>();
    for (int i = 0; i < peaks.size(); i++) {
      if (peaks.getMz(i) >= left && peaks.getMz(i) <= right) {
        isolated.add(new Peak(peaks.getMz(i), peaks.getIntensity(i), charges(peaks, i)));
      }
    }
    return isolated;
  }

  public static SortedMap<Double, Peak> toPrecursors(String comment) {
    String[] strs = Strs.split(comment, ';');
    if (Tools.isSet(strs)) {
      SortedMap<Double, Peak> precursors = new TreeMap<Double, Peak>();
      for (int i = 0; i < strs.length; i++) {
        String[] items = Strs.split(strs[i], '/');
        if (items != null && items.length > 2) {
          Peak pk = new Peak(Stats.toDouble(items[0]), Stats.toDouble(items[1]), Stats.toInt(items[2]));
          precursors.put(pk.getMz(), pk);
        }
      }
      return precursors;
    }
    return null;
  }
  public static ConcurrentSkipListMap<Double, Peak> toPrecursorsSkipListMap(String comment)
  {
    String[] strs = Strs.split(comment, ';');
    if (Tools.isSet(strs)) {
      ConcurrentSkipListMap<Double, Peak> precursors = new ConcurrentSkipListMap<>();
      for (int i = 0; i < strs.length; i++) {
        String[] items = Strs.split(strs[i], '/');
        if (items != null && items.length > 2) {
          Peak pk = new Peak(Stats.toDouble(items[0]), Stats.toDouble(items[1]), Stats.toInt(items[2]));
          precursors.put(pk.getMz(), pk);
        }
      }
      return precursors;
    }
    return null;
  }

  // is the putative MH supported by the precursors detected in the isolation window?
  public static AnnotatedPeak verifyCalcMH(int max_z, AnnotatedPeak calc, Tolerance tol, SortedMap<Double, Peak>... isolated_precursors)
  {
    // n opinion without additional information
    if (!Tools.isSet(isolated_precursors)) return calc;

    // now look for the best charge state as supported by the observed precursors
    int best = 0, isotopes = 0;
    calc.setIntensity(0d);
    calc.setVerifiedCharge(0);
    calc.setOriginalMz(0d);
    for (double z = max_z; z >= 1d; z--) {
      double mz = Peaks.MH2Mz(calc.getMz(), (int) z), ai = 0d, c12 = 0d;
      isotopes = 0;
      for (int c13 = 0; c13 <= z; c13++) {
        Map<Double, Peak> c = null;
        for (SortedMap<Double, Peak> precursors : isolated_precursors)
        {
          c = precursors.subMap(tol.getMin(mz + c13 * 1.003355d / z), tol.getMax(mz + c13 * 1.003355d / z));
          if (Tools.isSet(c)) break;
        }
        if (!Tools.isSet(c)) {
          isotopes = c13 - 1;
          break;
        }
        // accumulate the intensities
        ai += AbsIntensitySum(c.values());
        isotopes = c13;
        if (c12 == 0) c12 = Peaks.centroid(c.values());
      }
      if (isotopes > best || (isotopes == best && isotopes > 0 && ai > calc.getIntensity())) {
        best = isotopes;
        calc.setIntensity(ai);
        calc.setOriginalMz(c12).setProperty("isotopes", best).setVerifiedCharge((int) z);
      }
    }

    // no need to look for the true c12 since we starts with the calculated MH
    return calc;
  }
  public static AnnotatedPeak verifyCalcMH(int min_z, int max_z, AnnotatedPeak calc, Tolerance tol, ImmutableNavigableMap<Peak>... isolated_precursors)
  {
    // n opinion without additional information
    if (!Tools.isSet(isolated_precursors)) return calc;

    // now look for the best charge state as supported by the observed precursors
    int best = 0, isotopes = 0;
    calc.setIntensity(0d);
    calc.setVerifiedCharge(0);
    calc.setOriginalMz(0d);
    for (double z = max_z; z >= min_z; z--)
    {
      double mz = Peaks.MH2Mz(calc.getMz(), (int) z), ai = 0d, c12 = 0d;
      isotopes = 0;
      for (int c13 = 0; c13 <= z; c13++)
      {
        Peak pk = null; ImmutableNavigableMap<Peak> prec=null;
        for (ImmutableNavigableMap<Peak> precursors : isolated_precursors)
        {
          pk = precursors!=null?Peaks.query4peak(precursors, tol.getMin(mz+c13*1.003355d/z), tol.getMax(mz+c13*1.003355d/z)):null;
          if (pk!=null) { prec=precursors; break; }
        }
        if (pk==null)
        {
          isotopes = c13-1;
          break;
        }
        // accumulate the intensities
        ai += pk.getIntensity();
        isotopes = c13;
        if (c12 == 0) c12 = pk.getMz();
//        int[] c = null; ImmutableNavigableMap<Peak> prec=null;
//        for (ImmutableNavigableMap<Peak> precursors : isolated_precursors)
//        {
//          c = precursors!=null?precursors.query(tol.getMin(mz + c13 * 1.003355d / z), tol.getMax(mz + c13 * 1.003355d / z)):null;
//          if (c!=null) { prec=precursors; break; }
//        }
//        if (c==null)
//        {
//          isotopes = c13-1;
//          break;
//        }
//        // accumulate the intensities
//        ai += AbsIntensitySum(prec.fetchVals(c));
//        isotopes = c13;
//        if (c12 == 0) c12 = Peaks.centroid(prec.fetchVals(c), null,null);
//        c=null;
      }
      if (isotopes > best || (isotopes == best && isotopes > 0 && ai > calc.getIntensity()))
      {
        best = isotopes;
        calc.setIntensity(ai);
        calc.setOriginalMz(c12).setProperty("isotopes", best).setVerifiedCharge((int) z);
      }
    }

    // no need to look for the true c12 since we starts with the calculated MH
    return calc;
  }
  public static boolean query4index(ImmutableNavigableMap<PeakMatch> peaks, double k0, double k1, Collection<Integer> index)
  {
    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);

    if (start>=0)
    {
      int j0=-1, j1=-1;
      for (int k=start; k<peaks.getKeys().length; k++)
      {
        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
      }
      if (j0>=0 && j1>j0)
      {
        for (int j=j0; j<j1; j++) index.add(j);
        return true;
      }
    }

    return false;
  }

  public static Peak query4peak(ImmutableNavigableMap<Peak> peaks, double k0, double k1)
  {
    double ai=0, mz=0;
    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);

    if (start>=0)
    {
      int j0=-1, j1=-1;
      for (int k=start; k<peaks.getKeys().length; k++)
      {
        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
      }
      if (j0>=0 && j1>j0)
      {
        for (int j=j0; j<j1; j++)
        {
          ai+=Math.abs(peaks.getVals()[j].getIntensity());
          mz+=peaks.getKeys()[j];
        }
        mz/=(double )(j1-j0);
      }
    }

    return ai>0?(new Peak(mz, ai, 0)):null;
  }
  // doo we have good evidence for such charge state in our precursors
  synchronized public static boolean hasIsotopeEnvelop(
      Table<Integer, Integer, IsoEnvelope> sIsoEnveloped, int charge, double minDP, Tolerance tol,
      ImmutableNavigableMap<Peak>... isolated_precursors)
  {
    // n opinion without additional information
    if (!Tools.isSet(isolated_precursors)) return false;

    // now look for the best charge state as supported by the observed precursors
    for (ImmutableNavigableMap<Peak> precursors : isolated_precursors)
      if (precursors!=null && precursors.size()>0)
        for (Peak p : precursors.getVals())
        {
          double dp = scoreIsotope(sIsoEnveloped, p.getMz(), charge, tol, isolated_precursors);
          if (dp>minDP) return true;
        }

    return false;
  }
  // giving the m/z and charge, check the evidence
  public static double scoreIsotope(Table<Integer, Integer, IsoEnvelope> sIsoEnveloped, double mz, int charge, Tolerance tol, ImmutableNavigableMap<Peak>... isolated_precursors)
  {
    IsoEnvelope ev = Isotopes.estIsotopesByMz(mz, charge, 10d, 100d, sIsoEnveloped); int order=0;
    for (Peak iso : ev.getPredicted())
    {
      double ai=0d; // gather the ions from all precursors
      for (ImmutableNavigableMap<Peak> precursors : isolated_precursors)
        if (precursors!=null)
        {
          ai+=Peaks.query4ai(precursors, tol.getMin(iso.getMz()), tol.getMax(iso.getMz()));
//          int[] found = precursors.query(tol.getMin(iso.getMz()), tol.getMax(iso.getMz()));
//          if (Tools.isSet(found))
//            for (Peak p : precursors.fetchVals(found)) ai+=p.getIntensity();
        }

      if (ai>0) { ev.addIsotope(new Peak(iso.getMz(), ai)); order++; } else break;
    }
    return order>=charge?Similarity.dp(ev.getPredicted().subList(0,order), ev.getIsotopes().subList(0,order)):0d;
  }
  @Deprecated
  public static AnnotatedPeak verifyCalcMH(int max_z, AnnotatedPeak calc, Tolerance tol, ConcurrentSkipListMap<Double, Peak>... isolated_precursors)
  {
    // n opinion without additional information
    if (!Tools.isSet(isolated_precursors)) return calc;

    // now look for the best charge state as supported by the observed precursors
    int best = 0, isotopes = 0;
    calc.setIntensity(0d);
    calc.setVerifiedCharge(0);
    calc.setOriginalMz(0d);
    for (double z = max_z; z >= 1d; z--) {
      double mz = Peaks.MH2Mz(calc.getMz(), (int) z), ai = 0d, c12 = 0d;
      isotopes = 0;
      for (int c13 = 0; c13 <= z; c13++) {
        ConcurrentNavigableMap<Double, Peak> c = null;
        for (ConcurrentSkipListMap<Double, Peak> precursors : isolated_precursors)
        {
          c = precursors.subMap(tol.getMin(mz + c13 * 1.003355d / z), tol.getMax(mz + c13 * 1.003355d / z));
          if (c!=null && c.firstEntry()!=null) break;
        }
        if (c==null || c.firstEntry()==null || c.size()==0)
        {
          isotopes = c13-1;
          break;
        }
        // accumulate the intensities
        ai += AbsIntensitySum(c.values());
        isotopes = c13;
        if (c12 == 0) c12 = Peaks.centroid(c.values());
      }
      if (isotopes > best || (isotopes == best && isotopes > 0 && ai > calc.getIntensity())) {
        best = isotopes;
        calc.setIntensity(ai);
        calc.setOriginalMz(c12).setProperty("isotopes", best).setVerifiedCharge((int) z);
      }
    }

    // no need to look for the true c12 since we starts with the calculated MH
    return calc;
  }

  // is the putative MH supported by the precursors detected in the isolation window?
  public static double hasPrecursor(int z, Tolerance tol, Range<Double> isolation, ImmutableNavigableMap<Peak>... isolated_precursors)
  {
    double spacing = 1.00238d;

    // n opinion without additional information
    if (Tools.isSet(isolated_precursors))
      for (ImmutableNavigableMap<Peak> precursors : isolated_precursors)
        if (precursors!=null)
          for (double p0 : precursors.getKeys()) {
            if (p0 > isolation.upperEndpoint()) return 0;
            double maxmz = 0, c13 = 0;
            for (c13 = 1; c13 <= z; c13++)
            {
              if (precursors.query4counts(tol.getMin(p0+c13*spacing/(double) z), tol.getMax(p0+c13*spacing/(double) z))>0) break;
              maxmz = p0 + c13 * spacing / (double) z;
            }
            if (c13 > (z < 2 ? 1 : 2) && maxmz > isolation.lowerEndpoint()) return Peaks.toMH(p0,z);
          }

    // no need to look for the true c12 since we starts with the calculated MH
    return 0;
  }

  @Deprecated
  // is the putative MH supported by the precursors detected in the isolation window?
  public static double hasPrecursor(int z, Tolerance tol, Range<Double> isolation, SortedMap<Double, Peak>... isolated_precursors) {
    double spacing = 1.00238d;

    // n opinion without additional information
    if (Tools.isSet(isolated_precursors))
      for (SortedMap<Double, Peak> precursors : isolated_precursors)
        for (Double p0 : precursors.keySet()) {
          if (p0 > isolation.upperEndpoint()) return 0;
          double maxmz = 0, c13 = 0;
          for (c13 = 1; c13 <= z; c13++) {
            Map<Double, Peak> c = precursors.subMap(tol.getMin(p0 + c13 * spacing / (double) z), tol.getMax(p0 + c13 * spacing / (double) z));
            if (!Tools.isSet(c)) break;
            maxmz = p0 + c13 * spacing / (double) z;
          }
          if (c13 > (z < 2 ? 1 : 2) && maxmz > isolation.lowerEndpoint()) return p0;
        }

    // no need to look for the true c12 since we starts with the calculated MH
    return 0;
  }

  public static int getMostIntense(PeakList ms, Collection<Integer> peaks)
  {
    if (ms!=null && peaks!=null)
    {
      double base = 0;
      int pk = -1;
      for (Integer i : peaks)
        if (ms.getIntensity(i) > base)
        {
          base = ms.getIntensity(i);
          pk = i;
        }

      return pk;
    }
    return -1;
  }
  public static int getMostIntense(PeakList ms, Integer... peaks)
  {
    if (ms!=null && peaks!=null)
    {
      double base = 0;
      int pk = -1;
      for (Integer i : peaks)
        if (ms.getIntensity(i) > base)
        {
          base = ms.getIntensity(i);
          pk = i;
        }

      return pk;
    }
    return -1;
  }

  public static int find(SortedMap<Double, Integer> peaks, PeakList ms, double m, Tolerance tol)
  {
    int found = -1;
    if (m>0)
    {
      Map<Double, Integer> pk = peaks.subMap(tol.getMin(m), tol.getMax(m));
      // locate the most intense peak
      if (pk.size()>0) found = Peaks.getMostIntense(ms, pk.values());
      //Tools.dispose(  pk);
    }
    return found;
  }
  public static int find(ImmutableNavigableMap<Integer> peaks, PeakList ms, double m, Tolerance tol)
  {
    int found = -1;
    if (m>0)
    {
      int[] pks = peaks.query(tol.getMin(m), tol.getMax(m));
      // locate the most intense peak
      if (Tools.isSet(pks)) found = Peaks.getMostIntense(ms, peaks.fetchVals(pks));
      //Tools.dispose(  pk);
    }
    return found;
  }
  public static List<AnnotatedPeak> accumulate(List<AnnotatedPeak> A, List<AnnotatedPeak> B)
  {
    if (A==null) return B;
    if (B==null) return A;

    for (AnnotatedPeak pk : B)
      if (!A.contains(pk)) A.add(pk);

    Collections.sort(A);
    return A;
  }
  public static double[] MH2Mzs(double mh, int zlower, int zupper)
  {
    if (zupper<zlower) return new double[] {mh};

    double mzs[] = new double[zupper-zlower+1];
    for (int z=zlower; z<=zupper; z++)
      mzs[z-zlower] = MnH2MnH(mh, 1, z);

    return mzs;
  }
  public static AnnotatedPeak setProperties(AnnotatedPeak pk, TabFile f, String... keys) throws IOException
  {
    if (f!=null && Tools.isSet(keys))
      for (String key : keys)
        pk.setProperty(key, f.getDouble(key));

    return pk;
  }
  public static double query4ai(ImmutableNavigableMap<Peak> peaks, double k0, double k1)
  {
    double ai=0;
    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);

    if (start>=0)
    {
      int j0=-1, j1=-1;
      for (int k=start; k<peaks.getKeys().length; k++)
      {
        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
      }
      if (j0>=0 && j1>j0)
        for (int j=j0; j<j1; j++)
          ai+=Math.abs(peaks.getVals()[j].getIntensity());
    }

    return ai;
  }
}
