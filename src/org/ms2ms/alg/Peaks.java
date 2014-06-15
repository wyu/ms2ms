package org.ms2ms.alg;

import com.google.common.collect.Lists;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.DoublePeakList;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.core.ms.spectrum.Peak;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.util.*;

/** Collection of algorithms pertaining to the MS Peak
 *
 * Created by wyu on 4/23/14.
 */
public class Peaks
{
  public static final String OBT_CID    = "OrbitrapCID";
  public static final String LTQ_CID    = "LTQCID";
  public static final String OBT_HR_CID = "OrbitrapCIDHr";
  public static final String OBT_HCD    = "OrbitrapHCD";
  public static final String QTOF       = "QTOF";

  public static final String CID       = "cid";
  public static final String HCD       = "hcd";
  public static final String ETD       = "etd";

  static class IntensityDesendComparator implements Comparator<Peak> { public int compare(Peak o1, Peak o2) { return o1!=null && o2!=null ? Double.compare(o2.getIntensity(), o1.getIntensity()):0; } }
  static class IntensityAscendComparator implements Comparator<Peak> { public int compare(Peak o1, Peak o2) { return o1!=null && o2!=null ? Double.compare(o1.getIntensity(), o2.getIntensity()):0; } }

  public static boolean isType(PepLibPeakAnnotation s, IonType... types)
  {
    try
    {
      IonType ion = s.getOptFragmentAnnotation().get().getIonType();
      if (s!=null && Tools.isSet(types))
        for (IonType t : types)
          if (t.equals(ion)) return true;
    }
    catch (IllegalStateException e)
    {

    }

    return false;
  }
  public static boolean hasType(Collection<PepLibPeakAnnotation> s, IonType... types)
  {
    if (Tools.isSet(s))
      for (PepLibPeakAnnotation A : s)
        if (isType(A, types)) return true;

    return false;
  }

  /** Count the number of valid peaks between x1 and x2
   *
   * @param ms is a PeakList object. zero or negative intensity indicates invalid peak
   * @param x1, x2 define the m/z bound where the counting is to take place
   * @return the counts of peaks bound by (x1, x2) that's above 0
   */
  public static int countValid(PeakList ms, double x1, double x2)
  {
    if (ms==null || ms.size()==0) return 0;

    int counts=0;
    for (int i=0; i<ms.size(); i++)
      if (ms.getMz(i)>=x1 && ms.getMz(i)<=x2 && isValidIntensity(ms.getIntensity(i))) counts++;

    return counts;
  }
  public static double getMinIntensity(PeakList msms) { return getMinIntensity(msms, Double.MIN_VALUE, Double.MAX_VALUE); }
  public static double getMinIntensity(PeakList msms, double x0, double x1)
  {
    if (msms==null || msms.size()==0) return 0;

    double baseline = Double.MAX_VALUE;
    for (int i=0; i<msms.size(); i++)
      if (msms.getMz(i)>=x0 && msms.getMz(i)<=x1 && msms.getIntensity(i)<baseline) baseline=msms.getIntensity(i);

    return baseline;
  }
  public static double getmeanIntensity(PeakList msms) { return getmeanIntensity(msms, Double.MIN_VALUE, Double.MAX_VALUE); }
  public static double getmeanIntensity(PeakList msms, double x0, double x1)
  {
    if (msms==null || msms.size()==0) return 0;

    double sum = 0d;
    for (int i=0; i<msms.size(); i++)
      if (msms.getMz(i)>=x0 && msms.getMz(i)<=x1) sum+=msms.getIntensity(i);

    return sum/(double )msms.size();
  }
  public static double getMeanIntensity(Collection<Peak> msms)
  {
    if (msms==null || msms.size()==0) return 0;

    double sum = 0d;
    for (Peak p : msms) sum+=p.getIntensity();

    return sum/(double )msms.size();
  }
  public static double getBaseline(PeakList msms, double x0, double x1, int top, boolean exclude_precursor)
  {
    List<Peak> baselines = new ArrayList<Peak>();
    for (int i=0; i<msms.size(); i++)
      if (isValidIntensity(msms.getIntensity(i)) &&
        (!exclude_precursor || !hasType(msms.getAnnotations(i), IonType.p)) &&
          msms.getMz(i)>=x0 && msms.getMz(i)<=x1)
        baselines.add(new Peak(msms.getMz(i), msms.getIntensity(i), 1));

    Collections.sort(baselines, new IntensityAscendComparator());
    if (baselines.size()>top)
    {
      Collection<Double> pts = new ArrayList<Double>();
      for (int i=0; i<baselines.size()-top; i++) pts.add(baselines.get(i).getIntensity());
      return Stats.mean(pts) + Stats.stdev(pts);
    }
    return Tools.isSet(baselines) ? -1d* getMeanIntensity(Lists.partition(baselines, 3).get(0)) : 0;
  }
  public static boolean isValidIntensity(double s) { return s>0; }

  public static StringBuffer printIon(StringBuffer buf, double mz, double ai, int z)
  {
    if (buf==null) buf = new StringBuffer();
    buf.append(Tools.d2s(mz, 4) + ", " + Tools.d2s(ai, 1) + ", +" + z);
    return buf;
  }
  public static StringBuffer print(StringBuffer buf, Peak p)
  {
    return printIon(buf, p.getMz(), p.getIntensity(), p.getCharge());
  }
  public static StringBuffer print(StringBuffer buf, AnnotatedPeak p)
  {
    buf = printIon(buf, p.getMz(), p.getIntensity(), p.getCharge());
    buf.append(", " + Tools.d2s(p.getSNR(), 2));
    return buf;
  }
  public static StringBuffer printAnnot(StringBuffer buf, Collection<PepLibPeakAnnotation> annotations)
  {
    if (annotations!=null)
    {
      if (buf==null) buf = new StringBuffer();
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
  public static StringBuffer print(StringBuffer buf, PeakList<PepLibPeakAnnotation> msms, boolean annot)
  {
    if (msms == null || msms.size()==0) return null;

    if (buf==null) buf = new StringBuffer();
    buf.append("Precursor: ");
    buf = printIon(buf, msms.getPrecursor().getMz(),msms.getPrecursor().getIntensity(), msms.getPrecursor().getCharge());
    buf.append("\n");
    for (int i=0; i<msms.size(); i++)
    {
      // skip the peak if annotation is required
      if (annot && msms.getAnnotations(i)==null) continue;

      buf = printIon(buf, msms.getMz(i),msms.getIntensity(i), 0); buf.append(", ");
      if (annot) buf = printAnnot(buf, msms.getAnnotations(i));
//
//      if (msms.getAnnotations(i)!=null)
//        for (PepLibPeakAnnotation anno : msms.getAnnotations(i))
//        {
//          PepFragAnnotation f = anno.getOptFragmentAnnotation().get();
//          double loss = f.getNeutralLoss().getMolecularMass();
//          buf.append("+" + f.getCharge() +
//            "," + f.getIonType() + (loss!=0?(loss>0?"+":"")+Math.round(loss):"") + ", " + Tools.d2s(f.getTheoreticalMz(), 4) + ", " + f.getFragment() + ", " + f.getIsotopeComposition() + "; ");
//        }
      buf.append("\n");
    }

    return buf;
  }
  public static boolean isCtermFrag(IonType s)
  {
    return s!=null?(s.equals(IonType.y) || s.equals(IonType.x) || s.equals(IonType.z) || s.equals(IonType.w)):false;
  }
  public static boolean isNtermFrag(IonType s)
  {
    return s!=null?(s.equals(IonType.b) || s.equals(IonType.a) || s.equals(IonType.c) || s.equals(IonType.d)):false;
  }
  public static double toMass(double mz, int z) { return (mz-1.00078)*z; }
  public static double toMass(Peak p)           { return toMass(p.getMz(), p.getCharge()); }
  public static double toPPM(double m0, double m1) { return 1E6*(m1-m0)/m0; }

  /** makeup a peaklist using string shorthand
   *
   * @param mz and z are the m/z value and charge state of the precursor.
   * @param frags variable number of fragment ions. e.g. "334.5", "562/23", "mz/ai". Only the mz is required
   * @return
   */
  public static PeakList<PepLibPeakAnnotation> newPeakList(double mz, int z, String... frags)
  {
    PeakList<PepLibPeakAnnotation> spec = new DoublePeakList<PepLibPeakAnnotation>();
    spec.setPrecursor(new Peak(mz, 0d, z));
    // go thro the fragment ions if any
    for (String f : frags)
    {
      try
      {
        String[] strs = f.split("/");
        spec.add(Double.valueOf(strs[0]), strs.length>2?Double.valueOf(strs[2]):0d);
      }
      catch (Exception e) { e.printStackTrace(); }
    }
    return spec;
  }

  /** Parse the precursor specification from the web form
   *
   * @param line : "501.1/2 5051/3 505.6"
   * @param z : the default charge if not specified with the m/z above
   * @return array of precursor peaks
   */
  public static Peak[] newPeaks(String line, int z)
  {
    if (!Tools.isSet(line)) return null;
    String[] strs = line.split(";|,|\\s");
    Peak[]  peaks = new Peak[strs.length];
    for (int i=0; i<strs.length; i++)
    {
      String[] items = strs[i].split("/");
      peaks[i] = new Peak(Double.valueOf(items[0]), 0d, items.length>1?Integer.valueOf(items[1]):z);
    }
    return peaks;
  }
//  public static String toString(Peak... peaks)
//  {
//    String out=null;
//    if (peaks!=null)
//      for (Peak p : peaks)
//        out = Tools.extend(out, Tools.d2s(p.getMz(), 5)p.toString(), "; ");
//
//    return out;
//  }
/*
  public static AnnotatedSpectrum newPeakList(String frags)
  {
    PeakList<PepLibPeakAnnotation> spec = new DoublePeakList<PepLibPeakAnnotation>();
    try
    {
      ions.setPrecursor(new Peak(Double.valueOf(pmz), 0, pz));
      String[] fs=frags.split(";|,|\\s");
      for (String f : fs)
        ions.add(Double.valueOf(f), 0d);
    }
    catch (Exception e) {}
    return ions;
  }
*/
}

