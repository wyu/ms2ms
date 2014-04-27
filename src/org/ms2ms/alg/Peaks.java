package org.ms2ms.alg;

import com.google.common.collect.Lists;
import org.expasy.mzjava.core.mol.AtomicSymbol;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.core.ms.spectrum.Peak;
import org.expasy.mzjava.proteomics.mol.AAMassCalculator;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.util.*;

/** Collection of algorithms pertaining to the MS Peak
 *
 * Created by wyu on 4/23/14.
 */
public class Peaks
{
  static class IntensityDesendComparator implements Comparator<Peak> { public int compare(Peak o1, Peak o2) { return o1!=null && o2!=null ? Double.compare(o2.getIntensity(), o1.getIntensity()):0; } }
  static class IntensityAscendComparator implements Comparator<Peak> { public int compare(Peak o1, Peak o2) { return o1!=null && o2!=null ? Double.compare(o1.getIntensity(), o2.getIntensity()):0; } }

  public static boolean isType(PepLibPeakAnnotation s, IonType... types)
  {
    IonType ion = s.getOptFragmentAnnotation().get().getIonType();
    if (s!=null && Tools.isSet(types))
      for (IonType t : types)
        if (t.equals(ion)) return true;

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
  public static double getmeanIntensity(Collection<Peak> msms)
  {
    if (msms==null || msms.size()==0) return 0;

    double sum = 0d;
    for (Peak p : msms) sum+=p.getIntensity();

    return sum/(double )msms.size();
  }
  public static double getBaseline(PeakList msms, double x0, double x1, int top, double val)
  {
    List<Peak> baselines = new ArrayList<Peak>();
    for (int i=0; i<msms.size(); i++)
      if (isValidIntensity(msms.getIntensity(i)) &&
          msms.getMz(i)>=x0 && msms.getMz(i)<=x1)
        baselines.add(new Peak(msms.getMz(i), msms.getIntensity(i), 1));

    Collections.sort(baselines, new IntensityDesendComparator());
    if (baselines.size()>top)
    {
      Collection<Double> pts = new ArrayList<Double>();
      for (int i=0; i<baselines.size()-top; i++) pts.add(baselines.get(i).getIntensity());
      return Stats.mean(pts) + Stats.stdev(pts)*val;
    }
    return -1d*getmeanIntensity(Lists.partition(baselines, 3).get(0));
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
  public static StringBuffer print(StringBuffer buf, PeakList<PepLibPeakAnnotation> msms)
  {
    if (msms == null || msms.size()==0) return null;

    if (buf==null) buf = new StringBuffer();
    buf.append("Precursor: ");
    buf = printIon(buf, msms.getPrecursor().getMz(),msms.getPrecursor().getIntensity(), msms.getPrecursor().getCharge());
    buf.append("\n");
    for (int i=0; i<msms.size(); i++)
    {
      buf = printIon(buf, msms.getMz(i),msms.getIntensity(i), 0); buf.append(", ");

      if (msms.getAnnotations(i)!=null)
        for (PepLibPeakAnnotation anno : msms.getAnnotations(i))
        {
          PepFragAnnotation f = anno.getOptFragmentAnnotation().get();
          double loss = f.getNeutralLoss().getMolecularMass();
          buf.append("+" + f.getCharge() +
            "," + f.getIonType() + (loss!=0?(loss>0?"+":"")+Math.round(loss):"") + ", " + Tools.d2s(f.getTheoreticalMz(), 4) + ", " + f.getFragment() + "; ");
        }
      buf.append("\n");
    }

    return buf;
  }
}
