package org.ms2ms.test;

import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import org.junit.Test;
import org.ms2ms.algo.PeptideLC;
import org.ms2ms.data.Features;
import org.ms2ms.data.ms.LcMsMsFeatures;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.utils.Strs;

public class RtPredictTest extends TestAbstract
{
//  String[] mods = Strs.toStringArray("Oxidation (M)","Unmodified","Deamidation (NQ)");
  String[] mods = Strs.toStringArray("Unmodified");

  @Test
  public void LsRPbyDevi()
  {
    Multimap<String, Double> devis = example2D().toResidualNET(null);
    PeptideLC.screenSVR(PeptideLC.SVRprobs(devis, 1f, false),1, 5,
        Range.closed(4d,6d), Range.closed(-5d,-3d), Range.closed(-5d,-3d));
  }
  @Test
  public void FitDeviByLsRP()
  {
    Multimap<String, Double> devis = example1D().toResidualNET(null);
    PeptideLC.fitSVR(PeptideLC.SVRprobs(devis, 0.8f, false),5,-4,-3);
  }
  @Test
  public void SSRcalc()
  {
    // first test of LsRT port
    Multimap<String, Double> peptides = example1D().toPeptideRt(10d, mods);
    PeptideLC.SSRCalc(peptides.keySet());
  }
  @Test
  public void LsRPv1()
  {
    // first test of LsRT port
    Multimap<String, Double> peptides = example1D().toPeptideRt(10d, mods);

//    PeptideLC.screenSVR(PeptideLC.SVRprobs(peptides, 1f, true),
//        1, Range.closed(2d,8d), Range.closed(-8d,-2d), Range.closed(-3d,-1d));
    PeptideLC.fitSVR(PeptideLC.SVRprobs(peptides, 0.8f, true),8, -6, -1);
  }
  @Test
  public void LCwidth()
  {
    LcMsMsFeatures plasma1d = example1D();

    System.out.println("Peptide\tRT\tLeft\tRight");
    for (Features F : plasma1d.getIons().values())
    {
      if (F.get(LcMsMsFeatures.COL_PEPTIDE)!=null)
      {
        Range rt = (Range)F.get(LcMsMsFeatures.COL_RT_BOUND);
        System.out.println(F.getStr(LcMsMsFeatures.COL_SEQ)+"\t"+F.getFloat(LcMsMsFeatures.COL_RT)+"\t"+
            rt.lowerEndpoint()+"\t"+rt.upperEndpoint());
      }
    }
  }
  private LcMsMsFeatures example2D()
  {
    return new LcMsMsFeatures(
        "/Users/kfvf960/OneDrive - AZCollaboration/contrib/data/Plasma.Mann/txt.20171226",
        new OffsetPpmTolerance(15), "20150227_QEp1_LC7_GaPi_SA_PlasmaDepleted_.*","20141029_QEp1_LC7_NiKu_SA_Plasma_male1_.*");
  }
  private LcMsMsFeatures example1D()
  {
    return new LcMsMsFeatures(
        "/Users/kfvf960/OneDrive - AZCollaboration/contrib/data/Plasma.Mann/txt.20171226",
        new OffsetPpmTolerance(15), "20150805_QEp1_LC7_PhGe_SA_100min_2.*");
  }
//  // Multimap<String, Double> peptides = example1D();
//  private Multimap<String, Double> deviFromSSRcalc(Multimap<String, Double> peptides, Multimap<String, Double> H)
//  {
//    // first test oRf LsRT port
//    Map<String, Double> predicted = PeptideLC.SSRCalc(peptides.keySet());
//
//    // no more than 15, up to the size-2
//    List<WeightedObservedPoint> Rs = new ArrayList<>();
//    for (String key : predicted.keySet())
//      for (Double rt : peptides.get(key))
//        Rs.add(new WeightedObservedPoint(1, rt, predicted.get(key)));
//
//    // quadratic fit is not suitable because of the possibility of curving back up at higher score. Even though the fit is often better.
//    Fitted fit = new Fitted().fit(1, Rs);
//
//    // convert the observed RT to predicted H
//    Multimap<String, Double> devi = HashMultimap.create();
//    for (String key : predicted.keySet())
//      for (Double rt : peptides.get(key))
//      {
//        double Hi = fit.polynomial(rt);
//        if (H!=null) H.put(key, Hi);
//        devi.put(key,Hi-predicted.get(key));
//      }
//
//    return devi;
//  }
//  private Multimap<String, Double> example1D()
//  {
//    LcMsMsFeatures plasma1d = new LcMsMsFeatures(
//        "/Users/kfvf960/OneDrive - AZCollaboration/contrib/data/Plasma.Mann/txt.20171226",
//        new OffsetPpmTolerance(15), "20150805_QEp1_LC7_PhGe_SA_100min_2.*");
//
//    Multimap<String, Double> peptide_rt = HashMultimap.create();
//    for (Features F : plasma1d.getIons().values())
//    {
//      if (F.get(LcMsMsFeatures.COL_PEPTIDE)!=null)
//        peptide_rt.put(F.getStr(LcMsMsFeatures.COL_SEQ), F.getDouble(LcMsMsFeatures.COL_RT));
//    }
//
//    return peptide_rt;
//  }
}
