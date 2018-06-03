package org.ms2ms.test;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.junit.Test;
import org.ms2ms.algo.PeptideLC;
import org.ms2ms.data.Features;
import org.ms2ms.data.ms.LcMsMsFeatures;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.math.Fitted;
import org.ms2ms.utils.Tools;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RtPredictTest extends TestAbstract
{
  @Test
  public void LsRPbyDevi()
  {
    Multimap<String, Double> devis = example1D().toResidualNET(null);
    PeptideLC.screenSVR(PeptideLC.SVRprobs(devis, 1f, false),
        Range.closed(-8d,8d), Range.closed(-8d,8d), Range.closed(-3d,-1d));
  }
  @Test
  public void FitDeviByLsRP()
  {
    Multimap<String, Double> devis = example1D().toResidualNET(null);
    PeptideLC.fitSVR(PeptideLC.SVRprobs(devis, 0.8f, false),6,-5,-3);
  }
  @Test
  public void SSRcalc()
  {
    // first test of LsRT port
    Multimap<String, Double> peptides = example1D().toPeptideRt();
    PeptideLC.SSRCalc3(peptides.keySet());
  }
  @Test
  public void LsRPv1()
  {
    // first test of LsRT port
    Multimap<String, Double> peptides = example1D().toPeptideRt();

    PeptideLC.screenSVR(PeptideLC.SVRprobs(peptides, 1f, true),
        Range.closed(-8d,8d), Range.closed(-8d,8d), Range.closed(-3d,-1d));
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
//    Map<String, Double> predicted = PeptideLC.SSRCalc3(peptides.keySet());
//
//    // no more than 15, up to the size-2
//    List<WeightedObservedPoint> Rs = new ArrayList<>();
//    for (String peptide : predicted.keySet())
//      for (Double rt : peptides.get(peptide))
//        Rs.add(new WeightedObservedPoint(1, rt, predicted.get(peptide)));
//
//    // quadratic fit is not suitable because of the possibility of curving back up at higher score. Even though the fit is often better.
//    Fitted fit = new Fitted().fit(1, Rs);
//
//    // convert the observed RT to predicted H
//    Multimap<String, Double> devi = HashMultimap.create();
//    for (String peptide : predicted.keySet())
//      for (Double rt : peptides.get(peptide))
//      {
//        double Hi = fit.polynomial(rt);
//        if (H!=null) H.put(peptide, Hi);
//        devi.put(peptide,Hi-predicted.get(peptide));
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
