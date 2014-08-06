package org.ms2ms.test;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.data.Dataframe;
import org.ms2ms.runner.Aligner;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.util.Collections;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/17/14
 * Time: 6:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class DataframeTest extends TestAbstract
{
  Tolerance rttol=new PpmTolerance(1.5E5), mztol=new PpmTolerance(15);
  String root = "/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/";
  Dataframe evidences = null, peptides=null;

  @Before
  public void setUp()
  {
    //evidences = new Dataframe(root+"evidence1k.txt", '\t');
    evidences = new Dataframe(root+"evidence.txt", '\t');
    peptides  = new Dataframe(root+"allPeptides.txt", '\t');
  }
  @Test
  public void pivoting() throws Exception
  {
    // Dataframe pivot(String col, String val, Stats.Aggregator func, String... rows)
    Dataframe out = evidences.pivot("Raw file", "Retention time calibration", Stats.Aggregator.COUNT, "Modified sequence");
    System.out.println("\n" + out.display());
  }
  @Test
  public void splitting() throws Exception
  {
    // Dataframe pivot(String col, String val, Stats.Aggregator func, String... rows)
    Map<Object, Dataframe> outs = evidences.split("Raw file");

    for (Object obj : outs.keySet())
      System.out.println(obj.toString() + "\n" + outs.get(obj).display());
  }
  @Test
  public void interpolating() throws Exception
  {
    // Dataframe pivot(String col, String val, Stats.Aggregator func, String... rows)
    double xs[] = evidences.getDoubleCol("Retention time");
    double ys[] = evidences.getDoubleCol("Retention time calibration");
    double Xs[] = evidences.getDoubleCol("Calibrated retention time start");

    double[] Ys = Stats.interpolate(xs, ys, 0.3, Xs);

    evidences.addVar("interpolated", Ys);
    System.out.println("\n" + evidences.display());

    evidences.addVar("Calibrated RT", Stats.sum(evidences.getDoubleCol("Retention time"), evidences.getDoubleCol("Retention time calibration")));
  }
  @Test
  public void aligning() throws Exception
  {
    //Map<Object, Dataframe> outs = evidences.split("Raw file");
    Map<Object, Dataframe> outs = peptides.split("Raw file");

    Aligner aligner = new Aligner(new Tolerance[] {mztol, rttol}, "m/z", "Retention time");
    aligner.run(outs.values().toArray(new Dataframe[] {}));

    System.out.println("\n" + aligner.print());

    for (Object obj : outs.keySet())
      System.out.println(obj.toString() + "\n" + outs.get(obj).display());
  }
}
