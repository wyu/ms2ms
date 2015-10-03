package org.ms2ms.test.maxquant;

import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.alg.MsStats;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.runner.Aligner;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;

import java.util.Map;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   11/3/14
 */
public class OutputTest extends TestAbstract
{
  Tolerance rttol=new PpmTolerance(1.5E5), mztol=new PpmTolerance(15);
  String root = "/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/";
  Dataframe evidences = null, peptides=null, msms=null, scans=null;
  MaxQuant mq = null;

  @Before
  public void setUp()
  {
/*
    // allPeptides100.txt  evidence100.txt  msms100.txt  msmsScans100.txt
    evidences = Dataframe.readtable(root+"evidence1k.txt",    '\t', false).setTitle("evidence");
    peptides  = Dataframe.readtable(root+"allPeptides1k.txt", '\t', false).setTitle("peptide");
    msms      = Dataframe.readtable(root+"msms1k.txt",        '\t', false).setTitle("msms");
    scans     = Dataframe.readtable(root+"msmsScans1k.txt",   '\t', false).setTitle("scan");

    evidences.removeCols(MaxQuant.sRmdEvidence);
    msms.removeCols(MaxQuant.sRmdMsms);
    scans.removeCols(MaxQuant.sRmdScan);

    evidences.renameCol("id", "Evidence ID");
    evidences.renameCol("m/z", "m/z-calibrated");
    evidences.renameCol("Mass", "Mass-predicted");
*/
  }
  @Test
  public void prepareMSnBin() throws Exception
  {
    // prepare the binary dump of MS/MS and composite dataframe
    mq = new MaxQuant(root, "/media/data/test/mzXML/");

//    Dataframe msms = Dataframe.readtable(root+"mergedScans.txt",    '\t').setTitle("msms+"),
    Dataframe msms = mq.readMsMsWithAnnotations(),
        survey = Dataframe.readtable(root+"scan_survey.txt",    '\t').setTitle("survey"),
        offsets = Dataframe.merge(msms, survey, true, true, "Raw file", "Scan number").setTitle("offsets");

    IOs.write(root+"composite_scans.txt", offsets.display().toString());
  }
  @Test
  public void RtCalibration() throws Exception
  {
    prepare();

    // Dataframe pivot(String col, String val, MsStats.Aggregator func, String... rows)
    double xs[] = evidences.getDoubleCol("Retention time", true);
    double ys[] = evidences.getDoubleCol("Retention time calibration", true);
    double Xs[] = evidences.getDoubleCol("Calibrated retention time start", true);

    double[] Ys = MsStats.interpolate(xs, ys, 0.3, Xs);

    evidences.addVar("interpolated", Ys);
    System.out.println("\n" + evidences.display());

    evidences.addVar("Calibrated RT", MsStats.sum(evidences.getDoubleCol("Retention time", true), evidences.getDoubleCol("Retention time calibration", true)));
  }
  @Test
  public void mergeMQnSurvey() throws Exception
  {
    prepare();

    mq = new MaxQuant(root, "/media/data/test/mzXML/");

//    Dataframe msms = Dataframe.readtable(root+"mergedScans.txt",    '\t').setTitle("msms+"),
    Dataframe msms = mq.readMsMsWithAnnotations(),
            survey = Dataframe.readtable(root+"scan_survey.txt",    '\t').setTitle("survey"),
           offsets = Dataframe.merge(msms, survey, true, true, "Raw file", "Scan number").setTitle("offsets");

    IOs.write(root+"composite_scans.txt", offsets.display().toString());
  }
  @Test
  public void newMQ() throws Exception
  {
    prepare();

    mq = new MaxQuant(root, "/media/data/test/mzXML/");
    mq.init();
  }
  @Test
  public void summaryFile() throws Exception
  {
    prepare();

    Dataframe summary = Dataframe.readtable(root+"summary.txt",    '\t').setTitle("summary");
    summary = summary.subset("Raw file!='Total'");
    System.out.println(summary.asVar("Raw file").getFactors().size());
  }
  @Test
  // create a composite view of MSMS scans with peptide annotation if available
  public void toMsMsComposite() throws Exception
  {
    prepare();

    // force the ID columns to be categorical
         msms.init(Var.VarType.CATEGORICAL, "Evidence ID");
    evidences.init(Var.VarType.CATEGORICAL, "Evidence ID");
    // join the msms table with the evidence table, which contains the annotated LCMS features by the Evidence ID
    Dataframe annotations = Dataframe.merge(msms, evidences, true, true, "Evidence ID").setTitle("annot");
    // msms.txt contains only the annotated scans, msmsScans.txt hasAccession all the scans
    Dataframe out = Dataframe.merge(annotations, scans, true, true, "Raw file","Scan number");

    IOs.write(root + "mergedAnnotAll.txt", annotations.display("\t", "").toString());
    IOs.write(root + "mergedScansAll.txt", out.display("\t", "").toString());
  }
  @Test
  public void interpolating() throws Exception
  {
    prepare();

    // Dataframe pivot(String col, String val, MsStats.Aggregator func, String... rows)
    double xs[] = evidences.getDoubleCol("Retention time", true);
    double ys[] = evidences.getDoubleCol("Retention time calibration", true);
    double Xs[] = evidences.getDoubleCol("Calibrated retention time start", true);

    double[] Ys = MsStats.interpolate(xs, ys, 0.3, Xs);

    evidences.addVar("interpolated", Ys);
    System.out.println("\n" + evidences.display());

    evidences.addVar("Calibrated RT", MsStats.sum(evidences.getDoubleCol("Retention time", true), evidences.getDoubleCol("Retention time calibration", true)));
  }
  @Test
  public void aligning() throws Exception
  {
    prepare();

    //Map<Object, Dataframe> outs = evidences.split("Raw file");
    Map<Object, Dataframe> outs = peptides.split("Raw file");

    Aligner aligner = new Aligner(new Tolerance[] {mztol, rttol}, "m/z", "Retention time");
    aligner.run(outs.values().toArray(new Dataframe[] {}));

    System.out.println("\n" + aligner.print());

    for (Object obj : outs.keySet())
      System.out.println(obj.toString() + "\n" + outs.get(obj).display());
  }
  protected void prepare()
  {
    // allPeptides100.txt  evidence100.txt  msms100.txt  msmsScans100.txt
    evidences = Dataframe.readtable(root+"evidence1k.txt",    '\t', false).setTitle("evidence");
    peptides  = Dataframe.readtable(root+"allPeptides1k.txt", '\t', false).setTitle("peptide");
    msms      = Dataframe.readtable(root+"msms1k.txt",        '\t', false).setTitle("msms");
    scans     = Dataframe.readtable(root+"msmsScans1k.txt",   '\t', false).setTitle("scan");

    evidences.removeCols(MaxQuant.sRmdEvidence);
    msms.removeCols(MaxQuant.sRmdMsms);
    scans.removeCols(MaxQuant.sRmdScan);

    evidences.renameCol("id", "Evidence ID");
    evidences.renameCol("m/z", "m/z-calibrated");
    evidences.renameCol("Mass", "Mass-predicted");
  }
}
