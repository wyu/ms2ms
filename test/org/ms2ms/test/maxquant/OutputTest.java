package org.ms2ms.test.maxquant;

import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.io.MaxQuant;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.runner.Aligner;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Stats;

import java.io.FileWriter;
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

  @Before
  public void setUp()
  {
    // allPeptides100.txt  evidence100.txt  msms100.txt  msmsScans100.txt
//    evidences = Dataframe.readtable(root+"evidence1k.txt",    '\t').setTitle("evidence");
//    peptides  = Dataframe.readtable(root+"allPeptides1k.txt", '\t').setTitle("peptide");
//    msms      = Dataframe.readtable(root+"msms1k.txt",        '\t').setTitle("msms");
//    scans     = Dataframe.readtable(root+"msmsScans1k.txt",   '\t').setTitle("scan");

    evidences = Dataframe.readtable(root+"evidence.txt",    '\t');
    msms      = Dataframe.readtable(root+"msms.txt", '\t');
    scans     = Dataframe.readtable(root+"msmsScans.txt", '\t');

    evidences.removeCols(MaxQuant.sRmdEvidence);
    msms.removeCols(MaxQuant.sRmdMsms);
    scans.removeCols(MaxQuant.sRmdScan);
//    evidences.removeCols("Sequence","Length","Modifications","Oxidation (M) Probabilities","Oxidation (M) Score Diffs","Acetyl (Protein N-term)","Oxidation (M)","Missed cleavages","Proteins","Leading Proteins","Leading Razor Protein","TypeMatch time difference","Match m/z difference","Match q-value","Match score","PIF","Fraction of total spectrum","Base peak fraction","Reverse","Contaminant","Combinatorics","AIF MS/MS IDs","Oxidation (M) site IDs","Type");
//    msms.removeCols("Precursor","Sequence","Length","Missed cleavages","Modifications","Oxidation (M) Probabilities","Oxidation (M) Score Diffs","Acetyl (Protein N-term)","Oxidation (M)","Proteins","Fragmentation","Mass analyzer","Type","Score diff","Localization prob","Combinatorics","PIF","Fraction of total spectrum","Base peak fraction","Precursor Full ScanNumber","Precursor Intensity","Precursor Apex Fraction","Precursor Apex Offset","Precursor Apex Offset Time","Matches","Intensities","Mass Deviations [Da]","Mass Deviations [ppm]","Masses","Neutral loss level","ETD identification type","Reverse","All scores","All sequences","All modified sequences","Oxidation (M) site IDs","Scan type","Modified sequence","PEP","Score","Delta score","Protein group IDs");
//    scans.removeCols("Collision energy","Summations","Identified","MS/MS IDs","Sequence","Length","Mass analyzer","Parent intensity fraction","Fraction of total spectrum","Base peak fraction","Precursor full scan number","Precursor intensity","Precursor apex fraction","Precursor apex offset","Precursor apex offset time","Proteins","Score","Intens Comp Factor","CTCD Comp","RawOvFtT","AGC Fill","Modified sequence","PEP","Score","Delta score","Protein group IDs");
    // replace the "id" variable with "Evidence ID" so they can be joint
    evidences.renameCol("id", "Evidence ID");
    evidences.renameCol("m/z", "m/z-calibrated");
    evidences.renameCol("Mass", "Mass-predicted");
  }
  @Test
  // create a composite view of MSMS scans with peptide annotation if available
  public void toMsMsComposite() throws Exception
  {
    // force the ID columns to be categorical
         msms.init(Var.VarType.CATEGORICAL, "Evidence ID");
    evidences.init(Var.VarType.CATEGORICAL, "Evidence ID");
    // join the msms table with the evidence table, which contains the annotated LCMS features by the Evidence ID
    Dataframe annotations = Dataframe.merge(msms, evidences, true, "Evidence ID").setTitle("annot");
    // msms.txt contains only the annotated scans, msmsScans.txt has all the scans
    Dataframe out = Dataframe.merge(annotations, scans, true, "Raw file","Scan number");

    IOs.write(root + "mergedAnnotAll.txt", annotations.display("\t").toString());
    IOs.write(root + "mergedScansAll.txt", out.display("\t").toString());
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
