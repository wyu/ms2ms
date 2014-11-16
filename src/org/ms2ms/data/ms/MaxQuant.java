package org.ms2ms.data.ms;

import org.ms2ms.io.MsReaders;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Stats;

import java.util.Map;

/** Information and algorithms regarding Maxquant output
 *
 * ** Copyright 2014-2015 ms2ms.org
 * Author: wyu
 * Date:   11/6/14
 */
public class MaxQuant
{
  public static final String V_RT     = "RT-MS2.scan"; // not calibrated
  public static final String V_MZ     = "m/z"; // not calibrated
  public static final String V_TIC    = "Total ion current";
  public static final String V_OFFSET = "FilePointer";

  public static String[] sRmdEvidence = {"Sequence","Length","Modifications","Oxidation (M) Probabilities","Oxidation (M) Score Diffs","Acetyl (Protein N-term)","Oxidation (M)","Missed cleavages","Proteins","Leading Proteins","Leading Razor Protein","TypeMatch time difference","Match m/z difference","Match q-value","Match score","PIF","Fraction of total spectrum","Base peak fraction","Reverse","Contaminant","Combinatorics","AIF MS/MS IDs","Oxidation (M) site IDs","Type","Protein group IDs","Match time difference","Scan event number","Scan index"};
  public static String[] sRmdMsms     = {"Precursor","Sequence","Length","Missed cleavages","Modifications","Oxidation (M) Probabilities","Oxidation (M) Score Diffs","Acetyl (Protein N-term)","Oxidation (M)","Proteins","Fragmentation","Mass analyzer","Type","Score diff","Localization prob","Combinatorics","PIF","Fraction of total spectrum","Base peak fraction","Precursor Full ScanNumber","Precursor Intensity","Precursor Apex Fraction","Precursor Apex Offset","Precursor Apex Offset Time","Matches","Intensities","Mass Deviations [Da]","Mass Deviations [ppm]","Masses","Neutral loss level","ETD identification type","Reverse","All scores","All sequences","All modified sequences","Oxidation (M) site IDs","Scan type","Modified sequence","PEP","Score","Delta score","Protein group IDs","Scan event number","Scan index"};
  public static String[] sRmdScan     = {"Collision energy","Summations","Identified","MS/MS IDs","Sequence","Length","Mass analyzer","Parent intensity fraction","Fraction of total spectrum","Base peak fraction","Precursor full scan number","Precursor intensity","Precursor apex fraction","Precursor apex offset","Precursor apex offset time","Proteins","Score","Intens Comp Factor","CTCD Comp","RawOvFtT","AGC Fill","Modified sequence","PEP","Score","Delta score","Protein group IDs","Scan event number","Scan index"};

  private String    mResultDir, mRawDir, mSpCache;
  private Dataframe mSummary, mMsMs;

  public MaxQuant() { super(); }
  public MaxQuant(Dataframe msms) { super(); mMsMs=msms; }
  public MaxQuant(String result, String raw) { mResultDir=result; mRawDir=raw; }

  public void init()
  {
    // grab the summary first
    mSummary = Dataframe.readtable(mResultDir+"summary.txt",'\t').setTitle("summary");
    mSummary = mSummary.subset("Raw file!='Total'");
    mSpCache = "cache"+System.nanoTime();

    Dataframe survey = MsReaders.surveyMzXML(mSummary.getStrCol("Raw file"), mRawDir, mSpCache, null, 2);

    IOs.write(mResultDir+"scan_survey.txt", survey.display("\t", "").toString());

    // read the tables of MS/MS scans
    mMsMs = readMsMsWithAnnotations();
    Dataframe offsets = Dataframe.merge(mMsMs, survey, true, true, "Raw file", "Scan number").setTitle("offsets");

    IOs.write(mResultDir+"/composite_scans.txt", offsets.display("\t", "").toString());
  }
  public Dataframe readMsMsWithAnnotations() { return readMsMsWithAnnotations(mResultDir); }
  public Dataframe readMsMsWithAnnotations(String root)
  {
    mResultDir = root;

    Dataframe evidences = Dataframe.readtable(root+"evidence.txt",  '\t').removeCols(MaxQuant.sRmdEvidence).setTitle("evidence");
    Dataframe msms      = Dataframe.readtable(root+"msms.txt",      '\t').removeCols(MaxQuant.sRmdMsms).setTitle("msms");
    Dataframe scans     = Dataframe.readtable(root+"msmsScans.txt", '\t').removeCols(MaxQuant.sRmdScan).setTitle("scan");

    // replace the "id" variable with "Evidence ID" so they can be joint
    evidences.renameCol("id", "Evidence ID");
    evidences.renameCol("m/z", "m/z-calibrated");
    evidences.renameCol("Mass", "Mass-predicted");
    evidences.renameCol("Retention time", "RT-feature-calibrated");

     msms.renameCol("Retention time", "RT-MS2");
    scans.renameCol("Retention time", "RT-MS2");

    // force the ID columns to be categorical
    msms.init(     Var.VarType.CATEGORICAL, "Evidence ID");
    evidences.init(Var.VarType.CATEGORICAL, "Evidence ID");
    // join the msms table with the evidence table, which contains the annotated LCMS features by the Evidence ID
    Dataframe annotations = Dataframe.merge(msms, evidences, true, true, "Evidence ID").setTitle("annot");
    // msms.txt contains only the annotated scans, msmsScans.txt has all the scans
    mMsMs = Dataframe.merge(annotations, scans, true, true, "Raw file","Scan number").setTitle("msms+");

    IOs.write(root + "mergedAnnots.txt", annotations.display().toString());
    IOs.write(root + "mergedScans.txt", mMsMs.display().toString());
    return mMsMs;
  }
  public void imputeRT(Dataframe data, float rttol) throws Exception
  {
    if (data==null || !data.hasVars("Raw file","Retention time","Retention time calibration","Calibrated retention time start")) return;

    Map<Object, Dataframe> raw_data = data.split("Raw file");
    // go thro the run separately
    for (Dataframe d : raw_data.values())
    {
      double xs[] = d.getDoubleCol("Retention time");
      double ys[] = d.getDoubleCol("Retention time calibration");
      double Xs[] = d.getDoubleCol("Calibrated retention time start");

      double[] Ys = Stats.interpolate(xs, ys, rttol, Xs);
      d.addVar("interpolated", Ys);
      d.addVar("Calibrated RT", Stats.sum(d.getDoubleCol("Retention time"), d.getDoubleCol("Retention time calibration")));
    }
  }
}
