package org.ms2ms.data.ms;

import org.ms2ms.algo.MsStats;
import org.ms2ms.io.MsReaders;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.utils.IOs;

import java.util.Map;

/** Information and algorithms regarding Maxquant output
 *
 * ** Copyright 2014-2015 ms2ms.org
 * Author: wyu
 * Date:   11/6/14
 */
public class MaxQuant extends LcMsMsDataset
{
  public static final String V_RT      = "RT-MS2.scan"; // not calibrated
  public static final String V_MZ      = "m/z"; // not calibrated
  public static final String V_TIC     = "Total ion current";
  public static final String V_OFFSET  = "FilePointer";
  public static final String V_MSEQ    = "Modified sequence";
  public static final String V_SEQ     = "Sequence";
  public static final String V_Z       = "Charge";
  public static final String V_CLUSTER = "ClusterID";
  public static final String V_RAWFILE = "Raw file";
  public static final String V_SCAN    = "Scan number";
  public static final String V_ID_EVI  = "Evidence ID";

  public static String[] sRmdEvidence = {"Length","Modifications","Oxidation (M) Probabilities","Oxidation (M) Score Diffs","Acetyl (Protein N-term)","Oxidation (M)","Missed cleavages","Proteins","Leading Proteins","Leading Razor Protein","TypeMatch time difference","Match m/z difference","Match q-value","Match score","PIF","Fraction of total spectrum","Base peak fraction","Reverse","Contaminant","Combinatorics","AIF MS/MS IDs","Oxidation (M) site IDs","Type","Protein group IDs","Match time difference","Scan event number","Scan index"};
  public static String[] sRmdMsms     = {"Precursor","Sequence","Length","Missed cleavages","Modifications","Oxidation (M) Probabilities","Oxidation (M) Score Diffs","Acetyl (Protein N-term)","Oxidation (M)","Proteins","Fragmentation","Mass analyzer","Type","Score diff","Localization prob","Combinatorics","PIF","Fraction of total spectrum","Base peak fraction","Precursor Full ScanNumber","Precursor Intensity","Precursor Apex Fraction","Precursor Apex Offset","Precursor Apex Offset Time","Matches","Intensities","Mass Deviations [Da]","Mass Deviations [ppm]","Masses","Neutral loss level","ETD identification type","Reverse","All scores","All sequences","All modified sequences","Oxidation (M) site IDs","Scan type","Modified sequence","PEP","Score","Delta score","Protein group IDs","Scan event number","Scan index"};
  public static String[] sRmdScan     = {"Collision energy","Summations","Identified","MS/MS IDs","Sequence","Length","Mass analyzer","Parent intensity fraction","Fraction of total spectrum","Base peak fraction","Precursor full scan number","Precursor intensity","Precursor apex fraction","Precursor apex offset","Precursor apex offset time","Proteins","Score","Intens Comp Factor","CTCD Comp","RawOvFtT","AGC Fill","Modified sequence","PEP","Score","Delta score","Protein group IDs","Scan event number","Scan index"};

//  private String    mResultDir/*, mRawDir, mSpCache*/;
  private Dataframe mMsMs;

  public MaxQuant() { super(); }
  public MaxQuant(String s) { super(s); }
  public MaxQuant(Dataframe msms) { super(); mMsMs=msms; }
  public MaxQuant(String result, String raw) { mResultRoot=result; mRawfileRoot=raw; }

  public void init()
  {
    // grab the summary first
    mSummary = Dataframe.readtable(mResultRoot+"/summary.txt",'\t').setTitle("summary");
    mSummary = mSummary.subset("Raw file!='Total'");
    mSpCacheName = "cache"+System.nanoTime();

    Dataframe survey = MsReaders.surveyMzXML(mSummary.getStrCol("Raw file"), mRawfileRoot, mSpCacheName, null, 2);

    IOs.write(mResultRoot + "scan_survey.txt", survey.display("\t", "").toString());

    // readSpectrumIdentifier the tables of MS/MS scans
    mMsMs = readMsMsWithAnnotations();
    Dataframe offsets = Dataframe.merge(mMsMs, survey, true, true, "Raw file", "Scan number").setTitle("offsets");

    IOs.write(mResultRoot+"/composite_scans.txt", offsets.display("\t", "").toString());
  }

  public Dataframe readMsMsWithAnnotations() { return readMsMsWithAnnotations(mResultRoot); }
  public Dataframe readMsMsWithAnnotations(String root)
  {
    mResultRoot = root;

    Dataframe evidences = Dataframe.readtable(root+"/evidence.txt",  '\t', false).removeCols(MaxQuant.sRmdEvidence).setTitle("evidence");
    Dataframe msms      = Dataframe.readtable(root+"/msms.txt",      '\t', false).removeCols(MaxQuant.sRmdMsms).setTitle("msms");
    Dataframe scans     = Dataframe.readtable(root+"/msmsScans.txt", '\t', false).removeCols(MaxQuant.sRmdScan).setTitle("scan");

    // replace the "id" variable with "Evidence ID" so they can be joint
    evidences.renameCol("id", V_ID_EVI);
    evidences.renameCol("m/z", "m/z-calibrated");
    evidences.renameCol("Mass", "Mass-predicted");
    evidences.renameCol("Retention time", "RT-feature-calibrated");

     msms.renameCol("Retention time", "RT-MS2");
    scans.renameCol("Retention time", "RT-MS2");

    // force the ID columns to be categorical
    msms.init(     Var.VarType.CATEGORICAL, V_ID_EVI);
    evidences.init(Var.VarType.CATEGORICAL, V_ID_EVI);
    // join the msms table with the evidence table, which contains the annotated LCMS features by the Evidence ID
    Dataframe annotations = Dataframe.merge(msms, evidences, true, true, V_ID_EVI).setTitle("annot");
    // msms.txt contains only the annotated scans, msmsScans.txt hasAccession all the scans
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
      double xs[] = d.getDoubleCol("Retention time", true);
      double ys[] = d.getDoubleCol("Retention time calibration", true);
      double Xs[] = d.getDoubleCol("Calibrated retention time start", true);

      double[] Ys = MsStats.interpolate(xs, ys, rttol, Xs);
      d.setVar("interpolated", Ys);
      d.setVar("Calibrated RT", MsStats.sum(d.getDoubleCol("Retention time", true), d.getDoubleCol("Retention time calibration", true)));
    }
  }
}
