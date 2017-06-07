package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.algo.MsStats;
import org.ms2ms.algo.PSMs;
import org.ms2ms.data.Features;
import org.ms2ms.io.MsReaders;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.*;

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
  public static String[] sProteinInfo = {"Peptides","Razor + unique peptides","Unique peptides","Identification type","Sequence coverage","Intensity","iBAQ","LFQ intensity","MS/MS count"};

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
      d.setVar("Calibrated RT", MsStats.matrix_sum(d.getDoubleCol("Retention time", true), d.getDoubleCol("Retention time calibration", true)));
    }
  }
  public static PeptideMatch newPeptideMatchFromEvidence(Map<String, Object> info)
  {
    PeptideMatch M = new PeptideMatch(info.get("Sequence").toString());

    if (Stats.toDouble(info.get("Mass error [Da]")) !=null) M.setMassDiff(Stats.toDouble(info.get("Mass error [Da]")));
    if (Stats.toDouble(info.get("Mass"))            !=null) M.setNeutralPeptideMass(Stats.toDouble(info.get("Mass")));
    if (Stats.toInt(   info.get("Missed cleavages"))!=null) M.setNumMissedCleavages(Stats.toInt(info.get("Missed cleavages")));
    M.setRejected(info.get("Leading proteins").toString().indexOf("REV_")==0 || info.get("Leading proteins").toString().indexOf("CON_")==0);
    PSMs.addScore(M, info, "Score","Match q-value");

    // TODO add the mod
//    [1] "Sequence"                            "Length"                              "Modifications"
//      [4] "Modified.sequence"                   "Deamidation..NQ..Probabilities"      "Oxidation..M..Probabilities"
//      [7] "Deamidation..NQ..Score.Diffs"        "Oxidation..M..Score.Diffs"           "Acetyl..Protein.N.term."
//      [10] "Deamidation..NQ."                    "Gln..pyro.Glu"                       "Oxidation..M."
//      [13] "Missed.cleavages"                    "Proteins"                            "Leading.proteins"
//      [16] "Leading.razor.protein"               "Gene.names"                          "Protein.names"
//      [19] "Type"                                "Raw.file"                            "Experiment"
//      [22] "MS.MS.m.z"                           "Charge"                              "m.z"
//      [25] "Mass"                                "Resolution"                          "Uncalibrated...Calibrated.m.z..ppm."
//      [28] "Uncalibrated...Calibrated.m.z..Da."  "Mass.error..ppm."                    "Mass.error..Da."
//      [31] "Uncalibrated.mass.error..ppm."       "Uncalibrated.mass.error..Da."        "Max.intensity.m.z.0"
//      [34] "Retention.time"                      "Retention.length"                    "Calibrated.retention.time"
//      [37] "Calibrated.retention.time.start"     "Calibrated.retention.time.finish"    "Retention.time.calibration"
//      [40] "Match.time.difference"               "Match.m.z.difference"                "Match.q.value"
//      [43] "Match.score"                         "Number.of.data.points"               "Number.of.scans"
//      [46] "Number.of.isotopic.peaks"            "PIF"                                 "Fraction.of.total.spectrum"
//      [49] "Base.peak.fraction"                  "PEP"                                 "MS.MS.count"
//      [52] "MS.MS.scan.number"                   "Score"                               "Delta.score"
//      [55] "Combinatorics"                       "Intensity"                           "Reverse"
//      [58] "Potential.contaminant"               "id"                                  "Protein.group.IDs"
//      [61] "Peptide.ID"                          "Mod..peptide.ID"                     "MS.MS.IDs"
//      [64] "Best.MS.MS"                          "AIF.MS.MS.IDs"                       "Deamidation..NQ..site.IDs"
//      [67] "Oxidation..M..site.IDs"
    return M;
  }
  public static LcMsMsFeatures newMaxquant(String folder) throws IOException
  {
    LcMsMsFeatures lcmsms = new LcMsMsFeatures();

    Map<Long, Features> /*distincts = readMaxquant(new TabFile(folder+"/peptides.txt", TabFile.tabb), "id"),*/
        evidences = readMaxquant(new TabFile(folder+"/evidence.txt", TabFile.tabb), "id"),
        proteins = readMaxquant(new TabFile(folder+"/proteinGroups.txt", TabFile.tabb), "id");

    // pull the list of experiments
    Dataframe       summary = new Dataframe(folder+"/summary.txt", '\t');
    Set<String> experiments = new TreeSet<>(Strs.toStrings(summary.col("Experiment").values())),
                       runs = new TreeSet<>(Strs.toStrings(summary.col("Raw file").values()));

    // capture the protein, peptides and alignment
    Table<ProteinID, String, Map<String, Float>> protein_info = HashBasedTable.create();
    Table<PeptideFeature, String, Float>      peptide_expt_ai = HashBasedTable.create();
    Multimap<String, String>                         expt_run = HashMultimap.create();
    for (Long pid : proteins.keySet())
    {
      ProteinID protein = new ProteinID(pid); Features P = proteins.get(pid);
      protein.setName(Strs.split(P.get("Protein names").toString(), ';')[0]
      ).setGene(Strs.split(P.get("Gene names"   ).toString(), ';')[0]).setAccession(P.get("Majority protein IDs").toString());
      // gather the distinct peptides and their charge instances
      for (String expt : experiments)
        for (String info : sProteinInfo)
        {
          // fill-in the reported information about the protein in each experiment
          if (protein_info.get(protein, expt)==null) protein_info.put(protein, expt, Maps.newHashMap());
          protein_info.get(protein, expt).put(info, P.getFloat(info+" "+expt));
        }
      // now focus on the peptide and MS/MS
      Long[] ev_ids = Stats.toLongArray(P.get("Evidence IDs"), ';');
      if (ev_ids!=null)
        for (Long ev : ev_ids)
        {
          // deposit the peptide match
          PeptideMatch M = protein.put( newPeptideMatchFromEvidence(evidences.get(ev).getProperties()));

          Features F = evidences.get(ev);
          String run = F.get("Raw file").toString(), expt=F.get("Experiment").toString();
          Double  mz = Stats.toDouble(F.get("m/z")),
              rt = Stats.toDouble(F.get("Calibrated retention time")),
              ai = Stats.toDouble(F.get("Intensity"));

          // setup the entry
          PeptideFeature pf = protein.put(M, F.get("Modified sequence").toString(), Stats.toInt(F.get("Charge")));
          if (pf==null && expt==null && ai==null) peptide_expt_ai.put(pf.setProteinID(protein), expt, ai.floatValue());
          expt_run.put(expt, run);
        }
    }
    return lcmsms.setProteinInfo(protein_info).setPeptideExptAI(peptide_expt_ai);
  }
  public static Map<String, TreeBasedTable<Double, Double, Float>> newIntensityDistributions(String folder) throws IOException
  {
    System.out.println("Reading the features from "+folder);

    Map<String, TreeBasedTable<Double, Double, Float>> run_mz_rt_ai = new HashMap<>();

    int counts=0; TabFile features = new TabFile(folder+"/allPeptides.txt", TabFile.tabb);
    while (features.hasNext())
    {
      // Raw file        Type    Charge  m/z     Mass    Uncalibrated m/z        Resolution      Number of data points
      // Number of scans Number of isotopic peaks        PIF     Mass fractional part    Mass deficit    Mass precision [ppm]
      // Max intensity m/z 0     Retention time  Retention length        Retention length (FWHM) Min scan number
      // Max scan number Identified      MS/MS IDs       Sequence        Length  Modifications   Modified sequence
      // Proteins        Score   Intensity       Intensities     Isotope pattern MS/MS Count     MSMS Scan Numbers       MSMS Isotope Indices
      String run = features.get("Raw file");
      Double mz = features.getDouble("m/z"), rt = features.getDouble("Retention time"), ai = features.getDouble("Intensity");

      if (!run_mz_rt_ai.containsKey(run)) run_mz_rt_ai.put(run, TreeBasedTable.<Double, Double, Float>create());
      if ( run_mz_rt_ai.get(run).get(mz,rt)==null) run_mz_rt_ai.get(run).put(mz,rt, 0f);
      run_mz_rt_ai.get(run).put(mz,rt, ai.floatValue()+run_mz_rt_ai.get(run).get(mz,rt));

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();

    return run_mz_rt_ai;
  }
  public static Map<Long, Features> readMaxquant(TabFile features, String id) throws IOException
  {
    System.out.println("Reading the results from "+features.getFileName() + " anchored to " + id);
    int counts=0; Map<Long, Features> results = new HashMap<>();
    while (features.hasNext())
    {
      if (features.get(id)!=null)
        results.put(features.getLong(id), new Features().setProperties(features.getMappedRow()));

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();

    return results;
  }
  public static void fillBaseline(String folder, Tolerance RT, Tolerance M) throws IOException
  {
    Map<String, TreeBasedTable<Double, Double, Float>> dist = newIntensityDistributions(folder);
    LcMsMsFeatures lcmsms = newMaxquant(folder);

    lcmsms.toPeptideExptMatrix().csv("/Users/yuw/Downloads/before.csv", 3, "\t");

    // estimate the baseline
    for (PeptideFeature P : lcmsms.getPeptideExptAI().rowKeySet())
      for (String expt : lcmsms.getPeptideExptAI().row(P).keySet())
        if (lcmsms.getPeptideExptAI().get(P, expt)==null)
        {
          // come up with an estimate
          Float sum=0f;
          for (String run : lcmsms.getExptRun().get(expt))
            sum+=estBaseline(dist.get(run), Tools.toBound(M, P.getNeutralPeptideMass()/P.getCharge()), Tools.toBound(RT, P.getRT()));

          lcmsms.getPeptideExptAI().put(P, expt, sum);
        }

    lcmsms.toPeptideExptMatrix().csv("/Users/yuw/Downloads/after.csv", 3, "\t");
  }
  public static Float estBaseline(TreeBasedTable<Double, Double, Float> mz_rt_ai, Range<Double> mz, Range<Double> rt)
  {
    Map<Double, Map<Double, Float>> mz0 = mz_rt_ai.rowMap().subMap(mz.lowerEndpoint(), mz.upperEndpoint());
    Float baseline=0f;
    if (mz0!=null)
      for (Map<Double, Float> R : mz0.values())
        for (Double r : R.keySet())
          if (rt.contains(r) && (baseline==0 || baseline>R.get(r))) baseline=R.get(r);

    return baseline;
  }
  public static void readProteins(TabFile features) throws IOException
  {
    System.out.println("Reading the protein groups from "+features.getFileName());
    int counts=0;
    while (features.hasNext())
    {
      // Protein IDs     Majority protein IDs    Peptide counts (all)    Peptide counts (razor+unique)   Peptide counts (unique)
      // Protein names   Gene names      Fasta headers   Number of proteins      Peptides        Razor + unique peptides  Unique peptides
      // Peptides A1     Peptides A2     Peptides A3     Peptides A4     Peptides A5     Peptides N1     Peptides N2     Peptides N3     Peptides N4     Peptides N5
      // Razor + unique peptides A1      Razor + unique peptides A2      Razor + unique peptides A3      Razor + unique peptides A4      Razor + unique peptides A5
      // Razor + unique peptides N1      Razor + unique peptides N2      Razor + unique peptides N3      Razor + unique peptides N4      Razor + unique peptides N5
      // Unique peptides A1      Unique peptides A2      Unique peptides A3      Unique peptides A4      Unique peptides A5
      // Unique peptides N1      Unique peptides N2      Unique peptides N3      Unique peptides N4      Unique peptides N5
      // Sequence coverage [%]   Unique + razor sequence coverage [%]    Unique sequence coverage [%]    Mol. weight [kDa]
      // Sequence length Sequence lengths        Q-value Score   Identification type A1  Identification type A2  Identification type A3  Identification type A4  Identification type A5
      // Identification type N1  Identification type N2  Identification type N3  Identification type N4  Identification type N5
      // Sequence coverage A1 [%]        Sequence coverage A2 [%]        Sequence coverage A3 [%]        Sequence coverage A4 [%]        Sequence coverage A5 [%]
      // Sequence coverage N1 [%]        Sequence coverage N2 [%]        Sequence coverage N3 [%]        Sequence coverage N4 [%]        Sequence coverage N5 [%]
      // Intensity       Intensity A1    Intensity A2    Intensity A3    Intensity A4    Intensity A5    Intensity N1    Intensity N2    Intensity N3    Intensity N4    Intensity N5
      // iBAQ    iBAQ A1 iBAQ A2 iBAQ A3 iBAQ A4 iBAQ A5 iBAQ N1 iBAQ N2 iBAQ N3 iBAQ N4 iBAQ N5
      // LFQ intensity A1        LFQ intensity A2        LFQ intensity A3        LFQ intensity A4        LFQ intensity A5
      // LFQ intensity N1        LFQ intensity N2        LFQ intensity N3        LFQ intensity N4        LFQ intensity N5
      // MS/MS count A1  MS/MS count A2  MS/MS count A3  MS/MS count A4  MS/MS count A5  MS/MS count N1  MS/MS count N2  MS/MS count N3  MS/MS count N4  MS/MS count N5  MS/MS count
      // Only identified by site Reverse Potential contaminant   id      Peptide IDs     Peptide is razor
      // Mod. peptide IDs        Evidence IDs    MS/MS IDs       Best MS/MS      Deamidation (NQ) site IDs
      // Oxidation (M) site IDs  Deamidation (NQ) site positions Oxidation (M) site positions
      Double mz = features.getDouble("m/z"), rt = features.getDouble("Retention time");
      // add the features

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();
  }
}
