package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.ms2ms.data.Binary;
import org.ms2ms.data.Features;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.util.*;
import java.util.regex.Pattern;

/** starting with LCMS features from Dinosauer, combine with the MS2 scans from Maxquant, look for the utilization, etc
 *
 * Created by Wen Yu on 4/7/17.
 */
public class LcMsMsFeatures implements Binary
{
  public static final String COL_RUN        = "run";
  public static final String COL_MZ         = "mz";
  public static final String COL_Z          = "charge";
  public static final String COL_RT         = "rt";
  public static final String COL_RT_BOUND   = "rt range";
  public static final String COL_MASS       = "mass";
  public static final String COL_MASS_CAL   = "calc mass";
  public static final String COL_MASS_OBS   = "obs mass";
  public static final String COL_PPM_PREC   = "precision in ppm";
  public static final String COL_PPM        = "ppm";
  public static final String COL_SCAN       = "Scan number";
  public static final String COL_SCAN_BOUND = "Scan range";
  public static final String COL_AI_APEX    = "intensityApex";
  public static final String COL_AI_AREA    = "intensitySum";
  public static final String COL_SEQ        = "sequence";
  public static final String COL_PEPTIDE    = "peptide";
  public static final String COL_PROT       = "protein";
  public static final String COL_MOD        = "mods";
  public static final String COL_PSMID      = "psmID";
  public static final String COL_FDR        = "fdr";
  public static final String COL_SCORE      = "score";
  public static final String COL_PROPERTY   = "Properties";

  static
  {
    Features.setOrderBy(COL_RUN, COL_RT, COL_MZ);
  }

  // RT, m
  private TreeBasedTable<Integer, String, Double>   mScan_RT;
  private MultiTreeTable<Integer, String, Features> mScan_MS2;
  private Map<Long, Features>                       mPSMs;
  private TreeBasedTable<Double, Double, Features>  mIons, mMsMs; // mz,rt, feature
  private TreeMultimap<Double, Features>            mPeptides;

  private Table<ProteinID, String, Map<String, Float>> mProteinInfo;
  private Table<PeptideFeature, String, Float>      mPeptideExptAI;;
  private Multimap<String, String>                  mExptRun;

  private Collection<Pattern> mRunFilters;

  // assemble the protein and peptide matches
  private Map<String, ProteinID> mProteinIDs;
  private Table<PeptideMatch, String, Features> mPeptideRunFeatures;
  private Table<String, String, Features> mSampleProperties;

  public LcMsMsFeatures() { super(); }
  public LcMsMsFeatures(TabFile lcms, String mq, OffsetPpmTolerance mz, String... runs)
  {
    init(mq, lcms, mz, runs);
  }
  public LcMsMsFeatures(String mq, OffsetPpmTolerance mz, String... runs)
  {
    init(mq, null, mz, runs);
  }

  private void init(String mq, TabFile lcms, OffsetPpmTolerance mz, String... runs)
  {
    try
    {
      initRunFilters(runs);

      readMaxquantMS1(     new TabFile(mq+"/msScans.txt", "\t"));
      readMaxquantMsMs(    new TabFile(mq+"/msmsScans.txt", "\t"));
      if (IOs.exists(mq+"/allPeptides.txt")) readMaxquantFeatures(new TabFile(mq+"/allPeptides.txt", "\t"));
      readMaxquantEvidence(new TabFile(mq+"/evidence.txt", "\t"));

      Peptides2Ions();

      readDinosauer(lcms);
    }
    catch (IOException ie) { ie.printStackTrace(); }
  }

  private void initRunFilters(String... runs)
  {
    if (Tools.isSet(runs))
    {
      mRunFilters = new ArrayList<>();
      for (String run : runs)
        mRunFilters.add(Pattern.compile(run));
    }
  }
  public Table<ProteinID, String, Map<String, Float>> getProteinInfo()   { return mProteinInfo; }
  public Table<PeptideFeature, String, Float>         getPeptideExptAI() { return mPeptideExptAI; }
  public Multimap<String, String>                     getExptRun()       { return mExptRun; }
  public Collection<String>                           getRuns4Expt(String s) { return mExptRun!=null?mExptRun.get(s):null; }
  public RowSortedTable<Double, Double, Features>     getIons()          { return mIons; }

  public LcMsMsFeatures setProteinInfo(  Table<ProteinID, String, Map<String, Float>> s) { mProteinInfo  =s; return this; }
  public LcMsMsFeatures setPeptideExptAI(Table<PeptideFeature, String, Float>         s) { mPeptideExptAI=s; return this; }
  public LcMsMsFeatures setExptRuns(     Multimap<String, String>                     s) { mExptRun      =s; return this; }

  private boolean okRun(String s)
  {
    if (!Tools.isSet(mRunFilters)) return true;
    for (Pattern filter : mRunFilters)
      if (filter.matcher(s).matches()) return true;

    return false;
  }
  public void readMaxquantMsMs(TabFile features) throws IOException
  {
    if (mMsMs==null)  mMsMs = TreeBasedTable.create();
    int counts=0; mScan_MS2 = MultiTreeTable.create();

    System.out.println("Reading the MS/MS features from "+features.getFileName());
    while (features.hasNext())
    {
      if (!okRun(features.get("Raw file"))) continue;

      // Raw file        Scan number     Retention time  Ion injection time      Total ion current       Collision energy        Summations      Base peak intensity     Elapsed time
      // Identified      MS/MS IDs       Sequence        Length  Filtered peaks  m/z     Mass    Charge  Type    Fragmentation   Mass analyzer   Parent intensity fraction
      // Fraction of total spectrum      Base peak fraction      Precursor full scan number      Precursor intensity     Precursor apex fraction Precursor apex offset
      // Precursor apex offset time      Scan event number       Modifications   Modified sequence       Proteins        Score   Intens Comp Factor      CTCD Comp       RawOvFtT
      // AGC Fill        Scan index      MS scan index   MS scan number

      Features F = addMaxquantPSM(addMaxquantCore(new Features(), features), features);

      mMsMs.put(    F.getDouble(COL_MZ),   F.getDouble(COL_RT), F);
      mScan_MS2.put(F.getInt(   COL_SCAN), F.getStr(  COL_RUN), F);

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();
  }
  private Features addMaxquantCore(Features F, TabFile features)
  {
    // add the features
    F.add(COL_MZ,       features.getDouble("m/z"));
    F.add(COL_RT,       features.getDouble("Retention time"));
    F.add(COL_Z,        features.getInt("Charge"));
    F.add(COL_MASS,     features.getDouble("Mass"));
    F.add(COL_PPM_PREC, features.getDouble("Mass precision [ppm]"));

    F.add(COL_SCAN,     features.getInt("Scan number"));
    F.add(COL_Z,        features.getInt("Charge"));
    F.add(COL_MASS,     features.getDouble("Mass"));
    F.add(COL_AI_APEX,  features.getDouble("Precursor intensity"));
    F.add(COL_RUN,      features.get("Raw file"));

    if (features.getInt("Min scan number")!=null && features.getInt("Max scan number")!=null)
      F.add(COL_RT_BOUND, Range.closed(
        mScan_RT.get(features.getInt("Min scan number"), features.get("Raw file")),
        mScan_RT.get(features.getInt("Max scan number"), features.get("Raw file"))));

    if (features.getInt("Min scan number")!=null && features.getInt("Max scan number")!=null)
      F.add(COL_SCAN_BOUND, Range.closed(features.getInt("Min scan number"), features.getInt("Max scan number")));

    return F;
  }
  private Features addMaxquantPSM(Features F, TabFile features)
  {
    F.add(COL_PSMID,   features.get("MS/MS IDs"));
    if (Strs.isSet(features.get("Sequence")))
    {
      F.add(COL_SEQ,     features.get("Sequence"));
      F.add(COL_PEPTIDE, features.get("Modified sequence"));
      F.add(COL_PROT,    features.get("Proteins"));
      if (!"Unmodified".equals(features.get("Modifications"))) F.add(COL_MOD, features.get("Modifications"));

      F.add(COL_FDR,     features.getDouble("PEP"));
      F.add(COL_SCORE,   features.getDouble("Score"));
    }

    return F;
  }
  public void readMaxquantEvidence(TabFile features) throws IOException
  {
    if (mPSMs==null) mPSMs = new TreeMap<>();

    int counts=0;

    System.out.println("Reading the evidence features from "+features.getFileName());
    while (features.hasNext())
    {
      if (!okRun(features.get("Raw file"))) continue;

      // Sequence        Length  Modifications   Modified sequence       Deamidation (NQ) Probabilities  Oxidation (M) Probabilities     Deamidation (NQ) Score Diffs
      // Oxidation (M) Score Diffs       Acetyl (N-term) Deamidation (NQ)        Oxidation (M)   Missed cleavages        Proteins        Leading proteins
      // Leading razor protein   Gene names      Protein names   Type    Raw file        MS/MS m/z       Charge  m/z     Mass    Resolution      Uncalibrated - Calibrated m/z [ppm]
      // Uncalibrated - Calibrated m/z [Da]      Mass error [ppm]        Mass error [Da] Uncalibrated mass error [ppm]   Uncalibrated mass error [Da]    Max intensity m/z 0
      // Retention time  Retention length        Calibrated retention time       Calibrated retention time start Calibrated retention time finish        Retention time calibration
      // Match time difference  Match m/z difference    Match q-value   Match score     Number of data points   Number of scans Number of isotopic peaks        PIF
      // Fraction of total spectrum      Base peak fraction      PEP     MS/MS count     MS/MS scan number       Score   Delta score     Combinatorics   Intensity
      // Reverse Potential contaminant   id      Protein group IDs       Peptide ID      Mod. peptide ID MS/MS IDs       Best MS/MS      AIF MS/MS IDs   Deamidation (NQ) site IDs
      // Oxidation (M) site IDs

      // add the features
      Features F = addMaxquantPSM(addMaxquantCore(new Features(), features), features);
      F.add(COL_PPM,  features.getDouble("Mass error [ppm]"));

      if (F.get(COL_PSMID)!=null)
        for (String id : Strs.split(F.getStr(COL_PSMID), ';'))
          if (!"".equals(id)) mPSMs = Tools.put(mPSMs, Long.parseLong(id), F);

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();
  }

  public void readMaxquantMS1(TabFile features) throws IOException
  {
    System.out.println("Reading the MS1 scans from "+features.getFileName()); int counts=0;

    mScan_RT = TreeBasedTable.create();
    while (features.hasNext())
    {
      if (!okRun(features.get("Raw file"))) continue;

      // Raw file        Scan number     Scan index      Retention time  Cycle time      Ion injection time
      // Base peak intensity     Total ion current       MS/MS count     Mass calibration        Peak length
      // Isotope pattern length  Multiplet length        Peaks / s       Single peaks / s        Isotope patterns / s
      // Single isotope patterns / s     Multiplets / s  Identified multiplets / s       Multiplet identification rate [%]
      // MS/MS / s       Identified MS/MS / s    MS/MS identification rate [%]   Intens Comp Factor
      // CTCD Comp       RawOvFtT        AGC Fill
      mScan_RT.put(features.getInt("Scan number"), features.get("Raw file"), features.getDouble("Retention time"));

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();
  }
  public void readMaxquantFeatures(TabFile features) throws IOException
  {
    // initiate the collections
    if (mIons==null) mIons = TreeBasedTable.create();

    System.out.println("Reading the features from "+features.getFileName()); int counts=0;
    while (features.hasNext())
    {
      if (!okRun(features.get("Raw file"))) continue;

      // Raw file        Type    Charge  m/z     Mass    Uncalibrated m/z        Resolution      Number of data points
      // Number of scans Number of isotopic peaks        PIF     Mass fractional part    Mass deficit    Mass precision [ppm]
      // Max intensity m/z 0     Retention time  Retention length        Retention length (FWHM) Min scan number
      // Max scan number Identified      MS/MS IDs       Sequence        Length  Modifications   Modified sequence
      // Proteins        Score   Intensity       Intensities     Isotope pattern MS/MS Count     MSMS Scan Numbers       MSMS Isotope Indices
      Features F = addMaxquantPSM(addMaxquantCore(new Features(), features), features);
      F.add(COL_PSMID, Long.parseLong(F.getStr(COL_PSMID)));
      F.add(COL_AI_APEX,  features.getDouble("Intensity"));

      mIons.put(F.getDouble(COL_RT), F.getDouble(COL_MZ), F);

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts); features.close();
  }

  public LcMsMsFeatures readDinosauer(String features)
  {
    try
    {
      readDinosauer(new TabFile(features, TabFile.tabb));
    }
    catch (IOException ie)
    {
      ie.printStackTrace();
    }
    return this;
  }
  public void readDinosauer(TabFile features) throws IOException
  {
    if (features==null) return;

    System.out.println("Reading the features from "+features.getFileName()); int counts=0;
    while (features.hasNext())
    {
      // mz      mostAbundantMz  charge  rtStart rtApex  rtEnd   fwhm    nIsotopes       nScans  averagineCorr   mass    massCalib       intensityApex   intensitySum
      //Map<String, String> row  = features.nextRow();

      Double mz = features.getDouble("mz"), rt = features.getDouble("rtApex");
      // add the features
      addFeature(rt, mz, COL_MZ,       mz);
      addFeature(rt, mz, COL_RT,       rt);
      addFeature(rt, mz, COL_Z,        features.getDouble("charge"));
      addFeature(rt, mz, COL_MASS,     features.getDouble("mass"));
      addFeature(rt, mz, COL_AI_APEX,  features.getDouble("intensityApex"));
      addFeature(rt, mz, COL_RT_BOUND, Range.closed(features.getDouble("rtStart"), features.getDouble("rtEnd")));
      addFeature(rt, mz, COL_PROPERTY, features.getMappedRow());

      if (++counts%10000==0) System.out.print(".");
      if (counts%1000000==0) System.out.println();
    }
    System.out.println(counts);
  }
  private TreeMultimap<Double, Features> Ions2Peptides(OffsetPpmTolerance mz)
  {
    System.out.println("Rolling up the isotope clusters...");

    // arrange the ions by their intensities
    TreeMultimap<Double, Features>       ai_pts = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
    TreeMultimap<Double, Features> mass_feature = TreeMultimap.create(), rt_features = TreeMultimap.create();
    for (Table.Cell<Double, Double, Features> cell : mIons.cellSet())
    {
      ai_pts.put(cell.getValue().getDouble(COL_AI_APEX), cell.getValue());
      mass_feature.put(cell.getValue().getDouble(COL_MASS), cell.getValue());
    }
    // roll-up the charge clusters
    mPeptides = TreeMultimap.create();
    for (Double ai : ai_pts.keySet())
      for (Features pk : ai_pts.get(ai))
      {
        if (!pk.isValid()) continue;
        Range<Double> rts  = (Range<Double>) pk.get(COL_RT_BOUND);
        Double mass        = pk.getDouble(COL_MASS);
        Double mass_unique = Tools.unique(mPeptides.keySet(), mass);

        Map<Double, Collection<Features>> slice = mass_feature.asMap().subMap(mz.getMin(mass), mz.getMax(mass));

        if (Tools.isSet(slice))
        {
          rt_features.clear();
          for (Collection<Features> ss : slice.values())
            for (Features F : ss)
              if (F.getDouble(COL_RT)!=null)
                rt_features.put(F.getDouble(COL_RT), F);

          Map<Double, Collection<Features>> sect = rt_features.asMap().subMap(rts.lowerEndpoint()-1, rts.upperEndpoint()+1);

          if (sect!=null)
            for (Collection<Features> ss : sect.values())
              for (Features F : ss)
              {
                // need to make sure the RT ranges are compatible
                if (!Tools.isA4B((Range<Double>) F.get(COL_RT_BOUND), rts, 75d)) continue;
                mPeptides.put(mass_unique, F.clone()); F.invalidate();
              }
        }
      }

    ai_pts.clear(); ai_pts=null; mass_feature.clear(); mass_feature=null; rt_features.clear(); rt_features=null;

    System.out.println(mPeptides.keySet().size()+" Peptide clusters based on " + mIons.size() + " charge clusters");
    return mPeptides;
  }
  public void Peptides2Ions()
  {
    if (Tools.isSet(mIons) && Tools.isSet(mPSMs))
      for (Features F : mIons.values())
      {
        Long id = F.getLong(COL_PSMID);
        if (id!=null && mPSMs.containsKey(id))
        {
          F.add(COL_MASS_OBS,         F.getDouble(COL_MASS));
          F.add(COL_MASS, mPSMs.get(id).getDouble(COL_MASS));
          F.add(COL_PPM,  mPSMs.get(id).getDouble(COL_PPM));
        }
      }
  }

  public void Peptides2MsMs(OffsetPpmTolerance mz)
  {
    int ms1=0, ms2=0, pep=0;
    for (Double mass : mPeptides.keySet())
    {
      ms1++; int found = 0;
      for (Features F : mPeptides.get(mass))
      {
        // do we have any MS/MS matching any one of the features?
        List<Features> msms = hasMsMs(F, mz);
        if (msms!=null) found += msms.size();
      }
      if (found>0) { pep++; ms2+=found; }
    }
    System.out.println("MS2/Peptides/MS1: " + ms2+"/"+pep+"/"+ms1);
  }
  // do we have a corresponding MsMs for the MS1 feature?
  private List<Features> hasMsMs(Features lcms, OffsetPpmTolerance mz)
  {
    Range<Double> m = mz.getBoundary(lcms.getDouble(COL_MZ));
    SortedMap<Double, Map<Double, Features>> slice  = mMsMs.rowMap().subMap(m.lowerEndpoint(), m.upperEndpoint());
    // any MS2?
    if (Tools.isSet(slice))
    {
      Range<Double> rt = (Range<Double> )lcms.get(COL_RT_BOUND);

      List<Features> found = new ArrayList<>();
      for (Map<Double, Features> fs : slice.values())
        for (Features ff : fs.values())
          if (rt.contains(ff.getDouble(COL_RT))) { found.add(ff); ff.add("LCMS", lcms); }

      return found;
    }
    return null;
  }
  private LcMsMsFeatures addFeature(Double rt, Double mz, String key, Object val)
  {
    if (mIons==null) mIons = TreeBasedTable.create();
    // add the feature if not present
    Features F = mIons.get(rt, mz);
    if (F==null)
    {
      F = new Features(); mIons.put(rt, mz, F);
    }
    F.add(key, val);

    return this;
  }
//  // remember m/z first
//  private LcMsMsFeatures addMsMs(Double rt, Double mz, String key, Object val)
//  {
//    if (mMsMs==null) mMsMs = TreeBasedTable.create();
//    // add the feature if not present
//    Features F = mMsMs.get(mz, rt);
//    if (F==null)
//    {
//      F = new Features(); mMsMs.put(mz, rt, F);
//    }
//    F.add(key, val);
//
//    return this;
//  }

  public TreeMultimap<Double, Features> getFeaturesByIntensity()
  {
    if (!Tools.isSet(mIons)) return null;

    TreeMultimap<Double, Features> F = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
    for (Features ff : mIons.values())
      F.put(ff.getDouble("intensityApex"), ff);

    return F;
  }
  public LcMsMsFeatures addProteinInfo(ProteinID protein, String expt, String info, Float val)
  {
    if (mProteinInfo==null) mProteinInfo=HashBasedTable.create();
    if (mProteinInfo.get(protein, expt)==null) mProteinInfo.put(protein, expt, Maps.newHashMap());
    if (val      !=null) mProteinInfo.get(protein, expt).put(info, val);

    return this;
  }
  public LcMsMsFeatures addPeptideExptIntensity(PeptideFeature pf, String expt, Double ai)
  {
    if (mPeptideExptAI==null) mPeptideExptAI=HashBasedTable.create();
    if (pf!=null && expt!=null && ai!=null)
      mPeptideExptAI.put(pf, expt, Stats.sumFloats(mPeptideExptAI.get(pf, expt), ai.floatValue()));

    return this;
  }
  public LcMsMsFeatures addExptRun(String expt, String run)
  {
    if (mExptRun==null) mExptRun=HashMultimap.create();
    if (expt!=null && run!=null) mExptRun.put(expt, run);

    return this;
  }

  public Dataframe toPeptideExptMatrix()
  {
    Dataframe df = new Dataframe("peptide-charge vs experiment matrix");
    // going thro each peptide features
    for (PeptideFeature F : getPeptideExptAI().rowKeySet())
    {
      String row = F.getTitle()+"z"+F.getCharge()+","+Tools.d2s(F.getRT(), 3)+"sec";
      // put in the peptide cols
      df.put(row, "Peptide",   F.getTitle());
      df.put(row, "Sequence",  F.toSymbolString());
      df.put(row, "PSMs",      F.getMatches().size());
      df.put(row, "m/z",       F.getMz());
      df.put(row, "z",         F.getCharge());
      df.put(row, "Protein",   F.getProteinID().getName());
      df.put(row, "Accession", F.getProteinID().getAccession());
      df.put(row, "Gene",      F.getProteinID().getGene());

      for (String expt : getPeptideExptAI().row(F).keySet())
        df.put(row, expt, getPeptideExptAI().get(F, expt));
    }

    return df.init(false);
  }
  public Dataframe toProteinExptMatrix()
  {
    // let's build the protein-peptide groups
    Multimap<ProteinID, PeptideFeature> protein_peptide = HashMultimap.create();
    for (PeptideFeature F : getPeptideExptAI().rowKeySet()) protein_peptide.put(F.getProteinID(), F);

    Dataframe df = new Dataframe("protein vs experiment matrix");
    // going thro each peptide features
    for (ProteinID pid : protein_peptide.keySet())
    {
      String row = pid.toString();
      df.put(row, "Protein",   pid.getName());
      df.put(row, "Accession", pid.getAccession());
      df.put(row, "Gene", pid.getGene());
      df.put(row, "Peptides", protein_peptide.get(pid).size());

      int psm=0;
      for (PeptideFeature F : protein_peptide.get(pid)) psm+=F.getMatches().size();
      df.put(row, "PSMs",     psm);

      for (String expt : getExptRun().keySet())
      {
        double ai=0;
        for (PeptideFeature F : protein_peptide.get(pid))
          if (getPeptideExptAI().get(F, expt)!=null) ai+=getPeptideExptAI().get(F, expt);

        df.put(row, expt, ai);
      }
    }

    return df.init(false);
  }

  public Map<String, TreeBasedTable<Double, Double, Features>> indexing()
  {
    Map<String, TreeBasedTable<Double, Double, Features>>  RunRtAiIons = new TreeMap<>();
    for (Features F : mIons.values())
    {
      TreeBasedTable<Double, Double, Features> RtAiIons = RunRtAiIons.get(F.get(COL_RUN));
      if (RtAiIons==null)
      {
        RtAiIons = TreeBasedTable.create(); RunRtAiIons.put((String )F.get(COL_RUN), RtAiIons);
      }
      if (F.get(COL_PEPTIDE)!=null) RtAiIons.put(F.getDouble(COL_RT), F.getDouble(COL_AI_APEX), F);
    }

    return RunRtAiIons;
  }
  public SortedMap<Double, Features> peering(TreeBasedTable<Double, Double, Features>  RtAiIons, Features F,
                                              double dRT, SortedMap<Double, Features> AiIons)
  {
    if (AiIons==null) AiIons = new TreeMap<>(Ordering.natural().reversed()); else AiIons.clear();
    // gather the peers
    SortedMap<Double, Map<Double, Features>> slice = RtAiIons.rowMap().subMap(F.getDouble(COL_RT)-dRT, F.getDouble(COL_RT)+dRT);
    if (slice!=null)
      for (Map<Double, Features> s : slice.values())
        if (!s.containsValue(F)) AiIons.putAll(s);

    return AiIons;
  }
  public Multimap<String, Features> peers(TreeBasedTable<Double, Double, Features>  RtAiIons, Features F,
                                          double dRT, Multimap<String, Features> PepIons,
                                          Multimap<Double, Features>  ppmIons)
  {
    if (PepIons==null) PepIons = TreeMultimap.create(); else PepIons.clear();
    if (ppmIons==null) ppmIons = TreeMultimap.create(); else ppmIons.clear();
    // gather the peers
    SortedMap<Double, Map<Double, Features>> slice = RtAiIons.rowMap().subMap(F.getDouble(COL_RT)-dRT, F.getDouble(COL_RT)+dRT);
    if (slice!=null)
      for (Map<Double, Features> s : slice.values())
        for (Features f : s.values())
          if (!s.containsValue(F)/* && f.get(LcMsMsFeatures.COL_PEPTIDE)!=null*/)
          {
            PepIons.put(f.getStr(LcMsMsFeatures.COL_PEPTIDE), f);
            Tools.put(ppmIons, f.getDouble(LcMsMsFeatures.COL_PPM),  f);
          }

    return PepIons;
  }
  public SortedMap<Double, Features> peer2(TreeBasedTable<Double, Double, Features>  RtAiIons, Features F,
                                             double dRT, SortedMap<Double, Features> AiIons)
  {
    if (AiIons==null) AiIons = new TreeMap<>(Ordering.natural().reversed()); else AiIons.clear();
    // gather the peers
    SortedMap<Double, Map<Double, Features>> slice = RtAiIons.rowMap().subMap(F.getDouble(COL_RT)-dRT, F.getDouble(COL_RT)+dRT);
    if (slice!=null)
      for (Map<Double, Features> s : slice.values())
        if (!s.containsValue(F)) AiIons.putAll(s);

    return AiIons;
  }

  public static Double mean(Collection<Features> fs, String col)
  {
    if (!Tools.isSet(fs)) return null;

    Collection<Double> nums = new ArrayList<>();
    for (Features f : fs)
      Tools.add(nums, f.getDouble(col));

    Double avg = Stats.mean(nums); nums=Tools.dispose(nums);
    return avg;
  }
  @Override
  public void write(DataOutput ds) throws IOException
  {
//    IOs.writeIntStrDouble(   ds, mScan_RT);
//    MsIO.writeIntStrFeatures(ds, mScan_MS2);
//    IOs.writeLongMap(ds, mDistinctPeptides);
//    IOs.writeDoubleDoubleBin(ds, mIons);
//    IOs.writeDoubleDoubleBin(ds, mMsMs); // mz,rt, feature
//    IOs.writeDoubleMultimap(ds, mPeptides);
//
//    private Table<ProteinID, String, Map<String, Float>> mProteinInfo;
//    private Table<PeptideFeature, String, Float>      mPeptideExptAI;;
//    private Multimap<String, String>                  mExptRun;
//
//    // assemble the protein and peptide matches
//    private Map<String, ProteinID> mProteinIDs;
//    private Table<PeptideMatch, String, Features> mPeptideRunFeatures;
//    private Table<String, String, Features> mSampleProperties;
  }

  @Override
  public void read(DataInput ds) throws IOException
  {

  }
}
