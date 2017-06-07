package org.ms2ms.data.ms;

import com.compomics.util.experiment.identification.matches.PeptideMatch;
import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.Features;
import org.ms2ms.data.Point;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.util.*;

/** starting with LCMS features from Dinosauer, combine with the MS2 scans from Maxquant, look for the untilization, etc
 *
 * Created by Wen Yu on 4/7/17.
 */
public class LcMsMsFeatures
{
  public static final String COL_MZ         = "mz";
  public static final String COL_Z          = "charge";
  public static final String COL_RT         = "rt";
  public static final String COL_RT_BOUND   = "rt range";
  public static final String COL_MASS       = "mass";
  public static final String COL_SCAN       = "Scan number";
  public static final String COL_SCAN_BOUND = "Scan range";
  public static final String COL_AI_APEX    = "intensityApex";
  public static final String COL_AI_AREA    = "intensitySum";
  public static final String COL_PROPERTY   = "Properties";

  // RT, m
  private TreeBasedTable<Integer, String, Double>   mScan_RT;
  private MultiTreeTable<Integer, String, Features> mScan_MS2;
  private Map<Long, Features>                       mDistinctPeptides;
  private TreeBasedTable<Double, Double, Features>  mIons, mMsMs; // mz,rt, feature
  private TreeMultimap<Double, Features>            mPeptides;

  private Table<ProteinID, String, Map<String, Float>> mProteinInfo;
  private Table<PeptideFeature, String, Float>      mPeptideExptAI;;
  private Multimap<String, String>                  mExptRun;

  // assemble the protein and peptide matches
  private Map<String, ProteinID> mProteinIDs;
  private Table<PeptideMatch, String, Features> mPeptideRunFeatures;
  private Table<String, String, Features> mSampleProperties;

  public LcMsMsFeatures() { super(); }
  public LcMsMsFeatures(TabFile lcms, String mq, OffsetPpmTolerance mz)
  {
    try
    {
      readMaxquantMS1(     new TabFile(mq+"/msScans.txt", "\t"));
      readMaxquantMsMs(new TabFile(mq+"/msmsScans.txt", "\t"));
      readMaxquantFeatures(new TabFile(mq+"/allPeptides.txt", "\t"));

      readDinosauer(lcms);

      System.out.println(mIons.size());

//      Ions2Peptides(mz);
//      Peptides2MsMs(mz);
    }
    catch (IOException ie)
    {

    }
  }
  public LcMsMsFeatures(String mq, OffsetPpmTolerance mz)
  {
    try
    {
      readMaxquantMS1(     new TabFile(mq+"/msScans.txt", "\t"));
      readMaxquantMsMs(    new TabFile(mq+"/msmsScans.txt", "\t"));
      readMaxquantFeatures(new TabFile(mq+"/allPeptides.txt", "\t"));

      System.out.println(mIons.size());

//      Ions2Peptides(mz);
//      Peptides2MsMs(mz);
    }
    catch (IOException ie)
    {

    }
  }

  private void init()
  {

  }

  public Table<ProteinID, String, Map<String, Float>> getProteinInfo()   { return mProteinInfo; }
  public Table<PeptideFeature, String, Float>         getPeptideExptAI() { return mPeptideExptAI; }
  public Multimap<String, String>                     getExptRun()       { return mExptRun; }
  public Collection<String>                           getRuns4Expt(String s) { return mExptRun!=null?mExptRun.get(s):null; }

  public LcMsMsFeatures setProteinInfo(  Table<ProteinID, String, Map<String, Float>> s) { mProteinInfo  =s; return this; }
  public LcMsMsFeatures setPeptideExptAI(Table<PeptideFeature, String, Float>         s) { mPeptideExptAI=s; return this; }
  public LcMsMsFeatures setExptRuns(     Multimap<String, String>                     s) { mExptRun      =s; return this; }

  public void readMaxquantMsMs(TabFile features) throws IOException
  {
    System.out.println("Reading the features from "+features.getFileName());
    int counts=0; mScan_MS2 = MultiTreeTable.create();
    while (features.hasNext())
    {
      // Raw file        Scan number     Retention time  Ion injection time      Total ion current       Collision energy        Summations      Base peak intensity     Elapsed time
      // Identified      MS/MS IDs       Sequence        Length  Filtered peaks  m/z     Mass    Charge  Type    Fragmentation   Mass analyzer   Parent intensity fraction
      // Fraction of total spectrum      Base peak fraction      Precursor full scan number      Precursor intensity     Precursor apex fraction Precursor apex offset
      // Precursor apex offset time      Scan event number       Modifications   Modified sequence       Proteins        Score   Intens Comp Factor      CTCD Comp       RawOvFtT
      // AGC Fill        Scan index      MS scan index   MS scan number
      //Map<String, String> row  = features.nextRow();

      Double mz = features.getDouble("m/z"), rt = features.getDouble("Retention time");
      // add the features
      addMsMs(rt, mz, COL_MZ,       mz);
      addMsMs(rt, mz, COL_RT,       rt);
      addMsMs(rt, mz, COL_SCAN,     features.getInt("Scan number"));
      addMsMs(rt, mz, COL_Z,        features.getInt("Charge"));
      addMsMs(rt, mz, COL_MASS,     features.getDouble("Mass"));
      addMsMs(rt, mz, COL_AI_APEX,  features.getDouble("Precursor intensity"));
//      addMsMs(rt, mz, COL_PROPERTY, new HashMap<>(features.getMappedRow()));

      mScan_MS2.put(features.getInt("Scan number"), features.get("Raw file"), mMsMs.get(mz,rt));

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
    System.out.println("Reading the features from "+features.getFileName()); int counts=0;
    while (features.hasNext())
    {
      // Raw file        Type    Charge  m/z     Mass    Uncalibrated m/z        Resolution      Number of data points
      // Number of scans Number of isotopic peaks        PIF     Mass fractional part    Mass deficit    Mass precision [ppm]
      // Max intensity m/z 0     Retention time  Retention length        Retention length (FWHM) Min scan number
      // Max scan number Identified      MS/MS IDs       Sequence        Length  Modifications   Modified sequence
      // Proteins        Score   Intensity       Intensities     Isotope pattern MS/MS Count     MSMS Scan Numbers       MSMS Isotope Indices
      Double mz = features.getDouble("m/z"), rt = features.getDouble("Retention time");
      // add the features
      addFeature(rt, mz, COL_MZ,       mz);
      addFeature(rt, mz, COL_RT,       rt);
      addFeature(rt, mz, COL_Z,        features.getInt("Charge"));
      addFeature(rt, mz, COL_MASS,     features.getDouble("Mass"));
      addFeature(rt, mz, COL_AI_APEX,  features.getDouble("Intensity"));
      addFeature(rt, mz, COL_RT_BOUND, Range.closed(
          mScan_RT.get(features.getInt("Min scan number"), features.get("Raw file")),
          mScan_RT.get(features.getInt("Max scan number"), features.get("Raw file"))));
      addFeature(rt, mz, COL_SCAN_BOUND, Range.closed(features.getInt("Min scan number"), features.getInt("Max scan number")));
//      addFeature(rt, mz, COL_PROPERTY, new HashMap<>(features.getMappedRow()));

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

//        if (Math.abs(1251.636-mass)<0.01 && Math.abs(77-pk.getDouble(COL_RT))<2)
//          System.out.print("");

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
  // remember m/z first
  private LcMsMsFeatures addMsMs(Double rt, Double mz, String key, Object val)
  {
    if (mMsMs==null) mMsMs = TreeBasedTable.create();
    // add the feature if not present
    Features F = mMsMs.get(mz, rt);
    if (F==null)
    {
      F = new Features(); mMsMs.put(mz, rt, F);
    }
    F.add(key, val);

    return this;
  }

  public TreeMultimap<Double, Features> getFeaturesByIntensity()
  {
    if (!Tools.isSet(mIons)) return null;

    TreeMultimap<Double, Features> F = TreeMultimap.create(Ordering.natural().reverse(), Ordering.natural());
    for (Features ff : mIons.values())
      F.put(ff.getDouble("intensityApex"), ff);

    return F;
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
}
