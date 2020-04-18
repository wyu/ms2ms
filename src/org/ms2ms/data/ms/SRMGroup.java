package org.ms2ms.data.ms;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Range;
import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Similarity;
import org.ms2ms.data.Point;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.math.Points;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;
import toools.collections.Lists;

import java.io.IOException;
import java.io.Writer;
import java.util.*;

public class SRMGroup implements Ion, Comparable<SRMGroup>, Cloneable
{
  private String mPeptideSequence, mProteinId;
  private float mRT, mPrecursorMz, mDpSimilarity;
  private int mCharge;

//  private TreeMap<Float, Float> mTransitions;
//  private TreeMultimap<Float, Point> mXIC;
//  private TreeMap<Float, Point> mFeatures;
  private TreeMap<Float, SRM> mSRMs;

  public SRMGroup() { super(); }
  public SRMGroup(String peptide)
  {
    super();
    mPeptideSequence = peptide;
    mSRMs            = new TreeMap<>();
  }
  public SRMGroup(String peptide, float rt, float mz, int z)
  {
    super();
    mPeptideSequence = peptide;
    mSRMs            = new TreeMap<>();

    mRT=rt; mPrecursorMz=mz; mCharge=z;
  }

  @Override public int compareTo(SRMGroup o) { return Double.compare(getMz(), o.getMz()); }

  @Override public float getMz()     { return mPrecursorMz; }
  @Override public float getMH()     { return 0; }
  @Override public int   getCharge() { return mCharge; }

  public float  getRT()           { return mRT; }
  public float  getSimilarity()   { return mDpSimilarity; }
  public String getSequence()     { return mPeptideSequence; }
  public String getProteinId()    { return mProteinId; }
  public SRM    getCompositeSRM() { return mSRMs!=null?mSRMs.get(0f):null; }

  public SRM getCompositeProfile(boolean round2sec)
  {
    SRM cpo = getCompositeSRM();

    if (cpo!=null && Tools.isSet(cpo.getXIC()) && cpo.getFeature()!=null)
    {
      SRM profile = new SRM(cpo.getFragmentMz(), 0f);
      double rt0 = cpo.getFeature().getRT();
      for (LcMsPoint p : cpo.getXIC())
      {
        profile.addXIC(round2sec?(Math.round(60d*(p.getRT()-rt0))):(float )(p.getRT()-rt0), (float )p.getIntensity());
      }
      return profile;
    }
    return mSRMs!=null?mSRMs.get(0f):null;
  }

  public Map<Float, SRM> getSRMs() { return mSRMs; }

  public SRMGroup setMz(  float s) { mPrecursorMz=s; return this; }
  public SRMGroup setCharge(int s) { mCharge=s; return this; }
  public SRMGroup setSimilarity(float s) { mDpSimilarity=s; return this; }
  public SRMGroup setProteinId(String s) { mProteinId=s; return this; }

  private SRMGroup addTransition(float frag, float intensity)
  {
    mSRMs.put(frag, new SRM(frag, intensity));
    return this;
  }
  private LcMsPoint addXICPoint(float frag, float rt, Double intensity, Double mz, int scan)
  {
    if (mSRMs.get(frag)==null) mSRMs.put(frag, new SRM());
    return mSRMs.get(frag).addXIC(rt, intensity!=null? (float )intensity.doubleValue():0f, mz!=null? (float)mz.doubleValue():0f, scan);
  }

  public SRMGroup composite()
  {
    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    HashMap<Float, Integer>     rt_sn = new HashMap<>();
    for (Float frag : mSRMs.keySet())
      if (frag>0) // only the MS2 XIC
        for (LcMsPoint pk : mSRMs.get(frag).getXIC())
          if (pk.getY()>0)
          {
            rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));
            rt_sn.put((float )pk.getX(), pk.getScan());
          }

    // create the composite trace
    double n = mSRMs.size(); int ms1_start=0;
    SRM ms1 = mSRMs.get(-1f);
    mSRMs.put(0f, new SRM(0f, 0f));
    for (Float rt : rt_ai.keySet())
    {
      // look for the ms1 trace
      if (ms1!=null && Tools.isSet(ms1.getXIC()))
        for (int i=ms1_start; i<ms1.getXIC().size()-1; i++)
          if (rt>=ms1.getXIC().get(i).getX() && rt<=ms1.getXIC().get(i+1).getX())
          {
            rt_ai.put(rt, Math.log10(0.5*(ms1.getXIC().get(i).getY()+ms1.getXIC().get(i+1).getY()))); ms1_start=i+1; break;
          }

      float v = (float )Math.pow(10d, Stats.sum(rt_ai.get(rt))/n);
      mSRMs.get(0f).addXIC(rt, v, 0f, rt_sn.get(rt)*-1);
    }

    return this;
  }
  public SRMGroup impute(float gap)
  {
    for (SRM srm : mSRMs.values()) srm.impute(gap);
    return this;
  }
  public SRMGroup centroid(float min_ri, float rt_span)
  {
    Point cpo = Points.centroid(mSRMs.get(0f).getXIC(), 5d, Range.closed(0d, 1000d));
    if (cpo!=null)
    {
      Range<Double> rt_range = Range.closed(cpo.getX()-rt_span, cpo.getX()+rt_span);
      List<Point> injects = new ArrayList();
      for (Float frag : mSRMs.keySet())
      {
        SRM srm = mSRMs.get(frag);
        srm.setFeature(Points.centroid(srm.getXIC(), 5d, rt_range));

        if (frag>0)
        {
          injects.clear();
          for (LcMsPoint p : srm.getXIC())
            if (p.getIntensity()>0)
              injects.add(new Point(p.getFillTime(), p.getIntensity()));

          if (Tools.isSet(injects) && srm.getFeature()!=null)
            srm.getFeature().setFillTime(Points.centroid(injects));
        }
      }
    }
    return this;
  }
  public SRMGroup scoreSimillarity()
  {
    List<Float> lib = new ArrayList<>(), obs = new ArrayList<>();
    for (SRM srm : mSRMs.values())
      if (srm.getFeature()!=null && srm.getFragmentMz()>0 && srm.getLibraryIntensity()>0)
      {
        lib.add((float )Math.sqrt(srm.getLibraryIntensity()));
        obs.add((float )Math.sqrt(srm.getFeature().getY()));
      }

    mDpSimilarity = Similarity.dp(lib, obs);

    return this;
  }
  public SRMGroup calcFeatureExclusivity(float span, int apex_pts)
  {
    if (Tools.isSet(mSRMs) && mSRMs.get(0f)!=null && mSRMs.get(0f).getFeature()!=null)
    {
      double rt = mSRMs.get(0f).getFeature().getX();
      for (SRM srm : mSRMs.values()) srm.calPeakPct(rt, span, apex_pts);
    }
    return this;
  }
  public SRMGroup scanMS2(SortedMap<Double, Peak> peaks, float rt, int scan, float fill_time, Tolerance tol)
  {
    for (Float k : mSRMs.keySet())
    {
      if (k<=0) continue;

      SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(k), tol.getMax(k));
      if (pks!=null && pks.size()>0)
        addXICPoint(k, rt, Peaks.IntensitySum(pks.values()), Peaks.centroid(pks.values()), scan).setFillTime(fill_time);
      else
        addXICPoint(k, rt, 0d, 0d, scan).setFillTime(fill_time);
    }
    return this;
  }
  public SRMGroup scanMS1(SortedMap<Double, Peak> peaks, float rt, int scan, float fill_time, Tolerance tol)
  {
    SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(getMz()), tol.getMax(getMz()));
    if (pks!=null && pks.size()>0)
      addXICPoint(-1f, rt, Peaks.IntensitySum(pks.values()), Peaks.centroid(pks.values()), scan).setFillTime(fill_time);
    return this;
  }
  public static void headerXIC(Writer w) throws IOException
  {
    headerGroup(w);
    w.write("FragMz\tscan\txic.rt\txic.ai\tinjection\timputed\txic.mz\n");
  }
  public SRMGroup printXIC(Writer w) throws IOException
  {
    for (Float frag : mSRMs.keySet())
      for (LcMsPoint pk : mSRMs.get(frag).getXIC())
      {
        printGroup(w);
        w.write(Tools.d2s(frag,4)             +"\t");
        w.write(Tools.d2s(pk.getScan(),2)     +"\t");
        w.write(Tools.d2s(pk.getRT(),3)       +"\t");
        w.write(Tools.d2s(pk.getIntensity(),2)+"\t");
        w.write(Tools.d2s(pk.getFillTime(),2)+"\t");
        w.write(pk.isImputed()                   +"\t");
        w.write(Tools.d2s(pk.getMz(),3)       +"\n");
      }

    return this;
  }
  public static void headerFeatures(Writer w) throws IOException
  {
    headerGroup(w);
    w.write("FragMz\tNumPts\txicL\txicR\tPkEx\tPkExAll\tfeature.rt\tfeature.ai\tinjection\tfeature.apex\n");
  }
  public void printFeatures(Writer w) throws IOException
  {
    for (Float frag : mSRMs.keySet())
      if (mSRMs.get(frag).getFeature()!=null)
      {
        SRM srm = mSRMs.get(frag);
        printGroup(w);
        w.write(Tools.d2s(frag,4)+"\t");
        w.write(srm.getXIC().size()+"\t");
        w.write(Tools.d2s(srm.getXIC().get(0).getX(),3)+"\t");
        w.write(Tools.d2s(srm.getXIC().get(srm.getXIC().size()-1).getX(),3)+"\t");
        w.write(Tools.d2s(srm.getPeakPct(),2)+"\t");
        w.write(Tools.d2s(srm.getPeakPctAll(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getX(),3)+"\t");
        w.write(Tools.d2s(srm.getFeature().getY(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getFillTime(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getApex(),2)+"\n");
      }
  }
  public static void headerGroup(Writer w) throws IOException
  {
    w.write("Peptide\tz\tPrecMz\tRT\tSimilarity\t");
  }
  private void printGroup(Writer w) throws IOException
  {
    w.write(getSequence()+"\t");
    w.write(getCharge()+"\t");
    w.write(Tools.d2s(getMz(),4)+"\t");
    w.write(Tools.d2s(getRT(),3)+"\t");
    w.write(Tools.d2s(getSimilarity(),3)+"\t");
  }
  // replace the fragment mz of the transitions with random pick from all fragments
  public SRMGroup mutate(ListMultimap<Integer, Float> frag_bank, Random rnd)
  {
    mPeptideSequence = getSequence()+"_DECOY";
    Collection<Float> frags = new ArrayList<>(mSRMs.keySet());

    mSRMs.clear(); mSRMs = new TreeMap<>();
    for (Float frag : frags)
    {
      addTransition(Lists.getRandomSubset(frag_bank.get((int )Math.round(frag*0.01)),1, rnd).get(0), 0f);
    }
    return this;
  }

  @Override
  public SRMGroup clone()
  {
    SRMGroup cloned = new SRMGroup(getSequence(), getRT(), getMz(), getCharge());

    cloned.setSimilarity(getSimilarity());

    cloned.mSRMs = new TreeMap<>();
    for (Float frag : mSRMs.keySet()) cloned.mSRMs.put(frag, mSRMs.get(frag));

    return cloned;
  }

  @Override
  public String toString()
  {
    String out = (Strs.isSet(getProteinId())?(getProteinId()+"::"):"")+getSequence() + ", m/z" + Tools.d2s(getMz(), 4) + ", " +
        Tools.d2s(getRT(), 2) + "min, dp=" + Tools.d2s(getSimilarity(),2) + " ";

    if (Tools.isSet(mSRMs))
      for (Float frag : mSRMs.keySet())
      {
        if (!Tools.isSet(mSRMs.get(frag).getXIC())) continue;

        if (frag==-1f) out += ("ms1$"); else if (frag==0f) out += ("cpo$"); else  out += (Tools.d2s(frag, 3)+"$");

        out += (mSRMs.get(frag).getXIC().size()+";");
      }

    return out;
  }

  public static MultiTreeTable<Float, Float, SRMGroup> readTransitions(String trfile)
  {
    try
    {
      MultiTreeTable<Float, Float, SRMGroup> groups = new MultiTreeTable<Float, Float, SRMGroup>();

      TabFile                          tr = new TabFile(trfile, ",");
      Map<String, SRMGroup> peptide_group = new HashMap<>();

      // going thro the rows
      long row_counts = 0;
      while (tr.hasNext())
      {
        if (++row_counts % 10000  ==0) System.out.print(".");

        // "","PrecursorMz","ProductMz","LibraryIntensity","ProteinId","PeptideSequence","ModifiedPeptideSequence","PrecursorCharge","ProductCharge","FragmentType","FragmentSeriesNumber","NormalizedRetentionTime"
        // "2261",1044.48640687972,405.176849365234,33.29616,"P62258","AAFDDAIAELDTLSEESYK","AAFDDAIAELDTLSEESYK",2,1,"b",4,92.1534957885742
        int z = tr.get("PrecursorCharge", 0);
        String peptide = tr.getStr("ModifiedSequence", "ModifiedPeptideSequence")+(z!=0?("#"+z):"");

        SRMGroup group = peptide_group.get(peptide);

        if (group==null) {
          try
          {
            group = new SRMGroup(peptide, tr.get("NormalizedRetentionTime", 0f), tr.getFloat("PrecursorMz"), tr.get("PrecursorCharge", 0));
            if (tr.get("ProteinId")!=null) group.setProteinId(tr.get("ProteinId"));

            groups.put(group.getMz(), group.getRT(), group);
            peptide_group.put(peptide, group);
          }
          catch (NullPointerException e) {
            System.out.println();
          }
        }
        group.addTransition(tr.getFloat("ProductMz"), tr.get("LibraryIntensity", 0f));
      }
      tr.close();

      return groups;
    }
    catch (IOException e) {}

    return null;
  }
  public static SRMGroup buildProteinProfile(Collection<SRMGroup> groups, String proteinid, boolean round2sec, String... peptides)
  {
    SRMGroup protein = new SRMGroup("Summary: " + proteinid);

    // create the composite trace
    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    for (SRMGroup group : groups)
      if (group.getCompositeSRM()!=null && Strs.isA(group.getSequence(), peptides))
      {
        SRM pr = group.getCompositeProfile(round2sec);
        protein.getSRMs().put(group.getMz(), pr);
        for (LcMsPoint pk : pr.getXIC())
          if (pk.getY()>0)
          {
            // convert to seconds
            rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));
          }
      }

    double n = 0;
    for (Float rt : rt_ai.keySet())
      if (rt_ai.get(rt).size()>n) n = rt_ai.get(rt).size();

    // create the composite trace
    SRM profile = new SRM(0f, 0f);
    for (Float rt : rt_ai.keySet())
    {
      // look for the ms1 trace
      float v = (float )Math.pow(10d, Stats.sum(rt_ai.get(rt))/n);
      profile.addXIC(rt, v);
    }
    protein.getSRMs().put(0f, profile);

    return protein;
  }
}
