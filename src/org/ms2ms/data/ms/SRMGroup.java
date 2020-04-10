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
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;
import toools.collections.Lists;

import java.io.IOException;
import java.io.Writer;
import java.util.*;

public class SRMGroup implements Ion, Comparable<SRMGroup>, Cloneable
{
  private String mPeptideSequence;
  private float mRT, mPrecursorMz, mDpSimilarity;
  private int mCharge;

//  private TreeMap<Float, Float> mTransitions;
//  private TreeMultimap<Float, Point> mXIC;
//  private TreeMap<Float, Point> mFeatures;
  private TreeMap<Float, SRM> mSRMs;

  public SRMGroup() { super(); }
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
  public float getRT() { return mRT; }
  public float getSimilarity() { return mDpSimilarity; }
  public String getSequence() { return mPeptideSequence; }

  public Map<Float, SRM> getSRMs() { return mSRMs; }

  public SRMGroup setMz(  float s) { mPrecursorMz=s; return this; }
  public SRMGroup setCharge(int s) { mCharge=s; return this; }
  public SRMGroup setSimilarity(float s) { mDpSimilarity=s; return this; }

  private SRMGroup addTransition(float frag, float intensity)
  {
    mSRMs.put(frag, new SRM(frag, intensity));
    return this;
  }
  private SRMGroup addXICPoint(float frag, float rt, float intensity)
  {
    if (intensity>0)
    {
      if (mSRMs.get(frag)==null) mSRMs.put(frag, new SRM());
      mSRMs.get(frag).addXIC(rt, intensity);
    }
    return this;
  }

  public SRMGroup composite()
  {
    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    for (Float frag : mSRMs.keySet())
      if (frag>0) // only the MS2 XIC
        for (Point pk : mSRMs.get(frag).getXIC())
          rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));

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
      mSRMs.get(0f).addXIC(rt, v);
    }

    return this;
  }
  public SRMGroup centroid(float min_ri, float rt_span)
  {
    Point cpo = Points.centroid(mSRMs.get(0f).getXIC(), 5d, Range.closed(0d, 1000d));
    if (cpo!=null)
    {
      Range<Double> rt_range = Range.closed(cpo.getX()-rt_span, cpo.getX()+rt_span);
      for (Float frag : mSRMs.keySet())
      {
        SRM srm = mSRMs.get(frag);
        srm.setFeature(Points.centroid(srm.getXIC(), 5d, rt_range));
      }
    }
    return this;
  }
  public SRMGroup scoreSimillarity()
  {
    List<Float> lib = new ArrayList<>(), obs = new ArrayList<>();
    for (SRM srm : mSRMs.values())
      if (srm.getFeature()!=null && srm.getFragmentMz()>0)
      {
        lib.add(        srm.getLibraryIntensity());
        obs.add((float )srm.getFeature().getY());
      }

    mDpSimilarity = Similarity.dp(lib, obs);

    return this;
  }
  public SRMGroup calcFeatureExclusivity(float span)
  {
    if (Tools.isSet(mSRMs) && mSRMs.get(0f)!=null && mSRMs.get(0f).getFeature()!=null)
    {
      double rt = mSRMs.get(0f).getFeature().getX();
      for (SRM srm : mSRMs.values()) srm.calPeakPct(rt, span);
    }
    return this;
  }
  public SRMGroup scanMS2(SortedMap<Double, Peak> peaks, float rt, Tolerance tol)
  {
    for (Float k : mSRMs.keySet())
    {
      if (k<=0) continue;

      SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(k), tol.getMax(k));
      if (pks!=null && pks.size()>0)
        addXICPoint(k, rt, (float )Peaks.IntensitySum(pks.values()));
    }
    return this;
  }
  public SRMGroup scanMS1(SortedMap<Double, Peak> peaks, float rt, Tolerance tol)
  {
    SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(getMz()), tol.getMax(getMz()));
    if (pks!=null && pks.size()>0)
      addXICPoint(-1f, rt, (float )Peaks.IntensitySum(pks.values()));
    return this;
  }
  public static void headerXIC(Writer w) throws IOException
  {
    headerGroup(w);
    w.write("FragMz\txic.rt\txic.ai\n");
  }
  public SRMGroup printXIC(Writer w) throws IOException
  {
    for (Float frag : mSRMs.keySet())
      for (Point pk : mSRMs.get(frag).getXIC())
      {
        printGroup(w);
        w.write(Tools.d2s(frag,4)+"\t");
        w.write(Tools.d2s(pk.getX(),3)+"\t");
        w.write(Tools.d2s(pk.getY(),2)+"\n");
      }

    return this;
  }
  public static void headerFeatures(Writer w) throws IOException
  {
    headerGroup(w);
    w.write("FragMz\tNumPts\txicL\txicR\tPkEx\tPkExAll\tfeature.rt\tfeature.ai\n");
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
        w.write(Tools.d2s(srm.getFeature().getY(),2)+"\n");
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

  class SRM implements Cloneable
  {
    private float mFragmentMz, mLibraryIntensity, mApex, mArea, mPkPct=0, mPkPctAll=0;
    private List<Point> mXIC;
    private Point mFeature;

    SRM()
    {
      super();
      mFragmentMz=mLibraryIntensity=0;
      mXIC = new ArrayList<>();
    }
    SRM(float frag, float ai)
    {
      mFragmentMz=frag; mLibraryIntensity=ai;
      mXIC = new ArrayList<>();
    }

    public float getFragmentMz() { return mFragmentMz; }
    public float getLibraryIntensity() { return mLibraryIntensity; }
    public float getApex()       { return mApex; }
    public float getArea()       { return mArea; }
    public float getPeakPct()    { return mPkPct; }
    public float getPeakPctAll() { return mPkPctAll; }

    public List<Point> getXIC() { return mXIC; }
    public Point getFeature() { return mFeature; }

    public SRM addXIC(float rt, float ai) { if (ai>0) mXIC.add(new Point(rt,ai)); return this; }

    public SRM setFeature(Point s) { mFeature=s; return this; }

    public SRM calPeakPct(double rt, double span)
    {
      if (!Tools.isSet(getXIC())) return this;

      double inside=0, outside=0, all=0;
      Range<Double> inner = Range.closed(rt-span, rt+span), outer = Range.closed(rt-2d*span, rt+2d*span);

      for (Point p : getXIC())
      {
        if       (inner.contains(p.getX())) inside+=p.getY();
        else if (!outer.contains(p.getX())) outside+=p.getY();
        all += p.getY();
      }
      mPkPct    = (float )(100f*inside/(inside+outside));
      mPkPctAll = (float )(100f*inside/all);
      return this;
    }

    public SRM clone()
    {
      SRM cloned      = new SRM(getFragmentMz(), getLibraryIntensity());
      cloned.mApex    = mApex; cloned.mArea = mArea; cloned.mPkPct = mPkPct; cloned.mPkPctAll = mPkPctAll;
      cloned.mFeature = mFeature;

      cloned.mXIC     = new ArrayList<>();
      for (Point p : mXIC) cloned.mXIC.add(p);

      return cloned;
    }
    public SRM mutate(ListMultimap<Integer, Float> frag_bank, Random rnd)
    {
      mApex=mArea=mPkPct=mPkPctAll = 0f;
      mXIC.clear();

      Integer idx = (int )Math.round(getFragmentMz()*0.01);
      mFragmentMz = Lists.getRandomSubset(frag_bank.get(idx),1, rnd).get(0);

      return this;
    }
    @Override
    public String toString()
    {
      String out = "";
      if (Tools.isSet(getXIC())) out = "xic="+getXIC().size();
      if (mFeature!=null)        out += ", " + mFeature.toString();

      return out;
    }
  }

  @Override
  public String toString()
  {
    String out = getSequence() + ", m/z" + Tools.d2s(getMz(), 4) + ", " + Tools.d2s(getRT(), 2) + "min, dp=" +
        Tools.d2s(getSimilarity(),2) + " ";

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
        String peptide = tr.get("ModifiedPeptideSequence")+"#"+tr.get("PrecursorCharge");
        SRMGroup group = peptide_group.get(peptide);

        if (group==null) {
          try
          {
            group = new SRMGroup(peptide, tr.getFloat("NormalizedRetentionTime"), tr.getFloat("PrecursorMz"), tr.getInt("PrecursorCharge"));
            groups.put(group.getMz(), group.getRT(), group);
            peptide_group.put(peptide, group);
          }
          catch (NullPointerException e) {
            System.out.println();
          }
        }
        group.addTransition(tr.getFloat("ProductMz"), tr.getFloat("LibraryIntensity"));
      }
      tr.close();

      return groups;
    }
    catch (IOException e) {}

    return null;
  }
}
