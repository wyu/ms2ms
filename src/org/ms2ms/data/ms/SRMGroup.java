package org.ms2ms.data.ms;

import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.data.Point;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.Writer;
import java.util.*;

public class SRMGroup implements Ion, Comparable<SRMGroup>
{
  private String mPeptideSequence;
  private float mRT, mPrecursorMz;
  private int mCharge;

  private TreeMap<Float, Float> mTransitions;
  private TreeMultimap<Float, Point> mXIC;

  public SRMGroup() { super(); }
  public SRMGroup(String peptide, float rt, float mz, int z)
  {
    super();
    mPeptideSequence = peptide;
    mTransitions     = new TreeMap<>();
    mXIC             = TreeMultimap.create();

    mRT=rt; mPrecursorMz=mz; mCharge=z;
  }

  @Override public int compareTo(SRMGroup o) { return Double.compare(getMz(), o.getMz()); }

  @Override public float getMz()     { return mPrecursorMz; }
  @Override public float getMH()     { return 0; }
  @Override public int   getCharge() { return mCharge; }
  public float getRT() { return mRT; }
  public String getSequence() { return mPeptideSequence; }

  public TreeMap<Float, Float> getTransitions() { return mTransitions; }
  public Collection<Float> getFragmentMzs(Float f0, Float f1)
  {
    SortedMap<Float, Float> g = mTransitions.subMap(f0,f1);
    return (g!=null && g.size()>0) ? g.values() : null;
  }
  public SRMGroup setMz(  float s) { mPrecursorMz=s; return this; }
  public SRMGroup setCharge(int s) { mCharge=s; return this; }

  private SRMGroup addTransition(float frag, float intensity)
  {
    mTransitions.put(frag, intensity);
    return this;
  }
  private SRMGroup addXICPoint(float frag, float rt, float intensity)
  {
    mXIC.put(frag, new Point(rt, intensity));
    return this;
  }

  public SRMGroup composite()
  {
    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    for (Float frag : mXIC.keySet())
      for (Point pk : mXIC.get(frag))
      {
        rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));
      }

    // create the composite trace
    for (Float rt : rt_ai.keySet())
      mXIC.put(0f, new Point(rt, Math.pow(10d, Stats.sum(rt_ai.get(rt)))));

    return this;
  }
  public SRMGroup scanMS2(SortedMap<Double, Peak> peaks, float rt, Tolerance tol)
  {
    for (Float k : mTransitions.keySet())
    {
      SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(k), tol.getMax(k));
      if (pks!=null && pks.size()>0)
        addXICPoint(k, rt, (float )Peaks.IntensitySum(pks.values()));
    }
    return this;
  }
  public static void headerXIC(Writer w) throws IOException
  {
    w.write("Peptide\tz\tPrecMz\tRT\tFragMz\txic.rt\txic.ai\n");
  }
  public void printXIC(Writer w) throws IOException
  {
    for (Float frag : mXIC.keySet())
      for (Point pk : mXIC.get(frag))
      {
        w.write(getSequence()+"\t");
        w.write(getCharge()+"\t");
        w.write(Tools.d2s(getMz(),4)+"\t");
        w.write(getRT()+"\t");

        w.write(Tools.d2s(frag,4)+"\t");
        w.write(Tools.d2s(pk.getX(),3)+"\t");
        w.write(Tools.d2s(pk.getY(),2)+"\n");
      }
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
          group = new SRMGroup(peptide, tr.getFloat("NormalizedRetentionTime"), tr.getFloat("PrecursorMz"), tr.getInt("PrecursorCharge"));
          groups.put(group.getMz(), group.getRT(), group);
          peptide_group.put(peptide, group);
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
