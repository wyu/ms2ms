package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;
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

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.*;

public class SRMGroup implements Ion, Comparable<SRMGroup>, Cloneable
{
  private String mPeptideSequence, mProteinId;
  private float mRT, mPrecursorMz, mDpSimilarity, mIRT, mReportedRT, mNetworkNiche;
  private int mCharge, mQualifiedSRMs;
  private Map<String, Object> mNetworkStats;

  private Point mCompositeFeature0=null;
  private List<Peak> mCompositePeaks=null;

//  private TreeMap<Float, Float> mTransitions;
//  private TreeMultimap<Float, Point> mXIC;
//  private TreeMap<Float, Point> mFeatures;
  private TreeMap<Float, SRM> mSRMs;
  private SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> mNetwork;

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

  @Override public int compareTo(SRMGroup o)
  {
    int c = Double.compare(getMz(), o.getMz());
    if (c==0) c = Float.compare(getRT(), o.getRT());
    if (c==0) c = mPeptideSequence.compareTo(o.getSequence());

    if (c==0 && Tools.isSet(mSRMs) && Tools.isSet(o.getSRMs()))
    {
      c = Integer.compare(getSRMs().size(), o.getSRMs().size());
    }
    return c;
  }

  @Override public float getMz()     { return mPrecursorMz; }
  @Override public float getMH()     { return 0; }
  @Override public int   getCharge() { return mCharge; }

  public int    getNumQualifiedSRMs() { return mQualifiedSRMs; }

  public float  getRT(boolean iRT) { return iRT?mIRT:mRT; }

  public float  getRT()           { return mRT; }
  public float  getIRT()          { return mIRT; }
  public float  getReportedRT()   { return mReportedRT; }
  public float  getSimilarity()   { return mDpSimilarity; }
  public String getSequence()     { return mPeptideSequence; }
  public String getProteinId()    { return mProteinId; }
  public SRM    getCompositeSRM() { return mSRMs!=null?mSRMs.get(0f):null; }

  public SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> getNetwork() { return mNetwork; }

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
  public SRM getSRM(Float s) { return mSRMs!=null?mSRMs.get(s):null; }

  public SRMGroup setNumQualifiedSRMs(int s) { mQualifiedSRMs=s; return this; }

  public SRMGroup setMz(  float s) { mPrecursorMz=s; return this; }
  public SRMGroup setRT(  float s) { mRT=s; return this; }
  public SRMGroup setRT( Double s) { if (s!=null) mRT=(float )s.doubleValue(); return this; }

  public SRMGroup setIRT( float s) { mIRT=s; return this; }
  public SRMGroup setReportedRT( float s) { mReportedRT=s; return this; }

  public SRMGroup setCharge(int s) { mCharge=s; return this; }
  public SRMGroup setSimilarity(float s) { mDpSimilarity=s; return this; }
  public SRMGroup setProteinId(String s) { mProteinId=s; return this; }

  public float calcSrmYield(float span)
  {
    if (Tools.isSet(getSRMs()))
    {
      float expected = getSRMs().keySet().size(), found=0f, rt0=getRT();

      if (getSRMs().containsKey(0f) && getSRMs().get(0f).getFeature()!=null)
      {
        rt0 = (float )getSRMs().get(0f).getFeature().getRT();
        for (SRM srm : getSRMs().values())
          if (srm.getFragmentMz()<=0f) expected--;
          else
          {
            if (srm.getFeature()!=null &&
                ((srm.getPeakBoundary()!=null && srm.getPeakBoundary().contains(srm.getFeature().getRT())) ||
                    Math.abs(srm.getFeature().getRT()-rt0)<=span)) found++;
          }
      }
      mQualifiedSRMs = Math.round(found);

      return (expected!=0?(100f*found/expected):0f);
    }
    // how much of the SRM were productive
    return 0f;
  }
  private SRMGroup addTransition(float frag, float intensity)
  {
    mSRMs.put(frag, new SRM(frag, intensity));
    return this;
  }
  private LcMsPoint addXICPoint(float frag, float rt, Double intensity, Double mz, int scan, Double fill, boolean keep_zero)
  {
    if (keep_zero || intensity>0)
    {
      if (mSRMs.get(frag)==null) mSRMs.put(frag, new SRM());
      return mSRMs.get(frag).addXIC(rt, intensity!=null? (float )intensity.doubleValue():0f, mz!=null? (float)mz.doubleValue():0f, scan).setFillTime(fill);
    }
    return null;
  }
  public SRMGroup disposeSRMs()
  {
    if (Tools.isSet(mSRMs))
      for (Float frag  : mSRMs.keySet())
        if (frag!=0) getSRM(frag).disposeXIC(); // 20200605, keep the composite alive for protein IDs

    return this;
  }
  // shift the RT of ms1 XIC so they line up with those of ms2
  public SRMGroup shift_ms1()
  {
    if (!Tools.isSet(mSRMs) || mSRMs.get(-1f)==null) return this;

    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    for (Float frag : mSRMs.keySet())
      if (frag>0) // only the MS2 XIC
        for (LcMsPoint pk : mSRMs.get(frag).getXIC())
          if (pk.getY()>0)
          {
            rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));
          }

    // create the composite trace
    int ms1_start=0;
    SRM ms1 = mSRMs.get(-1f), ms1a = ms1.clone(); ms1a.getXIC().clear();
    for (Float rt : rt_ai.keySet())
    {
      // look for the ms1 trace
      if (ms1!=null && Tools.isSet(ms1.getXIC()))
        for (int i=ms1_start; i<ms1.getXIC().size()-1; i++)
          if (rt>=ms1.getXIC().get(i).getX() && rt<=ms1.getXIC().get(i+1).getX())
          {
            ms1a.addXIC(rt.doubleValue(), 0.5*(ms1.getXIC().get(i).getY()+ms1.getXIC().get(i+1).getY()), ms1.getXIC().get(i).getMz(),
                ms1.getXIC().get(i).getScan());
            break;
          }
    }
    mSRMs.put(-10f, ms1a);

    return this;
  }
  public SRMGroup networking(float min_dp)
  {
    if (!Tools.isSet(getSRMs())) return this;

    // create a new network
    mNetwork = new SimpleDirectedWeightedGraph<>(DefaultWeightedEdge.class);
    List<Float> traces = new ArrayList<>(getSRMs().keySet());
    // remove the default ms1 sinceit's not sync with the ms2 on RT. We should have called shift_ms1 first to create a sync'd version @-10
    traces.remove(-1f);
    // remove the composite trace from the network since it's redundant for this purpose
    traces.remove(0f);

    for (int i=0; i<traces.size(); i++)
    {
      SRM x = getSRM(traces.get(i));
      if (!Tools.isSet(x.getXIC()) || x.getXIC().size()<5) continue;

        mNetwork.addVertex(x);
      for (int j=0; j<traces.size(); j++)
      {
        SRM y = getSRM(traces.get(j));
        if (!Tools.isSet(y.getXIC()) || y.getXIC().size()<5) continue;

        mNetwork.addVertex(y);
        if (i!=j)
        {
          // calc the similarity
          double dp = Similarity.dp_points(x.getXIC(), y.getXIC());
          if (dp>=min_dp) mNetwork.setEdgeWeight(mNetwork.addEdge(x,y), dp);
        }
      }
    }
    if (mNetwork!=null && Tools.isSet(mNetwork.edgeSet())) mNetworkStats = SRM.inspectNetwork(mNetwork);

    return this;
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
  public SRMGroup fill(float baseline)
  {
    Multimap<Float, Float> Xs = TreeMultimap.create();
    for (SRM srm : mSRMs.values())
      if (srm.getFragmentMz()>0 && Tools.isSet(srm.getXIC()))
        for (LcMsPoint pt : srm.getXIC())
        {
          Xs.put(0f,              (float )pt.getRT());
          Xs.put(srm.getFragmentMz(), (float )pt.getRT());
        }

    for (SRM srm : mSRMs.values()) srm.fill(0f, Xs);

    Xs = Tools.dispose(Xs);

    return this;
  }
  public SRMGroup impute(float gap)
  {
    for (SRM srm : mSRMs.values()) srm.impute(gap);
    return this;
  }
  public SRMGroup centroid(LcSettings settings)
  {
    return centroid(settings.getBaseRI(), settings.getPeakWidth()*settings.getOuterMultiple(), settings.getPeakWidth());
  }
  public SRMGroup centroid(double min_ri, float rt_span, float rt_width)
  {
    Point cpo = Points.centroid(mSRMs.get(0f).getXIC(), min_ri, Range.closed(getRT()-rt_span*2d, getRT()+2d*rt_span));
    if (cpo!=null) {
      // check the alternative peaks by 1st derivatives
      mCompositePeaks = mSRMs.get(0f).detectPeak();
      if (Tools.isSet(mCompositePeaks)) {
        double best = Double.MAX_VALUE;
        Peak closed2expected = null, centroid_found = null;
        for (Peak pk : mCompositePeaks) {
          if (Math.abs(pk.getMz() - cpo.getX()) <= rt_span / 4d) centroid_found = pk;
          if (Math.abs(pk.getMz() - getRT()) < best) {
            closed2expected = pk;
            best = Math.abs(pk.getMz() - getRT());
          }
        }
        if (closed2expected != null && centroid_found != closed2expected)
        {
          mCompositeFeature0 = cpo;
          cpo = Points.centroid(mSRMs.get(0f).getXIC(), min_ri, Range.closed(closed2expected.getMz() - rt_width / 2d, closed2expected.getMz() + rt_width / 2d));
        }
      }
//      if (cpo==null)
//      {
//        System.out.println();
//      }

    }
    if (cpo!=null)
    {
      Range<Double> rt_range = Range.closed(cpo.getX()-rt_span, cpo.getX()+rt_span);
      List<Point> injects = new ArrayList(), mzs = new ArrayList<>();
      for (Float frag : mSRMs.keySet())
      {
        SRM srm = mSRMs.get(frag);
        srm.setFeature(Points.centroid(srm.getXIC(), min_ri, rt_range));

        if (frag>0)
        {
          injects.clear();  mzs.clear();
          for (LcMsPoint p : srm.getXIC())
            if (p.getIntensity()>0)
            {
              injects.add(new Point(p.getFillTime(), p.getIntensity()));
              mzs.add(    new Point(p.getMz(), p.getIntensity()));
            }

          if (Tools.isSet(injects) && srm.getFeature()!=null)
            srm.getFeature().setFillTime(Points.centroid(injects));
          if (Tools.isSet(mzs) && srm.getFeature()!=null)
            srm.getFeature().setMz(      Points.centroid(mzs));
        }
      }
    }
    return this;
  }
  public SRMGroup centroidByPeakBoundry()
  {
    List<Point> pts = new ArrayList(), mzs = new ArrayList();
    for (Float frag : mSRMs.keySet())
    {
      SRM srm = mSRMs.get(frag);
      if (srm.getFeature()!=null && Tools.isSet(srm.getPeakBoundary()))
      {
        pts.clear();
        for (LcMsPoint p : srm.getXIC())
          if (p.getIntensity()>0 && srm.getPeakBoundary().contains(p.getX()))
          {
            pts.add(new Point(p.getX(), p.getIntensity()));
            mzs.add(new Point(p.getMz(), p.getIntensity()));
          }

        if (Tools.isSet(pts))
        {
          srm.getFeature().setInitialCentroidRt(srm.getFeature().getX());
          srm.getFeature().setX( Points.centroid(pts));
          srm.getFeature().setMz(Points.centroid(mzs));
          // calc the stdev
          List<Double> devi = new ArrayList<>(mzs.size());
          for (Point p : mzs) devi.add(1E6*(p.getX()-srm.getFeature().getMz())/srm.getFeature().getMz());
          srm.getFeature().setMzStdev(Stats.stdev(devi));

          devi = (List )Tools.dispose(devi);
        }
      }
    }
    pts = (List )Tools.dispose(pts);
    mzs = (List )Tools.dispose(mzs);

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
  public SRMGroup calcFeatureExclusivity(LcSettings settings)
  {
    return calcFeatureExclusivity(settings.getPeakWidth()/2f, settings.getApexPts(), settings.getBaseRI()/100d);
  }
  public SRMGroup calcFeatureExclusivity(float span, int apex_pts, double peak_base)
  {
    if (Tools.isSet(mSRMs) && mSRMs.get(0f)!=null && mSRMs.get(0f).getFeature()!=null)
    {
      double rt = mSRMs.get(0f).getFeature().getX();
      for (SRM srm : mSRMs.values()) srm.calPeakPct(rt, span, apex_pts, peak_base);
    }
    return this;
  }
  public SRMGroup scanMS2(SortedMap<Double, Peak> peaks, float rt, int scan, Double fill_time, Tolerance tol, boolean keep_zero)
  {
    for (Float k : mSRMs.keySet())
    {
      if (k<=0) continue;

      SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(k), tol.getMax(k));
      if (pks!=null && pks.size()>0)
        addXICPoint(k, rt, Peaks.IntensitySum(pks.values()), Peaks.centroid(pks.values()), scan, fill_time, keep_zero);
      else
        addXICPoint(k, rt, 0d, 0d, scan, fill_time, keep_zero);
    }
    return this;
  }
  public SRMGroup scanMS1(SortedMap<Double, Peak> peaks, float rt, int scan, Double fill_time, Tolerance tol, boolean keep_zero)
  {
    SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(getMz()), tol.getMax(getMz()));
    if (pks!=null && pks.size()>0)
    {
      addXICPoint(-1f, rt, Peaks.IntensitySum(pks.values()), Peaks.centroid(pks.values()), scan, fill_time, keep_zero);
    }
    return this;
  }
  public static void headerXIC(Writer w) throws IOException
  {
    headerGroup(w);
    w.write("FragMz\tscan\txic.rt\txic.ai\tinjection\timputed\txic.mz\n");
  }
  public static void header_xic(Writer w) throws IOException
  {
    // xic.ai~xic.rt|FragMz
    w.write("\"Peptide\\tFragMz\txic.rt\txic.ai\n");
  }
  public SRMGroup print_xic(Writer w, float... frags) throws IOException
  {
    if (!Tools.isSet(frags)) frags = Tools.toFloatArray(mSRMs.keySet());
    for (Float frag : frags)
      for (LcMsPoint pk : mSRMs.get(frag).getXIC())
      {
        if (pk.getIntensity()<=0) continue;

        w.write(getSequence()+"\t");
        w.write(Tools.d2s(frag,4)             +"\t");
        w.write(Tools.d2s(pk.getRT(),3)       +"\t");
        w.write(Tools.d2s(pk.getIntensity(),2)+"\n");
      }

    return this;
  }

  public SRMGroup printXIC(Writer w, float... frags) throws IOException
  {
    if (!Tools.isSet(mSRMs)) return this;

    if (!Tools.isSet(frags)) frags = Tools.toFloatArray(mSRMs.keySet());
    for (Float frag : frags)
      if (mSRMs.get(frag)!=null && Tools.isSet(mSRMs.get(frag).getXIC()))
        for (LcMsPoint pk : mSRMs.get(frag).getXIC())
        {
          if (pk.getIntensity()<=0) continue;

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
    w.write("FragMz\tNumPts\txicL\txicR\tPkEx\tPkExFull\tPkExAll\tfeature.rt\tfeature.ai\tfeature.rt0\tinjection\tfeature.mz\tfeature.apex\tfeature.area\tppm.stdev\tlower\tupper\n");
  }
  public SRMGroup printFeatures(Writer w) throws IOException
  {
    for (Float frag : mSRMs.keySet())
      if (mSRMs.get(frag).getFeature()!=null && mSRMs.get(frag)!=null && Tools.isSet(mSRMs.get(frag).getXIC()))
      {
        SRM srm = mSRMs.get(frag);

        printGroup(w);
        w.write(Tools.d2s(frag,4)+"\t");
        w.write(srm.getXIC().size()+"\t");
        w.write(Tools.d2s(srm.getXIC().get(0).getX(),3)+"\t");
        w.write(Tools.d2s(srm.getXIC().get(srm.getXIC().size()-1).getX(),3)+"\t");
        w.write(Tools.d2s(srm.getPeakPct(),2)+"\t");
        w.write(Tools.d2s(srm.getPeakPctFull(),2)+"\t");
        w.write(Tools.d2s(srm.getPeakPctAll(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getX(),3)+"\t");
        w.write(Tools.d2s(srm.getFeature().getY(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getInitialCentroidRt(),3)+"\t");
        w.write(Tools.d2s(srm.getFeature().getFillTime(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getMz(),4)+"\t");
        w.write(Tools.d2s(srm.getFeature().getApex(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getArea(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getMzStdev(),2)+"\t");
        w.write(Tools.d2s(srm.getPeakBoundary()!=null?srm.getPeakBoundary().lowerEndpoint():0,2)+"\t");
        w.write(Tools.d2s(srm.getPeakBoundary()!=null?srm.getPeakBoundary().upperEndpoint():0,2)+"\n");
      }
    return this;
  }
  public static void headerGroup(Writer w) throws IOException
  {
    w.write("Peptide\tz\tPrecMz\tRT\tSimilarity\tYield\tNetV\tNetE\tSCC\tSCC1\t");
  }
  private void printGroup(Writer w) throws IOException
  {
    w.write(getSequence()+"\t");
    w.write(getCharge()+"\t");
    w.write(Tools.d2s(getMz(),4)+"\t");
    w.write(Tools.d2s(getRT(),3)+"\t");
    w.write(Tools.d2s(getSimilarity(),3)+"\t");
    w.write(Tools.d2s(calcSrmYield(0.5f), 2)+"\t");
    w.write((mNetwork!=null?mNetwork.vertexSet().size():0)+"\t");
    w.write((mNetwork!=null?mNetwork.edgeSet().size():0)+"\t");
    w.write((mNetworkStats!=null?((int )mNetworkStats.get("SCC: length")):0)+"\t");
    w.write((mNetworkStats!=null && mNetworkStats.get("SCC: node size")!=null? ((int[] )mNetworkStats.get("SCC: node size"))[0]:0)+"\t");
  }
  // replace the fragment mz of the transitions with random pick from all fragments
  public SRMGroup mutate(ListMultimap<Integer, Float> frag_bank, Random rnd)
  {
    mPeptideSequence = getSequence()+"_DECOY";
    Collection<Float> frags = new ArrayList<>(mSRMs.keySet());

    mSRMs = (TreeMap )Tools.dispose(mSRMs); mSRMs = new TreeMap<>();
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
        Tools.d2s(getRT(), 2) + "min";

    int N=0;
    if (Tools.isSet(mSRMs))
    {
      for (Float frag : mSRMs.keySet())
        if (Tools.isSet(mSRMs.get(frag).getXIC())) N++;

//      for (Float frag : mSRMs.keySet())
//      {
//        if (!Tools.isSet(mSRMs.get(frag).getXIC())) continue;
//
//        if (frag==-1f) out += ("ms1$"); else if (frag==0f) out += ("cpo$"); else  out += (Tools.d2s(frag, 3)+"$");
//
//        out += (mSRMs.get(frag).getXIC().size()+";");
//      }
      return out+", N=" + N + " ";
    }

    return out;
  }

  public static MultiTreeTable<Float, Float, SRMGroup> readTransitions(String trfile, String delim, Map<String, String> cols, boolean use_iRT, String... protein_ids)
  {
    try
    {
      MultiTreeTable<Float, Float, SRMGroup> groups = new MultiTreeTable<>();

      TabFile                          tr = new TabFile(trfile, delim).setColMappings(cols);
      Map<String, SRMGroup> peptide_group = new HashMap<>();

      // going thro the rows
      long row_counts = 0;
      while (tr.hasNext())
      {
        if (++row_counts % 10000  ==0) System.out.print(".");

        // "","PrecursorMz","ProductMz","LibraryIntensity","ProteinId","PeptideSequence","ModifiedPeptideSequence","PrecursorCharge","ProductCharge","FragmentType","FragmentSeriesNumber","NormalizedRetentionTime"
        // "2261",1044.48640687972,405.176849365234,33.29616,"P62258","AAFDDAIAELDTLSEESYK","AAFDDAIAELDTLSEESYK",2,1,"b",4,92.1534957885742
        int z = tr.get("PrecursorCharge", 0);
        String seq = tr.getStr("ModifiedSequence", "ModifiedPeptideSequence"), peptide;

//        if (!(seq.indexOf("ELGTVM[Oxidation (M)]R#2")>=0 && z==2)) continue;
        if (Tools.isSet(protein_ids) && !Strs.hasA(tr.get("ProteinId"), 0, protein_ids)) continue;
        if (Strs.isSet(seq))
        {
          peptide = seq.replaceAll("_","");
          if (seq.indexOf('#')<0) peptide = peptide+(z!=0?("#"+z):"");
        }
        else
        {
          throw new RuntimeException("Peptide sequence not found in the transition list!");
        }
        SRMGroup group = peptide_group.get(peptide);

        if (group==null) {
          try
          {
            group = new SRMGroup(peptide, tr.get("NormalizedRetentionTime", 0f), tr.getFloat("PrecursorMz"), tr.get("PrecursorCharge", 0));
            if (tr.get("ProteinId")!=null) group.setProteinId(tr.get("ProteinId"));
            if (tr.get("iRT")!=null) group.setIRT(tr.get("iRT", 0f));
            if (tr.get("ReportedRT")!=null) group.setReportedRT(tr.get("ReportedRT", 0f));

            groups.put(group.getMz(), group.getRT(use_iRT), group);
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
    catch (IOException e) { e.printStackTrace(); }

    return null;
  }
  public static ProteinID buildProteinProfile(LcSettings settings, Collection<SRMGroup> groups, String proteinid, String... peptides)
  {
    ProteinID protein_id = new ProteinID(null, proteinid);
    SRMGroup protein     = new SRMGroup("Summary");

    // create the composite trace
    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    Range<Float> window = Range.closed(settings.getSpan()*-1, settings.getSpan());
    for (SRMGroup group : groups)
      if (Strs.isA(group.getSequence(), peptides))
      {
        protein_id.addSRMGroup(group, group.getSequence());
        if (group.getCompositeSRM()!=null && Tools.isSet(group.getCompositeSRM().getXIC()))
        {
          protein.getSRMs().put(group.getMz(), group.getCompositeSRM().shift(window, settings.getGridSize(),
              (float )group.getCompositeSRM().getFeature().getRT()));
//          SRM pr = group.getCompositeProfile(round2sec);
//          if (pr!=null && Tools.isSet(pr.getXIC()))
//          {
//            protein.getSRMs().put(group.getMz(), pr);
////            protein_id.addSRMGroup(group, group.getSequence());
//            for (LcMsPoint pk : pr.getXIC())
//              if (pk.getY()>0)
//              {
//                // convert to seconds
//                rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));
//              }
//          }
        }
      }

    if (Tools.isSet(protein.getSRMs()))
    {
      protein.composite();
      protein_id.setCompositeSRMGroup(protein);
    }
//    if (Tools.isSet(rt_ai))
//    {
//      double n = 0;
//      for (Float rt : rt_ai.keySet())
//        if (rt_ai.get(rt).size()>n) n = rt_ai.get(rt).size();
//
//      // create the composite trace
//      SRM profile = new SRM(0f, 0f);
//      for (Float rt : rt_ai.keySet())
//      {
//        // look for the ms1 trace
//        float v = (float )Math.pow(10d, Stats.sum(rt_ai.get(rt))/n);
//        profile.addXIC(rt, v);
//      }
//      protein.getSRMs().put(0f, profile);
//      protein_id.setCompositeSRMGroup(protein);
//    }
    return protein_id;
  }
  public static List<Peak> extractRtCal(MultiTreeTable<Float, Float, SRMGroup> landmarks,
                                        MultiTreeTable<Float, Float, SRMGroup> lib,
                                        boolean use_iRT, float min_peak_exclusivity, float min_apex)
  {
    if (!Tools.isSet(landmarks) || !Tools.isSet(lib)) return null;

//    System.out.println("Protein\tSequence\tRT\tfeature.rt\tPkEx\tPkExFront\tfeature.apex");

    Multimap<String, SRMGroup> pep_lib = HashMultimap.create();
    for (SRMGroup grp : lib.values())
      pep_lib.put(grp.getSequence(), grp);

    List<Peak> cals = new ArrayList<>();
    for (SRMGroup grp : landmarks.values())
    {
      SRM cpo = grp.getSRMs().get(0f);
      Collection<SRMGroup> libs = pep_lib.get(grp.getSequence());
//      if (cpo!=null && Tools.isSet(libs))
//      {
//        System.out.println(grp.getProteinId()+"\t"+grp.getSequence()+"\t"+ grp.getRT()+"\t"+cpo.getFeature().getRT() + "\t" + cpo.getPeakPct() + "\t" + cpo.getPeakPctFull()+"\t"+cpo.getFeature().getIntensity());
//      }
      if (cpo!=null && Tools.isSet(libs) && cpo.getPeakPct()>=min_peak_exclusivity && cpo.getFeature().getIntensity()>min_apex)
      {
        double iRTs=0, reportedRT=0;
        for (SRMGroup g : libs) { iRTs+=g.getIRT(); reportedRT+=g.getReportedRT(); }
        iRTs/=libs.size(); reportedRT/=libs.size();

        cals.add(new Peak((use_iRT?iRTs:reportedRT), cpo.getFeature().getRT()));
      }
    }
    pep_lib = Tools.dispose(pep_lib);

    Collections.sort(cals);
    return cals;
  }
  public static MultiTreeTable<Float, Float, SRMGroup> interpolate(MultiTreeTable<Float, Float, SRMGroup> lib, List<Peak> cals, boolean use_iRT)
  {
    MultiTreeTable<Float, Float, SRMGroup> mapped = MultiTreeTable.create();
    for (SRMGroup grp : lib.values())
    {
//      if (grp.getSequence().indexOf("AAEAAINILK")>=0)
//        System.out.println();
//      float rt0 = (use_iRT?grp.getIRT():grp.getReportedRT());

      grp.setRT(Peaks.interpolateForY(cals, (use_iRT?grp.getIRT():grp.getReportedRT())));
      if (!Double.isInfinite(grp.getRT())) mapped.put(grp.getMz(), grp.getRT(), grp);

//      Peak pk   = new Peak(rt0, 0d);
//      int index = Collections.binarySearch(cals, pk), left, right;
//      if (index >= 0)
//      {
//        left = index; right = index;
//      }
//      else  // (-(insertion point) - 1)
//      {
//        index = -1 * index - 1;
//        left = (index > 0 ? index-1 : -1);
//        right = (index < cals.size() ? index : -1);
//      }
//      if (left>=0 && right>=left && right<cals.size())
//      {
//        if (left==right) grp.setRT((float )cals.get(left).getIntensity());
//        else
//        {
//          grp.setRT(Peaks.interpolateForY(cals.get(left), cals.get(right), (double )rt0));
//        }
//
//        // add the transition with the mapped RT
//        mapped.put(grp.getMz(), grp.getRT(), grp);
//      } else
//      {
//        grp.setRT(0f);
//      }
    }
    // remove the old one to save the memory
    lib = Tools.dispose(lib);

    return mapped;
  }
//  public static MultiTreeTable<Float, Float, SRMGroup> fit(MultiTreeTable<Float, Float, SRMGroup> lib, List<Peak> cals, boolean use_iRT)
//  {
//    MultiTreeTable<Float, Float, SRMGroup> mapped = MultiTreeTable.create();
//    for (SRMGroup grp : lib.values())
//    {
////      if (grp.getSequence().indexOf("AAEAAINILK")>=0)
////        System.out.println();
//      float rt0 = (use_iRT?grp.getIRT():grp.getReportedRT());
//
//      Peak pk   = new Peak(rt0, 0d);
//      int index = Collections.binarySearch(cals, pk), left, right;
//      if (index >= 0)
//      {
//        left = index; right = index;
//      }
//      else  // (-(insertion point) - 1)
//      {
//        index = -1 * index - 1;
//        left = (index > 0 ? index-1 : -1);
//        right = (index < cals.size() ? index : -1);
//      }
//      if (left>=0 && right>=left && right<cals.size())
//      {
//        if (left==right) grp.setRT((float )cals.get(left).getIntensity());
//        else
//        {
//          grp.setRT(Peaks.interpolateForY(cals.get(left), cals.get(right), (double )rt0));
//        }
//
//        // add the transition with the mapped RT
//        mapped.put(grp.getMz(), grp.getRT(), grp);
//      } else
//      {
//        grp.setRT(0f);
//      }
//    }
//    // remove the old one to save the memory
//    lib = Tools.dispose(lib);
//
//    return mapped;
//  }
  // write the updated transitions as a new library
  public static void writeTransitions(Collection<SRMGroup> groups, String filename)
  {
    try
    {
      FileWriter tr = new FileWriter(filename);

      tr.write("ProteinId\tModifiedPeptideSequence\tNormalizedRetentionTime\tPrecursorMz\tPrecursorCharge\tiRT\tProductMz\tLibraryIntensity\tReportedRT\tReportedApex\n");
      for (SRMGroup grp : groups)
        for (SRM srm : grp.getSRMs().values())
        {
          // "","PrecursorMz","ProductMz","LibraryIntensity","ProteinId","PeptideSequence","ModifiedPeptideSequence","PrecursorCharge","ProductCharge","FragmentType","FragmentSeriesNumber","NormalizedRetentionTime"
          // "2261",1044.48640687972,405.176849365234,33.29616,"P62258","AAFDDAIAELDTLSEESYK","AAFDDAIAELDTLSEESYK",2,1,"b",4,92.1534957885742
          tr.write(grp.getProteinId()        +"\t"); // ProteinId
          tr.write(grp.getSequence()         +"\t"); // ModifiedPeptideSequence
          tr.write(grp.getRT()               +"\t"); // NormalizedRetentionTime
          tr.write(grp.getMz()               +"\t"); // PrecursorMz
          tr.write(grp.getCharge()           +"\t"); // PrecursorCharge
          tr.write(grp.getIRT()              +"\t"); // iRT
          tr.write(srm.getFragmentMz()       +"\t"); // ProductMz
          tr.write(srm.getLibraryIntensity() +"\t"); // LibraryIntensity
          tr.write(grp.getRT()               +"\t"); // ReportedRT
          tr.write(srm.getApex()             +"\n"); // ReportedRT
        }

      tr.close();
    }
    catch (IOException e) {}
  }
}
