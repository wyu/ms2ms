package org.ms2ms.data.ms;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.jgrapht.Graph;
import org.jgrapht.alg.connectivity.KosarajuStrongConnectivityInspector;
import org.jgrapht.alg.interfaces.StrongConnectivityAlgorithm;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;
import org.ms2ms.Disposable;
import org.ms2ms.data.Point;
import org.ms2ms.math.Points;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;
import toools.collections.Lists;

import java.io.IOException;
import java.io.Writer;
import java.util.*;

public class SRM implements Cloneable, Disposable, Comparable<SRM>
{
  private float mFragmentMz, mLibraryIntensity, mApex, mArea, mPkPct=0, mPkPctFull =0, mPkPctAll=0, mFillTime=0;
  Range<Double> mPeakBoundary;

  private List<LcMsPoint> mXIC;
  private LcMsFeature mFeature;

  public SRM()
  {
    super();
    mFragmentMz=mLibraryIntensity=0;
    mXIC = new ArrayList<>();
  }
  public SRM(float frag, float ai)
  {
    mFragmentMz=frag; mLibraryIntensity=ai;
    mXIC = new ArrayList<>();
  }

  public float getFragmentMz() { return mFragmentMz; }
  public float getLibraryIntensity() { return mLibraryIntensity; }
  public float getApex()       { return mApex; }
  public float getArea()       { return mArea; }
  public float getPeakPct()    { return mPkPct; }
  public float getPeakPctFull() { return mPkPctFull; }
  public float getPeakPctAll() { return mPkPctAll; }
  public float getFillTime()   { return mFillTime; }
  public Range<Double> getPeakBoundary() { return mPeakBoundary; }

  public List<LcMsPoint> getXIC() { return mXIC; }
  public LcMsFeature getFeature() { return mFeature; }
  public LcMsPoint get(int i) { return mXIC.get(i); }

  public LcMsPoint addXIC(float rt, float ai)
  {
    if (ai>0) mXIC.add(new LcMsPoint(rt,ai));
    return Tools.isSet(mXIC)?mXIC.get(mXIC.size()-1):null;
  }
  public LcMsPoint addXIC(float rt, float ai, float mz, int scan, double ppm)
  {
    if (mXIC==null) mXIC = new ArrayList<>();
    mXIC.add(new LcMsPoint(rt,ai,mz,scan, (float )ppm)); return mXIC.get(mXIC.size()-1);
  }
  public LcMsPoint addXIC(double rt, double ai, double mz, int scan)
  {
    return addXIC((float )rt, (float )ai, (float )mz, scan, Float.NaN);
  }

  public SRM setFeature(Point s) { if (s!=null) mFeature = new LcMsFeature(s); return this; }
  public SRM setFeature(LcMsFeature s) { mFeature=s; return this; }
  public SRM setFillTime(float s) { mFillTime=s; return this; }

  public SRM calPeakPct(double rt, double span, int apex_pts, double peak_base)
  {
    if (!Tools.isSet(getXIC())) return this;

    // change the definition on May 16, 2020. outside is now 3x of the LC peak span, instead of 2x
    double inside=0, outside=0, inside_front=0, outside_front=0, all=0;
    Range<Double> inner = Range.closed(rt-span, rt+span), outer = Range.closed(rt-3d*span, rt+3d*span),
            inner_front = Range.closed(rt-span, rt),        outer_front = Range.closed(rt-3d*span, rt);

    List<Double> tops = new ArrayList<>();
    int apex_i=-1; double best=Float.MAX_VALUE;
    for (int i=0; i<getXIC().size(); i++)
    {
      Point p = get(i);
      if (inner.contains(p.getX()))
      {
        inside+=p.getY();
        if (p.getY()>0) tops.add(p.getY());
      }
      if (inner_front.contains(p.getX())) inside_front+=p.getY();
      // change the definition on May 16, 2020. outside is sum of intensity within 'outside', instead of intensities outside of 'outside'
      if (outer.contains(      p.getX()))  outside      +=p.getY();
      if (outer_front.contains(p.getX()))  outside_front+=p.getY();

      all += p.getY();

      if (Math.abs(p.getX()-rt)<best) { apex_i=i; best=Math.abs(p.getX()-rt); }
    }
    if (getFeature()!=null && Tools.isSet(tops)) {
      Collections.sort(tops, Ordering.natural().reversed());
      // OK with SRM with fewer points
      getFeature().setApex(Stats.mean(tops.subList(0, Math.min(tops.size(), apex_pts))));
      mApex = (float )getFeature().getApex();
    }
    if (apex_i>0 && getFeature()!=null && getFeature().getApex()>0)
    {
      double cut = getFeature().getApex()*peak_base, left=0, right=0, area=0;
      for (int i=apex_i; i>0; i--)
        if (get(i).getIntensity()<cut && (i<=0 || get(i-1).getIntensity()<cut)) { left=get(i).getRT(); break; } else area+=get(i).getIntensity();

      for (int i=apex_i+1; i<getXIC().size()-1; i++)
        if (get(i).getIntensity()<cut && (i>=getXIC().size() || get(i+1).getIntensity()<cut)) { right=get(i).getRT(); break; } else area+=get(i).getIntensity();

      if (left>0 && right>left) { mPeakBoundary = Range.closed(left, right); mArea=(float )area; getFeature().setArea(area); }
    }
    // change the definition on May 16, 2020. exclusivity is now inside/outside, for local exclusivity, instead of inside_intensity / (intensity outside of 2x span)
    mPkPctFull = (float )(100f*inside/outside); // a local exclusivity (+-3x LC peak span)
    mPkPct     = (float )(100f*inside_front/outside_front); // a local exclusivity (+-3x LC peak span)
    mPkPctAll  = (float )(100f*mArea/all);      // this is now a global (+-5min) exclusivity
    tops = (List )Tools.dispose(tops);

    return this;
  }
  // center = the targeted RT,
  // span   = the span of the RT window where the peaks are expected
  public List<Peak> detectPeak(float center, float span)
  {
    double base_deri = 0d;
    for (Point p : getXIC())
      if (Math.abs(p.getX()-center)<=span && p.getY()>base_deri) { base_deri=p.getY(); }

    // require the peak to be at least 10% of the base
    base_deri /= 10d;

    List<Point> deri = Points.deriv1stBySG5(getXIC());
    List<Peak> pks = new ArrayList<>();

    double deri_max=0;
    if (Tools.isSet(deri))
      for (int i=0; i<deri.size()-1; i++)
      {
        if (deri.get(i).getY()>deri_max) deri_max=deri.get(i).getY();

//        if (deri.get(i+1).getY()>=0 && deri.get(i  ).getY()>deri.get(i+1).getY() &&
//            deri.get(i+2).getY()<=0 && deri.get(i+2).getY()>deri.get(i+3).getY())
        if (deri.get(i).getY()>=0 && deri.get(i+1).getY()<=0)
        {
          // the peak top is at the zero-intercept
          Point top = Points.interpolateByY(deri.get(i), deri.get(i+1), 0d);
          if (deri_max>base_deri)
          {
            pks.add(new Peak(top.getX(), deri_max)); deri_max=0;
          }
        }
      }

    deri = (List )Tools.dispose(deri);
    return pks;
 }
  public SRM clone()
  {
    SRM cloned      = new SRM(getFragmentMz(), getLibraryIntensity());
    cloned.mApex    = mApex; cloned.mArea = mArea; cloned.mPkPct = mPkPct; cloned.mPkPctAll = mPkPctAll;
    cloned.mFeature = mFeature;

    cloned.mXIC     = new ArrayList<>();
    for (LcMsPoint p : mXIC) cloned.mXIC.add(p);

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
  public SRM fill(float baseline, Multimap<Float, Float> xs)
  {
    if (Tools.isSet(xs))
      for (Float x : xs.get(0f))
        // need to bypass the non-zero check in addXIC
        if (!xs.get(getFragmentMz()).contains(x)) mXIC.add(new LcMsPoint(x,baseline));

    Collections.sort(getXIC());

    return this;
  }
  public SRM impute(float gap)
  {
//      if (Math.abs(mFragmentMz-445.2769)<0.01)
//        System.out.println();
    int last=-1, ahead=mXIC.size();
    if (Tools.isSet(mXIC) && mXIC.size()>5)
      for (int i=0; i<mXIC.size(); i++)
      {
        if (get(i).getIntensity()>0) last=i;
        else if (last>=0 && get(i).getRT()-get(last).getRT()<=gap)
        {
          // a missing point
          for (int j=i+1; j<mXIC.size(); j++)
            if (mXIC.get(j).getY()>0) { ahead=j; break; }

          if (ahead<mXIC.size() && get(ahead).getRT()-get(i).getRT()<=gap)
          {
            if (ahead==last) get(i).setIntensity(get(last).getIntensity()/2d);
            else             get(i).setIntensity(get(last).getIntensity()+(get(i).getRT()-get(last).getRT()) * (get(ahead).getIntensity()-get(last).getIntensity()) / (get(ahead).getRT()-get(last).getRT()));
            get(i).isImputed(true);
          }
        }
      }

    return this;
  }
  // shift the RT axis to +- around the center, then interpolate them to the new RT points
  public SRM shift(Range<Float> xs, int steps, Float center)
  {
    if (center!=null)
      for (LcMsPoint p : getXIC()) p.setRT(p.getRT()-center);

    if (Tools.isSet(xs))
    {
      List<LcMsPoint> pts = new ArrayList<>(steps);
      float step = (xs.upperEndpoint()-xs.lowerEndpoint())/steps;
      for (double x=xs.lowerEndpoint(); x<=xs.upperEndpoint(); x+=step)
      {
        // interpolate from the existing array
        Point p = Points.interpolate(getXIC(), x, false);
        if (p==null) p = new Point(x, 0);

        pts.add(new LcMsPoint(p));
      }
      mXIC = (List )Tools.dispose(mXIC);
      mXIC = pts;
    }
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
  public String getUID(String s)
  {
    return s+"_"+Tools.d2s(getFragmentMz(), 4);
  }

  public SRM disposeXIC()
  {
    mXIC = (List )Tools.dispose(mXIC);
    return this;
  }
  @Override
  public void dispose()
  {
    disposeXIC();
    mFeature = null;
  }

  @Override
  public int compareTo(SRM o)
  {
    int c = Float.compare(mFragmentMz, o.mFragmentMz);
    if (c==0 && mFeature!=null)
    {
      c = Double.compare(mFeature.getApex(), o.mFeature.getApex());
    }
    return (c);
  }

  public static void headerEdges(Writer w) throws IOException
  {
    w.write("Source\tTarget\tdp\n");
  }
  public static void printEdges(Writer w, SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> net, String tag) throws IOException
  {
    if (net!=null && Tools.isSet(net.edgeSet()))
      for (DefaultWeightedEdge e : net.edgeSet())
      {
        w.write(net.getEdgeSource(e).getUID(tag)+"\t");
        w.write(net.getEdgeTarget(e).getUID(tag)+"\t");
        w.write(Tools.d2s(net.getEdgeWeight(e),2)+"\n");
      }
  }
  public static void headerNodes(Writer w) throws IOException
  {
    w.write("UID\tFragMz\n");
  }
  public static void printNodes(Writer w, SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> net, String tag) throws IOException
  {
    if (net!=null && Tools.isSet(net.vertexSet()))
      for (SRM v : net.vertexSet())
      {
        w.write(v.getUID(tag)+"\t");
        w.write(Tools.d2s(v.getFragmentMz(),2)+"\n");
      }
  }
  public static Map<String, Object> inspectNetwork(SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> net)
  {
    StrongConnectivityAlgorithm<String, DefaultWeightedEdge> scAlg =
        new KosarajuStrongConnectivityInspector(net);

    // a graph is said to be strongly connected if every vertex is reachable from every other vertex
    List<Graph<String, DefaultWeightedEdge>> stronglyConnectedSubgraphs = scAlg.getStronglyConnectedComponents();

    List<Integer> sc = new ArrayList<>(stronglyConnectedSubgraphs.size());

    // prints the strongly connected components
    for (int i = 0; i < stronglyConnectedSubgraphs.size(); i++)
      if (stronglyConnectedSubgraphs.get(i).vertexSet().size()>1)
        sc.add(stronglyConnectedSubgraphs.get(i).vertexSet().size());

    Collections.sort(sc, Ordering.natural().reversed());
    if (sc.size()>2)
      System.out.print("");

    Map<String, Object> props = new HashMap<>();
    props.put("SCC: length", stronglyConnectedSubgraphs.size());
    props.put("SCC: node size", sc.toArray(new Integer[sc.size()]));
    props.put("Network: node size", net.vertexSet().size());
    props.put("Network: edge size", net.edgeSet().size());

//    BronKerboschCliqueFinder<SRM, DefaultWeightedEdge> clique = new BronKerboschCliqueFinder(net);
    return props;
  }
}

