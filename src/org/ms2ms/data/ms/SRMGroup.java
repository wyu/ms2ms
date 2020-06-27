package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;
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

import static org.ms2ms.math.Stats.closed;

public class SRMGroup implements Ion, Comparable<SRMGroup>, Cloneable
{
  private String mPeptideSequence, mProteinId;
  private float mRT, mRtOffset, mPrecursorMz, mDpSimilarity, mIRT, mReportedRT, mNetworkNiche;
  private int mCharge, mQualifiedSRMs;
  private Map<String, Object> mNetworkStats;

  private LcMsFeature mCompositePeak=null, mCompositeFeature=null;

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
    mRtOffset=0;
  }
  public SRMGroup(String peptide, float rt, float mz, int z)
  {
    super();
    mPeptideSequence = peptide;
    mSRMs            = new TreeMap<>();

    mRT=rt; mPrecursorMz=mz; mCharge=z; mRtOffset=0;
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
  public float  getRtOffset()     { return mRtOffset; }
  public float  getReportedRT()   { return mReportedRT; }
  public float  getSimilarity()   { return mDpSimilarity; }
  public String getSequence()     { return mPeptideSequence; }
  public String getProteinId()    { return mProteinId; }
  public SRM getComposite() { return mSRMs!=null?mSRMs.get(0f):null; }

  public SimpleDirectedWeightedGraph<SRM, DefaultWeightedEdge> getNetwork() { return mNetwork; }

  public SRM getCompositeProfile(boolean round2sec)
  {
    SRM cpo = getComposite();

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
  public SRMGroup setRtOffset(float s) { mRtOffset=s; return this; }
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
      return mSRMs.get(frag).addXIC(rt, intensity!=null? (float )intensity.doubleValue():0f,
          mz!=null? (float)mz.doubleValue():0f, scan, 1E6*(mz-frag)/frag).setFillTime(fill);
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
      if (frag>0 && mSRMs.get(frag)!=null && mSRMs.get(frag).getXIC()!=null) // only the MS2 XIC
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
    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create(), rt_ai0 = TreeMultimap.create(), rt_pm = TreeMultimap.create();
    HashMap<Float, Integer>     rt_sn = new HashMap<>();

    for (Float frag : mSRMs.keySet())
      if (frag>0 && mSRMs.get(frag)!=null && mSRMs.get(frag).getXIC()!=null) // only the MS2 XIC
        for (LcMsPoint pk : mSRMs.get(frag).getXIC())
        {
          if (pk.getY()>0)
          {
            rt_ai.put((float )pk.getX(), Math.log10(pk.getY()));
            rt_sn.put((float )pk.getX(), pk.getScan());
            if (pk.getPPM()!=Float.NaN) rt_pm.put((float )pk.getX(), pk.getPPM());
          }
          rt_ai0.put((float )pk.getX(), pk.getY());
        }

    // create the composite trace
    double n = mSRMs.size(); int ms1_start=0;
    SRM ms1 = mSRMs.get(-1f);
    mSRMs.put(0f, new SRM(0f, 0f));
    for (Float rt : rt_ai0.keySet())
    {
      float v = 0f; Integer scan = 0; Double ppm=0d;
      if (rt_ai.containsKey(rt))
      {
        // look for the ms1 trace
        if (ms1!=null && Tools.isSet(ms1.getXIC()))
          for (int i=ms1_start; i<ms1.getXIC().size()-1; i++)
            if (rt>=ms1.getXIC().get(i).getX() && rt<=ms1.getXIC().get(i+1).getX())
            {
              rt_ai.put(rt, Math.log10(0.5*(ms1.getXIC().get(i).getY()+ms1.getXIC().get(i+1).getY()))); ms1_start=i+1; break;
            }

        v    = (float )Math.pow(10d, Stats.sum(rt_ai.get(rt))/n);
        scan = rt_sn.get(rt)*-1;
        ppm  = Stats.sum(rt_pm.get(rt))/n;
      }
      mSRMs.get(0f).addXIC(rt, v, 0f, scan, ppm);
    }
    rt_ai = (TreeMultimap )Tools.dispose(rt_ai);
    rt_ai0= (TreeMultimap )Tools.dispose(rt_ai0);
    rt_pm = (TreeMultimap )Tools.dispose(rt_pm);
    rt_sn = (HashMap )Tools.dispose(rt_sn);

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

    // need to add the missing points for very sparse data
    if (Tools.isSet(Xs.get(0f)) && Xs.get(0f).size()>10)
    {
      List<Float> xx = new ArrayList<>(Xs.get(0f)), dx = new ArrayList<>();

      Collections.sort(xx);
      for (int i=1; i<xx.size(); i++) dx.add(xx.get(i)-xx.get(i-1));
      // sort the steps
      Collections.sort(dx);
      float step = Stats.meanFloats(dx.subList(0, 5));
      for (int i=1; i<xx.size(); i++)
        if (xx.get(i)-xx.get(i-1)>step*1.5)
          for (float x=xx.get(i-1)+step; xx.get(i)-x>step/2; x+=step)
            Xs.get(0f).add(x);

      xx = (List )Tools.dispose(xx);
      dx = (List )Tools.dispose(dx);
    }

    for (SRM srm : mSRMs.values()) srm.fill(baseline, Xs);

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
    return centroid(settings.getPeakWidth(), settings.getQuanSpan(), settings.getQuanOffset(), settings.getMinSNR(), settings.getBaseRI());
  }
  public SRMGroup centroid(float rt_width, float quan_span, float quant_offset, float minSNR, float min_ri)
  {
    // nothing to do without a composite!
    if (getComposite()==null || !Tools.isSet(getComposite().getXIC())) return this;

//    if (Strs.isA(getSequence(), "WPPEDEISKPEVPEDVDLDLKK#4","LMYEHELEQLR#2","VLVLGDSGVGK#2","FEFFEGLLFEGR#2"))
//    if (Strs.isA(getSequence(), "FLNSFASMHR#3","TNPPLIQEKPAK#3","TQSPC[Carbamidomethyl (C)]FGDDDPAKK#3"))
//    if (Strs.isA(getSequence(), "EFSHIAFLTIK#3","DLFDSM[Oxidation (M)]DK#2","VADEGSFTC[Carbamidomethyl (C)]FVSIR#2","DGLLPTGLGQR#2","LIPGC[Carbamidomethyl (C)]EVILATPYGR#2"))
    if (Strs.isA(getSequence(), "﻿NSLIQMTILNLLPR#2","HSTSGTDEGEDGDEPDDGSNDVVDLLPR#3","M[Oxidation (M)]GPGAASGGERPNLK#2","SYNPFDDDGEDEGAR#2","M[Oxidation (M)]QHNLEQQIQAR#2","LADVDKDGLLDDEEFALANHLIK#3"))
      System.out.print("");
    // check the alternative peaks by 1st derivatives
    getComposite().detectPeak(getRT(), quan_span);
    if (mSRMs.get(-1f)!=null) mSRMs.get(-1f).detectPeak(getRT(), quan_span);

    double outer_span   = ((Tools.back(getComposite().getXIC()).getRT()-Tools.front(getComposite().getXIC()).getRT())-rt_width)/2f;

    if (outer_span<=0) return this;

    double rt0 = getRT()+quant_offset;
    Range<Double> inner = Range.closed(rt0-quan_span,  rt0+quan_span),
                 suburb = Range.closed(rt0-quan_span*1.5,  rt0+quan_span*1.5),
                 outer  = Range.closed(rt0-outer_span, rt0+outer_span);

    // setup the MS1 list
    SortedMap<Double, LcMsFeature> rt_ms1 = new TreeMap<>();
    if (getSRM(-1f)!=null && Tools.isSet(getSRM(-1f).getPeaks()))
      for (LcMsFeature p : getSRM(-1f).getPeaks())
        if (inner.contains(p.getX())) rt_ms1.put(p.getX(), p);

    // determine the best peak from this XIC trace
    LcMsFeature selected=null;
    if (!Tools.isSet(getComposite().getPeaks()))
    {
      // nothing!!
      System.out.print("");
    }
    else if (getComposite().getPeaks().size()==1) selected = getComposite().getPeaks().get(0);
    else
    {
      LcMsFeature highest_unbound = null, highest_burb=null;
      List<LcMsFeature> inbounds = new ArrayList<>(), tops = new ArrayList<>();
      for (LcMsFeature pk : getComposite().getPeaks())
      {
        if (highest_unbound==null && !suburb.contains(pk.getX()) && outer.contains(pk.getX())) highest_unbound = pk;
        if (highest_burb   ==null && !inner.contains(pk.getX()) && suburb.contains(pk.getX())) highest_burb = pk;
        if (inner.contains(pk.getX())) inbounds.add(pk);
      }

      for (int i=0; i<inbounds.size()-1; i++)
      {
        LcMsFeature pk = getComposite().getPeaks().get(i);
        if (inner.contains(pk.getX())) tops.add(pk);
        if (pk.getArea()/getComposite().getPeaks().get(i+1).getArea()>=minSNR) break;
      }
      // a single good hit is a good one!
      if (tops.size()==0 && inbounds.size()>0) tops.add(inbounds.get(0));

      if (Tools.isSet(tops))
      {
        // only check the MS1 if the top two are too close in snr
        if (tops.size()>1 && tops.get(0).getSNR()<tops.get(1).getSNR()*1.5f && Tools.isSet(rt_ms1))
        {
          for (LcMsFeature pk : tops)
            if (hasMS1(rt_ms1, pk.getX(), rt_width/2f, 0f))
            { selected=pk; selected.wasBasedOn(LcMsFeature.criteria.ms1); break; }
        }
        if (selected==null) { selected = tops.get(0); selected.wasBasedOn(LcMsFeature.criteria.snr); }
      }
      // don;t worry about the burb if the unbound one is obviously better
      if (highest_burb!=null && highest_unbound!=null &&
          highest_burb.getSNR()<highest_unbound.getSNR()/10) highest_burb=null;

      // only expand outside of the quan window if supported by the MS1
      if (highest_unbound!=null && !hasMS1(rt_ms1, highest_unbound.getX(), rt_width/2f, 3f)) highest_unbound=null;
      if (highest_burb   !=null && !hasMS1(rt_ms1, highest_burb.getX(),    rt_width/2f, 3f)) highest_burb   =null;

      if (highest_burb!=null && (selected==null || selected.getY()*3f<highest_burb.getY()))
      {
        selected = highest_burb; selected.wasBasedOn(LcMsFeature.criteria.suburb);
      }
      else if (highest_unbound!=null && (selected==null || selected.getY()*10f<highest_unbound.getY()))
      {
        selected = highest_unbound; selected.wasBasedOn(LcMsFeature.criteria.outer);
      }
      // grabt the top one from within the quan span. Have to qualify the feature later on
      if (selected==null && Tools.isSet(tops)) selected = tops.get(0);

      // check for likely problem
      if (Tools.isSet(inbounds) && selected!=inbounds.get(0) &&
          inbounds.get(0).getSNR()>selected.getSNR()*10)
        System.out.print("");

      inbounds = (List )Tools.dispose(inbounds); tops = (List )Tools.dispose(tops);
    }
    if (selected!=null)
    {
      selected.setApexRT(selected.getX()).setAbundance(selected.getY());
      getComposite().setFeature(selected);
      // determine the peak boundary
      getComposite().setPeakBoundary(getComposite().getFeature(), min_ri);
      // propagate to the individual transition using the same peak boundary
      centroidTransitions(getComposite().getFeature(), getComposite().getPeakBoundary());
    }

    return this;
  }
  private boolean hasMS1(SortedMap<Double, LcMsFeature> rt_ms1, double rt, float rt_tol, float snr)
  {
    Map<Double, LcMsFeature> slice = rt_ms1.subMap(rt-rt_tol/2d, rt+rt_tol/2d);
    if (snr>0 && Tools.isSet(slice))
    {
      int n=0;
      for (LcMsFeature F : slice.values()) if (F.getSNR()>=snr) n++;
      return (n>0);
    }
    return Tools.isSet(slice);
  }
  private SRMGroup centroidTransitions(LcMsFeature cpo, Range<Double> rt_range)
  {
    if (cpo==null || !Tools.isSet(rt_range))
      return this;
    // propagate the features to each transition once the composite feature is determined
    List<Point> injects = new ArrayList(), mzs = new ArrayList<>();
    for (Float frag : mSRMs.keySet())
    {
      SRM srm = mSRMs.get(frag);

      srm.updateFeature(Points.centroid(srm.getXIC(), 0d, rt_range));

      if (frag>0)
      {
        injects.clear();  mzs.clear();
        for (LcMsPoint p : srm.getXIC())
          if (p.getIntensity()>0 && rt_range.contains(p.getRT()))
          {
            injects.add(new Point(p.getFillTime(), p.getIntensity()));
            mzs.add(    new Point(p.getMz(), p.getIntensity()));
          }

        if (Tools.isSet(injects) && srm.getFeature()!=null)
          srm.getFeature().setFillTime(Points.centroid(injects));
        if (Tools.isSet(mzs) && srm.getFeature()!=null)
          srm.getFeature().setMz(Points.centroid(mzs), srm.getFragmentMz());
      }
    }
    injects = (List )Tools.dispose(injects); mzs = (List )Tools.dispose(mzs);

    return this;
  }
  public SRMGroup centroidByPeakBoundry(int apex_pts)
  {
//    if (getSequence().equals("TITLEVEPSDTIENVK#2"))
//      System.out.print("");
    Range<Double> boundary = getComposite().getPeakBoundary();

    if (!Tools.isSet(boundary))
      return this;

    List<Point>  pts  = new ArrayList(),   mzs = new ArrayList();
    List<Double> devi = new ArrayList<>(), mex = new ArrayList(), devi0 = new ArrayList<>(), mex0 = new ArrayList();
    Map<Double, Double> rt_ij = new TreeMap<>();
    List<Double> tops = new ArrayList<>();

    for (Float frag : mSRMs.keySet())
    {
      SRM srm = mSRMs.get(frag);
      if (srm.getFeature()==null) srm.setFeature(new LcMsFeature());
      // look for the peak properties
      pts.clear(); mzs.clear(); mex.clear(); devi.clear();
      for (LcMsPoint p : srm.getXIC())
        if (p.getIntensity()>0 && boundary.contains(p.getX()))
        {
          pts.add(new Point(p.getX(),  p.getIntensity()));
          mzs.add(new Point(p.getMz(), p.getIntensity()));
          if (frag==0f) tops.add(p.getIntensity());
          if (frag>0) devi.add(p.getPPM());
          if (frag>0 && p.getFillTime()>0) rt_ij.put(p.getX(), p.getFillTime());
        } else
        {
          // outside of the peak boundary
          if (p.getMz()>0 && frag>0) mex.add(p.getPPM());
        }

      if (Tools.isSet(pts))
      {
        srm.getFeature().setInitialCentroidRt(srm.getFeature().getX());
        srm.getFeature().setX( Points.centroid(pts));
        srm.getFeature().setMz(Points.centroid(mzs), srm.getFragmentMz());
        srm.getFeature().setArea(Points.sumY(pts));
      }
      // on the individual SRM
      srm.getFeature().setMzStdev(Tools.isSet(devi) && devi.size()>2?Stats.stdev(devi):Double.NaN);
      srm.getFeature().setMzStdevEx(Tools.isSet(mex) && mex.size()>2?Stats.stdev(mex ):Double.NaN);
      // add them to the composite
      devi0.addAll(devi); mex0.addAll(mex);
    }
    // set the composite
    if (getComposite()!=null && getComposite().getFeature()!=null)
    {
      getComposite().getFeature().setPPM(Tools.isSet(devi0) && devi0.size()>2?(float )Stats.mean(devi0):Float.NaN);
      getComposite().getFeature().setMzStdev(  Tools.isSet(devi0) && devi0.size()>2?Stats.stdev(devi0):Double.NaN);
      getComposite().getFeature().setMzStdevEx(Tools.isSet(mex0) && mex0.size()>2?Stats.stdev(mex0 ):Double.NaN);

      if (Tools.isSet(tops)) {
        Collections.sort(tops, Ordering.natural().reversed());
        // OK with SRM with fewer points
        getComposite().getFeature().setApex(Stats.mean(tops.subList(0, Math.min(tops.size(), apex_pts))));
//      mApex = (float )getFeature().getApex();
      }


      List<Point>  ijs = new ArrayList<>();
      for (LcMsPoint p : getComposite().getXIC())
        if (rt_ij.containsKey(p.getRT())) ijs.add(new Point(rt_ij.get(p.getRT()), p.getIntensity()));

      getComposite().getFeature().setFillTime(Tools.isSet(ijs) ? Points.centroid(ijs):Double.NaN);
      ijs = (List )Tools.dispose(ijs);
    }
    else
    {
      System.out.print("");
    }

    pts = (List )Tools.dispose(pts);  mzs  = (List )Tools.dispose(mzs);
    devi= (List )Tools.dispose(devi); devi0= (List )Tools.dispose(devi0);
    mex = (List )Tools.dispose(mex);  mex0 = (List )Tools.dispose(mex0);
    rt_ij = (Map )Tools.dispose(rt_ij);

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
    return calcFeatureExclusivity(settings.getQuanSpan(), settings.getApexPts());
  }
  public SRMGroup calcFeatureExclusivity(float quan_span, int apex_pts)
  {
    if (Tools.isSet(mSRMs) && mSRMs.get(0f)!=null && mSRMs.get(0f).getFeature()!=null)
    {
//      if (getSequence().equals("NFDVGHVPIR#2"))
//        System.out.print("");

      double rt = mSRMs.get(0f).getFeature().getX();
      for (SRM srm : mSRMs.values())
        srm.calPeakExclusivity(rt, quan_span, getComposite().getPeakBoundary());
    }
    return this;
  }
//  public SRMGroup setPeakBoundary(double peak_base)
//  {
//    if (Tools.isSet(mSRMs) && getComposite()!=null && getComposite().getFeature()!=null)
//      for (SRM srm : mSRMs.values())
//        srm.setPeakBoundary(srm.getFeature(), peak_base);
//
//    return this;
//  }
  public SRMGroup scanMS2(SortedMap<Double, Peak> peaks, float rt, int scan, Double fill_time, Tolerance tol, boolean keep_zero)
  {
    for (Float k : mSRMs.keySet())
    {
      if (k<=0) continue;

      SortedMap<Double, Peak> pks = peaks.subMap(tol.getMin(k), tol.getMax(k));
      if (pks!=null && !pks.isEmpty())
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
//    headerGroup(w);
    w.write("Peptide\tz\tPrecMz\tRT\tFragMz\tscan\txic.rt\txic.ai\tinjection\timputed\txic.ppm\n");
  }
  public static void header_xic(Writer w) throws IOException
  {
    // xic.ai~xic.rt|FragMz
    w.write("Peptide\tFragMz\txic.rt\txic.ai\txic.ppm\n");
  }
  public SRMGroup print_xic(Writer w, boolean keep_zero, float... frags) throws IOException
  {
    if (!Tools.isSet(frags)) frags = Tools.toFloatArray(mSRMs.keySet());
    for (Float frag : frags)
      if (mSRMs!=null && mSRMs.get(frag)!=null && mSRMs.get(frag).getXIC()!=null)
        for (LcMsPoint pk : mSRMs.get(frag).getXIC())
          if (keep_zero || pk.getIntensity()>0)
          {
            w.write(getSequence()+"\t");
            w.write(frag             +"\t");
            w.write(pk.getRT()       +"\t");
            w.write(pk.getIntensity()+"\t");
            w.write(pk.getPPM()      +"\n");
          }

    return this;
  }

  public SRMGroup printXIC(Writer w, boolean keep_zero, float... frags) throws IOException
  {
    if (!Tools.isSet(mSRMs)) return this;

    if (!Tools.isSet(frags)) frags = Tools.toFloatArray(mSRMs.keySet());
    for (Float frag : frags)
      if (mSRMs.get(frag)!=null && Tools.isSet(mSRMs.get(frag).getXIC()))
        for (LcMsPoint pk : mSRMs.get(frag).getXIC())
          if (keep_zero || pk.getIntensity()>0)
          {
            w.write(getSequence()+"\t");
            w.write(getCharge()+"\t");
            w.write(Tools.d2s(getMz(),4)+"\t");
            w.write(Tools.d2s(getRT(),3)+"\t");
            w.write(Tools.d2s(frag,4)             +"\t");
            w.write(Tools.d2s(pk.getScan(),2)     +"\t");
            w.write(Tools.d2s(pk.getRT(),3)       +"\t");
            w.write(Tools.d2s(pk.getIntensity(),2)+"\t");
            w.write(Tools.d2s(pk.getFillTime(),2)+"\t");
            w.write(pk.isImputed()                   +"\t");
            w.write(Tools.d2s(pk.getPPM(),2)      +"\n");
          }

    return this;
  }
  public static void headerFeatures(Writer w) throws IOException
  {
    headerGroup(w);
    w.write("FragMz\tNumPts\tPkEx\tPkExAll\tfeature.rt\tfeature.ai\tinjection\tfeature.apex\tfeature.area\tfeature.ppm\tfeature.snr\tppm.stdev\tppm.stdev.ex\trule\tinbound.snr\tlower\tupper\n");
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
        w.write(Tools.d2s(srm.getPeakPct(),2)+"\t");
        w.write(Tools.d2s(srm.getPeakPctAll(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getX(),3)+"\t");
        w.write(Tools.d2s(srm.getFeature().getY(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getFillTime(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getApex(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getArea(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getPPM(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getSNR(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getMzStdev(),2)+"\t");
        w.write(Tools.d2s(srm.getFeature().getMzStdevEx(),2)+"\t");
        w.write(srm.getFeature().wasBasedOn()+"\t");
        w.write((Tools.isSet(srm.getPeaks())?Tools.d2s(srm.getPeaks().get(0).getSNR(),2):"")+"\t");
        w.write(Tools.d2s(srm.getPeakBoundary()!=null?srm.getPeakBoundary().lowerEndpoint():0,2)+"\t");
        w.write(Tools.d2s(srm.getPeakBoundary()!=null?srm.getPeakBoundary().upperEndpoint():0,2)+"\n");
      }
    return this;
  }
  public static void headerGroup(Writer w) throws IOException
  {
    w.write("Peptide\tz\tPrecMz\tRT\toffset\tSimilarity\tYield\tNetV\tNetE\tSCC\tSCC1\t");
  }
  private void printGroup(Writer w) throws IOException
  {
    w.write(getSequence()+"\t");
    w.write(getCharge()+"\t");
    w.write(Tools.d2s(getMz(),4)+"\t");
    w.write(Tools.d2s(getRT(),3)+"\t");
    w.write(Tools.d2s(getRtOffset(),2)+"\t");
    w.write(Tools.d2s(getSimilarity(),3)+"\t");
    w.write(Tools.d2s(calcSrmYield(0.5f), 2)+"\t");
    w.write((mNetwork!=null?mNetwork.vertexSet().size():0)+"\t");
    w.write((mNetwork!=null?mNetwork.edgeSet().size():0)+"\t");
    w.write((mNetworkStats!=null?((int )mNetworkStats.get("SCC: length")):0)+"\t");
    w.write((mNetworkStats!=null && mNetworkStats.get("SCC: node size")!=null? ((Integer[] )mNetworkStats.get("SCC: node size"))[0]:0)+"\t");
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

    cloned.setRtOffset(getRtOffset());
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
            System.out.print("!");
          }
        }
        group.addTransition(tr.getFloat("ProductMz"), tr.get("LibraryIntensity", 0f));
      }
      tr.close();
      System.out.println();

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
//    TreeMultimap<Float, Double> rt_ai = TreeMultimap.create();
    Range<Float> window = Range.closed(settings.getSpan()*-1, settings.getSpan());
    for (SRMGroup group : groups)
      if (Strs.isA(group.getSequence(), peptides))
      {
        protein_id.addSRMGroup(group, group.getSequence());
        if (group.getComposite()!=null && Tools.isSet(group.getComposite().getXIC()) && group.getComposite().getFeature()!=null)
        {
          protein.getSRMs().put(group.getMz(), group.getComposite().shift(window, settings.getGridSize(),
              (float )group.getComposite().getFeature().getRT()));
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
                                        boolean use_iRT, float min_peak_exclusivity, float min_apex, Writer w) throws IOException
  {
    if (!Tools.isSet(landmarks) || !Tools.isSet(lib)) return null;

//    System.out.println("Protein\tSequence\tRT\tfeature.rt\tPkEx\tPkExFront\tfeature.apex");
    if (w!=null) w.write("Protein\tPeptide\trtA\trtB\n");

    Multimap<String, SRMGroup> pep_lib = HashMultimap.create();
    for (SRMGroup grp : lib.values())
      pep_lib.put(grp.getSequence(), grp);

    List<Peak> cals = new ArrayList<>();
    for (SRMGroup grp : landmarks.values())
    {
      SRM cpo = grp.getSRMs().get(0f);
      Collection<SRMGroup> libs = pep_lib.get(grp.getSequence());
      if (cpo!=null && Tools.isSet(libs) &&
          // need both the local and global check to ensure the quality of the calibrants. WYU 20200625
          cpo.getPeakPct()>=min_peak_exclusivity && cpo.getPeakPctAll()>=min_peak_exclusivity &&
          cpo.getFeature().getIntensity()>min_apex)
      {
        double iRTs=0, reportedRT=0;
        for (SRMGroup g : libs) { iRTs+=g.getIRT(); reportedRT+=g.getReportedRT(); }
        iRTs/=libs.size(); reportedRT/=libs.size();

        Peak p = new Peak((use_iRT?iRTs:reportedRT), cpo.getFeature().getRT());
        cals.add(p);
        if (w!=null)
        {
          w.write(grp.getProteinId()+"\t"+grp.getSequence()+"\t" + p.getMz()+"\t"+p.getIntensity()+"\n");
        }
      }
    }
    pep_lib = Tools.dispose(pep_lib);
    if (w!=null) w.close();

    Collections.sort(cals);
    return cals;
  }
  public static MultiTreeTable<Float, Float, SRMGroup> interpolate(MultiTreeTable<Float, Float, SRMGroup> lib, List<Peak> cals, boolean use_iRT)
  {
    MultiTreeTable<Float, Float, SRMGroup> mapped = MultiTreeTable.create();
    for (SRMGroup grp : lib.values())
    {
      grp.setRT(Peaks.estimateForY(cals, (use_iRT?grp.getIRT():grp.getReportedRT())));
      if (!Double.isInfinite(grp.getRT())) mapped.put(grp.getMz(), grp.getRT(), grp);
    }
    // remove the old one to save the memory
    lib = Tools.dispose(lib);

    return mapped;
  }
  public static MultiTreeTable<Float, Float, SRMGroup> calibrate(LcSettings setting, MultiTreeTable<Float, Float, SRMGroup> lib, List<Peak> cals, boolean use_iRT)
  {
    MultiTreeTable<Float, Float, SRMGroup> mapped = MultiTreeTable.create();
    if (setting.isCalMethod(LcSettings.calibration.pt2, LcSettings.calibration.SG5))
    {
      // smooth the dRT/RT curve if asked
      if (setting.toSmoothRT()) cals = Peaks.smoothBySG5(cals);

      for (SRMGroup grp : lib.values())
      {
        grp.setRT(Peaks.estimateForY(cals, (use_iRT?grp.getIRT():grp.getReportedRT())));
        if (!Double.isInfinite(grp.getRT())) mapped.put(grp.getMz(), grp.getRT(), grp);
      }
    } else if (setting.isCalMethod(LcSettings.calibration.loess))
    {
      double[] xs = new double[cals.size()], ys = new double[cals.size()];
      for (int i=0; i<cals.size(); i++) { xs[i] = cals.get(i).getMz(); ys[i] = cals.get(i).getIntensity(); }

      try
      {
        PolynomialSplineFunction poly = new LoessInterpolator(setting.getBandwidth(), 2).interpolate(xs, ys);
        // compute the interpolated value
        Range<Double> bound = closed(xs);
        for (SRMGroup grp : lib.values())
        {
          double x = use_iRT?grp.getIRT():grp.getReportedRT();
          grp.setRT(bound.contains(x)?poly.value(x):0);
          if (grp.getRT()>0) mapped.put(grp.getMz(), grp.getRT(), grp);
        }
      }
      catch (MathException e)
      {
        throw new RuntimeException("Something is wrong with loess regression!", e);
      }
    }
    // remove the old one to save the memory
    lib = Tools.dispose(lib);

    return mapped;
  }
  // write the updated transitions as a new library
  public static void writeTransitions(Collection<SRMGroup> groups, String filename)
  {
    try
    {
      FileWriter tr = new FileWriter(filename);

      tr.write("ProteinId\tModifiedPeptideSequence\tNormalizedRetentionTime\tPrecursorMz\tPrecursorCharge\tiRT\tProductMz\tLibraryIntensity\tReportedRT\n");
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
          tr.write(grp.getRT()               +"\n"); // ReportedRT
        }

      tr.close();
    }
    catch (IOException e) {}
  }
}
