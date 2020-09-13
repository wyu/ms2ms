package org.ms2ms.data.ms;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.jgrapht.Graph;
import org.jgrapht.alg.connectivity.KosarajuStrongConnectivityInspector;
import org.jgrapht.alg.interfaces.StrongConnectivityAlgorithm;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;
import org.ms2ms.Disposable;
import org.ms2ms.data.Point;
import org.ms2ms.math.Fitted;
import org.ms2ms.math.Points;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import toools.collections.Lists;

import java.io.IOException;
import java.io.Writer;
import java.nio.DoubleBuffer;
import java.util.*;

public class SRM implements Cloneable, Disposable, Comparable<SRM>
{
  private SRMGroup.IsoLable mIsotopeLabel = SRMGroup.IsoLable.L;

  private boolean mIsStronglyConnected=false;

  private int mIsotope=0, mSizeNonzero=0, mNumEdges=0, mCharge;
  private float mFragmentMz, mLibraryIntensity, mPkPct=0, mPkPctAll=0, mBackground=0f, mPeaksArea=0f, mPrecursorMz=0f;
  private String mFragmentType;
  Range<Double> mPeakBoundary=null;

  private List<LcMsPoint> mXIC;
  private LcMsFeature mFeature=null, mMostSimilar=null;
  private List<LcMsFeature> mPeaks=null;

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

  public SRMGroup.IsoLable getIsotopeLabel() { return mIsotopeLabel; }

  public boolean isStronglyConnected() { return mIsStronglyConnected; }
  public boolean isIsotopeLabel(SRMGroup.IsoLable s) { return Tools.equals(s, mIsotopeLabel); }
  public int getIsotope() { return mIsotope; }
  public int getXicSize() { return mSizeNonzero; }
  public int getNumEdges() { return mNumEdges; }
  public int getCharge() { return mCharge; }

  public String getFragmentType() { return mFragmentType; }

  public float getFragmentMz() { return mFragmentMz; }
  public float getLibraryIntensity() { return mLibraryIntensity; }
  public float getPrecursorMz() { return mPrecursorMz; }

//  public float getApex()       { return mApex; }
  public float getPeaksArea()       { return mPeaksArea; }
  public float getBackground()        { return mBackground; }
//  public float getSNR(float apex)        { return mBackground!=0?apex/mBackground:0f; }

  public float getPeakPct()    { return mPkPct; }
  public float getPeakPctAll() { return mPkPctAll; }
//  public float getFillTime()   { return mFillTime; }
//  public float getSNR() { return getBackground()>0 && getFeature()!=null ? (float )getFeature().getApex()/getBackground():0f; }
  public Range<Double> getPeakBoundary() { return mPeakBoundary; }
  public List<LcMsFeature> getPeaks() { return mPeaks; }

  public List<LcMsPoint> getXIC() { return mXIC; }
  public LcMsFeature getFeature() { return mFeature; }
  public LcMsFeature getMostSimilar() { return mMostSimilar; }

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
  public SRM updateFeature(Point s)
  {
    if      (mFeature==null && s!=null) mFeature = new LcMsFeature(s);
    else if (mFeature!=null && s!=null) mFeature.set(s);
    return this;
  }

//  public SRM setFeature(Point s) { if (s!=null) mFeature = new LcMsFeature(s); return this; }
  public SRM setIsotopeLabel(SRMGroup.IsoLable s) { mIsotopeLabel=s; return this; }
  public SRM setFeature(LcMsFeature s) { mFeature=s; return this; }
  public SRM setPeakBoundary(Range<Double> s) { mPeakBoundary=s!=null?Range.closed(s.lowerEndpoint(), s.upperEndpoint()):null; return this; }

  public SRM setFragmentType(String s) { if (Strs.isSet(s)) { mFragmentType=s; } return this; }

  public SRM isStronglyConnected(boolean s) { mIsStronglyConnected = s; return this; }

  //  public SRM setFillTime(float s) { mFillTime=s; return this; }
  public SRM setBackground(float s) { mBackground=s; return this; }
  public SRM setNumEdges(int s) { mNumEdges=s; return this; }
  public SRM setCharge(int s) { mCharge=s; return this; }
  public SRM setIsotope(int s) { mIsotope=s; return this; }
  public SRM setPrecursorMz(float s) { mPrecursorMz =s; return this; }
  public SRM calPeakExclusivity(double rt, double quan_span, Range<Double> boundary)
  {
    if (!Tools.isSet(getXIC()) || !Tools.isSet(boundary)) return this;

    if (boundary.lowerEndpoint()-quan_span > rt)
      System.out.print("");
    // change the definition on May 16, 2020. outside is now 3x of the LC peak span, instead of 2x
    double inside_front=0, outside_front=0, all=0;
    Range<Double> inner_front = Range.closed(boundary.lowerEndpoint(), rt),
                  outer_front = Range.closed(boundary.lowerEndpoint()-quan_span, rt);

    double best=Float.MAX_VALUE;
    for (int i=0; i<getXIC().size(); i++)
    {
      Point p = get(i);
      if (inner_front.contains(p.getX())) inside_front+=p.getY();
      // change the definition on May 16, 2020. outside is sum of intensity within 'outside', instead of intensities outside of 'outside'
      if (outer_front.contains(p.getX()))  outside_front+=p.getY();

      all += p.getY();

      if (Math.abs(p.getX()-rt)<best) { best=Math.abs(p.getX()-rt); }
    }
//    if (apex_i>0 && getFeature()!=null && getFeature().getApex()>0)
//    {
//      double cut = getFeature().getApex()*peak_base, left=0, right=0, area=0;
//      for (int i=apex_i; i>0; i--)
//        if (get(i).getIntensity()<cut && (i<=0 || get(i-1).getIntensity()<cut)) { left=get(i).getRT(); break; } else area+=get(i).getIntensity();
//
//      for (int i=apex_i+1; i<getXIC().size()-1; i++)
//        if (get(i).getIntensity()<cut && (i>=getXIC().size() || get(i+1).getIntensity()<cut)) { right=get(i).getRT(); break; } else area+=get(i).getIntensity();
//
//      if (left>0 && right>left) { mPeakBoundary = Range.closed(left, right); getFeature().setArea(area); }
//    }
    // change the definition on May 16, 2020. exclusivity is now inside/outside, for local exclusivity, instead of inside_intensity / (intensity outside of 2x span)
    mPkPct     = (float )(100f*inside_front/outside_front); // a local exclusivity (+-3x LC peak span)
    mPkPctAll  = (float )(100f*getFeature().getArea()/all);      // this is now a global (+-5min) exclusivity
//    tops = (List )Tools.dispose(tops);

    return this;
  }
  public LcMsFeature setPeakBoundary(LcMsFeature lead, double peak_base, float lc_width)
  {
    if (lead!=null && lead.getY()>0 && lead.getApexPos()>=0)
    {
      double cut = getFeature().getY()*peak_base/100d, left=0, right=0, area=0, r0 = lead.getX();
      for (int i=lead.getApexPos(); i>0; i--)
      {
        // no more than a full expected LC peak width, effectively 2x the expected
        if (r0-get(i).getX()>lc_width || (get(i).getIntensity()<cut && (i<=0 || get(i-1).getIntensity()<cut)))
        {
          left=get(i).getRT(); break;
        }
        else if (i==1)
        {
          area+=get(i).getIntensity()+get(i-1).getIntensity(); left=get(0).getRT(); break;
        }
        else area+=get(i).getIntensity();
      }

      for (int i=lead.getApexPos()+1; i<getXIC().size()-1; i++)
      {
        // no more than a full expected LC peak width, effectively 2x the expected
        if (get(i).getX()-r0>lc_width || (get(i).getIntensity()<cut && (i>=getXIC().size() || get(i+1).getIntensity()<cut)))
        {
          right=get(i).getRT(); break;
        }
        else if (i==getXIC().size()-2)
        {
          area+=get(i).getIntensity()+get(i+1).getIntensity(); right=get(i+1).getRT(); break;
        }
        else area+=get(i).getIntensity();
      }

      if (left>0 && right>left) { mPeakBoundary = Range.closed(left, right); lead.setArea(area); }
    }
    return lead;
  }
  // center = the targeted RT,
  // span   = the span of the RT window where the peaks are expected
  public List<LcMsFeature> detectPeak(float center, float span)
  {
    if (!Tools.isSet(getXIC()) || getXIC().size()<=5) return null;

    List<Point>      deri = Points.deriv1stBySG5(getXIC());
    mPeaks = new ArrayList<>();

    // get the start and end first
    mPeaks.add(new LcMsFeature(getXIC().subList(0, 5)));
    mPeaks.add(new LcMsFeature(getXIC().subList(getXIC().size()-5, getXIC().size())));

    if (Tools.isSet(deri) && deri.size()>2)
    {
      // figure out the background level
      List<Double> ys = Points.toYs(deri);

      Collections.sort(ys, Ordering.natural().reversed());
      double cutoff = ys.get((int )(ys.size()*0.25));

      double deri_max=0, deri_min= Double.MAX_VALUE; int pos_deri_max=-1; LcMsFeature lastF=null;
      if (Tools.isSet(deri))
        for (int i=0; i<deri.size()-1; i++)
        {
          if      (deri.get(i).getY()>deri_max) { deri_max=deri.get(i).getY(); pos_deri_max=i; }
          else if (deri.get(i).getY()<deri_min)
          {
            deri_min=deri.get(i).getY();
            if (lastF!=null && i+2<mXIC.size()) lastF.setUpper(get(i+2).getX());
          }

          if (deri.get(i).getY()>=0 && deri.get(i+1).getY()<=0)
          {
            // the peak top is at the zero-intercept
            Point top = Points.interpolateByY(deri.get(i), deri.get(i+1), 0d);
            Point apx = Points.interpolate(getXIC().get(i+2), getXIC().get(i+3), top.getX());

            if (deri_max>cutoff)
            {
              double area = LcMsPoint.sumY(getXIC().subList(pos_deri_max+2, i+3));
              if (area>0)
              {
                lastF = new LcMsFeature(top.getX(), apx.getY(), deri_max).setApexPos(i+2
                ).setPointWidth(i+1-pos_deri_max).setArea(area).setLower(get(pos_deri_max+2).getX());
                mPeaks.add(lastF);
              }
              deri_max=0;
            }
          }
        }
      ys = (List )Tools.dispose(ys);
      // calculate the SNR
      if (mPeaks.size()>1)
      {
        int mid = Math.round(mPeaks.size()/2);
        Collections.sort(mPeaks, new LcMsFeature.SimilarityDesendComparator());
        mMostSimilar = mPeaks.get(0);
        Collections.sort(mPeaks, new LcMsFeature.AreaDesendComparator());
        // figor out the background level
        double background = LcMsFeature.sumArea(mPeaks.subList(mid, mPeaks.size()))/(mPeaks.size()-mid);
        // set the total peak area
        mPeaksArea = 0f;
        for (LcMsFeature p : mPeaks)
        {
          p.setSNR(p.getArea()/background);
          mPeaksArea+=p.getArea();
        }
        for (LcMsFeature p : mPeaks) p.setExclusivity(100f*p.getArea()/mPeaksArea);
      }
    }
//    else
//      System.out.println();

    deri = (List )Tools.dispose(deri);

    return mPeaks;
 }
  public SRM clone()
  {
    SRM cloned = new SRM(getFragmentMz(), getLibraryIntensity());

    cloned.mIsotopeLabel = mIsotopeLabel;
    cloned.mIsotope = mIsotope; cloned.mSizeNonzero = mSizeNonzero;
    cloned.mBackground = mBackground; cloned.mPeaksArea = mPeaksArea; cloned.mPrecursorMz = mPrecursorMz;
    cloned.mPkPct = mPkPct; cloned.mPkPctAll = mPkPctAll;

    if (mPeakBoundary!=null) cloned.mPeakBoundary = Range.closed(mPeakBoundary.lowerEndpoint(), mPeakBoundary.upperEndpoint());
    if (mFeature!=null) cloned.mFeature = new LcMsFeature(mFeature);
    if (mMostSimilar!=null) cloned.mMostSimilar = new LcMsFeature(mMostSimilar);

    if (Tools.isSet(mPeaks)) cloned.mPeaks = new ArrayList<>(mPeaks);

    if (Tools.isSet(mXIC))
    {
      cloned.mXIC     = new ArrayList<>();
      for (LcMsPoint p : mXIC)
        if (p!=null) cloned.mXIC.add(p);
    }

    return cloned;
  }
  public SRM mutate(ListMultimap<Integer, Float> frag_bank, Random rnd)
  {
//    mApex=mArea=0f;
    mPkPct=mPkPctAll = 0f;
    mXIC.clear();

    Integer idx = (int )Math.round(getFragmentMz()*0.01);
    mFragmentMz = Lists.getRandomSubset(frag_bank.get(idx),1, rnd).get(0);

    return this;
  }
  public SRM fill(float baseline, Multimap<Float, Float> xs)
  {
    if (Tools.isSet(xs))
    {
      if (mXIC==null) mXIC = new ArrayList<>();
      mSizeNonzero = mXIC.size();
      for (Float x : xs.get(0f))
        // need to bypass the non-zero check in addXIC
        if (!xs.get(getFragmentMz()).contains(x)) mXIC.add(new LcMsPoint(x,baseline));
    }

    if (Tools.isSet(getXIC())) Collections.sort(getXIC());

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
  // return the estimate ratio to the control chennel by template fitting
  public double[] ratio2control(SRM ctrl, float edge)
  {
    if (ctrl==null || !Tools.isSet(ctrl.getXIC()) || !Tools.isSet(getXIC())) return new double[] {0d, 0d};

//    List<WeightedObservedPoint> pts = new ArrayList<>();
    List<Point> pts = new ArrayList<>(), ptt = new ArrayList<>();

//    System.out.println("\nlogAssay\tlogCtrl\tweight\tRT\tAssay\tCtrl");
//    if (ctrl.getPeakBoundary()==null)
//      System.out.println("");
//
    double r1=0d, r2=0, diff=0d, n=0d,
        rtL=ctrl.getPeakBoundary()!=null?ctrl.getPeakBoundary().lowerEndpoint()-edge:0d,
        rtR=ctrl.getPeakBoundary()!=null?ctrl.getPeakBoundary().upperEndpoint()+edge:Double.MAX_VALUE;
    for (int i=0; i<getXIC().size(); i++)
      if (ctrl.get(i).getX()>=rtL && ctrl.get(i).getX()<=rtR && get(i).getIntensity()>0 && ctrl.get(i).getIntensity()>0)
      {
        r1 +=      get(i).getIntensity();
        r2 += ctrl.get(i).getIntensity();

        diff += Math.log10(get(i).getIntensity())-Math.log10(ctrl.get(i).getIntensity()); n++;
        pts.add(new Point(Math.log10(get(i).getIntensity()), Math.log10(ctrl.get(i).getIntensity())));
        ptt.add(new Point(Math.log10(get(i).getIntensity())-Math.log10(ctrl.get(i).getIntensity()), Math.sqrt(ctrl.get(i).getIntensity())));

//        WeightedObservedPoint pt = new WeightedObservedPoint(Math.sqrt(ctrl.get(i).getIntensity()), Math.log10(get(i).getIntensity()), Math.log10(ctrl.get(i).getIntensity()));
////        pts.add(pt);
//        System.out.println(pt.getX()+"\t"+pt.getY()+"\t"+pt.getWeight()+"\t"+get(i).getX()+"\t"+get(i).getIntensity()+"\t"+ctrl.get(i).getIntensity());
      }

//    System.out.println("\nRT\tAssay\tCtrl");
//    for (int i=0; i<getXIC().size(); i++)
//      if (get(i).getIntensity()>0 || ctrl.get(i).getIntensity()>0)
//      {
//        System.out.println(get(i).getX()+"\t"+get(i).getIntensity()+"\t"+ctrl.get(i).getIntensity());
//      }

    if (n>3)
    {
      diff /=n;
      List<Point> errs = new ArrayList<>();
      for (double k=-2; k<2; k+=0.1)
      {
        double err2=0;
        for (Point pt : pts)
        {
          err2 += Math.abs(pt.getY()+diff+k-pt.getX());
        }
        errs.add(new Point(diff-k, 100d/err2));
      }
      Collections.sort(ptt);
      double dif = Points.centroid(errs), difw = Points.centroid(ptt);
      return (new double[] {
          ctrl.getFeature().getArea()*Math.pow(10, dif),
          ctrl.getFeature().getArea()*Math.pow(10, diff),
          ctrl.getFeature().getArea()*Math.pow(10, difw)});
    }
    else if (n==0) return new double[] {0d, 0d};

    return new double[] {r1/r2, r1/r2};
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
//        Point p = Points.interpolate(getXIC(), x, false);
        LcMsPoint p = LcMsPoint.interpolate(getXIC(), x, true);
        if (p==null) p = new LcMsPoint(x, 0);

        pts.add(p);
      }
      mXIC = (List )Tools.dispose(mXIC);
      mXIC = pts;
    }
    return this;
  }
  @Override
  public String toString()
  {
    String out = getIsotopeLabel().toString();

    out = Strs.extend(out, Strs.isSet(getFragmentType())?getFragmentType():("m/z"+Tools.d2s(getFragmentMz(), 3)), ",");
//    if (getIsotope()>0)        Strs.extend(out, getIsotope()+"", "iso ");
    if (Tools.isSet(getXIC())) out = Strs.extend(out,getXIC().size()+"", ", xic ");
    if (mFeature!=null)        out = Strs.extend(out,  mFeature.toString(), ", ");

    return out;
  }
  public String getUID(String s)
  {
    return s+"_"+Tools.d2s(getFragmentMz(), 4)+"_"+getFragmentType();
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
  public static void headerFeatures(Writer w) throws IOException
  {
    w.write("isoL\tiso\tNumPts\tPkEx\tPkExAll\tsc\tedge\tfeature.rt\tfeature.ai\tinjection\tfeature.apex\tfeature.area\tfeature.ppm\tfeature.snr\tppm.stdev\tppm.stdev.ex\trule\tinbound.snr\tlower\tupper");
  }
  public static void headerAssayFeatures(Writer w, SRMGroup.IsoLable assay) throws IOException
  {
    String H = assay.toString();
    w.write(H+".PkExAll\t"+H+".sc\t"+H+".edge\t"+H+".rt\t"+H+".ai\t"+H+".inj\t"+H+".apex\t"+H+".area\t"+H+".ppm\t"+H+".snr\t"+H+".r\t"+H+".delta\t"+H+".fmol");
  }
  public void printFeature(Writer w) throws IOException
  {
    w.write(getIsotopeLabel()+"\t");
    w.write(getIsotope()+"\t");
    w.write(getXicSize()+"\t");
    w.write(Tools.d2s(getPeakPct(),2)+"\t");
    w.write(Tools.d2s(getPeakPctAll(),2)+"\t");
    w.write(isStronglyConnected()+"\t");
    w.write(getNumEdges()+"\t");

    if (getFeature()!=null)
    {
      w.write(Tools.d2s(getFeature().getX(),3)+"\t");
      w.write(Tools.d2s(getFeature().getY(),2)+"\t");
      w.write(Tools.d2s(getFeature().getFillTime(),2)+"\t");
      w.write(Tools.d2s(getFeature().getApex(),2)+"\t");
      w.write(Tools.d2s(getFeature().getArea(),2)+"\t");
      w.write(Tools.d2s(getFeature().getPPM(),2)+"\t");
      w.write(Tools.d2s(getFeature().getSNR(),2)+"\t");
      w.write(Tools.d2s(getFeature().getMzStdev(),2)+"\t");
      w.write(Tools.d2s(getFeature().getMzStdevEx(),2)+"\t");
      w.write(getFeature().wasBasedOn()+"\t");
    } else
    {
      w.write("\t\t\t\t\t\t\t\t\t\t");
    }
    w.write((Tools.isSet(getPeaks())?Tools.d2s(getPeaks().get(0).getSNR(),2):"")+"\t");
    w.write(Tools.d2s(getPeakBoundary()!=null?getPeakBoundary().lowerEndpoint():0,2)+"\t");
    w.write(Tools.d2s(getPeakBoundary()!=null?getPeakBoundary().upperEndpoint():0,2));
  }
  public void printKeyFeatures(Writer w) throws IOException
  {
    w.write(Tools.d2s(getPeakPctAll(),2)+"\t");
    w.write(isStronglyConnected()+"\t");
    w.write(getNumEdges()+"\t");

    if (getFeature()!=null)
    {
      w.write(Tools.d2s(getFeature().getX(),3)+"\t");
      w.write(Tools.d2s(getFeature().getY(),2)+"\t");
      w.write(Tools.d2s(getFeature().getFillTime(),2)+"\t");
      w.write(Tools.d2s(getFeature().getApex(),2)+"\t");
      w.write(Tools.d2s(getFeature().getArea(),2)+"\t");
      w.write(Tools.d2s(getFeature().getPPM(),2)+"\t");
      w.write(Tools.d2s(getFeature().getSNR(),2)+"\t");
      w.write(getFeature().getRatio2Ctrl()+"\t");
      w.write(getFeature().getDelta2Ctrl()+"\t");
      w.write(getFeature().getFmols()+"");
    }
    else
    {
      w.write("\t\t\t\t\t\t\t\t\t");
    }
  }
  public static void printKeyBlanks(Writer w) throws IOException
  {
    w.write("\t\t\t\t\t\t\t\t\t\t");
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
    Map<String, Object> props = new HashMap<>();

    // count the edges from the nodes
    for (SRM v : net.vertexSet())
      v.setNumEdges(net.edgesOf(v).size());

    StrongConnectivityAlgorithm<SRM, DefaultWeightedEdge> scAlg =
        new KosarajuStrongConnectivityInspector(net);

    // a graph is said to be strongly connected if every vertex is reachable from every other vertex
    List<Graph<SRM, DefaultWeightedEdge>> stronglyConnectedSubgraphs = scAlg.getStronglyConnectedComponents();

    List<Integer> sc = new ArrayList<>(stronglyConnectedSubgraphs.size());

    // prints the strongly connected components
    int n=0;
    for (int i = 0; i < stronglyConnectedSubgraphs.size(); i++)
      if (stronglyConnectedSubgraphs.get(i).vertexSet().size()>1)
      {
        sc.add(stronglyConnectedSubgraphs.get(i).vertexSet().size());
        if (stronglyConnectedSubgraphs.get(i).vertexSet().size()>n)
        {
          props.put("SCG", stronglyConnectedSubgraphs.get(i));
          n = stronglyConnectedSubgraphs.get(i).vertexSet().size();
        }
      }

    if (props.get("SCG")!=null)
      for (SRM srm : ((Graph<SRM, DefaultWeightedEdge> )props.get("SCG")).vertexSet())
        srm.isStronglyConnected(true);

    Collections.sort(sc, Ordering.natural().reversed());
//    if (sc.size()>2)
//      System.out.print("");

    props.put("SCC: length", stronglyConnectedSubgraphs.size());
    props.put("SCC: node size", sc.toArray(new Integer[sc.size()]));
    props.put("Network: node size", net.vertexSet().size());
    props.put("Network: edge size", net.edgeSet().size());

//    BronKerboschCliqueFinder<SRM, DefaultWeightedEdge> clique = new BronKerboschCliqueFinder(net);
    return props;
  }
}

