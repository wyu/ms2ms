package org.ms2ms.data.ms;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;
import org.ms2ms.Disposable;
import org.ms2ms.data.Point;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;
import toools.collections.Lists;

import java.util.*;

public class SRM implements Cloneable, Disposable
{
  private float mFragmentMz, mLibraryIntensity, mApex, mArea, mPkPct=0, mPkPctAll=0, mFillTime=0;
  Range<Double> mPeakBoundary;

  private List<LcMsPoint> mXIC;
  private LcMsFeature mFeature;

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
  public LcMsPoint addXIC(float rt, float ai, float mz, int scan)
  {
    if (mXIC==null) mXIC = new ArrayList<>();
    mXIC.add(new LcMsPoint(rt,ai,mz,scan)); return mXIC.get(mXIC.size()-1);
  }

  public SRM setFeature(Point s) { if (s!=null) mFeature = new LcMsFeature(s); return this; }
  public SRM setFeature(LcMsFeature s) { mFeature=s; return this; }
  public SRM setFillTime(float s) { mFillTime=s; return this; }

  public SRM calPeakPct(double rt, double span, int apex_pts, double peak_base)
  {
    if (!Tools.isSet(getXIC())) return this;

    double inside=0, outside=0, all=0;
    Range<Double> inner = Range.closed(rt-span, rt+span), outer = Range.closed(rt-2d*span, rt+2d*span);

    List<Double> tops = new ArrayList<>();
    int apex_i=-1; double best=Float.MAX_VALUE;
    for (int i=0; i<getXIC().size(); i++)
    {
      Point p = get(i);
      if       (inner.contains(p.getX())) { inside+=p.getY(); tops.add(p.getY()); }
      else if (!outer.contains(p.getX())) outside+=p.getY();
      all += p.getY();

      if (Math.abs(p.getX()-rt)<best) { apex_i=i; best=Math.abs(p.getX()-rt); }
    }
    if (getFeature()!=null && tops.size()>apex_pts) {
      Collections.sort(tops, Ordering.natural().reversed());
      getFeature().setApex(Stats.mean(tops.subList(0, apex_pts)));
      mApex = (float )getFeature().getApex();
    }
    if (apex_i>0 && getFeature()!=null)
    {
      double cut = getFeature().getApex()*peak_base, left=0, right=0, area=0;
      for (int i=apex_i; i>0; i--)
        if (get(i).getIntensity()<cut && (i<=0 || get(i-1).getIntensity()<cut)) { left=get(i).getRT(); break; } else area+=get(i).getIntensity();

      for (int i=apex_i+1; i<getXIC().size()-1; i++)
        if (get(i).getIntensity()<cut && (i>=getXIC().size() || get(i+1).getIntensity()<cut)) { right=get(i).getRT(); break; } else area+=get(i).getIntensity();

      if (left>0 && right>left) { mPeakBoundary = Range.closed(left, right); mArea=(float )area; getFeature().setArea(area); }
    }
    mPkPct    = (float )(100f*inside/(inside+outside));
    mPkPctAll = (float )(100f*inside/all);

    tops = (List )Tools.dispose(tops);

    return this;
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
        if (!xs.get(getFragmentMz()).contains(x)) addXIC(x, baseline);

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
  @Override
  public String toString()
  {
    String out = "";
    if (Tools.isSet(getXIC())) out = "xic="+getXIC().size();
    if (mFeature!=null)        out += ", " + mFeature.toString();

    return out;
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
}

