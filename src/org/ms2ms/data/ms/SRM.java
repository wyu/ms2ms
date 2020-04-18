package org.ms2ms.data.ms;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.Range;
import org.ms2ms.data.Point;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;
import toools.collections.Lists;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class SRM implements Cloneable
{
  private float mFragmentMz, mLibraryIntensity, mApex, mArea, mPkPct=0, mPkPctAll=0, mFillTime=0;
  private List<LcMsPoint> mXIC;
  private LcMsPoint mFeature;

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

  public List<LcMsPoint> getXIC() { return mXIC; }
  public LcMsPoint getFeature() { return mFeature; }
  public LcMsPoint get(int i) { return mXIC.get(i); }

  public LcMsPoint addXIC(float rt, float ai)
  {
    if (ai>0) mXIC.add(new LcMsPoint(rt,ai));
    return mXIC.get(mXIC.size()-1);
  }
  public LcMsPoint addXIC(float rt, float ai, float mz, int scan)
  {
    mXIC.add(new LcMsPoint(rt,ai,mz,scan)); return mXIC.get(mXIC.size()-1);
  }

  public SRM setFeature(Point s) { if (s!=null) mFeature = new LcMsPoint(s); return this; }
  public SRM setFeature(LcMsPoint s) { mFeature=s; return this; }
  public SRM setFillTime(float s) { mFillTime=s; return this; }

  public SRM calPeakPct(double rt, double span, int apex_pts)
  {
    if (!Tools.isSet(getXIC())) return this;

    double inside=0, outside=0, all=0;
    Range<Double> inner = Range.closed(rt-span, rt+span), outer = Range.closed(rt-2d*span, rt+2d*span);

    List<Double> tops = new ArrayList<>();
    for (Point p : getXIC())
    {
      if       (inner.contains(p.getX())) { inside+=p.getY(); tops.add(p.getY()); }
      else if (!outer.contains(p.getX())) outside+=p.getY();
      all += p.getY();
    }
    if (getFeature()!=null && tops.size()>apex_pts) {
      Collections.sort(tops, Ordering.natural().reversed());
      getFeature().setApex(Stats.mean(tops.subList(0, apex_pts)));
//      if (Double.isNaN(getFeature().getApex()))
//        System.out.println();
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
}

