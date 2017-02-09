package org.ms2ms.data.ms;

import org.ms2ms.math.Histogram;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.util.EnumMap;
import java.util.List;

/** Representation of an individual component score in a composite scheme
 *
 * Created by yuw on 10/10/16.
 */
public class ScoreModel
{
  public enum eType
  {
    exact("Exact"), open("Open"), all("ALL"), sim("Simulated");

    private final String name;
    eType(String n) { this.name = n; }

    public String getName() { return name; }
  };

  private String mName;
  private EnumMap<eType, Histogram> mDecoys;
  private EnumMap<eType, Double> mOffsets;
  //private Double mCenter, mSigma, mWeight=1d, mFactor=1d, mQvalSlope, mQvalIntercept, mBaseline;
  private Double mBaseline,mCenter;

  public ScoreModel() { super(); }
  public ScoreModel(String s) { super(); mName=s; }

//  public Double getQvalSlope()     { return mQvalSlope; }
//  public Double getQvalIntercept() { return mQvalIntercept; }
  public Double getBaseline()      { return mBaseline; }
  public Double getCenter()        { return mCenter; }

//  public double getWeight() { return mWeight; }
  public Double getOffset(eType t) { return mOffsets.get(t); }
  public String getName() { return mName; }

//  public Double calcScoreAtQval(double q) { return getQvalIntercept()!=null && getQvalSlope()!=null?q*getQvalSlope()+getQvalIntercept():null; }

  public ScoreModel setCenter(double s) { mCenter=s; return this; }
//  public ScoreModel setSigma(double s) { mSigma=s; return this; }
//  public ScoreModel setQvalCoeffs(double slope, double intercept)
//  {
//    mQvalSlope=slope; mQvalIntercept=intercept;
//    return this;
//  }
//  public ScoreModel setFactor(double s) { mFactor=s; return this; }
//  public ScoreModel setWeigth(double s) { mWeight=s; return this; }
  public ScoreModel setOffset(eType t, double s)
  {
    if (mOffsets==null) mOffsets = new EnumMap<>(eType.class);
    mOffsets.put(t, s);
    return this;
  }

  public ScoreModel addExactDecoy(double s) { return add(s, eType.exact); }
  public ScoreModel addOpenDecoy( double s) { return add(s, eType.open); }

  public ScoreModel add(double s, eType t)
  {
    if (mDecoys==null) mDecoys = new EnumMap<>(eType.class);
    if (mDecoys.get(t)==null) mDecoys.put(t, new Histogram(t.getName()));

    mDecoys.get(t).add(s);
    return this;
  }

  public void generate(eType... types)
  {
    if (mDecoys!=null && Tools.isSet(types))
      for (eType type : types)
      {
        Histogram H = mDecoys.get(type);
        if (H!=null) H=H.generate2pts(15, 0.5);
        if (H!=null) H=H.assessTruncated(0);
      }
  }
  public void generate(int steps, eType... types)
  {
    if (mDecoys!=null && Tools.isSet(types))
      for (eType type : types)
      {
        Histogram H = mDecoys.get(type);
        if (H!=null) H=H.generate(steps);
      }
  }
  public double score(double s, eType type)
  {
    return s-getOffset(type);
  }
  public ScoreModel model(eType main, eType open)
  {
    if (mDecoys!=null &&
        mDecoys.get(main)!=null && mDecoys.get(main).getData()!=null && mDecoys.get(main).getData().size()>1 &&
        mDecoys.get(open)!=null && mDecoys.get(open).getData()!=null && mDecoys.get(open).getData().size()>1)
    {
      setOffset(main, 0d).setOffset(open, 0d);
      setOffset(eType.open, 10d * Math.log10((double )mDecoys.get(open).getData().size()/(double )mDecoys.get(main).getData().size()));
    }
    else
    {
      setOffset(eType.open, 10d * Math.log10(getOffset(open)/Math.max(0.5, getOffset(main))));
      setOffset(main, 0d);
    }

//    // use the centroid and upper quartile for normalization. Gaussian fit is not robust enough
//    generate(eType.exact,eType.open);
//
//    if (Tools.isSet(mDecoys))
//    {
//      if (mDecoys.get(main)!=null && mDecoys.get(main).getData()!=null && mDecoys.get(main).getData().size()>2)
//      {
////        mOffsets = new EnumMap<>(eType.class);
//        double count0 = Math.log10(mDecoys.get(main).getData().size())*10;
//        // setup the offsets using the counts
//        for (eType t : mDecoys.keySet())
//          setOffset(t, 10d * Math.log10(mDecoys.get(t).getData().size()) - count0);
//      }
//      else if (getOffset(main)>0)
//      {
//        double count0 = Math.log10(getOffset(main))*10;
//        // setup the offsets using the counts
//        for (eType t : mDecoys.keySet())
//          setOffset(t, 10d * Math.log10(getOffset(t)) - count0);
//      }
//      else
//      {
//        // disable the offsets
//        for (eType t : mDecoys.keySet()) setOffset(t, 0d);
//      }

//      Histogram all = new Histogram(eType.all.getName());
//      all.addAll(mDecoys.get(main).getData());
//      for (Double x : mDecoys.get(open).getData()) all.add(x-getOffset(open));

////      mOffsets = new EnumMap<>(eType.class);
//      for (eType t : mDecoys.keySet())
//      {
//        setOffset(t, 0d);
//        try
//        {
//          setOffset(t, mDecoys.get(t).getCenter() - mDecoys.get(main).getCenter());
//          for (Double x : mDecoys.get(t).getData()) all.add((x-getOffset(t)));
//        }
//        // skip the part if we run into some NULL pointer
//        catch (NullPointerException e)
//        {
//          // deposit the uncalibrated points
//          for (Double x : mDecoys.get(t).getData()) all.add(x);
//        }
//      }
      // avoid the first 2 points to reduce the effect of truncated distribution
//      Tools.put(mDecoys, eType.all, all.generate(all.getData().size()>100?25:15).assessTruncated(2));

//      List<Double> baseline = all.getData().subList(Math.max(all.getData().size() - 5, 0),all.getData().size());
//      if (Tools.isSet(baseline))
//      {
//        if (baseline.size()<2) mBaseline=baseline.get(0);
//        else mBaseline = Stats.mean(baseline)+Stats.stdev(baseline)*2d;
//      }
//      mCenter=all.getCenter();
//      mCenter=all.getCenter(); mSigma=all.getSigma()!=null?Math.abs(all.getSigma()):null;
//    }

    return this;
  }

}
