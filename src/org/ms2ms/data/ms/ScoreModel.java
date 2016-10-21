package org.ms2ms.data.ms;

import org.ms2ms.math.Histogram;
import org.ms2ms.utils.Tools;

import java.util.EnumMap;

/** Representation of an individual component score in a composite scheme
 *
 * Created by yuw on 10/10/16.
 */
public class ScoreModel
{
  public enum eType
  {
    exact("Exact"), open("Open"), all("ALL");

    private final String name;
    eType(String n) { this.name = n; }

    public String getName() { return name; }
  };

  private String mName;
  private EnumMap<eType, Histogram> mDecoys;
  private EnumMap<eType, Double> mOffsets;
  private Double mCenter, mSigma, mWeight=1d, mFactor=1d;

  public ScoreModel() { super(); }
  public ScoreModel(String s) { super(); mName=s; }

  public double getWeight() { return mWeight; }
  public String getName() { return mName; }

  public ScoreModel setFactor(double s) { mFactor=s; return this; }
  public ScoreModel setWeigth(double s) { mWeight=s; return this; }

  public ScoreModel addExactDecoy(double s) { return add(s, eType.exact); }
  public ScoreModel addOpenDecoy( double s) { return add(s, eType.open); }

  private ScoreModel add(double s, eType t)
  {
    if (mDecoys==null) mDecoys = new EnumMap<>(eType.class);
    if (mDecoys.get(t)==null) mDecoys.put(t, new Histogram(t.getName()));

    mDecoys.get(t).add(s);
    return this;
  }

  private void generate(eType... types)
  {
    if (mDecoys!=null && Tools.isSet(types))
      for (eType type : types)
      {
        Histogram H = mDecoys.get(type);
        if (H!=null) H=H.generate2pts(15, 0.5);
        if (H!=null) H=H.assessTruncated();
      }
  }
  public double score(double s, eType type)
  {
    return mWeight*(s-mCenter-mOffsets.get(type))/mSigma;
  }
  public ScoreModel model(eType main)
  {
    // use the centroid and upper quartile for normalization. Gaussian fit is not robust enough
    generate(eType.exact,eType.open);

//    mDecoys.get(eType.exact).printHistogram();
//    mDecoys.get(eType.open ).printHistogram();

    if (Tools.isSet(mDecoys))
    {
      Histogram all = new Histogram(eType.all.getName()); mOffsets = new EnumMap<>(eType.class);
      for (eType t : mDecoys.keySet())
      {
        mOffsets.put(t, 0d);
        try
        {
          mOffsets.put(t,mDecoys.get(t).getCenter()-mDecoys.get(main).getCenter())/*, y0=mDecoys.get(t).getSigma()/mDecoys.get(main).getSigma()*/;
          //for (Double x : mDecoys.get(t).getData()) all.add((x-x0)/y0);
          for (Double x : mDecoys.get(t).getData()) all.add((x-mOffsets.get(t)));
        }
        // skip the part if we run into some NULL pointer
        catch (NullPointerException e)
        {
          // deposit the uncalibrated points
          for (Double x : mDecoys.get(t).getData()) all.add(x);
        }
      }
      Tools.put(mDecoys,eType.all, all.generate2pts(25, 0.5).assessTruncated());

      mCenter=all.getCenter(); mSigma=all.getSigma();
    }

    return this;
  }

}
