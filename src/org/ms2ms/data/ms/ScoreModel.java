package org.ms2ms.data.ms;

import org.ms2ms.data.Point;
import org.ms2ms.math.Histogram;
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
  private Double mCenter, mSigma, mWeight, mScoreOffset=0d, mFactor=1d;

  public ScoreModel() { super(); }
  public ScoreModel(String s) { super(); mName=s; }

  public String getName() { return mName; }

  public ScoreModel setFactor(double s) { mFactor=s; return this; }

  public ScoreModel addExactDecoy(double s) { return add(s, eType.exact); }
  public ScoreModel addOpenDecoy( double s) { return add(s, eType.open); }

  private ScoreModel add(double s, eType t)
  {
    if (mDecoys==null) mDecoys = new EnumMap<>(eType.class);
    if (mDecoys.get(t)==null) mDecoys.put(t, new Histogram(t.getName()));

    mDecoys.get(t).add(s);
    return this;
  }

  public ScoreModel model(eType main)
  {
    // use the centroid and upper quartile for normalization. Gaussian fit is not robust enough
    mDecoys.get(eType.exact).generate2pts(15, 0.5).assessTruncated();
    mDecoys.get(eType.open ).generate2pts(25, 0.5).assessTruncated();

    mDecoys.get(eType.exact).printHistogram();
    mDecoys.get(eType.open ).printHistogram();

    Histogram all = new Histogram(eType.all.getName());
    for (eType t : mDecoys.keySet())
    {
      double x0=mDecoys.get(t).getCenter()-mDecoys.get(main).getCenter(), y0=mDecoys.get(t).getSigma()/mDecoys.get(main).getSigma();
      for (Double x : mDecoys.get(t).getData()) all.add((x-x0)/y0);
    }
    mDecoys.put(eType.all, all.generate2pts(25, 0.5).assessTruncated());

    return this;
  }

}
