package org.ms2ms.data.ms;

import com.google.common.collect.Ordering;
import org.ms2ms.math.Histogram;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.util.*;

/** Representation of an individual component score in a composite scheme
 *
 * Created by yuw on 10/10/16.
 */
public class ScoreModel
{
  public enum eType
  {
    exact("Exact"), open("Open"), all("ALL"), sim("Simulated"), y("y"), b("b"), bs("bootstrapped");

    private final String name;
    eType(String n) { this.name = n; }

    public String getName() { return name; }
  };

  private String mName;
  private EnumMap<eType, Histogram> mDecoys, mNorms;
  private EnumMap<eType, Double> mOffsets;
  private Double mBaseline,mCenter;

  public ScoreModel() { super(); }
  public ScoreModel(String s) { super(); mName=s; }

  public Double getBaseline()      { return mBaseline; }
  public Double getCenter()        { return mCenter; }

  public Double getOffset(eType t) { return mOffsets.get(t); }
  public String getName() { return mName; }

  public ScoreModel setCenter(double s) { mCenter=s; return this; }
  public ScoreModel setOffset(eType t, double s)
  {
    if (mOffsets==null) mOffsets = new EnumMap<>(eType.class);
    mOffsets.put(t, s);
    return this;
  }
  public void clear(eType... types)
  {
    if (Tools.isSet(types))
      for (eType t : types)
      {
        if (mDecoys!=null && mDecoys.get(t)!=null) mDecoys.get(t).clear();
        if (mNorms !=null && mNorms.get( t)!=null) mNorms.get( t).clear();
      }
  }
  public ScoreModel add(Histogram hist, boolean decoy, eType t)
  {
    if      ( decoy && mDecoys==null) mDecoys = new EnumMap<>(eType.class);
    else if (!decoy && mNorms ==null) mNorms  = new EnumMap<>(eType.class);

    if (decoy) mDecoys.put(t, hist); else mNorms.put(t, hist);
    return this;
  }
  public ScoreModel add(boolean decoy, eType t, double s, boolean distinct)
  {
    if      ( decoy && mDecoys==null) mDecoys = new EnumMap<>(eType.class);
    else if (!decoy && mNorms ==null) mNorms  = new EnumMap<>(eType.class);

    if      ( decoy && mDecoys.get(t)==null) mDecoys.put(t, new Histogram("decoy_"+t.toString()));
    else if (!decoy &&  mNorms.get(t)==null)  mNorms.put(t, new Histogram( "norm_"+t.toString()));

    if (decoy) mDecoys.get(t).add(s, distinct); else mNorms.get(t).add(s, distinct);

    return this;
  }

  public ScoreModel bootstrapping(int samples)
  {
    add(bootstrapping(samples, true,  eType.y, eType.b), true,  eType.bs);
    add(bootstrapping(samples, false, eType.y, eType.b), false, eType.bs);

    return this;
  }
  public Histogram bootstrapping(int samples, Boolean decoy, eType C, eType N)
  {
    List<Double> Nc=new ArrayList<>(), Cc=new ArrayList<>();

    if ((decoy==null || decoy) && mDecoys!=null)
    {
      if (mDecoys.get(N)!=null) Nc.addAll(mDecoys.get(N).getData());
      if (mDecoys.get(C)!=null) Cc.addAll(mDecoys.get(C).getData());
    }
    if ((decoy==null || !decoy) && mNorms!=null)
    {
      if (mNorms.get(N)!=null) Nc.addAll(mNorms.get(N).getData());
      if (mNorms.get(C)!=null) Cc.addAll(mNorms.get(C).getData());
    }
    // sort the scores
    Collections.sort(Nc); Collections.sort(Cc);

    // boot strap for the combined distribution
    Histogram combo = new Histogram("Bootstrapping estimation of the score distribution: "+decoy+"/"+N+"/"+C);
    Random      RND = new Random(System.nanoTime());
    for (int i=0; i<samples; i++)
    {
      combo.add(Nc.get(RND.nextInt(Nc.size()))+Cc.get(RND.nextInt(Cc.size())));
    }
    combo.generate(25);
    Collections.sort(combo.getData(), Ordering.natural().reverse());
    Tools.dispose(Nc, Cc);

    return combo;
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

    return this;
  }

  public StringBuffer toWiki(StringBuffer buf)
  {
    if (buf==null) buf = new StringBuffer();

    buf.append(getName()+"\n");
    for (eType t : mDecoys.keySet())
    {
      buf.append(t.toString()+"\n");
      buf = mDecoys.get(t).wikiHistogram(buf);
    }

    return buf;
  }
}
