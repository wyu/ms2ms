package org.ms2ms.data.ms;

import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.ms2ms.data.Binary;
import org.ms2ms.math.Histogram;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.*;

/** Representation of an individual component score in a composite scheme
 *
 * Created by yuw on 10/10/16.
 */
public class ScoreModel implements Binary
{
  public enum eType
  {
    exact("Exact"), open("Open"), all("ALL"), sim("Simulated"), y("y"), b("b"), bs("bootstrapped");

    private final String name;
    eType(String n) { this.name = n; }

    public String getName() { return name; }
  };

  private String                    mName;
  private EnumMap<eType, Histogram> mDecoys, mNorms, mBoth;
  private EnumMap<eType, Double>    mOffsets;
  private Double                    mBaseline,mCenter;
  private Integer                   mBootstrapSize=1000;
  private SortedSetMultimap<String, Double> mBootsrapped;

  public static float[] sTrypicity = {0f,0f,0f,0f,0.5f,1f,2f};

  public ScoreModel() { super(); }
  public ScoreModel(String s) { super(); mName=s; }

  public Double getBaseline()      { return mBaseline; }
  public Double getCenter()        { return mCenter; }

  public Double getOffset(eType t) { return mOffsets.get(t); }
  public String getName() { return mName; }
  public Integer getBootstrapSize() { return mBootstrapSize; }

  public int size(boolean decoy, eType t)
  {
    if ( decoy && mDecoys!=null && mDecoys.get(t)!=null && mDecoys.get(t).getData()!=null) return mDecoys.get(t).getData().size();
    if (!decoy && mNorms !=null &&  mNorms.get(t)!=null &&  mNorms.get(t).getData()!=null) return  mNorms.get(t).getData().size();
    return 0;
  }

  public Histogram get(Boolean decoy, eType t)
  {
    if (decoy==null &&          mBoth  !=null) return   mBoth.get(t);
    if (decoy!=null && decoy && mDecoys!=null) return mDecoys.get(t);
    if (decoy!=null &&!decoy && mNorms !=null) return  mNorms.get(t);
    return null;
  }
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
  public ScoreModel add(Histogram hist, Boolean decoy, eType t)
  {
    if (decoy==null)
    {
      if (mBoth==null) mBoth = new EnumMap<>(eType.class);
      mBoth.put(t, hist);
    }
    else
    {
      if      ( decoy && mDecoys==null) mDecoys = new EnumMap<>(eType.class);
      else if (!decoy && mNorms ==null) mNorms  = new EnumMap<>(eType.class);

      if (decoy) mDecoys.put(t, hist); else mNorms.put(t, hist);
    }
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

  public ScoreModel fitEval(int samples)
  {
    // prepare for the e-val calculation
    if (mDecoys!=null && mDecoys.get(eType.bs)!=null) mDecoys.get(eType.bs).fitEval(samples,2,0.95);
    if (mNorms !=null &&  mNorms.get(eType.bs)!=null)  mNorms.get(eType.bs).fitEval(samples,2,0.95);
    if (mBoth  !=null &&   mBoth.get(eType.bs)!=null)   mBoth.get(eType.bs).fitEval(samples*2,3, 0.95);

    return this;
  }
  // if y or b is absence, we'll allow it to extend up to 'extension' as in the algorithm itself
  public ScoreModel bootstrapping(int samples, double jitter, double extension, double trypicity)
  {
    mBootstrapSize=samples;
    add(bootstrapping(samples, jitter, extension, trypicity, true,  eType.y, eType.b), true,  eType.bs);
    add(bootstrapping(samples, jitter, extension, trypicity, false, eType.y, eType.b), false, eType.bs);
    // now do the combined dist
    add(bootstrapping(samples, jitter, extension, trypicity, null,  eType.y, eType.b), null,  eType.bs);

    return this;
  }
  public Histogram bootstrapping(int samples, double jitter, double extension, double trypicity, Boolean decoy, eType C, eType N)
  {
    List<Double> Nc=new ArrayList<>(), Cc=new ArrayList<>();

    if (decoy==null&& mDecoys!=null && mNorms!=null)
    {
      // want both types together
      if (mDecoys.get(N)!=null) Nc.addAll(mDecoys.get(N).getData());
      if (mDecoys.get(C)!=null) Cc.addAll(mDecoys.get(C).getData());

      if (mNorms.get(N)!=null) Nc.addAll(mNorms.get(N).getData());
      if (mNorms.get(C)!=null) Cc.addAll(mNorms.get(C).getData());
    }
    if ((decoy!=null &&  decoy) && mDecoys!=null)
    {
      if (mDecoys.get(N)!=null) Nc.addAll(mDecoys.get(N).getData());
      if (mDecoys.get(C)!=null) Cc.addAll(mDecoys.get(C).getData());
    }
    if ((decoy!=null && !decoy) && mNorms!=null)
    {
      if (mNorms.get(N)!=null) Nc.addAll(mNorms.get(N).getData());
      if (mNorms.get(C)!=null) Cc.addAll(mNorms.get(C).getData());
    }
    // must have at least a value in each type
    if (Nc.size()==0) Nc.add(0d);
    if (Cc.size()==0) Cc.add(0d);

    // sort the scores
    Collections.sort(Nc); Collections.sort(Cc);

    // boot strap for the combined distribution
    Histogram combo = new Histogram("Bootstrapping estimation of the score distribution: "+decoy+"/"+N+"/"+C);
    Random      RND = new Random(System.nanoTime());
    for (int i=0; i<samples; i++)
    {
      double nt = Nc.get(RND.nextInt(Nc.size())), ct = Cc.get(RND.nextInt(Cc.size())),
            enz = sTrypicity[RND.nextInt(sTrypicity.length)]*trypicity;

//      if (nt<extension/10) nt += RND.nextDouble()*extension;
//      if (ct<extension/10) ct += RND.nextDouble()*extension;
      combo.add(nt+ct+RND.nextDouble()*jitter+enz);
    }
    combo.generate(25);
    Collections.sort(combo.getData(), Ordering.natural().reverse());
    Tools.dispose(Nc, Cc);

    return combo;
//    // prepare for the e-val calculation
//    return combo.fitEval((int )(samples*0.05d));
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
  public void dispose_intermediates(EnumMap<eType, Histogram> data)
  {
    if (Tools.isSet(data))
      for (eType t : data.keySet()) if (!t.equals(eType.bs)) data.get(t).dispose_data();
  }
  public void dispose(EnumMap<eType, Histogram>... data)
  {
    if (Tools.isSet(data))
      for (EnumMap<eType, Histogram> d : data)
        if (d!=null)
          for (eType t : d.keySet())
            if (d.get(t)!=null) d.get(t).dispose();
  }
  public void retainTops(String tag, EnumMap<eType, Histogram> data, eType t, int tops)
  {
    if (Tools.isSet(data) && data.get(t)!=null && data.get(t).getData()!=null && data.get(t).getData().size()>=tops)
      mBootsrapped.putAll(tag, data.get(t).getData().subList(0, tops));
  }
  public void dispose_intermediates()
  {
    // copy the necessary data before disposing the histogram
    mBootsrapped = TreeMultimap.create(Ordering.natural(), Ordering.natural().reverse());

    retainTops("Decoys", mDecoys, eType.bs, 25);
    retainTops("Norms",  mNorms,  eType.bs, 25);
    retainTops("Both",   mBoth,   eType.bs, 25);

    dispose(mDecoys, mNorms, mBoth);
  }
  public ScoreModel clone_intermediates(int tops)
  {
    ScoreModel cloned = new ScoreModel(getName());

    cloned.mBootsrapped = TreeMultimap.create(Ordering.natural(), Ordering.natural().reverse());

    cloned.retainTops("Decoys", mDecoys, eType.bs, tops);
    cloned.retainTops("Norms",  mNorms,  eType.bs, tops);
    cloned.retainTops("Both", mBoth, eType.bs, tops);

    return cloned;
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, mName);
    IOs.write(ds, mBaseline);
    IOs.write(ds, mCenter);
    IOs.write(ds, mBootstrapSize);
    IOs.writeStringDoubles(ds, mBootsrapped);

    write(ds, mDecoys);
    write(ds, mNorms);
    write(ds, mBoth);
    writeEnumDouble(ds, mOffsets);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mName          = IOs.read(ds, "");
    mBaseline      = IOs.read(ds, mBaseline);
    mCenter        = IOs.read(ds, 0d);
    mBootstrapSize = IOs.read(ds, 0);

    mBootsrapped   = (SortedSetMultimap )IOs.readStringDoubles(ds, TreeMultimap.create(Ordering.natural(), Ordering.natural().reverse()));

    mDecoys        = readEnumHistogram(ds);
    mNorms         = readEnumHistogram(ds);
    mBoth          = readEnumHistogram(ds);
    mOffsets       = readEnumDouble(   ds);
  }

  private void writeEnumDouble(DataOutput ds, EnumMap<eType, Double> data) throws IOException
  {
    IOs.write(ds, Tools.isSet(data) ? data.size() : 0);

    if (Tools.isSet(data))
      for (eType key : data.keySet())
      {
        IOs.write(ds, key.toString());
        IOs.write(ds, data.get(key));
      }
  }
  private EnumMap<eType, Double> readEnumDouble(DataInput ds) throws IOException
  {
    int n = IOs.read(ds, 0);
    if (n>0)
    {
      EnumMap<eType, Double> data = new EnumMap<>(eType.class);
      for (int i=0; i<n; i++)
        data.put(eType.valueOf(IOs.read(ds, "")), IOs.read(ds,0d));

      return data;
    }
    return null;
  }
  private <T extends Binary> void write(DataOutput ds, EnumMap<eType, T> data) throws IOException
  {
    IOs.write(ds, Tools.isSet(data) ? data.size() : 0);

    if (Tools.isSet(data))
      for (eType key : data.keySet())
      {
        IOs.write(ds, key.toString());
        IOs.write(ds, data.get(key));
      }
  }
  private EnumMap<eType, Histogram> readEnumHistogram(DataInput ds) throws IOException
  {
    int n = IOs.read(ds, 0);
    if (n>0)
    {
      EnumMap<eType, Histogram> data = new EnumMap<>(eType.class);
      for (int i=0; i<n; i++)
      {
        eType   key = eType.valueOf(IOs.read(ds, ""));
        Histogram H = new Histogram(); H.read(ds);
        data.put(key, H);
      }

      return data;
    }
    return null;
  }
}
