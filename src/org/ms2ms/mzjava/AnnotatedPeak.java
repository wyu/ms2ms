package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.Polarity;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.data.Binary;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/** A more elaborate class to describe the MS peak
 *
 * Created by wyu on 4/27/14.
 */
public class AnnotatedPeak extends Peak implements Binary
{
//  public static final String FREQUENCY = "Frequency";
//  public static final String SNR       = "SNR";
//  public static final String MZ        = "mz";
//  public static final String COUNTS    = "Counts";
//  public static final String OUTLIER   = "outlier";

  private boolean mIsOutlier=false;
  private int mVerifiedCharge=0, mIsotopes=1;
  private long mCounts;
  private double mSNR, mFreq, mOrigMz;
  private Map<String, Double> mAnnotations = null;
  private IonType ionType = IonType.unknown;

  public AnnotatedPeak() { super(); }
  public AnnotatedPeak(double mz, double ai) { super(mz, ai); }
  public AnnotatedPeak(double mz, double ai, int z) { super(mz, ai, z); }
  public AnnotatedPeak(double mz, double ai, int z, double snr)
  {
    super(mz, ai, z);
    setSNR(snr);
  }
  public AnnotatedPeak(double mz, double ai, int z, double snr, double freq)
  {
    super(mz, ai, z);
    setSNR(snr).setFrequency(freq);
  }
  public AnnotatedPeak(double mz, double ai, int z, float snr, IonType type)
  {
    super(mz, ai, z);
    setSNR(snr); ionType = type;
  }
  public AnnotatedPeak(AnnotatedPeak s)
  {
    super(s.getMz(), s.getIntensity(), s.getCharge());
    mAnnotations=s.mAnnotations; ionType=s.getIonType();
  }
  public boolean hasProperty(String k, Double s) { return mAnnotations!=null && Tools.equals(s, mAnnotations.get(k)); }
//  public boolean hasProperty(String k) { return mAnnotations!=null && mAnnotations.get(k)!=null; }

  public int     getIsotopes()         { return mIsotopes; }
  public int     getVerifiedCharge()   { return mVerifiedCharge; }
  public boolean isOutlier()           { return mIsOutlier; }
  public Double  getProperty(String k) { return mAnnotations!=null?mAnnotations.get(k):null; }
  public long    getCounts()           { return mCounts; }
  public double  getSNR()              { return mSNR; }
  public double  getFrequency()        { return mFreq; }
  public double  getOriginalMz()       { return mOrigMz; }
  public IonType getIonType()          { return ionType; }

  public Map<String, Double> getProperties() { return mAnnotations; }

  public AnnotatedPeak removeProperty(String s) { if (mAnnotations!=null) mAnnotations.remove(s); return this; }

  public AnnotatedPeak setCharge(        int z) { super.setMzAndCharge(getMz(), z); return this; }
  public AnnotatedPeak setVerifiedCharge(int z) { mVerifiedCharge=z; return this; }
  public AnnotatedPeak setIsotopes(      int s) { mIsotopes      =s; return this; }
  public AnnotatedPeak increIsotope()           { mIsotopes++;       return this; }
  public AnnotatedPeak isOutlier(    boolean s) { mIsOutlier     =s; return this; }
  public AnnotatedPeak setSNR(        double s) { mSNR           =s; return this; }
  public AnnotatedPeak setFrequency(  double s) { mFreq          =s; return this; }
  public AnnotatedPeak setOriginalMz( double s) { mOrigMz        =s; return this; }
  public AnnotatedPeak setCounts(       long s) { mCounts        =s; return this; }

  public AnnotatedPeak setProperty(String k, double s)
  {
    if (mAnnotations==null) mAnnotations=new HashMap<>(24);

    mAnnotations.put(k,s);
    return this;
  }
  public AnnotatedPeak setProperty(AnnotatedPeak s, String... ks)
  {
    // nothing to copy over
    if (s.mAnnotations==null) return this;
    // with an initial capacity
    if (mAnnotations==null) mAnnotations=new HashMap<>(24);

    if (!Tools.isSet(ks)) mAnnotations.putAll(s.mAnnotations);
    else
    {
      for (String k : ks) mAnnotations.put(k,s.getProperty(k));;
    }
    return this;
  }

  @Override
  public String toString()
  {
    String line = "m/z"+Tools.d2s(getMz(),2)+", %"+Tools.d2s(getIntensity(),4)+", z"+getCharge()+
                ", S/N"+Tools.d2s(getSNR(),1)+(getIsotopes()>1?"$"+getIsotopes():"");

    if (Tools.isSet(mAnnotations))
      for (String key : mAnnotations.keySet())
        if (mAnnotations.get(key)!=null)
        {
          int deci=1, t = key.indexOf('#'); Double val = mAnnotations.get(key);
          if (t>0)
          {
            deci = Integer.valueOf(key.substring(t+1)); key = key.substring(0, t);
          }
          line=Strs.extend(line, key+":"+Tools.d2s(val, deci), ", ");
        }

    return line;
  }
  public AnnotatedPeak clone()
  {
    AnnotatedPeak p = new AnnotatedPeak();

    p.setMzAndCharge(getMz(), getCharge());
    p.setIntensity(getIntensity());

    p.isOutlier(isOutlier()).setVerifiedCharge(getVerifiedCharge()).setCounts(getCounts());
    p.setSNR(getSNR()).setFrequency(getFrequency()).setOriginalMz(getOriginalMz());
    if (mAnnotations!=null) p.mAnnotations = new HashMap<>(mAnnotations);
    p.ionType=ionType;

    return p;
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds, getMz());IOs.write(ds, getCharge()); IOs.write(ds, getIntensity());

    IOs.write(ds, mIsOutlier);
    IOs.write(ds, mVerifiedCharge);
    IOs.write(ds, mIsotopes);
    IOs.write(ds, mCounts);
    IOs.write(ds, mSNR);
    IOs.write(ds, mFreq);
    IOs.write(ds, mOrigMz);
    IOs.writeStrDouble(ds, mAnnotations);
    IOs.write(ds, ionType.toString());
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    setMzAndCharge(IOs.read(ds, getMz()), IOs.read(ds, getCharge())); setIntensity(IOs.read(ds, getIntensity()));

    isOutlier(        IOs.read(ds, mIsOutlier));
    setVerifiedCharge(IOs.read(ds, mVerifiedCharge));
    setIsotopes(      IOs.read(ds, mIsotopes));
    setCounts(        IOs.read(ds, mCounts));
    setSNR(           IOs.read(ds, mSNR));
    setFrequency(     IOs.read(ds, mFreq));
    setOriginalMz(    IOs.read(ds, mOrigMz));
    mAnnotations =    IOs.readStrDouble(ds);
    ionType.valueOf(  IOs.read(ds,""));
  }
}
