package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.HashMap;
import java.util.Map;

/** A more elaborate class to describe the MS peak
 *
 * Created by wyu on 4/27/14.
 */
public class AnnotatedPeak extends Peak
{
//  public static final String FREQUENCY = "Frequency";
//  public static final String SNR       = "SNR";
//  public static final String MZ        = "mz";
//  public static final String COUNTS    = "Counts";
//  public static final String OUTLIER   = "outlier";

  private boolean mIsOutlier=false;
  private int mVerifiedCharge=0;
  private long mCounts;
  private double mSNR, mFreq, mOrigMz;
  private Map<String, Double> mAnnotations = null;
  private IonType ionType ;

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

  public int getVerifiedCharge() { return mVerifiedCharge; }
  public boolean isOutlier() { return mIsOutlier; }
  public Double getProperty(String k) { return mAnnotations.get(k); }
  public long getCounts() { return mCounts; }
  public double getSNR() { return mSNR; }
  public double getFrequency() { return mFreq; }
  public double getOriginalMz() { return mOrigMz; }
  public IonType getIonType() { return ionType; }

  public AnnotatedPeak removeProperty(String s) { if (mAnnotations!=null) mAnnotations.remove(s); return this; }

  public AnnotatedPeak setCharge(int z) { super.setMzAndCharge(getMz(), z); return this; }
  public AnnotatedPeak setVerifiedCharge(int z) { mVerifiedCharge=z; return this; }
  public AnnotatedPeak isOutlier(boolean s) { mIsOutlier=s; return this; }
  public AnnotatedPeak setSNR(        double s) { mSNR=s; return this; }
  public AnnotatedPeak setFrequency(  double s) { mFreq=s; return this; }
  public AnnotatedPeak setOriginalMz(double s) { mOrigMz=s; return this; }
  public AnnotatedPeak setCounts(long s) { mCounts=s; return this; }
  public AnnotatedPeak setProperty(   String k, double s)
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
    String line = "m/z"+ Tools.d2s(getMz(),2)+", %"+Tools.d2s(getIntensity(),4)+", z"+getCharge()+", S/N"+Tools.d2s(getSNR(),1);

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
    return new AnnotatedPeak(this);
  }
}
