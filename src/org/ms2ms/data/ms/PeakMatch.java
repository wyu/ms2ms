package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.peaklist.Polarity;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.utils.Copyable;
import org.ms2ms.Disposable;
import org.ms2ms.algo.Peaks;
import org.ms2ms.data.Binary;
import org.ms2ms.data.collect.ImmutableNavigableMap;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.*;

/** A simplified, riskier version Peak
 * Created by yuw on 3/26/17.
 */
public class PeakMatch extends PeakFragmentMatch implements Copyable<PeakMatch>, Comparable<PeakMatch>
{
  private double mz=0.0d, intensity=0.0d, mass=0.0d, mSNR=0d, mFreq=0d, mOrigMz=0d, mCalcMz=0d, mScore=0d, mRawAI=0d;
  private Float mz_low, mz_high;
//  private Polarity polarity;
//  private int[] chargeList=new int[0];
  private int charge=0;

  private boolean mIsOutlier=false;
  private int mVerifiedCharge=0, mIsotopes=1, mIndex;
  private long mCounts=1;
//  private Map<String, Double> mAnnotations = null;
  private IonType ionType = IonType.unknown ;

  // temp objects
//  static int[] tAccumSamples = Stats.fillIntArray(new int[255], 0);

  @Override
  public void dispose()
  {
//    polarity  =null;
//    chargeList=null;
    ionType   =null;
  }

  public static class IntensityDesendComparator implements Comparator<PeakMatch>
  {
    public int compare(PeakMatch o1, PeakMatch o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getIntensity(), o1.getIntensity()) : 0;
    }
  }

//  public PeakMatch() { this.polarity=Polarity.UNKNOWN; }
  public PeakMatch() { }
  public PeakMatch(double mz, double intensity)
  {
    this.setValues(mz, intensity, 0);
  }

  public PeakMatch(double mz, double intensity, int z)
  {
//    this.polarity=Polarity.UNKNOWN;
    this.setValues(mz, intensity, z);
  }

//  public static PeakMatch noIntensity(double mz, int z)
//  {
//    return new PeakMatch(mz, 0.0D, z);
//  }

//  public PeakMatch(double mz, double ai, int z, double snr)
//  {
//    this.setValues(mz, ai, z);
//    setSNR(snr);
//  }
  public PeakMatch(double mz, double ai, int z, double snr, double freq)
  {
    this.setValues(mz, ai, z);
    setSNR(snr).setFrequency(freq);
  }
//  public PeakMatch(double mz, double ai, int z, float snr, IonType type)
//  {
//    this.setValues(mz, ai, z);
//    setSNR(snr); ionType = type;
//  }
  public PeakMatch(PeakMatch peak)
  {
    this.mz       =peak.mz;
    this.mass     =peak.mass;
    this.intensity=peak.intensity;
    this.charge   =peak.charge;
    this.mass     =peak.mass;
    this.mFreq    =peak.mFreq;
    this.mSNR     =peak.mSNR;

    ionType=peak.getIonType();
  }
  public PeakMatch(PeakEntry peak)
  {
    this.mz       =peak.getMz();
    this.intensity=peak.getIntensity();
    this.charge   =peak.getCharge();
    this.mFreq    =peak.getFrequency();
    this.mSNR     =peak.getSNR();
    this.mCalcMz  =peak.getCalcMz();
  }
  public void setValues(double mz, double intensity, int charge)
  {
    this.setMzAndCharge(mz, charge);
    this.intensity=intensity;
  }
  public void invalidate()
  {
    mz=mass=intensity=mass=mFreq=mSNR=0;
    charge=0;

    ionType=IonType.unknown;
  }
  public boolean isValid() { return mz!=0 && intensity>0; }

  public double getMz()        { return this.mz; }
  public Float getMzLow()     { return this.mz_low; }
  public Float getMzHigh()    { return this.mz_high; }
  public int getCharge()       { return this.charge; }
  public double getIntensity() { return this.intensity; }
  public double getRawAI() { return mRawAI; }

  public void setMz(double mz) { this.mz=mz; }
  public PeakMatch setRawAI(double s) { mRawAI=s; return this; }
  public PeakMatch setMzExpectedBound(OffsetPpmTolerance tol)
  {
    double err=tol.calcError(mz), offset=tol.calcOffset(mz);
    this.mz_low=(float )(mz-err+offset); this.mz_high=(float )(mz+err+offset);
    return this;
  }

  public void setMzAndCharge(double mz, int z)
  {
    setMz(mz); setCharge(z);
  }

  public void setIntensity(double intensity)
  {
    this.intensity=intensity;
  }

  public PeakMatch copy() { return new PeakMatch(this); }

  public int compareTo(PeakMatch peak)
  {
    int cmp=Double.compare((double) this.getCharge(), (double) peak.getCharge());
    if (cmp!=0) {
      return cmp;
    } else {
      cmp=Double.compare(this.mz, peak.mz);
      return cmp!=0 ? cmp : Double.compare(this.intensity, peak.intensity);
    }
  }

  public int     getIndex()            { return mIndex; }
  public int     getIsotopes()         { return mIsotopes; }
  public int     getVerifiedCharge()   { return mVerifiedCharge; }
  public boolean isOutlier()           { return mIsOutlier; }
  public long    getCounts()           { return mCounts; }
  public double  getSNR()              { return mSNR; }
  public double  getFrequency()        { return mFreq; }
  public double  getOriginalMz()       { return mOrigMz; }
  public double  getCalcMz()           { return mCalcMz; }
  public double  getScore()            { return mScore; }
  public IonType getIonType()          { return ionType; }
  public boolean isIonType(IonType s)  { return Tools.equals(s, ionType); }

  public PeakMatch setIndex(         int z) { mIndex=z; return this; }
  public PeakMatch setCharge(        int z) { charge=z; return this; }
  public PeakMatch setVerifiedCharge(int z) { mVerifiedCharge=z; return this; }
  public PeakMatch setIsotopes(      int s) { mIsotopes      =s; return this; }
  public PeakMatch increIsotope()           { mIsotopes++;       return this; }
  public PeakMatch isOutlier(    boolean s) { mIsOutlier     =s; return this; }
  public PeakMatch setSNR(        double s) { mSNR           =s; return this; }
  public PeakMatch setFrequency(  double s)
  {
    if (s==0)
      System.out.println();
    mFreq          =s; return this;
  }
  public PeakMatch setOriginalMz( double s) { mOrigMz        =s; return this; }
  public PeakMatch setCalcMz(     double s) { mCalcMz        =s; return this; }
  public PeakMatch setScore(      double s) { mScore         =s; return this; }
  public PeakMatch setCounts(       long s) { mCounts        =s; return this; }
  public PeakMatch setIonType(   IonType s) { ionType        =s; return this; }

//  public PeakMatch calcGapScore(int gap, double min_ppm)
//  {
//    // no point to proceed...
//    if (gap<=0 || gap>250) return setScore(0);
//
//    // average AA mass: 115, largest - smallest: 129
//    // match.getMz() is the actual mass deviation of the fragment
//    int     bins = (int )(Math.log(2)/Math.log(1d+1E-6*Math.max(Math.abs(getMz()), min_ppm))),
//        nsamples = 0, ntrials=Math.max(1,(int )Math.round(((gap-1)*115d+129d)*getFrequency()));
//
//    // cumulative numbers of gaps
//    if (tAccumSamples[gap]==0)
//    {
//      for (int i=1; i<=gap; i++)
//        if (i<19 && nsamples<bins) nsamples+=Math.exp(Stats.ln_combination(19, i)); else break;
//
//      tAccumSamples[gap]=nsamples;
//    }
//    else nsamples = tAccumSamples[gap];
//
//    double score0=0;
//    if (nsamples<bins/2)
//    {
//      score0 = -0.07491232 + 0.41163668*Math.log((double )nsamples/(double )bins) + 0.40504996*Math.log((double )ntrials);
//      if (score0 > -0.1) score0=0;
//    }
//
//    return setScore(score0);
//  }

  @Override
  public String toString()
  {
    return isValid()? ("m/z"+Tools.d2s(getMz(),2)+", %"+Tools.d2s(getIntensity(),4)+", z"+getCharge()+
        ", S/N"+Tools.d2s(getSNR(),1)+(getIsotopes()>1?"$"+getIsotopes():"")+(", scr="+Tools.d2s(getScore()*-10d,1))) : "---";
  }
  public PeakMatch clone()
  {
    PeakMatch p = new PeakMatch();

    p.setMzAndCharge(getMz(), getCharge());
    p.setIntensity(getIntensity());

    p.isOutlier(isOutlier()).setVerifiedCharge(getVerifiedCharge()).setCounts(getCounts());
    p.setSNR(getSNR()).setFrequency(getFrequency()).setOriginalMz(getOriginalMz());
    p.ionType=ionType;
    p.mCalcMz=mCalcMz; p.mScore=mScore;

    return p;
  }

  public boolean equals(Object o)
  {
    if (this==o) {
      return true;
    } else if (o!=null&&this.getClass()==o.getClass())
    {
      PeakMatch peak=(PeakMatch) o;
      return Double.compare(peak.intensity, this.intensity)==0 &&
          Double.compare(peak.mz, this.mz)==0 && this.charge==peak.charge &&
      mIsOutlier==peak.mIsOutlier && mVerifiedCharge==peak.mVerifiedCharge && mIsotopes==peak.mIsotopes &&
      mCounts==peak.mCounts && mSNR==peak.mSNR && mFreq==peak.mFreq && mOrigMz==peak.mOrigMz &&
          ionType==peak.ionType/* && Tools.equals(mAnnotations, peak.mAnnotations)*/;
    } else {
      return false;
    }
  }

  public int hashCode()
  {
    long temp=this.mz!=0.0D ? Double.doubleToLongBits(this.mz) : 0L;
    int result=(int) (temp^temp>>>32);
    temp=this.intensity!=0.0D ? Double.doubleToLongBits(this.intensity) : 0L;
    result=31*result+(int) (temp^temp>>>32);
    result=31*result+this.charge;
//    result=31*result+this.polarity.hashCode();

//    result=31*result+mIsOutlier==peak.mIsOutlier && mVerifiedCharge==peak.mVerifiedCharge && mIsotopes==peak.mIsotopes &&
//        mCounts==peak.mCounts && mSNR==peak.mSNR && mFreq==peak.mFreq && mOrigMz==peak.mOrigMz &&
//        ionType==peak.ionType && Tools.equals(mAnnotations, peak.mAnnotations);

    return result;
  }

  /** static helper **/
  public static <T extends PeakMatch> double centroid(Collection<T> points)
  {
    if (!Tools.isSet(points)) return 0;

    double sumXY = 0, sumY = 0;
    for (PeakMatch xy : points)
    {
      sumXY += xy.getMz() * xy.getIntensity();
      sumY  += xy.getIntensity();
    }
    return sumY != 0 ? sumXY / sumY : 0;
  }
  public static <T extends PeakMatch> boolean hasNegativeIntensity(T... peaks)
  {
    if (Tools.isSet(peaks))
      for (T peak : peaks) if (peak.getIntensity()<0) return true;

    return false;
  }
  public static <T extends PeakMatch> double AbsIntensitySum(T... peaks)
  {
    if (!Tools.isSet(peaks)) return 0;

    double sum = 0d;
    for (T peak : peaks) sum += Math.abs(peak.getIntensity());

    return sum;
  }
  public static SortedMap<Double, PeakMatch> toPeaksWithExclusion(PeakList ms, Range<Double>... exclusion)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double, PeakMatch> peaks = new TreeMap<>();

    int i=0, left,right;
    while (i<ms.size())
    {
      double raw_ai = Tools.isSet(ms.getAnnotations(i)) ? ((IsotopePeakAnnotation)Tools.front(ms.getAnnotations(i))).getIntensity():0d;
      // save the c13 with negative intensity, WYU 20170318
      if (Peaks.hasC13(ms.getAnnotations(i)))
      {
        peaks.put(ms.getMz(i), new PeakMatch(ms.getMz(i), ms.getIntensity(i)*-1d).setRawAI(raw_ai));
        i++; continue;
      }

      double mz=ms.getMz(i);

      // check for exclusion
      if (Tools.isSet(exclusion))
      {
        boolean found=false;
        for (Range<Double> excluded : exclusion)
          if (excluded.contains(mz)) { found=true; break; }
        if (found) { i++; continue; }
      }

      // estimates the local frequency
      left=i-10;right=i+10;
      if      (left <0)           { left =0;          right=Math.min(left+20,ms.size()-1); }
      else if (right>ms.size()-1) { right=ms.size()-1; left=Math.max(right-20, 0); }

      peaks.put(mz, new PeakMatch(ms.getMz(i), ms.getIntensity(i), 0,0,
          Math.max(1,Peaks.countC12(ms, left, right))/(ms.getMz(right)-ms.getMz(left))).setRawAI(raw_ai));
      // advance the pointer
      i++;
    }
    return peaks;
  }

  public static double query4ai(ImmutableNavigableMap<PeakMatch> peaks, double m, OffsetPpmTolerance tol)
  {
    double err=tol.calcError(m), offset=tol.calcOffset(m), k0=m-err-offset, k1=m+err+offset, ai=0;
    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);

    if (start>=0)
    {
      int j0=-1, j1=-1;
      for (int k=start; k<peaks.getKeys().length; k++)
      {
        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
      }
      if (j0>=0 && j1>j0)
        for (int j=j0; j<j1; j++)
          ai+=Math.abs(peaks.getVals()[j].getIntensity());
    }

    return ai;
  }
  public static Collection<Double> query4keys(ImmutableNavigableMap<PeakMatch> peaks, double m, OffsetPpmTolerance tol, Collection<Double> keys)
  {
    double k0=tol.getMin(m), k1=tol.getMax(m);
    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);

    if (start>=0)
    {
      int j0=-1, j1=-1;
      for (int k=start; k<peaks.getKeys().length; k++)
      {
        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
      }
      if (j0>=0 && j1>j0)
        for (int j=j0; j<j1; j++)
          keys.add(peaks.getKeys()[j]);
    }

    return keys;
  }

  // negative: OK with negative intensity?
  public static IsoEnvelope query4isotope(ImmutableNavigableMap<PeakMatch> peaks, double m, OffsetPpmTolerance tol, boolean negative)
  {
    double err=tol.calcError(m), offset=tol.getOffset(m), k0=m-err-offset, k1=m+err+offset;
    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);

    IsoEnvelope pk = null;
    if (start>=0)
    {
      int j0=-1, j1=-1;
      for (int k=start; k<peaks.getKeys().length; k++)
      {
        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
      }
      if (j0>=0 && j1>j0)
      {
        double m0=0, ai=0, fr=0d;

        for (int j=j0; j<j1; j++)
        {
          if (!negative && peaks.getVals()[j].getIntensity()<0) return null;

          m0+=peaks.getVals()[j].getMz();
          fr+=peaks.getVals()[j].getFrequency();
          ai+=Math.abs(peaks.getVals()[j].getIntensity());
        }
        pk = new IsoEnvelope(m0/(double )(j1-j0),ai, 0);
        pk.setScore(fr/(double) (j1-j0));
      }
    }

    return pk;
  }
  public static int query4counts(ImmutableNavigableMap<PeakMatch> peaks, OffsetPpmTolerance tol, double m)
  {
    double err=tol.calcError(m), offset=tol.calcOffset(m);
    return peaks.query4counts(m-err-offset, m+err+offset);
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds,mz);
    IOs.write(ds,intensity);
    IOs.write(ds,mass);
    IOs.write(ds,mSNR);
    IOs.write(ds,mFreq);
    IOs.write(ds,mOrigMz);
    IOs.write(ds,mCalcMz);
    IOs.write(ds,mScore);
    IOs.write(ds,mz_low);
    IOs.write(ds,mz_high);
    IOs.write(ds,charge);

    IOs.write(ds,mIsOutlier);
    IOs.write(ds,mVerifiedCharge);
    IOs.write(ds,mIsotopes);
    IOs.write(ds,mIndex);
    IOs.write(ds,mCounts);
    IOs.write(ds,ionType.toString());
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mz       =IOs.read(ds, 0d);
    intensity=IOs.read(ds, 0d);
    mass     =IOs.read(ds, 0d);
    mSNR     =IOs.read(ds, 0d);
    mFreq    =IOs.read(ds, 0d);
    mOrigMz  =IOs.read(ds, 0d);
    mCalcMz  =IOs.read(ds, 0d);
    mScore   =IOs.read(ds, 0d);
    mz_low   =IOs.read(ds, 0F);
    mz_high  =IOs.read(ds, 0F);
    charge   =IOs.read(ds, 0);

    mIsOutlier     =IOs.read(ds, false);
    mVerifiedCharge=IOs.read(ds, 0);
    mIsotopes      =IOs.read(ds, 0);
    mIndex         =IOs.read(ds, 0);
    mCounts        =IOs.read(ds, 0L);
    ionType.valueOf(IOs.read(ds,""));
  }
}