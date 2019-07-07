package org.ms2ms.data.ms;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.utils.Copyable;
import org.ms2ms.Disposable;
import org.ms2ms.algo.Peaks;
import org.ms2ms.data.Binary;
import org.ms2ms.data.collect.ImmutableNavigableMap;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.IsotopePeakAnnotation;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;
import java.util.SortedMap;
import java.util.TreeMap;

/** A slim-down version PeakMatch
 * Created by yuw on 7/6/2019.
 */
public class PeakEntry implements Copyable<PeakEntry>, Comparable<PeakEntry>, Disposable, Binary
{
  private float mz=0.0f, intensity=0.0f, mSNR=0f, mFreq=0f, mScore=0f, mCalcMz=0f;
  private char charge=0;

//  // temp objects
//  static int[] tAccumSamples = Stats.fillIntArray(new int[255], 0);

  @Override
  public void dispose()
  {
  }

  public static class IntensityDesendComparator implements Comparator<PeakEntry>
  {
    public int compare(PeakEntry o1, PeakEntry o2) {
      return o1 != null && o2 != null ? Double.compare(o2.getIntensity(), o1.getIntensity()) : 0;
    }
  }

//  public PeakMatch() { this.polarity=Polarity.UNKNOWN; }
  public PeakEntry() { }
  public PeakEntry(double mz, double intensity)
  {
    this.setValues((float )mz, (float )intensity, 0);
  }

  public PeakEntry(float mz, float intensity, int z)
  {
//    this.polarity=Polarity.UNKNOWN;
    this.setValues(mz, intensity, z);
  }

  public PeakEntry(double mz, double ai, int z, float snr, float freq)
  {
    this.setValues((float )mz, (float )ai, z);
    setSNR(snr).setFrequency(freq);
  }
  public PeakEntry(PeakEntry peak)
  {
    this.mz       =peak.mz;
    this.intensity=peak.intensity;
    this.charge   =peak.charge;
    this.mFreq    =peak.mFreq;
    this.mSNR     =peak.mSNR;
    this.mScore   =peak.mScore;
    this.mCalcMz  =peak.mCalcMz;
  }
  public void setValues(float mz, float intensity, int charge)
  {
    this.setMzAndCharge(mz, charge);
    this.intensity=intensity;
  }
  public void invalidate()
  {
    mz=intensity=mFreq=mSNR=0;
    charge=0;
  }
  public boolean isValid() { return mz!=0 && intensity>0; }

  public float getMz()        { return this.mz; }
  public int getCharge()       { return this.charge; }
  public float getIntensity() { return this.intensity; }

  public void setMz(float mz) { this.mz=mz; }
  public void setCalcMz(float mz) { this.mCalcMz=mz; }
  public void setMzAndCharge(float mz, int z)
  {
    setMz(mz); setCharge(z);
  }

  public void setIntensity(float intensity)
  {
    this.intensity=intensity;
  }

  public PeakEntry copy() { return new PeakEntry(this); }

  public int compareTo(PeakEntry peak)
  {
    int cmp=Double.compare((float) this.getCharge(), (float) peak.getCharge());
    if (cmp!=0) {
      return cmp;
    } else {
      cmp=Double.compare(this.mz, peak.mz);
      return cmp!=0 ? cmp : Double.compare(this.intensity, peak.intensity);
    }
  }

  public float  getSNR()              { return mSNR; }
  public float  getFrequency()        { return mFreq; }
  public float  getScore()            { return mScore; }
  public float  getCalcMz()           { return mCalcMz; }

  public PeakEntry setCharge(int z) { charge=(char )z; return this; }
  public PeakEntry setSNR(float s) { mSNR           =s; return this; }
  public PeakEntry setFrequency(float s) { mFreq =s; return this; }
  public PeakEntry setScore(float s) { mScore         =s; return this; }

  public PeakEntry calcGapScore(int gap, float min_ppm, int[] tAccumSamples)
  {
    // no point to proceed...
    if (gap<=0 || gap>250) return setScore(0);

    // average AA mass: 115, largest - smallest: 129
    // match.getMz() is the actual mass deviation of the fragment
    int     bins = (int )(Math.log(2)/Math.log(1d+1E-6*Math.max(Math.abs(getMz()), min_ppm))),
        nsamples = 0, ntrials=Math.max(1,(int )Math.round(((gap-1)*115d+129d)*getFrequency()));

    // cumulative numbers of gaps
    if (tAccumSamples[gap]==0)
    {
      for (int i=1; i<=gap; i++)
        if (i<19 && nsamples<bins) nsamples+=Math.exp(Stats.ln_combination(19, i)); else break;

      tAccumSamples[gap]=nsamples;
    }
    else nsamples = tAccumSamples[gap];

    float score0=0;
    if (nsamples<bins/2)
    {
      score0 = (float )(-0.07491232f + 0.41163668f*Math.log((float )nsamples/(float )bins) + 0.40504996*Math.log(ntrials));
      if (score0 > -0.1) score0=0;
    }

    return setScore(score0);
  }

  @Override
  public String toString()
  {
    return isValid()? ("m/z"+Tools.d2s(getMz(),2)+", %"+Tools.d2s(getIntensity(),4)+", z"+getCharge()+
        ", S/N"+Tools.d2s(getSNR(),1)+(", scr="+Tools.d2s(getScore()*-10d,1))) : "---";
  }
  public PeakEntry clone()
  {
    PeakEntry p = new PeakEntry();

    p.setMzAndCharge(getMz(), getCharge());
    p.setIntensity(getIntensity());

    p.setSNR(getSNR()).setFrequency(getFrequency());
    p.mScore=mScore;

    return p;
  }

  public boolean equals(Object o)
  {
    if (this==o) {
      return true;
    } else if (o!=null&&this.getClass()==o.getClass())
    {
      PeakEntry peak=(PeakEntry) o;
      return Double.compare(peak.intensity, this.intensity)==0 &&
          Double.compare(peak.mz, this.mz)==0 && this.charge==peak.charge &&
       mSNR==peak.mSNR && mFreq==peak.mFreq ;
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

    return result;
  }

  /** static helper **/
//  public static <T extends PeakEntry> float centroid(Collection<T> points)
//  {
//    if (!Tools.isSet(points)) return 0;
//
//    float sumXY = 0, sumY = 0;
//    for (PeakEntry xy : points)
//    {
//      sumXY += xy.getMz() * xy.getIntensity();
//      sumY  += xy.getIntensity();
//    }
//    return sumY != 0 ? sumXY / sumY : 0;
//  }
//  public static <T extends PeakEntry> boolean hasNegativeIntensity(T... peaks)
//  {
//    if (Tools.isSet(peaks))
//      for (T peak : peaks) if (peak.getIntensity()<0) return true;
//
//    return false;
//  }
//  public static <T extends PeakEntry> float AbsIntensitySum(T... peaks)
//  {
//    if (!Tools.isSet(peaks)) return 0;
//
//    float sum = 0f;
//    for (T peak : peaks) sum += Math.abs(peak.getIntensity());
//
//    return sum;
//  }
  public static SortedMap<Double, PeakEntry> toPeaksWithExclusion(PeakList ms, Range<Double>... exclusion)
  {
    // walking thro the peaks and recording the matching peaks
    SortedMap<Double, PeakEntry> peaks = new TreeMap<>();

    int i=0, left,right;
    while (i<ms.size())
    {
      // save the c13 with negative intensity, WYU 20170318
      if (Peaks.hasC13(ms.getAnnotations(i)))
      {
        peaks.put(ms.getMz(i), new PeakEntry(ms.getMz(i), ms.getIntensity(i)*-1d));
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

      peaks.put(mz, new PeakEntry((float )ms.getMz(i), (float )ms.getIntensity(i), 0,0,
          (float )(Math.max(1,Peaks.countC12(ms, left, right))/(ms.getMz(right)-ms.getMz(left)))));
      // advance the pointer
      i++;
    }
    return peaks;
  }

//  public static float query4ai(ImmutableNavigableMap<PeakEntry> peaks, float m, OffsetPpmTolerance tol)
//  {
//    double err=tol.calcError(m), offset=tol.calcOffset(m), k0=m-err-offset, k1=m+err+offset;
//    float ai=0;
//    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);
//
//    if (start>=0)
//    {
//      int j0=-1, j1=-1;
//      for (int k=start; k<peaks.getKeys().length; k++)
//      {
//        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
//        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
//      }
//      if (j0>=0 && j1>j0)
//        for (int j=j0; j<j1; j++)
//          ai+=Math.abs(peaks.getVals()[j].getIntensity());
//    }
//
//    return ai;
//  }
//  public static Collection<Double> query4keys(ImmutableNavigableMap<PeakEntry> peaks, float m, OffsetPpmTolerance tol, Collection<Double> keys)
//  {
//    double k0=tol.getMin(m), k1=tol.getMax(m);
//    int i0=Math.max(0, peaks.index(k0)), start=peaks.start(i0);
//
//    if (start>=0)
//    {
//      int j0=-1, j1=-1;
//      for (int k=start; k<peaks.getKeys().length; k++)
//      {
//        if (j0==-1 && peaks.getKeys()[k]>=k0)   j0=k;
//        if (          peaks.getKeys()[k]> k1) { j1=k; break; }
//      }
//      if (j0>=0 && j1>j0)
//        for (int j=j0; j<j1; j++)
//          keys.add(peaks.getKeys()[j]);
//    }
//
//    return keys;
//  }

  public static boolean query4index(ImmutableNavigableMap<PeakEntry> peaks, double k0, double k1, Collection<Integer> index)
  {
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
      {
        for (int j=j0; j<j1; j++) index.add(j);
        return true;
      }
    }

    return false;
  }

  // negative: OK with negative intensity?
  public static IsoEnvelope query4isotope(ImmutableNavigableMap<PeakEntry> peaks, float m, OffsetPpmTolerance tol, boolean negative)
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
        pk = new IsoEnvelope(m0/(float )(j1-j0),ai, 0);
        pk.setScore(fr/(float) (j1-j0));
      }
    }

    return pk;
  }
//  public static int query4counts(ImmutableNavigableMap<PeakEntry> peaks, OffsetPpmTolerance tol, float m)
//  {
//    double err=tol.calcError(m), offset=tol.calcOffset(m);
//    return peaks.query4counts(m-err-offset, m+err+offset);
//  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.write(ds,mz);
    IOs.write(ds,intensity);
    IOs.write(ds,mSNR);
    IOs.write(ds,mFreq);
    IOs.write(ds,mScore);
    IOs.write(ds,charge);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mz       =IOs.read(ds, 0f);
    intensity=IOs.read(ds, 0f);
    mSNR     =IOs.read(ds, 0f);
    mFreq    =IOs.read(ds, 0f);
    mScore   =IOs.read(ds, 0f);
    charge   =IOs.read(ds, '0');
  }
}