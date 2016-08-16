package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.utils.Tools;

import java.util.HashMap;
import java.util.Map;

/** A more elaborate class to describe the MS peak
 *
 * Created by wyu on 4/27/14.
 */
public class AnnotatedPeak extends Peak
{
  public static final String FREQUENCY = "Frequency";
  public static final String SNR       = "SNR";
  public static final String MZ        = "mz";

//  private double SNR;
  private Map<String, Double> mAnnotations = new HashMap<>();
  private IonType ionType ;

  public AnnotatedPeak() { super(); }
  public AnnotatedPeak(double mz, double ai)
  {
    super(mz, ai);
  }
  public AnnotatedPeak(double mz, double ai, int z, double snr)
  {
    super(mz, ai, z);
    setSNR(snr);
  }
  public AnnotatedPeak(double mz, double ai, int z, double snr, double freq)
  {
    super(mz, ai, z);
    setSNR(snr).setProperty(FREQUENCY, freq);
  }
  public AnnotatedPeak(double mz, double ai, int z, double snr, IonType type)
  {
    super(mz, ai, z);
    setSNR(snr); ionType = type;
  }
  public AnnotatedPeak(AnnotatedPeak s)
  {
    super(s.getMz(), s.getIntensity(), s.getCharge());
    mAnnotations=s.mAnnotations; ionType=s.getIonType();
  }
  public Double getProperty(String k) { return mAnnotations.get(k); }
  public Double getSNR() { return getProperty(SNR); }
  public IonType getIonType() { return ionType; }

  public AnnotatedPeak setSNR(Double s) { return setProperty(SNR, s); }
  public AnnotatedPeak setProperty(String k, double s) { mAnnotations.put(k,s); return this; }

  @Override
  public String toString()
  {
    return "m/z"+ Tools.d2s(getMz(),2)+", %"+Tools.d2s(getIntensity(),4)+", z"+getCharge()+", S/N"+Tools.d2s(getSNR(),1);
  }
  public AnnotatedPeak clone()
  {
    return new AnnotatedPeak(this);
  }
}
