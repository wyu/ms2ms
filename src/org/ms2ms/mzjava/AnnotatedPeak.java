package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.utils.Tools;

/** A more elaborate class to describe the MS peak
 *
 * Created by wyu on 4/27/14.
 */
public class AnnotatedPeak extends Peak
{
  private double SNR;
  private IonType ionType ;

  public AnnotatedPeak() { super(); }
  public AnnotatedPeak(double mz, double ai, int z, double snr)
  {
    super(mz, ai, z);
    SNR=snr;
  }
  public AnnotatedPeak(double mz, double ai, int z, double snr, IonType type)
  {
    super(mz, ai, z);
    SNR=snr; ionType = type;
  }
  public AnnotatedPeak(AnnotatedPeak s)
  {
    super(s.getMz(), s.getIntensity(), s.getCharge());
    SNR=s.getSNR(); ionType=s.getIonType();
  }
  public double getSNR() { return SNR; }
  public IonType getIonType() { return ionType; }

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
