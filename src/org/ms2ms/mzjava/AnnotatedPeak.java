package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.core.ms.spectrum.Peak;

/** A more elaborate class to descrbe the MS peak
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

  public AnnotatedPeak clone()
  {
    return new AnnotatedPeak(this);
  }
}
