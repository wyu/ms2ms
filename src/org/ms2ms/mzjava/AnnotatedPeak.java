package org.ms2ms.mzjava;

import org.expasy.mzjava.core.ms.spectrum.Peak;

/** A more elaborate class to descrbe the MS peak
 *
 * Created by wyu on 4/27/14.
 */
public class AnnotatedPeak extends Peak
{
  private double SNR;
  public AnnotatedPeak() { super(); }
  public AnnotatedPeak(double mz, double ai, int z, double snr)
  {
    super(mz, ai, z);
    SNR=snr;
  }

  public double getSNR() { return SNR; }
}
