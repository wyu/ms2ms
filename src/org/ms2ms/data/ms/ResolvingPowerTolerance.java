//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by Fernflower decompiler)
//

package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Writer;

/**
 * Steve: Further to our meeting yesterday, I was checking some of the mass windows Skyline uses. When I tell Skyline we have 30,000 resolution at 200 m/z, Skyline uses that information to define ppm windows of increasing size across the m/z range. E.g.:
 * ~30 ppm window at 440 m/z
 * ~46 ppm window at 764 m/z
 * ~65 ppm window at 1528 m/z
 *
 * I think this is best practice, as the resolution and peak width will vary across the m/z range. I’m not sure exactly how Skyline calculates the window size- haven’t thought about this for years! Do any of you know the right equation for this?
 *
 Wen: The discussion thread below seems to be quite relevant to our question. The discussion also touched on the “High selectivity extraction” option in Skyline, where only 1x FWHM is used.
 https://skyline.ms/announcements/home/support/thread.view?rowId=26327
 https://skyline.ms/announcements/home/support/thread.view?entityId=a6ebc95f-05ec-1033-bad2-da202582f4ca&_anchor=22458#row:22458

 The full width of the extraction window of a particular mz is: 2 * mz * sqrt(mz) / sqrt(referenceMz) / resolution

 We can convert the formula to: ppm = 2*10E6 * sqrt(r) / R; where “R” is the resolving power, and “r” is the ratio of the m/z in consideration to the m/z where the resolving power is defined.  With this formula, the m vs ppm in the example Steve described would be:

 The 'resolving power' in Orbitrap is proportional to square root of [(kq.m-1)], which is the kinetic energy of already formed ion along the axis of the analyzer. This 'frequency resolving power' is equal to twice the mass resolving power. The mass resolving power is determined by the FWHM of the analyte MS peak. When the FWHM decreases, i.e. increase the resolution, that the inaccuracy of the determination of the MS peak 'area' is reduced. This is especially important for overlaping MS peaks, i.e. very close m/z. Therefore with high-to-ultra-high resolving power which is in the capacity of Orbitrap is achieved not only accurate identification of analytes in a complex multicomponent mixture, but as well a highly accurate quantitation (i.e intensities determination of closely disposed MS peaks).
 This matter is very well explained with corresponding equations and graphs in the following paper:
 1. R. Perry, R. Cooks, R. Noll, Orbitrap mass spectrometry: instrumentation, ion motion and application, Mass Spectrometry Reviews, 2008, 27, 661– 699
 */
public class ResolvingPowerTolerance extends Tolerance {
  private double mResolvingPower=30000d;
  private double mNominalMz     =200d;

  public ResolvingPowerTolerance(double rp) {
    this.mResolvingPower = rp; this.mNominalMz=200d;
  }
  public ResolvingPowerTolerance(double rp, double mz) {
    this.mResolvingPower = rp; this.mNominalMz=mz;
  }

  public boolean withinTolerance(double expected, double actual) {
    double error = this.calcError(expected);
    double min = expected - error;
    double max = expected + error;
    return actual >= min && actual <= max;
  }

  public double getError(double expectedMr, double actualMr) {
    double delta = actualMr / expectedMr;
    return actualMr > expectedMr ? delta / expectedMr * 1000000.0D : -(delta / expectedMr) * 1000000.0D;
  }

  public Location check(double expected, double actual) {
    double error = this.calcError(expected);
    double min = expected - error;
    double max = expected + error;
    if (actual < min) {
      return Location.SMALLER;
    } else {
      return actual > max ? Location.LARGER : Location.WITHIN;
    }
  }

  public double getMin(double mz) {
    double error = this.calcError(mz);
    return mz - error;
  }

  public double getMax(double mz) {
    double error = this.calcError(mz);
    return mz + error;
  }

  private double calcError(double expectedMass) {
    // The full width of the extraction window of a particular mz is: 2 * mz * sqrt(mz) / sqrt(referenceMz) / resolution
    // for example, RP=30000 @ m/z 200
    // ~30 ppm window at 440 m/z
    // ~46 ppm window at 764 m/z
    // ~65 ppm window at 1528 m/z
    // to capture the obs RP to ppm relationship, the formula will be changed to:
    return (expectedMass/mResolvingPower) * Math.sqrt(expectedMass/(2d*mNominalMz));
  }
  public void print(PrintStream w, double... mzs) throws IOException
  {
    w.print("Resolving Power: "+mResolvingPower+" @ m/z"+mNominalMz+"\n");
    w.print("mz\tmzL\tmzU\tErr\tppm\n");
    if (Tools.isSet(mzs))
       for (double mz : mzs)
       {
         w.print(mz+"\t"+Tools.d2s(getMin(mz),4)+"\t"+Tools.d2s(getMax(mz),4)+"\t"+Tools.d2s(calcError(mz),4)+"\t"+Tools.d2s(1E6*calcError(mz)/mz, 1)+"\n");
       }
  }
}
