package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.utils.Settings;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 10/5/14
 * Time: 3:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class LcSettings extends Settings
{
  public static final String KY_RT_TOL = "RTTol";      // likely variation
  public static final String KY_LC_WIDTH = "Expected LC peak width";
  public static final String KY_RT_SPAN = "RT Window";  // max variation
  public static final String KY_SMOOTHRT  = "RT smoothed";
  public static final String KY_USE_IRT  = "Use iRT";
  public static final String KY_PEAK_EXCL  = "Min peak exclusivity";
  public static final String KY_PEAK_APEX  = "Min peak apex";

  public static LcSettings nLC = null;
  static
  {
    nLC = new LcSettings().setMinApex(10000).setMinPeakExclusivity(50f).setPeakWidth(0.15f).setSpan(5f);
  }


  public Tolerance getTolerance()          { return properties!=null?(Tolerance )properties.get(KY_RT_TOL):null; }
  public float     getSpan()               { return getFloat(KY_RT_SPAN); }
  public float     getPeakWidth()          { return getFloat(KY_LC_WIDTH); }
  public float     getMinPeakExclusivity() { return getFloat(KY_PEAK_EXCL); }
  public float     getMinApex()            { return getFloat(KY_PEAK_APEX); }

  public boolean toSmoothRT() { return getBoolean(KY_SMOOTHRT); }
  public boolean useiRT() { return getBoolean(KY_USE_IRT); }

  public LcSettings setSpan(              float s) { set(KY_RT_SPAN,   s); return this; }
  public LcSettings setPeakWidth(         float s) { set(KY_LC_WIDTH,  s); return this; }
  public LcSettings setMinPeakExclusivity(float s) { set(KY_PEAK_EXCL, s); return this; }
  public LcSettings setMinApex(           float s) { set(KY_PEAK_APEX, s); return this; }

  public LcSettings toSmoothRT(         boolean s) { set(KY_SMOOTHRT, s); return this; }
  public LcSettings useiRT(             boolean s) { set(KY_USE_IRT, s); return this; }

  public Settings setRt(float tol, float window)
  {
    properties.put(KY_RT_TOL,  new AbsoluteTolerance(tol));
    properties.put(KY_RT_SPAN, new AbsoluteTolerance(window));
    return this;
  }

  @Override
  public LcSettings clone() throws CloneNotSupportedException
  {
    return (LcSettings )super.clone();
  }
}
