package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.utils.Settings;
import org.ms2ms.utils.Strs;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 10/5/14
 * Time: 3:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class LcSettings extends Settings
{
  public enum calibration  { loess, SG5, pt2 }

  public static final String KY_RT_TOL = "RTTol";      // likely variation
  public static final String KY_LC_WIDTH = "Expected LC peak width at half-height";
  public static final String KY_RT_SPAN = "RT Window";  // max variation
  public static final String KY_SMOOTHRT  = "RT smoothed";
  public static final String KY_USE_IRT  = "Use iRT";
  public static final String KY_PEAK_EXCL  = "Min peak exclusivity";
  public static final String KY_PEAK_APEX  = "Min peak apex";
  public static final String KY_APEX_PTS  = "Min peak apex points";
  public static final String KY_PEAK_BASE  = "Centroid bound by % apex";
  public static final String KY_PEAK_MULTIPLE  = "outer bound as multiple of LC width";
  public static final String KY_GRID_SIZE  = "grid size";
  public static final String KY_RTCAL  = "RT calibration method";
  public static final String KY_BANDWIDTH  = "Loess Bandwidth";
  public static final String KY_QUAN_SPAN  = "RT Window during quantitation";
  public static final String KY_QUAN_OFFSET  = "RT offset during quantitation";
  public static final String KY_SNR  = "min SNR to qualify";

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
  public float     getBaseRI()             { return getFloat(KY_PEAK_BASE); }
  public float     getOuterMultiple()      { return getFloat(KY_PEAK_MULTIPLE); }
  public float     getBandwidth()          { return getFloat(KY_BANDWIDTH); }
  public float     getQuanSpan()           { return getFloat(KY_QUAN_SPAN); }
  public float     getQuanOffset()         { return getFloat(KY_QUAN_OFFSET, 0f); }
  public float     getMinSNR()             { return getFloat(KY_SNR); }

  public int       getApexPts()            { return getInteger(KY_APEX_PTS); }
  public int       getGridSize()           { return getInteger(KY_GRID_SIZE); }

  public boolean   isCalMethod(     calibration... s)
  {
    for (calibration c : s)
      if (Strs.equals(c.name(), getString(KY_RTCAL))) return true;

    return false;
  }

  public boolean toSmoothRT() { return getBoolean(KY_SMOOTHRT); }
  public boolean useiRT() { return getBoolean(KY_USE_IRT); }

  public LcSettings setSpan(              float s) { set(KY_RT_SPAN,   s); return this; }
  public LcSettings setPeakWidth(         float s) { set(KY_LC_WIDTH,  s); return this; }
  public LcSettings setMinPeakExclusivity(float s) { set(KY_PEAK_EXCL, s); return this; }
  public LcSettings setMinApex(           float s) { set(KY_PEAK_APEX, s); return this; }
  public LcSettings setBaseRI(            float s) { set(KY_PEAK_BASE, s); return this; }
  public LcSettings setOuterMultiple(     float s) { set(KY_PEAK_MULTIPLE, s); return this; }
  public LcSettings setBandwidth(         float s) { set(KY_BANDWIDTH, s); return this; }
  public LcSettings setQuanSpan(          float s) { set(KY_QUAN_SPAN, s); return this; }
  public LcSettings setQuanOffset(        float s) { set(KY_QUAN_OFFSET, s); return this; }
  public LcSettings setMinSNR(            float s) { set(KY_SNR, s); return this; }

  public LcSettings setApexPts(             int s) { set(KY_APEX_PTS, s); return this; }
  public LcSettings setGridSize(            int s) { set(KY_GRID_SIZE, s); return this; }
  public LcSettings setCalMethod(   calibration s) { set(KY_RTCAL, s.name()); return this; }

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
