package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.utils.Settings;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 10/4/14
 * Time: 7:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class ClusteringSettings extends Settings
{
  public static final String TOL_RT     = "RTTol";      // likely variation
  public static final String WINDOW_RT  = "RT Window";  // max variation
  public static final String SPEC_TYPE  = "SpecType";

  private MsSettings mInstrument;

  public ClusteringSettings() { super(); }
  public ClusteringSettings(MsSettings ms) { super(); mInstrument=ms; }

  public MsSettings getInstrument() { return mInstrument; }
  public Tolerance  getRtTol()    { return properties!=null?(Tolerance )properties.get(TOL_RT):null; }
  public Tolerance  getRtWindow() { return properties!=null?(Tolerance )properties.get(WINDOW_RT):null; }

  public byte[]    getSpecType()     { return getBytes(SPEC_TYPE); }
//  public int       getZFloat()       { return getInteger(Z_FLOAT); }

  public ClusteringSettings setRt(float tol, float window)
  {
    properties.put(   TOL_RT, new AbsoluteTolerance(tol));
    properties.put(WINDOW_RT, new AbsoluteTolerance(window));

    return this;
  }
}
