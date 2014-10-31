package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.alg.Peaks;
import org.ms2ms.nosql.HBasePeakList;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 10/4/14
 * Time: 7:53 PM
 * To change this template use File | Settings | File Templates.
 */
public class ClusteringSettings extends MsSettings
{
  public static final String TOL_RT     = "RTTol";      // likely variation
  public static final String WINDOW_RT  = "RT Window";  // max variation
  public static final String SPEC_TYPE = "SpecType";
  public static final String Z_FLOAT   = "Zfloat";

  public static final MsSettings ORBITRAP = new MsSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.5d), Peaks.CID);

  static
  {

  }


  public ClusteringSettings() { super(); }

  public Tolerance getRtTol()    { return properties!=null?(Tolerance )properties.get(TOL_RT):null; }
  public Tolerance getRtWindow() { return properties!=null?(Tolerance )properties.get(WINDOW_RT):null; }

  public byte[]    getSpecType()     { return getBytes(SPEC_TYPE); }

  public int       getZFloat()       { return getInteger(Z_FLOAT); }
}
