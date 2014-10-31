package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.alg.Peaks;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.utils.Settings;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/24/14
 * Time: 11:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class MsSettings extends Settings
{
  public static final String TOL_PREC  = "PrecTol";
  public static final String TOL_FRAG  = "FragTol";
  public static final String FRAG_MODE = "FragMode";

  public MsSettings() { super(); properties.put(ClusteringSettings.Z_FLOAT, 0); }
  public MsSettings(byte[] type, Tolerance precursor, Tolerance frag, String mode)
  {
    super();
    properties.put(TOL_PREC,  precursor);
    properties.put(TOL_FRAG,  frag);
    properties.put(FRAG_MODE, mode);
  }

  public Tolerance getPrecursorTol() { return properties!=null?(Tolerance )properties.get(TOL_PREC):null; }
  public Tolerance getFragmentTol()  { return properties!=null?(Tolerance )properties.get(TOL_FRAG):null; }
  public String    getFragMode()     { return getString(FRAG_MODE); }
}
