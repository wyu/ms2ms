package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.algo.Peaks;
import org.ms2ms.nosql.ms.HBasePeakList;
import org.ms2ms.utils.Settings;

/** Parameters associated with a single stage of MS analyzer
 *
 * User: wyu
 * Date: 5/24/14
 */
public class MsSettings extends Settings
{
  public static final String TOL_PREC  = "PrecTol";
  public static final String TOL_FRAG  = "FragTol";
  public static final String ISOLATION_PREC  = "Precursor isolation";
  public static final String FRAG_MODE = "FragMode";
  public static final String Z_FLOAT    = "Zfloat";

  public static final MsSettings LTQ          = new MsSettings(HBasePeakList.SPEC_TRAP_CID, new AbsoluteTolerance(0.5d), new AbsoluteTolerance(0.5d),  Peaks.CID);
  public static final MsSettings ORBITRAP     = new MsSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.5d),  Peaks.CID);
  public static final MsSettings ORBITRAP_HR  = new MsSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.05d), Peaks.CID);
  public static final MsSettings ORBITRAP_HCD = new MsSettings(HBasePeakList.SPEC_TRAP_HCD, new PpmTolerance(15d), new AbsoluteTolerance(0.05d), Peaks.HCD);

  public static MsSettings    ORBI = null;
  static
  {
    ORBI = new MsSettings(
    ).setFragmentTol( new OffsetPpmTolerance(15d, 0d)
    ).setPrecursorTol(new OffsetPpmTolerance(15d, 0d));
  }

  public MsSettings() { super(); properties.put(Z_FLOAT, 0); }
  public MsSettings(MsSettings s) { super(s);  }
  public MsSettings(byte[] type, Tolerance precursor, Tolerance frag, String mode)
  {
    super();
    properties.put(TOL_PREC,  precursor);
    properties.put(TOL_FRAG,  frag);
    properties.put(FRAG_MODE, mode);
  }

  public Tolerance getPrecursorTol() { return properties!=null?(Tolerance )properties.get(TOL_PREC):null; }
  public Tolerance getPrecursorIsolation()  { return properties!=null?(Tolerance )properties.get(ISOLATION_PREC):null; }
  public Tolerance getFragmentTol()  { return properties!=null?(Tolerance )properties.get(TOL_FRAG):null; }
  public String    getFragMode()     { return getString(FRAG_MODE); }

  public MsSettings setFragmentTol( Tolerance s) { properties.put(TOL_FRAG, s); return this; }
  public MsSettings setPrecursorTol(Tolerance s) { properties.put(TOL_PREC, s); return this; }

  @Override
  public MsSettings clone() throws CloneNotSupportedException
  {
    return (MsSettings )super.clone();
  }
}
