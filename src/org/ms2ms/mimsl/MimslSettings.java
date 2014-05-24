package org.ms2ms.mimsl;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.alg.Peaks;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.utils.Settings;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/20/14
 * Time: 10:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class MimslSettings extends Settings
{
  public static final String SPEC_TYPE = "SpecType";
  public static final String TOL_PREC  = "PrecTol";
  public static final String TOL_FRAG  = "FragTol";
  public static final String FRAG_MODE = "FragMode";
  public static final String Z_FLOAT   = "Zfloat";

  public static MimslSettings ORBI_HL_CID  =null;
  public static MimslSettings ORBI_HL_ETD  =null;
  public static MimslSettings ORBI_HH_CID  =null;
  public static MimslSettings ORBI_HH_HCD  =null;
  public static MimslSettings QEXACT_HH_HCD=null;
  public static MimslSettings QTOF_MM_CID  =null;
  public static MimslSettings LTQ_LL_CID   =null;

  static
  {
    ORBI_HL_CID   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.5d), Peaks.CID);
    ORBI_HL_ETD   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.5d), Peaks.ETD);
    ORBI_HH_CID   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new PpmTolerance(15d), Peaks.CID);
    ORBI_HH_HCD   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new PpmTolerance(15d), Peaks.HCD);
    QEXACT_HH_HCD = new MimslSettings(HBasePeakList.SPEC_TRAP_HCD, new PpmTolerance(15d), new PpmTolerance(15d), Peaks.HCD);
    QTOF_MM_CID   = new MimslSettings(HBasePeakList.SPEC_QTOF,     new PpmTolerance(250d),new PpmTolerance(250d), Peaks.CID);
  }
  public MimslSettings() { super(); properties.put(Z_FLOAT, 0); }
  public MimslSettings(byte[] type, Tolerance precursor, Tolerance frag, String mode)
  {
    super();
    properties.put(SPEC_TYPE, type);
    properties.put(TOL_PREC,  precursor);
    properties.put(TOL_FRAG,  frag);
    properties.put(FRAG_MODE, mode);
    properties.put(Z_FLOAT,   0);
  }
/*
  public MimslSettings(char type, Tolerance precursor, Tolerance frag)
  {
    super();
    properties.put(SPEC_TYPE, type);
    properties.put(TOL_PREC,  precursor);
    properties.put(TOL_FRAG,  frag);
  }
*/

  public byte[]    getSpecType()     { return getBytes(SPEC_TYPE); }
  public Tolerance getPrecursorTol() { return properties!=null?(Tolerance )properties.get(TOL_PREC):null; }
  public Tolerance getFragmentTol()  { return properties!=null?(Tolerance )properties.get(TOL_FRAG):null; }
  public String    getFragMode()     { return getString(FRAG_MODE); }
  public int       getZFloat()       { return getInteger(Z_FLOAT); }

  public static MimslSettings valueOf(String s)
  {
    switch (s)
    {
      case Peaks.OBT_HR_CID: return ORBI_HH_CID;
      case Peaks.LTQ_CID:    return LTQ_LL_CID;
      case Peaks.OBT_HCD:    return ORBI_HH_HCD;
      case Peaks.QTOF:       return QTOF_MM_CID;
    }
    return ORBI_HL_CID;
  }
}
