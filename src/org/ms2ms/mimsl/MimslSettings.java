package org.ms2ms.mimsl;

import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.ms.ClusteringSettings;
import org.ms2ms.data.ms.MsSettings;
import org.ms2ms.algo.Peaks;
import org.ms2ms.nosql.ms.HBasePeakList;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/20/14
 * Time: 10:42 PM
 * To change this template use File | Settings | File Templates.
 */
public class MimslSettings extends MsSettings
{
  public static final String SPEC_TYPE = "SpecType";
  public static final String Z_FLOAT   = "Zfloat";

  public static MimslSettings ORBI_HL_CID  =null;
  public static MimslSettings ORBI_HL_ETD  =null;
  public static MimslSettings ORBI_HH_CID  =null;
  public static MimslSettings ORBI_HH_HCD  =null;
  public static MimslSettings QEXACT_HH_HCD=null;
  public static MimslSettings QTOF_MM_CID  =null;
  public static MimslSettings LTQ_LL_CID   =null;

  //HBaseProteomics.prepareLib("/media/data/splib/2013", lib, 50d, 450d, 7, 4d) + " entries prepared");
  private double half_width=50d, min_mz=450d, min_snr=4d;
  private int min_pk=4;

  static
  {
    ORBI_HL_CID   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.5d), Peaks.CID);
    ORBI_HL_ETD   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new AbsoluteTolerance(0.5d), Peaks.ETD);
    ORBI_HH_CID   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new PpmTolerance(15d), Peaks.CID);
    ORBI_HH_HCD   = new MimslSettings(HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(15d), new PpmTolerance(15d), Peaks.HCD);
    QEXACT_HH_HCD = new MimslSettings(HBasePeakList.SPEC_TRAP_HCD, new PpmTolerance(15d), new PpmTolerance(15d), Peaks.HCD);
    QTOF_MM_CID   = new MimslSettings(HBasePeakList.SPEC_QTOF,     new PpmTolerance(250d),new PpmTolerance(250d), Peaks.CID);
  }
  public MimslSettings() { super(); }
  public MimslSettings(byte[] type, Tolerance precursor, Tolerance frag, String mode)
  {
    super(type, precursor, frag, mode);
    properties.put(ClusteringSettings.SPEC_TYPE, type);
    properties.put(MsSettings.Z_FLOAT,   0);
  }
  public static MimslSettings valueOf(String s)
  {
    if      (s.equals(Peaks.OBT_HR_CID)) return ORBI_HH_CID;
    else if (s.equals(Peaks.LTQ_CID))    return LTQ_LL_CID;
    else if (s.equals(Peaks.OBT_HCD))    return ORBI_HH_HCD;
    else if (s.equals(Peaks.QTOF))       return QTOF_MM_CID;

    return ORBI_HL_CID;
  }
  public MimslSettings set(double halfw, double minmz, double minsnr, int minpk)
  {
    half_width=halfw; min_mz=minmz; min_snr=minsnr; min_pk=minpk;
    return this;
  }
  public double getHalfWidth() { return half_width; }
  public double getMinMz()     { return min_mz; }
  public double getMinSNR()    { return min_snr; }
  public int    getMinPeaks()  { return min_pk; }
  public byte[] getSpecType()  { return getBytes(ClusteringSettings.SPEC_TYPE); }
  public int    getZFloat()    { return getInteger(MsSettings.Z_FLOAT); }
}
