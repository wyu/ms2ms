package org.ms2ms.data.ms;

import org.ms2ms.algo.Peptides;

public class PsmSettings extends MsSettings
{
  private float[] mAAs;

  public static final String MIN_MH = "min mh";
  public static final String MOD_NT = "nt mod mh";
  public static final String MOD_CT = "ct mod mh";
  public static final String MAXZ_FRAG = "max fragment charge";
  public static final String MIN_MH_FRAGZ = "min frag mh to consider the higher charge state";
  public static final String MIN_FRAG_NL = "min length of fragment to consider neutral loss";


  public static PsmSettings ORBI_HCD=null, ORBI_HCD_TMT=null;

    static
    {
      ORBI_HCD     = new PsmSettings(MsSettings.ORBITRAP_HCD).setAAs(ResidueBase.CIM.getAAs());
      ORBI_HCD_TMT = new PsmSettings(MsSettings.ORBITRAP_HCD).setAAs(ResidueBase.TMT10.getAAs());
    }

    public PsmSettings() { super(); }
    public PsmSettings(MsSettings s) { super(s); }

    public double getMinMH() { return getDouble(MIN_MH, 0d); }
    public double getNtMod() { return getDouble(MOD_NT, 0d); }
    public double getCtMod() { return getDouble(MOD_CT, 0d); }

    public float[] getAAs() { return mAAs; }
    public PsmSettings setAAs(float[] s) { mAAs=s; return this; }

    public double getNt0() { return (getAAs()['^'] +   Peptides.H)+getNtMod(); }
    public double getCt0() { return (getAAs()['$'] + 2*Peptides.H)+getCtMod(); }

    public int    getMinFrag4NL()      { return getDouble(MIN_FRAG_NL, 2d).intValue(); }
    public int    getMaxFragCharge()   { return getDouble(MAXZ_FRAG, 1d).intValue(); }
    public double getMinMHFragCharge() { return getDouble(MIN_MH_FRAGZ, 350d); }
}
