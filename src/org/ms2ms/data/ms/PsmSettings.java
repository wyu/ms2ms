package org.ms2ms.data.ms;

import org.ms2ms.algo.Peptides;

public class PsmSettings extends MsSettings
{
    private float[] mAAs;

    public static final String MIN_MH = "min mh";
    public static final String MOD_NT = "nt mod mh";
    public static final String MOD_CT = "ct mod mh";

    public static PsmSettings ORBI_HCD  =null;

    static
    {
        ORBI_HCD = new PsmSettings(MsSettings.ORBITRAP_HCD);
    }

    public PsmSettings() { super(); }
    public PsmSettings(MsSettings s) { super(s); }

    public double getMinMH() { return getDouble(MIN_MH, 0d); }
    public double getNtMod() { return getDouble(MOD_NT, 0d); }
    public double getCtMod() { return getDouble(MOD_CT, 0d); }

    public float[] getAAs() { return mAAs; }

    public double getNt0() { return (getAAs()['^'] +   Peptides.H)+getNtMod(); }
    public double getCt0() { return (getAAs()['$'] + 2*Peptides.H)+getCtMod(); }
}
