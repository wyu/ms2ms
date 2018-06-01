package org.ms2ms.test;

import org.junit.Test;
import org.ms2ms.algo.PeptideLC;
import org.ms2ms.data.Features;
import org.ms2ms.data.ms.LcMsMsFeatures;
import org.ms2ms.data.ms.OffsetPpmTolerance;

import java.util.HashMap;
import java.util.Map;

public class RtPredictTest extends TestAbstract
{
  @Test
  public void LsRPv1()
  {
    // first test of LsRT port
    Map<String, Float> peptides = example1D();

    PeptideLC.SVR(peptides, 1f);
  }
  private Map<String, Float> example1D()
  {
    LcMsMsFeatures plasma1d = new LcMsMsFeatures(
        "/Users/kfvf960/OneDrive - AZCollaboration/contrib/data/Plasma.Mann/txt.20171226",
        new OffsetPpmTolerance(15), "20150805_QEp1_LC7_PhGe_SA_100min_2.*");

    Map<String, Float> peptide_rt = new HashMap<>();
    for (Features F : plasma1d.getIons().values())
    {
      if (F.get(LcMsMsFeatures.COL_PEPTIDE)!=null)
        peptide_rt.put(F.getStr(LcMsMsFeatures.COL_SEQ), F.getFloat(LcMsMsFeatures.COL_RT));
    }

    return peptide_rt;
  }
}
