package org.ms2ms.data.ms;

import com.google.common.collect.Range;

/** A representation of LcMs feature from MaxQuant or similar program.
 *
 */
public class LcMsFeature
{
  private Measurable    mMz, mRt;
  private String        mPeptide;
  private Range<Double> mRtBound;

  public LcMsFeature() {};
  public LcMsFeature(String peptide, double mz, double rt)
  {
    mPeptide=peptide; mMz = Measurable.newMz(mz); mRt=Measurable.newRt(rt);
  }


}
