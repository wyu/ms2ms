package org.ms2ms.data.ms;

import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 8/4/16.
 */
public class FragmentEntry implements Comparable<FragmentEntry>
{
  private int mLength;
  private Long mPeptideKey;
  private Double mMH;
  private FragmentEntry mPrev;

  public FragmentEntry() { super(); }
  public FragmentEntry(Double mh, Long peptide, FragmentEntry prev, int len)
  {
    super();
    mMH=mh; mPeptideKey=peptide; mPrev=prev; mLength=len;
  }

  public int           getLength()     { return mLength; }
  public Double        getMH()         { return mMH; }
  public Long          getPeptideKey() { return mPeptideKey; }
  public FragmentEntry getPrev()       { return mPrev; }

  public FragmentEntry setLen(int s)   { mLength=s; return this; }

  @Override
  public int compareTo(FragmentEntry o)
  {
    return Long.compare(mPeptideKey, o.getPeptideKey());
  }

  @Override
  public String toString()
  {
    return mPeptideKey+":"+mLength+"#"+ (mMH!=null?Tools.d2s(mMH, 4):"NUL")+"<-"+(mPrev!=null?Tools.d2s(mPrev.getMH(),4):"NUL");
  }
}
