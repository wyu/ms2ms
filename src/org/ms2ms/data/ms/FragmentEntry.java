package org.ms2ms.data.ms;

/**
 * Created by yuw on 8/4/16.
 */
public class FragmentEntry implements Comparable<FragmentEntry>
{
  private Long mProteinKey, mPeptideKey;
  private Double mNextFragment;

  public FragmentEntry() { super(); }
  public FragmentEntry(Long protein, Long peptide, Double frag)
  {
    super();
    mProteinKey=protein; mPeptideKey=peptide; mNextFragment=frag;
  }

  public Long   getProteinKey() { return mProteinKey; }
  public Long   getPeptideKey() { return mPeptideKey; }
  public Double next()          { return mNextFragment; }

  @Override
  public int compareTo(FragmentEntry o)
  {
    return Long.compare(mProteinKey, o.getProteinKey());
  }
}
