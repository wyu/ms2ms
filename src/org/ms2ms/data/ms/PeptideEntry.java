package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;

/** An entry of peptide without the actual sequence or mass information
 *
 * Created by yuw on 8/3/16.
 */
public class PeptideEntry
{
  private IonType mIonType; // enum
  private int  mTerminal;
  private long mProteinKey;

  public PeptideEntry() { super(); }
  public PeptideEntry(long pointer, int terminal, IonType type)
  {
    super();
    mProteinKey=pointer; mTerminal=terminal; mIonType=type;
  }

  public long    getProteinKey() { return mProteinKey; }
  public int     getTerminal()   { return mTerminal; }

  public IonType getIonType()    { return mIonType; }
  public boolean isY()           { return IonType.y==mIonType; }
  public boolean isB()           { return IonType.b==mIonType; }
}
