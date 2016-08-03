package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;

/** An entry of peptide without the actual sequence or mass information
 *
 * Created by yuw on 8/3/16.
 */
public class PeptideEntry
{
  private IonType mIonType; // enum
  private int mSequencePointer, mTerminal;

  public PeptideEntry() { super(); }
  public PeptideEntry(int pointer, int terminal, IonType type)
  {
    super();
    mSequencePointer=pointer; mTerminal=terminal; mIonType=type;
  }

  public int     getSequencePointer() { return mSequencePointer; }
  public int     getTerminal()        { return mTerminal; }

  public IonType getIonType()         { return mIonType; }
  public boolean isY()                { return IonType.y==mIonType; }
  public boolean isB()                { return IonType.b==mIonType; }
}
