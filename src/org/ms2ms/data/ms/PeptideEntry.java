package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.algo.MsStats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.pride.jmztab.model.Mod;

import java.util.List;

/** An entry of peptide without the actual sequence or mass information
 *
 * Created by yuw on 8/3/16.
 */
public class PeptideEntry
{
  private IonType mIonType; // enum
  private int  mTerminal;
  private long mProteinKey;
  private List<ModLocation> mMods;

  public PeptideEntry() { super(); }
  public PeptideEntry(long pointer, int terminal, IonType type, List<ModLocation> loc)
  {
    super();
    mProteinKey=pointer; mTerminal=terminal; mIonType=type; mMods=loc;
  }

  public long    getProteinKey() { return mProteinKey; }
  public int     getTerminal()   { return mTerminal; }
  public List<ModLocation> getMods() { return mMods; }

  public IonType getIonType()    { return mIonType; }
  public boolean isY()           { return IonType.y==mIonType; }
  public boolean isB()           { return IonType.b==mIonType; }

  @Override
  public String toString()
  {
    return getProteinKey()+":"+getTerminal()+"%"+getIonType()+(Tools.isSet(mMods) ?(", $"+ Strs.toString(mMods,";")):"");
  }
}
