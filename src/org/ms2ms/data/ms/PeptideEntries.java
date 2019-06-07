package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/** An entry of key pointing to multiple proteins, without the actual sequence or mass information
 *
 * Created by yuw on 8/3/16.
 */
public class PeptideEntries
{
  private IonType mIonType; // enum
  private int  mTerminal;
  private Set<Long> mProteinKeys;
  private List<ModLocation> mMods;

  public PeptideEntries() { super(); }
  public PeptideEntries(long pointer, int terminal, IonType type, List<ModLocation> loc)
  {
    super();
    mTerminal=terminal; mIonType=type; mMods=loc;
    mProteinKeys = new HashSet<>(); mProteinKeys.add(pointer);
  }

  public Set<Long>         getProteinKeys() { return mProteinKeys; }
  public int               getTerminal()   { return mTerminal; }
  public List<ModLocation> getMods() { return mMods; }

  public IonType getIonType()    { return mIonType; }
  public boolean isY()           { return IonType.y==mIonType; }
  public boolean isB()           { return IonType.b==mIonType; }

  @Override
  public String toString()
  {
    return getProteinKeys().size()+":"+getTerminal()+"%"+getIonType()+(Tools.isSet(mMods) ?(", $"+ Strs.toString(mMods,";")):"");
  }
}
