package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

/** An entry of peptide without the actual sequence or mass information
 *
 * Created by yuw on 8/3/16.
 */
public class PeptideEntry
{
  private IonType mIonType; // enum
  private int  mTerminal;
//  private Long mProteinKey;
  private int mPeptSeqKey;
//  private List<ModLocation> mMods;
  private ModLocation[] mMods;

  public PeptideEntry() { super(); }
  public PeptideEntry(int pointer, int terminal, IonType type, List<ModLocation> loc)
  {
    super();
    setPeptSeqKey(pointer); mTerminal=terminal; mIonType=type; mMods=Tools.isSet(loc)?loc.toArray(new ModLocation[] {}):null;
  }

  public PeptideEntry setPeptSeqKey(int s)
  {
    mPeptSeqKey=s;
    return this;
  }
  public int    getPeptSeqKey() { return mPeptSeqKey; }
  public int     getTerminal()   { return mTerminal; }
  public List<ModLocation> getMods() { return mMods!=null?Arrays.asList(mMods):null; }
  public Map<Integer, Double> getLocMods()
  {
    if (Tools.isSet(mMods))
    {
      Map<Integer, Double> mods = new TreeMap<>();
      for (ModLocation mod : mMods)
        mods.put(mod.locations, mod.mods);

      return mods;
    }

    return null;
  }

  public IonType getIonType()    { return mIonType; }
  public boolean isY()           { return IonType.y==mIonType; }
  public boolean isB()           { return IonType.b==mIonType; }

  @Override
  public String toString()
  {
    return getPeptSeqKey()+":"+getTerminal()+"%"+getIonType()+(Tools.isSet(mMods) ?(", $"+ Strs.toString(mMods,";")):"");
  }
}
