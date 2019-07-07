package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Tools;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class ProteinSegments
{
  private IonType[]    mIonTypes;
  private int          mSegmentCapacity = 10000;
  private int          mProteinKey;
  private FpmSlot[][] mIonTypeSegments = new FpmSlot[IonType.values().length][];  // maximum of all ion types
  private Integer[]    mSegmentEnds    = new Integer[ IonType.values().length];

  public ProteinSegments() { super(); }
  public ProteinSegments(IonType... ions)
  {
    super();
    if (Tools.isSet(ions))
    {
      mIonTypes = ions;
      for (IonType ion : ions)
      {
        mIonTypeSegments[ion.ordinal()] = new FpmSlot[mSegmentCapacity];
        mSegmentEnds[    ion.ordinal()] = 0;
      }
    }
  }

  public int getProteinKey() { return mProteinKey; }
  public FpmSlot   getIon(IonType s, int p) { return mIonTypeSegments[s.ordinal()][p]; }
  public FpmSlot[] getIons(IonType s) { return mIonTypeSegments[s.ordinal()]; }

  public int getSegmentEnd(IonType ion) { return mSegmentEnds[ion.ordinal()]; }
  public int size()
  {
    return Tools.isSet(mSegmentEnds)?Stats.sum(mSegmentEnds):0;
  }
  public ProteinSegments start()
  {
    if (Tools.isSet(mSegmentEnds))
      Arrays.fill(mSegmentEnds,0);

    return this;
  }
  public ProteinSegments ensureCapacity(int n, double multiple)
  {
    if (Tools.isSet(mIonTypeSegments))
      for (IonType ion : mIonTypes)
        ensureCapacity(ion, n, multiple);

    return this;
  }

  public ProteinSegments ensureCapacity(IonType ion, int n, double multiple)
  {
    if (getSegmentEnd(ion)+n>=mIonTypeSegments[ion.ordinal()].length)
    {
      // extends the array
      int new_capacity = (int )(mIonTypeSegments[ion.ordinal()].length*multiple);
      mIonTypeSegments[ion.ordinal()] = Arrays.copyOf(mIonTypeSegments[ion.ordinal()], new_capacity);
    }
    return this;
  }
  public ProteinSegments addSegment(FpmMatch s)
  {
    // add the entry by ion type
    int p = s.getIonType().ordinal(); mProteinKey = s.getProtein();
    mIonTypeSegments[p][mSegmentEnds[p]] = s.getEntry();
    mSegmentEnds[p]++;

    return this;
  }

  public void sortByPosition()
  {
    for (IonType ion : mIonTypes)
    {
      Arrays.sort(mIonTypeSegments[ion.ordinal()], 0, mSegmentEnds[ion.ordinal()], new FpmSlot.PositionComparator());
    }
  }
  public void sortByGapScore(IonType ion)
  {
    Arrays.sort(mIonTypeSegments[ion.ordinal()], 0, mSegmentEnds[ion.ordinal()], new FpmSlot.GapScoreDecendComparator());
  }
}
