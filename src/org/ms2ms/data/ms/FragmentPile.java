package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.ms2ms.algo.Peaks;
import org.ms2ms.math.Stats;

import java.util.*;

/** AbstractPile of FragmentEntry from the screen process. We need to store and sort them  **/
/* it consists of 1-N trunks of long arrays of FragmentEntry  */

public class FragmentPile extends AbstractPile<FragmentMatch>
{
  // zero or more trunks of fragment entry vectors to hold the bulk data
  private PeakMatch[] track = new PeakMatch[255]; // no more than 255 residue long for now.
  // holding the query results
//  private FragmentMatch[] mFragSeries = new FragmentMatch[255];

  public FragmentPile()       { super(); }
  public FragmentPile(int s)  { super(s); init(); }

  public PeakMatch[]   getTrack() { return track; }

  @Override
  public int getKeyAt(int pile, int idx) { return get(pile,idx).getEntry().getPeptideKey(); }
  @Override
  public FragmentMatch[] newPile(int s) { return new FragmentMatch[s]; }

  public int toTrack(PeakList ms, Map<Integer, Double> ion_freq, OffsetPpmTolerance tol)
  {
    // once a block of FragmentEntry was located in the 'pile', we can turn them into a 'track'
    int track_size = 0;
    for (int i=0; i<mSeriesEnd; i++)
    {
      FragmentEntry E = at(i).getEntry();
      int         ion = at(i).getObsIndex();
      // starting from the lower mass end
      double mz0 = Peaks.MnH2MnH(ms.getMz(ion), (int )E.getCharge(), 1);

      track[track_size++] = new PeakMatch(
          Stats.ppm(mz0, E.getMH())+tol.getOffset(mz0),   // mz deviation in ppm
          mz0,                         // mz of the observed fragment
          E.getLength(),               // order of the fragment in the series
          Math.abs(ms.getIntensity(ion)),               // the intensity
          ion_freq.get(ion)
      ).setOriginalMz(mz0).setCalcMz(E.getMH()).setIndex(ion
      ).setIonType(E.isProDirected()? IonType.p:IonType.unknown);
    }

    return track_size;
  }
}
