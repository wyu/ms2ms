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
  private PeakEntry[] track = new PeakEntry[255]; // no more than 255 residue long for now.
  private double[]   scores = new double[255];

  public FragmentPile()       { super(); }
  public FragmentPile(int s)  { super(s); init(); }

  public PeakEntry[]   getTrack() { return track; }

  @Override
  public FragmentPile init()
  {
    super.init();
    for (int i=0; i<255; i++) track[i] = new PeakEntry();
    return this;
  }
  @Override
  public int getKeyAt(int pile, int idx) { return mDataPiles.get(pile)[idx].getEntry().getPeptideKey(); }
  @Override
  public FragmentMatch[] newPile(int s) { return new FragmentMatch[s]; }

  // combine the toTrack and inspec4screen without actually creating the FpmSlot object
  public double[] screenTrack(PeakList ms, double[] ion_freq, OffsetPpmTolerance tol, int[] tAccumSamples)
  {
    int gap=0, contig_start=0, contig_last=0, best=0;
    double percentile=0, score=0, sumAI=0d;

    int track_size = 0;
    for (int i=0; i<mSeriesEnd; i++)
    {
      FragmentSlot E = at(i).getEntry();
      int        ion = at(i).getObsIndex();

      // starting from the lower mass end
      float mz0 = (float )Peaks.MnH2MnH(ms.getMz(ion), (int )E.getCharge(), 1);

      PeakEntry M = track[track_size++];

      // mz deviation in ppm, mz of the observed fragment, order of the fragment in the series
      M.setValues((float )(Stats.ppm(mz0, E.getMH())+tol.getOffset(mz0)), mz0, E.getLength());
      M.setSNR((float )Math.abs(ms.getIntensity(ion)));               // the intensity
      M.setFrequency((float )ion_freq[ion]);

      gap          =  M.getCharge()-contig_last;
      percentile   = (M.getIntensity()*0.01); // set a minimum
      // accumualte the gap score
      if (gap>0)
      {
        M.calcGapScore(gap, 1f, tAccumSamples);
        score+=M.getScore() - Math.log10(percentile);
      }
      // move the contig start if gap>1
      if (contig_start==0 || M.getCharge()-contig_last>1)  contig_start=M.getCharge();
      // always update to the preceeding ion
      contig_last=M.getCharge(); sumAI+=percentile;
      // update the longest contig if qualified
      if (contig_last-contig_start>best) best=contig_last-contig_start;
    }
    scores[0]=track_size;
    scores[1]=best;
    scores[2]=-10d*score;

    return scores;
  }
  public int toTrack(PeakList ms, double[] ion_freq, OffsetPpmTolerance tol)
  {
    // once a block of FragmentEntry was located in the 'pile', we can turn them into a 'track'
    int track_size = 0;
    for (int i=0; i<mSeriesEnd; i++)
    {
      FragmentSlot  E = at(i).getEntry();
      int         ion = at(i).getObsIndex();

      // starting from the lower mass end
      float mz0 = (float )Peaks.MnH2MnH(ms.getMz(ion), (int )E.getCharge(), 1);

      PeakEntry M = track[track_size++];

      // mz deviation in ppm, mz of the observed fragment, order of the fragment in the series
      M.setValues((float )(Stats.ppm(mz0, E.getMH())+tol.getOffset(mz0)), mz0, E.getLength());
      M.setSNR((float )Math.abs(ms.getIntensity(ion)));               // the intensity
      M.setFrequency((float )ion_freq[ion]);

//      M.setOriginalMz(mz0);
//      M.setCalcMz(E.getMH());
//      M.setIndex(ion);

//      track[track_size++] = new PeakEntry(
//          Stats.ppm(mz0, E.getMH())+tol.getOffset(mz0),   // mz deviation in ppm
//          mz0,                         // mz of the observed fragment
//          E.getLength(),               // order of the fragment in the series
//          Math.abs(ms.getIntensity(ion)),               // the intensity
//          ion_freq.get(ion)
//      ).setOriginalMz(mz0).setCalcMz(E.getMH()).setIndex(ion
//      ).setIonType(E.isProDirected()? IonType.p:IonType.unknown);
    }

    return track_size;
  }
}
