package org.ms2ms.data.ms;

import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.data.Dataset;
import org.ms2ms.data.MultiTreeTable;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 10/3/14
 * Time: 10:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class LcMsMsDataset implements Dataset
{
  private String mName;
  private MultiTreeTable<Double, Double, Long> mMzRtFileOffset;
  private TreeMultimap<  Double, Long>         mTicFileOffset; // negative TIC to put the most intense in the front

  public LcMsMsDataset()         { super(); }
  public LcMsMsDataset(String s) { mName=s; }

  public String getName() { return mName; }
  public MultiTreeTable<Double, Double, Long> getMzRtFileOffset() { return mMzRtFileOffset; }
  public TreeMultimap<  Double, Long>         getTicFileOffset( ) { return mTicFileOffset; }

  public LcMsMsDataset add(MsnSpectrum ms, long offset)
  {
    if (ms!=null)
    {
      if (mMzRtFileOffset==null) mMzRtFileOffset=MultiTreeTable.create();
      if ( mTicFileOffset==null)  mTicFileOffset=TreeMultimap.create();

      mMzRtFileOffset.put(ms.getPrecursor().getMz(), ms.getRetentionTimes().getFirst().getTime(), offset);
       mTicFileOffset.put(ms.getTotalIonCurrent()*-1d, offset);
    }
    return this;
  }
}
