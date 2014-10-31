package org.ms2ms.data.ms;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.data.Dataset;
import org.ms2ms.data.HData;
import org.ms2ms.data.MultiTreeTable;
import org.ms2ms.utils.Tools;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Collection;

/** The logical unit of LC-MS/MS runs from a given study. Implicitly, they should have the same experimental
 *  characteristics and similar protein composition that are amendable for clustering
 *
 * User: wyu
 * Date: 10/3/14
 * To change this template use File | Settings | File Templates.
 */
public class LcMsMsDataset implements Dataset
{
  public static String sTempDir = "/tmp/";
  private String             mName, mSpCacheName;
  private Collection<String> mRawFilenames;
  private HData              mData; // contains the sample, variable and intensities data for downstream analysis
  private ClusteringSettings mSpSettings;

  private MultiTreeTable<Double, Double, Long> mMzRtFileOffset; // mz and RT index to the file offset for each MS/MS
  private TreeMultimap<  Double, Long>         mTicFileOffset;  // negative TIC to put the most intense in the front
  private TreeMultimap<  Double, Double>       mTicMz;          // negative TIC to put the most intense in the front

  /** the constructors **/
  public LcMsMsDataset()         { super(); }
  public LcMsMsDataset(String s) { mName=s; }

  /** Simple getter and setters  **/
  public String                               getName()               { return mName; }
  public HData                                getData()               { return mData; }
  public MultiTreeTable<Double, Double, Long> getMzRtFileOffset()     { return mMzRtFileOffset; }
  public TreeMultimap<  Double, Long>         getTicFileOffset( )     { return mTicFileOffset; }
  public ClusteringSettings                   getClusteringSettings() { return mSpSettings; }
  public Collection<String>                   getRawFilenames()       { return mRawFilenames; }

  public LcMsMsDataset setName(String    s) { mName = s; return this; }
  public LcMsMsDataset setRawFilename(String... s) { mRawFilenames = Lists.newArrayList(s); return this; }

  /** more complex access functions **/
  public LcMsMsDataset add(MsnSpectrum ms, long offset)
  {
    if (ms!=null)
    {
      if (mMzRtFileOffset==null) mMzRtFileOffset=MultiTreeTable.create();
      if ( mTicFileOffset==null)  mTicFileOffset=TreeMultimap.create();
      if (         mTicMz==null)          mTicMz=TreeMultimap.create();

      mMzRtFileOffset.put(ms.getPrecursor().getMz(), ms.getRetentionTimes().getFirst().getTime(), offset);
       mTicFileOffset.put(ms.getTotalIonCurrent()*-1d, offset);
               mTicMz.put(ms.getTotalIonCurrent()*-1d, ms.getPrecursor().getMz());
    }
    return this;
  }
  public RandomAccessFile getSpCacheFile(int mslevel) throws IOException
  {
    if (!Tools.isSet(mSpCacheName))
      mSpCacheName = (Tools.isSet(getName())?getName():System.nanoTime()+"");

    return new RandomAccessFile(sTempDir + mSpCacheName+".ms"+mslevel, "rw");
  }
  public Collection<Long> getFileoffsetsByMzRt(double mz, double rt)
  {
    return getMzRtFileOffset()!=null ? getMzRtFileOffset().subset(
        getClusteringSettings().getPrecursorTol().getMin(mz),
        getClusteringSettings().getPrecursorTol().getMax(mz),
        getClusteringSettings().getRtWindow().getMin(rt),
        getClusteringSettings().getRtWindow().getMax(rt)) : null;

  }

  /*** Methods that apply directly to LcMsDataSet   *****
   *
   */

  /** Determine the appropriate parameters in situ.
   *
   * @param settings
   * @return
   */
  public ClusteringSettings learn(ClusteringSettings settings)
  {

    return settings;
  }
  /** perform self-clustering among the MS/MS in the dataset.
   *
   * @return
   */
  public LcMsMsDataset cluster()
  {


    return this;
  }
}
