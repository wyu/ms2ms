package org.ms2ms.data.ms;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.data.Dataset;
import org.ms2ms.data.HData;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.io.MsReaders;
import org.ms2ms.nosql.ms.HBaseProteomics;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
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
abstract public class LcMsMsDataset implements Dataset
{
  public static String sTempDir = "/tmp/";
  public static String COL_MZ = "mz";
  public static String COL_RT = "rt";
  public static String COL_OFFSET = "offset";
  public static String COL_TIC = "TIC";

  protected String             mName, mSpCacheName, mRawfileRoot, mResultRoot;
  private Collection<String> mRawFilenames;
  protected Dataframe          mSummary, mSurvey;
  private HData              mData; // contains the sample, variable and intensities data for downstream analysis
  private ClusteringSettings mSpSettings;

//  private MultiTreeTable<Double, Double, Long> mMzRtFileOffset; // mz and RT index to the file offset for each MS/MS
//  private TreeMultimap<  Double, Long>         mTicFileOffset;  // negative TIC to put the most intense in the front
//  private TreeMultimap<  Double, Double>       mTicMz;          // negative TIC to put the most intense in the front

  /** the constructors **/
  public LcMsMsDataset()         { super(); }
  public LcMsMsDataset(String s) { mName=s; }

  /** Simple getter and setters  **/
  public String                               getName()               { return mName; }
  public HData                                getData()               { return mData; }
  public MultiTreeTable<Double, Double, Long> getMzRtFileOffset()     { return null; }
  public TreeMultimap<  Double, Long>         getTicFileOffset( )     { return null; }
  public ClusteringSettings                   getClusteringSettings() { return mSpSettings; }
  public Collection<String>                   getRawFilenames()       { return mRawFilenames; }
  public String                               getRawfileRoot()        { return mRawfileRoot; }
  public String                               getResultRoot()         { return mResultRoot; }
  public Dataframe                            getSummary()            { return mSummary; }

  public LcMsMsDataset setName(String    s) { mName = s; return this; }
  public LcMsMsDataset setRawFilename(String... s) { mRawFilenames = Lists.newArrayList(s); return this; }
  public LcMsMsDataset setRawFileRoot(String s)    { mRawfileRoot  = s; return this; }
  public LcMsMsDataset setResultsRoot(String s)    { mResultRoot  = s; return this; }
  public LcMsMsDataset setSummary(Dataframe s) { mSummary = s; return this; }

  public static int seekRow(Dataframe data, String run, String scan)
  {
    return 0;
  }
  abstract public Dataframe readMsMsWithAnnotations();
  abstract public void init();

  /** more complex access functions **/
  public LcMsMsDataset add(MsnSpectrum ms, long offset)
  {
    if (ms!=null)
    {
      if (mSurvey==null) mSurvey = new Dataframe("MS/MS Survey");
      int row = mSurvey.rows().size()+1;
      mSurvey.put(row, MaxQuant.V_MZ, ms.getPrecursor().getMz()
            ).put(row, MaxQuant.V_RT, ms.getRetentionTimes().getFirst().getTime()
            ).put(row, MaxQuant.V_OFFSET, offset
            ).put(row, MaxQuant.V_TIC, ms.getTotalIonCurrent());
    }
    return this;
  }
  public RandomAccessFile getSpCacheFile(int mslevel) throws IOException
  {
    if (!Strs.isSet(mSpCacheName))
      mSpCacheName = (Strs.isSet(getName())?getName():System.nanoTime()+"");

    return new RandomAccessFile(sTempDir + mSpCacheName+".ms"+mslevel, "rw");
  }
  public Collection<Long> getFileoffsetsByMzRt(double mz, double rt)
  {
//    return getMzRtFileOffset()!=null ? getMzRtFileOffset().subset(
//        getClusteringSettings().getInstrument().getPrecursorTol().getMin(mz),
//        getClusteringSettings().getInstrument().getPrecursorTol().getMax(mz),
//        getClusteringSettings().getRtWindow().getMin(rt),
//        getClusteringSettings().getRtWindow().getMax(rt)) : null;

    return null;
  }

  /*** Methods that apply directly to LcMsDataSet   *****
   *
   */

  public void importMaxquant(String root, String mzXML)
  {

  }
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
