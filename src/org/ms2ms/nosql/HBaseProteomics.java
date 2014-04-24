package org.ms2ms.nosql;

import org.apache.hadoop.hbase.client.HConnection;
import org.apache.hadoop.hbase.client.HConnectionManager;
import org.apache.hadoop.hbase.client.HTableInterface;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;

import java.io.IOException;
import java.util.Collection;

/**
 * Created by wyu on 4/20/14.
 */
public class HBaseProteomics extends HBaseAbstract
{
  public static void ensurePeakListTable() throws IOException
  {
      createTable(HBasePeakList.TBL_PEAKLIST, new String[] {
      HBasePeakList.FAM_FLAG.toString(), HBasePeakList.FAM_PRECURSOR.toString(),
      HBasePeakList.FAM_STAT.toString()});
  }
  public static void save(Collection<LibrarySpectrum> spectra) throws IOException
  {
    // ensure that the table has been created
    ensurePeakListTable();
    // connection to the cluster
    HConnection conn = HConnectionManager.createConnection(conf);

    // When the cluster connection is established get an HTableInterface for each operation or thread.
    // HConnection.getTable(...) is lightweight. The table is really just a convenient place to call
    // table method and for a temporary batch cache.
    // It is in fact less overhead than HTablePool had when retrieving a cached HTable.
    // The HTableInterface returned is not thread safe as before.
    // It's fine to get 1000's of these.
    // Don't cache the longer than the lifetime of the HConnection
    HTableInterface table = conn.getTable(HBasePeakList.TBL_PEAKLIST);

    // save the spectra to the table
    for (LibrarySpectrum spec : spectra)
      HBasePeakList.save(table, spec);

    // just flushes outstanding commit, no futher cleanup needed, can be omitted.
    // HConnection holds no references to the returned HTable objects, they can be GC'd as soon as they leave scope.
    table.close(); conn.close(); // done with the cluster, release resources
  }

}
