package org.ms2ms.nosql;

import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.HColumnDescriptor;
import org.apache.hadoop.hbase.HTableDescriptor;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.ms2ms.alg.Peaks;
import org.ms2ms.mimsl.MIMSL;
import org.ms2ms.mzjava.AnnotatedPeak;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.UUID;

/**
 * Created by wyu on 4/20/14.
 */
public class HBaseProteomics extends HBaseAbstract
{
  public static void ensureTables() throws IOException
  {
    // make sures the tables were properly created
     createTable(HBasePeakList.TBL_PEAKLIST, new String[] {
      HBasePeakList.FAM_FLAG, HBasePeakList.FAM_PRECURSOR, HBasePeakList.FAM_PROP});
    createTable(HBasePeakList.TBL_MSMSINDEX, new String[] {
      HBasePeakList.FAM_ID, HBasePeakList.FAM_PROP});
  }
  public static void listTables() throws IOException
  {
    HConnection conn = HConnectionManager.createConnection(HBaseConfiguration.create());
    for (HTableDescriptor table : conn.listTables())
    {
      // get the number of row. Can be very expansive for a large table!!
      HTableInterface tbl = conn.getTable(table.getTableName());
      ResultScanner scanner = tbl.getScanner(new Scan());
      long counts=0;
      for (Result rs = scanner.next(); rs != null; rs = scanner.next()) { counts++; }
      tbl.close();

      System.out.println(table.getTableName() + " with " + table.getColumnFamilies().length + " Col families and " + counts + " Rows");
      for (HColumnDescriptor col : table.getFamilies())
      {
        System.out.println("  " + col.getNameAsString());
      }
    }
    conn.close();
  }
  public static void save(Collection<LibrarySpectrum> spectra) throws IOException
  {
    // ensure that the table has been created
    ensureTables();
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
  // MIMSL
  public static void index(Collection<LibrarySpectrum> spectra, double half_width, double min_mz, int min_pk, double min_snr) throws IOException
  {
    // ensure that the table has been created
    ensureTables();
    // connection to the cluster
    HConnection conn = HConnectionManager.createConnection(conf);

    HTableInterface peaklist = conn.getTable(HBasePeakList.TBL_PEAKLIST),
                      indice = conn.getTable(HBasePeakList.TBL_MSMSINDEX);

    // save the spectra to the table
    long counts=0;
    for (LibrarySpectrum spec : spectra)
    {
      if (++counts%100 ==0) System.out.print(".");
      if (++counts%5000==0) System.out.print("\n");
//      System.out.println(Peaks.print(null, spec, false));

      spec.setId(UUID.randomUUID());
      HBasePeakList.save(peaklist, spec);
      List<AnnotatedPeak> sigs = MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr);

      for (AnnotatedPeak sig : sigs)
      {
//        System.out.println("Sig: " + Peaks.print(null, sig));
        // TODO need to work out the composite key incoporating the n/c and mod flag
        Put row = new Put(Bytes.toBytes(sig.getMz()));
        // byte[] family, byte[] qualifier, byte[] value
        row.add(Bytes.toBytes(HBasePeakList.FAM_ID),   Bytes.toBytes(HBasePeakList.COL_UUID), Bytes.toBytes(spec.getId().toString()));
        row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_SIG),  Bytes.toBytes(sig.getMz()));
        row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_SNR),  Bytes.toBytes(sig.getSNR()));
        row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_Z),    Bytes.toBytes(sig.getCharge()));
        // TODO row.add(Bytes.toBytes(HBasePeakList.FAM_MZ), Bytes.toBytes(HBasePeakList.COL_MMOD), Bytes.toBytes(spec.getPrecursor().getCharge()));
        indice.put(row);
      }
    }
    peaklist.close(); indice.close(); conn.close(); // release resources

    System.out.println("\nPreparaing the MsMs Index.");
  }
}
