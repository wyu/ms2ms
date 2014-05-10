package org.ms2ms.nosql;

import com.google.common.collect.Range;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.HColumnDescriptor;
import org.apache.hadoop.hbase.HTableDescriptor;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.filter.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.DoublePeakList;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.peakprocessor.PeakProcessorChain;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.core.ms.spectrum.Peak;
import org.expasy.mzjava.proteomics.io.ms.spectrum.MsLibReader;
import org.expasy.mzjava.proteomics.io.ms.spectrum.msp.MspCommentParser;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.expasy.mzjava.utils.URIBuilder;
import org.ms2ms.mimsl.MIMSL;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.mzjava.MspAnnotationResolver2;
import org.ms2ms.utils.Tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

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
  public static void index(Collection<LibrarySpectrum> spectra, char spec_type, double half_width, double min_mz, int min_pk, double min_snr) throws IOException
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

      index(peaklist, indice, spec_type, spec, MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr));
    }
    peaklist.close(); indice.close(); conn.close(); // release resources
  }
  public static void index(HTableInterface peaklist, HTableInterface indice, char spec_type,
                           LibrarySpectrum spec, List<AnnotatedPeak> sigs) throws IOException
  {
    spec.setId(UUID.randomUUID());
    HBasePeakList.save(peaklist, spec);

    boolean verbose = spec.getPrecursor().getMz()>=500.7 && spec.getPrecursor().getMz()<=500.8;

    if (verbose) System.out.print(Tools.d2s(spec.getPrecursor().getMz(), 3) + ", +" + spec.getPrecursor().getCharge() + "(" + sigs.size() + "): ");
    for (AnnotatedPeak sig : sigs)
    {
      if (verbose) System.out.print(Tools.d2s(sig.getMz(), 2) + ",");

      // TODO need to work out the composite key incoporating the n/c and mod flag
      Put row = new Put(HBasePeakList.row4MsMsIndex(spec_type, (float )spec.getPrecursor().getMz(), (byte )spec.getPrecursor().getCharge()));
      // byte[] family, byte[] qualifier, byte[] value
      row.add(Bytes.toBytes(HBasePeakList.FAM_ID),   Bytes.toBytes(HBasePeakList.COL_UUID), Bytes.toBytes(spec.getId().toString()));
      row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_MZ),   Bytes.toBytes(sig.getMz()));
      row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_SNR),  Bytes.toBytes(sig.getSNR()));
      row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_Z),    Bytes.toBytes(sig.getCharge()));
      row.add(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_IONS), Bytes.toBytes(sigs.size()));
      // TODO row.add(Bytes.toBytes(HBasePeakList.FAM_MZ), Bytes.toBytes(HBasePeakList.COL_MMOD), Bytes.toBytes(spec.getPrecursor().getCharge()));
      indice.put(row);
    }
    if (verbose && Tools.isSet(sigs)) System.out.println();
  }

  /** Query the HBase for the candidates. No scoring in this call
   *
   * @param signatures
   */
  public static Collection<AnnotatedSpectrum> query(PeakList<PepLibPeakAnnotation> signatures, char spec_type, Tolerance tol, double offset) throws IOException
  {
    if (signatures == null || signatures.getPrecursor()==null) return null;

    // setup the HTable for query
    Configuration conf = HBaseConfiguration.create();
    HConnection conn = null;
    try
    {
      conn = HConnectionManager.createConnection(conf);
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
    // get the number of row. Can be very expansive for a large table!!
    HTableInterface tbl = conn.getTable(HBasePeakList.TBL_MSMSINDEX);

    Map<UUID, AnnotatedSpectrum> id_candidate = new HashMap<UUID, AnnotatedSpectrum>();
    Collection<Range<Peak>> slices = MIMSL.enumeratePrecursors(tol, 0, signatures.getPrecursor());
    for (Range<Peak> range : slices)
    {
      Scan scan = new Scan(HBasePeakList.row4MsMsIndex(spec_type, range.lowerEndpoint().getMz()+offset, range.lowerEndpoint().getCharge()),
                           HBasePeakList.row4MsMsIndex(spec_type, range.upperEndpoint().getMz()+offset, range.lowerEndpoint().getCharge()));
      if (signatures.size()>0)
      {
        List<Filter> filters = new ArrayList<Filter>(signatures.size());
        for (int k=0; k<signatures.size(); k++)
        {
          // set the range of row keys
          filters.add(new FilterList(FilterList.Operator.MUST_PASS_ALL,
            new SingleColumnValueFilter(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_MZ), CompareFilter.CompareOp.GREATER_OR_EQUAL, new BinaryComparator(Bytes.toBytes(tol.getMin(signatures.getMz(k))))),
            new SingleColumnValueFilter(Bytes.toBytes(HBasePeakList.FAM_PROP), Bytes.toBytes(HBasePeakList.COL_MZ), CompareFilter.CompareOp.LESS_OR_EQUAL, new BinaryComparator(Bytes.toBytes(tol.getMax(signatures.getMz(k)))))));
        }
        FilterList filterList = new FilterList(FilterList.Operator.MUST_PASS_ONE, filters);
        scan.setFilter(filterList);
      }
      // go thro the peaks
      ResultScanner scanner = tbl.getScanner(scan);
      for (Result rs = scanner.next(); rs != null; rs = scanner.next()) {
        UUID id = UUID.fromString(
          HBaseAbstract.getString(rs, HBasePeakList.FAM_ID, HBasePeakList.COL_UUID));
        if (id_candidate.get(id) == null) id_candidate.put(id, new AnnotatedSpectrum());
        id_candidate.get(id).add(
          HBaseAbstract.getDouble(rs, HBasePeakList.FAM_PROP, HBasePeakList.COL_MZ),
          HBaseAbstract.getDouble(rs, HBasePeakList.FAM_PROP, HBasePeakList.COL_SNR));
        id_candidate.get(id).setId(id);
        // record the number the indice
        id_candidate.get(id).setIonIndexed(
          HBaseAbstract.getInt(rs, HBasePeakList.FAM_PROP, HBasePeakList.COL_IONS));
        id_candidate.get(id).setIonQueried(signatures.size());
        id_candidate.get(id).setIonMatched(id_candidate.get(id).size());
        id_candidate.get(id).setPrecursor(signatures.getPrecursor());
        id_candidate.get(id).setMzQueried(signatures.getPrecursor().getMz());
      }
      System.out.println("Candidates with the precursor m/z range of " + range.toString());
      for (UUID id : id_candidate.keySet())
      {
        System.out.print(
          Tools.d2s(id_candidate.get(id).getPrecursor().getMz(), 3) + ", +" + id_candidate.get(
            id).getPrecursor().getCharge() + "(" + id_candidate.get(id).size() + "): ");
        for (int i = 0; i < id_candidate.get(id).size(); i++) {
          System.out.print(Tools.d2s(id_candidate.get(id).getMz(i), 2) + ",");
        }
        System.out.println();
      }
    }
    return id_candidate.values();
  }
  public static long prepareMsp(File src, char spec_type, double half_width, double min_mz, int min_pk, double min_snr) throws IOException
  {
    // connection to the cluster
    HConnection conn = HConnectionManager.createConnection(conf);

    HTableInterface peaklist = conn.getTable(HBasePeakList.TBL_PEAKLIST),
                      indice = conn.getTable(HBasePeakList.TBL_MSMSINDEX);

    System.out.println("Preparing " + src.getAbsolutePath());

    BufferedReader reader = new BufferedReader(new FileReader(src));
    MsLibReader msp = new MsLibReader(reader, URIBuilder.UNDEFINED_URI,
      PeakList.Precision.FLOAT,
      new MspCommentParser(), new MspAnnotationResolver2(),
      Pattern.compile(
        "^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+)\\s+\"([^\"]+)\"$"),
      new PeakProcessorChain<PepLibPeakAnnotation>());

    long counts=0; LibrarySpectrum spec = null;
    while (msp.hasNext())
    {
      if (++counts%100 ==0) System.out.print(".");
      if (++counts%5000==0) System.out.print("\n");
      spec = msp.next();
      index(peaklist, indice, spec_type, spec, MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr));
      spec = null;
    }
    peaklist.close(); indice.close(); conn.close(); // release resources

    return counts;
  }
  public static long prepareMsps(String root, char spec_type, double half_width, double min_mz, int min_pk, double min_snr, String... msps) throws IOException
  {
    long counts = 0;
    for (String msp : msps)
    {
      counts += prepareMsp(new File(root + "/" + msp), spec_type, half_width, min_mz, min_pk, min_snr);
    }
    return counts;
  }
  public static Collection<AnnotatedSpectrum> loadPeakLists(Collection<AnnotatedSpectrum> spectra) throws IOException
  {
    HConnection    conn = HConnectionManager.createConnection(conf);
    HTableInterface tbl = conn.getTable(HBasePeakList.TBL_PEAKLIST);
    for (AnnotatedSpectrum spec : spectra)
    {
      HBasePeakList.getPeakList(tbl, spec.getId()).toPeakList(spec);
    }
    tbl.close();

    return spectra;
  }
}
