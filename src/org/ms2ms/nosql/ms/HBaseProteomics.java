package org.ms2ms.nosql.ms;

import com.google.common.collect.Range;
import org.apache.hadoop.hbase.HColumnDescriptor;
import org.apache.hadoop.hbase.HTableDescriptor;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.filter.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.peaklist.PeakProcessorChain;
import org.expasy.mzjava.proteomics.io.ms.spectrum.MsLibReader;
import org.expasy.mzjava.proteomics.io.ms.spectrum.SptxtReader;
import org.expasy.mzjava.proteomics.io.ms.spectrum.msp.MspCommentParser;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodModificationResolver;
import org.expasy.mzjava.proteomics.ms.consensus.PeptideConsensusSpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.expasy.mzjava.utils.URIBuilder;
import org.ms2ms.mimsl.MIMSL;
import org.ms2ms.mimsl.MimslSettings;
import org.ms2ms.mzjava.*;
import org.ms2ms.nosql.HBases;
import org.ms2ms.utils.Tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.*;

/**
 * Created by wyu on 4/20/14.
 */
public class HBaseProteomics extends HBases
{
  static public String HUMAN  = "Homo Sapians";
  static public String MOUSE  = "Mouse";
  static public String ECOLI  = "EColi";
  static public String YEAST  = "Yeast";
  static public String MM     = "Msmegmatis";
  static public String WORM   = "Ce";
  static public String DM     = "Dm";
  static public String BOVINE = "Bovine";
  static public String CHICK  = "Chicken";
  static public String DRADI  = "Dradiodurans";
  static public String DROSO  = "Drosophila";
  static public String RAT    = "Rat";
  static public String MIXED  = "Mixed";
  static public String POMBE  = "S.Pombe";

  public static void ensureTables() throws IOException
  {
    // make sures the tables were properly created
    createTable(HBasePeakList.TBL_PEAKLIST, HBases.FAM_FLAG, HBasePeakList.FAM_PRECURSOR, HBases.FAM_PROP);
    createTable(HBases.TBL_MSMSINDEX, HBases.FAM_ID, HBases.FAM_PROP);
    createTable(HBases.TBL_TABLE,    HBases.FAM_ID, HBases.FAM_PROP);
    createTable(Bytes.toString(HBaseSpLib.TBL_SPLIB), HBases.FAM_ID, HBases.FAM_PROP);
  }
  public static void listTables() throws IOException
  {
    HConnection conn = getConnection();
    for (HTableDescriptor table : conn.listTables())
    {
      // cells the number of row. Can be very expansive for a large table!!
      //HTableInterface tbl = conn.getTable(table.getTableName());
      HTableInterface tbl = conn.getTable(table.getName());
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
  public static void save(Collection<PeptideConsensusSpectrum> spectra) throws IOException
  {
    // ensure that the table contains been created
    ensureTables();
    // connection to the cluster
    HConnection conn = getConnection();

    // When the cluster connection is established cells an HTableInterface for each operation or thread.
    // HConnection.getTable(...) is lightweight. The table is really just a convenient place to call
    // table method and for a temporary batch cache.
    // It is in fact less overhead than HTablePool had when retrieving a cached HTable.
    // The HTableInterface returned is not thread safe as before.
    // It's fine to cells 1000's of these.
    // Don't cache the longer than the lifetime of the HConnection
    HTableInterface table = conn.getTable(HBasePeakList.TBL_PEAKLIST);

    // save the spectra to the table
    for (PeptideConsensusSpectrum spec : spectra)
      HBasePeakList.save(table, spec);

    // just flushes outstanding commit, no futher cleanup needed, can be omitted.
    // HConnection holds no references to the returned HTable objects, they can be GC'd as soon as they leave scope.
    table.close(); conn.close(); // done with the cluster, release resources
  }
  // MIMSL
  public static void index(Collection<PeptideConsensusSpectrum> spectra, byte[] spec_type, double half_width, double min_mz, int min_pk, double min_snr) throws IOException
  {
    // ensure that the table contains been created
    ensureTables();
    // connection to the cluster
    HConnection conn = getConnection();

    HTableInterface peaklist = conn.getTable(HBasePeakList.TBL_PEAKLIST),
                      indice = conn.getTable(HBases.TBL_MSMSINDEX);

    // save the spectra to the table
    long counts=0;
    for (PeptideConsensusSpectrum spec : spectra)
    {
      if (++counts%100 ==0) System.out.print(".");
      if (++counts%5000==0) System.out.print("\n");

      index(peaklist, indice, spec_type, spec, MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr));
    }
    peaklist.close(); indice.close(); conn.close(); // release resources
  }
  public static void index(HTableInterface peaklist, HTableInterface indice, byte[] spec_type,
                           PeptideConsensusSpectrum spec, List<AnnotatedPeak> sigs) throws IOException
  {
    spec.setId(UUID.randomUUID());
    HBasePeakList.save(peaklist, spec);

    //boolean verbose = spec.getPrecursor().getMz()>=500.7 && spec.getPrecursor().getMz()<=500.8;

    //if (verbose) System.out.print(Tools.d2s(spec.getPrecursor().getMz(), 3) + ", +" + spec.getPrecursor().getCharge() + "(" + sigs.size() + "): ");
    for (AnnotatedPeak sig : sigs)
    {
      //if (verbose) System.out.print(Tools.d2s(sig.getMz(), 2) + ",");

      // TODO need to work out the composite key incoporating the n/c and mod flag
      Put row = new Put(HBasePeakList.row4MsMsIndex(spec_type, (float )spec.getPrecursor().getMz(), (byte )spec.getPrecursor().getCharge()));
      // byte[] family, byte[] qualifier, byte[] value
      row.add(HBases.FAM_ID,   HBasePeakList.COL_UUID, Bytes.toBytes(spec.getId().toString()));
      row.add(HBases.FAM_PROP, HBasePeakList.COL_MZ,   Bytes.toBytes(sig.getMz()));
      row.add(HBases.FAM_PROP, HBasePeakList.COL_SNR,  Bytes.toBytes(sig.getSNR()));
      row.add(HBases.FAM_PROP, HBasePeakList.COL_Z,    Bytes.toBytes(sig.getCharge()));
      row.add(HBases.FAM_PROP, HBasePeakList.COL_IONS, Bytes.toBytes(sigs.size()));
      // TODO row.add(Bytes.toBytes(HBasePeakList.FAM_MZ), Bytes.toBytes(HBasePeakList.COL_MMOD), Bytes.toBytes(spec.getPrecursor().getCharge()));
      indice.put(row);
    }
    //if (verbose && Tools.isSet(sigs)) System.out.println();
  }

  /** Query the HBases for the candidates. No scoring in this call
   *
   * @param signatures
   */
  @Deprecated
  public static Collection<AnnotatedSpectrum> query(PeakList<PepLibPeakAnnotation> signatures, byte[] spec_type, Tolerance tol, double offset) throws IOException
  {
    if (signatures == null || signatures.getPrecursor()==null) return null;

    // setup the HTable for query
    HConnection conn = getConnection();
    // cells the number of row. Can be very expansive for a large table!!
    HTableInterface tbl = conn.getTable(HBases.TBL_MSMSINDEX);

    Map<UUID, AnnotatedSpectrum> id_candidate = new HashMap<UUID, AnnotatedSpectrum>();
    Collection<Range<Peak>> slices = MIMSL.enumeratePrecursors(tol, 0, signatures.getPrecursor());
    for (Range<Peak> range : slices)
    {
      Scan scan = new Scan(HBasePeakList.row4MsMsIndex(spec_type, (float )(range.lowerEndpoint().getMz()+offset), (byte )range.lowerEndpoint().getCharge()),
                           HBasePeakList.row4MsMsIndex(spec_type, (float )(range.upperEndpoint().getMz()+offset), (byte )range.lowerEndpoint().getCharge()));
      if (signatures.size()>0)
      {
        List<Filter> filters = new ArrayList<Filter>(signatures.size());
        for (int k=0; k<signatures.size(); k++)
        {
          // set the range of row keys
          filters.add(new FilterList(FilterList.Operator.MUST_PASS_ALL,
            new SingleColumnValueFilter(HBases.FAM_PROP, HBasePeakList.COL_MZ, CompareFilter.CompareOp.GREATER_OR_EQUAL, new BinaryComparator(Bytes.toBytes(tol.getMin(signatures.getMz(k))))),
            new SingleColumnValueFilter(HBases.FAM_PROP, HBasePeakList.COL_MZ, CompareFilter.CompareOp.LESS_OR_EQUAL,    new BinaryComparator(Bytes.toBytes(tol.getMax(signatures.getMz(k)))))));
        }
        FilterList filterList = new FilterList(FilterList.Operator.MUST_PASS_ONE, filters);
        scan.setFilter(filterList);
      }
      // go thro the peaks
      ResultScanner scanner = tbl.getScanner(scan);
      for (Result rs = scanner.next(); rs != null; rs = scanner.next()) {
        UUID id = UUID.fromString(
          HBases.getString(rs, HBases.FAM_ID, HBasePeakList.COL_UUID));
        if (id_candidate.get(id) == null) id_candidate.put(id, new AnnotatedSpectrum());
        id_candidate.get(id).add(
          HBases.getDouble(rs, HBases.FAM_PROP, HBasePeakList.COL_MZ),
          HBases.getDouble(rs, HBases.FAM_PROP, HBasePeakList.COL_SNR));
        id_candidate.get(id).setId(id);
        // record the number the indice
        id_candidate.get(id).setIonIndexed(
          HBases.getInt(rs, HBases.FAM_PROP, HBasePeakList.COL_IONS));
        id_candidate.get(id).setIonQueried(signatures.size());
        id_candidate.get(id).setIonMatched(id_candidate.get(id).size());
        id_candidate.get(id).setPrecursor(signatures.getPrecursor());
        id_candidate.get(id).setMzQueried((float )signatures.getPrecursor().getMz());
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
  public static Collection<AnnotatedSpectrum> query(Peak[] precursors, MimslSettings settings, double offset, Peak... frags) throws IOException
  {
    if (frags==null || precursors==null) return null;

    // setup the HTable for query
    HConnection conn = getConnection();
    // cells the number of row. Can be very expansive for a large table!!
    HTableInterface tbl = conn.getTable(HBases.TBL_MSMSINDEX);

    Map<UUID, AnnotatedSpectrum> id_candidate = new HashMap<UUID, AnnotatedSpectrum>();
    Collection<Range<Peak>> slices = MIMSL.enumeratePrecursors(settings.getPrecursorTol(), settings.getZFloat(), precursors);
    for (Range<Peak> range : slices)
    {
      Tools.putAll(id_candidate, fetch(tbl, settings.getSpecType(),
          (float )(range.lowerEndpoint().getMz()+offset), (float )(range.upperEndpoint().getMz()+offset),
          (byte )range.lowerEndpoint().getCharge(), settings.getFragmentTol(), frags));
/*

      Scan scan = new Scan(HBasePeakList.row4MsMsIndex(settings.getSpecType(), range.lowerEndpoint().getMz()+offset, range.lowerEndpoint().getCharge()),
                           HBasePeakList.row4MsMsIndex(settings.getSpecType(), range.upperEndpoint().getMz()+offset, range.lowerEndpoint().getCharge()));
      if (frags.length>0)
      {
        List<Filter> filters = new ArrayList<Filter>(frags.length);
        for (int k=0; k<frags.length; k++)
        {
          // set the range of row keys
          filters.add(new FilterList(FilterList.Operator.MUST_PASS_ALL,
              new SingleColumnValueFilter(HBases.FAM_PROP, HBasePeakList.COL_MZ, CompareFilter.CompareOp.GREATER_OR_EQUAL,
                  new BinaryComparator(Bytes.toBytes(settings.getFragmentTol().getMin(frags[k].getMz())))),
              new SingleColumnValueFilter(HBases.FAM_PROP, HBasePeakList.COL_MZ, CompareFilter.CompareOp.LESS_OR_EQUAL,
                  new BinaryComparator(Bytes.toBytes(settings.getFragmentTol().getMax(frags[k].getMz()))))));
        }
        FilterList filterList = new FilterList(FilterList.Operator.MUST_PASS_ONE, filters);
        scan.setFilter(filterList);
      }
      // go thro the peaks
      ResultScanner scanner = tbl.getScanner(scan);
      for (Result rs = scanner.next(); rs != null; rs = scanner.next()) {
        UUID id = UUID.fromString(
            HBases.getString(rs, HBases.FAM_ID, HBasePeakList.COL_UUID));
        if (id_candidate.cells(id) == null) id_candidate.put(id, new AnnotatedSpectrum());
        id_candidate.cells(id).add(
            HBases.getDouble(rs, HBases.FAM_PROP, HBasePeakList.COL_MZ),
            HBases.getDouble(rs, HBases.FAM_PROP, HBasePeakList.COL_SNR));
        id_candidate.cells(id).setId(id);
        // record the number the indice
        id_candidate.cells(id).setIonIndexed(
            HBases.getInt(rs, HBases.FAM_PROP, HBasePeakList.COL_IONS));
        id_candidate.cells(id).setIonQueried(frags.length);
        id_candidate.cells(id).setIonMatched(id_candidate.cells(id).size());
        id_candidate.cells(id).setPrecursor(null);
        id_candidate.cells(id).setPrecursors(precursors);
//        id_candidate.cells(id).setMzQueried((float )signatures.getPrecursor().getMz());
      }
*/
      System.out.println("Candidates with the precursor m/z range of " + range.toString());
      for (UUID id : id_candidate.keySet())
      {
        id_candidate.get(id).setPrecursors(precursors);
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
  public static Map<UUID, AnnotatedSpectrum> fetch(HTableInterface tbl,
     byte[] stype, float mzlow, float mzhigh, byte z, Tolerance frag, Peak... frags) throws IOException
  {
    Scan scan = new Scan(HBasePeakList.query4MsMsIndex(stype, mzlow, z), HBasePeakList.query4MsMsIndex(stype, mzhigh, z));
    if (frags!=null && frags.length>0)
    {
      List<Filter> filters = new ArrayList<Filter>(frags.length);
      for (int k=0; k<frags.length; k++)
      {
        // set the range of row keys
        filters.add(new FilterList(FilterList.Operator.MUST_PASS_ALL,
            new SingleColumnValueFilter(HBases.FAM_PROP, HBasePeakList.COL_MZ, CompareFilter.CompareOp.GREATER_OR_EQUAL,
                new BinaryComparator(Bytes.toBytes(frag.getMin(frags[k].getMz())))),
            new SingleColumnValueFilter(HBases.FAM_PROP, HBasePeakList.COL_MZ, CompareFilter.CompareOp.LESS_OR_EQUAL,
                new BinaryComparator(Bytes.toBytes(frag.getMax(frags[k].getMz()))))));
      }
      FilterList filterList = new FilterList(FilterList.Operator.MUST_PASS_ONE, filters);
      scan.setFilter(filterList);
    }
    // go thro the peaks
    ResultScanner scanner = tbl.getScanner(scan);
    Map<UUID, AnnotatedSpectrum> id_candidate = new HashMap<UUID, AnnotatedSpectrum>();
    for (Result rs = scanner.next(); rs != null; rs = scanner.next()) {
      UUID id = UUID.fromString(
          HBases.getString(rs, HBases.FAM_ID, HBasePeakList.COL_UUID));
      if (id_candidate.get(id) == null) id_candidate.put(id, new AnnotatedSpectrum());
      id_candidate.get(id).add(
          HBases.getDouble(rs, HBases.FAM_PROP, HBasePeakList.COL_MZ),
          HBases.getDouble(rs, HBases.FAM_PROP, HBasePeakList.COL_SNR));
      id_candidate.get(id).setId(id);
      // record the number the indice
      id_candidate.get(id).setIonIndexed(
          HBases.getInt(rs, HBases.FAM_PROP, HBasePeakList.COL_IONS));
      id_candidate.get(id).setIonQueried(frags.length);
      id_candidate.get(id).setIonMatched(id_candidate.get(id).size());
//      id_candidate.cells(id).setPrecursor(null);
      // assume the range is symmetrical for now
      id_candidate.get(id).setMzQueried(0.5f*(mzlow+mzhigh));
    }
/*
    System.out.println("Candidates with the precursor m/z range of " + range.toString());
    for (UUID id : id_candidate.keySet())
    {
      System.out.print(
          Tools.d2s(id_candidate.cells(id).getPrecursor().getMz(), 3) + ", +" + id_candidate.cells(
              id).getPrecursor().getCharge() + "(" + id_candidate.cells(id).size() + "): ");
      for (int i = 0; i < id_candidate.cells(id).size(); i++) {
        System.out.print(Tools.d2s(id_candidate.cells(id).getMz(i), 2) + ",");
      }
      System.out.println();
    }
*/
    return id_candidate;
  }
/*
  public static long prepareMsp(File src, byte[] spec_type, double half_width, double min_mz, int min_pk, double min_snr) throws IOException, URISyntaxException
  {
    // connection to the cluster
    HConnection conn = HConnectionManager.createConnection(conf);

    HTableInterface peaklist = conn.getTable(HBasePeakList.TBL_PEAKLIST),
                      indice = conn.getTable(HBases.TBL_MSMSINDEX);

    System.out.println("Preparing " + src.getAbsolutePath());

    BufferedReader reader = new BufferedReader(new FileReader(src));
    MsLibReader msp = new MsLibReader(reader, URIBuilder.UNDEFINED_URI,
      PeakList.Precision.FLOAT,
      new MspCommentParser(), new MspAnnotationResolver2(),
      Pattern.compile("^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+\"([^\"]+)\"$"),
      new PeakProcessorChain<PepLibPeakAnnotation>());

    long counts=0; PeptideConsensusSpectrum spec = null;
    while (msp.hasNext())
    {
      if (++counts%1000  ==0) System.out.print(".");
      if (  counts%200000==0) System.out.print("\n");
      spec = msp.next(); spec.setSpectrumSource(new URI(src.getName()));
      index(peaklist, indice, spec_type, spec, MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr));
      spec = null;
    }
    peaklist.close(); indice.close(); conn.close(); // release resources

    return counts;
  }
*/
  public static MsLibReader newLibReader(String src)
  {
    MsLibReader splib=null;
    try
    {
      if (src.endsWith(".sptxt"))
      {
        splib = SptxtReader.newBuilder(new BufferedReader(new FileReader(src)),
            URIBuilder.UNDEFINED_URI, PeakList.Precision.FLOAT).build();

        // work around the unknown mod per Oliver's suggestion
        Field f = MsLibReader.class.getDeclaredField("modResolver");
        f.setAccessible(true);
        UnimodModificationResolver modResolver = (UnimodModificationResolver) f.get(splib);
        modResolver.putOverrideUnimod("USM_C", Modification.parseModification("H"));
        modResolver.putOverrideUnimod("Propionamide:13C(3)", Modification.parseModification("H"));
        modResolver.putOverrideUnimod("USM_n_230.170762", Modification.parseModification("H"));
        modResolver.putOverrideUnimod("USM_K_357.257892", Modification.parseModification("H"));
//      modResolver.putTranslate("USM_C", "Unimod version of USM_C");
//      modResolver.putTranslate("Propionamide:13C(3)", "Unknown");
      }
      else if (src.endsWith(".msp"))
      {
        splib = new MsLibReader(new BufferedReader(new FileReader(src)), URIBuilder.UNDEFINED_URI,
            PeakList.Precision.FLOAT,
            new MspCommentParser(), new MspAnnotationResolver2(),
            new PeakProcessorChain<PepLibPeakAnnotation>(),
            new NumModResolver("^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+\"([^\"]+)\"$"));
      }
    }
    catch (IOException ioe)
    {
      throw new RuntimeException(ioe);
    }
    catch (NoSuchFieldException nse)
    {
      throw new RuntimeException(nse);
    }
    catch (IllegalAccessException iae)
    {
      throw new RuntimeException(iae);
    }
    return splib;
  }
/*
  public static long prepareLib(File src, byte[] spec_type, double half_width, double min_mz, int min_pk, double min_snr) throws IOException
  {
    // connection to the cluster
    HConnection conn = HConnectionManager.createConnection(conf);

    System.out.println("Preparing " + src.getAbsolutePath());
    HTableInterface peaklist = conn.getTable(HBasePeakList.TBL_PEAKLIST),
                      indice = conn.getTable(HBases.TBL_MSMSINDEX);
    MsLibReader        splib = newLibReader(src);
*/
/*
    BufferedReader reader = new BufferedReader(new FileReader(src));
    MsLibReader     splib = new SptxtReader(reader, URIBuilder.UNDEFINED_URI, PeakList.Precision.FLOAT);

    try
    {
      // work around the unknown mod per Oliver's suggestion
      Field f = MsLibReader.class.getDeclaredField("modResolver");
      f.setAccessible(true);
      UnimodModificationResolver modResolver = (UnimodModificationResolver) f.cells(splib);
      modResolver.putOverrideUnimod("USM_C", Modification.parseModification("H"));
      modResolver.putOverrideUnimod("Propionamide:13C(3)", Modification.parseModification("H"));
      modResolver.putOverrideUnimod("USM_n_230.170762", Modification.parseModification("H"));
      modResolver.putOverrideUnimod("USM_K_357.257892", Modification.parseModification("H"));
//      modResolver.putTranslate("USM_C", "Unimod version of USM_C");
//      modResolver.putTranslate("Propionamide:13C(3)", "Unknown");
    }
    catch (NoSuchFieldException nse)
    {
      throw new RuntimeException(nse);
    }
    catch (IllegalAccessException iae)
    {
      throw new RuntimeException(iae);
    }
*//*


    long counts=0; PeptideConsensusSpectrum spec = null;
    while (splib.hasNext())
    {
      if (++counts%1000  ==0) System.out.print(".");
      if (  counts%200000==0) System.out.print("\n");
      try
      {
        spec = splib.next(); spec.setSpectrumSource(new URI(src.getName()));
        index(peaklist, indice, spec_type, spec, MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr));
        spec = null;
      }
      catch (URISyntaxException ue)
      {
        ue.printStackTrace();
      }
      catch (IllegalStateException e)
      {
        System.out.print("undefined residue!");
        //e.printStackTrace();
      }
    }
    peaklist.close(); indice.close(); conn.close(); splib.close(); // release resources

    return counts;
  }
*/
/*
  public static long prepareLibs(String root, byte[] spec_type, double half_width, double min_mz, int min_pk, double min_snr, String... msps) throws IOException, URISyntaxException
  {
    long counts = 0;
    for (String msp : msps)
    {
      if      (msp.endsWith(".msp"))   counts += prepareMsp(  new File(root + "/" + msp), spec_type, half_width, min_mz, min_pk, min_snr);
      else if (msp.endsWith(".sptxt")) counts += prepareSpLib(new File(root + "/" + msp), spec_type, half_width, min_mz, min_pk, min_snr);
//
//      counts += prepareMsp(new File(root + "/" + msp), spec_type, half_width, min_mz, min_pk, min_snr);
    }
    return counts;
  }
*/
  public static long prepareLib(String root, HBaseSpLib lib, double half_width, double min_mz, int min_pk, double min_snr) throws IOException, URISyntaxException
  {
    // connection to the cluster
    HConnection         conn = getConnection();
    String               src = root+"/"+lib.getName()+(lib.isFormat(HBaseSpLib.LIB_MSP)?".msp":(lib.isFormat(HBaseSpLib.LIB_SPTXT)?".sptxt":""));
    HTableInterface peaklist = conn.getTable(HBasePeakList.TBL_PEAKLIST),
                      indice = conn.getTable(HBases.TBL_MSMSINDEX);
    MsLibReader        splib = newLibReader(src);

    System.out.println("Preparing " + src);

    long counts=0; PeptideConsensusSpectrum spec = null;
    while (splib.hasNext())
    {
      if (++counts%1000  ==0) System.out.print(".");
      if (  counts%200000==0) System.out.print("\n");
      try
      {
        spec = splib.next(); spec.setSpectrumSource(new URI(src));
        index(peaklist, indice, lib.getSpecType(), spec, MIMSL.toSignature(spec, half_width, min_mz, min_pk, min_snr));
        spec = null;
      }
      catch (URISyntaxException ue)
      {
        ue.printStackTrace();
      }
      catch (IllegalStateException e)
      {
        System.out.print("undefined residue!");
        //e.printStackTrace();
      }
    }
    peaklist.close(); indice.close(); conn.close(); splib.close(); // release resources
    // cells the entries counts
    HBaseSpLib.increEntries(Bytes.toBytes(lib.getName()), counts);
    return counts;
  }
  public static Collection<AnnotatedSpectrum> loadPeakLists(Collection<AnnotatedSpectrum> spectra) throws IOException
  {
    HConnection    conn = getConnection();
    HTableInterface tbl = conn.getTable(HBasePeakList.TBL_PEAKLIST);
    for (AnnotatedSpectrum spec : spectra)
    {
      HBasePeakList.getPeakList(tbl, spec.getId()).toPeakList(spec);
    }
    tbl.close();

    return spectra;
  }
  public static UUID randUUID(HTableInterface tbl, byte[] spectype, Random rand) throws IOException
  {
    if (rand==null) rand=new Random(System.nanoTime());

    float mz = rand.nextFloat() * 1000f+350f; int z = rand.nextInt(3)+1;
    Map<UUID, AnnotatedSpectrum> id_candidate = fetch(tbl, spectype, mz-0.015f, mz+0.015f, (byte )z, null);

    return Tools.front(id_candidate.keySet());
  }
  public static Map<PeptideConsensusSpectrum.Status, PeptideConsensusSpectrum> sampleRecovery(String src, int nsig, int interval, MimslSettings settings) throws IOException
  {
    // connection to the cluster
    HConnection          conn = HBases.getConnection();
    MsLibReader         splib = HBaseProteomics.newLibReader(src);
    Random               rand = new Random(System.nanoTime());
    PeptideConsensusSpectrum      spec = null;
    List<AnnotatedPeak> peaks = null;
    Peak[]            sampled = new Peak[nsig];

    System.out.println("Sampling " + src);

    long counts=0;
    Map<PeptideConsensusSpectrum.Status, PeptideConsensusSpectrum> findings = new HashMap<PeptideConsensusSpectrum.Status, PeptideConsensusSpectrum>();
    while (splib.hasNext())
    {
      if (++counts%1000  ==0) System.out.print(".");
      if (  counts%200000==0) System.out.print("\n");
      try
      {
        spec = splib.next();
        if (counts%interval==0 && spec.getPeptide()!=null && spec.getPeptide().getModificationCount()==0)
        {
          peaks = MIMSL.toSignature(spec, settings.getHalfWidth(), settings.getMinMz(), settings.getMinPeaks(), settings.getMinSNR());
          // extract a subset of signature peaks
          int n=0;
          for (int i=0; i<peaks.size(); i++)
          {
            int k=rand.nextInt(peaks.size());
            if (!Tools.contains(sampled, peaks.get(k), nsig)) sampled[n++] = peaks.get(k);
            if (n>=nsig) break;
          }
          // issue the MIMSL query using the sampled signature peaks
          List<AnnotatedSpectrum> candidates = MIMSL.run(new Peak[] {spec.getPrecursor()}, settings, n<nsig?Arrays.copyOfRange(sampled, 0, nsig):sampled);
          // locate the candidate that agree with the starter
          int rank=-1;
          if (Tools.isSet(candidates))
            for (int i=0; i<candidates.size(); i++)
              if (candidates.get(i).getPeptide().equals(spec.getPeptide())) rank=i;

          spec.setMsLevel(rank); // save the ranking
          spec.setStatus(rank==0?PeptideConsensusSpectrum.Status.SINGLETON:(rank<0? PeptideConsensusSpectrum.Status.UNKNOWN: PeptideConsensusSpectrum.Status.CONFLICTING_ID));
          findings.put(spec.getStatus(), spec);
        }
        spec=null; peaks=null;
      }
      catch (IllegalStateException e)
      {
        System.out.print("undefined residue!");
        //e.printStackTrace();
      }
    }
    conn.close(); splib.close(); // release resources

    return findings;
  }
}
