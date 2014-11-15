package org.ms2ms.nosql;

import org.apache.hadoop.hbase.client.Get;
import org.apache.hadoop.hbase.client.HTableInterface;
import org.apache.hadoop.hbase.client.Put;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.peaklist.DoublePeakList;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.UUID;

/** An Entity class for I/O with HBase database. Not to be used directly for data analysis
 *  Only expose the PeakList for public consumption
 *
 * Created by wyu on 4/18/14.
 */
public final class HBasePeakList extends MsSpectrum
{
  private static final long serialVersionUID = 8472732522296541667L;

  // Try to keep the ColumnFamily names as small as possible, preferably one character (e.g. "d" for data/default).
  // https://hbase.apache.org/book/rowkey.design.html
  static public byte[] FAM_PRECURSOR = Bytes.toBytes("p");
  static public String TBL_PEAKLIST  = "PeakList";

  static public byte[] COL_MZ        = Bytes.toBytes("mz");
  static public byte[] COL_Z         = Bytes.toBytes( "z");
  //static final byte[] COL_AI        = "ai".getBytes();
  static public byte[] COL_IONS      = Bytes.toBytes("io");
  static public byte[] COL_UUID      = Bytes.toBytes("id");
  static public byte[] COL_MMOD      = Bytes.toBytes("mm"); // mass of the modification on the fragment
  static public byte[] COL_SIG       = Bytes.toBytes("sg"); // m/z of the signature fragment
  static public byte[] COL_SNR       = Bytes.toBytes("sr"); // m/z of the signature fragment

  static public byte[] SPEC_TRAP_CID = Bytes.toBytes("c");
  static public byte[] SPEC_TRAP_HCD = Bytes.toBytes("h");
  static public byte[] SPEC_TRAP_ETD = Bytes.toBytes("e");
  static public byte[] SPEC_QTOF     = Bytes.toBytes("q");

  protected int cursor;            // current position of the peak
  private   byte[] ppmList, flucList; // the variances of the peaks in m/z and intensity
  private   String peptide, protein;

  protected float maxPPM;
  protected float maxFluc;

  public <A extends PeakAnnotation> HBasePeakList(PeakList<A> src)
  {
    super(src);
    // the stats

    if (src instanceof LibrarySpectrum)
    {
      peptide = ((LibrarySpectrum )src).getPeptide().toString().replace("(Carbamidomethyl)", "");
      protein = Tools.front(((LibrarySpectrum) src).getProteinAccessionNumbers());
    }
/*

    // TODO copy the data from the source
    precursorMz =        src.getPrecursor().getMz();
    precursorAi =(float )src.getPrecursor().getIntensity();
    precursorZ  =        src.getPrecursor().getCharge();
    // the stats
    maxIntensity=(float )src.getBasePeakIntensity();
    size        =        src.size();

    if (src instanceof LibrarySpectrum)
    {
      peptide = ((LibrarySpectrum )src).getPeptide().toString().replace("(Carbamidomethyl)", "");
      protein = Tools.front(((LibrarySpectrum) src).getProteinAccessionNumbers());
    }

    mzList        = new double[src.size()];
    intensityList = new short[ src.size()];
    for (int i=0; i<src.size(); i++)
    {
      mzList[       i]=src.getMz(i);
      intensityList[i]=(short )(src.getIntensity(i)*Short.MAX_VALUE/maxIntensity+Short.MIN_VALUE);
    }
*/
  }

  public <A extends PeakAnnotation> DoublePeakList<A> toPeakList()
  {
    DoublePeakList<A> peaks = new DoublePeakList<A>(size);
    // initialize the content of the new list
    peaks.setPrecursor(new Peak(precursorMz, precursorAi, precursorZ));
    for (int i=0; i<size; i++)
    {
      peaks.add(getMz(i), getIntensity(i));
    }

    return peaks;
  }
  public <A extends PeakAnnotation> PeakList<A> toPeakList(PeakList<A> peaks)
  {
    peaks.clear();
    // initialize the content of the new list
    peaks.setPrecursor(new Peak(precursorMz, precursorAi, precursorZ));
    for (int i=0; i<size; i++)
      peaks.add(getMz(i), getIntensity(i));

    if (peaks instanceof AnnotatedSpectrum)
    {
      // can't set the peptide so we have to put it into the comment
      // TODO need to setup a proper graph to store the relationship
      ((AnnotatedSpectrum)peaks).setComment(peptide + "^" + protein);
    }
    return peaks;
  }

  public static HBasePeakList fromBytes(byte[] s)
  {
    try
    {
      ByteArrayInputStream bai = new ByteArrayInputStream(s);
      ObjectInputStream     in = new ObjectInputStream(bai);
      HBasePeakList          e = (HBasePeakList) in.readObject();
      in.close(); bai.close();
      return e;
    }
    catch(IOException i)
    {
      i.printStackTrace();
    }
    catch(ClassNotFoundException c)
    {
      System.out.println("Employee class not found");
      c.printStackTrace();
    }
    return null;
  }
  public static byte[] toBytes(HBasePeakList s)
  {
    ByteArrayOutputStream bao = new ByteArrayOutputStream();
    try
    {
      ObjectOutputStream out = new ObjectOutputStream(bao);
      out.writeObject(s);
      out.close(); bao.close();
    }
    catch (IOException e)
    { throw new RuntimeException("Error during persistence", e); }

    return bao.toByteArray();
  }

  public static void save(HTableInterface tbl, PeakList spec) throws IOException
  {
    if (spec==null) return;

    Put row = new Put(row4PeakList(spec.getId()));
    // byte[] family, byte[] qualifier, byte[] value
    row.add(FAM_PRECURSOR, COL_MZ, Bytes.toBytes(spec.getPrecursor().getMz()));
    row.add(FAM_PRECURSOR, COL_Z,  Bytes.toBytes(spec.getPrecursor().getCharge()));
    // create the peaks bytes
    row.add(HBase.FAM_PROP, COL_IONS, toBytes(new HBasePeakList(spec)));

    tbl.put(row);
  }
  public static HBasePeakList getPeakList(HTableInterface tbl, UUID id) throws IOException
  {
    Get g = new Get(row4PeakList(id));
    byte[] value = tbl.get(g).getValue(HBase.FAM_PROP, COL_IONS);
    HBasePeakList peaks = HBasePeakList.fromBytes(value);

    return peaks;
  }
  public static byte[] row4PeakList(UUID id) { return Bytes.toBytes(id.toString()); }
  public static byte[] row4MsMsIndex(byte[] spec_type, float mz, byte z)
  {
    // tag the system time in nanosec to ensure unique row key
    return Bytes.add(query4MsMsIndex(spec_type, mz,z), Bytes.toBytes(System.nanoTime()));
  }
  public static byte[] query4MsMsIndex(byte[] spec_type, float mz, byte z)
  {
    // tag the system time in nanosec to ensure unique row key
    return Bytes.add(spec_type, Bytes.add(Bytes.toBytes(mz), new byte[] {z}));
  }
//  public static byte[] row4MsMsIndex(byte[] spec_type, double mz, int z)
//  {
//    return row4MsMsIndex(spec_type, (float )mz, (byte )z);
//  }
}
