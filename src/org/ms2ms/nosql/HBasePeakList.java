package org.ms2ms.nosql;

import org.apache.hadoop.hbase.client.HTableInterface;
import org.apache.hadoop.hbase.client.Put;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.peaklist.DoublePeakList;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.Peak;
import org.expasy.mzjava.core.ms.spectrum.PeakAnnotation;

import java.io.*;
import java.util.UUID;

/** An Entity class for I/O with HBase database. Not to be used directly for data analysis
 *  Only expose the PeakList for public consumption
 *
 * Created by wyu on 4/18/14.
 */
public final class HBasePeakList implements Serializable
{
  static final String TBL_LIBSPECTRUM  = "Spectrum";

  static final byte[] FAM_PRECURSOR = "precursor".getBytes();
  static final byte[] FAM_STAT      = "stat".getBytes();
  static final byte[] FAM_IONS      = "ions".getBytes();
  static final byte[] FAM_FLAG      = "flag".getBytes();

  static final byte[] COL_MZ        = "mz".getBytes();
  static final byte[] COL_Z         =  "z".getBytes();
  static final byte[] COL_AI        = "ai".getBytes();

  private      int size,              // length of the peak list
                   cursor,            // current position of the peak
                   precursorZ;        // the precursor charge
  private     UUID id;
  private double[] mzList;            // m/z of the peaks
  private  short[] intensityList;     // relative intensities of the peaks
  private   byte[] ppmList, flucList; // the variances of the peaks in m/z and intensity

  // the upper bound of the peak intensity and variance in m/z and intensity
  private float maxIntensity, maxPPM, maxFluc, precursorAi;
  private double precursorMz;         // the precursor charge

  public <A extends PeakAnnotation> HBasePeakList(PeakList<A> src)
  {
    super();
    // TODO copy the data from the source
    precursorMz =        src.getPrecursor().getMz();
    precursorAi =(float )src.getPrecursor().getIntensity();
    precursorZ  =        src.getPrecursor().getCharge();
    // the stats
    maxIntensity=(float )src.getBasePeakIntensity();
    size        =        src.size();

    mzList        = new double[src.size()];
    intensityList = new short[ src.size()];
    for (int i=0; i<src.size(); i++)
    {
      mzList[       i]=src.getMz(i);
      intensityList[i]=(short )(src.getIntensity(i)*Short.MAX_VALUE/maxIntensity+Short.MIN_VALUE);
    }
  }

  /** compression and uncompression routines  **/
  private double getMz(int i)
  {
    rangeCheck(i); return mzList[i];
  }
  private double getIntensity(int i)
  {
    rangeCheck(i);
    return maxIntensity * (intensityList[i]-Short.MIN_VALUE)/(2d*Short.MAX_VALUE);
  }

  private void rangeCheck(int index)
  {
    if (index >= size) throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size);
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

    Put row = new Put(Bytes.toBytes(spec.getId().toString()));
    // byte[] family, byte[] qualifier, byte[] value
    row.addImmutable(FAM_PRECURSOR, COL_MZ, Bytes.toBytes(spec.getPrecursor().getMz()));
    row.addImmutable(FAM_PRECURSOR, COL_Z,  Bytes.toBytes(spec.getPrecursor().getCharge()));

    // create the peaks bytes
    row.addImmutable(FAM_IONS, null, toBytes(new HBasePeakList(spec)));

    tbl.put(row);
  }
}
