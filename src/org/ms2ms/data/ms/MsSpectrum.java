package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeDiscrete;
import org.expasy.mzjava.core.ms.spectrum.ScanNumberDiscrete;
import org.expasy.mzjava.core.ms.spectrum.TimeUnit;
import org.ms2ms.Disposable;

import java.io.*;

/** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:  Dedicated class to serve the binary I/O
 * <p/>
 * Author: wyu
 * Date:   11/13/14
 */
public class MsSpectrum  implements Serializable, Disposable
{
  private static final long serialVersionUID = 8472752523296641667L;

  protected int precursorZ, msLevel;   // the precursor charge
  protected float[] mzList;            // m/z of the peaks
  protected float[] intensityList;     // relative intensities of the peaks
  protected float precursorAi;
  protected double precursorMz;         // the precursor charge
  protected int size, scan;              // length of the peak list
  // the upper bound of the peak intensity and variance in m/z and intensity
  protected float maxIntensity, rt;
//  protected String scan;

  public <A extends PeakAnnotation> MsSpectrum(PeakList<A> src)
  {
    super(); init(src);
//    // TODO copy the data from the source
//    precursorMz =        src.getPrecursor().getMz();
//    precursorAi =(float )src.getPrecursor().getIntensity();
//    precursorZ  =        src.getPrecursor().getCharge();
//    maxIntensity=(float )src.getBasePeakIntensity();
//
//    mzList        = new float[src.size()];
//    intensityList = new float[ src.size()];
//    for (int i=0; i<src.size(); i++)
//    {
//      mzList[       i]=(float )src.getMz(i);
//      intensityList[i]=(float )src.getIntensity(i);
//    }
//    size = src.size();
  }
  public MsSpectrum(MsnSpectrum src)
  {
    super(); init(src);
    scan    = src.getScanNumbers().getFirst().getValue();
    msLevel = src.getMsLevel();
    // assume the time is in seconds
    if (src.getRetentionTimes()!=null)
      rt    = (float )src.getRetentionTimes().getFirst().getTime() / 60f;
  }

  private <A extends PeakAnnotation> void init(PeakList<A> src)
  {
    // TODO copy the data from the source
    precursorMz =        src.getPrecursor().getMz();
    precursorAi =(float )src.getPrecursor().getIntensity();
    precursorZ  =        src.getPrecursor().getCharge();

    if (src.size()>0) maxIntensity=(float )src.getBasePeakIntensity();

    mzList        = new float[src.size()];
    intensityList = new float[ src.size()];
    for (int i=0; i<src.size(); i++)
    {
      mzList[       i]=(float )src.getMz(i);
      intensityList[i]=(float )src.getIntensity(i);
    }
    size = src.size();
  }
  /** compression and uncompression routines  **/
  protected double getMz(int i)
  {
    rangeCheck(i); return mzList[i];
  }

  protected double getIntensity(int i)
  {
    rangeCheck(i); return intensityList[i];
//    return maxIntensity * (intensityList[i]-Short.MIN_VALUE)/(2d*Short.MAX_VALUE);
  }

  private void rangeCheck(int index)
  {
    if (index >= size) throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size);
  }

  @Override
  public void dispose()
  {
    intensityList=null; mzList=null;
  }
  public static MsSpectrum adopt(MsnSpectrum s) { return new MsSpectrum(s); }
  public MsnSpectrum toMsnSpectrum()
  {
    MsnSpectrum out = new MsnSpectrum();
    out.setMsLevel(msLevel);
    out.setPrecursor(new Peak(precursorMz, precursorAi, precursorZ));
    for (int i=0; i<mzList.length; i++)
    {
      out.add(mzList[i], intensityList[i]);
    }
    out.addScanNumber(scan);
    out.addRetentionTime(new RetentionTimeDiscrete(rt, TimeUnit.MINUTE));

    return out;
  }
  public static MsSpectrum fromBytes(byte[] s)
  {
    try
    {
      ByteArrayInputStream bai = new ByteArrayInputStream(s);
      ObjectInputStream in = new ObjectInputStream(bai);
      MsSpectrum          e = (MsSpectrum ) in.readObject();
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
  public static byte[] toBytes(MsSpectrum s)
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
}
