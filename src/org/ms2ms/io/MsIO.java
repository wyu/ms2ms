package org.ms2ms.io;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.*;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/** Responsible for writing out the MS objects in binary format. Not to be confused with import duty
 *
 * User: wyu
 * Date: 9/29/14
 */
public class MsIO extends IOs
{
/*
  public static void read(RandomAccessFile w)
  {
    try
    {
      MsSpectrum
      ByteArrayInputStream bai = new ByteArrayInputStream(s);
      ObjectInputStream     in = new ObjectInputStream(w);
      MsSpectrum          e = (MsSpectrum ) in.readObject();
      in.close(); bai.close();
      return;
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
    return;
  }
  public static long write(ObjectOutputStream w, MsnSpectrum ms)
  {
    try
    {
      long pos = w.;
      MsSpectrum m = new MsSpectrum(ms);
      // Write object out to disk
      w.writeObject (m);
      Tools.dispose( m);
    }
    catch (IOException e)
    { throw new RuntimeException("Error during persistence", e); }
  }
*/

  // BufferedWriter
  public static void write(DataOutput w, PeakList ms) throws IOException
  {
    w.writeDouble(ms.getPrecursor().getMz());
    w.writeDouble(ms.getPrecursor().getIntensity());
    w.writeInt(   ms.getPrecursor().getCharge());

    w.writeInt(   ms.size());
    if (ms.size()>0)
    {
      for (int i=0; i<ms.size(); i++)
      {
        w.writeDouble(ms.getMz(i));
        w.writeDouble(ms.getIntensity(i));
      }
    }
  }
  public static void write(DataOutput w, MsnSpectrum ms) throws IOException
  {
    write(w, (PeakList)ms);

    int counts = ms.getRetentionTimes()!=null&&ms.getRetentionTimes().size()>0?ms.getRetentionTimes().size():0;
    w.writeInt(counts);
    if (counts>0)
      for (RetentionTime rt : ms.getRetentionTimes()) w.writeDouble(rt.getTime());
  }
  public static PeakList read(DataInput w, PeakList ms) throws IOException
  {
    if (ms!=null)
    {
      ms.setPrecursor(new Peak(w.readDouble(), w.readDouble(), w.readInt()));

      int counts = w.readInt(); ms.clear();
      if (counts>0)
        for (int i=0; i<counts; i++)
          ms.add(w.readDouble(), w.readDouble());
    }
    return ms;
  }
  public static MsnSpectrum read(DataInput w, MsnSpectrum ms) throws IOException
  {
    if (ms!=null)
    {
      ms = (MsnSpectrum )read(w, (PeakList)ms);

      int counts = w.readInt();
      // clear the content of the retention times if necessary
      if (ms.getRetentionTimes()!=null) ms.getRetentionTimes().clear();
      if (counts>0)
        for (int i=0; i<counts; i++)
          ms.addRetentionTime(new RetentionTimeDiscrete(w.readDouble(), TimeUnit.SECOND));
    }
    return ms;
  }
  public static List<MsnSpectrum> readSpectra(String s, long[] offsets)
  {
    RandomAccessFile F = null;
    try
    {
      try
      {
        F = new RandomAccessFile(s, "r");
        List<MsnSpectrum> spec = readSpectra(F, offsets);
        F.close();
        return spec;
      }
      catch (FileNotFoundException fne)
      {
        throw new RuntimeException("Not able to location the file: " + s, fne);
      }
      finally
      {
        if (F!=null) F.close();
      }
    }
    catch (IOException ie) { throw new RuntimeException("Error while reading the spectra", ie); }
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s, long[] offsets) throws IOException
  {
    if (s==null || !Tools.isSet(offsets)) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>(offsets.length);
    for (long offset : offsets)
    {
      s.seek(offset);
      spectra.add(read(s, new MsnSpectrum()));
    }
    return spectra;
  }
  public static void writeSpectra(String s, List<MsnSpectrum> spectra) throws IOException
  {
    if (s==null || !Tools.isSet(spectra)) return;
    // the output
    RandomAccessFile F = null;
    try
    {
      try
      {
        F = new RandomAccessFile(s, "rw");
        for (MsnSpectrum m : spectra)
        {
          write(F, m);
        }
        F.close();
        return;
      }
      catch (FileNotFoundException fne)
      {
        throw new RuntimeException("Not able to locate the file: " + s, fne);
      }
      finally
      {
        if (F!=null) F.close();
      }
    }
    catch (IOException ie) { throw new RuntimeException("Error while writing the spectra", ie); }
  }
}
