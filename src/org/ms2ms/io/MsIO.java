package org.ms2ms.io;

import org.expasy.mzjava.core.io.ms.spectrum.MgfWriter;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.*;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.model.mzml.*;

import java.io.*;
import java.util.*;

/** Responsible for writing out the MS objects in binary format. Not to be confused with import duty
 *
 * User: wyu
 * Date: 9/29/14
 */
public class MsIO extends IOs
{
  public static MsSpectrum read(RandomAccessFile w, long offset)
  {
    try
    {
      w.seek(offset);
      return read(w);
    }
    catch(IOException i)
    {
      i.printStackTrace();
    }
    return null;
  }

  public static MsSpectrum read(RandomAccessFile w)
  {
    try
    {
      byte[] bs = new byte[w.readInt()]; w.read(bs);
      MsSpectrum spec = MsSpectrum.fromBytes(bs); bs=null;
      return spec;
    }
    catch(IOException i)
    {
//      i.printStackTrace();
    }
    return null;
//    catch(ClassNotFoundException c)
//    {
//      System.out.println("Employee class not found");
//      c.printStackTrace();
//    }
  }
  public static long write(RandomAccessFile w, MsSpectrum ms)
  {
    try
    {
      long p1 = w.getFilePointer();
      byte[] bs = MsSpectrum.toBytes(ms);
      w.writeInt(bs.length); w.write(bs); bs=null;
      return p1;
    }
    catch (IOException e)
    { throw new RuntimeException("Error during persistence", e); }
  }

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
    write(w, (PeakList) ms);

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
  public static List<MsnSpectrum> readSpectra(String s) { return readSpectra(s, null); }
  public static List<MsnSpectrum> readSpectra(String s, long[] offsets)
  {
    RandomAccessFile F = null;
    try
    {
      try
      {
        F = new RandomAccessFile(s, "r");
        List<MsnSpectrum> spec = Tools.isSet(offsets) ? readSpectra(F, offsets) : readSpectra(F);
        F.close();
        return spec;
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
    catch (IOException ie) { throw new RuntimeException("Error while reading the spectra", ie); }
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s) throws IOException
  {
    if (s==null) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>();
    try
    {
      while (1==1)
      {
        spectra.add(read(s).toMsnSpectrum());
      }
    }
    catch (Exception ie) {}
    return spectra;
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s, long[] offsets) throws IOException
  {
    if (s==null || !Tools.isSet(offsets)) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>(offsets.length);
    for (long offset : offsets)
    {
      s.seek(offset);
      spectra.add(read(s).toMsnSpectrum());
    }
    return spectra;
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s, Map<Long, String> offset_row) throws IOException
  {
    if (s==null || !Tools.isSet(offset_row)) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>(offset_row.size());
    for (Long offset : offset_row.keySet())
    {
      s.seek(offset);
      MsnSpectrum spec = read(s).toMsnSpectrum();
      spec.setComment(offset_row.get(offset)); spectra.add(spec);
    }
    return spectra;
  }
  public static void writeSpectra(String s, Collection<MsnSpectrum> spectra) throws IOException
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
          write(F, MsSpectrum.adopt(m));
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
  public static String get(Collection<CVParam> params, String key)
  {
    if (Tools.isSet(params) && Strs.isSet(key))
      for (Iterator lCVParamIterator = params.iterator(); lCVParamIterator.hasNext();)
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals(key))
        {
          return lCVParam.getValue().trim();
        }
      }
    return null;
  }
  public static Integer getInt(   Collection<CVParam> params, String key) { return Stats.toInt(   get(params, key)); }
  public static Double  getDouble(Collection<CVParam> params, String key) { return Stats.toDouble(get(params, key)); }

  public static String getIf(Collection<CVParam> params, String tag_required, String val_required)
  {
    if (Tools.isSet(params) && Strs.isSet(tag_required) && Strs.isSet(val_required))
      for (Iterator lCVParamIterator = params.iterator(); lCVParamIterator.hasNext();)
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals(tag_required))
        {
          String unitRT = lCVParam.getUnitAccession().trim();
          if (unitRT.equals(val_required)) //... Validating rt unit (mins or secs) ...//
          {
            return lCVParam.getValue().trim();
          }
        }
      }
    return null;
  }
  public static boolean hasAccession(Collection<CVParam> params, String tag_required)
  {
    if (Tools.isSet(params) && Strs.isSet(tag_required))
      for (Iterator lCVParamIterator = params.iterator(); lCVParamIterator.hasNext();)
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals(tag_required)) return true;
      }
    return false;
  }
  public static Integer ScanNumberFromSpectrumRef(String ref)
  {
    String[] ids = Strs.split(ref, ' ', true);
    for (String id : ids)
      if (id.indexOf("scan=")==0) return Stats.toInt(id.substring(5));

    return null;
  }
  // write the content of the spectrum from an mzML file to another MGF file
  public static void writeMGF(MgfWriter w, uk.ac.ebi.jmzml.model.mzml.Spectrum ms) throws IOException
  {
    if (w==null || ms==null) return;




  }
}
