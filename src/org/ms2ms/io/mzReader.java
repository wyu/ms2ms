package org.ms2ms.io;

import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Tools;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.Writer;
import java.util.Collection;
import java.util.SortedMap;

/** probably not required. MzxmlReader from mzJava.core.io.ms.spectrum should be preferred
 *
 * User: hliu
 * Date: 8/20/14
 */
public class mzReader
{
  public  static String SCAN                = "scan";
  public  static String NUMBER                = "num";
  public  static String RT                = "retentionTime";

  // Terminology as defined in mzXML,
  public  static String PEAKS               = "peaks";
  public  static String PRECISION           = "precision";
  public  static String BYTE_ORDER          = "byteOrder";
  public  static String PAIR_ORDER          = "pairOrder";
  public  static String PRECURSOR_MZ        = "precursorMz";
  public  static String PRECURSOR_INTENSITY = "precursorIntensity";
  public  static String SCAN_COUNT          = "scanCount";
  public  static String START_TIME          = "startTime";
  public  static String END_TIME            = "endTime";
  public  static String PARENT_FILE         = "parentFile";
  public  static String FILE_NAME           = "fileName";
  public  static String FILE_TYPE           = "fileType";
  public  static String FILE_SHA1           = "fileSha1";
  public  static String CATEGORY            = "category";
  public  static String VALUE               = "value";
  public  static String MS_RUN              = "msRun";
  public  static String MS_INSTRUMENT       = "msInstrument";
  public  static String MS_MANUFACTURER     = "msManufacturer";
  public  static String MS_MODEL            = "msModel";
  public  static String MS_LEVEL            = "msLevel";
  public  static String MS_IONISATION       = "msIonisation";
  public  static String MS_MASS_ANALYZER    = "msMassAnalyzer";
  public  static String MS_DETECTOR         = "msDetector";

  protected XMLStreamReader mReader;

  public mzReader() { super(); }
  public mzReader(XMLStreamReader s) { mReader=s; start(); }
  public mzReader(String s)
  {
    try
    {
      mReader = XMLInputFactory.newInstance().createXMLStreamReader(s, new FileInputStream(s));
      start();
    }
    catch (Exception e)
    {
      throw new RuntimeException("Not able to open file: " + s , e);
    }
  }

  // reading the header of the data file
  public void start()
  {
  }

  // reading the closing element if any
  public void end()
  {
  }

  protected static void XIC(MsnSpectrum ms, double ppm, double dRT,
                     MultiTreeTable<Double, Double, String> rt_mz_row,
                     MultiTreeTable<Double, Double, Peak>   rt_mz_ms1,
                     MultiTreeTable<Double, Double, Peak>   rt_mz_mz1)
  {
    double rt = ms.getRetentionTimes().getFirst().getTime()/60d, mz=ms.getPrecursor().getMz();
    // deposit the ms1 intensity if fell within the bound of a key precursor
    SortedMap<Double, TreeMultimap<Double, String>> slice = rt_mz_row.getData().subMap(rt-dRT,rt+dRT);
    if (Tools.isSet(slice))
    {
      SortedMap<Double, Peak> peaks = Spectra.toPeaks(ms);
      for (Double RT : slice.keySet())
        for (Double mz2 : slice.get(RT).keySet())
        {
          SortedMap<Double, Peak> pks = peaks.subMap(mz2-ppm*mz2/1E6, mz2+ppm*mz2/1E6);
          if (Tools.isSet(pks))
            for (Peak p : pks.values())
            {
              rt_mz_ms1.put(RT, mz2, new Peak(rt,        p.getIntensity()));
              rt_mz_mz1.put(RT, mz2, new Peak(p.getMz(), p.getIntensity()));
            }
        }

      peaks = (SortedMap )Tools.dispose(peaks);
    }

  }
  public void writeXIC(Writer w, MultiTreeTable<Double, Double, Peak> rt_mz_ms1) throws IOException
  {
    if (Tools.isSet(rt_mz_ms1))
      for (Double rt : rt_mz_ms1.getData().keySet())
        for (Double mz2 : rt_mz_ms1.row(rt).keySet())
          for (Peak p : rt_mz_ms1.get(rt, mz2))
            w.write(rt+"\t"+mz2+"\t"+p.getMz()+"\t"+p.getIntensity()+"\n");
  }
  protected static Dataframe detectXIC(
      Dataframe out, String row, String mz_col, String RT_col,
      MultiTreeTable<Double, Double, Peak> rt_mz_ms1, MultiTreeTable<Double, Double, Peak> rt_mz_mz1)
  {
    Collection<Peak> XIC = rt_mz_ms1.get(out.getDouble(row, RT_col), out.getDouble(row, mz_col)),
                  XIC_mz = rt_mz_mz1.get(out.getDouble(row, RT_col), out.getDouble(row, mz_col));
    if (Tools.isSet(XIC) && XIC.size()>2)
    {
      double rt_c = Peaks.centroid(XIC);
      out.put(row, "XIC.centroid", rt_c);
      out.put(row, "XIC.area",     Peaks.AbsIntensitySum(XIC));
      out.put(row, "XIC.mz",       Peaks.centroid(XIC_mz));
      out.put(row, "XIC.N",        XIC.size());
      // locate the boundary
      Peak p0=null; boolean found=false;
      for (Peak p : XIC)
      {
        if      (p!=null && p0!=null && p.getMz()<rt_c && !found) { out.put(row, "XIC.left",   p.getMz()); found=true; }
        else if (p==null && p0!=null && p.getMz()>rt_c)           { out.put(row, "XIC.right", p0.getMz()); found=false; }
        p0=p;
      }
      if (p0!=null && p0.getMz()>rt_c && found) out.put(row, "XIC.right", p0.getMz());
    }

    return out;
  }
}
