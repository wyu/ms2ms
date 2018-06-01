package org.ms2ms.io;

import com.google.common.collect.TreeMultimap;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.SortedMap;

/**
 * Deprecated! Use jmzml reader instead
 *
 * To change this template use File | Settings | File Templates.
 */
public class mzMLReader extends mzReader
{
  static
  {
    SCAN      = "scan";
    NUMBER    = "num";
    RT        = "retentionTime";
  }
  public static Dataframe readScanInfo(Dataframe data, double ppm, double dRT, String root, String pattern) throws IOException
  {
    String[] searches = IOs.listFiles(root, pattern);

    if (Tools.isSet(searches))
      for (String s : searches)
      {
//        data = mzMLReader.inferPrecursorsFromMS2(data, s).init(true);
        data = readScanInfo(data, s).init(true);
        data = mzMLReader.readPeptideFeatures(data, ppm, dRT, s).init(true);
      }

    return data;
  }
  public static Dataframe readScanInfo(Dataframe out, String filename) throws IOException
  {
    System.out.println("Reading "+filename+"...");

    File file = new File(filename); String run = file.getName().substring(0,file.getName().indexOf('.'));
    MzMLUnmarshaller mzml = new MzMLUnmarshaller(file, false, null);

    // looping through the scans
    MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> spectrumIterator = mzml.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);

    if (out==null) out = new Dataframe(filename);
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next(), true);
      String row = filename+"#"+ms.getScanNumbers().getFirst().getValue();
      out.put(row, "Scan",       ms.getScanNumbers().getFirst().getValue());
      out.put(row, "MsLevel",    ms.getMsLevel());
      out.put(row, "FragMethod", ms.getFragMethod());
      out.put(row, "RT",         ms.getRetentionTimes().getFirst().getTime()/60);
      out.put(row, "mz",         ms.getPrecursor().getMz());
      out.put(row, "z",          ms.getPrecursor().getCharge());
      out.put(row, "Intensity",  ms.getPrecursor().getIntensity());
      out.put(row, "Ions", ms.size());
      out.put(row, "TIC", ms.getTotalIonCurrent());
      out.put(row, "RunFile", filename);
      out.put(row, "Run", run);
    }
    return out;
  }
  public static Dataframe readPeptideFeatures(Dataframe out, double ppm, double dRT, String filename) throws IOException
  {
    System.out.println("Reading "+filename+"...");

    MultiTreeTable<Double, Double, String> rt_mz_row = MultiTreeTable.create();
    MultiTreeTable<Double, Double, Peak>   rt_mz_ms1 = MultiTreeTable.create(), rt_mz_mz1 = MultiTreeTable.create();

    for (String row : out.rows())
      if (out.cell(row,"RunFile").equals(filename) && out.getInteger(row, "MsLevel")>1)
        rt_mz_row.put(out.getDouble(row, "RT"), out.getDouble(row, "mz"), row);

    File file = new File(filename); String run = file.getName().substring(0,file.getName().indexOf('.'));
    MzMLUnmarshaller mzml = new MzMLUnmarshaller(file, false, null);

    // looping through the scans
    MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> spectrumIterator = mzml.unmarshalCollectionFromXpath(
        "/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);

    if (out==null) out = new Dataframe(filename);
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next(), false);
      if (ms.getMsLevel()==1)
      {
        double rt = ms.getRetentionTimes().getFirst().getTime()/60d, mz=ms.getPrecursor().getMz();
        // deposit the ms1 intensity if fell within the bound of a peptide precursor
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
    }
    // calculate the LC peak centroids
    for (String row : out.rows())
      if (out.cell(row,"RunFile").equals(filename) && out.getInteger(row, "MsLevel")>1)
      {
        Collection<Peak> XIC = rt_mz_ms1.get(out.getDouble(row, "RT"), out.getDouble(row, "mz")),
                      XIC_mz = rt_mz_mz1.get(out.getDouble(row, "RT"), out.getDouble(row, "mz"));
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
            else if (p==null && p0!=null && p.getMz()>rt_c) { out.put(row, "XIC.right", p0.getMz()); found=false; }
            p0=p;
          }
          if (p0!=null && p0.getMz()>rt_c && found) out.put(row, "XIC.right", p0.getMz());
        }
      }

    return out;
  }
  public static Dataframe inferPrecursorsFromMS2(Dataframe out, String filename) throws IOException
  {
    System.out.println("Reading "+filename+"...");

    File file = new File(filename); String run = file.getName().substring(0,file.getName().indexOf('.'));
    MzMLUnmarshaller mzml = new MzMLUnmarshaller(file, false, null);

    // looping through the scans
    MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> spectrumIterator = mzml.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);

    if (out==null) out = new Dataframe(filename);
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next(), false);
      if (ms.getMsLevel()==2)
      {
        Spectra.precursorByComplements(ms, null);

      }
      String row = filename+"#"+ms.getScanNumbers().getFirst().getValue();
      out.put(row, "Scan",       ms.getScanNumbers().getFirst().getValue());
    }
    return out;
  }
}
