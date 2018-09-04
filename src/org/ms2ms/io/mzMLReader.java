package org.ms2ms.io;

import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

import java.io.File;
import java.io.IOException;

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
        XIC(ms, ppm, dRT, rt_mz_row,rt_mz_ms1,rt_mz_mz1);
    }
    // calculate the LC peak centroids
    for (String row : out.rows())
      if (out.cell(row,"RunFile").equals(filename) && out.getInteger(row, "MsLevel")>1)
        out = detectXIC(out, row, "mz", "RT", rt_mz_ms1, rt_mz_mz1);

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
