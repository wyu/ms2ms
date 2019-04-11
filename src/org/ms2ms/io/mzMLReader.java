package org.ms2ms.io;

import com.google.common.collect.Range;
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
import java.util.Collection;

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
  public static Dataframe readScanInfo(Dataframe data, double ppm, double dRT, String root, String pattern, Collection<Double> neuloss, double... channels) throws IOException
  {
    String[] searches = IOs.listFiles(root, pattern);

    System.out.println("mzML files: "+ searches.length);

    if (Tools.isSet(searches))
      for (String s : searches)
      {
//        data = mzMLReader.inferPrecursorsFromMS2(data, s).init(true);
        data = readScanInfo(data, ppm, s, true, neuloss, channels).init(true);
        data = mzMLReader.readPeptideFeatures(data, ppm, dRT, s).init(true);
      }

    return data;
  }
  public static Dataframe readScanInfo(Dataframe out, double ppm, String filename, boolean loadIons, Collection<Double> neuloss, double... channels) throws IOException
  {
    System.out.println("Reading scans from "+filename+"...");

    boolean needIons = (loadIons || Tools.isSet(neuloss) || Tools.isSet(channels));

    File file = new File(filename); String run = file.getName().substring(0,file.getName().indexOf('.'));
    MzMLUnmarshaller mzml = new MzMLUnmarshaller(file, false, null);

    // looping through the scans
    MzMLObjectIterator<uk.ac.ebi.jmzml.model.mzml.Spectrum> spectrumIterator = mzml.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", uk.ac.ebi.jmzml.model.mzml.Spectrum.class);

    if (out==null) out = new Dataframe(filename);
    int rows=0;
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next(), !needIons);
      String row = filename+"#"+ms.getScanNumbers().getFirst().getValue();
      out.put(row, "Scan",       ms.getScanNumbers().getFirst().getValue());
      out.put(row, "MsLevel",    ms.getMsLevel());
      out.put(row, "FragMethod", ms.getFragMethod());
      out.put(row, "RT",         ms.getRetentionTimes().getFirst().getTime()/60);
      out.put(row, "mz",         ms.getPrecursor().getMz());
      out.put(row, "z",          ms.getPrecursor().getCharge());
      out.put(row, "Intensity",  ms.getPrecursor().getIntensity());
      out.put(row, "Ions",       ms.size());
      out.put(row, "TIC",        ms.getTotalIonCurrent());
      out.put(row, "BaseMz",     ms.getBasePeakMz());
      out.put(row, "BaseInt",    ms.getBasePeakIntensity());
      out.put(row, "RunFile",    filename);
      out.put(row, "Run",        run);
//      if (ms.getScanNumbers().getFirst().getValue()==7135)
//        System.out.println();
      if (ms!=null && ms.getMsLevel()>1 && ms.size()>0)
      {
        if (Tools.isSet(channels))
          for (double channel : channels)
          {
            double delta = channel*1E-6*ppm, sum=Spectra.sum(ms, Range.closed(channel-delta,channel+delta));
            if (sum>0)
              out.put(row, "LM"+Tools.d2s(channel,4), sum);
          }
        if (Tools.isSet(neuloss))
          for (double nl : neuloss)
            for (double iso=0d; iso<4d; iso=iso+1d)
            {
              // gather the precursor and its neutral loss
              double prec     = ms.getPrecursor().getMz()+(iso/ms.getPrecursor().getCharge()),
                     targeted = prec-(nl/ms.getPrecursor().getCharge()),
                     delta    = targeted*1E-6*ppm,
                     sumP     = Spectra.sum(ms, Range.closed(prec-delta,    prec+delta)),
                     sum      = Spectra.sum(ms, Range.closed(targeted-delta,targeted+delta));
              if (sum >0) out.put(row, "NL"+Tools.d2s(nl,4)+"iso"+((int )iso), sum);
              if (sumP>0) out.put(row, "Prec"+((int )iso), sumP);
            }
      }
      if (++rows%1000==0) System.out.print(".");
      // remove the spectrum from the memory
      ms.clear(); ms=null;
    }
    System.out.println(rows);

    return out;
  }
  public static Dataframe readPeptideFeatures(Dataframe out, double ppm, double dRT, String filename) throws IOException
  {
    System.out.println("Reading peptide features from "+filename+"...");

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
