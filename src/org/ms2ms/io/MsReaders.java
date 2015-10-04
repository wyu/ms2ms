package org.ms2ms.io;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.io.ms.spectrum.MzxmlReader;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.expasy.mzjava.core.ms.spectrum.RetentionTimeDiscrete;
import org.expasy.mzjava.core.ms.spectrum.ScanNumberDiscrete;
import org.expasy.mzjava.core.ms.spectrum.TimeUnit;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.data.HData;
import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.DataFormatException;
import java.util.zip.Inflater;

/**
 * User: wyu
 * Date: 8/22/14

 wyu@yubuntu:~/Projects/contrib/mzJava/bitbucket/mzjava/mzjava-core/mzjava-core/src/test/resources$ ls -l org/expasy/mzjava/core/io/ms/spectrum/
 total 3244
 -rw-r--r-- 1 wyu wyu     897 Sep 27 06:29 emptyscan.mzXML
 -rw-r--r-- 1 wyu wyu  949457 Sep 27 06:29 issue30.mzXML
 -rw-r--r-- 1 wyu wyu    2802 Sep 27 06:29 mgf_test2.mgf
 -rw-r--r-- 1 wyu wyu   38368 Sep 27 06:29 mgf_test.mgf
 -rw-r--r-- 1 wyu wyu 2205194 Sep 27 06:29 mzxml_test_2.mzXML
 -rw-r--r-- 1 wyu wyu   27744 Sep 27 06:29 mzxml_test-inconsistent.mzXML
 -rw-r--r-- 1 wyu wyu   27757 Sep 27 06:29 mzxml_test.mzXML
 -rw-r--r-- 1 wyu wyu   27744 Sep 27 06:29 mzxml-unsorted_test.mzXML
 -rw-r--r-- 1 wyu wyu   27744 Sep 27 06:29 xml_event_utils.mzXML
 */
public class MsReaders
{
  /** Initialize a HData where MS/MS scans are the rows.
   *
   * @param rootdir
   * @return
   * @throws IOException
   */
  public static HData readMaxquantScans(String rootdir) throws IOException
  {
    HData data = new HData();
    Dataframe scans = Dataframe.readtable(rootdir+"/msmsScans.txt", '\t', "Raw file","Scan number"),
              annot = Dataframe.readtable(rootdir+"/msms.txt", '\t', "Raw file","Scan number");


    return data;
  }
  public static List<MsnSpectrum> readMzXML(String data, Range<Double> rt, int... mslevel)
  {
    Logger.getLogger(MzxmlReader.class.getName()).setLevel(Level.SEVERE);

    List<MsnSpectrum> spectra = null;
    try
    {
      MzxmlReader    reader = MzxmlReader.newTolerantReader(new File(data), PeakList.Precision.FLOAT);
      int counts=0; spectra = new ArrayList<MsnSpectrum>();
      while (reader.hasNext())
      {
        MsnSpectrum spec = reader.next();
        if (++counts%100==0 && spectra.size()==0) System.out.print(".");
        if (counts%10000==0)                      System.out.println(counts);

        // move on before we reach the lower bound of the RT range if requested
        if (rt!=null && Spectra.before(spec.getRetentionTimes(), rt.lowerEndpoint())) continue;
        // quit if we move pass the limit
        if (rt!=null && Spectra.after( spec.getRetentionTimes(), rt.upperEndpoint())) break;
        if (Tools.contains(mslevel, spec.getMsLevel()))
        {
          spec.trimToSize();
          spectra.add(spec);
          if (spectra.size()%100==0) System.out.print("$");
        }
        else { spec.clear(); spec=null; }
      }
      reader.close();

      System.out.println("\n" + spectra.size() + " spectra imported");
    }
    catch (IOException ie)
    {
      throw new RuntimeException("Not able to open data file: " + data, ie);
    }
    return spectra;
  }

  /** Survey the mzXML file and return a table of m/z, rt and file pointers
   *
   * @param rt
   * @param mslevel
   * @return
   */
  public static Dataframe surveyMzXML(String[] rawfiles, String root, String cache, Range<Double> rt, int... mslevel)
  {
    Logger.getLogger(MzxmlReader.class.getName()).setLevel(Level.SEVERE);

    Dataframe stats=new Dataframe("stats");
    Map<Integer, RandomAccessFile> ms_bin = new HashMap<>();
    try
    {
      int tinrange=0;
      for (String rawfile : rawfiles)
      {
        System.out.println("surveying " + rawfile);

//        MzxmlReader reader = MzxmlReader.newTolerantReader(new File(root+"/"+rawfile+".mzXML"), PeakList.Precision.FLOAT);
        mzXMLReader reader = mzXMLReader.newTolerantReader(new File(root+"/"+rawfile+".mzXML"), PeakList.Precision.FLOAT);
        int counts=0, inrange=0;
        while (reader.hasNext())
        {
          MsnSpectrum spec = reader.next();
          if (++counts%1000==0) System.out.print(".");
          if (counts%100000==0) System.out.println(counts);

          // move on before we reach the lower bound of the RT range if requested
          if (rt!=null && Spectra.before(spec.getRetentionTimes(), rt.lowerEndpoint())) continue;
          // quit if we move pass the limit
          if (rt!=null && Spectra.after( spec.getRetentionTimes(), rt.upperEndpoint())) break;
          if (!Tools.isSet(mslevel) || Tools.contains(mslevel, spec.getMsLevel()))
          {
            // write the spectrum to the file
            if (!ms_bin.containsKey(spec.getMsLevel())) ms_bin.put(spec.getMsLevel(), new RandomAccessFile(root+"/"+cache+".ms"+spec.getMsLevel(), "rw"));

            RandomAccessFile ms = ms_bin.get(spec.getMsLevel());
            // populate the stats
            stats.put(tinrange,"TIC", Tools.d2s(spec.getTotalIonCurrent(), 1));
            stats.put(tinrange,"Raw file",rawfile);
            stats.put(tinrange,"Scan number",spec.getScanNumbers().getFirst().getValue());
//            stats.put(tinrange,"FilePointer",ms.getFilePointer());
            stats.put(tinrange,"FilePointer",MsIO.write(ms, MsSpectrum.adopt(spec)));
            stats.put(tinrange,"m/z",Tools.d2s(spec.getPrecursor().getMz(), 4));
            stats.put(tinrange,"Charge",spec.getPrecursor().getCharge());
            stats.put(tinrange,"RT",spec.getRetentionTimes().toString());
            inrange++; tinrange++;
//            MsIO.write(ms, spec);
          }
          else { spec.clear(); spec=null; }
        }
        reader.close();

        System.out.print(inrange+"/"+counts + " spectra recorded/surveyed" + "\n");
      }
      // close the cache files
      for (RandomAccessFile ms : ms_bin.values())
      {
        ms.close();
      }

      System.out.println("\n" + tinrange + " spectra recorded");
    }
    catch (IOException ie)
    {
      throw new RuntimeException("Not able to open data file", ie);
    }

    return stats.init();
  }
  public static LcMsMsDataset surveyMzXML(LcMsMsDataset data, Range<Double> rt, int... mslevel)
  {
    Logger.getLogger(MzxmlReader.class.getName()).setLevel(Level.SEVERE);

//    if (data==null) data=new LcMsMsDataset();
    try
    {
      RandomAccessFile ms2 = data.getSpCacheFile(2);
      for (String rawfile : data.getRawFilenames())
      {
        MzxmlReader reader = MzxmlReader.newTolerantReader(new File(data.getRawfileRoot()+"/"+rawfile), PeakList.Precision.FLOAT);
        int counts=0;
        while (reader.hasNext())
        {
          MsnSpectrum spec = reader.next();
          if (++counts%100==0) System.out.print(".");
          if (counts%10000==0) System.out.println(counts);

          // move on before we reach the lower bound of the RT range if requested
          if (rt!=null && Spectra.before(spec.getRetentionTimes(), rt.lowerEndpoint())) continue;
          // quit if we move pass the limit
          if (rt!=null && Spectra.after( spec.getRetentionTimes(), rt.upperEndpoint())) break;
          if (Tools.contains(mslevel, spec.getMsLevel()))
          {
            data.add(spec, ms2.getFilePointer());
            // write the spectrum to the file
            MsIO.write(ms2, spec);
          }
          else { spec.clear(); spec=null; }
        }
        reader.close();

//        System.out.println("\n" + data.getMzRtFileOffset().size() + " spectra imported");
      }
      ms2.close();
    }
    catch (IOException ie)
    {
      throw new RuntimeException("Not able to open data file", ie);
    }
    return data;
  }
  public static MsnSpectrum nextSpectrum(MzxmlReader reader) throws IOException
  {
    Logger.getLogger(MzxmlReader.class.getName()).setLevel(Level.SEVERE);
    return reader.hasNext()?reader.next():null;
  }
  public static MsnSpectrum from(Spectrum ms)
  {
    if (ms==null) return null;

    MsnSpectrum spec = new MsnSpectrum();

    // populate the new spectrum with the old
    spec.setMsLevel(MsIO.getInt(ms.getCvParam(), "MS:1000511"));
    spec.setSpectrumIndex(ms.getIndex());
    // set the scan filter:
    // 'ITMS + c NSI r d Full ms2 838.7747@cid35.00 [226.0000-2000.0000]'
    spec.setComment(MsIO.get(ms.getScanList().getScan().get(0).getCvParam(), "MS:1000512"));
    if (spec.getMsLevel()==3)
    {
      // 'FTMS + p NSI sps d Full ms3 838.77@cid35.00 432.89@hcd40.00 [120.00-1800.00]'
      String[] items = Strs.split(spec.getComment().substring(spec.getComment().indexOf("ms3"), spec.getComment().indexOf('[')), ' ', true);
      // grab the MS2 info and change the precursor to that of the MS2
      spec.setPrecursor(new Peak(Stats.toDouble(Strs.split(items[1], '@', true)[0]), 0));
      // set the activation mode
      spec.setFragMethod(Strs.split(items[2], '@', true)[1].substring(0, 3));
      // no spectrumref avail for MS3
    }
    else if (spec.getMsLevel()>1)
    {
      Precursor p = ms.getPrecursorList().getPrecursor().get(0);
      if (Strs.isSet(p.getSpectrumRef())) spec.setParentScanNumber(new ScanNumberDiscrete(MsIO.ScanNumberFromSpectrumRef(p.getSpectrumRef())));
      spec.setFragMethod(MsIO.hasAccession(p.getActivation().getCvParam(), "MS:1000133") ? "cid" : "");
      Double ai = MsIO.getDouble(p.getSelectedIonList().getSelectedIon().get(0).getCvParam(), "MS:1000042");
      Integer z = MsIO.getInt(p.getSelectedIonList().getSelectedIon().get(0).getCvParam(), "MS:1000041");
      spec.setPrecursor(new Peak(MsIO.getDouble(p.getSelectedIonList().getSelectedIon().get(0).getCvParam(), "MS:1000744"),
                                 ai!=null?ai:0, z!=null?z:0));
    }

    // set the scan number
    spec.addScanNumber(MsIO.ScanNumberFromSpectrumRef(ms.getId()));
    spec.addRetentionTime(new RetentionTimeDiscrete(MsIO.getDouble(ms.getScanList().getScan().get(0).getCvParam(), "MS:1000016"), TimeUnit.MINUTE));

    // read the ions
    Number[] mzs=null, ais=null;
    for (BinaryDataArray bin : ms.getBinaryDataArrayList().getBinaryDataArray())
    {
      int     precision = MsIO.hasAccession(bin.getCvParam(), "MS:1000521")?32:64;
      String  compressionType=MsIO.hasAccession(bin.getCvParam(), "MS:1000574")?"zlib":"none";
      try
      {
        if (MsIO.hasAccession(bin.getCvParam(), "MS:1000514"))
        {
          // decode the m/z string
          mzs = bin.getBinaryDataAsNumberArray();
        }
        else if (MsIO.hasAccession(bin.getCvParam(), "MS:1000515"))
        {
          // decode the m/z string
          ais = bin.getBinaryDataAsNumberArray();
        }
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }
    if (mzs!=null && ais!=null && mzs.length==ais.length)
      for (int i=0; i<mzs.length; i++)
      {
        spec.add(mzs[i].doubleValue(), ais[i].doubleValue());
      }

    return spec;
  }
  /*******************
   * Function returns an array of doubles that was decoded from the passed string
   * https://github.com/dfermin/lucXor/blob/master/src/lucxor/mzMLreader.java
   * @param peaks
   * @return
   * @throws IOException
   * @throws DataFormatException
   */
  private static double[] decode_string(String peaks, int precision, boolean isBase64Encoded, String compressionType) throws IOException, DataFormatException
  {
    byte[] decoded;

    if(isBase64Encoded) {
      decoded = org.apache.commons.codec.binary.Base64.decodeBase64(peaks);
    } else {
      decoded = peaks.getBytes();
    }

    if(compressionType.equals("zlib")) { // need to decompress the bytes
      decoded = zlibUncompressBuffer(decoded);
    }

    int decodedLen = decoded.length; // in bytes
    int chunkSize = precision / 8;
    int numPeaks = decodedLen / chunkSize;
    double[] retAry = new double[numPeaks];

    if (precision == 32) {

      int asInt;
      float asFloat = 0.0f;
      int offset;

      for (int i = 0; i < numPeaks; i++) {
        offset = i * chunkSize;

        asInt = (  (decoded[offset + 0] & 0xFF)) // zero shift
            | ((decoded[offset + 1] & 0xFF) << 8)
            | ((decoded[offset + 2] & 0xFF) << 16)
            | ((decoded[offset + 3] & 0xFF) << 24);
        asFloat = Float.intBitsToFloat(asInt);
        retAry[i] = asFloat;
      }
    } else if (precision == 64) {
      long asLong;
      double asDouble = 0.0d;
      int offset;

      for (int i = 0; i < numPeaks; i++) {
        offset = i * chunkSize;

        asLong = ( (long) (decoded[offset + 0] & 0xFF)) // zero shift
            | ((long) (decoded[offset + 1] & 0xFF) << 8)
            | ((long) (decoded[offset + 2] & 0xFF) << 16)
            | ((long) (decoded[offset + 3] & 0xFF) << 24)
            | ((long) (decoded[offset + 4] & 0xFF) << 32)
            | ((long) (decoded[offset + 5] & 0xFF) << 40)
            | ((long) (decoded[offset + 6] & 0xFF) << 48)
            | ((long) (decoded[offset + 7] & 0xFF) << 56);
        asDouble = Double.longBitsToDouble(asLong);

        retAry[i] = asDouble;
      }
    } else {
      throw new IllegalArgumentException("Precision can only be 32 or 64 bits.");
    }

    return retAry;
  }
  /** https://github.com/dfermin/lucXor/blob/master/src/lucxor/mzMLreader.java
   * Inflates zLib compressed byte[].
   * @param compressed zLib compressed bytes
   * @return inflated byte array
   * @throws IOException should never happen, ByteArrayOutputStream is in-memory
   * @throws DataFormatException in case of malformed input byte array
   */
  public static byte[] zlibUncompressBuffer(byte[] compressed) throws IOException, DataFormatException
  {

    Inflater decompressor = new Inflater();
    decompressor.setInput(compressed);

    ByteArrayOutputStream bos = new ByteArrayOutputStream(compressed.length);
    byte[] buf = new byte[decompressor.getRemaining() * 2];
    try {
      // Decompress the data
      while (decompressor.getRemaining() > 0) {
        int count = decompressor.inflate(buf);
        bos.write(buf, 0, count);
      }

    } finally {
      try {
        bos.close();
      } catch (IOException nope) {
        // This exception doesn't matter, but it totally should not happen
        throw nope;
      }
    }
    decompressor.end();
    byte[] result = bos.toByteArray();
    return result;
  }
}
