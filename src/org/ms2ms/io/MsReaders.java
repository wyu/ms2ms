package org.ms2ms.io;

import com.compomics.util.io.FilenameExtensionFilter;
import com.google.common.collect.Range;
import com.google.common.collect.RowSortedTable;
import com.google.common.collect.TreeBasedTable;
import org.expasy.mzjava.core.io.ms.spectrum.MgfWriter;
import org.expasy.mzjava.core.io.ms.spectrum.MzxmlReader;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.*;
import org.expasy.mzjava.proteomics.io.ms.ident.MzIdentMlReader;
import org.expasy.mzjava.proteomics.io.ms.ident.mzidentml.v110.AbstractParamType;
import org.expasy.mzjava.proteomics.io.ms.ident.mzidentml.v110.CVParamType;
import org.expasy.mzjava.proteomics.io.ms.ident.mzidentml.v110.UserParamType;
import org.ms2ms.algo.LCMSMS;
import org.ms2ms.algo.Peaks;
import org.ms2ms.algo.PurgingPeakProcessor;
import org.ms2ms.algo.Spectra;
import org.ms2ms.data.collect.ImmutableNavigableMap;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.data.HData;
import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.Precursor;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

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
  public static Dataframe surveyMzXML(String root, String cache, Range<Double> rt, int... mslevel)
  {
    String[] rawfiles = new File(root).list(new FilenameExtensionFilter("mzXML"));

    return surveyMzXML(rawfiles, root, cache, rt, mslevel);
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
        mzXMLReader reader = mzXMLReader.newTolerantReader(new File(root+"/"+rawfile), PeakList.Precision.FLOAT);
        int counts=0, inrange=0;
        TreeMap<Integer, Integer> scan_tin = new TreeMap<>();
        while (reader.hasNext())
        {
          MsnSpectrum spec = reader.next();
          if (++counts%1000==0) System.out.print(".");
          if (counts%100000==0) System.out.println(counts);

          int scan = spec.getScanNumbers().getFirst().getValue();

//          if (scan==5648)
//          {
//            System.out.println();
//          }
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
            if (spec.getMsLevel()==3)
            {
              //System.out.println();
              int tin_ms2 = scan_tin.get(spec.getParentScanNumber().getValue());
              stats.put(tin_ms2,"TIC.ms3", Tools.d2s(spec.getTotalIonCurrent(), 1));
              stats.put(tin_ms2,"Scan.ms3",scan);
              stats.put(tin_ms2,"m/z.ms3",Tools.d2s(spec.getPrecursor().getMz(), 4));
              stats.put(tin_ms2,"Charge.ms3",spec.getPrecursor().getCharge());
              stats.put(tin_ms2,"RT.ms3",spec.getRetentionTimes().getFirst().getTime()/60d);
              stats.put(tin_ms2,"MS.ms3",spec.getMsLevel());
            }
            else
            {
              stats.put(tinrange,"TIC", Tools.d2s(spec.getTotalIonCurrent(), 1));
              stats.put(tinrange,"Raw file",rawfile);
              stats.put(tinrange,"Scan",spec.getScanNumbers().getFirst().getValue());
//            stats.put(tinrange,"FilePointer",ms.getFilePointer());
              stats.put(tinrange,"FilePointer",MsIO.write(ms, MsSpectrum.adopt(spec)));
              stats.put(tinrange,"m/z",Tools.d2s(spec.getPrecursor().getMz(), 4));
              stats.put(tinrange,"Charge",spec.getPrecursor().getCharge());
              stats.put(tinrange,"RT",spec.getRetentionTimes().getFirst().getTime()/60d);
              stats.put(tinrange,"MS",spec.getMsLevel());

              Map<String, Double> ss = Spectra.survey(spec, Range.closed(0d, 250d), 5, 10.0);
              for (String tag : ss.keySet()) stats.put(tinrange,tag,ss.get(tag));

              scan_tin.put(spec.getScanNumbers().getFirst().getValue(), tinrange);
            }

            inrange++; tinrange++;
//            MsIO.write(ms, spec);
          }
          else { spec.clear(); spec=null; }
        }
        reader.close();

        System.out.print(inrange + "/" + counts + " spectra recorded/surveyed" + "\n");
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

    return stats.init(true);
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
      spec.setPrecursor(new Peak(LCMSMS.parseNominalPrecursorMz(spec.getComment()), 0));
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

    // readSpectrumIdentifier the ions
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
  public static <P extends AbstractParamType> Map<String, P> toCVMap(List<P> paramGroup)
  {
    if (paramGroup == null || paramGroup.isEmpty()) return null;

    Map<String, P> paramMap = new HashMap<>();
    for (P paramType : paramGroup)
      if (paramType instanceof CVParamType)
      {
        paramMap.put(((CVParamType )paramType).getAccession(), paramType);
      }
      else if (paramType instanceof UserParamType)
      {
        paramMap.put("user:"+paramType.getName(), paramType);
      }

    return paramMap;
  }
  public static ScanNumberList newScanNumberList(Map<String, AbstractParamType> cvParamTypes)
  {
    AbstractParamType cvParamType = cvParamTypes.get(MzIdentMlReader.SCAN_NUMBER_CV);
    if (cvParamType == null) return new ScanNumberList();

    try
    {
      int scanNumber = Integer.parseInt(cvParamType.getValue());
      return new ScanNumberList(scanNumber);
    } catch (NumberFormatException e) {

      throw new UnsupportedOperationException("Cannot yet parse scan number like " + cvParamType.getValue(), e);
    }
  }
  public static RetentionTimeList newRetentionTimeList(Map<String, AbstractParamType> cvParamTypes)
  {
    AbstractParamType cvParamType = cvParamTypes.get(MzIdentMlReader.RETENTION_TIME_CV);
    if (cvParamType == null) return new RetentionTimeList();

    String unit = cvParamType.getUnitName();
    TimeUnit timeUnit;
    if ("second".equals(unit)) {

      timeUnit = TimeUnit.SECOND;
    } else {

      throw new UnsupportedOperationException("Cannot parse units = " + unit);
    }

    double time;
    try {
      time = Double.parseDouble(cvParamType.getValue());
    } catch (NumberFormatException e) {

      throw new UnsupportedOperationException("Time values of " + cvParamType.getValue() + " are not yet suported", e);
    }

    RetentionTimeList retentionTimes = new RetentionTimeList();
    retentionTimes.add(time, timeUnit);

    return retentionTimes;
  }
  public static Dataframe splitMs2MS3(String mzml_root) throws Exception
  {
    File  xmlFile = new File(mzml_root+".mzML");
    String    run = xmlFile.getName().substring(0, xmlFile.getName().lastIndexOf("."));
    MgfWriter ms2 = new MgfWriter(new File(mzml_root+".ms2.mgf"), PeakList.Precision.FLOAT),
              ms3 = new MgfWriter(new File(mzml_root+".ms3.mgf"), PeakList.Precision.FLOAT),
             ms3y = new MgfWriter(new File(mzml_root+".ms3y.mgf"),PeakList.Precision.FLOAT),
             ms3b = new MgfWriter(new File(mzml_root+".ms3b.mgf"),PeakList.Precision.FLOAT);

    MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(xmlFile);

    // looping through the scans
    RowSortedTable<Double, Integer, MsnSpectrum> cache = TreeBasedTable.create();
    MzMLObjectIterator<Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);

    System.out.println("Splitting the MS2 and MS3 contents from:" +mzml_root);
    int counts=0;
    Dataframe xref = new Dataframe("MS2-MS3 xref: " + mzml_root);
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next());
      if (++counts % 1000 == 0) System.out.print(".");
      if      (ms.getMsLevel()==2)
      {
        ms2.write(ms);
        cache.put(LCMSMS.parseNominalPrecursorMz(ms.getComment()), ms.getScanNumbers().getFirst().getValue(), ms);
        String scan=ms.getScanNumbers().getFirst().getValue()+"", row=run+"#"+scan;
        xref.put(row,"Run", run);
        xref.put(row,"Scan", scan);
        xref.put(row,"RT2", ms.getRetentionTimes().getFirst().getTime());
        xref.put(row,"mz2", ms.getPrecursor().getMz());
        xref.put(row,"z2",  ms.getPrecursor().getCharge());
        xref.put(row,"MsLevel", 2);
      }
      else if (ms.getMsLevel()==3)
      {
        // figure out the parent MS2 scan within the isolation window
        Integer parent=null, pz=null; Double pmz=null;
        for (Double mz : cache.rowKeySet().subSet(ms.getPrecursor().getMz()-0.02, ms.getPrecursor().getMz()+0.02))
          if (parent==null || Collections.max(cache.row(mz).keySet())>parent)
          {
            parent = Collections.max(cache.row(mz).keySet());
            pmz    = cache.get(mz, parent).getPrecursor().getMz();
            pz     = cache.get(mz, parent).getPrecursor().getCharge();
          }
        if (parent!=null)
        {
          String row=run+"#"+parent;
          xref.put(row,"MS3", ms.getScanNumbers().getFirst().getValue());
          xref.put(row,"RT3", ms.getRetentionTimes().getFirst().getTime());
          xref.put(row,"mz3", ms.getPrecursor().getMz());
          xref.put(row,"z3",  ms.getPrecursor().getCharge());

          ms.setParentScanNumber(ms.getScanNumbers().getFirst());
          ms.getScanNumbers().clear();
          ms.addScanNumber(parent);
          // copy the charge state from the parent
          if (pmz!=null) ms.getPrecursor().setMzAndCharge(pmz, pz!=null?pz:ms.getPrecursor().getCharge());
        }
        else
        {
          System.out.println();
        }

        // remove the reporter ions to avoid search interference
        Spectra.notchUpto(ms, 150d);

        ms3.write(ms.copy(new PurgingPeakProcessor<>()));
        // change the procursor to the SPS MASS
        try
        {
          pmz = Double.parseDouble(Strs.split(Strs.split(ms.getComment(), '@', true)[1], ' ')[1]);
          // educated guess of the charge state
          pz  = (pmz>ms.getPrecursor().getMz()?ms.getPrecursor().getCharge()-1:ms.getPrecursor().getCharge());

          // assume this is a y-ion
          ms.getPrecursor().setMzAndCharge(pmz, pz);
          ms3y.write(ms);
          // or a b-ion
          ms.getPrecursor().setMzAndCharge(pmz-28, pz);
          ms3b.write(ms);

        } catch (Exception e) {}
      }
    }
    System.out.println();

    ms2.close(); ms3.close();

    FileWriter xw = new FileWriter(mzml_root+".xref");
    xref.init(true).write(xw, "\t");
    xw.close();

    return xref;
  }
  public static MsnSpectrum newQualSpec(String datafile) throws IOException
  {
    MsnSpectrum ms2 = new MsnSpectrum();

//    m/z     Intensity       Relative        Charge  Noise
//    120.08127       10492.3           2.37            0.00         2681.32
//    121.08474        5012.8           1.13            0.00         2695.37
    TabFile tab = new TabFile(datafile, TabFile.tabb);
    while (tab.hasNext())
    {
      double mz=tab.getDouble("m/z"), z=tab.getDouble("Charge");
      ms2.add((z>1?(mz*z-(z-1)*1.00783d):mz), tab.getDouble("Intensity")/tab.getDouble("Noise"));
    }
    return ms2;
  }

  public static void mzML2MGF(String mzml, double left, double right, String mgf, int... scans) throws IOException
  {
    System.out.println("Saving the scans to " + mgf);

//    File mzmlFile = new File("/Users/yuw/Apps/pipeline/Chorus/data/151107_Ecoli_12protein_Spike_MS2_R2_Fr6.mzML");
//    int[] scans = new int[] {20647, 21623, 22456};
    MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(new File(mzml));
    // looping through the scans
    MzMLObjectIterator<Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);

    System.out.println("Going through the MS2 scans from:"+mzml);
    int counts=0; MsnSpectrum ms1=null;
    Collection<MsnSpectrum> msms = new ArrayList<>();
    MgfWriter MGF = new MgfWriter(new File(mgf), PeakList.Precision.DOUBLE);
    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next());
      if (++counts % 1000 == 0) System.out.print(".");
      if      (ms.getMsLevel()==1 && ms1==null) ms1=ms;
      else if (ms.getMsLevel()==1 && ms1!=null)
      {
        writeMGF(MGF, ms1,ms, left, right, msms, scans);
        ms1=ms; msms.clear();
      }
      else if (ms.getMsLevel()==2) msms.add(ms);
    }
    if (Tools.isSet(msms) && ms1!=null) writeMGF(MGF, ms1,null, left, right, msms, scans);

    MGF.close();
  }
  public static SortedMap<Integer, MsnSpectrum> read(String ms1)
  {
    SortedMap<Integer, MsnSpectrum> spectra = new TreeMap<>();
    RandomAccessFile bin = null;

    try
    {
      try
      {
        bin = new RandomAccessFile(ms1, "r");
        while (1==1)
        {
          byte[] bs = new byte[bin.readInt()]; bin.read(bs);
          MsnSpectrum spec = MsSpectrum.fromBytes(bs).toMsnSpectrum(); bs=null;
          spectra.put(spec.getScanNumbers().getFirst().getValue(), spec);
        }
      }
      finally {
        if (bin!=null) bin.close();
      }
    }
    catch (IOException e)
    {

    }
    return spectra;
  }
  // read the MS1 scan and save them as the binary objects
  public static void mzML2MS1(String mzml, String peaks) throws IOException
  {
    System.out.println("Saving the scans to " + peaks);

    MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(new File(mzml));
    // looping through the scans
    MzMLObjectIterator<Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);

    System.out.println("Going through the MS1 scans from:"+mzml);
    int counts=0;
    RandomAccessFile bin = new RandomAccessFile(peaks, "rw"), scans = new RandomAccessFile(peaks+".scans", "rw");
    Map<Integer, Long> scan_fseek = new HashMap<>();

    while (spectrumIterator.hasNext())
    {
      MsnSpectrum ms = MsReaders.from(spectrumIterator.next());
      if (++counts % 1000 == 0) System.out.print(".");
      if (ms.getMsLevel()==1)
        scan_fseek.put(ms.getScanNumbers().getFirst().getValue(), MsIO.write(bin, MsSpectrum.adopt(ms)));
    }
    IOs.writeIntLongMap(scans, scan_fseek);
    bin.close(); scans.close();
  }

  private static void writeMGF(MgfWriter MGF, MsnSpectrum ms10, MsnSpectrum ms11,
                               double left, double right, Collection<MsnSpectrum> spectra, int... scans) throws IOException
  {
    for (MsnSpectrum ms : spectra)
    for (int scan : scans)
      if (scan==-1 || scan==ms.getScanNumbers().getFirst().getValue())
      {
        // parse the isolated m/z. It maybe different from the precursor m/z
        String[] strs = ms.getComment().split("@")[0].split(" ");
        double center = (Tools.isSet(strs)? Stats.toDouble(strs[strs.length-1]):ms.getPrecursor().getMz());

        // save the precursor isolation region
        String isolation = "isolation/"+ms10.getScanNumbers().getFirst().getValue()+"/"+center+
                       (ms11!=null?("/"+ms11.getScanNumbers().getFirst().getValue()):""),
               line = Peaks2Str(Peaks.isolate(ms10, center+left, center+right)) +
                   (ms11!=null?(","+Peaks2Str(Peaks.isolate(ms11, center+left, center+right))):"");

        if (Strs.isSet(line)) ms.setComment(isolation+","+line);
        MGF.write(ms);
      }
  }
  public static String Peaks2Str(List<Peak> isolated)
  {
    String line = null;
    if (Tools.isSet(isolated))
      for (int i=0; i<isolated.size(); i++)
        line = Strs.extend(line, isolated.get(i).getMz()+"/"+isolated.get(i).getIntensity()+"/"+isolated.get(i).getCharge(), ";");

    return line;
  }
}
