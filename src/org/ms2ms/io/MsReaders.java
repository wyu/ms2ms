package org.ms2ms.io;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.io.ms.spectrum.MzxmlReader;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.alg.Spectra;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.r.Dataframe;
import org.ms2ms.data.HData;
import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

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
        MzxmlReader reader = MzxmlReader.newTolerantReader(new File(root+"/"+rawfile+".mzXML"), PeakList.Precision.FLOAT);
        int counts=0, inrange=0;
        while (reader.hasNext())
        {
          MsnSpectrum spec = reader.next();
          if (++counts%100==0) System.out.print(".");
          if (counts%10000==0) System.out.println(counts);

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

        System.out.print("\n" + inrange+"/"+counts + " spectra recorded/surveyed");
      }
      // close the cache files
      for (RandomAccessFile ms : ms_bin.values()) ms.close();

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

    if (data==null) data=new LcMsMsDataset();
    try
    {
      RandomAccessFile ms2 = data.getSpCacheFile(2);
      for (String rawfile : data.getRawFilenames())
      {
        MzxmlReader reader = MzxmlReader.newTolerantReader(new File(data.getmRawfileRoot()+"/"+rawfile), PeakList.Precision.FLOAT);
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

        System.out.println("\n" + data.getMzRtFileOffset().size() + " spectra imported");
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
}
