package org.ms2ms.test.ms;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.io.ms.spectrum.MgfWriter;
import org.expasy.mzjava.core.io.ms.spectrum.MzxmlReader;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.junit.Test;
import org.ms2ms.data.ms.LcMsMsDataset;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.io.MsIO;
import org.ms2ms.io.MsReaders;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.junit.Assert.assertEquals;

public class MsReaderTest extends TestAbstract
{
//  String root = "/Users/hliu/Desktop/App/2014/data/mzXML-centroid/";
  String root = "/media/data/test/mzXML/";

  // TODO: Need to modify or extend MzxmlReader to read only selected msLevel or RT range, etc
  // peak processing takes lots of time!

  @Test
  public void prepApp() throws Exception
  {

  }

  @Test
  public void serializingViaBytes() throws Exception
  {
    // grab a new scan
    MsnSpectrum[] scans  = firstScan(root+"20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML", 2, 2);
    MsnSpectrum scan=scans[0], scan2 = scans[1];

    RandomAccessFile f = new RandomAccessFile("/tmp/myobject2.data", "rw");

    // Write object out to disk
    long p1 = MsIO.write(f, MsSpectrum.adopt(scan));
    long p2 = MsIO.write(f, MsSpectrum.adopt(scan2));

    f.close();

    f = new RandomAccessFile("/tmp/myobject2.data", "r");
    MsSpectrum m1 = MsIO.read(f);
    f.close();

    f = new RandomAccessFile("/tmp/myobject2.data", "r");
    MsSpectrum m3 = MsIO.read(f, p1);
    MsSpectrum m4 = MsIO.read(f, p2);

    f.close();
  }
  @Test
  public void serializing() throws Exception
  {
    // grab a new scan
    MsnSpectrum[] scans  = firstScan(root+"20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML", 2, 2);
    MsnSpectrum scan=scans[0], scan2 = scans[1];

    FileOutputStream f_out = new FileOutputStream("/tmp/myobject.data");

    // Write object with ObjectOutputStream
    ObjectOutputStream obj_out = new ObjectOutputStream (f_out);

    // Write object out to disk
    long p1 = f_out.getChannel().position();
    obj_out.writeObject(MsSpectrum.adopt(scan));
    long p2 = f_out.getChannel().position();
    obj_out.writeObject(MsSpectrum.adopt(scan2));

    obj_out.close(); f_out.close();

    FileInputStream f_in = new FileInputStream("/tmp/myobject.data");

    ObjectInputStream obj_in = new ObjectInputStream(f_in);
    MsSpectrum m1 = (MsSpectrum )obj_in.readObject();
    obj_in.close(); f_in.close();

    f_in = new FileInputStream("/tmp/myobject.data");
    obj_in = new ObjectInputStream(f_in);

    f_in.getChannel().position(p1);
//    MsSpectrum m3 = (MsSpectrum )obj_in.readObject();
    f_in.getChannel().position(p2);
    long p3 = f_in.getChannel().position();

    MsSpectrum m4 = (MsSpectrum )obj_in.readObject();
    obj_in.close(); f_in.close();
  }
  @Test
  public void statMzXMLs() throws IOException
  {
    String[] rawfiles = {"20081129_Orbi6_NaNa_SA_FASP_blacktips_01","20081129_Orbi6_NaNa_SA_FASP_blacktips_02"};
//    Dataframe test = MsReaders.surveyMzXML(rawfiles, root, "/tmp/survey01", Range.openClosed(20d, 21d), 2);
    Dataframe test = MsReaders.surveyMzXML(rawfiles, root, "/tmp/survey01", null, 2);

    IOs.write("/tmp/surveys.txt", test.display().toString());
  }
  @Test
  public void surveyMzXMLs() throws IOException
  {
    LcMsMsDataset test = new MaxQuant("survey");
    test.setRawFilename(root+"20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML");

    test = MsReaders.surveyMzXML(test, null, 2);

    // write the examples out in MGF format
    Collection<Long> ids = test.getMzRtFileOffset().subset(495.2d, 495.35d, 0d, Double.MAX_VALUE);
    RandomAccessFile bin = test.getSpCacheFile(2);
    MgfWriter mgf = new MgfWriter(new File("/tmp/examples495_3.mgf"), PeakList.Precision.DOUBLE);
    for (Long id : ids)
    {
      bin.seek(id);
      MsnSpectrum ms = MsIO.read(bin, new MsnSpectrum());
      mgf.write(ms);
    }
    mgf.close(); bin.close();

    assertEquals(test.getMzRtFileOffset().size(), 36831);
  }
  @Test
  public void nextSpec() throws IOException
  {
    Logger.getLogger(MzxmlReader.class.getName()).setLevel(Level.SEVERE);

    File          data = new File(root+"20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML");
    MzxmlReader reader = MzxmlReader.newTolerantReader(data, PeakList.Precision.FLOAT);

    int counts=0;
    while (reader.hasNext())
    {
      MsnSpectrum spec = reader.next();
      if (++counts%100==0) System.out.print(".");
      if (counts%10000==0) System.out.println(counts);
    }
    reader.close();

    System.out.println("\n" + counts + " spectra imported");
  }
  @Test
  public void byRT()
  {
    List<MsnSpectrum> spectra = MsReaders.readMzXML(root + "20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML", Range.closed(20d, 21d), 2);

    System.out.println(spectra.size());
  }
  private MsnSpectrum[] firstScan(String file, int msLevel, int len) throws IOException
  {
    Logger.getLogger(MzxmlReader.class.getName()).setLevel(Level.SEVERE);

//    File          data = new File(root+"20081129_Orbi6_NaNa_SA_FASP_blacktips_01.mzXML");
    int order = 0;
    File        data = new File(file);
    MsnSpectrum[] spectra = new MsnSpectrum[len];
    MzxmlReader reader = MzxmlReader.newTolerantReader(data, PeakList.Precision.FLOAT);
    while (reader.hasNext())
    {
      if (order>=len) break;

      MsnSpectrum spec = reader.next();
      if (spec!=null && spec.getMsLevel()==msLevel)
      {
        spectra[order] = spec;
        order++;
      }
    }
    reader.close();

    return spectra;
  }
}
