package org.ms2ms.io;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamReader;
import java.io.FileInputStream;

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

}
