package org.ms2ms.io;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
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
  public static Dataframe readScanInfo(Dataframe data, String root, String pattern) throws IOException
  {
    String[] searches = IOs.listFiles(root, pattern);

    if (Tools.isSet(searches))
      for (String s : searches)
        data = readScanInfo(data, s);

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
//
//  public Spectrum nextSpectrum()
//  {
//    return null;
//  }
//  private void parse(XMLStreamReader parser) throws IOException, XMLStreamException
//  {
// 		while(parser.hasNext())
//    {
// 			if(parser.isStartElement())
//      {
// 				if(parser.getLocalName().equalsIgnoreCase(SCAN))
//        { //root node
// 					if(parser.getAttributeCount() > 0)
//          {
// 						Integer scanNr = null;
// 						Integer msLevel = null;
// 						Double retentionTime = null;
// 						for(int a = 0; a < parser.getAttributeCount(); a++)
//            {
// 							if(XMLs.isA(parser, a, NUMBER))
//              {
// 								scanNr = XMLs.getInt(parser, a);
// 							}
//              else if (XMLs.isA(parser, a, RT)) {
// 								retentionTime = parseRetentionTime(parser.getAttributeValue(a));
// 							} else if (parser.getAttributeLocalName(a).equalsIgnoreCase("msLevel")) {
// 								msLevel = Integer.parseInt(parser.getAttributeValue(a));
// 							}
// 						}
//
// 						if(msLevel != null && scanNr != null && retentionTime != null && msLevel == 1) {
//// 							scanNrToTime.put(scanNr, retentionTime);
//// 							timeToScanNr.put(retentionTime,scanNr);
// 						}
// 					}
// 				}
// 			}
// 			parser.next();
// 		}
// 	}
//  /**
//  	 * Some manual string parsing due to the fact of some
//  	 * unexpected characters in the retentiontime field.
//  	 *
//  	 * + Conversion of retention time to minutes!!
//  	 *
//  	 * @param rt
//  	 * @return
//  	 */
//  	private Double parseRetentionTime(String rt) {
//  		Double result = null;
//
//  		boolean seconds = false;
//
//  		if(rt.startsWith("PT")) {
//  			rt = rt.substring(2);
//  		}
//  		if(rt.contains("S")) {
//  			rt = rt.replace("S", "");
//  			seconds = true; //Guess that the S stands for seconds
//  		}
//
//  		result = Double.parseDouble(rt);
//
//  		if(seconds) {//Guess that the S stands for seconds
//  			result = result / 60;
//  		}
//
//  		return result;
//  	}

}
