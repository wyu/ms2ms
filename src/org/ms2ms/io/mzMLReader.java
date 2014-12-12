package org.ms2ms.io;

import org.expasy.mzjava.core.ms.spectrum.Spectrum;
import org.ms2ms.data.XMLs;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/19/14
 * Time: 5:10 PM
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

  public Spectrum nextSpectrum()
  {
    return null;
  }
  private void parse(XMLStreamReader parser) throws IOException, XMLStreamException
  {
 		while(parser.hasNext())
    {
 			if(parser.isStartElement())
      {
 				if(parser.getLocalName().equalsIgnoreCase(SCAN))
        { //root node
 					if(parser.getAttributeCount() > 0)
          {
 						Integer scanNr = null;
 						Integer msLevel = null;
 						Double retentionTime = null;
 						for(int a = 0; a < parser.getAttributeCount(); a++)
            {
 							if(XMLs.isA(parser, a, NUMBER))
              {
 								scanNr = XMLs.getInt(parser, a);
 							}
              else if (XMLs.isA(parser, a, RT)) {
 								retentionTime = parseRetentionTime(parser.getAttributeValue(a));
 							} else if (parser.getAttributeLocalName(a).equalsIgnoreCase("msLevel")) {
 								msLevel = Integer.parseInt(parser.getAttributeValue(a));
 							}
 						}

 						if(msLevel != null && scanNr != null && retentionTime != null && msLevel == 1) {
// 							scanNrToTime.put(scanNr, retentionTime);
// 							timeToScanNr.put(retentionTime,scanNr);
 						}
 					}
 				}
 			}
 			parser.next();
 		}
 	}
  /**
  	 * Some manual string parsing due to the fact of some
  	 * unexpected characters in the retentiontime field.
  	 *
  	 * + Conversion of retention time to minutes!!
  	 *
  	 * @param rt
  	 * @return
  	 */
  	private Double parseRetentionTime(String rt) {
  		Double result = null;

  		boolean seconds = false;

  		if(rt.startsWith("PT")) {
  			rt = rt.substring(2);
  		}
  		if(rt.contains("S")) {
  			rt = rt.replace("S", "");
  			seconds = true; //Guess that the S stands for seconds
  		}

  		result = Double.parseDouble(rt);

  		if(seconds) {//Guess that the S stands for seconds
  			result = result / 60;
  		}

  		return result;
  	}
}
