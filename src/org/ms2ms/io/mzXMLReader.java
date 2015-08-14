package org.ms2ms.io;

import org.expasy.mzjava.core.io.ms.spectrum.MzxmlReader;
import org.expasy.mzjava.core.ms.peaklist.PeakAnnotation;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.EnumSet;

/** PROTOTYPE ONLY
 *
 * Created with IntelliJ IDEA.
 * User: hliu
 * Date: 8/21/14
 * Time: 3:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class mzXMLReader<A extends PeakAnnotation> extends MzxmlReader<A>
{
  // initialization if required
  static
  {

  }

  public mzXMLReader(File file, PeakList.Precision precision) throws IOException
  {
    super(file, precision);
  }
  public static <PA extends PeakAnnotation> mzXMLReader<PA> newTolerantReader(File file, PeakList.Precision precision) throws IOException {

    mzXMLReader<PA> reader = new mzXMLReader<PA>(file, precision);

    reader.removeConsistencyChecks(EnumSet.allOf(ConsistencyCheck.class));
    reader.acceptUnsortedSpectra();

    return reader;
  }

  @Override
  protected void consistencyChecks(MsnSpectrum spectrum) { }
}
