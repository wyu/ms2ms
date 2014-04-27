package org.ms2ms.splib;

import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.peakprocessor.PeakProcessorChain;
import org.expasy.mzjava.proteomics.io.ms.spectrum.MsLibReader;
import org.expasy.mzjava.proteomics.io.ms.spectrum.msp.MspCommentParser;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.expasy.mzjava.utils.URIBuilder;
import org.ms2ms.mzjava.MspAnnotationResolver2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.UUID;
import java.util.regex.Pattern;

/** collection of algorithms and utilities related to spectral libraries
 *
 * Created by wyu on 4/20/14.
 */
public class SpLibs
{
  public static Collection<LibrarySpectrum> readMsp(File src) throws IOException
  {
    //BufferedReader reader = new BufferedReader(new FileReader("/media/data/splib/human_crp_consensus_final_true_lib.msp"));
    Collection<LibrarySpectrum> spectra = new ArrayList<LibrarySpectrum>();
    BufferedReader               reader = new BufferedReader(new FileReader(src));
    MsLibReader                     msp = new MsLibReader(reader, URIBuilder.UNDEFINED_URI,
      PeakList.Precision.FLOAT,
      new MspCommentParser(), new MspAnnotationResolver2(),
      Pattern.compile(
        "^([+-]?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?)\\s+([+-]?\\d+)\\s+\"([^\"]+)\"$"),
      new PeakProcessorChain<PepLibPeakAnnotation>());

    while (msp.hasNext())
    {
      LibrarySpectrum spec = msp.next();
      spec.setId(UUID.randomUUID());
      spectra.add(spec);
/*
      LibrarySpectrum spec = msp.next();
      HBasePeakList test = new HBasePeakList(spec);
      byte[] b = HBasePeakList.toBytes(test);
      HBasePeakList test2 = HBasePeakList.fromBytes(b);
      PeakList peak2 = test2.toPeakList();
      Peptide peptide = spec.getPeptide();
*/
    }
    return spectra;
  }

}
