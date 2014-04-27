package org.ms2ms.test;

import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.junit.Test;
import org.ms2ms.nosql.HBaseProteomics;
import org.ms2ms.splib.SpLibs;

import java.io.*;
import java.util.Collection;

/** Reading the content of splib
 *
 * Created by wyu on 4/13/14.
 */
public class SpLibsTest extends TestAbstract
{
  @Test
  public void readMsp() throws IOException
  {
    Collection<LibrarySpectrum> spectra = SpLibs.readMsp(new File("/media/data/splib/human_crp_consensus_final_true_lib.msp"));
    // save the spectrum and indice to HBase
    HBaseProteomics.index(spectra, 50d, 450d, 7, 4d);
    HBaseProteomics.listTables();

    //assert spectra.size()==92;
  }
}
