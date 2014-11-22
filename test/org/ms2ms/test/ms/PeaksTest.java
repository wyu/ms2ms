package org.ms2ms.test.ms;

import hep.aida.ref.Test;
import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.io.MsIO;
import org.ms2ms.test.TestAbstract;

import java.util.List;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   11/21/14
 */
public class PeaksTest extends TestAbstract
{
  List<MsnSpectrum> spectra = null;

  public void setUp()
  {
    spectra = MsIO.readSpectra("data/examples_16__AAELIANSLATAGDGLIELR__z2.ms2");
  }

  @Test
  public void baseline() throws Exception
  {
     assert(1==1);
  }

}
