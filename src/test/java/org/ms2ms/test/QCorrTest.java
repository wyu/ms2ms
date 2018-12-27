package org.ms2ms.test;

import org.junit.Test;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;

import java.io.BufferedReader;
import java.io.FileReader;

public class QCorrTest extends TestAbstract
{
  @Test
  public void pk() throws Exception
  {
    String filename = "/Users/kfvf960/Apps/contrib/MsApps/Yates/160921_48Protein_TMT3plex_test.ms2_CCArrays.txt";
    BufferedReader CCA = new BufferedReader(new FileReader(filename));
    while (CCA.ready())
    {
//      S       1138    3
//      1543.7655       0
//      1543.75549      3324254.5
//      1543.74548      0

      String[] items = Strs.split(CCA.readLine(), '\t');
    }
  }
}
