package org.ms2ms.test.io;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.junit.Test;
import org.ms2ms.algo.PSMs;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;

import java.util.List;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   2/10/15
 */
public class IOsTest extends TestAbstract
{
  @Test
  public void StrsSplit() throws Exception
  {
    String p1 = "IETLMRNLM[15.9949]PWRK",
           p2 = "+229.163HMK+229.163K+229.163HAK+229.163K+229.163MK+229.163K+229.163QMK+229.163K+229.163";
    String r1 = "(^.*?\\[|\\]\\s*$)", r2 = "\\]\\s*,\\s*\\[";

    PeptideMatch m1 = PSMs.fromNumModSequence(p1),
                 m2 = PSMs.fromNumModSequence(p2);
    List<String> item01 = Strs.splits(p1, "[+-.\\d]+"),
    item02 = Strs.splits(p1, r1),
    item03 = Strs.splits(p1, r2);

    System.out.println(item01.size());
  }
  @Test
  public void recursiveListing()
  {
    IOs.listFiles("/tmp", new WildcardFileFilter("*.tmp"));
  }
}
