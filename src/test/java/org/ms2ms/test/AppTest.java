package org.ms2ms.test;

import org.junit.Ignore;
import org.junit.Test;
import org.ms2ms.apps.Prep;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   8/7/15
 */
@Ignore
public class AppTest extends TestAbstract
{
  String root = "/media/data/test/mzXML/";

  @Test
  public void prepApp() throws Exception
  {
    String[] args = new String[] {"-x", "mq", "-w", "/media/data/test/mzXML/", "-r", "/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt"};
    Prep.main(args);

  }

}
