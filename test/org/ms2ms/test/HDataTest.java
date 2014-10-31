package org.ms2ms.test;

import org.junit.Before;
import org.junit.Test;
import org.ms2ms.r.Dataframe;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/26/14
 * Time: 11:33 AM
 * To change this template use File | Settings | File Templates.
 */
public class HDataTest extends TestAbstract
{
  String root = "/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/";
  Dataframe evidences = null;

  @Before
  public void setUp()
  {
    evidences = new Dataframe(root+"evidence1k.txt", '\t');
  }

  @Test
  public void dataframe2hdata() throws Exception
  {

  }
}
