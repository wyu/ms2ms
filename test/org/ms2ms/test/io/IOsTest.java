package org.ms2ms.test.io;

import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.junit.Test;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.IOs;

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
  public void recursiveListing()
  {
    IOs.listFiles("/tmp", new WildcardFileFilter("*.tmp"));
  }
}
