package org.ms2ms.test.nosql;

import org.junit.Test;
import org.ms2ms.nosql.ms.HBaseProteomics;
import org.ms2ms.test.TestAbstract;

import java.io.IOException;

/**
 * Created by wyu on 4/21/14.
 */
public class HbaseTableSetup extends TestAbstract
{
  @Test
  public void createTables() throws IOException
  {
    HBaseProteomics.ensureTables();
  }
}
