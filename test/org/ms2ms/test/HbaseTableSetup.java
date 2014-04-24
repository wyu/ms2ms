package org.ms2ms.test;

import org.junit.Test;
import org.ms2ms.nosql.HBaseProteomics;

import java.io.IOException;

/**
 * Created by wyu on 4/21/14.
 */
public class HbaseTableSetup extends TestAbstract
{
  @Test
  public void createTables() throws IOException
  {
    HBaseProteomics.ensurePeakListTable();
  }
}
