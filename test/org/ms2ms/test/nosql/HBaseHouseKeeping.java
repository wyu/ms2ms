package org.ms2ms.test.nosql;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.nosql.HBase;
import org.ms2ms.nosql.ms.HBasePeakList;
import org.ms2ms.test.TestAbstract;

public class HBaseHouseKeeping extends TestAbstract
{
  private Configuration conf = null;

  @Before
  public void setUp()
  {
    conf = HBaseConfiguration.create();
  }
  @Test
  public void listTable()
  {
    System.out.println("===========show all record========");
    HBase.getAllRecord(conf, HBasePeakList.TBL_PEAKLIST);
  }
  @Test
  public void deleteTable() throws Exception
  {
    System.out.println("===========clear all record========");
    HBase.deleteTable(conf, HBasePeakList.TBL_PEAKLIST);
  }
  @Test
  public void test1()
  {
    try
    {
      String tablename = "scores";
      String[] familys = { "grade", "course" };
      HBase.createTable(conf, tablename, familys);

      // add record zkb
      HBase.addRecord(conf, tablename, "zkb", "grade", "", "5");
      HBase.addRecord(conf, tablename, "zkb", "course", "", "90");
      HBase.addRecord(conf, tablename, "zkb", "course", "math", "97");
      HBase.addRecord(conf, tablename, "zkb", "course", "art", "87");
      // add record baoniu
      HBase.addRecord(conf, tablename, "baoniu", "grade", "", "4");
      HBase.addRecord(conf, tablename, "baoniu", "course", "math", "89");

      System.out.println("===========cells one record========");
      HBase.getOneRecord(conf, tablename, "zkb");

      System.out.println("===========show all record========");
      HBase.getAllRecord(conf, tablename);

      System.out.println("===========del one record========");
      HBase.delRecord(conf, tablename, "baoniu");
      HBase.getAllRecord(conf, tablename);

      System.out.println("===========show all record========");
      HBase.getAllRecord(conf, tablename);
    } catch (Exception e)
    {
      e.printStackTrace();
    }
  }
}
