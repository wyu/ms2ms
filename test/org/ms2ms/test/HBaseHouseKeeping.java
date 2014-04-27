package org.ms2ms.test;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.util.Bytes;
import org.junit.Before;
import org.junit.Test;
import org.ms2ms.nosql.HBaseAbstract;
import org.ms2ms.nosql.HBasePeakList;

import java.io.IOException;

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
    HBaseAbstract.getAllRecord(conf, HBasePeakList.TBL_PEAKLIST);
  }
  @Test
  public void deleteTable() throws Exception
  {
    System.out.println("===========clear all record========");
    HBaseAbstract.deleteTable(conf, HBasePeakList.TBL_PEAKLIST);
  }
  @Test
  public void test1()
  {
    try
    {
      String tablename = "scores";
      String[] familys = { "grade", "course" };
      HBaseAbstract.createTable(conf, tablename, familys);

      // add record zkb
      HBaseAbstract.addRecord(conf, tablename, "zkb", "grade", "", "5");
      HBaseAbstract.addRecord(conf, tablename, "zkb", "course", "", "90");
      HBaseAbstract.addRecord(conf, tablename, "zkb", "course", "math", "97");
      HBaseAbstract.addRecord(conf, tablename, "zkb", "course", "art", "87");
      // add record baoniu
      HBaseAbstract.addRecord(conf, tablename, "baoniu", "grade", "", "4");
      HBaseAbstract.addRecord(conf, tablename, "baoniu", "course", "math", "89");

      System.out.println("===========get one record========");
      HBaseAbstract.getOneRecord(conf, tablename, "zkb");

      System.out.println("===========show all record========");
      HBaseAbstract.getAllRecord(conf, tablename);

      System.out.println("===========del one record========");
      HBaseAbstract.delRecord(conf, tablename, "baoniu");
      HBaseAbstract.getAllRecord(conf, tablename);

      System.out.println("===========show all record========");
      HBaseAbstract.getAllRecord(conf, tablename);
    } catch (Exception e)
    {
      e.printStackTrace();
    }
  }
}
