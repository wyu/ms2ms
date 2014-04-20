package org.ms2ms.nosql;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.hbase.*;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

/** the base class to support I/O routine against the HBASE
 *
 * Created by wyu on 4/16/14.
 */
abstract public class HbaseAbstract
{
  protected static Configuration conf;

  // initiate the configuratio object. Only one per JVM
  static
  {
    conf = HBaseConfiguration.create();
  }
  public static void execute(String tbl, Callable func) throws IOException
  {
    // connection to the cluster
    HConnection conn = HConnectionManager.createConnection(conf);

    // When the cluster connection is established get an HTableInterface for each operation or thread.
    // HConnection.getTable(...) is lightweight. The table is really just a convenient place to call
    // table method and for a temporary batch cache.
    // It is in fact less overhead than HTablePool had when retrieving a cached HTable.
    // The HTableInterface returned is not thread safe as before.
    // It's fine to get 1000's of these.
    // Don't cache the longer than the lifetime of the HConnection
    HTableInterface table = conn.getTable(tbl);

    // TODO do something with the table

    // just flushes outstanding commit, no futher cleanup needed, can be omitted.
    // HConnection holds no references to the returned HTable objects, they can be GC'd as soon as they leave scope.
    table.close();

    conn.close(); // done with the cluster, release resources
  }
  /*** Create a table */
  public static void creatTable(Configuration conf, String tableName, String[] familys) throws Exception
  {
    HBaseAdmin admin = new HBaseAdmin(conf);
    if (admin.tableExists(tableName))
    {
      System.out.println("table already exists!");
    } else
    {
      HTableDescriptor tableDesc = new HTableDescriptor(tableName);
      for (int i = 0; i < familys.length; i++) {
        tableDesc.addFamily(new HColumnDescriptor(familys[i]));
      }
      admin.createTable(tableDesc);
      System.out.println("create table " + tableName + " ok.");
    }
  }

  /*** Delete a table */
  public static void deleteTable(Configuration conf, String tableName) throws Exception
  {
    try
    {
      HBaseAdmin admin = new HBaseAdmin(conf);
      admin.disableTable(tableName);
      admin.deleteTable(tableName);
      System.out.println("delete table " + tableName + " ok.");
    } catch (MasterNotRunningException e)
    {
      e.printStackTrace();
    } catch (ZooKeeperConnectionException e)
    {
      e.printStackTrace();
    }
  }

  /*** Put (or insert) a row */
  public static void addRecord(Configuration conf, String tableName, String rowKey, String family, String qualifier, String value) throws Exception
  {
    try
    {
      HTable table = new HTable(conf, tableName);
      Put put = new Put(Bytes.toBytes(rowKey));
      put.add(Bytes.toBytes(family), Bytes.toBytes(qualifier), Bytes
        .toBytes(value));
      table.put(put);
      System.out.println("insert recored " + rowKey + " to table "
        + tableName + " ok.");
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  /*** Delete a row */
  public static void delRecord(Configuration conf, String tableName, String rowKey) throws IOException
  {
    HTable table = new HTable(conf, tableName);
    List<Delete> list = new ArrayList<Delete>();
    Delete del = new Delete(rowKey.getBytes());
    list.add(del);
    table.delete(list);
    System.out.println("del recored " + rowKey + " ok.");
  }

  /*** Get a row */
  public static void getOneRecord(Configuration conf, String tableName, String rowKey) throws IOException
  {
    HTable table = new HTable(conf, tableName);
    Get get = new Get(rowKey.getBytes());
    Result rs = table.get(get);
    for(KeyValue kv : rs.raw()){
      System.out.print(new String(kv.getRow()) + " " );
      System.out.print(new String(kv.getFamily()) + ":" );
      System.out.print(new String(kv.getQualifier()) + " " );
      System.out.print(kv.getTimestamp() + " " );
      System.out.println(new String(kv.getValue()));
    }
  }
  /*** Scan (or list) a table */
  public static void getAllRecord (Configuration conf, String tableName)
  {
    try{
      HTable table = new HTable(conf, tableName);
      Scan s = new Scan();
      ResultScanner ss = table.getScanner(s);
      for(Result r:ss){
        for(KeyValue kv : r.raw()){
          System.out.print(new String(kv.getRow()) + " ");
          System.out.print(new String(kv.getFamily()) + ":");
          System.out.print(new String(kv.getQualifier()) + " ");
          System.out.print(kv.getTimestamp() + " ");
          System.out.println(new String(kv.getValue()));
        }
      }
    } catch (IOException e){
      e.printStackTrace();
    }
  }
}
