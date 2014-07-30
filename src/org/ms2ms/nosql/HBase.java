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
abstract public class HBase
{
  static public String TBL_MSMSINDEX = "MsMsIndex";
  static public String TBL_TABLE     = "TableInfo";

  // Try to keep the ColumnFamily names as small as possible, preferably one character (e.g. "d" for data/default).
  // https://hbase.apache.org/book/rowkey.design.html
  static public byte[] FAM_FLAG      = Bytes.toBytes("f");
  static public byte[] FAM_PROP      = Bytes.toBytes("P");
  static public byte[] FAM_ID        = Bytes.toBytes("i");

  static public byte[] COL_ENDID     = Bytes.toBytes("id");
  static public byte[] COL_ENTRIES   = Bytes.toBytes("n");
  static public byte[] COL_NAME      = Bytes.toBytes("nm");
  static public byte[] COL_SOURCE    = Bytes.toBytes("s");
  static public byte[] COL_ORGANISM  = Bytes.toBytes("o");
  static public byte[] COL_VERSION   = Bytes.toBytes("v");
  static public byte[] COL_SPECTYPE  = Bytes.toBytes("t");
  static public byte[] COL_FORMAT    = Bytes.toBytes("f");

  static public String HUMAN  = "Homo Sapians";
  static public String MOUSE  = "Mouse";
  static public String ECOLI  = "EColi";
  static public String YEAST  = "Yeast";
  static public String MM     = "Msmegmatis";
  static public String WORM   = "Ce";
  static public String DM     = "Dm";
  static public String BOVINE = "Bovine";
  static public String CHICK  = "Chicken";
  static public String DRADI  = "Dradiodurans";
  static public String DROSO  = "Drosophila";
  static public String RAT    = "Rat";
  static public String MIXED  = "Mixed";
  static public String POMBE  = "S.Pombe";

  public static Configuration conf;

  // initiate the configuratio object. Only one per JVM
  static
  {
    conf = HBaseConfiguration.create();
  }
  public static HConnection getConnection() throws IOException { return HConnectionManager.createConnection(conf); }

/*
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
*/
  /*** Create a table */
  public static void createTable(String tableName, byte[]... familys) throws IOException
  {
    String[] strs = new String[familys.length];
    for (int i=0; i<familys.length; i++) strs[i] = Bytes.toString(familys[i]);
    createTable(conf, tableName, strs);
  }
  public static void createTable(Configuration conf, String tableName, String[] familys) throws IOException
  {
    HBaseAdmin admin = new HBaseAdmin(conf);
    if (admin.tableExists(tableName))
    {
      System.out.println(tableName + " table already exists!");
    } else
    {
      HTableDescriptor tableDesc = new HTableDescriptor(tableName);
      for (int i = 0; i < familys.length; i++) {
        tableDesc.addFamily(new HColumnDescriptor(familys[i]));
      }
      admin.createTable(tableDesc);
      System.out.println("create table " + tableName + " ok.");
    }
    if (admin!=null) admin.close();
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
/*
  protected static double getDouble(Result rs, String family, String col)
  {
    return getDouble(rs, Bytes.toBytes(family), Bytes.toBytes(col));
  }
  protected static int getInt(Result rs, String family, String col)
  {
    return getInt(rs, Bytes.toBytes(family), Bytes.toBytes(col));
  }
  protected static String getString(Result rs, String family, String col)
  {
    return getString(rs, Bytes.toBytes(family), Bytes.toBytes(col));
  }
*/
// Core utils that mimic the JDBC operations
//
  public static double getDouble(Result rs, byte[] family, byte[] col) { return Bytes.toDouble(rs.getValue(family, col)); }
  public static int    getInt(   Result rs, byte[] family, byte[] col) { return Bytes.toInt(   rs.getValue(family, col));  }
  public static String getString(Result rs, byte[] family, byte[] col) { return Bytes.toString(rs.getValue(family, col));  }

  public static long incre(byte[] table, byte[] rowkey, byte[] family, byte[] col, long increment) throws IOException
  {
    HConnection    conn = HBase.getConnection();
    HTableInterface tbl = conn.getTable(table);
    long next = incre(tbl, rowkey, HBase.FAM_ID, HBase.COL_ENTRIES, increment);
    tbl.close();conn.close();
    return next;
  }

  public static long incre(HTableInterface table, byte[] rowkey, byte[] family, byte[] col, long increment) throws IOException
  {
    // grab the current value
    byte[] id = select(table, rowkey, family, col);

    long next = (id==null?increment:Bytes.toLong(id)+increment);
    update(table, rowkey, family, col, Bytes.toBytes(next));

    return next;
  }
/*
  public static long incre(byte[] tbl, byte[] family, byte[] col, long increment) throws IOException
  {
    // connection to the cluster
    HConnection      conn = HConnectionManager.createConnection(conf);
    HTableInterface table = conn.getTable(tbl);

    long next = incre(table, tbl, )
    // grab the current value
    byte[] id = table.get(new Get(tbl)).getValue(family, col);

    long next = (id==null?increment:Bytes.toLong(id)+increment);
    Put row = new Put(Bytes.toBytes(tbl.toString()));
    // byte[] family, byte[] qualifier, byte[] value
    row.add(family, col, Bytes.toBytes(next));

    table.put(row);
    table.close(); conn.close(); // done with the cluster, release resources

    return next;
  }
*/
  public static void update(HTableInterface table, byte[] rowkey, byte[] family, byte[] col, byte[] val) throws IOException
  {
    Put row = new Put(rowkey);
    // byte[] family, byte[] qualifier, byte[] value
    row.add(family, col, val);
    table.put(row);
  }
  public static byte[] select(HTableInterface table, byte[] rowkey, byte[] family, byte[] col) throws IOException
  {
    return table.get(new Get(rowkey)).getValue(family, col);
  }
  // End of the core utils that mimic the JDBC operations //

/*
  public static void increSequence(HConnection conn, String table, long increment) throws IOException
  {
    HTableInterface tbl = conn.getTable(table);

    Get g = new Get(Bytes.toBytes(table));
    Result r = tbl.get(g);
    byte[] value = r.getValue(Bytes.toBytes("ID"), Bytes.toBytes("id"));

    UUID id = UUID.fromString(Bytes.toString(value)), next = id.;

  }
*/
  public static String get(Result row, byte[] family, byte[] col, String val)
  {
    return row!=null && row.getValue(family, col)!=null ? Bytes.toString(row.getValue(family, col)) : val;
  }
  public static long get(Result row, byte[] family, byte[] col, long val)
  {
    return row!=null && row.getValue(family, col)!=null ? Bytes.toLong(row.getValue(family, col)) : val;
  }
  public static char get(Result row, byte[] family, byte[] col, char val)
  {
    return row!=null && row.getValue(family, col)!=null ? Bytes.toString(row.getValue(family, col)).charAt(0) : val;
  }
  public static byte[] get(Result row, byte[] family, byte[] col, byte[] val)
  {
    return row!=null && row.getValue(family, col)!=null ? row.getValue(family, col) : val;
  }
  public static int  get(Result row, byte[] family, byte[] col, int val)
  {
    return row!=null && row.getValue(family, col)!=null ? Bytes.toInt(row.getValue(family, col)) : val;
  }

  public static boolean verifyConnection()
  {
    HConnection conn = null;
    try
    {
      try
      {
        HBaseAdmin admin = new HBaseAdmin(conf);
        HTableDescriptor tableDescriptor = new HTableDescriptor("dummy");
        tableDescriptor.addFamily(new HColumnDescriptor("personal"));
        tableDescriptor.addFamily(new HColumnDescriptor("contactinfo"));
        tableDescriptor.addFamily(new HColumnDescriptor("creditcard"));
        admin.createTable(tableDescriptor);
        admin.close();

        conn = HConnectionManager.createConnection(conf);
        if (!conn.isTableAvailable(Bytes.toBytes("dummy"))) throw new RuntimeException("Table creation failed!");
        System.out.println(conn.getRegionLocation(Bytes.toBytes("dummy"), Bytes.toBytes("test"), true).toString());
        if (conn!=null) conn.close();

        deleteTable(conf, "dummy");
      }
      catch (Exception e) { e.printStackTrace(); }
      finally
      {
        if (conn!=null && !conn.isClosed()) conn.close();
      }
    }
    catch (IOException ioe) { throw new RuntimeException(ioe); }

    return true;
  }
}
