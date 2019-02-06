package org.ms2ms.nosql.ms;

import org.apache.hadoop.hbase.HBaseConfiguration;
import org.apache.hadoop.hbase.client.*;
import org.apache.hadoop.hbase.util.Bytes;
import org.ms2ms.nosql.HBases;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/** Representation of spectral library in HBases
 *  User: wyu
 *  Date: 5/14/14
 */
public class HBaseSpLib implements Serializable
{
  static public byte[] TBL_SPLIB = Bytes.toBytes("SpLib");
  static public String LIB_MSP   = "msp";
  static public String LIB_SPTXT = "sptxt";

  private String name, path, source, organism, version, format;
  private long entries;
//  private int id;
  private byte[] spectype;

  public HBaseSpLib(String n, String s, String o, String v, String f, long e, byte[] stype)
  {
    spectype = stype;
    name     = n;
    source   = s;
    organism = o;
    version  = v;
    format   = f;
    entries  = e;
  }

  public HBaseSpLib(Result row)
  {
    if (row==null) return;
    // process the row
    spectype = HBases.get(row, HBases.FAM_PROP, HBases.COL_SPECTYPE, spectype);
    name     = HBases.get(row, HBases.FAM_PROP, HBases.COL_NAME, name);
    source   = HBases.get(row, HBases.FAM_PROP, HBases.COL_SOURCE, source);
    organism = HBases.get(row, HBases.FAM_PROP, HBases.COL_ORGANISM, organism);
    version  = HBases.get(row, HBases.FAM_PROP, HBases.COL_VERSION, version);
    format   = HBases.get(row, HBases.FAM_PROP, HBases.COL_FORMAT, format);
    entries  = HBases.get(row, HBases.FAM_PROP, HBases.COL_ENTRIES, entries);
//    id       = HBases.cells(row, HBases.FAM_ID,   HBases.COL_ENDID,   id);
  }
  public String getName()     { return name; }
  public byte[] getSpecType() { return spectype; }

  public boolean isFormat(String s) { return format!=null && format.equals(s); }

  @Override
  public String toString()
  {
    return Bytes.toString(spectype) + ": '" + name + "' (" + organism + ") by '" + source + "' #" + version + ": entries=" + entries;
  }
//  public static long increID(String tbl, long block) throws IOException
//  {
//    return incre(Bytes.toBytes(tbl), HBases.FAM_ID, HBases.COL_ENDID, block);
//  }
//  public static long increID(byte[] lib, long block) throws IOException
//  {
//    return HBases.incre(TBL_SPLIB, lib, HBases.FAM_ID, HBases.COL_ENDID, block);
//  }
  public static long increEntries(byte[] lib, long incre) throws IOException
  {
    return HBases.incre(TBL_SPLIB, lib, HBases.FAM_PROP, HBases.COL_ENTRIES, incre);
  }

  public static HBaseSpLib set(HTableInterface tbl, String source, String format,
                         String organism, String version, byte[] spectype, String name) throws IOException
  {
    // only need an integer
//    int id = (int )increID(tbl.getTableName(), 1l);
    Put row = new Put(Bytes.toBytes(name));
    // byte[] family, byte[] qualifier, byte[] value
//    row.add(HBases.FAM_ID,   HBases.COL_ENDID,   Bytes.toBytes(id));
    row.add(HBases.FAM_PROP, HBases.COL_NAME,    Bytes.toBytes(name));
    row.add(HBases.FAM_PROP, HBases.COL_SOURCE,  Bytes.toBytes(source));
    row.add(HBases.FAM_PROP, HBases.COL_ORGANISM,Bytes.toBytes(organism));
    row.add(HBases.FAM_PROP, HBases.COL_VERSION, Bytes.toBytes(version));
    row.add(HBases.FAM_PROP, HBases.COL_SPECTYPE,spectype);
    row.add(HBases.FAM_PROP, HBases.COL_FORMAT,  Bytes.toBytes(format));

    tbl.put(row);
    return new HBaseSpLib(name, source, organism, version, format, 0, spectype);

  }
  public static Collection<HBaseSpLib> getAll()
  {
    Collection<HBaseSpLib> libs = new ArrayList<HBaseSpLib>();
    try
    {
      HBaseProteomics.ensureTables();

      HConnection      conn = HConnectionManager.createConnection(HBaseConfiguration.create());
      HTableInterface table = conn.getTable(TBL_SPLIB);
      Scan             scan = new Scan();

      scan.setCaching(1); scan.setBatch(1);
      scan.addColumn(HBases.FAM_PROP, HBases.COL_NAME);

      int rows = 0;
      ResultScanner resultScanner = table.getScanner(scan);
      Iterator<Result> iterator = resultScanner.iterator();
      while (iterator.hasNext())
      {
        Result result = iterator.next();
        //if you want to cells the entire row
        Get get = new Get(result.getRow());
        Result entireRow = table.get(get);
        // process the row
        libs.add(new HBaseSpLib(entireRow));
      }
    }
    catch (IOException ie)
    {
      ie.printStackTrace();
    }
    return libs;
  }
}
