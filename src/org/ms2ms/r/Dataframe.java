package org.ms2ms.r;

import com.google.common.collect.*;
import org.ms2ms.Disposable;
import org.ms2ms.data.NameValue;
import org.ms2ms.data.collect.MapOfMultimap;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.utils.*;

import java.io.IOException;
import java.io.Writer;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/13/14
 * Time: 2:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class Dataframe implements Disposable
{
  private String                     mTitle;
  private boolean                    mKeepData = true;
  private List<String>               mRowIDs, mColIDs;
  private Map<String, Var>           mNameVar;
  private Table<String, String, Object> mData;
  private String[]                   mKeyCols;

  public Dataframe()                                         { super(); }
  public Dataframe(String s)                                 { super(); setTitle(s); }
  public Dataframe(String cvs, char delim, String... idcols) { super(); readTable(cvs, delim, idcols); setTitle(cvs); }

  public Dataframe setRowIds(List<String> s) { mRowIDs=s; return this; }
  public Dataframe setColIds(List<String> s) { mColIDs=s; return this; }

  //** factory method **//
  public static Dataframe csv(String csv, char delim, String... idcols)
  {
    Dataframe data = new Dataframe();
    data.readTable(csv, delim, idcols); data.setTitle(csv);
    return data;
  }

  //** Getters the Setters **//
  public int size() { return mData!=null?mData.rowKeySet().size():0; }
  public String       getTitle()      { return mTitle; }
//  public List<String> getColIds()     { return mColIDs; }
//  public List<String> getRowIds()     { return mRowIDs; }
  public String       getRowId(int i) { return mRowIDs!=null?mRowIDs.get(i):null; }
//  public Var          asVar(String s) { return mNameVar!=null?mNameVar.get(s):null; }
  public Var[]        asVars(String... s)
  {
    if (!hasVars(s)) return null;
    Var[] out = new Var[s.length];
    for (int i=0; i<s.length; i++) out[i] = asVar(s[i]);

    return out;
  }
  public Dataframe setTitle(String s) { mTitle=s; return this; }

  public Map<String, Object> row(int    i) { return mData!=null?mData.row(getRowId(i)):null; }
  public Map<String, Object> row(String s) { return mData!=null?mData.row(s):null; }
  public Map<String, Object> col(String s) { return mData!=null?mData.column(s):null; }
  public Var asVar(String s)
  {
    if (mNameVar==null)
    {
      mNameVar= new HashMap<String, Var>();
      if (Tools.isSet(mColIDs))
        for (String v : mColIDs) mNameVar.put(v, new Variable(v));
    }
    return mNameVar.get(s);
  }
//  public Dataframe put(String row, Var col, Object val)
//  {
//    if (mData==null) mData = HashBasedTable.create();
//    mData.put(row, col, val);
//    // update the variable cache. DON'T do it due to cost
////    addVar(col);
//
//    return this;
//  }
//  public Dataframe put(Var col, Object val)
//  {
//    if (mData  ==null) mData = HashBasedTable.create();
//    return put(mData.rowKeySet().size() + "", col, val);
//  }
  public Dataframe put(int row, String col, Object val) { return put(row+"", col, val); }
  public Dataframe put(String row, String col, Object val)
  {
    if (row!=null && col!=null && val!=null)
    {
      if (mData  ==null) mData = HashBasedTable.create();
      mData.put(row, col, val);
    }

    return this;
  }
//  public Dataframe put(String col, Object val)
//  {
//    if (mData  ==null) mData = HashBasedTable.create();
//    return put(mData.rowKeySet().size()+"", hasVars(col)?getVar(col):new Variable(col), val);
//  }
  public Dataframe addRowId(String row)
  {
    if (mRowIDs==null) mRowIDs = new ArrayList<String>();
    mRowIDs.add(row); return this;
  }
  public Dataframe addRow(String id, Map<String, Object> row)
  {
    if (mData==null) mData = HashBasedTable.create();
    for (String v : row.keySet())
      mData.put(id, v, row.get(v));

    return this;
  }
  public Var addVar(Var v)
  {
    if (mNameVar==null) mNameVar = new HashMap<String, Var>();
    if (mNameVar.put(v.toString(), v)==null)
    {
      // add to the var list if this is a new one
      if ( mColIDs==null) mColIDs = new ArrayList<String>();
      if (!mColIDs.contains(v.getName())) mColIDs.add(v.getName());
    }
    return v;
  }
  public boolean hasVars(String... vs)
  {
    if (!Tools.isSet(vs)) return false;
    for (String v : vs) if (asVar(v)==null) return false;

    return true;
  }
  public boolean hasVar(String s, boolean isCategorical)
  {
    if (Tools.isSet(s) && asVar(s)!=null && asVar(s).isCategorical()==isCategorical) return true;
    return false;
  }
  public boolean hasVar(Var s)
  {
    return (s!=null && mNameVar!=null && mNameVar.values().contains(s));
  }
//  public Object[] cells(String rowid, Var... vs)
//  {
//    Object[] lead = new Object[vs.length];
//    for (int i=0; i<vs.length; i++) lead[i]=row(rowid).get(vs[i]);
//
//    return lead;
//  }
  public Object[] cells(String rowid, String... vs)
  {
    Object[] lead = new Object[vs.length];
    for (int i=0; i<vs.length; i++) lead[i]=row(rowid).get(vs[i]);

    return lead;
  }
//  public Object get(String rowid, Var v) { return row(rowid).get(v); }
//  public Object get(String rowid, String v) { return row(rowid).get(getVar(v)); }
//  public Object cell(String rowid, Var v) { return mData!=null?mData.get(rowid,v):null; }
  public Object cell(String rowid, String s) { return mData!=null?mData.get(rowid,s):null; }
  public List<String> cols() { return mColIDs; }
  public List<String> rows() { return mRowIDs; }

  public Dataframe reorder(String... s)
  {
    if (!Tools.isSet(s) || mData==null || !Tools.isSet(mColIDs) || !Tools.isSet(mNameVar)) return this;

    mColIDs = new ArrayList<String>(s.length);
    for (String v : s) if (mNameVar.containsKey(v)) mColIDs.add(v);
    return this;
  }
  // simulate optional var so it's OK to call display(). Only the first element of the array is used
  public StringBuffer display() { return display("\t", ""); }
  public StringBuffer display(String delim, String empty)
  {
    StringBuffer buf = new StringBuffer();
    buf.append("rowid" + Strs.toString(cols(), delim) + "\n");
    for (String id : rows())
    {
      buf.append(id);
      for (String v : cols())
        buf.append(delim + (cells(id, v)!=null&& cells(id, v)[0]!=null? cells(id, v)[0]:empty));

      buf.append("\n");
    }
    return buf;
  }
  @Override
  public String toString()
  {
    return getTitle();
//    StringBuffer buf = new StringBuffer();
//    buf.append("rowid\t" + Strs.toString(getVars(), "\t") + "\n");
//    for (String id : getRowIds())
//    {
//      buf.append(id);
//      for (Var v : getVars())
//      {
//        Object t = cells(id, v)!=null? cells(id, v)[0]:null;
//        buf.append("\t" + (t!=null?(t instanceof String?(String)t:t.toString()):"--"));
//      }
//
//      buf.append("\n");
//    }
//    return buf.toString();
  }
  public void write(Writer writer, String delim)
  {
    try
    {
      writer.write(display(delim, "") + "\n");
    }
    catch (IOException io)
    {
      throw new RuntimeException("Failed to write the data frame to the output.", io);
    }
  }
  public SortedMap<Double, Double> getXY(String x, String y)
  {
    if (!hasVar(x,false) || hasVar(y,false)) return null;

    SortedMap<Double, Double> line = new TreeMap<Double, Double>();
//    Var vx=asVar(x), vy=getVar(y);
    for (String id : rows())
    {
      Tools.putNotNull(line, cell(id, x), cell(id, y));
    }
    return line;
  }
  public double[] getDoubleCol(String y)
  {
    if (!Tools.isSet(mData) || !hasVars(y)) return null;

    double[] ys = new double[rows().size()];
    for (int i=0; i<rows().size(); i++)
    {
      ys[i] = Stats.toDouble(cell(rows().get(i), y));
    }
    return ys;
  }
  public long[] getLongCol(String y)
  {
    if (!Tools.isSet(mData) || !hasVars(y)) return null;

    long[] ys = new long[rows().size()];
    for (int i=0; i<rows().size(); i++)
    {
      ys[i] = Stats.toLong(cell(rows().get(i), y));
    }
    return ys;
  }
  public long[] getLongCol(String y, Collection<String> rs)
  {
    if (!Tools.isSet(mData) || !hasVars(y)) return null;

    long[] ys = new long[rs.size()]; int i=0;
    for (String r : rs)
    {
      ys[i++] = Stats.toLong(cell(r, y));
    }
    return ys;
  }
  public String[] getStrCol(String y)
  {
    if (!Tools.isSet(mData) || !hasVars(y)) return null;

    String[] ys = new String[rows().size()];
    for (int i=0; i<rows().size(); i++)
    {
      ys[i] = cell(rows().get(i), y).toString();
    }
    return ys;
  }
  public Dataframe addVar(String v, double[] ys)
  {
    if (hasVars(v)) throw new RuntimeException("Variable " + v + " already exist!");

    if (ys!=null && ys.length==rows().size())
    {
      for (int i=0; i<ys.length; i++)
      {
        put(getRowId(i), v, ys[i]);
      }
    }
    return this;
  }
  public Dataframe renameCol(String from, String to)
  {
    if (!Tools.isSet(cols())) return this;
    // look
    int i = cols().indexOf(from);
    if (i>=0)
    {
      mColIDs.set(i, to);
      Var v = mNameVar.get(from);
      v.setName(to); mNameVar.remove(from); mNameVar.put(to, v);
      // move the actual column
      Map<String, Object> col = mData.column(from);
      if (Tools.isSet(col))
        for (String r : col.keySet())
        {
          mData.put(r, to, col.get(r));
          // remove the old column
          mData.remove(r, from);
        }
    }
    return this;
  }
  public Dataframe removeCols(String... cols)
  {
    if (Tools.isSet(cols) && mData!=null)
      for (String col : cols)
      {
        mColIDs.remove(col); mNameVar.remove(col);
        if (mData.column(col)!=null) mData.column(col).clear();
      }

    return this;
  }
  public Dataframe removeRows(String... rows)
  {
    if (Tools.isSet(rows) && mData!=null)
      for (String row : rows)
      {
        mRowIDs.remove(row);
        if (mData.row(row)!=null) mData.row(row).clear();
      }
    // re-init the columns since we removed some of the rows
    for (String col : cols()) init(asVar(col));

    return this;
  }
  //** builders **//
  public void readTable(String src, char delimiter, String... idcols)
  {
    if (!IOs.exists(src)) return;

    System.out.println("Reading the data table from " + src);
    TabFile csv=null;
    try
    {
      csv = new TabFile(src, delimiter);
      // convert the header to variables
      mColIDs = new ArrayList<>();
      mData    = HashBasedTable.create();
      for (String col : csv.getHeaders()) addVar(new Variable(col));
      // going thro the rows
      long row_counts = 0;
      while (csv.hasNext())
      {
        if (++row_counts % 10000  ==0) System.out.print(".");
//        if (  row_counts % 1000000==0) System.out.println();
        String id=null;
        if (Tools.isSet(idcols))
        {
          for (String col : idcols)
            id= Strs.extend(id, csv.get(col), "_");
        }
        else id = row_counts+"";

        addRowId(id);
        // deposit the cells
        for (String v : cols()) mData.put(id, v, csv.get(v));
      }
      csv.close();
      // setup the types
      init();
    }
    catch (IOException ioe)
    {
      throw new RuntimeException("Unable to access file: " + src, ioe);
    }
  }
  // go thro the table to determine the type of the variables. Convert them to number if necessary
  public Dataframe init()
  {
    if (!Tools.isSet(mData)) return this;
    if (mRowIDs ==null) { mRowIDs  = new ArrayList<>(mData.rowKeySet()); Collections.sort(mRowIDs); }
    if (mColIDs ==null)
    {
      mColIDs  = new ArrayList<>(mData.columnKeySet());
      mNameVar = new HashMap<>(mColIDs.size());
      for (String v : mColIDs) mNameVar.put(v, new Variable(v));
    }
    for (String v : mColIDs)
    {
      asVar(v).setFactors(null); init(asVar(v));
    }
//    System.out.println();
    return this;
  }
  public Dataframe initVars()
  {
    // re-init the columns since we removed some of the rows
    for (String col : cols()) init(asVar(col).getType(), col);
    return this;
  }
  public Dataframe init(Var v)
  {
    if (!Tools.isSet(mData)) return this;

    boolean       isNum=true;
    Set<Object> factors=new HashSet<Object>();
    for (String row : mRowIDs)
    {
      Object val = Stats.toNumber(cell(row, v.getName()));

      if (val instanceof String) isNum=false;
      if (val!=null && (!(val instanceof String) || ((String )val).length()>0)) factors.add(val);
      // put the cell back
      if (row!=null && v!=null && val!=null) mData.put(row, v.getName(), val);
    }
    if (v.isType(Var.VarType.UNKNOWN))
    {
      if (factors.size()>0 && (!isNum || factors.size()<Math.min(250, mRowIDs.size()*0.25)))
        v.setType(Var.VarType.CATEGORICAL);
      else v.setType(Var.VarType.CONTINOUOUS);
    }
    if ( v.isCategorical()) v.setFactors(factors);

    Tools.dispose(factors);
    return this;
  }
  public Dataframe init(Var.VarType type, String... vs)
  {
    if (!Tools.isSet(mData)) return this;

    for (String s : vs)
    {
      Var v = asVar(s);
      if (Tools.equals(type, Var.VarType.CATEGORICAL))
      {
        Set<Object> factors=new HashSet<>();
        for (String row : mRowIDs) factors.add(cell(row, s));
        v.setFactors(factors);
        Tools.dispose(factors);
      }
      else if (Tools.equals(type, Var.VarType.CONTINOUOUS))
      {
        for (String row : mRowIDs)
        {
          Object val = Stats.toNumber(cell(row, v.getName()));

          if (val instanceof String) break;
          // put the cell back
          if (row!=null && v!=null && val!=null) mData.put(row, v.getName(), val);
        }
      }
      v.setType(type);
    }
    return this;
  }

  /** Test whether the columns are identical
   *
   * @param A
   * @param B
   * @return
   */
  public boolean isEqualCols(String A, String B)
  {
    if (mColIDs==null || !mColIDs.contains(A) || !mColIDs.contains(B)) return false;

    MapDifference<String, Object> diff = Maps.difference(mData.column(A), mData.column(B));
    return diff.areEqual();
  }

  /** Produce a shallow copy of the self */
  @Override
  public Dataframe clone()
  {
    Dataframe out = new Dataframe();
    if (Tools.isSet(mData))
    {
      out.mData = TreeBasedTable.create(); out.mData.putAll(mData);
    }
    out.setTitle(getTitle());
    out.mKeepData = mKeepData;
    if (mRowIDs !=null) out.mRowIDs =new ArrayList<>(mRowIDs);
    if (mColIDs !=null) out.mColIDs =new ArrayList<>(mColIDs);
    if (mNameVar!=null) out.mNameVar=new HashMap<>(  mNameVar);

    return out;
  }

  /** Produce a view of the self. Any change to the Data object of the 'view' will alter the self, and vice versus
   *  Row and col IDs are free to change on their own
   *
   * @return
   */
  public Dataframe view()
  {
    Dataframe out = new Dataframe();

    out.setTitle(getTitle());
    out.mKeepData = mKeepData;
    if (mRowIDs !=null) out.mRowIDs =new ArrayList<>(mRowIDs);
    if (mColIDs !=null) out.mColIDs =new ArrayList<>(mColIDs);
    if (mNameVar!=null) out.mNameVar=new HashMap<>(  mNameVar);

    out.mData = mData;

    return out;
  }

  //********** factory methods ***************//
  // obs2 <- read.table(header=T, text='number  size type\n1   big  cat\n2 small  dog\n3 small  dog\n4   big  dog\n5   big  dog\n6   big  dog')
  public static Dataframe readtable(boolean header, String text)
  {
    Dataframe    f = new Dataframe();
    String[] lines = Strs.split(text, '\n');

    // test for the field delimiter. prefer tab if exist
    String token = lines[0].split("\t").length>1?"\t":"\\s+";
    // make the headers
    String[] headers = header&&lines.length>0?lines[0].split(token):Strs.toStringArray(Stats.newIntArray(0,lines[0].split(token).length));

    // fill out the data frame
    for (int i=(header?1:0); i<lines.length; i++)
    {
      f.addRowId(i+"");
      String[] fields = lines[i].split(token);
      for (int j=0; j<fields.length; j++)
        f.put(i+"", headers[j], fields[j]);
    }
    f.init();

    return f;
  }
  public static Dataframe readtable(String src, char delimiter, String... idcols)
  {
    Dataframe f = new Dataframe(src, delimiter, idcols);
    return f;
  }

  //********** R or Matlab style algorithms ***************//

  /** Split the data frame by the factors in variable 'v'
   *
   * @param v
   * @return
   */
  public Map<Object, Dataframe> split(String v)
  {
    if (v==null || !hasVar(v,true)) return null;

    Map<Object, Dataframe> outs = new HashMap<Object, Dataframe>();
    for (Object f : asVar(v).getFactors())
      outs.put(f, subset(v+"=="+(f instanceof String ? ("'"+f+"'"):f.toString())));

    return outs;
  }

  /** subset the rows according to the test conditions specified, similar to 'subset' function from R
   *
   * subset(animals, type=="cat")
   *
      size    type    name
   1  small   cat     lynx
   2  big     cat     tiger
   *
   * @return
   */
  public Dataframe subset(String test)
  {
    // taking
    Dataframe out = view();
    NameValue nv  = new NameValue();
    // parse the test and prepare the row subset
    List<String> rows = new ArrayList<>();
    String[]      ors = Strs.split(test, '|', true);
    if (Tools.isSet(ors))
      for (String or : ors)
      {
        String[]     ands = Strs.split(or, '&', true);
        List<String> subs = new ArrayList<>(rows());
        if (Tools.isSet(ands))
          for (String and : ands)
          {
            if (nv.parse(and, "==",">=","<=",">","<","!=","%in%") && out.hasVars(nv.name))
            {
              Iterator<String> itr = subs.iterator();
              while (itr.hasNext())
              {
                // remove the row if not meeting the test
                Object v=Stats.toNumber(cell(itr.next(), nv.name));
                if (v instanceof String)
                {
                  if ((Tools.equals(nv.token, "==") && !Tools.equals((String )v, nv.val)) ||
                      (Tools.equals(nv.token, "!=") &&  Tools.equals((String )v, nv.val))) itr.remove();
                }
                else if (v instanceof Double)
                {
                  if ((Tools.equals(nv.token, "==") && !Tools.equals((Double )v, nv.getNumber())) ||
                      (Tools.equals(nv.token, "!=") &&  Tools.equals((Double )v, nv.getNumber())) ||
                      (Tools.equals(nv.token, "<=") &&  (Double )v>nv.getNumber()) ||
                      (Tools.equals(nv.token, ">=") &&  (Double )v<nv.getNumber()) ||
                      (Tools.equals(nv.token, "<")  &&  (Double )v>=nv.getNumber()) ||
                      (Tools.equals(nv.token, ">")  &&  (Double )v<=nv.getNumber())) itr.remove();
                }
              }
            }
          }
        // combine the rows
        for (String r : subs)
          if (rows.size()==0 || !rows.contains(r)) rows.add(r);
      }

    // refresh the factors and update the rows
    return out.setRowIds(rows).initVars();
  }
  /** Partial implementation of R-aggregate
   *
   * @param by is a list of grouping categorical variables,
   * @return
   */
  public Dataframe aggregate(String... by)
  {
    Dataframe stats = new Dataframe();

    return stats;
  }
  public Dataframe melt(String... idvars)
  {
    return null;
  }

  /** generic transformation of a data frame in the style of Matlab
   *
   * pivot( [dose sbj], visit_name ) produces the following table

   []               []    'visit_name'    'visit_name'
   'dose'      'sbj'            'D0'            'D22'
   'dosed'    '1003'    [         1]    [         1]
   'dosed'    '1015'    [         1]    [         1]
   'dosed'    '1025'    [         1]    [         1]
   *
   * @param col is a categorical column whose factors will be used as the column header in the outgoing data frame
   * @param val is a numberic column whose values will be the cell in the outgoing data frame
   * @param func is the aggregate function if multiple values are found in a cell
   * @param rows are the columns that will transferred to the outgoing data frame
   * @return the outgoing data frame
   */
  public Dataframe pivot(String col, String val, Stats.Aggregator func, String... rows)
  {
    // make sure the column types are OK
    if ((Tools.isSet(rows) && !hasVars(rows)) || !hasVar(col, true) || !hasVar(val, true)) return null;
    // build the inventory
    ListMultimap<ArrayKey, Object> body = ArrayListMultimap.create();
    for (String rowid : rows())
    {
      body.put(new ArrayKey(ObjectArrays.concat(cell(rowid, col), cells(rowid, rows))), cell(rowid, val));
    }
    // construct the outgoing data frame
    Dataframe out = new Dataframe();
    for (ArrayKey keys : body.keySet())
    {
      String id = Strs.toString(Arrays.copyOfRange(keys.key, 1, keys.key.length), "");
      for (int i=1; i< keys.key.length; i++)
      {
        out.put(id, rows[i-1], keys.key[i]);
      }
      out.put(id, keys.key[0].toString(), Stats.aggregate(body.get(keys), func));
    }
    out.init(); Tools.dispose(body);
    out.reorder(ObjectArrays.concat(rows, Strs.toStringArray(asVar(col).getFactors()), String.class));

    return out;
  }
  public MultiTreeTable<Double, Double, String> index(String row, String col)
  {
    if (!hasVar(row,false) || !hasVar(col,false)) return null;

    MultiTreeTable<Double, Double, String> indice = MultiTreeTable.create();
//    Var vrow=getVar(row), vcol=getVar(col);
    for (String rowid : rows())
      indice.put(Stats.toDouble(cell(rowid, row)), Stats.toDouble(cell(rowid, col)), rowid);

    return indice;
  }
  public SortedSetMultimap<Double, String> index(String row)
  {
    if (!hasVar(row,false)) return null;

    TreeMultimap<Double, String> indice = TreeMultimap.create();
//    Var vrow=getVar(row), vcol=getVar(col);
    for (String rowid : rows())
      indice.put(Stats.toDouble(cell(rowid, row)), rowid);

    return indice;
  }

  public static MultiTreeTable<Double, Double, String>[] indice(String row, String col, Dataframe... frames)
  {
    if (!Tools.isSet(frames) || !Tools.isSet(row) || !Tools.isSet(col)) return null;

    MultiTreeTable<Double, Double, String>[] indices = new MultiTreeTable[frames.length];
    for (int i=0; i<frames.length; i++)
    {
      indices[i] = frames[i].index(row, col);
    }
    return indices;
  }

  /** ftable(animals), same as ftable(animals[,c("size","type","name")])
   *
              name   chihuahua   greatdane   lynx  tiger
   size   type
   big    cat               0          0      0     1
          dog               0          1      0     0
   small  cat               0          0      1     0
          dog               1          0      0     0
   *
   > ftable(animals[,c("size","type")])
          type cat dog
   size
   big          1   1
   small        1   1
   *
   * @param cols
   * @return
   */
  public Dataframe ftable(String... cols)
  {
    Dataframe out = new Dataframe();

    return out;
  }
  //** algorithms **//

  /**The returned data frame will contain:

   columns: all columns present in any provided data frame
   rows:    a set of rows from each provided data frame, with values in columns not present in the given data frame
            filled with missing (NA) values.

   The data type of columns will be preserved, as long as all data frames with a given column name agree on the
   data type of that column. If the data frames disagree, the column will be converted into a character strings.
   The user will need to coerce such character columns into an appropriate type.
   *
   * @param frames
   * @return
   */
  public static Dataframe smartbind(Dataframe... frames)
  {
    // prepare the merged columns
    Set<String> cols = new TreeSet<String>(); int order=0;
    for (Dataframe F : frames)
    {
      if (!Tools.isSet(F.getTitle())) F.setTitle(""+order++);
      cols.addAll(F.cols());
    }
    // the resulting dataframe
    Dataframe output = new Dataframe();
    for (Dataframe frame : frames)
    {
      for (String v : cols)
        for (String r : frame.rows())
          if (frame.cols().contains(v)) output.put(frame.getTitle()+"::"+r, v, frame.cell(r,v));
          // no value set if col didn;t exist for this dataframe. In R-routine, NA would the be the default
    }
    return output;
  }

  /** Merge two data frames by common columns or row names, or do other versions of database join operations.

   *
   * @param x and y : the data frames to be merged
//   * @param allx and ally : TRUE if the rows from x/y will be added to the output that contains no matching in the other.
   * @param all is true is the unmatched rows are to be placed in the merged data frame
   * @return the dataframe with the merge data
   */
  public static Dataframe merge(Dataframe x, Dataframe y, boolean combineSharedCol, boolean all, String... by)
  {
    if (x==null || y==null) return null;
    // cells the shared cols
    String[] shared=Strs.toStringArray(Tools.overlap(x.cols(), y.cols()));
    // set the by cols to the common if not specified
    if (!Tools.isSet(by)) by=shared;
    if (!Tools.isSet(by)) return null;
    // pool the matching rows
    MapOfMultimap<String, Integer, String> id_x_y = MapOfMultimap.create();
    for (String r : x.rows())
      id_x_y.put(Strs.toString(x.cells(r, by),"^"), 1, r);
    for (String r : y.rows())
      id_x_y.put(Strs.toString(y.cells(r, by),"^"), 2, r);

    // create the merged cols
    Table<Integer, String, String> xy_var_col = HashBasedTable.create();
    for (String v : x.cols())
      xy_var_col.put(1, v, !Tools.contains(shared, v) || Tools.contains(by, v) ? v : v + (Tools.isSet(x.getTitle())?"."+x.getTitle():".x"));
    for (String v : y.cols())
      xy_var_col.put(2, v, !Tools.contains(shared, v) || Tools.contains(by, v) ? v : v + (Tools.isSet(y.getTitle())?"."+y.getTitle():".y"));

    // create the merged results
    Dataframe out = new Dataframe();
    // deposite the columns
    for (String v : xy_var_col.values()) out.addVar(new Variable(v));

    for (String id : id_x_y.keySet())
    {
      if (id_x_y.get(id).keySet().size()>1)
        for (String xrow : id_x_y.get(id, 1))
          for (String yrow : id_x_y.get(id, 2))
          {
            String row = xrow+"."+yrow;
            out.addRowId(row);
            // make up the unique row id
            // deposit the A first
            for (String v : x.cols())
              out.put(row, xy_var_col.get(1, v), x.cell(xrow, v));
            for (String v : y.cols())
              out.put(row, xy_var_col.get(2, v), y.cell(yrow, v));
          }
      else if (all)
      {
        // singleton, not matched between x and y
        for (Integer xy : id_x_y.get(id).keySet())
        {
          for (String xyrow : id_x_y.get(id, xy))
          {
            String row = (xy==1?(xyrow+"."):("."+xyrow));
            out.addRowId(row);
            if (xy==1)
            {
              for (String v : x.cols())
                out.put(row, xy_var_col.get(xy, v), x.cell(xyrow, v));
            }
            else if (xy==2)
            {
              for (String v : y.cols())
                out.put(row, xy_var_col.get(xy, v), y.cell(xyrow, v));
            }
          }
        }

      }
    }
    // combined the shared cols bot in 'by' if asked
    if (combineSharedCol && shared!=null && shared.length>by.length)
    {
      for (String s : shared)
      {
        // skip the column if already in 'by'
        if (Tools.contains(by, s)) continue;
        if (out.isEqualCols(xy_var_col.get(1, s), xy_var_col.get(2, s)))
        {
          out.removeCols(xy_var_col.get(1, s));
          out.renameCol(xy_var_col.get(2, s), s);
        }
      }
    }
    Collections.sort(out.cols());
    return out;
  }

  @Override
  public void dispose()
  {
    mTitle=null;
    Tools.dispose(mRowIDs, mColIDs);
    Tools.dispose(mNameVar);
    Tools.dispose(mData);
  }
}

