package org.ms2ms.data;

import com.google.common.collect.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;
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
public class Dataframe
{
  private boolean                    mKeepData = true;
  private List<String>               mRowIDs;
  private List<Var>                  mColVars;
  private Map<String, Var>           mNameVar;
  private Table<String, Var, Object> mData;

  public Dataframe()                                         { super(); }
  public Dataframe(String cvs, char delim, String... idcols) { super(); readTable(cvs, delim, idcols);}

  //** Getters the Settings **//
  public int size() { return mData!=null?mData.rowKeySet().size():0; }
  public List<Var>    getVars()   { return mColVars; }
  public List<String> getRowIds() { return mRowIDs; }
  public String       getRowId(int i) { return mRowIDs!=null?mRowIDs.get(i):null; }
  public Var[]        toVars(String... s)
  {
    if (!hasVars(s)) return null;
    Var[] out = new Var[s.length];
    for (int i=0; i<s.length; i++) out[i] = getVar(s[i]);

    return out;
  }

  public Map<Var,    Object> row(int    i) { return mData!=null?mData.row(getRowIds().get(i)):null; }
  public Map<Var,    Object> row(String s) { return mData!=null?mData.row(s):null; }
  public Map<String, Object> col(Var    s) { return mData!=null?mData.column(s):null; }
  public Var getVar(String s)
  {
    if (mNameVar==null)
    {
      mNameVar= new HashMap<>();
      for (Var v : mColVars) mNameVar.put(v.toString(), v);
    }
    return mNameVar.get(s);
  }
  public Dataframe put(String row, Var col, Object val)
  {
    if (mData==null) mData = HashBasedTable.create();
    mData.put(row, col, val);
    // update the variable cache
    addVar(col);

    return this;
  }
  public Dataframe put(Var col, Object val)
  {
    if (mData  ==null) mData = HashBasedTable.create();
    return put(mData.rowKeySet().size() + "", col, val);
  }
  public Dataframe put(String row, String col, Object val)
  {
    return put(row, hasVars(col)?getVar(col):new Variable(col), val);
  }
  public Dataframe put(String col, Object val)
  {
    if (mData  ==null) mData = HashBasedTable.create();
    return put(mData.rowKeySet().size()+"", hasVars(col)?getVar(col):new Variable(col), val);
  }
  public Dataframe addRowId(String row)
  {
    if (mRowIDs==null) mRowIDs = new ArrayList<>();
    mRowIDs.add(row); return this;
  }
  public Dataframe addRow(String id, Map<Var, Object> row)
  {
    if (mData==null) mData = HashBasedTable.create();
    for (Var v : row.keySet())
      mData.put(id, v, row.get(v));

    return this;
  }
  public Var addVar(Var v)
  {
    if (mNameVar==null) mNameVar = new HashMap<>();
    if (mNameVar.put(v.toString(), v)==null)
    {
      // add to the var list if this is a new one
      if (mColVars==null) mColVars = new ArrayList<>();
      mColVars.add(v);
    }
    return v;
  }
  public boolean hasVars(String... vs)
  {
    if (!Tools.isSet(vs)) return false;
    for (String v : vs) if (getVar(v)==null) return false;

    return true;
  }
  public boolean hasVar(String s, boolean isCategorical)
  {
    if (Tools.isSet(s) && getVar(s)!=null && getVar(s).isCategorical()==isCategorical) return true;
    return false;
  }
  public boolean hasVar(Var s)
  {
    return (s!=null && mColVars!=null && mColVars.contains(s));
  }
  public Object[] get(String rowid, Var... vs)
  {
    Object[] lead = new Object[vs.length];
    for (int i=0; i<vs.length; i++) lead[i]=row(rowid).get(vs[i]);

    return lead;
  }
  public Object cell(String rowid, Var v) { return mData!=null?mData.get(rowid,v):null; }

  public Dataframe reorder(String... s)
  {
    if (!Tools.isSet(s) || mData==null || !Tools.isSet(mColVars) || !Tools.isSet(mNameVar)) return this;

    mColVars = new ArrayList<>(s.length);
    for (String v : s) if (mNameVar.containsKey(v)) mColVars.add(mNameVar.get(v));
    return this;
  }
  // simulate optional var so it's OK to call display(). Only the first element of the array is used
  public StringBuffer display(StringBuffer... bufs)
  {
    StringBuffer buf = (Tools.isSet(bufs) && bufs[0]!=null) ? bufs[0] : new StringBuffer();
    buf.append("rowid\t" + Strs.toString(getVars(), "\t") + "\n");
    for (String id : getRowIds())
    {
      buf.append(id);
      for (Var v : getVars())
        buf.append("\t" + (get(id, v)!=null&&get(id, v)[0]!=null?get(id, v)[0]:"--"));

      buf.append("\n");
    }
    return buf;
  }
  public void write(Writer writer, String delim)
  {
    try
    {
      writer.write("rowid" + delim + Strs.toString(getVars(), delim) + "\n");
      for (String id : getRowIds())
      {
        writer.write(id);
        for (Var v : getVars())
          writer.write(delim + (get(id, v)!=null&&get(id, v)[0]!=null?get(id, v)[0]:""));

        writer.write("\n");
      }
    }
    catch (IOException io)
    {
      throw new RuntimeException("Failed to write the data frame to the output.", io);
    }
  }
  public SortedMap<Double, Double> getXY(String x, String y)
  {
    if (!hasVar(x,false) || hasVar(y,false)) return null;

    SortedMap<Double, Double> line = new TreeMap<>();
    Var vx=getVar(x), vy=getVar(y);
    for (String id : getRowIds())
    {
      Tools.putNotNull(line, cell(id, vx), cell(id, vy));
    }
    return line;
  }
  public double[] getDoubleCol(String y)
  {
    if (!Tools.isSet(mData) || !hasVars(y)) return null;

    double[] ys = new double[getRowIds().size()];
    Var      vy = getVar(y);
    for (int i=0; i<getRowIds().size(); i++)
    {
      ys[i] = Stats.toDouble(cell(getRowIds().get(i), vy));
    }
    return ys;
  }
  public Dataframe addVar(String v, double[] ys)
  {
    if (hasVars(v)) throw new RuntimeException("Variable " + v + " already exist!");

    if (ys!=null && ys.length==getRowIds().size())
    {
      Var vv = addVar(new Variable(v));
      for (int i=0; i<ys.length; i++)
      {
        put(getRowId(i), vv, ys[i]);
      }
    }
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
      mColVars = new ArrayList<>();
      mData    = HashBasedTable.create();
      for (String col : csv.getHeaders()) mColVars.add(new Variable(col));
      // going thro the rows
      long row_counts = 0;
      while (csv.hasNext())
      {
        if (++row_counts % 10000  ==0) System.out.print(".");
        if (  row_counts % 1000000==0) System.out.println();
        String id=null;
        if (Tools.isSet(idcols))
        {
          for (String col : idcols)
            id= Strs.extend(id, csv.get(col), "_");
        }
        else id = row_counts+"";

        addRowId(id);
        // deposite the cells
        for (Var v : mColVars)
          mData.put(id, v, csv.get(v.toString()));
      }
      csv.close();
      System.out.println();
      // setup the types
      init();
    }
    catch (IOException ioe)
    {
      throw new RuntimeException("Unable to access file: " + src, ioe);
    }
  }
  // go thro the table to determine the type of the variables. Convert them to number if necessary
  protected void init()
  {
    if (!Tools.isSet(mData)) return;
    if (mRowIDs ==null) { mRowIDs  = new ArrayList<>(mData.rowKeySet()); Collections.sort(mRowIDs); }
    if (mColVars==null)
    {
      mColVars = new ArrayList<>(mData.columnKeySet());
      mNameVar = new HashMap<>(mColVars.size());
      for (Var v : mColVars) mNameVar.put(v.toString(), v);
    }

    //int counts=0;
    for (Var v : mColVars)
    {
      //System.out.print(v+"...");
      //if (++counts%10==0) System.out.println();
      init(v);
    }
    System.out.println();
  }
  protected void init(Var v)
  {
    if (!Tools.isSet(mData)) return;

    boolean       isNum=true;
    Set<Object> factors=new HashSet<>();
    for (String row : mRowIDs)
    {
      Object val = Stats.toNumber(mData.get(row, v));

      if (val instanceof String) isNum=false;
      if (val!=null && (!(val instanceof String) || ((String )val).length()>0)) factors.add(val);
      // put the cell back
      if (row!=null && v!=null && val!=null) mData.put(row, v, val);
    }
    if (v.isType(Var.VarType.UNKNOWN))
    {
      if (factors.size()>0 && (!isNum || factors.size()<Math.min(250, mRowIDs.size()*0.25)))
        v.setType(Var.VarType.CATEGORICAL);
      else v.setType(Var.VarType.CONTINOUOUS);
    }
    if ( v.isCategorical())
    {
      v.setFactors(factors);
      //Collections.sort(v.getFactors());
    }

    Tools.dispose(factors);
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

    Map<Object, Dataframe> outs = new HashMap<>();
    Var vv = getVar(v);
    for (String r : getRowIds())
    {
      Object key = row(r).get(vv);
      Dataframe F = outs.get(key);
      if (F==null) F = new Dataframe();
      F.addRow(r, mData.row(r));
      outs.put(key, F);
    }
    for (Dataframe d : outs.values()) d.init();

    return outs;
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
    if (!hasVars(rows) || !hasVar(col, true)) return null;
    // looping thro the rows
    Var vcol=getVar(col), vval=getVar(val); Var[] vrows=Variable.toVars(this,rows);
    // build the inventory
    ListMultimap<ArrayKey, Object> body = ArrayListMultimap.create();
    for (String rowid : getRowIds())
    {
      body.put(new ArrayKey(ObjectArrays.concat(get(rowid, vcol)[0], get(rowid, vrows))), get(rowid, vval)[0]);
    }
    // construct the outgoing data frame
    Dataframe out = new Dataframe();
    for (ArrayKey keys : body.keySet())
    {
      String id = Strs.toString(Arrays.copyOfRange(keys.key, 1, keys.key.length), "");
      for (int i=1; i< keys.key.length; i++)
      {
        out.put(id, vrows[i-1], keys.key[i]);
      }
      out.put(id, keys.key[0].toString(), Stats.aggregate(body.get(keys), func));
    }
    out.init(); Tools.dispose(body);
    out.reorder(ObjectArrays.concat(rows, Strs.toStringArray(vcol.getFactors()), String.class));

    return out;
  }
  public TreeBasedTable<Double, Double, String> index(String row, String col)
  {
    if (!hasVar(row,false) || !hasVar(col,false)) return null;

    TreeBasedTable<Double, Double, String> indice = TreeBasedTable.create();
    Var vrow=getVar(row), vcol=getVar(col);
    for (String rowid : getRowIds())
      indice.put(Stats.toDouble(cell(rowid, vrow)), Stats.toDouble(cell(rowid, vcol)), rowid);

    return indice;
  }

  public static TreeBasedTable<Double, Double, String>[] indice(String row, String col, Dataframe... frames)
  {
    if (!Tools.isSet(frames) || !Tools.isSet(row) || !Tools.isSet(col)) return null;

    TreeBasedTable<Double, Double, String>[] indices = new TreeBasedTable[frames.length];
    for (int i=0; i<frames.length; i++)
    {
      indices[i] = frames[i].index(row, col);
    }
    return indices;
  }
}

