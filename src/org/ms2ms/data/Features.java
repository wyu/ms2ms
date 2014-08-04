package org.ms2ms.data;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.utils.Stats;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/** A generic class to model rows of data in a data table. Used in an alignment routine
 *
 * User: wyu
 * Date: 7/13/14
 */
public class Features
{
/*
  private Dataframe mData;
  private Var[]            mVariable;
  private Tolerance[]      mTols;

  public Features() { super(); }
  public Features(Tolerance tol, Var[] vars)
  {
    super();
    mVariable   = vars;
    mTols       = Arrays.copyOf(new Tolerance[] {tol}, vars.length); // one relative tolerance for all variables
  }
  public Features(Tolerance[] tol, Var... vars)
  {
    super();
    if (tol.length!=vars.length) throw new RuntimeException("Equal number of tolerances and variables required!");

    mVariable  = vars;
    mTols      = tol; // one relative tolerance for all variables
  }
  public Features(Dataframe d, Tolerance[] tol, Var... vars)
  {
    super();

    if (tol.length!=vars.length) throw new RuntimeException("Equal number of tolerances and variables required!");

    mVariable  = vars;
    mTols      = tol; // one relative tolerance for all variables
    mData      = d;
  }

  public Var       var(int i) { return mVariable!=null?mVariable[i]:null; }
  public Tolerance tol(int i) { return mTols!=null?mTols[i]:null; }

  public String  getID()            { return null; }
  public Object  get(int i,Var key) { return mData!=null?mData.row(i).get(key):null; }
*/
/*
  public String  getString(Var key) { return (String )get(key); }
  public Double  getDouble(Var key) { return (Double )get(key); }
  public Float   getFloat( Var key) { return (Float  )get(key); }
  public Long    getLong(  Var key) { return (Long   )get(key); }
  public Integer getInt(   Var key) { return (Integer)get(key); }
*//*


  public Features setVars(Var... vs) { mVariable = vs; return this; }
  public Features setTols(Tolerance... ts)
  {
    mTols = ts;
    return this;
  }

  public boolean match(Features s)
  {
    for (int i=0; i<mVariable.length; i++)
      if (!mTols[i].withinTolerance(Stats.toDouble(get(var(i))), Stats.toDouble(s.get(s.var(i))))) return false;

    return true;
  }
  public double score(Features s)
  {
    double score = 1d;
    for (int i=0; i<mVariable.length; i++)
    {
      double x0 = Stats.toDouble(get(var(i))), delta = x0 - Stats.toDouble(s.get(s.var(i)));
      NormalDistribution norm = new NormalDistribution(0, (tol(i).getMax(x0)-tol(i).getMin(x0))/1.77d);
      score *= norm.density(delta); norm=null;
    }
    return score;
  }
*/
}
