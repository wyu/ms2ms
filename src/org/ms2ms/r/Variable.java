package org.ms2ms.r;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/13/14
 * Time: 1:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class Variable implements Var
{
  private String       mName;
  private List<Object> mFactors;
  private VarType      eType = VarType.UNKNOWN;


  public Variable(String s)            { setName(s); }
  public Variable(String s, VarType t) { setName(s); eType=t; }
  @Override
  public boolean isCategorical() { return isType(VarType.CATEGORICAL); }
  @Override
  public boolean isContinuous() { return isType(VarType.CONTINOUOUS); }
  @Override
  public boolean isType(VarType s) { return eType.equals(s);}
  @Override
  public String toString() { return mName; }

  public Var setType(VarType s) { eType=s; return this; }
  public Var setFactors(Collection s)
  {
    if (s==null) { mFactors=null; return this; }

    if (mFactors==null) mFactors = new ArrayList();
    mFactors.addAll(s);
    return this;
  }
  public Variable setName(String s) { mName=s; return this; }
  public List getFactors() { return mFactors; }

  public boolean equals(Variable s)
  {
    return toString().equals(s.toString()) && eType.equals(s.eType);
  }

  //**************  Utils  ***********************//

  public static Var[] toVars(Dataframe data, String... vs)
  {
    Var[] vrows = new Var[vs.length];
    for (int i=0; i<vs.length; i++) vrows[i]=data.getVar(vs[i]);

    return vrows;
  }
}
