package org.ms2ms.data;

import java.util.Collection;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/13/14
 * Time: 12:08 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Var
{
  public enum VarType { CONTINOUOUS, CATEGORICAL, UNKNOWN }

  public boolean isCategorical();
  public boolean isContinuous();
  public Var setType(VarType s);
  public Var setFactors(Collection s);
  public List getFactors();
  public boolean isType(VarType s);
}
