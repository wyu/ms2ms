package org.ms2ms.data.ms;

import java.util.ArrayList;
import java.util.Collection;

/** Measurement and its transformation
 *
 */
public class Measurable implements Comparable<Measurable>, Cloneable
{
  @Override
  public int compareTo(Measurable o) { return Double.compare(value, o.value); }

  enum eType {RAW, OBS, CALC, CALIB, RESIDUAL, UNKNOWN };
  enum eUnit {DALTON, SEC, MIN, NET, UNKNOWN };

  // prev step of the transformation
  public Measurable prev=null;

  public double value;
  public eType  type = eType.UNKNOWN;
  public eUnit  unit = eUnit.UNKNOWN;

  private Collection<Measurable> transformed = null;

  public Measurable() {}
  public Measurable(double s, eType t, eUnit u) { value=s; type=t; unit=u; }

  public Measurable transform(double val, eType t, eUnit u)
  {
    if (transformed==null) transformed = new ArrayList<>();

    prev = clone(); transformed.add(prev); value=val; type=t; unit=u;
    return this;
  }

  public Measurable clone()
  {
    Measurable cloned = new Measurable(value, type, unit);
    cloned.transformed=transformed; cloned.prev=prev;

    return cloned;
  }
  public static Measurable newMz(double s) { return new Measurable(s, eType.RAW, eUnit.DALTON); }
  public static Measurable newRt(double s) { return new Measurable(s, eType.RAW, eUnit.MIN); }
}
