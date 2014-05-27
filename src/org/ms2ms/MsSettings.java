package org.ms2ms;

import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.utils.Settings;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/24/14
 * Time: 11:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class MsSettings extends Settings
{
  public static final String SPEC_TYPE = "SpecType";
  public static final String TOL_PREC  = "PrecTol";
  public static final String TOL_FRAG  = "FragTol";
  public static final String FRAG_MODE = "FragMode";
  public static final String Z_FLOAT   = "Zfloat";

  public MsSettings() { super(); properties.put(Z_FLOAT, 0); }
  public MsSettings(byte[] type, Tolerance precursor, Tolerance frag, String mode)
  {
    super();
    properties.put(SPEC_TYPE, type);
    properties.put(TOL_PREC,  precursor);
    properties.put(TOL_FRAG,  frag);
    properties.put(FRAG_MODE, mode);
    properties.put(Z_FLOAT,   0);
  }

  public byte[]    getSpecType()     { return getBytes(SPEC_TYPE); }
  public Tolerance getPrecursorTol() { return properties!=null?(Tolerance )properties.get(TOL_PREC):null; }
  public Tolerance getFragmentTol()  { return properties!=null?(Tolerance )properties.get(TOL_FRAG):null; }
  public String    getFragMode()     { return getString(FRAG_MODE); }
  public int       getZFloat()       { return getInteger(Z_FLOAT); }
}
