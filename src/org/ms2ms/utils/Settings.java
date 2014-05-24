package org.ms2ms.utils;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 5/20/14
 * Time: 10:36 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class Settings implements Serializable
{
  protected Map<String, Object> properties;

  public Settings() { super(); properties = new HashMap<String, Object>();  }

  protected Character getChar(   String s) { return properties!=null?(Character )properties.get(s):0; }
  protected Double    getDouble( String s) { return properties!=null?(Double    )properties.get(s):null; }
  protected String    getString( String s) { return properties!=null?(String    )properties.get(s):null; }
  protected Long      getLong(   String s) { return properties!=null?(Long      )properties.get(s):null; }
  protected Integer   getInteger(String s) { return properties!=null?(Integer   )properties.get(s):null; }
  protected byte[]    getBytes(  String s) { return properties!=null?(byte[]    )properties.get(s):null; }

}
