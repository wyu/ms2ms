package org.ms2ms.data.ms;

import org.ms2ms.Disposable;
import org.ms2ms.data.Binary;

// A contract representing an actual PSM during the search. Not using 'interface' since we need to keep some variables here
//
abstract public class PSM implements Disposable, Binary
{
  public abstract String getPeptide();
  public abstract String getSequence();
}
