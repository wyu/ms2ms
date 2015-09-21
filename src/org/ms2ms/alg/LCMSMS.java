package org.ms2ms.alg;

import org.expasy.mzjava.core.ms.spectrum.ScanNumber;
import org.expasy.mzjava.core.ms.spectrum.ScanNumberList;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

/**
 * Created by yuw on 9/20/2015.
 */
public class LCMSMS
{
  public static String toScanStr(ScanNumberList scans)
  {
    String scan=null;
    if (Tools.isSet(scans))
    {
      for (ScanNumber s : scans)
        scan = Strs.extend(scan, s.getValue()+"", "+");
    }
    return scan;
  }
}
