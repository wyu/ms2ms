package org.ms2ms.algo;

import org.ms2ms.data.ms.PSM;

public interface Calib
{
  Calib addCalibrant(String tag, PSM T, Double rt, Double ai);
  Calib addPositive(PSM T, Double rt, Double ai);
  Calib addNegative(PSM T, Double rt);

}
