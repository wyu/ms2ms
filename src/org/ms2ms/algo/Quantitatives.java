package org.ms2ms.algo;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.Tools;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/15/14
 * Time: 9:27 PM
 * To change this template use File | Settings | File Templates.
 */
public class Quantitatives
{
  public static void align()
  {

  }
  /** Extract the RT calibration curves from the key alignment
   *
   * @param evidences is a data frame from Maxquant containing the key alignment
   * @param cols are the column headers for key, Rt and RT deviation in order
   * @return a new data frame containing the RT calibration curves
   */
  public static Dataframe toRtCalibrationCurves(Dataframe evidences, String... cols)
  {
    if (evidences==null || !Tools.isSet(cols) || cols.length<3) return null;

    String peptide=cols[0], rt=cols[1], delta=cols[2];
    Multimap<String, Object> peptides = HashMultimap.create();
    Dataframe deviations = new Dataframe();



    return deviations;
  }
}
