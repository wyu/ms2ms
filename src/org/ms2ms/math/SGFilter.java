package org.ms2ms.math;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * Authos: wyu
 * Date: Jan 31, 2007
 * Time: 9:32:19 PM
 * Version:
 * <p/>
 * Synposis:
 * <p/>
 * <p/>
 * Comments:
 */
public class SGFilter extends Object
{
  // key 1 = (data window + 1)/2, key 2 = order of polynomial
  private Table<Integer, Integer, Double[]> mCoeffs_0, mCoeffs_1, mCoeffs_2;

  SGFilter() { super(); init(); }

  private void init()
  {
    mCoeffs_0 = HashBasedTable.create();

    mCoeffs_0.put(1, 1, new Double[] {0.3333333,0.3333333,0.3333333});
    mCoeffs_0.put(2, 2, new Double[] {-0.08571429,0.34285714,0.48571429,0.34285714,-0.08571429});
    mCoeffs_0.put(3, 2, new Double[] {-0.0952381,0.1428571,0.2857143,0.3333333,0.2857143,0.1428571,-0.0952381});
    mCoeffs_0.put(4, 2, new Double[] {-0.09090909,0.06060606,0.16883117,0.23376623,0.25541126,0.23376623,0.16883117,0.06060606,-0.09090909});
    mCoeffs_0.put(5, 2, new Double[] {-0.08391608,0.02097902,0.10256410,0.16083916,0.19580420,0.20745921,0.19580420,0.16083916,0.10256410,0.02097902,-0.08391608});
    mCoeffs_0.put(6, 2, new Double[] {-7.692308e-02,2.953670e-17,6.293706e-02,1.118881e-01,1.468531e-01,1.678322e-01,1.748252e-01,1.678322e-01,1.468531e-01,1.118881e-01,6.293706e-02,2.953670e-17,-7.692308e-02});
    mCoeffs_0.put(7, 2, new Double[] {-0.07058824,-0.01176471,0.03800905,0.07873303,0.11040724,0.13303167,0.14660633,0.15113122,0.14660633,0.13303167,0.11040724,0.07873303,0.03800905,-0.01176471,-0.07058824});
    mCoeffs_0.put(8, 2, new Double[] {-0.06501548,-0.01857585,0.02167183,0.05572755,0.08359133,0.10526316,0.12074303,0.13003096,0.13312693,0.13003096,0.12074303,0.10526316,0.08359133,0.05572755,0.02167183,-0.01857585,-0.06501548});
    mCoeffs_0.put(3, 3, new Double[] {-0.0952381,0.1428571,0.2857143,0.3333333,0.2857143,0.1428571,-0.0952381});
    mCoeffs_0.put(4, 3, new Double[] { -0.09090909,0.06060606,0.16883117,0.23376623,0.25541126,0.23376623,0.16883117,0.06060606,-0.09090909});
    mCoeffs_0.put(5, 3, new Double[] {-0.08391608,0.02097902,0.10256410,0.16083916,0.19580420,0.20745921,0.19580420,0.16083916,0.10256410,0.02097902,-0.08391608});
    mCoeffs_0.put(6, 3, new Double[] {-7.692308e-02,6.900373e-16,6.293706e-02,1.118881e-01,1.468531e-01,1.678322e-01,1.748252e-01,1.678322e-01,1.468531e-01,1.118881e-01,6.293706e-02,-6.081433e-16,-7.692308e-02});
    mCoeffs_0.put(7, 3, new Double[] { -0.07058824,-0.01176471,0.03800905,0.07873303,0.11040724,0.13303167,0.14660633,0.15113122,0.14660633,0.13303167,0.11040724,0.07873303,0.03800905,-0.01176471,-0.07058824});
    mCoeffs_0.put(8, 3, new Double[] {-0.06501548,-0.01857585,0.02167183,0.05572755,0.08359133,0.10526316,0.12074303,0.13003096,0.13312693,0.13003096,0.12074303,0.10526316,0.08359133,0.05572755,0.02167183,-0.01857585,-0.06501548});
    mCoeffs_0.put(9, 3, new Double[] {-0.06015038,-0.02255639,0.01061477,0.03936311,0.06368863,0.08359133,0.09907121,0.11012826,0.11676249,0.11897391,0.11676249,0.11012826,0.09907121,0.08359133,0.06368863,0.03936311,0.01061477,-0.02255639,-0.06015038});
    mCoeffs_1.put(2, 2, new Double[] {-8.333333e-02,6.666667e-01,8.146443e-17,-6.666667e-01,8.333333e-02});
    mCoeffs_1.put(3, 2, new Double[] {-8.730159e-02,2.658730e-01,2.301587e-01,-4.578921e-17,-2.301587e-01,-2.658730e-01,8.730159e-02});
    mCoeffs_1.put(4, 2, new Double[] {});
    mCoeffs_1.put(5, 2, new Double[] {});
    mCoeffs_1.put(6, 2, new Double[] {});
    mCoeffs_1.put(7, 2, new Double[] {});
    mCoeffs_1.put(3, 3, new Double[] {});
    mCoeffs_1.put(4, 3, new Double[] {});
    mCoeffs_1.put(5, 3, new Double[] {});
    mCoeffs_1.put(6, 3, new Double[] {});
    mCoeffs_1.put(7, 3, new Double[] {});
    mCoeffs_1.put(8, 3, new Double[] {});
    mCoeffs_2.put(3, 2, new Double[] {});
    mCoeffs_2.put(4, 2, new Double[] {});
    mCoeffs_2.put(5, 2, new Double[] {});
    mCoeffs_2.put(6, 2, new Double[] {});
    mCoeffs_2.put(7, 2, new Double[] {});
    mCoeffs_2.put(8, 2, new Double[] {});
    mCoeffs_2.put(4, 2, new Double[] {});
  }

  public Double[] getCoeffs(int left, int right, int poly_order, int derivative)
  {

    return null;
  }
  /** Ported from MatchMaker.SavitzkyGolayFilter::setCoefficient
   *
   * @param window
   * @param derivative
   */
  private Double[] calcCoefficient(int window, int derivative)
  {
    double c, i, s, m = (window - 1) / 2;
    List<Double> coeffs = new ArrayList<Double>();

    switch (derivative)
    {
      case 0:
        for (i = 1; i <= window; i++)
        {
          s = i - 1 - m;
          coeffs.add(3 * (3*m*m + 3*m - 1 - 5*s*s)/
                        ((2*m + 3) * (2*m + 1) * (2*m - 1)));
        }
        break;
      case 1:
        for (i = 1; i <= window; i++)
        {
          s = i - 1 - m;
          c  = 5 * (5*s * (3*m*m*m*m +
                           6*m*m*m - 3*m + 1)-
                      7 * (3*m*m + 3*m - 1) * s*s*s);
          c /= ((2*m + 3) * (2*m + 1) * (2*m - 1) *
                  (m + 2) * (m + 1) * m * (m - 1));
          coeffs.add(c);
        }
        break;
      case 2:
        for (i = 1; i <= window; i++)
        {
          s = i - 1 - m;
          c  = 30 * (3*s*s - m * (m + 1));
          c /= ((2*m + 3) * (2*m + 1) *
                (2*m - 1) * (m + 1)*(m));
          coeffs.add(c);
        }
        break;
      case 3:
        for (i = 1; i <= window; i++)
        {
          s = i - 1 - m;
          c  = 210 * ((5*s*s*s) - (3*m*m + 3*m - 1) * s);
          c /= ((2*m + 3) * (2*m + 1) * (2*m - 1)*(m + 2) *
                  (m + 1) * (m) * (m - 1));
          coeffs.add(c);
        }
        break;
    }
    return coeffs.toArray(new Double[]{});
  }

}
