package org.ms2ms.alg;

import com.google.common.collect.Range;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * ** Copyright 2014-2015 ms2ms.org
 * <p/>
 * Description:
 * <p/>
 * Author: wyu
 * Date:   11/21/14
 */
public class Similarity
{
  /**
      * Check for the positive evidence of duplicated spectrum using dot-product
      * algorithm as outlined by
      *     Stein, SE, Scott DR. "Optimization and Testing of Mass Spectral Library
      *     Search Algorithms for Compound Identification", JASMS 1994, 5, 859-866
      *
      *  porting notes: the minor peaks are expected to be tagged via 'invalidate()'
      *                 prior to the call.
      *
      * @param             A = Profile A
      * @param             B = Profile B
      * @param           tol = m/z tolerance in dalton
      * @param        sqrted = transform the Y by sqrt() if TRUE
      *
      * @return dot-product from A to B, 0 if not enough matching points or empty profile in A or B
      */
  public static double dp(List<? extends Peak> A, List<? extends Peak> B, Tolerance tol, boolean highest, boolean sqrted)
  {
    boolean verbose = false;

    // false if one of the spectra is empty
    if (A == null || A.isEmpty() || B == null || B.isEmpty()) return 0;
    // the variables Peak           a, b;
    int    match_cnt = 0, ia     = 0, ib         = 0, start_b;
    double dp        = 0, best_Y = 0, best_error = 1E23, abs_err,
           sum_a     = 0, sum_b  = 0, sum_ab     = 0;

    while (ia < A.size())
    {
      // ignore any invalidated point
      // if (!A.cells(ia).isValid()) { ++ia; continue; }
      if (verbose) System.out.print(">>" + A.get(ia).getMz() + ", " + A.get(ia).getIntensity());
      // reset the goalposts
      best_Y = 0; best_error = 1E23; start_b = -1;
      // recycle to avoid access violation
      if (ib >= B.size()) ib = 0;
      while (ib < B.size())
      {
        // skip the minor points that's labled
        // if (!B.cells(ib).isValid())                               { ++ib; continue; }
        if ( B.get(ib).getMz() > tol.getMax(A.get(ia).getMz()) + 0.1) break;
        if (!tol.withinTolerance(A.get(ia).getMz(), B.get(ib).getMz())) { ++ib; continue; }
        if (start_b == -1) start_b = ib;
        // find the highest peak in the matching window
        if   ( highest && (B.get(ib).getIntensity() > best_Y))
          best_Y = B.get(ib).getIntensity();
          // WYU 10-13-2000, option for matching the closest
        else
        {
          abs_err = Math.abs(B.get(ib).getMz() - A.get(ia).getMz());
          if (!highest && (abs_err  < best_error))
          { best_Y = B.get(ib).getIntensity(); best_error  = abs_err; }
        }
        ++ib;
      }
      if (verbose) System.out.println("    !!" + best_Y);
      if (best_Y > 0) match_cnt++;
      // calculate the summation terms with sqrt() scaling of abundance
      if (sqrted)
      {
        sum_a  +=   A.get(ia).getIntensity();
        sum_b  +=           best_Y;
        sum_ab += Math.sqrt(best_Y * A.get(ia).getIntensity());
      }
      else
      {
        sum_a  +=   A.get(ia).getIntensity() * A.get(ia).getIntensity();
        sum_b  +=           best_Y   * best_Y;
        sum_ab +=   A.get(ia).getIntensity() * best_Y;
      }
      // increment the outer loop
      ++ia;
      // move the b iterator back a little
      if (start_b > 0) { ib = start_b - 1; start_b = -1; } else ib = 0;
    }
    // calculating the dot-product
    if (sum_a > 0 && sum_b > 0)
      dp = sum_ab * sum_ab / (sum_a * sum_b);

    return dp;
  }

  /** Assuming that the registration is already done.
   *
   * @param A
   * @param B
   * @return dot-product
   */
  public static double dp(List<? extends Peak> A, List<? extends Peak> B)
  {
    // false if one of the spectra is empty or with diff dimension
    if (A.isEmpty() || B.isEmpty() || A.size() != B.size()) return 0;

    double dp = 0, sum_a = 0, sum_b = 0, sum_ab = 0;
    for (int ia = 0; ia < A.size(); ia++) {
       sum_a  +=   A.get(ia).getIntensity() * A.get(ia).getIntensity();
       sum_b  +=   B.get(ia).getIntensity() * B.get(ia).getIntensity();
       sum_ab +=   A.get(ia).getIntensity() * B.get(ia).getIntensity();
    }
    // calculating the dot-product
    return sum_ab * sum_ab / (sum_a * sum_b);
  }
  public static double bonanza(List<? extends Peak> A,
                               List<? extends Peak> B,
                               double               delta,
                               double               precursor_AB,
                               double               min_x)
  {
    boolean verbose = false;

    // false if one of the spectra is empty
    if (A.isEmpty() || B.isEmpty()) return 0;

    double mab = 0, uab = 0;
    for (Peak a : A)
    {
      if (a.getMz() < min_x) continue;
      for (Peak b : B)
      {
        if (b.getMz() < min_x) continue;
        // test the match betwwwn the point
        if (Math.abs(a.getMz() - b.getMz()) < delta ||
            Math.abs(precursor_AB - a.getMz() + b.getMz()) < delta)
        {
          // tag them as matched
          Peaks.invalidate(a, b);
          mab += a.getIntensity() * b.getIntensity();
//          if (verbose) System.out.println(">>" + a.getX() + " <--> " + b.getX());
        }
      }
    }
    // account for the unmatched
    for (Peak a : A) if (a.getMz() > min_x && Peaks.isValid(a)) uab += a.getIntensity() * a.getIntensity();
    for (Peak b : B) if (b.getMz() > min_x && Peaks.isValid(b)) uab += b.getIntensity() * b.getIntensity();

    Peaks.validate(A, B);

    return mab + uab != 0 ? (mab / (mab + uab)) : 0;
  }
//  // to allow the offset be adjusted by the possible charge state
//  public static double bonanza(List<? extends XYPoint> A,
//                               List<? extends XYPoint> B,
//                               double                  delta,
//                               double                  precursor_AB,
//                               double                  min_x,
//                               int                     max_frag_z)
//  {
//    boolean verbose = false;
//
//    // false if one of the spectra is empty
//    if (A.isEmpty() || B.isEmpty()) return 0;
//
//    Collection<Double> ladder = new ArrayList<Double>();
//    ladder.add(0d);
//    for (int z = max_frag_z; z > 0; z--)
//      ladder.add(precursor_AB / (double )z);
//
//    double mab = 0, uab = 0;
//    for (XYPoint a : A)
//    {
//      if (a.getX() < min_x) continue;
//      for (XYPoint b : B)
//      {
//        if (b.getX() < min_x) continue;
//        // test the match betwwwn the point
//        //if (Math.abs(a.getX() - b.getX()) < delta ||
//        //    Math.abs(precursor_AB - a.getX() + b.getX()) < delta)
//        if (Toolbox.nearA(Math.abs(a.getX() - b.getX()), delta, ladder))
//        {
//          // tag them as matched
//          a.invalidate(); b.invalidate();
//          mab += a.getY() * b.getY();
//          if (verbose) System.out.println(">>" + a.getX() + " <--> " + b.getX());
//        }
//      }
//    }
//    // account for the unmatched
//    for (XYPoint a : A) if (a.getX() > min_x && a.isValid()) uab += a.getY() * a.getY();
//    for (XYPoint b : B) if (b.getX() > min_x && b.isValid()) uab += b.getY() * b.getY();
//
//    validate(A); validate(B);
//
//    return mab + uab != 0 ? (mab / (mab + uab)) : 0;
//  }

  //--------------------------------------------------------------------------
  /** Produce a list of the signature peaks that can be used to screen candidate profile
   *  quickly. It's also important to check the integrity of the signature aside from dot-product
   *  criteria so small, yet significant changes due to residue swapping, etc will be noticed.
   *
   * @param              A = Profile
   * @param      n_regions = # of sub-regions into which the profie would be divided.
   *                         The base point from each sub-region is selected as the index point.
   * @param         n_tops = # of the most intense points that's not already selected from
   *                         the sub-region, to be added to the indices.
   * @param min_separation = min separation between the index points in dalton
   * @return The index points as a List
   */
  public static <T extends Peak> List<Peak> index(List<T> A, int n_regions, int n_tops, double min_separation, double min_ai)
  {
    Collections.sort(A, new Peaks.MzAscendComparator());

    // initiating a new List
    List<Peak> indices = new ArrayList<>();
    // required variables
    double range = A.get(A.size()-1).getMz() - A.get(0).getMz(), last_X = 0;

    // setup the regions by mz width
    if      (n_regions  < 0) n_regions = (int )Math.round(range / Math.abs(n_regions));
      // default to 5 as min# of regions
    else if (n_regions == 0 || n_regions  < 5) n_regions = 5;
    // default to 3 as min# of tops
    if      (n_tops    <= 0) n_tops    = 3;
    // not enough points to make a meaningful indices
    if (A.size() < n_regions * 2 + n_tops) return indices;
    // inspect the regions
    int left = 0; double step = range / n_regions;
    for (double bound  = A.get(0).getMz() + step;
         bound <= A.get(A.size()-1).getMz() + step * 0.1; bound += step)
    {
      double top_Y = 0; int top_i = -1;
      while (left < A.size() && A.get(left).getMz() <= bound)
      {
        // skip the invalidated point, WYU 111607
        if (!Peaks.isValid(A.get(left))) { ++left; continue; }
        if (         A.get(left).getIntensity()    > top_Y &&
            Math.abs(A.get(left).getMz() - last_X) > min_separation)
        { top_i = left; top_Y = A.get(left).getIntensity(); }
        ++left;
      }
      if (top_i >= 0 && top_Y >= min_ai)
      {
        indices.add(A.get(top_i)); last_X = A.get(top_i).getMz();
      }
    }
    // prepare the points so the most intense points is at the front
    Collections.sort(A, new Peaks.IntensityDesendComparator());
    // grab from the most intense points that are not already in
    int top_cnt = 0;
    for (left = 0; left < A.size(); left++)
    {
      // skip the invalidated point, WYU 111607
      if (!Peaks.isValid(A.get(left))) continue;
      // make sure it's not in the indices already
      boolean matched = false;
      for (Peak p : indices)
        if (Math.abs(p.getMz() - A.get(left).getMz()) <= min_separation) { matched = true; break; }
      if (!matched) { indices.add(A.get(left)); top_cnt++; }

      if (top_cnt >= n_tops) break;
    }
    // putting the points back to the natural order
    Collections.sort(A,       new Peaks.MzAscendComparator());
    Collections.sort(indices, new Peaks.MzAscendComparator());

    return indices;
  }
  // hypergeometric distribution: scores overlap between the position of predicted and observed peaks
  public static double hg_similarity(List<? extends Peak> A, List<? extends Peak> B,
                                     double delta, long n_bins, boolean matchHighest, Map<Peak, Peak> outcomes)
  {
    outcomes = Peaks.overlap(A, B, delta, matchHighest, false, outcomes);
    int matched = 0, na = 0, nb = 0;
    for (Peak p : A) na += Peaks.isC12(p) ? 1 : 0;
    for (Peak p : B) nb += Peaks.isC12(p) ? 1 : 0;

    if (outcomes != null && outcomes.size() > 0)
      for (Map.Entry<Peak, Peak> E : outcomes.entrySet())
        if (Peaks.isC12(E.getKey()) && Peaks.isC12(E.getValue())) matched++;

    return -1 * Stats.hypergeometricPval1(matched, nb, nb, n_bins);
  }
  public static double hg_similarity(List<? extends Peak> A, List<? extends Peak> B, double delta)
  {
    Range<Double> mz_range = Range.closed(Math.min(Tools.front(A).getMz(), Tools.front(B).getMz()),
                                          Math.max(Tools.back( A).getMz(), Tools.back( B).getMz()));
    Map<Peak, Peak> outcomes = Peaks.overlap(A, B, delta, true, false, null);
    return -1 * Stats.hypergeometricPval1(outcomes.size(), A.size(), B.size(),
        (long )((mz_range.upperEndpoint() - mz_range.lowerEndpoint()) / delta));
  }
}
