package org.ms2ms.runner;

import com.google.common.collect.TreeBasedTable;
import net.sf.mzmine.data.IsotopePattern;
import net.sf.mzmine.data.PeakListRow;
import net.sf.mzmine.data.RawDataFile;
import net.sf.mzmine.data.impl.SimplePeakListRow;
import net.sf.mzmine.modules.peaklistmethods.alignment.join.JoinAlignerParameters;
import net.sf.mzmine.modules.peaklistmethods.alignment.join.RowVsRowScore;
import net.sf.mzmine.modules.peaklistmethods.isotopes.isotopepatternscore.IsotopePatternScoreCalculator;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.util.PeakUtils;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.Dataframe;
import org.ms2ms.data.Features;
import org.ms2ms.data.Var;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * User: wyu
 * Date: 7/23/14
 */
public class Aligner
{
  protected Tolerance[] mTols;
  protected String[]    mCols;

  class F2F implements Comparable<F2F>
  {
    private Features mRow1, mRow2;
    private double mScore;

    protected F2F(Features from, Features to)
    {
      mRow1=from; mRow2=to;
      mScore=mRow1.score(mRow2);
    }

    /*** This method returns the peak list row which is being aligned */
    public Features getRow1() { return mRow1; }
    /*** This method returns the row of aligned peak list */
    public Features getRow2() { return mRow2; }

    /*** This method returns score between the these two peaks (the lower score, the better match) */
    double getScore() { return mScore; }

    /*** @see java.lang.Comparable#compareTo(java.lang.Object) */
    public int compareTo(F2F object)
    {
      // We must never return 0, because the TreeSet in JoinAlignerTask would treat such elements as equal
      return mScore < object.getScore() ? 1 : -1;
    }
  }
  static public double match()
  {

  }

  public static <T extends Features> void align(String r1, String r2, Tolerance t1, Tolerance t2, Dataframe... traces)
  {
    Dataframe                            aligned = new Dataframe();
    TreeBasedTable<Double, Double, String> index = aligned.index(r1, r2);

    // Iterate source peak lists
    for (int i=0; i<traces.length; i++)
    {
      // the variables
      Var v1=traces[i].getVar(r1), v2=traces[i].getVar(r2);
      // Create a sorted set of scores matching
      TreeSet<F2F> scoreSet = new TreeSet<F2F>();
      // Calculate scores for all possible alignments of this row
      for (String rowid : traces[i].getRowIds())
      {
        Double x1= Stats.toDouble(traces[i].cell(rowid,v1)), x2= Stats.toDouble(traces[i].cell(rowid, v2));
        Collection<String> candidates = Tools.slice(index, t1.getMin(x1), t1.getMax(x1), t2.getMin(x2), t2.getMax(x2));

        // Calculate scores and store them
        for (String candidate : candidates)
        {
          if (sameChargeRequired) {
            if (!PeakUtils.compareChargeState(row, candidate))
              continue;
          }

          if (sameIDRequired) {
            if (!PeakUtils.compareIdentities(row, candidate))
              continue;
          }

          if (compareIsotopePattern) {
            IsotopePattern ip1 = row.getBestIsotopePattern();
            IsotopePattern ip2 = candidate.getBestIsotopePattern();

            if ((ip1 != null) && (ip2 != null)) {
              ParameterSet isotopeParams = parameters
                  .getParameter(
                      JoinAlignerParameters.compareIsotopePattern)
                  .getEmbeddedParameters();

              if (!IsotopePatternScoreCalculator.checkMatch(ip1,
                  ip2, isotopeParams)) {
                continue;
              }
            }
          }

          RowVsRowScore score = new RowVsRowScore(row, candidate,
              mzRange.getSize() / 2, mzWeight,
              rtRange.getSize() / 2, rtWeight);

          scoreSet.add(score);

        }

        processedRows++;

      }

      // Create a table of mappings for best scores
      Hashtable<PeakListRow, PeakListRow> alignmentMapping = new Hashtable<PeakListRow, PeakListRow>();

      // Iterate scores by descending order
      Iterator<RowVsRowScore> scoreIterator = scoreSet.iterator();
      while (scoreIterator.hasNext()) {

        RowVsRowScore score = scoreIterator.next();

        // Check if the row is already mapped
        if (alignmentMapping.containsKey(score.getPeakListRow()))
          continue;

        // Check if the aligned row is already filled
        if (alignmentMapping.containsValue(score.getAlignedRow()))
          continue;

        alignmentMapping.put(score.getPeakListRow(),
            score.getAlignedRow());

      }

      // Align all rows using mapping
      for (PeakListRow row : allRows) {

        PeakListRow targetRow = alignmentMapping.get(row);

        // If we have no mapping for this row, add a new one
        if (targetRow == null) {
          targetRow = new SimplePeakListRow(newRowID);
          newRowID++;
          alignedPeakList.addRow(targetRow);
        }

        // Add all peaks from the original row to the aligned row
        for (RawDataFile file : row.getRawDataFiles()) {
          targetRow.addPeak(file, row.getPeak(file));
        }

        // Add all non-existing identities from the original row to the
        // aligned row
        PeakUtils.copyPeakListRowProperties(row, targetRow);

        processedRows++;

      }

    } // Next peak list

  }
}
