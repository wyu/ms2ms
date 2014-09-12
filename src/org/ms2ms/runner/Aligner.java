package org.ms2ms.runner;

import com.google.common.collect.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.expasy.mzjava.core.ms.Tolerance;
import org.ms2ms.data.Dataframe;
import org.ms2ms.data.Features;
import org.ms2ms.data.Var;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Strs;
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
//  protected Table<Dataframe, String, Var> mStr2Var;

  protected Table<Long, Dataframe, String> mAlignment;      // id, frame and rowid to the starting data frame
  protected List<TreeMap<Double, Long>>    mIndice;

  public Aligner(Tolerance[] tols, String... cols)
  {
    mCols=cols; mTols=tols;
  }

  public int getNumVars() { return mTols!=null?mTols.length:0; }

  public String print()
  {
    StringBuffer buf = new StringBuffer();

    // pre-fetch the cols and rows to save time
    Collection<Dataframe> cols = mAlignment.columnKeySet();
    Collection<Long>      rows = mAlignment.rowKeySet();

    // print the headers
    for (Dataframe frm : cols)
    {
      buf.append(Strs.rtuncate(frm.getTitle(), 4) + "\t");
    }
    buf.append("\n");

    for (Long row : rows)
    {
      for (Dataframe frm : cols)
      {
        Object val = frm.cell(mAlignment.get(row, frm), "m/z");
        buf.append((val!=null?val.toString():"--") + "\t");
      }
      buf.append("\n");
    }

    return buf.toString();
  }
  class F2F implements Comparable<F2F>
  {
    private Long mID1, mID2;
    private String mRow1, mRow2;
    private double mScore;

    protected F2F(String from, String to, Long id2, double s)
    {
      mRow1=from; mRow2=to; mID2=id2; mScore=s;
    }

    /*** This method returns the peak list row which is being aligned */
    public String getRow1() { return mRow1; }
    /*** This method returns the row of aligned peak list */
    public String getRow2() { return mRow2; }
    public Long   getID2()  { return mID2; }
    /*** This method returns score between the these two peaks (the lower score, the better match) */
    double getScore() { return mScore; }

    F2F setRow1(String s) { mRow1=s; return this; }
    F2F setRow2(String s) { mRow2=s; return this; }
    F2F setID1( Long   s) { mID1 =s; return this; }
    F2F setID2( Long   s) { mID2 =s; return this; }

    /*** @see java.lang.Comparable#compareTo(java.lang.Object) */
    public int compareTo(F2F object)
    {
      // We must never return 0, because the TreeSet in JoinAlignerTask would treat such elements as equal
      return mScore < object.getScore() ? 1 : -1;
    }
  }
  public double match(Dataframe A, String rowidA, Dataframe B, String rowidB)
  {
    double score = 1d;
    for (int i=0; i<mCols.length; i++)
    {
      double x0 = Stats.toDouble(A.cell(rowidA, mCols[i])), delta = Math.abs(x0 - Stats.toDouble(B.cell(rowidB, mCols[i])));
      NormalDistribution norm = new NormalDistribution(0, (mTols[i].getMax(x0)-mTols[i].getMin(x0))/1.77d);
      score *= norm.density(delta); norm=null;
    }
    return score;
  }
  public <T extends Features> void run(Dataframe... traces)
  {
    mAlignment = HashBasedTable.create(); // id, frame and rowid
    mIndice    = new ArrayList<TreeMap<Double, Long>>(getNumVars());

    // map the col header and variables
//    mStr2Var = HashBasedTable.create();
//    for (Dataframe d : traces)
//        for (String m : mCols)
//            mStr2Var.put(d, m, d.getVar(m));

    // Iterate source peak lists
    int processedRows = 0; long maxid=0; boolean fresh=true; Collection<Long> candidates = new HashSet<Long>();
    for (int i=0; i<traces.length; i++)
    {
      System.out.print("Frame " + (i + 1) + ": " + traces[i].size());
      // the variables
      Var[] vs=traces[i].toVars(mCols);
      // Create a sorted set of scores matching
      TreeSet<F2F> scoreSet = new TreeSet<F2F>();
      // Calculate scores for all possible alignments of this row
      for (String rowid : traces[i].getRowIds())
      {
        // collect the rows from the current alignment that match to the row in consideration within the tolerances
        candidates.clear(); fresh=true;
        if (mIndice.size()>0)
          for (int k=0; k<getNumVars(); k++)
          {
            Double                x = Stats.toDouble(traces[i].cell(rowid,vs[k]));
            Map<Double, Long> slice = mIndice.get(k).subMap(mTols[k].getMin(x), mTols[k].getMax(x));
            if (slice==null)
            {
              // setup the new candidates
              if  (fresh) candidates.addAll(slice.values());
                // retain the shared rows only
              else candidates=Tools.common(candidates, slice.values());
            }
            fresh=false;
          }

        // Calculate scores and store them
        if (Tools.isSet(candidates))
          for (Long candidate : candidates)
            for (Dataframe d : mAlignment.row(candidate).keySet())
            {
              double score = match(traces[i], rowid, d, mAlignment.get(candidate, d));
              if (score!=0) scoreSet.add(new F2F(rowid, mAlignment.get(candidate, d), candidate, score));
            }

        processedRows++;
      }

      // Create a table of mappings for best scores
      BiMap<String, F2F> mapping = HashBiMap.create();

      // Iterate scores by descending order
      for (F2F score : scoreSet)
        mapping = Tools.putNew(mapping, score.getRow1(), score);

      // Align all rows using mapping
      for (String rowid : traces[i].getRowIds())
      {
        F2F target = mapping.get(rowid);

        // If we have no mapping for this row, add a new one
        if (target == null)
        {
          target = new F2F(rowid, null, maxid++, 0);
        }
        target.setRow2(rowid);
        // Add all peaks from the original row to the aligned row
        mAlignment.put(target.getID2(), traces[i], rowid);
        // update the indice
        for (int k=0; k<getNumVars(); k++)
        {
          if (mIndice.size()-1<k) mIndice.add(new TreeMap<Double, Long>());
          mIndice.get(k).put(Stats.toDouble(traces[i].cell(rowid,vs[k])), target.getID2());
        }

        processedRows++;
      }
      System.out.println("\t" + mAlignment.rowKeySet().size());
    } // Next peak list
  }
}
