package org.ms2ms.data;

import com.google.common.collect.*;
import org.ms2ms.utils.ArrayKey;
import org.ms2ms.utils.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.Arrays;

/** Representation of high-dimensional omic dataset. L-shape struct as described in ArrayStudio doc.
 *
 *       Sample (X)
 *
 *   V   Measurements (Y)       Y is a continuous number for expression/CNV-log2/mutation-freq/
 *   a                            imputed-SNP:allele-dose-or-3-genotype-probs/methylation:freqs, or
 *   r                            genotype for SNP/genotyping, integer (code index) for CNV allele.
 *  (A)
 *
 * User: wyu
 * Date: 7/13/14
 */
public class HData
{
  private Dataframe mSamples, mVariables;
  private Table<String, String, Object> mMeasurements;

  public HData() { super(); }
  public HData(Dataframe data, String result, String sample, String variable, Range<Integer>... samples)
  { init(data, result, sample, variable, samples); }

  // col2idx returns index in variable matching 'col'
  // row2idx returns index in samples  matching 'row'


  /** From a long-skinny data frame with one observation per row. Duplicate observations should be coded
   *  with a separate 'duplicate' column in sample-description section to avoid exception
   *
   *  Date      Sample  Test  Visit Subject Result
   *  8/21/2003 Z702    CD3   D0    23566   700
   *  8/21/2003 X780    CD3   D0    23566   650
   *  8/20/2003 X780    CD8   D100  23566   1090
   *  8/20/2003 Z702    CD8   D100  23566   1001
   *
   *  HData(data, "Result", "Sample", "Test", new Range<Integer>(1,5))
   * @param data is the source data frame
   * @param sample
   * @param variable
   * @param samples
   * @return
   */
  public HData init(Dataframe data, String result, String sample, String variable, Range<Integer>... samples)
  {
    if (data==null || !data.hasVar(result,false) || !data.hasVar(sample,true) || !data.hasVar(variable,true)) return this;

/*
    Var vrslt=data.getVar(result), vsmpl=data.getVar(sample), vvar=data.getVar(variable);
    // build the inventory
    ListMultimap<ArrayKey, Object> body = ArrayListMultimap.create();
    for (String rowid : getRowIds())
    {
      body.put(new ArrayKey(ObjectArrays.concat(get(rowid, vcol)[0], get(rowid, vrows))), get(rowid, vval)[0]);
    }
    // construct the outgoing data frame
    Dataframe out = new Dataframe();
    for (ArrayKey keys : body.keySet())
    {
      String id = Strs.toString(Arrays.copyOfRange(keys.key, 1, keys.key.length), "");
      for (int i=1; i< keys.key.length; i++)
      {
        out.put(id, vrows[i-1], keys.key[i]);
      }
      out.put(id, keys.key[0].toString(), Stats.aggregate(body.get(keys), func));
    }
    out.init(); Tools.dispose(body);
    out.reorder(ObjectArrays.concat(rows, Strs.toStringArray(vcol.getFactors()), String.class));
*/
    return this;
  }
}
