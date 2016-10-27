package org.ms2ms.algo;

import com.google.common.collect.Range;
import com.hfg.bio.Element;
import com.hfg.bio.ms.Peptide;
import com.hfg.bio.ms.Protein;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.data.ms.IsoEnvelope;
import org.ms2ms.math.Histogram;
import org.ms2ms.utils.Tools;

import java.util.*;
import java.io.*;

/** Adaptation of emass program from Alan Rockwood
 * User: wyu
 * Date: Apr 8, 2010
 * Time: 1:13:41 PM
 *
 emass.cpp: Calculation of accurate masses and intensities of isotopic peaks

 Based on an algorithm developed by Alan L. Rockwood.

 Published in
 Rockwood, A.L. and Haimi, P.: "Efficent calculation of Accurate Masses of Isotopic Peaks",
 Journal of The American Society for Mass Spectrometry
 JASMS 03-2263, in press

 Copyright (c) 2005 Perttu Haimi and Alan L. Rockwood

 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided
 that the following conditions are met:

     * Redistributions of source code must retain the above copyright notice, this list of conditions
       and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions
       and the following disclaimer in the documentation and/or other materials provided with the distribution.
     * Neither the author nor the names of any contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
public class Isotopes
{
  // some inner classes
  class Pattern
  {
    public List<Peak> pattern = new ArrayList<Peak>();
  }
  class SuperAtomList
  {
    public List<Pattern> super_atom_list = new ArrayList<Pattern>();
  }
  //class SuperAtomData
  //{
  //  public List<SuperAtomList> super_atom_data;
  //}

  //typedef std::vector<peak> Pattern;              // index: peak_number
  //typedef std::vector<Pattern> SuperAtomList;        // index: bit_number
  //typedef std::vector<SuperAtomList> SuperAtomData;  // index: element_number

  static final double   ELECTRON_MASS = 0.00054858;
  static final double      DUMMY_MASS = -10000000;
  static final double  AVERAGINE_MASS = 111.0543052;
  public static final double       DELTA_C13 = 1.00335d;
  static Map<String, Long>         EM = new HashMap<String, Long>();
  static List<List<List<Peak>>>   SAD = new ArrayList<List<List<Peak>>>();

  public static Histogram dp_prediction = new Histogram();

  static
  {
    try
    {
      //init_data("/home/wyu/contrib/Rockwood/emass/ISOTOPE.DAT");
      init_data("/Users/yuw/Documents/Apps/contrib/ms2ms/data/ISOTOPE.DAT");
    }
    catch (Exception e) { throw new RuntimeException("Not able to initialize the isotope util: ", e); }
  }
  // The averagine formula
  // C4.9384 H7.7583 N1.3577 O1.4773 S0.0417
  // M(mono) = 111.0543052
  // M(avg)  = 111.1237368
  
  /** suck in the parameters
X  2
1  0.9
2  0.1

H  2
1.0078246  0.99985
2.0141021  0.00015
  */
  public static boolean init_data(String filename) throws IOException
  {
    return init_data(filename, SAD, EM);
  }

  private static boolean init_data(String filename, List<List<List<Peak>>> sad, Map<String, Long> em) throws IOException
  {
    BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename)));

    sad.clear();
     em.clear();

    long elemindex = 0; int state = 0;
    List<List<Peak>> current = null;
    while(reader.ready())
    {
      String[] strs = reader.readLine().split("\\s+");
      String element;
      if (state == 0)
      {
        element = strs[0];
        em.put(element, elemindex);
        current = new ArrayList<List<Peak>>();
        sad.add(current);
        elemindex++;
        state = 1;
      }
      else if (state == 1) // isotope
      {
        try
        {
          Peak p = new Peak(new Double(strs[0]), new Double(strs[1]));
          if (current.size() == 0) current.add(new ArrayList<Peak>());

          List<Peak> idist = current.get(0);
          // fill the gaps in the patterns with zero abundancy peaks
          if(idist.size() > 0)
          {
            double prevmass = Tools.back(idist).getMz();
            for(int i = 0; i < (int )(p.getMz() - prevmass - 0.5); i++) idist.add(new Peak(DUMMY_MASS, 0));
          }
          // insert the peak
          idist.add(p);
        }
        catch (NumberFormatException e)
        {
          //e.printStackTrace();
          state = 0; // no more isotope data
        }
      }
    }
    reader.close();
    EM = em;
    
    return true;
  }

// Merge two patterns to one.
  public static void convolute_basic(List<Peak> h, final List<Peak> g, final List<Peak> f)
  {
    h.clear();
    int g_n = g.size(), f_n = f.size();
    if(g_n == 0 || f_n == 0) return;

    for (int k = 0; k < g_n + f_n - 1; k++)
    {
      double sumweight = 0, summass = 0;
      int start = k < (f_n - 1) ? 0 : k - f_n + 1, // max(0, k-f_n+1)
            end = k < (g_n - 1) ? k : g_n - 1;       // min(g_n - 1, k)
      for(int i = start; i <= end; i++)
      {
        double weight = g.get(i).getIntensity() * f.get(k - i).getIntensity();
        double   mass = g.get(i).getMz()        + f.get(k - i).getMz();
        sumweight += weight;
        summass   += weight * mass;
      }
      h.add(new Peak(sumweight == 0 ? DUMMY_MASS : summass / sumweight, sumweight));
    }
  }

  // Prune the small peaks from both sides but leave them within the pattern.
  public static void prune(List<Peak> f, double limit)
  {
    // prune the front
    ListIterator<Peak> itr = f.listIterator();
    while (itr.hasNext()) if(itr.next().getIntensity() <= limit) itr.remove();

    // prune the end
    itr = f.listIterator();
    while(itr.hasPrevious())
    {
      if(f.size() == 0) break;
      if(itr.previous().getIntensity() <= limit) itr.remove();
    }
  }

  // script: RMS (weighted) (ATGC)1000
  public static void print_pattern(List<Peak> result, int digits)
  {
    // find the maximum
    double max_area = 0;
    for(Peak pk : result)
    {
      if(max_area < pk.getIntensity()) max_area = pk.getIntensity();
    }
    if(max_area == 0) return; // empty pattern

    double print_limit = Math.pow(10.0, -digits) / 2;

    for(Peak pk : result)
    {
      double mass = pk.getMz();
      double rel_area = pk.getIntensity();
      double val_perc = rel_area / max_area * 100;

      if(mass != DUMMY_MASS && val_perc >= print_limit)
        System.out.println(Tools.d2s(mass, digits) + " " + Tools.d2s(val_perc, digits));
    }
  }

  synchronized public static List<Peak> calculate(List<Peak> tmp, List<Peak> result, Map<Integer, Long> fm, double limit, long charge)
  {
    for (Map.Entry<Integer, Long> i : fm.entrySet())
    {
      Integer atom_index = i.getKey();
      List<List<Peak>> sal = SAD.get(atom_index);
      Long n = i.getValue();
      int  j = 0;
      while(n > 0)
      {
        int sz = sal.size();
        if (j == sz)
        {
          // expand the list
          sal.add(new ArrayList<Peak>());
          // Make new superatom from previous largest superatom. We are trying to avoid copying on assignment here.
          convolute_basic(sal.get(j), sal.get(j - 1), sal.get(j - 1));
          prune(sal.get(j), limit);
        }
        if((n & 1) != 0L) // digit is 1, convolute result
        {
          convolute_basic(tmp , result, sal.get(j));
          prune(tmp, limit);

          // Hopefully the swap implementation will not copy all elements.
          List<Peak> intermediate = new ArrayList<Peak>(tmp);
             tmp.clear();    tmp.addAll(result);
          result.clear(); result.addAll(intermediate);
        }
        n >>>= 1;
        j++;
      }
    }

    // take charge into account
    for(Peak pk : result)
    {
      if      (charge > 0) pk.setMzAndCharge(pk.getMz() / Math.abs(charge) - ELECTRON_MASS, (int) charge);
      else if (charge < 0) pk.setMzAndCharge(pk.getMz() / Math.abs(charge) + ELECTRON_MASS, (int )charge);
    }
    return result;
  }

  /** Not protonated, neutral molecule only
   *
   * @param peptide
   * @return
   */
  public static Map<Integer, Long> newFormulaMap(Protein peptide, int charge_h)
  {
    Map<Integer, Long>          formula = new HashMap<Integer, Long>();
    Map<Element, Float> comp = peptide.getElementalComposition();
    for (Element aa : comp.keySet())
      if (comp.get(aa) > 0)
        formula = increFormulaMap(formula, EM, aa, comp.get(aa).intValue());

    // add the teminal groups
    formula = increFormulaMap(formula, EM, Element.HYDROGEN, charge_h);
    //formula = increFormulaMap(formula, EM, "O", 1);

    return formula;
  }
  public static Map<Integer, Long> increFormulaMap(Map<Integer, Long> formula, Map<String, Long> em, Element el, long atoms)
  {
    // cells the element index
    int  index = em.get(el.getSymbol()).intValue();
    formula.put(index, atoms + (formula.get(index) != null ? formula.get(index) : 0));

    return formula;
  }
  public static Map<Integer, Long> newFormulaMap(String sequence, int charge_h)
  {
    Protein peptide = new Protein(); peptide.addChain(new Peptide("id", sequence));
    return newFormulaMap(peptide, charge_h);
  }
  public static Map<Integer, Long> newFormulaMapByAveragine(double c12)
  {
    // C4.9384 H7.7583 N1.3577 O1.4773 S0.0417
    double multiple = c12 / AVERAGINE_MASS;
    Map<Integer, Long>          formula = new HashMap<Integer, Long>();

    formula = increFormulaMap(formula, EM, Element.CARBON,   Math.round(4.9384 * multiple));
    formula = increFormulaMap(formula, EM, Element.HYDROGEN, Math.round(7.7583 * multiple));
    formula = increFormulaMap(formula, EM, Element.NITROGEN, Math.round(1.3577 * multiple));
    formula = increFormulaMap(formula, EM, Element.OXYGEN,   Math.round(1.4773 * multiple));
    formula = increFormulaMap(formula, EM, Element.SULFUR,   Math.round(0.0417 * multiple));

    return formula;
  }
  public static IsoEnvelope calcIsotopesForPeptide(String peptide, int charge)
  {
    List<Peak>   result = new ArrayList<Peak>();
    Map<Integer, Long> fm = newFormulaMap(peptide, charge);

    // initialize the result
    result.add(new Peak(0.0, 1.0));

    calculate(new ArrayList<Peak>(), result, fm, 0, charge);
    return new IsoEnvelope(result, charge);
  }
  public static IsoEnvelope calcIsotopesForPeptide(Protein peptide, int charge)
  {
    List<Peak>     result = new ArrayList<Peak>();
    Map<Integer, Long> fm = newFormulaMap(peptide, charge);

    // initialize the result
    result.add(new Peak(0.0, 1.0));

    calculate(new ArrayList<Peak>(), result, fm, 0, charge);
    return new IsoEnvelope(result, charge);
  }
  public static IsoEnvelope calcIsotopesByMz(double c12, int charge, double minri, double ai)
  {
    double               limit = 0;
    List<Peak>          result = new ArrayList<Peak>(), tmp = new ArrayList<Peak>();
    Map<Integer, Long> formula = newFormulaMapByAveragine(c12 * charge);

    // initialize the result
    result.add(new Peak(0.0, 1.0));

    calculate(tmp, result, formula, limit, charge);

    // move the predictions to the c12
    double offset = c12 - result.get(0).getMz(), min_ai = result.get(0).getIntensity()*minri*0.01;
    // scale the ai to the 1st c12
    ai /= result.get(0).getIntensity();

    Iterator<Peak> itr = result.iterator();
    while (itr.hasNext())
    {
      Peak pk = itr.next();
      pk.setMzAndCharge(pk.getMz()+offset, pk.getChargeList());
      pk.setIntensity(pk.getIntensity()*ai);
      if (pk.getIntensity() < min_ai) itr.remove();
    }

    // scale the isotopes
//    double scale = c12 / result.get(0).getMz(), min_ai = result.get(0).getIntensity() * 0.001;
//    Iterator<Peak> itr = result.iterator();
//    while (itr.hasNext())
//    {
//      Peak pk = itr.next();
//      pk.setMzAndCharge(pk.getMz() * scale, pk.getChargeList());
//      if (pk.getIntensity() < min_ai) itr.remove();
//    }
    return new IsoEnvelope(result, charge);
  }
  public static IsoEnvelope subtract(List<Peak> isolation, IsoEnvelope iso, Tolerance tol, boolean purge)
  {
    boolean          verbose = false;
    double               sum = 0d, sumiso = 0d;
    Map<Peak, Peak> outcomes = new TreeMap<Peak, Peak>();

    if (verbose) System.out.print(">>" + Tools.d2s(iso.getMz(), 4) + "/" + Tools.d2s(iso.getIntensity() / 1000, 0) + "(+" + iso.getCharge() + "): ");

    // figuring out the scaling factor
    double multiplcity = 0d, c12 = 0d;
    for (int i = 0; i < iso.getPredicted().size(); i++)
    {
      if (iso.getPredicted(i).getIntensity() > 0.05) multiplcity++;
      Peak pt = new Peak(iso.getPredicted(i).getMz(), c12 == 0 ? 0d : iso.getPredicted(i).getIntensity() * c12 / iso.getPredicted(0).getIntensity());
      //MsIon pos = XYPoint_Util.findByPPM(pt, isolation, ppm * (Math.pow(i, 2) + 0.1));
      Peak pos = Peaks.find(isolation, pt, tol);
      if (pos != null)
      {
        if (verbose) System.out.print(Tools.d2s(pos.getMz(), 4) + "; ");
        outcomes.put(iso.getPredicted(i), pos);
        iso.addIsotope(pos);
        sum += pos.getIntensity(); sumiso += iso.getPredicted(i).getIntensity();
        if (i == 0)
        {
          iso.setMzAndCharge(pos.getMz(), pos.getChargeList());
          iso.setIntensity(pos.getIntensity()); c12 = pos.getIntensity();
        }
      }
      else break;
    }

    double sst = 0d, highest = Peaks.getBasePeak(isolation).getIntensity();
    Peak scale = optimum_scale(outcomes, Range.closed(0.1d, highest * 5d), 50);
    sum  /= (double )outcomes.size();
    for (Peak prdt : outcomes.keySet())
      sst += (outcomes.get(prdt).getIntensity() - sum) * (outcomes.get(prdt).getIntensity() - sum);

    iso.setScore(1 - scale.getIntensity() / sst);
    iso.setChargeScore(outcomes.size() >= iso.getCharge() ? -10d * Math.log(Math.pow(1d / multiplcity, outcomes.size())) : 0);
    if (verbose) System.out.println(", score = " + Tools.d2s(iso.getScore(), 2) + ", zscore = " + Tools.d2s(iso.getChargeScore(), 1) + " ***");

    if (!purge) return iso;

//    for (XYPoint E : outcomes.values()) E.invalidate();
/*
    // to subtract or remove the observed peaks
    double ln_fold_tol = ri_tol > 0 ? Math.abs(Math.log(ri_tol)) : -1;
    for (Map.Entry<XYPoint, XYPoint> E : outcomes.entrySet())
    {
      double expected_ri = E.getKey().getY() * scale.getX(),
                    fold = Math.abs(Math.log(E.getValue().getY() / expected_ri));
      if (fold <= ln_fold_tol)
      {
        E.getValue().setY(E.getValue().getY() - expected_ri);
      }
      else
      {
        isolation.remove(E.getValue());
      }
    } */
    return iso;
  }
  public static IsoEnvelope subtract(double[] mzs, double[] ais, Tolerance tol, IsoEnvelope iso)
  {
    boolean          verbose = false;
    double               sum = 0d, sumiso = 0d;
    Map<Peak, Peak> outcomes = new TreeMap<Peak, Peak>();

    if (verbose) System.out.print(">>" + Tools.d2s(iso.getMz(), 4) + "/" + Tools.d2s(iso.getIntensity() / 1000, 0) + "(+" + iso.getCharge() + "): ");

    // figuring out the scaling factor
    double multiplcity = 0d, c12 = 0d;
    for (int i = 0; i < iso.getPredicted().size(); i++)
    {
      if (iso.getPredicted(i).getIntensity() > 0.05) multiplcity++;
      Peak pt = new Peak(iso.getPredicted(i).getMz(), c12 == 0 ? 0d :
          iso.getPredicted(i).getIntensity() * c12 / iso.getPredicted(0).getIntensity());
      //MsIon pos = XYPoint_Util.findByPPM(pt, isolation, ppm * (Math.pow(i, 2) + 0.1));
      int pos = Peaks.find(mzs, pt.getMz(), tol);
      if (pos>= 0)
      {
        if (verbose) System.out.print(Tools.d2s(mzs[pos], 4) + "; ");
        outcomes.put(iso.getPredicted(i), new Peak(mzs[pos], ais[pos], 0));
        iso.addIsotope(new Peak(mzs[pos], ais[pos], 1));
        sum += ais[pos]; sumiso += iso.getPredicted(i).getIntensity();
        if (i == 0)
        {
          iso.setMzAndCharge(mzs[pos], iso.getChargeList());
          iso.setIntensity(ais[pos]); c12 = ais[pos];
        }
      }
      else break;
    }

    double sst = 0d, highest = ais[Peaks.getBasePeak(ais)];
    Peak scale = optimum_scale(outcomes, Range.closed(0.1d, highest * 5d), 50);
    sum  /= (double )outcomes.size();
    for (Peak prdt : outcomes.keySet())
      sst += (outcomes.get(prdt).getIntensity() - sum) * (outcomes.get(prdt).getIntensity() - sum);

    iso.setScore(1 - scale.getIntensity() / sst);
    iso.setChargeScore(outcomes.size() >= iso.getCharge() ? -10d * Math.log(Math.pow(1d / multiplcity, outcomes.size())) : 0);
    if (verbose) System.out.println(
        ", score = "+Tools.d2s(iso.getScore(),2)+", zscore = "+Tools.d2s(iso.getChargeScore(),1)+" ***");

    return iso;
  }
  public static IsoEnvelope subtract(double[] mzs, double[] ais, Tolerance tol, String peptide, int charge)
  {
    return subtract(mzs, ais, tol, calcIsotopesForPeptide(peptide, charge));
  }
  public static IsoEnvelope subtract(double[] mzs, double[] ais, Tolerance tol, Protein peptide, int charge)
  {
    return subtract(mzs, ais, tol, calcIsotopesForPeptide(peptide, charge));
  }
  public static IsoEnvelope subtract(List<Peak> isolation, Tolerance tol, double ri_tol,String peptide, int charge)
  {
    return subtract(isolation, calcIsotopesForPeptide(peptide, charge), tol, true);
  }
  public static IsoEnvelope subtract(List<Peak> isolation, Tolerance tol, double ri_tol, Protein peptide, int charge)
  {
    return subtract(isolation, calcIsotopesForPeptide(peptide, charge), tol, true);
  }
  protected static Peak optimum_scale(Map<Peak, Peak> outcomes, Range<Double> bound, int steps)
  {
    boolean verbose = false;
    double delta = (bound.upperEndpoint() - bound.lowerEndpoint()) / (double )steps,
          delta2 = delta / (double )steps, ss2 = 0d, best = Double.MAX_VALUE;
    Peak optimum = new Peak(0d, 0d);
    for (double x = bound.lowerEndpoint(); x <= bound.upperEndpoint(); x += delta)
    {
      ss2 = 0d;
      for (Map.Entry<Peak, Peak> E : outcomes.entrySet())
        ss2 += (E.getKey().getIntensity() * x - E.getValue().getIntensity()) *
               (E.getKey().getIntensity() * x - E.getValue().getIntensity());

      if (ss2 < best)
      {
        best = ss2; optimum.setMzAndCharge(x); optimum.setIntensity(best);
        if    (verbose) System.out.print(" **");
      }
      else if (verbose) System.out.print("   ");

      if (verbose) System.out.println(Tools.d2s(x, 1)+" --> "+ss2+": "+Tools.d2s(optimum.getMz(), 1)+" <-- "+best);
    }
    if (best < Double.MAX_VALUE)
      for (double x = optimum.getMz() - delta; x <= optimum.getMz() + delta; x += delta2)
      {
        ss2 = 0d;
        for (Map.Entry<Peak, Peak> E : outcomes.entrySet())
          ss2 += (E.getKey().getIntensity() * x - E.getValue().getIntensity()) *
                 (E.getKey().getIntensity() * x - E.getValue().getIntensity());

        if (ss2 < best)
        {
          best = ss2; optimum.setMzAndCharge(x); optimum.setIntensity(best);
          if    (verbose) System.out.print(" $$");
        }
        else if (verbose) System.out.print("   ");

        if (verbose) System.out.println(Tools.d2s(x,1)+" --> "+ss2+": "+Tools.d2s(optimum.getMz(), 1)+" <-- "+best);
      }

    return optimum;
  }
  public static boolean isEmbedd(IsoEnvelope A, Tolerance tol, Collection<IsoEnvelope> others)
  {
    if (!Tools.isSet(others) || A == null) return false;
    for (IsoEnvelope other : others)
      if (other != null && other.contains(A, tol)) return true;

    return false;
  }
//  public static <T extends Peak> List<IsoEnvelope> decodePrecursors(double[] mzs, double[] ais, int z, Tolerance tol, int max_z, Protein... peptides)
//  {
//    List<Peak>           peaks = new ArrayList<Peak>(isolation.getRawData());
//    List<IsoEnvelope> precursors = new ArrayList<IsoEnvelope>();
//
//    // remove the peaks corresponding to the peptides if required
//    if (Tools.isSet(peptides))
//      for (Protein peptide : peptides)
//      {
//        if (peptide == null) continue;
//        // force a component for each peptide
//        IsoEnvelope iso = subtract(mzs, ais, tol, peptide, z);
//        precursors.add(iso);
//      }
//
//    Collections.sort(isolation.getRawData());
//    Collections.sort(peaks, MsIon.INTENSITY_DESCEND);
//
//    for (MsIon ion : peaks)
//    {
//      // treating every peak as a possible c12
//      MsIsotopes best = null;
//      for (int z = 1; z <= max_z; z++)
//      {
//        MsIsotopes expected = Isotope_Util.subtract(isolation.getRawData(), Isotope_Util.calcIsotopesByMz(ion.getMz(), z), ppm, false);
//        if (best == null || expected.getChargeScore() > best.getChargeScore()) best = expected;
//      }
//      if (best != null && !Isotope_Util.isEmbedd(best, ppm, precursors) && best.getScore() > 0)
//      {
//        best.getIsotopes().clear();
//        precursors.add(Isotope_Util.subtract(isolation.getRawData(), best, ppm, true));
//      }
//    }
//    // figure out the noise level
//    double baseline = 0d;
//    Collection<T> remainders = new ArrayList<T>();
//    for (T pt : (List<T> )isolation.getRawData())
//      if (pt.isValid()) { remainders.add(pt); baseline += pt.getY(); }
//
//    if (remainders.size() > 2) baseline /= (double )remainders.size(); else baseline = 0d;
//
//    // reconcile the minimum representation of the precursors
//    Iterator<MsIsotopes> itr = precursors.iterator();
//    while (itr.hasNext())
//    {
//      MsIsotopes item = itr.next();
//      if (item.getY() < baseline * 1.64d) itr.remove(); else System.out.println(item.toString());
//    }
//
//    Collections.sort(precursors, XYPoint.Y_DESCEND);
//
//    return precursors;
//  }
}
