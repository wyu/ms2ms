package org.ms2ms.algo;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.hfg.bio.Element;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.ms2ms.data.collect.ImmutableNavigableMap;
import org.ms2ms.data.ms.IsoEnvelope;
import org.ms2ms.data.ms.OffsetPpmTolerance;
import org.ms2ms.data.ms.PeakMatch;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.*;

/** Non-static version of Isotopes
 *
 * Created by yuw on 3/6/17.
 */
public class Isotopics
{
  // initialization methods
  final double   ELECTRON_MASS = 0.00054858;
  final double      DUMMY_MASS = -10000000;
  final double  AVERAGINE_MASS = 111.0543052;
  final double       DELTA_C13 = 1.00335d;

  Map<String, Long>         EM = new HashMap<String, Long>();
  List<List<List<Peak>>>   SAD = new ArrayList<List<List<Peak>>>();

  Table<Integer, Integer, List<Peak>> mMassIso = HashBasedTable.create();
//  public static Histogram dp_prediction = new Histogram();

  public Isotopics()
  {
    try
    {
      init_data(Isotopes.class.getClassLoader().getResourceAsStream("org/ms2ms/algo/ISOTOPE.DAT"), SAD, EM);
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }

  public List<Peak> predictEnvelope(double c12, int charge, double minri)
  {
    Integer      mass = (int )c12 * charge;
    List<Peak> result = mMassIso.get(charge, mass);

    if (result==null)
    {
      // initialize the result
      result = new ArrayList<>();
      result.add(new Peak(0.0, 1.0));
      result = calculate(new ArrayList<>(), result, newFormulaMapByAveragine(mass), 0, charge);

      Iterator<Peak> itr = result.iterator();
      while (itr.hasNext())
        if (itr.next().getIntensity() < minri) itr.remove();

      mMassIso.put(charge, mass, result);
    }
//    // a nasty bug, WYU 20170625
//    else if (result.get(0).getCharge()!=charge)
//    {
//      // make sure the charge states are right
//      List<Peak> results = new ArrayList<>(result.size());
//      for (Peak p : result)
//        results.add(new Peak(Peaks.MnH2MnH(p.getMz(), p.getCharge(), charge), p.getIntensity(), charge));
//
//      return results;
//    }

    return result;
  }
  public IsoEnvelope calcIsotopesByMz(double c12, int charge, double minri, double ai)
  {
    List<Peak> result = predictEnvelope(c12, charge, minri);

    // move the predictions to the c12
    double offset = c12 - result.get(0).getMz();
    // scale the ai to the 1st c12
    ai /= result.get(0).getIntensity();

    List<Peak> isotopes = new ArrayList<>();
    for (Peak pk : result)
      isotopes.add(new Peak(pk.getMz()+offset, ai*pk.getIntensity(), pk.getChargeList()));

    IsoEnvelope iso = new IsoEnvelope(isotopes, charge);
    isotopes=(List )Tools.dispose(isotopes);
    return iso;
  }

  public List<Peak> calculate(List<Peak> tmp, List<Peak> result, Map<Integer, Long> fm, double limit, long charge)
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

  // Merge two patterns to one.
  public void convolute_basic(List<Peak> h, final List<Peak> g, final List<Peak> f)
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
  public void prune(List<Peak> f, double limit)
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
  public void print_pattern(List<Peak> result, int digits)
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

  public Map<Integer, Long> newFormulaMapByAveragine(double c12)
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

  public Map<Integer, Long> increFormulaMap(Map<Integer, Long> formula, Map<String, Long> em, Element el, long atoms)
  {
    // cells the element index
    int  index = em.get(el.getSymbol()).intValue();
    formula.put(index, atoms + (formula.get(index) != null ? formula.get(index) : 0));

    return formula;
  }

  private boolean init_data(InputStream dat, List<List<List<Peak>>> sad, Map<String, Long> em) throws IOException
  {
    BufferedReader reader = new BufferedReader(new InputStreamReader(dat));

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

//  public IsoEnvelope gatherIsotopes(ImmutableNavigableMap<PeakMatch> peaks, double mz, double mh, OffsetPpmTolerance tol)
//  {
//    // return null, if has negative intensity
//    IsoEnvelope iso = PeakMatch.query4isotope(peaks, mz, tol, false);
//
//    // for large peptide, also check the 1st c13
//    if (iso==null && mz>1000)
//    {
//      iso = PeakMatch.query4isotope(peaks, mz+Peptides.C13, tol, false);
//      if (iso!=null)
//        iso = new IsoEnvelope(iso.getMz()-Peptides.C13, iso.getIntensity(), 0).setScore(iso.getScore());
//    }
//
//    // quit if not matching to the first mz!
//    if (iso!=null)
//    {
//      iso.setMz(iso.getMz()).setCharge((int) Math.round(mh/iso.getMz()));
//      // check for any below m/z
//      if (PeakMatch.query4counts(peaks, tol, mz-(1.0025/(double) iso.getCharge()))>0)
//      {
//        // got an incorrect c12!
//        iso.setChargeScore(-1);
//      }
//      else
//      {
//        iso.setChargeScore(1d);
//        // looking forward
//        for (int i=1; i<iso.getCharge()+1; i++)
//        {
//          double ai = PeakMatch.query4ai(peaks, mz+(i*1.0025/(double) iso.getCharge()), tol);
//          if (ai<=0) break;
//
//          iso.setChargeScore(iso.getChargeScore()+1d);
//          iso.setIntensity(iso.getIntensity()+ai);
//        }
//      }
//
//      return iso;
//    }
//
//    return null;
//  }

  public IsoEnvelope gatherIsotopes(ImmutableNavigableMap<PeakMatch> peaks, double mz, double mh, OffsetPpmTolerance tol)
  {
    // return null, if has negative intensity
    IsoEnvelope iso = PeakMatch.query4isotope(peaks, mz, tol, false);

    // for large peptide, also check the 1st c13
    iso = gatherFirstIsotope(iso, peaks, mz, tol);
    return gather(iso, peaks, mz, mh, tol);
  }

  public IsoEnvelope gatherFirstIsotope(IsoEnvelope iso, ImmutableNavigableMap<PeakMatch> peaks, double mz, OffsetPpmTolerance tol)
  {
    // for large peptide, also check the 1st c13
    if (iso==null && mz>1000)
    {
      iso = PeakMatch.query4isotope(peaks, mz+Peptides.C13, tol, false);
      if (iso!=null)
        iso = new IsoEnvelope(iso.getMz()-Peptides.C13, iso.getIntensity(), 0).setScore(iso.getScore());
    }
    return iso;
  }
  public IsoEnvelope gather(IsoEnvelope iso, ImmutableNavigableMap<PeakMatch> peaks, double mz, double mh, OffsetPpmTolerance tol)
  {
    // quit if not matching to the first mz!
    if (iso!=null)
    {
      iso.setMz(iso.getMz()).setCharge((int) Math.round(mh/iso.getMz()));
      // check for any below m/z
      if (PeakMatch.query4counts(peaks, tol, mz-(1.0025/(double) iso.getCharge()))>0)
      {
        // got an incorrect c12!
        iso.setChargeScore(-1);
      }
      else
      {
        iso.setChargeScore(1d);
        // looking forward
        for (int i=1; i<iso.getCharge()+1; i++)
        {
          double ai = PeakMatch.query4ai(peaks, mz+(i*1.0025/(double) iso.getCharge()), tol);
          if (ai<=0) break;

          iso.setChargeScore(iso.getChargeScore()+1d);
          iso.setIntensity(iso.getIntensity()+ai);
        }
      }

      return iso;
    }
    return null;
  }

//  // deprecated
//  public IsoEnvelope gatherIsotopes(ImmutableNavigableMap<PeakMatch> peaks, double mz, double mh, OffsetPpmTolerance tol)
//  {
//    int[]         mz_ai = peaks.query(tol.toActualBoundary(mz));
//    PeakMatch[] matched = peaks.fetchVals(mz_ai);
//
//    // quit if not matching to the first mz!
//    if (matched!=null && !PeakMatch.hasNegativeIntensity(matched))
//    {
//      IsoEnvelope iso = new IsoEnvelope();
//      iso.setMzAndCharge(Stats.mean0(peaks.fetchKeys(mz_ai)), 0);
//      iso.setMzAndCharge(iso.getMz(), (int) Math.round(mh/iso.getMz()));
//      iso.setIntensity(PeakMatch.AbsIntensitySum(matched));
//      iso.setScore(matched[0].getFrequency());
//
//      // check for any below m/z
//      double[] bound = tol.toActualBoundary(mz-(1.0025/(double) iso.getCharge()));
//      mz_ai=peaks.query(bound);
//      if (mz_ai!=null)
//      {
//        // got an incorrect c12!
//        iso.setChargeScore(-1);
//      }
//      else
//      {
//        iso.setChargeScore(1d);
//        // looking forward
//        for (int i=1; i<iso.getCharge()+1; i++)
//        {
//          bound = tol.toActualBoundary(mz+(i*1.0025/(double) iso.getCharge()));
//          mz_ai = peaks.query(bound);
//          if (mz_ai==null) break;
//
//          matched = peaks.fetchVals(mz_ai);
//          iso.setChargeScore(iso.getChargeScore()+1d);
//          iso.setIntensity(iso.getIntensity()+PeakMatch.AbsIntensitySum(matched));
//          matched=null; mz_ai=null; bound=null;
//        }
//      }
//      bound=null; mz_ai=null; matched=null;
//
//      return iso;
//    }
//    mz_ai=null; matched=null;
//
//    return null;
//  }
}
