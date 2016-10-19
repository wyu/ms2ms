package org.ms2ms.algo;

import com.google.common.collect.*;
import org.apache.commons.collections.map.HashedMap;
import org.expasy.mzjava.core.mol.SymbolSequence;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.Tolerance;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodManager;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodMod;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmentAnnotator;
import org.expasy.mzjava.proteomics.ms.fragment.PeptideFragmenter;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

/**
 * Created by yuw on 9/14/2015.
 */
public class Peptides
{
  static final ImmutableMap.Builder<Character, Float> AAsBuilder =
      new ImmutableMap.Builder<Character, Float>()
          .put('G',57.02146f  ).put('A', 71.03711f).put('S', 87.03203f).put('P', 97.05276f).put('V', 99.06841f)
          .put('T', 101.04768f).put('I', 113.08406f).put('L', 113.08406f).put('N', 114.04293f).put('D', 115.02694f)
          .put('Q', 128.05858f).put('E', 129.04259f).put('M', 131.04049f).put('H', 137.05891f).put('F', 147.06841f)
          .put('R', 156.10111f).put('C', 160.03065f).put('Y', 163.06333f).put('W', 186.07931f).put('K', 128.09496f)
          .put('U', 150.95363f).put('^', 1.00783f).put('$', 17.00273f); // N/C-terminal mod

  public static final String sAA = "ACDEFGHIKLMNPQRSTVWY";
  public static <P extends SymbolSequence> Multimap<String, P> toSequenceMap(Collection<P>... proteins)
  {
    if (!Tools.isSet(proteins)) return null;
    // create a map
    Multimap<String, P> seq_prots = HashMultimap.create();
    for (Collection<P> ps : proteins)
      for (P p : ps) seq_prots.put(p.toSymbolString(), p);

    return seq_prots;
  }
  public static PeptideFragmentAnnotator newHCDAnnotator()
  {
    // the fragmenter is responsible for generating a theoretical spectrum from the peptide
    PeptideFragmenter fragmenter =
      new PeptideFragmenter(EnumSet.of(IonType.b, IonType.y, IonType.a), PeakList.Precision.DOUBLE);

    // the annotator needs to delegate to the fragmenter the fragmentation process and to
    // the internal aligner the tolerance for aligning peaks
    return new PeptideFragmentAnnotator(fragmenter, new AbsoluteTolerance(0.1));
  }
  public static String keepAAs(String seq)
  {
    String peptide=null;
    if (Strs.isSet(seq))
    {
      for (char c : seq.toCharArray())
        if (sAA.indexOf(c)>=0) peptide=Strs.extend(peptide, c+"","");
    }
    return peptide;
  }
  public static Multimap<Float, Character> newMassAAs(String fixed)
  {
    return toMassAAs(newAAsMass(fixed));
  }
  public static TreeMultimap<Float, Character> toMassAAs(Map<Character, Float> AAs)
  {
    if (Tools.isSet(AAs))
    {
      TreeMultimap<Float, Character> m = TreeMultimap.create();
      for (Character AA : AAs.keySet()) m.put(AAs.get(AA), AA);
      return m;
    }
    return null;
  }
  public static Map<String, Float> modAAsMass(String... unimods)
  {
    Map<String, Float> AAs = new HashMap<>(); Map<Character, Float> As = AAsBuilder.build();
    for (Character c : As.keySet()) AAs.put(c+"", As.get(c));

    Map<String, String[]> mod_site = new HashMap<>();
    for (String unimod : unimods)
    {
      String[] strs = unimod.split(":");
      mod_site.put(strs[0], strs[1].split("/"));
    }
    UnimodManager manager = UnimodManager.getInstance();
    List<UnimodMod>    mm = manager.getModificationList();
    for (UnimodMod m : mm)
      for (String unimod : mod_site.keySet())
        if (Strs.isA(unimod, m.getFullName(), m.getInterimName(), m.getPsiMsMsName(), m.getLabel()))
          for (String site : mod_site.get(unimod))
            if (AAs.containsKey(site))
              AAs.put(site+"^"+m.getLabel(), (float )(m.getMolecularMass()+AAs.get(site)));

    return AAs;
  }
  public static Map<String, Float> x2ModAAsMass(Map<String, Float> As)
  {
    As.remove("^"); As.remove("$");

    Map<String, Float> added = new HashMap<>();
    for (String x1 : As.keySet())
      for (String x2 : As.keySet())
      {
        Float sum = As.get(x1)+As.get(x2);
        if (!As.containsKey(sum) && !added.containsValue(sum)) added.put(x1+"-"+x2, sum);
      }

    // combine the entries
    added.putAll(As);

    return added;
  }
  // setup a simple mapping of the AA symbol and their incremental masses
  public static ImmutableMap<Character, Float> newAAsMass(String fixed)
  {
//    Map<Character, Double> AAs = new HashMap<>();
//    AAs.put('G', 57.02146d );AAs.put('A', 71.03711d );AAs.put('S', 87.03203d );AAs.put('P', 97.05276d );AAs.put('V', 99.06841d);
//    AAs.put('T', 101.04768d);AAs.put('I', 113.08406d);AAs.put('L', 113.08406d);AAs.put('N', 114.04293d);AAs.put('D', 115.02694d);
//    AAs.put('Q', 128.05858d);AAs.put('E', 129.04259d);AAs.put('M', 131.04049d);AAs.put('H', 137.05891d);AAs.put('F', 147.06841d);
//    AAs.put('R', 156.10111d);AAs.put('C', 160.03065d);AAs.put('Y', 163.06333d);AAs.put('W', 186.07931d);AAs.put('K', 128.09496d); AAs.put('U',150.95363d);
//    AAs.put('^', 1.00783d  ); // N-terminal mod
//    AAs.put('$', 17.00273d ); // C-terminal mod

    if ("tmt10".equalsIgnoreCase(fixed) || "tmt6".equalsIgnoreCase(fixed))
    {
      TreeMap<Character, Float> map = new TreeMap<>(AAsBuilder.build());
      map.put('K', 357.25789f); map.put('^', 229.162932f); // TMT-10

      return ImmutableMap.copyOf(map);
//      AAs.put('K', 357.25789d);
//      AAs.put('^', 229.162932d); // TMT-10
    }
    return null;
  }
  // calculate the MH value of the peptide defined by the positions (left, right, inclusive) onto the 'sequence'
  public static double calcMH(char[] sequence, int left, int right, float[] AAs)
  {
//    double y=AAs.get('$')+2* 1.00783;
    double y=AAs['$']+AAs['^']+2*1.007825;
    for (int c=right; c>=left; c--) y+=AAs[sequence[c]];

    return y;
  }
  public static double calcMH(char[] sequence, int left, int right, Map<Character, Float> AAs)
  {
//    double y=AAs.get('$')+2* 1.00783;
    double y=AAs.get('$')+AAs.get('^')+2*1.007825;
    for (int c=right; c>=left; c--) y+=AAs.get(sequence[c]);

    return y;
  }
  public static boolean isTryptic(char[] sequence, int n0, int n1)
  {
    // not checking for sequence validity to save time. Do it before call this function!

    // should the start or end of the sequence or conforming to the tryptic specificity
    return n0==0 || n0==sequence.length-1 ||
        (n0>=0 && n1<sequence.length && (sequence[n0]=='K' || sequence[n0]=='R') && (sequence[n1]!='P'));
  }
  public static boolean isTryptic(char n0, char n1)  { return (n0=='K' || n0=='R') && n1!='P'; }

  public static Collection<Integer> seekRemoval(String sequence, boolean reverse, double calcM, double deltaM, Tolerance tol, Range<Integer> isoErr, float[] AAs)
  {
    if (!Strs.isSet(sequence) || "-".equals(sequence)) return null;

    List<Integer> removed = new ArrayList<>();

    float ct=(float )deltaM;
    int start=reverse?sequence.length()-1:0, stop=reverse?0:sequence.length()-1,step=reverse?-1:1;
    for (int i=start; i<=stop; i+=step)
    {
      ct+=AAs[sequence.charAt(i)]; removed.add(i);
      for (int iso=isoErr.lowerEndpoint(); iso<=isoErr.upperEndpoint(); iso++)
      {
        if (tol.withinTolerance(calcM+iso*Isotopes.DELTA_C13, ct+calcM)) return removed;
      }
    }
    return null;
  }
}
