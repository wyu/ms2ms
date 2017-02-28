package org.ms2ms.algo;

import com.google.common.collect.*;
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
  public static final double C   = 12.00;
  public static final double H   = 1.007825;
  public static final double O   = 15.994915;
  public static final double N   = 14.003074;
  public static final double S   = 31.972072;
  public static final double H2O = 2d*H+O;
  public static final double NH3 = N+3d*H;
  public static final double CO  = C+O;
  public static final double OH  = H+O;

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
  public static TreeMultimap<Float, String> toMassAAs2x(Map<Character, Float> AAs)
  {
    if (Tools.isSet(AAs))
    {
      TreeMultimap<Float, String> m = TreeMultimap.create();
      for (Character AA : AAs.keySet())
      {
        if (sAA.indexOf(AA)<0) continue;
//        m.put(AAs.get(AA), AA+""); // take care of the single residue first
        for (Character aa : AAs.keySet())
        {
          if (sAA.indexOf(aa)<0) continue;
          String newAA = AA.toString()+aa.toString();
          if (!m.containsValue(newAA)) m.put(AAs.get(AA)+AAs.get(aa), newAA);
        }
      }
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
    if ("tmt10".equalsIgnoreCase(fixed) || "tmt6".equalsIgnoreCase(fixed))
    {
      TreeMap<Character, Float> map = new TreeMap<>(AAsBuilder.build());
      map.put('K', 357.25789f); map.put('^', 229.162932f); // TMT-10

      return ImmutableMap.copyOf(map);
    }
    return ImmutableMap.copyOf(new TreeMap<>(AAsBuilder.build()));
  }
  // calculate the MH value of the peptide defined by the positions (left, right, inclusive) onto the 'sequence'
  public static double calcMH(char[] sequence, int left, int right, float[] AAs)
  {
//    double y=AAs.get('$')+2* 1.00783;
    double y=AAs['$']+AAs['^']+2*1.007825;
    for (int c=right; c>=left; c--)
      if (c>=0 && c<sequence.length)
        y+=AAs[sequence[c]];
      else
        System.out.print("");

    return y;
  }
  public static double calcMH(char[] sequence, int left, int right, Map<Character, Float> AAs)
  {
//    double y=AAs.get('$')+2* 1.00783;
    double y=AAs.get('$')+AAs.get('^')+2*1.007825;
    for (int c=right; c>=left; c--) y+=AAs.get(sequence[c]);

    return y;
  }
  public static int numTryptic(String peptide)
  {
    // check the uncleaved site
    int sites=0;
    for (int i=0; i<peptide.length()-1; i++)
      if (isTryptic(peptide.charAt(i), peptide.charAt(i+1))) sites++;

    return sites;
  }
  public static boolean isTryptic(char[] sequence, int n0, int n1)
  {
    // not checking for sequence validity to save time. Do it before call this function!

    // should the start or end of the sequence or conforming to the tryptic specificity
    return n0==0 || n0==sequence.length-1 ||
        (n0>=0 && n1<sequence.length && (sequence[n0]=='K' || sequence[n0]=='R') && (sequence[n1]!='P'));
  }
  public static boolean isTryptic(char n0, char n1)  { return (n0=='K' || n0=='R') && n1!='P'; }

  public static List<Integer> seekRemoval(String sequence, boolean reverse, double calcM, double deltaM, Tolerance tol, Range<Integer> isoErr, float[] AAs)
  {
    if (!Strs.isSet(sequence) || "-".equals(sequence)) return null;

    List<Integer> removed = new ArrayList<>();

    float ct=(float )deltaM;
    int start=reverse?sequence.length()-1:0, stop=reverse?0:sequence.length()-1,step=reverse?-1:1, i=start;
    while (1==1)
    {
      ct-=AAs[sequence.charAt(i)]*step; removed.add(i);
      for (int iso=isoErr.lowerEndpoint(); iso<=isoErr.upperEndpoint(); iso++)
      {
        if (tol.withinTolerance(calcM+iso*Isotopes.DELTA_C13, ct+calcM)) return removed;
      }
      i+=step;
      if (reverse?(i<0):(i>stop)) break;
    }
    return null;
  }
  public static Map<String, Double> newCtStarters(String tag)
  {
    Map<String, Double> Ct = new HashMap<>();

    Multimap<Float, Character> AAs = Peptides.newMassAAs(tag);
    for (Float f : AAs.keySet())
      Ct.put(Tools.front(AAs.get(f))+"", f+Peptides.H+Peptides.H2O);

    return Ct;
  }
  // generate the predicted fragments of the peptide.
  public static Map<Float, String> toFragments(char[] peptide, Map<Integer, Double> mods, float[] AAs)
  {
    if (peptide==null || peptide.length==0) return null;

//    Float H=1.007825f, O=15.994915f, C=12f, CO=C+O, OH=O+H, N=14.003074f, NH3=N+3*H, H2O=2*H+O,
//        b=(AAs['^']+1.007825f), y=(AAs['$']+2*1.007825f), a=b-CO,z=y-NH3;
    double b=(AAs['^']+H), y=(AAs['$']+2*H), a=b-CO,z=y-NH3;

    Map<Float, String> frags = new TreeMap<>(); int y_18=0, y_17=0, b_18=0, b_17=0;
    for (int i=0; i<peptide.length; i++)
    {
      int j=peptide.length-i-1;
      y+=AAs[peptide[j]]+(mods!=null&&mods.get(j)!=null?mods.get(j).floatValue():0f);
      b+=AAs[peptide[i]]+(mods!=null&&mods.get(i)!=null?mods.get(i).floatValue():0f);

      if ("RKQN".indexOf(peptide[j])>=0) y_17++;
      if ("STED".indexOf(peptide[j])>=0) y_18++;
      if ("RKQN".indexOf(peptide[i])>=0) b_17++;
      if ("STED".indexOf(peptide[i])>=0) b_18++;

      frags.put((float )y,      "y"+(i+1));
      frags.put((float)(y-NH3), "z"+(i+1));
      frags.put((float )b,      "b"+(i+1));
      frags.put((float )(b-CO), "a"+(i+1));

      // any neutral loss?
      if (y_17>0) frags.put((float )(y-NH3), "y"+(i+1)+"-17");
      if (b_18>0) frags.put((float )(b-H2O), "b"+(i+1)+"-18");

      // get the internal ions
      if (i>0)
      {
        Float internal=(float )H, intn28=null;
        for (int k=2; k<5; k++)
        {
          if (i+k>=peptide.length) break;
          // calculate the expected mass
          internal+=AAs[peptide[k+i]]; intn28=(float )(internal-CO);
          if (!frags.containsKey(internal)) frags.put(internal, Strs.toString(peptide, i,i+k));
          if (!frags.containsKey(intn28))   frags.put(intn28,   Strs.toString(peptide, i,i+k)+"-28");
        }
      }
    }
    return frags;
  }
  // generate the predicted fragments of the peptide.
  public static float[] toYs(char[] peptide, Map<Integer, Double> mods, float[] AAs)
  {
    if (peptide==null || peptide.length==0) return null;

    Float   y  = (AAs['$']+2*1.007825f);
    float[] ys = new float[peptide.length];
    for (int i=0; i<peptide.length; i++)
    {
      int j=peptide.length-i-1;
      y+=AAs[peptide[j]]+(mods!=null&&mods.get(j)!=null?mods.get(j).floatValue():0f);
      ys[i] = y;
    }
    return ys;
  }
  // generate the predicted fragments of the peptide.
  public static int[] toZYs(char[] peptide)
  {
    if (peptide==null || peptide.length==0) return null;

    float z=0; int[] ys = new int[peptide.length];
    for (int i=0; i<peptide.length; i++)
    {
      z    += ("KR".indexOf(peptide[peptide.length-i-1])>=0?1:("H".indexOf(peptide[i])>=0?0.5:0));
      ys[i] = Math.round(z);
    }
    return ys;
  }
  // generate the predicted fragments of the peptide.
  public static float[] toBs(char[] peptide, Map<Integer, Double> mods, float[] AAs)
  {
    if (peptide==null || peptide.length==0) return null;

    float   b  = (AAs['^']+1.007825f);
    float[] bs = new float[peptide.length];
    for (int i=0; i<peptide.length; i++)
    {
      b+=AAs[peptide[i]]+(mods!=null&&mods.get(i)!=null?mods.get(i).floatValue():0f);
      bs[i] = b;
    }
    return bs;
  }
  public static int[] toZBs(char[] peptide)
  {
    if (peptide==null || peptide.length==0) return null;

    float z=0.5f; int[] bs = new int[peptide.length];
    for (int i=0; i<peptide.length; i++)
    {
      z    += ("KR".indexOf(peptide[i])>=0?1:("H".indexOf(peptide[i])>=0?0.5:0));
      bs[i] = Math.round(z);
    }
    return bs;
  }
}
