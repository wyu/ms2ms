package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.ms2ms.algo.Peptides;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;

import java.util.*;

// holder of AA definition and modifications
//
public class ResidueBase implements Cloneable
{
  public static ResidueBase REDUCED, CIM, TMT10;
  public static Map<Character, Double> sVarCharMods=null;

  // user-defined residue bases
  public Map<String, ResidueBase> mUserBases;

  static
  {
    sVarCharMods = new HashMap<>();
    sVarCharMods.put('N',  0.984016d); // deamidation
//    sVarMods.put('Q',  0.984016d); // deamidation at a much slower rate
    sVarCharMods.put('M', 15.994915d);
    sVarCharMods.put('^', -229.162932d); // for testing only
    sVarCharMods.put('c',-17.026549d); // N-terminal cyclization
    sVarCharMods.put('q',-17.026549d); // N-terminal cyclization

    REDUCED = new ResidueBase(); REDUCED.init("reduced");
    CIM     = new ResidueBase();     CIM.init("LFQ");
    TMT10   = new ResidueBase();   TMT10.initTMT();
  }

  private Map<Character, Float> mAADic;
  private Map<String, Float> mDenovoModAAs;

  private float[] mAAs;

  private NumericalModVault mModVault;

  private Map<Double, Character> mAACt =null;
  private double[] mTrivialMods =null;
  private ImmutableSetMultimap<Character, Double> mVarAAMods4Exact=null;

//  private Multimap<String, Double> sVarMods=null;
  private MultiTreeTable<String, String, Double> sExactMods=null;
  private Map<Double, Character> sAACt =null;
  private SortedMap<Double, Float> sModCharge=null;
  private SortedMap<Double, Float> sFacileMods=null;

  // pull out the variable mods at the N/C-term
  private Collection<Double> mNTmods, mCTmods;


  public ResidueBase() { super(); }

  public static ResidueBase valueOf(String s)
  {
    if (Strs.isA(Strs.toupper(s),"TMT10","TMT6")) return TMT10;
    if (Strs.isA(Strs.toupper(s),"REDUCED"))      return REDUCED;

    return CIM;
  }
  public Map<Character, Float> getAADic()          { return mAADic; }
  public Map<String, Float>    getDenovoModAAs()
  {
    if (mDenovoModAAs!=null && mDenovoModAAs.containsKey("I"))
    {
      // no need to keep multiple AA of the same mass
      mDenovoModAAs.remove("L");
      mDenovoModAAs.remove("J");
    }
    return mDenovoModAAs;
  }
  public float[]               getAAs()            { return mAAs; }
  public double[]              getTrivialMods()    { return mTrivialMods; }
//  public Range<Double>         getModBound()       { return mModBound; }
  public NumericalModVault     getModVault()       { return mModVault; }
  public SortedMap<Double, Float> getModCharge()   { return sModCharge; }
  public SortedMap<Double, Float> getFacileMods()  { return sFacileMods; }
  public MultiTreeTable<String, String, Double> getExactMods() { return sExactMods; }
  public Map<Double, Character> getAACt()          { return mAACt; }
  public ImmutableSetMultimap<Character, Double> getVarAAMods4Exact() { return mVarAAMods4Exact; }
  public Collection<Double> getNTmods() { return mNTmods; }
  public Collection<Double> getCTmods() { return mCTmods; }

  public ResidueBase setAAs(Map<Character, Float> aas)
  {
    mAAs = Peptides.toAAs(aas);
    return this;
  }
  public ResidueBase setAAs(float[] aas)
  {
    mAAs = aas;
    return this;
  }
  public ResidueBase setModVault(             String s) { mModVault=new NumericalModVault(s); return this; }
  public ResidueBase setModVault(  NumericalModVault s) { mModVault     =s; return this; }
  public ResidueBase setAADic( Map<Character, Float> s) { mAADic=s; return this; }

  public void config(String cfg_file)
  {

  }

  private void init(String fixed)
  {
    mAADic = Peptides.newAAsMass(fixed);
    mAAs   = Peptides.toAAs(mAADic);

    mAACt = Maps.newHashMap();
    mAACt.put(175.1190d, 'R');
    mAACt.put(147.1128d, 'K');

    // set up an array for quick lookup
//    mVarAAMods4Exact['K'] = 27.994915d;
    ImmutableSetMultimap.Builder<Character, Double> builder = ImmutableSetMultimap.builder();
    builder.put('N', Peptides.N2D);
    builder.put('M', Peptides.O);
    builder.put('W', Peptides.O*2d);
    mVarAAMods4Exact= builder.build();

    sFacileMods=Maps.newTreeMap();
    //sVarMods.put( "*", -0.984016d); // C-terminal, any
    sFacileMods.put(162.052824d, 75f); // Hexose
    sFacileMods.put(192.063388d,  0f); // Heptose
    sFacileMods.put(79.966331d,  75f); // Phospho

//    // move the static vars here
//    sVarMods= HashMultimap.create();
//    //sVarMods.put( "*", -0.984016d); // C-terminal, any
//    sVarMods.put( "*",   1.003355d); // c12->c13
//    sVarMods.put( "*",   2.006710d); // c12->2xc13
//    sVarMods.put( "N",   0.984016d); // deamidation
//    sVarMods.put( "M",  15.994915d);
//    sVarMods.put( "K",  27.994915d); // Formylation @ Nt and K
//    sVarMods.put("^C", -17.026549d); // N-terminal cyclization
//    sVarMods.put("^Q", -17.026549d); // N-terminal cyclization
//    sVarMods.put("^S", -18.01056d);  // N-terminal cyclization

    sModCharge=Maps.newTreeMap();

    sAACt =Maps.newHashMap();
    sAACt.put(175.1190d, 'R');
    sAACt.put(147.1128d, 'K');

    mTrivialMods = new double[] {Peptides.C13, Peptides.C13*2};

    // mods that ought to be placed in the 'Exact' matches because their are expected at a particular residue type
    sExactMods = MultiTreeTable.create();
    sExactMods.put("*", "M", 15.994915d);
    sExactMods.put("*", "N", 0.984016d);
    sExactMods.put("^", "Q", -17.026549d);
    sExactMods.put("^", "C", -17.026549d);

    mDenovoModAAs = Peptides.modAAsMass("Oxidation:M");
    mDenovoModAAs.remove("^"); mDenovoModAAs.remove("$"); mDenovoModAAs.remove("K");

    // pull out the variable mods at the N/C-term
    mNTmods = new ArrayList<>(); mCTmods = new ArrayList<>();

    if (sVarCharMods.containsKey('^')) mNTmods.add(sVarCharMods.get('^')); else mNTmods.add(0d);
    if (sVarCharMods.containsKey('$')) mCTmods.add(sVarCharMods.get('$')); else mCTmods.add(0d);
  }
  private void initTMT()
  {
    mAADic = Peptides.newAAsMass("tmt10");
    mAAs   = Peptides.toAAs(mAADic);

    mAACt = Maps.newHashMap();
    mAACt.put(175.1190d, 'R');
    mAACt.put(376.2757d, 'K');

    // set up an array for quick lookup
    ImmutableSetMultimap.Builder<Character, Double> builder = ImmutableSetMultimap.builder();
    builder.put('N', Peptides.N2D);
    builder.put('M', Peptides.O);
    builder.put('W', Peptides.O*2d);
    builder.put('^', Peptides.TMT10);
    builder.put('^', Peptides.TMT10_LOSS);
    builder.put('S', Peptides.TMT10);
    builder.put('H', Peptides.TMT10);
    builder.put('T', Peptides.TMT10);
    builder.put('Y', Peptides.TMT10);
    builder.put('K', Peptides.TMT10);
    builder.put('R', Peptides.TMT10);
    // HSTY, TMT over-labeling, http://pubs.acs.org/doi/full/10.1021/acs.jproteome.5b00072?src=recsys
    mVarAAMods4Exact= builder.build();

    sFacileMods=Maps.newTreeMap();
    //sVarMods.put( "*", -0.984016d); // C-terminal, any
    sFacileMods.put(162.052824d, 75f); // Hexose
    sFacileMods.put(192.063388d,  0f); // Heptose
    sFacileMods.put(79.966331d,  75f); // Phospho

    // move the static vars here
//    sVarMods= HashMultimap.create();
//    //sVarMods.put( "*", -0.984016d); // C-terminal, any
//    sVarMods.put( "*",   1.003355d); // c12->c13
//    sVarMods.put( "*",   2.006710d); // c12->2xc13
//    sVarMods.put( "N",   0.984016d); // deamidation
//    sVarMods.put( "M",  15.994915d);
//    sVarMods.put( "K",  27.994915d); // Formylation @ Nt and K
//    sVarMods.put("^C", -17.026549d); // N-terminal cyclization
//    sVarMods.put("^Q", -17.026549d); // N-terminal cyclization
//    sVarMods.put("^S", -18.01056d);  // N-terminal cyclization

    sModCharge=Maps.newTreeMap();
    sModCharge.put(Peptides.TMT10, 1f);

    sAACt =Maps.newHashMap();
    sAACt.put(175.1190d, 'R');
    sAACt.put(376.2757d, 'K');

    mTrivialMods = new double[] {Peptides.C13, Peptides.C13*2};

    // mods that ought to be placed in the 'Exact' matches because their are expected at a particular residue type
    sExactMods = MultiTreeTable.create();
    sExactMods.put("*", "M", 15.994915d);
    sExactMods.put("*", "N", 0.984016d);
    sExactMods.put("*", "S", Peptides.TMT10);
    sExactMods.put("*", "H", Peptides.TMT10);
    sExactMods.put("*", "T", Peptides.TMT10);
    sExactMods.put("*", "Y", Peptides.TMT10);
    sExactMods.put("^", "Q", -17.026549d);
    sExactMods.put("^", "C", -17.026549d);
    sExactMods.put("*", "K", Peptides.TMT10*-1d); // unmodified K by TMT6/10

    mDenovoModAAs = Peptides.modAAsMass("TMT6plex:K/^","Oxidation:M");
    mDenovoModAAs.remove("^"); mDenovoModAAs.remove("$"); mDenovoModAAs.remove("^^TMT6plex"); mDenovoModAAs.remove("K");

    // pull out the variable mods at the N/C-term
    mNTmods = new ArrayList<>(); mCTmods = new ArrayList<>();

    // always start with 'no' mod. 20190228
    mNTmods.add(0d); mCTmods.add(0d);

    if (sVarCharMods.containsKey('^')) mNTmods.add(sVarCharMods.get('^'));
    if (sVarCharMods.containsKey('$')) mCTmods.add(sVarCharMods.get('$'));
  }
  public StringBuffer printParams(StringBuffer buf)
  {
    buf = Tools.printParam(buf, "AADic", mAADic, "tmt10");
    buf = Tools.printParam(buf, "DenovoModAAs", mDenovoModAAs);
    buf = Tools.printParam(buf, "AACt", mAACt);
    buf = Tools.printParam(buf, "sAACt", sAACt);
    buf = Tools.printParam(buf, "ModCharge", sModCharge);

    buf = Tools.printParam(buf, "AAs", Strs.toString(mAAs, ","), "");
    buf = Tools.printParam(buf, "TrivialMods", Strs.toString(mTrivialMods, ",", 2), "");
    buf = Tools.printParam(buf, "NTmods", Strs.toString(mNTmods, ","), "");
    buf = Tools.printParam(buf, "CTmods", Strs.toString(mCTmods, ","), "");

    buf = Tools.printParam(buf, "VarAAMods", mVarAAMods4Exact, "");
//    buf = Tools.printParam(buf, "sVarMods", sVarMods, "");

    buf.append("ExactMods\n");
    for (String row : sExactMods.keySet())
    {
      buf.append("|"+row);
      for (String col : sExactMods.columnSet())
        buf.append("|"+(sExactMods.get(row, col)!=null?Strs.toString(sExactMods.get(row, col), ";"):" "));
      buf.append("|\n");
    }

    return buf;
  }
  @Override
  public ResidueBase clone()
  {
    ResidueBase cloned = new ResidueBase(); // do we need to clone the indices?

    if (mAADic!=null) cloned.mAADic = new HashMap<>(mAADic);

    cloned.mAAs         = Tools.cloneFloatArray(mAAs);
    cloned.mTrivialMods = Tools.cloneDoubleArray(mTrivialMods);
    cloned.mVarAAMods4Exact = ImmutableSetMultimap.copyOf(mVarAAMods4Exact);
    cloned.mModVault    = (mModVault  !=null?mModVault.clone():null);
//    cloned.sVarMods     = (sVarMods   !=null?TreeMultimap.create(sVarMods):null);
    cloned.sExactMods   = (sExactMods !=null?sExactMods.clone():null);
    cloned.sModCharge   = (sModCharge !=null?new TreeMap<>(sModCharge):null);
    cloned.sFacileMods  = (sFacileMods!=null?new TreeMap<>(sFacileMods):null);
    cloned.mNTmods      = (mNTmods    !=null?new ArrayList<>(mNTmods):null);
    cloned.mCTmods      = (mCTmods    !=null?new ArrayList<>(mCTmods):null);
    cloned.sAACt        = (sAACt      !=null? new TreeMap<>(sAACt):null);

    if (mDenovoModAAs!=null) cloned.mDenovoModAAs = new HashMap<>(mDenovoModAAs);
    if (mAACt!=null)         cloned.mAACt         = new TreeMap<>(mAACt);

    return cloned;
  }
}
