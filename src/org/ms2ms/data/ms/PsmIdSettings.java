package org.ms2ms.data.ms;

import org.expasy.mzjava.core.ms.spectrum.MsnSpectrum;
import org.ms2ms.algo.Spectra;
import org.ms2ms.graph.Property;
import org.ms2ms.utils.Settings;

public class PsmIdSettings extends Settings
{
  public static final String KEY_ACS      = "ac";
  public static final String KEY_ASSIGN   = "assign";
  public static final String KEY_BUNDLE   = "bundle";
  public static final String KEY_C13_KEEP = "c13";
  public static final String KEY_CANDIDATE= "candidate";
  public static final String KEY_CENTROID = "centroid";
  public static final String KEY_CRAP     = "crap";
  public static final String KEY_DB       = "db";
  public static final String KEY_DUPL     = "dupl";
  public static final String KEY_EVAL_N   = "bootstrap";
  public static final String KEY_EXPAND   = "exx";
  public static final String KEY_FALSE    = "neg";
  public static final String KEY_FASTA    = "fasta";
  public static final String KEY_FDR      = "fdr";
  public static final String KEY_FILTER   = "filter";
  public static final String KEY_FILTER0  = "cfg.filter";
  public static final String KEY_FINAL    = "final";
  public static final String KEY_FRAG_PPM = "ppm.f";
  public static final String KEY_KNOWN    = "known";
  public static final String KEY_LC       = "lc";
  public static final String KEY_LM       = "lm";
  public static final String KEY_MGF      = "mgf";
  public static final String KEY_MATCHES  = "matches";
  public static final String KEY_MID      = "mid";
  public static final String KEY_MOTIF    = "motif";
  public static final String KEY_MULTIZ_DB= "zs";
  public static final String KEY_MZML     = "mzml";
  public static final String KEY_NLOSS    = "nl";
  public static final String KEY_OUT      = "o";
  public static final String KEY_PREQ     = "preq";
  public static final String KEY_PROC     = "proc";
  public static final String KEY_PREC_PPM = "ppm.p";
  public static final String KEY_QVAL     = "qval";
  public static final String KEY_ROUND    = "r";
  public static final String KEY_REFS     = "ref";
  public static final String KEY_RESIDUE  = "rb";
  public static final String KEY_RUN      = "run";
  public static final String KEY_SAMPLE   = "sample";
  public static final String KEY_SAMPLES  = "samples";
  public static final String KEY_SCAN     = "scan";
  public static final String KEY_SCANINFO = "scaninfo";
  public static final String KEY_TEST     = "test";
  public static final String KEY_TOP      = "top";
  public static final String KEY_THREAD   = "thread";
  public static final String KEY_TRACK    = "track";
  public static final String KEY_TRUE     = "pos";
  public static final String KEY_WORKING  = "w";
  public static final String KEY_VERIFIED = "verified";
  public static final String KEY_ZZ_KEEP  = "zz";

  public static PsmIdSettings ORBI_ORBI = null;

  static
  {
    ORBI_ORBI = new PsmIdSettings(
    ).setPrecision(   new OffsetPpmTolerance(2.5d, 0d)
    ).setAccuracy(    new OffsetPpmTolerance(10d, 0d).isIncremental(true)
    ).setPrecursorTol(new OffsetPpmTolerance(15d, 0d));
  }

  public PsmIdSettings() { super(); }
  public PsmIdSettings(Property s)
  {
    setDefaults();

    setDuplicateProb(s.getDouble(KEY_DUPL));

    setMinTrackSize(  s.getInt(KEY_TRACK));
    setMinMotifs(     s.getInt(KEY_MOTIF));
    keepC13(          s.getBoolean(KEY_C13_KEEP));
    keepZZ(           s.getBoolean(KEY_ZZ_KEEP));

    setDefaultTracks( s.getInt(KEY_TRACK));
    setDefaultMotifs( s.getInt(KEY_MOTIF));

    useMultiZFrags(   s.getBoolean(KEY_MULTIZ_DB));
    setEvalSamples(   s.getInt(KEY_EVAL_N));
    expandExacts(     s.getBoolean(KEY_EXPAND));
    setBootstrapMid(  s.getFloat(KEY_MID));
    setFalseTags(     s.getProperty(KEY_FALSE));
    setTrueTags(      s.getProperty(KEY_TRUE));
    setRefAssignments(s.getProperty(KEY_REFS));
    setModVault(      s.getProperty(KEY_KNOWN));
    // only a subset of scans if requested
    setSelectedMSMS(  s.getProperty(KEY_SCAN));

    setResidueBase(ResidueBase.valueOf(s.getProperty(KEY_RESIDUE)));

    // set the spectrum pre-processor
    if      ( "snr".equals(s.getProperty(KEY_PROC))) setProcessor(Spectra.PrepProcessor.LOCAL_SNR);
    else if ("norm".equals(s.getProperty(KEY_PROC))) setProcessor(Spectra.PrepProcessor.LOCAL_NORM);
  }
  public PsmIdSettings setDefaults()
  {
    setPrecision(   new OffsetPpmTolerance(2.5d, 0d));
    setAccuracy(    new OffsetPpmTolerance(10d, 0d).isIncremental(true));
    setPrecursorTol(new OffsetPpmTolerance(15d, 0d));

    setMaxNullMotifs(5).setPeakTransform(0.5d);

    setMinSeriesZval(1.25);
    setMinPeptideSize(4);
    setMinMatchProb(100d);

    // values taken from Bonfire.java
    setDuplicateProb(105d);

    setMinTrackSize(  3);
    setMinMotifs(     2);
    keepC13(          true);
    keepZZ(           false);
    setDefaultTracks( 3);
    setDefaultMotifs( 2);
    useMultiZFrags(   false);
    setEvalSamples(   1000);
    expandExacts(     true);
    setBootstrapMid(  0f);
    setResidueBase(ResidueBase.TMT10);

    // set the spectrum pre-processor
    setProcessor(Spectra.PrepProcessor.LOCAL_NORM);

    return this;
  }
  public PsmIdSettings updateSettings(MsnSpectrum ms)
  {
    if (getDefaultTracks()==null) setDefaultTracks(getMinTrackSize());
    if (getDefaultMotifs()==null) setDefaultMotifs(getMinMotifs());

    double  L = (ms.getPrecursor().getMass()-2*getBase().getAAs()['^'])/getBase().getAAs()['D']; // assumming 2x TMT10 tags
    int decre = (L<8?(L<7?2:1):0);

    // shorter key need lower threshold. Assume 2+ in case of altZ
    setMinTrackSize(Math.max(1, getDefaultTracks()-decre));
    setMinMotifs(   Math.max(0, getDefaultMotifs()-decre));

    return this;
  }

  /*** Getters and Setters ***/
  //
  public Float   getBootstrapMid()   { return (Float )properties.get(KEY_MID); }
  public Integer getDefaultMotifs()  { return (Integer )properties.get("DefaultMotifs"); }
  public Integer getDefaultTracks()  { return (Integer )properties.get("DefaultTracks"); }
  public Double  getDuplicateProb()  { return (Double  )properties.get(KEY_DUPL); }
  public Integer getEvalSamples()    { return (Integer )properties.get(KEY_EVAL_N); }
  public Integer getMaxNullMotifs()  { return (Integer )properties.get("MaxNullMotifs"); }
  public Double  getMinMatchProb()   { return (Double )properties.get("MinMatchProb"); }
  public Integer getMinMotifs()      { return (Integer )properties.get(KEY_MOTIF); }
  public Integer getMinPeptideSize() { return (Integer )properties.get("MinPeptideSize"); }
  public Double  getMinSeriesZval()  { return (Double  )properties.get("MinSeriesZval"); }
  public Integer getMinTrackSize()   { return (Integer )properties.get(KEY_TRACK); }
  public Integer getNumCandidates()  { return (Integer )properties.get(KEY_CANDIDATE); }
  public String  getOutfile()        { return (String )properties.get(KEY_OUT); }
  public Double  getPeakTransform()  { return (Double  )properties.get("PeakTransform"); }
  public Double  getTrypicity()      { return (Double  )properties.get("getTrypicity"); }
  public String  getTrueTags()       { return (String )properties.get(KEY_TRUE); }

  public ResidueBase           getBase()         { return (ResidueBase )properties.get(KEY_RESIDUE); }
  public OffsetPpmTolerance    getAccuracy()     { return (OffsetPpmTolerance )properties.get(KEY_FRAG_PPM); }
  public OffsetPpmTolerance    getPrecursorTol() { return (OffsetPpmTolerance )properties.get(KEY_PREC_PPM); }
  public OffsetPpmTolerance    getPrecision()    { return (OffsetPpmTolerance )properties.get("Precision"); }
  public Spectra.PrepProcessor getProcess()      { return (Spectra.PrepProcessor )properties.get(KEY_PROC); }

  public boolean expandExacts()   { return (Boolean )properties.get(KEY_EXPAND); }
  public boolean keepC13()        { return (Boolean )properties.get(KEY_C13_KEEP); }
  public boolean keepZZ()         { return (Boolean )properties.get(KEY_ZZ_KEEP); }
  public boolean useMultiZFrags() { return (Boolean )properties.get("useMultiZFrags"); }

  public PsmIdSettings expandExacts(     Boolean s) { properties.put(KEY_EXPAND, s); return this; }
  public PsmIdSettings keepC13(          Boolean s) { properties.put(KEY_C13_KEEP, s); return this; }
  public PsmIdSettings keepZZ(           Boolean s) { properties.put(KEY_ZZ_KEEP,  s); return this; }
  public PsmIdSettings setBootstrapMid(    Float s) { properties.put(KEY_MID, s); return this; }
  public PsmIdSettings setDecoyMultiple( Double  s) { properties.put("DecoyMultiple", s); return this; }
  public PsmIdSettings setDefaultMotifs( Integer s) { properties.put("DefaultMotifs", s); return this; }
  public PsmIdSettings setDefaultTracks( Integer s) { properties.put("DefaultTracks", s); return this; }
  public PsmIdSettings  setDuplicateProb( Double s) { properties.put(KEY_DUPL, s); return this; }
  public PsmIdSettings setEvalSamples(   Integer s) { properties.put(KEY_EVAL_N, s); return this; }
  public PsmIdSettings setFalseTags(      String s) { properties.put(KEY_FALSE, s); return this; }
  public PsmIdSettings setMaxNullMotifs( Integer s) { properties.put("MaxNullMotifs", s); return this; }
  public PsmIdSettings setMinMatchProb(   Double s) { properties.put("MinMatchProb", s); return this; }
  public PsmIdSettings setMinMotifs(     Integer s) { properties.put(KEY_MOTIF, s); return this; }
  public PsmIdSettings setMinPeptideSize(Integer s) { properties.put("MinPeptideSize", s); return this; }
  public PsmIdSettings setMinSeriesZval( Double  s) { properties.put("MinSeriesZval", s); return this; }
  public PsmIdSettings setMinTrackSize(  Integer s) { properties.put(KEY_TRACK, s); return this; }
  public PsmIdSettings setModVault(       String s) { properties.put(KEY_KNOWN, s); return this; }
  public PsmIdSettings setOutfile(        String s) { properties.put(KEY_OUT, s);         return this; }
  public PsmIdSettings setPeakTransform(  Double s) { properties.put("PeakTransform", s); return this; }
  public PsmIdSettings setRefAssignments( String s) { properties.put(KEY_REFS, s); return this; }
  public PsmIdSettings setSelectedMSMS(   String s) { properties.put(KEY_SCAN, s); return this; }
  public PsmIdSettings setTrueTags(       String s) { properties.put(KEY_TRUE, s); return this; }
  public PsmIdSettings useMultiZFrags(  Boolean  s) { properties.put("useMultiZFrags", s); return this; }

  public PsmIdSettings setAccuracy(    OffsetPpmTolerance s) { properties.put(KEY_FRAG_PPM, s); return this; }
  public PsmIdSettings setPrecision(   OffsetPpmTolerance s) { properties.put("Precision", s); return this; }
  public PsmIdSettings setPrecursorTol(OffsetPpmTolerance s) { properties.put(KEY_PREC_PPM, s); return this; }
  // set the spectrum pre-processor
  public PsmIdSettings setProcessor(Spectra.PrepProcessor s) { properties.put(KEY_PROC, s); return this; }
  public PsmIdSettings setResidueBase(        ResidueBase s) { properties.put(KEY_RESIDUE, s); return this; }
}
