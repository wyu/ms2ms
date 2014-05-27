package org.ms2ms.test;

import org.apache.hadoop.hbase.client.HConnection;
import org.apache.hadoop.hbase.client.HConnectionManager;
import org.apache.hadoop.hbase.client.HTableInterface;
import org.apache.hadoop.hbase.util.Bytes;
import org.expasy.mzjava.core.ms.AbsoluteTolerance;
import org.expasy.mzjava.core.ms.PpmTolerance;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.proteomics.io.ms.spectrum.MsLibReader;
import org.expasy.mzjava.proteomics.io.ms.spectrum.SptxtReader;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodModificationResolver;
import org.expasy.mzjava.proteomics.ms.spectrum.LibrarySpectrum;
import org.expasy.mzjava.proteomics.ms.spectrum.PepLibPeakAnnotation;
import org.expasy.mzjava.utils.URIBuilder;
import org.junit.Test;
import org.ms2ms.alg.Peaks;
import org.ms2ms.mimsl.MIMSL;
import org.ms2ms.mimsl.MimslSettings;
import org.ms2ms.mzjava.AnnotatedPeak;
import org.ms2ms.mzjava.AnnotatedSpectrum;
import org.ms2ms.nosql.HBase;
import org.ms2ms.nosql.HBasePeakList;
import org.ms2ms.nosql.HBaseProteomics;
import org.ms2ms.nosql.HBaseSpLib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.List;

/** Test the application logic
 *
 * Created by wyu on 4/22/14.
 */
public class MrmInspiredML extends TestAbstract
{
  @Test
  public void testRun() throws IOException
  {
  }

  /** suite of tests to back up the manuscript */
  @Test
  public void queryTest() throws IOException
  {
    // 500.730, +2(6): 318.20,520.19,568.30,683.25,782.32,869.35,
//    PeakList<PepLibPeakAnnotation> ions = Peaks.newPeakList(500.73d, 2, "318.2", "568.3"); //782.32,869.35
    // nist_nci_stdmix_consensus_final_true_lib.msp
    // C.AAADPHECYAK.V/2, m/z588.2642; 375.67/2, 424.19/2, 481.71/2, 517.22/2, 552.74/2, 613.27/1, 750.32/1, 795.31/1, 847.38/1, 958.37/1,962.40/1, 1029.41/1, 1033.44/1
    PeakList<PepLibPeakAnnotation> ions = Peaks.newPeakList(588.3d, 2, "517.2", "613.3");
    List<AnnotatedSpectrum>  candidates = MIMSL.run(ions, HBasePeakList.SPEC_TRAP_CID, new PpmTolerance(500d), new AbsoluteTolerance(0.5));

    System.out.println(MIMSL.printCandidates(null, candidates));
  }
  // test the recovery of known spectra at various parameter
  @Test
  public void sampleRecovery() throws IOException
  {
    HBaseProteomics.sampleRecovery("/media/data/splib/2013/HumanPlasma_2012-08_all.sptxt", 1, 1000, MimslSettings.ORBI_HL_CID);
  }


/*  String libs_dir = "/bioinfo/apps/linux/ms/twin/libs/nist/";
  Multimap<Double, PTM> ptms;
  static MsMsDictionary dictionary = null;

  static final String hsa = "hsa_consensus_final_true_lib", human = "human_consensus_final_true_lib",
    stdmix = "nist_nci_stdmix_consensus_final_true_lib",  hsa_selected = "hsa_selected_final_true_lib",
    crp = "human_crp_consensus_final_true_lib",        crp_selected = "human_crp_selected_final_true_lib",
    b2mg = "human_b2mg_consensus_final_true_lib",      sigma = "sigmaups1_consensus_final_true_lib";
  // combining the assignments from multiple libraries and indice it
  public void vtestMimslLibPreparation() throws Exception
  {
    SpectralLib_Util.prepareMimslLib("/data/library/human_cid_it_20130625", "/data/sequences/cp9606.fasta", 8d, 50d, 7,
      new BinPepLibFile("/data/library/human_plus.bin"),
      new BinPepLibFile("/data/library/human_20130623_it_cid.bin"));
  }
  public void vtestMimslIndivLibPreparation() throws Exception
  {
    SpectralLib_Util.prepareMimslLib("/data/library/human_20130625_test.bin", "/data/sequences/cp9606.fasta", 8d, 50d, 7,
      new BinPepLibFile("/data/library/human_test_it_cid.bin"));
  }
  public void vtestQueryBySignature() throws Exception
  {
    // YSVFFQSLAVVEQEMK, a high quality peptide ID with high background, in the library
    //List<MsIon> signatures = MsIon_Util.setCharge(newSignatureFromMsMsSpectrum(conn, 68568244L), 0);
    //List<MsIon> signatures = MsIon_Util.setCharge(newSignatureFromMsMsSpectrum(conn, 68539694L), 0);
    //List<MsIon> signatures = MsIon_Util.setCharge(newSignatureFromMsMsSpectrum(conn, 68563004L), 0);
    List<MsIon> signatures = MsIon_Util.setCharge(newSignatureFromMsMsSpectrum(conn, 68554927L), 0);

    System.out.println("The Signature ions generated from the PSM:");
    for (MsIon ion : signatures)
      System.out.println(ion.toString());

    Mimp mimp = new Mimp(MsMsDictionary.LIB_HUMAN);
    mimp.setSettings(MsSettings.ORBITRAP);

    mimp = randomQuery(mimp, signatures, 6, null, new Random(System.nanoTime()));
    mimp.announce("", 25, 10d);

    //mimp = randomQuery(mimp, signatures, 4, new MsMsPeak(610.362, 0), new Random(System.nanoTime()));
    //mimp.announce("", 25, 10d);
  }
  public void vtestModQueryBySignature() throws Exception
  {
    // IEDVGsDEEDDSGK, a high quality peptide ID with high background, in the library
    List<MsIon> signatures = MsIon_Util.setCharge(newSignatureFromMsMsSpectrum(conn, 76547545L), 0);

    for (MsIon ion : signatures)
      System.out.println(ion.toString());

    Mimp mimp = new Mimp(MsMsDictionary.LIB_HUMAN);
    mimp.setSettings(MsSettings.ORBITRAP);
    //mimp.toModsArtifact(true);

    mimp.toModsTranslational(true);

    mimp = randomQuery(mimp, signatures, 6, new MsMsPeak(787.8, 0d, 2), new Random(System.nanoTime()));
    mimp.announce("", 25, 10d);
  }
  // starting a list of peptide ids that are known to exist in the library, test for the recovery
  public void vtestBayesianProphetWithPeptideIDs() throws Exception
  {
    Collection<Long> msids = Spectre_MsMsSpectrum.getRowIdsAnnotatedForAnalysisRun(conn, 59885L, 2, 5d, new Range<Double>(0d, 670d));
    Multimap<String, MsMsAssignment> seq_assignment =
      Qualitative_Util.toMsMsAssignments(Spectre_MsMsSpectrum.getRowsAnnotated(conn, msids, MsMsSpectrumRetrievalFlags.ASSIGNMENT));
    // mark the known assignments as "ref"
    int known = MsMsDictionary.subset(MsMsDictionary.LIB_HUMAN, seq_assignment, true).size();

    System.out.println();
    System.out.println("Reference Assignments (in-library/totals): " + known + "/" + seq_assignment.size());
    Mimp mimp = new Mimp(MsMsDictionary.LIB_HUMAN);
    mimp.setSettings(MsSettings.ORBITRAP);
    Mimp.sVerbose = false;

    // deal with the best case scenario first
    BayesianProphet prophet = newBayesianProphet(conn, mimp, seq_assignment.values(), 6, 8d, 50d, 7, true);
    System.out.println(Wiki_Util.newChart(prophet, null));
  }

  public void vtestRecoveryWithPeptideIDs() throws Exception
  {
    Collection<Long> msids = Spectre_MsMsSpectrum.getRowIdsAnnotatedForAnalysisRun(conn, 59885L, 2, 5d, null);
    Multimap<String, MsMsAssignment> seq_assignment =
      Qualitative_Util.toMsMsAssignments(Spectre_MsMsSpectrum.getRowsAnnotated(conn, msids, MsMsSpectrumRetrievalFlags.ASSIGNMENT)),
      seq_assignment_pos = MsMsDictionary.subset(MsMsDictionary.LIB_HUMAN, seq_assignment, true);

    System.out.println();
    System.out.println("Reference Assignments (in-library/totals): " + seq_assignment_pos.size() + "/" + seq_assignment.size());
    Mimp mimp = new Mimp(MsMsDictionary.LIB_HUMAN);
    mimp.setSettings(MsSettings.LTQ);
    Mimp.sVerbose = false;

    // deal with the best case scenario first
    MatchStats stats6true = randomQuery(conn, mimp, seq_assignment_pos.values(), 6, 4d, true, true, new Random(System.nanoTime()));
  }
  public void vtestNegativeRecoveryWithPeptideIDs() throws Exception
  {
    Collection<Long> msids = Spectre_MsMsSpectrum.getRowIdsAnnotatedForAnalysisRun(conn, 59885L, 2, 5d, null);
    Multimap<String, MsMsAssignment> seq_assignment =
      Qualitative_Util.toMsMsAssignments(Spectre_MsMsSpectrum.getRowsAnnotated(conn, msids, MsMsSpectrumRetrievalFlags.ASSIGNMENT)),
      seq_assignment_neg = MsMsDictionary.subset(MsMsDictionary.LIB_HUMAN, seq_assignment, false);

    System.out.println();
    System.out.println("Reference Assignments (not-in-library/totals): " + seq_assignment_neg.size() + "/" + seq_assignment.size());
    Mimp mimp = new Mimp(MsMsDictionary.LIB_HUMAN);
    mimp.setSettings(MsSettings.LTQ);
    Mimp.sVerbose = false;

    // deal with the best case scenario first, expect non match
    MatchStats stats6true = randomQuery(conn, mimp, seq_assignment_neg.values(), 6, 4d, true, false, new Random(System.nanoTime()));
  }
  public void vtestPostives() throws Exception
  {
    List<ProteinId> pids = ProteinId_Util.readProteinIDs("/tmp/qualifiedID.pid");
    Collection<MsMsAssignment> assignments = new ArrayList<MsMsAssignment>();


    Collection<ProteinId> selected = new ArrayList<ProteinId>();
    for (int i = 0; i < 10; i++) selected.add(pids.get(i));
    Db_Util.populateMsMs(conn, selected, MsMsSpectrumRetrievalFlags.ION_ATTR, true);

    for (ProteinId p : selected)
    {
      for (MsMsAssignment A : Toolbox.values(p.getAssignments()))
        if (A.getAssignment().getErrorPct() < 5d && A.getConcensusSize() > 2)
        {
          assignments.add(Qualitative_Util.populateIonAnnotations(A, IonActivationType.CID, 0.5f));
        }
    }
  }

  public void vtestSaveHumanBinaryLib() throws Exception
  {
    String lib = libs_dir + "humans";
    SpectralLib_Util.Src2Bin(lib + ".bin",
      new NistPepLibFile(libs_dir + human        + ".msp"),
      new NistPepLibFile(libs_dir + hsa          + ".msp"),
      new NistPepLibFile(libs_dir + hsa_selected + ".msp"),
      new NistPepLibFile(libs_dir + stdmix       + ".msp"),
      new NistPepLibFile(libs_dir + crp          + ".msp"),
      new NistPepLibFile(libs_dir + crp_selected + ".msp"),
      new NistPepLibFile(libs_dir + b2mg         + ".msp"),
      new NistPepLibFile(libs_dir + sigma        + ".msp"));

    dictionary = new MsMsDictionary().index(new BinPepLibFile(lib + ".bin"), 10d,  50d, 7);

    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(lib + ".indices", "rw");
    dictionary.write(indices);
    indices.close();

    //testQueryBy3Frag();
  }
  public void vtestSaveMiscBinaryLib() throws Exception
  {
    String lib = libs_dir + "miscs";
    SpectralLib_Util.Src2Bin(lib + ".bin",
      new NistPepLibFile(libs_dir + hsa          + ".msp"),
      new NistPepLibFile(libs_dir + hsa_selected + ".msp"),
      new NistPepLibFile(libs_dir + stdmix       + ".msp"),
      new NistPepLibFile(libs_dir + crp          + ".msp"),
      new NistPepLibFile(libs_dir + crp_selected + ".msp"),
      new NistPepLibFile(libs_dir + b2mg         + ".msp"),
      new NistPepLibFile(libs_dir + sigma        + ".msp"));

    dictionary = new MsMsDictionary().index(new BinPepLibFile(lib + ".bin"), 10d, 50d, 7);
    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(lib + ".indices", "rw");

    dictionary.write(indices);
    indices.close();

    //testQueryBy3FragMIsc();
  }
  public void vtestReadCombinedBinaryLib() throws Exception
  {
    String lib = libs_dir + "combined";

    dictionary = new MsMsDictionary().index(new BinPepLibFile(lib + ".bin"), 10d, 50d, 7);
    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(lib + ".indices", "rw");

    dictionary.write(indices);
    indices.close();

    //testPubic();
  }
  public void vtestSaveISBLib() throws Exception
  {
    String lib = libs_dir + "combined",
      root_misc = "/bioinfo/apps/linux/ms/twin/libs/misc/",
      root_isb = root_misc + "ISB_Hs_plasma_consensus_20070706_PUBLIC/";
    SpectralLib_Util.Src2Bin(lib + ".bin",
      new NistPepLibFile(root_isb + "ISB_Hs_plasma_consensus_20070706_PUBLIC_Q1.sptxt"),
      new NistPepLibFile(root_isb + "ISB_Hs_plasma_consensus_20070706_PUBLIC_Q2.sptxt"),
      new NistPepLibFile(root_isb + "ISB_Hs_plasma_consensus_20070706_PUBLIC_Q0.sptxt"));

    dictionary = new MsMsDictionary().index(new BinPepLibFile(lib + ".bin"), 10d, 50d, 7);
    BufferedRandomAccessFile indices = new BufferedRandomAccessFile(lib + ".indices", "rw");

    dictionary.write(indices);
    indices.close();
  }
  public void vtestPrecursorEnumeration() throws Exception
  {
    MultiTreeMap<Float, PTM> unimods = Sequence_Util.newUnimodPTMs();
    Collection<Range<Float>> slices = Mimp.enumeratePrecursor(MsSettings.QTOF, unimods.keySet(), MsIon_Util.toIons(466.9d).values(), 3);
    for (Range<Float> slice : slices)
      System.out.println(slice.toString());
  }
  public void vtestHuman() throws Exception
  {
    //doQuery(466.9d, MsMsDictionary.LIB_HUMAN, 629.1d);
    doQuery(466.9d, MsMsDictionary.LIB_HUMAN, 629.1d, 549.0d, 664.6d);
    //doQuery(null, MsMsDictionary.LIB_HUMAN, 629.1d, 549.0d, 664.6d);
    //doQuery(null, MsMsDictionary.LIB_HUMAN, 629.1d, 549.0d, 664.6d, 492.5d);
  }
  public void vtestHumanMod() throws Exception
  {
    //doMod(679.819, 2, MsMsDictionary.LIB_MISC, 1041.3, 1188.4, 811.24);
    doMod(679.36, 3, MsMsDictionary.LIB_HUMAN, 795.8, 944.1, 982.4, 1054.5);
  }
  public void vtestQueryBy3FragMIsc() throws Exception
  {
    doQuery(466.9d, MsMsDictionary.LIB_MISC, 629.1d, 549.0d, 664.6d);
  }
  private void doQuery(Double precursor, String lib, Double... frags) throws Exception
  {
    Mimp mimp = new Mimp(lib);
    mimp.setSettings(MsSettings.LTQ);
    if (precursor != null) mimp.setPrecursors(precursor);
    mimp.setFrags(frags);
    mimp.run();

    mimp.announce("The final candidates", 1000, 10d);
  }
  private void doQuery(MsMsAssignment assignment, Mimp mimp, int picks)  throws Exception
  {
    List<MsIon> signatures = MimpCore.extractSignatureFragments(assignment, 10d, 50d, 7);

    mimp.setPrecursors(assignment.getMsMs().getPrecursorMz());

    // scramble the order of the fragments
    Collections.sort(signatures, MsPoint.Y_DESCEND);

    mimp.setFrags(signatures.subList(0, picks));

    mimp.run();

    mimp.announce("The final candidates", 1000, 10d);
  }
  private void doMod(Double precursormz, int precursorz, String libname, Double... frags) throws Exception
  {
    Mimp mimp = new Mimp(libname);
    mimp.setSettings(MsSettings.ORBITRAP);
    mimp.toModsArtifact(true);

    mimp.setPrecursor(precursormz, precursorz);
    mimp.setFrags(frags);
    mimp.run();

    mimp.announce("The final candidates", 1000, 10d);
  }
  private static boolean inLibrary(String libname, Spectre_MsMsSpectrum ms) throws Exception
  {
    Multimap<String, MsMsAssignment> seq_assignment = HashMultimap.create(), seq_assignment_pos = null;

    seq_assignment.put(ms.getPeptideAssignment().getSequence(), ms.toMsMsAssignment());
    seq_assignment_pos = MsMsDictionary.subset(libname, seq_assignment, true);

    return Toolbox.isSet(seq_assignment_pos);
  }
  private static List<MsIon> newSignatureFromMsMsSpectrum(Connection conn, Long msid) throws Exception
  {
    Spectre_MsMsSpectrum spec = Spectre_MsMsSpectrum.getRow(conn, msid, MsMsSpectrumRetrievalFlags.ASSIGNMENT_ION);
    spec.getPeptideAssignment().matchToSpectrum(conn, spec.getRawSpectrum(), 0.5f, true);

    return MimpCore.extractSignatureFragments(spec.getRawSpectrum(), 4d, 50d, 7);
  }
  private static Mimp randomQuery(Mimp mimp, List<MsIon> signatures, int sample, MsIon precursor, Random rand) throws Exception
  {
    if (precursor != null) mimp.setPrecursor(precursor.getMz(), precursor.getCharge());

    //if (rand == null) rand = new Random(System.nanoTime());
    // intensity_bias shuld already be set
    mimp.setFrags(rand != null ? Toolbox.sample(signatures, sample, rand) :
      Toolbox.sample(signatures, sample, MsIon.INTENSITY_DESCEND));
    mimp.run();

    //mimp.announce("The final candidates", 1000);
    return mimp;
  }
  private static MatchStats randomQuery(Connection conn, Mimp mimp,
    Collection<MsMsAssignment> peptides,
    int sample,
    double min_snr,
    boolean use_precursor,
    boolean match,
    Random rand) throws Exception
  {
    MatchStats stats = new MatchStats();
    //Random      rand = new Random(System.nanoTime());
    Histogram deltas = new Histogram();

    int order = 0;
    for (MsMsAssignment peptide : peptides)
    {
      //System.out.println(">>" + peptide.getAssignment().toString());
      if (!Toolbox.isSet(peptide.getMsMs().getRawData()))
      {
        Spectre_MsMsSpectrum row = Spectre_MsMsSpectrum.getRow(conn, peptide.getMsMs().getPrimaryKey(), MsMsSpectrumRetrievalFlags.ASSIGNMENT_ION);
        row.getPeptideAssignment().matchToSpectrum(conn, row.getRawSpectrum(), 0.5f, true);
        peptide.getMsMs().setData(row.getRawSpectrum().getData());
        //peptide.setMsMs(row.toMsMsAssignment().getMsMs());
      }
      List<MsIon> signatures = MimpCore.extractSignatureFragments(peptide, min_snr, 50d, 7);
      mimp = randomQuery(mimp, MsIon_Util.setCharge(signatures, 0), sample,
        use_precursor && peptide.getMsMs() != null ? peptide.getMsMs().getPrecursor() : null, rand);

      Map<PeptideHit, MsMsAssignment> candidates = mimp.fdr(10, min_snr);
      //boolean matched = mimp.isMatch(candidates, 10d, match ? peptide.getAssignment() : null, 1, Math.log(1));
      boolean matched = mimp.isMatch(candidates, 10d, 1, Math.log(1));
      if (matched) stats.A_cnt++; else stats.B_cnt++;

      if (match != matched)
      {
        //System.out.println(Mimp.sMessage.toString());
        System.out.println();
        System.out.println(Wiki_Util.newMsMsWiki(peptide.getMsMs().getName() + ": " + peptide.getConcensusSize() + " votes @ " + peptide.getAssignment().getErrorPct() + "%", null, peptide.getMsMs().getPrimaryKey()));
        for (MsIon ion : signatures)
        {
          boolean used = mimp.getFragmentMzs().values().contains(ion);
          System.out.print((used ? "*" : "") + Toolbox.d2s(ion.getMz(), 3) + "/" + Toolbox.d2s(ion.getIntensity(), 0) + (used ? "*" : "") + "; ");
        }
        System.out.println();
        mimp.announce(peptide.getAssignment().toString(), candidates);
        System.out.println("hit/miss: " + stats.A_cnt + "/" + stats.B_cnt);
      }
      if (Toolbox.isSet(candidates))
      {
        PeptideHit top = Toolbox.front(candidates.keySet());
        if (top != null && peptide != null)
          if (Toolbox.equals(top.getSequence(), peptide.getAssignment().getSequence()) == match) deltas.add(top.getScoreDelta());
      }
      Mimp.sMessage = new StringBuffer();

      if (++order % 50   == 0) System.out.print(".");
    }
    deltas.generate(50);
    System.out.println(Wiki_Util.newChart(deltas));
    System.out.println("hit/miss: " + stats.A_cnt + "/" + stats.B_cnt);

    return stats;
  }
  private static BayesianProphet newBayesianProphet(Connection conn, Mimp mimp,
    Collection<MsMsAssignment> peptides,
    int sample, double min_snr, double half_width, int tops, boolean use_precursor) throws Exception
  {
    MatchStats stat = new MatchStats();
    Histogram score_fp = new Histogram("Scores of the false positives");
    BayesianProphet prophet = new BayesianProphet("Score Differentials", new Range<Double>(0d, 25d), 50),
      scores = new BayesianProphet("Scores", new Range<Double>(0d, 50d), 50);
    //Random      rand = new Random(System.nanoTime());
    Random      rand = null;

    System.out.println("New Bayesian Prophet>>" + (use_precursor ? " with precursor" : "") + ", window/top/snr=" + half_width + "/" + tops + "/" + min_snr);
    System.out.println("with sample=" + sample + (rand != null ? " randomly selected fragments." : "most intensed fragments by SNR."));

    int order = 0, found = 0, unfound = 0;
    for (MsMsAssignment peptide : peptides)
    {
      //System.out.println(">>" + peptide.getAssignment().toString());
      Spectre_MsMsSpectrum row = null;
      if (!Toolbox.isSet(peptide.getMsMs().getRawData()))
      {
        row = Spectre_MsMsSpectrum.getRow(conn, peptide.getMsMs().getPrimaryKey(), MsMsSpectrumRetrievalFlags.ASSIGNMENT_ION);
        row.getPeptideAssignment().matchToSpectrum(conn, row.getRawSpectrum(), 0.5f, true);
        peptide.getMsMs().setData(row.getRawSpectrum().getData());
      }
      List<MsIon> signatures = MimpCore.extractSignatureFragments(peptide, min_snr, half_width, tops);
      mimp = randomQuery(mimp, MsIon_Util.setCharge(signatures, 0), sample,
        use_precursor && peptide.getMsMs() != null ? peptide.getMsMs().getPrecursor() : null, rand);

      Map<PeptideHit, MsMsAssignment> candidates = mimp.fdr(10, min_snr);

      // check by the score, rank and delta score
      boolean       matched = mimp.isMatch(candidates, 10d, 1, 3d);
      PeptideHit top_ranked = Toolbox.isSet(candidates) ? Toolbox.front(candidates.keySet()) : null;

      if (candidates != null && candidates.size() == 1)
        if (matched) stat.A_cnt++; else stat.B_cnt++;

      // record the counts
      if (peptide.is(Candidate.eVerdict.ref))
        if (matched) stat.ref_assigned_cnt++; else
        {
          stat.ref_missed_cnt++;
          if (top_ranked != null && top_ranked.getScoreDelta() != null && top_ranked.getScoreDelta() < 1 &&
            top_ranked.getScore() != null && top_ranked.getScore() > 15)
          {
            System.out.println("\nh4. False Negative with delta < 1, score > 15: ");
            message(mimp, signatures, peptide, candidates);
          }
        }
      else
      {
        if (matched)
        {
          stat.decoy_assigned_cnt++;
          if (top_ranked != null && top_ranked.getScoreDelta() != null && top_ranked.getScoreDelta() > 10)
          {
            System.out.println("\nh4. False Positive with delta > 10: ");
            message(mimp, signatures, peptide, candidates);
          }
        }
        else stat.decoy_missed_cnt++;
      }

      if (top_ranked != null)
      {
        if (peptide.is(Candidate.eVerdict.ref))
        {
          if (top_ranked.getScoreDelta() != null) prophet.addPositive(top_ranked.getScoreDelta());
          scores.addPositive(top_ranked.getScore());
        }
        else
        {
          if (top_ranked.getScoreDelta() != null)
          {
            prophet.addNegative(top_ranked.getScoreDelta());
          }
          if (matched) score_fp.add(top_ranked.getScore());
          scores.addNegative(top_ranked.getScore());
        }
      }
      Mimp.sMessage = new StringBuffer();

      Toolbox.dispose(row);
      Toolbox.dispose(candidates);
      Toolbox.dispose(signatures);
      //Toolbox.dispose(peptide.getMsMs().getRawData());

      if (++order % 100   == 0) System.out.print(".");
    }
    System.out.println("found/unfound/totals: " + found + "/" + unfound + "/" + order);
    System.out.println("score delta, positives/negatives: " + prophet.getPositives().getTotals() + "/" + prophet.getNegatives().getTotals());
    System.out.println("      score, positives/negatives: " +  scores.getPositives().getTotals() + "/" +  scores.getNegatives().getTotals());

    System.out.println("singleton match/not: " + stat.A_cnt + "/" + stat.B_cnt);
    System.out.println("  ref assign/miss: " + stat.ref_assigned_cnt + "/" + stat.ref_missed_cnt);
    System.out.println("decoy assign/miss: " + stat.decoy_assigned_cnt + "/" + stat.decoy_missed_cnt);
    System.out.println("false positive/false negative: " +
      Toolbox.d2s((double )stat.decoy_assigned_cnt / (stat.decoy_assigned_cnt + stat.decoy_missed_cnt), 2) + "/" +
      Toolbox.d2s((double )stat.ref_missed_cnt     / (stat.ref_assigned_cnt   + stat.ref_missed_cnt), 2));

    System.out.println(Wiki_Util.newChart(score_fp.generate(50)));
    System.out.println(Wiki_Util.newChart(scores, null));

    return prophet;
  }
  private static void message(Mimp mimp, List<MsIon> signatures, MsMsAssignment peptide, Map<PeptideHit, MsMsAssignment> candidates)
  {
    //System.out.println(Mimp.sMessage.toString());
    System.out.println("\n" + signatures.size());
    System.out.println(Wiki_Util.newMsMsWiki(peptide.getMsMs().getName() + ": " + peptide.getConcensusSize() + " votes @ " + peptide.getAssignment().getErrorPct() + "%", null, peptide.getMsMs().getPrimaryKey()));
    for (MsIon ion : signatures)
    {
      boolean used = mimp.getFragmentMzs().values().contains(ion);
      System.out.print((used ? "*" : "") + Toolbox.d2s(ion.getMz(), 3) + "/" + Toolbox.d2s(ion.getIntensity(), 0) + (used ? "*" : "") + "; ");
    }
    System.out.println();
    mimp.announce("+" + peptide.getAssignment().getCharge() + ", " + peptide.getAssignment().toString(), candidates);
  }
  public static void newRepoFromProteinMatrix(Connection conn, String matrix, String libname, int row_limit) throws Exception
  {
    String binlib = MsMsDictionary.libsdir + "/" + libname;
    //SpectralLib_Util.saveProteinIdsAsBinaryLib(conn, binlib, 0.15f, ProteinId_Util.newProteinMatrix_View("1335991687500.bin", "wyu"));

    ProteinIDCaller         caller = new ProteinIDCaller();
    ProteinIdReportMgr  reportMgr  = new ProteinIdReportMgr(Spectre.getUserSpectreDir("wyu"));
    ProteinIdReportInfo reportInfo = reportMgr.getReportInfo(matrix);
    BufferedRandomAccessFile  file = new BufferedRandomAccessFile(reportInfo.getReportFile(), "r");

    Toolbox.isnull(file);
    Map<String, Long> peptide_pointer = caller.readAssignmentIndex(file, true);

    SpectralLib_Util.saveProteinIdsAsBinaryLib(conn, caller, file, peptide_pointer, row_limit,
      SpectralLib_Util.newBinPepLibFilesByInstIonActivation(binlib));

    file.close();
  }
  public static void doProteinIds(String... binfiles) throws Exception
  {
    List<ProteinId> pids = Qualitative_Util.newProteinIDs(
      Sequence_Util.newProteins("/bioinfo/scratch/wyu/staging/sequences/cp_reviewed_9606.fasta"));

    for (String filename : binfiles)
    {
      BinPepLibFile bin = new BinPepLibFile(filename);
      pids = ProteinId_Util.doProteinIds(pids, bin);
    }
  }*/
}
