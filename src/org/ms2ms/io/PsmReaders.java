package org.ms2ms.io;

import com.google.common.base.Optional;
import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import info.monitorenter.cpdetector.io.FileFilterExtensions;
import org.expasy.mzjava.core.mol.NumericMass;
import org.expasy.mzjava.core.ms.spectrum.TimeUnit;
import org.expasy.mzjava.proteomics.io.ms.ident.MzIdentMlReader;
import org.expasy.mzjava.proteomics.io.ms.ident.PSMReaderCallback;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.alg.Peptides;
import org.ms2ms.mzjava.NumModResolver;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/**
 * Created by yuw on 9/14/2015.
 */
public class PsmReaders
{
  private final Map<SpectrumIdentifier, List<PeptideMatch>> searchResultMap = new HashMap<>();

  PSMReaderCallback insertIdResultCB = new PSMReaderCallback() {
    @Override
    public void resultRead(SpectrumIdentifier identifier, PeptideMatch searchResult) {

      List<PeptideMatch> results;
      if (searchResultMap.containsKey(identifier)) {
        results = searchResultMap.get(identifier);
      } else {
        results = new ArrayList<>();
        searchResultMap.put(identifier, results);
      }

      results.add(searchResult);
    }
  };

  public Map<SpectrumIdentifier, List<PeptideMatch>> getResultMap() { return searchResultMap; }

  public void parseMSGFplusMzID(String filename)
  {
    MzIdentMlReader reader = new MzIdentMlReader(new NumModResolver());

    FileInputStream file = null;
    try
    {
      try
      {
        file = new FileInputStream(new File(filename));
        reader.parse(new FileInputStream(new File(filename)), insertIdResultCB);
      }
      finally
      {
        if (file!=null) file.close();
      }
    }
    catch (FileNotFoundException e1)
    {
      throw new RuntimeException("File not found: " + filename);
    }
    catch (IOException e2)
    {
      throw new RuntimeException("File I/O exception!");
    }
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMascotPD(String filename) throws IOException
  {
    TabFile file = new TabFile(filename, TabFile.comma);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = file.get("Spectrum File"), scan = file.get("First Scan");
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run+"#"+scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("m/z [Da]"));
        id.setPrecursorIntensity(file.getDouble("Intensity"));
        id.setSpectrumFile(run);
        id.addRetentionTime(file.getDouble("RT [min]"), TimeUnit.MINUTE);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = new PeptideMatch(file.get("Sequence").toUpperCase());
      String[] mods = Strs.split(file.get("Modifications"), ';');
      if (Tools.isSet(mods))
        for (String mod : mods)
        {
          String tag = mod.substring(0, mod.indexOf('(')).trim(), m = mod.substring(mod.indexOf('(')+1, mod.indexOf(')'));
          if      (Strs.equals(tag, "N-Term")) match.addModificationMatch(ModAttachment.N_TERM, new Modification(m, new NumericMass(0d)));
          else if (Strs.equals(tag, "C-Term")) match.addModificationMatch(ModAttachment.C_TERM, new Modification(m, new NumericMass(0d)));
          else
          {
            match.addModificationMatch(new Integer(tag.substring(1, tag.length()))-1, new Modification(m, new NumericMass(0d)));
          }
        }
      match.addProteinMatch(new PeptideProteinMatch(file.get("Protein Group Accessions"),
          Optional.of(file.get("Protein Descriptions")), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("?M [ppm]")*id.getPrecursorMz().get()*id.getAssumedCharge().get()*1E-6);
      match.setNumMatchedIons(new Integer(matched[0]));
      match.setTotalNumIons(new Integer(matched[1]));
      match.setNumMissedCleavages(file.getInt("# Missed Cleavages"));
      match.setRank(file.getInt("Rank"));
      match.addScore("IonScore", file.getDouble("IonScore"));
      match.addScore("Exp Value", file.getDouble("Exp Value"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readSequestHT(String filename) throws IOException
  {
    TabFile file = new TabFile(filename, TabFile.tabb);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = file.get("Spectrum File"), scan = file.get("First Scan");
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("m/z [Da]"));
//        id.setPrecursorIntensity(file.getDouble("Intensity"));
        id.setSpectrumFile(run);
        id.addRetentionTime(file.getDouble("RT [min]"), TimeUnit.MINUTE);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = new PeptideMatch(file.get("Annotated Sequence").toUpperCase());
      String[] mods = Strs.split(file.get("Modifications"), ';');
      if (Tools.isSet(mods))
        for (String mod : mods)
        {
          String tag = mod.substring(0, mod.indexOf('(')).trim(), m = mod.substring(mod.indexOf('(')+1, mod.indexOf(')'));
          if      (Strs.equals(tag, "N-Term")) match.addModificationMatch(ModAttachment.N_TERM, new Modification(m, new NumericMass(0d)));
          else if (Strs.equals(tag, "C-Term")) match.addModificationMatch(ModAttachment.C_TERM, new Modification(m, new NumericMass(0d)));
          else
          {
            match.addModificationMatch(new Integer(tag.substring(1, tag.length()))-1, new Modification(m, new NumericMass(0d)));
          }
        }
      match.addProteinMatch(new PeptideProteinMatch(file.get("Master Protein Accessions"),
          Optional.of("-"), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("DeltaM [ppm]")*id.getPrecursorMz().get()*id.getAssumedCharge().get()*1E-6);
//      match.setNumMatchedIons(new Integer(matched[0]));
//      match.setTotalNumIons(new Integer(matched[1]));
      match.setNumMissedCleavages(file.getInt("# Missed Cleavages"));
      match.setRank(file.getInt("Rank"));
      match.addScore("XCorr", file.getDouble("XCorr"));
      match.addScore("Percolator q-Value", file.getDouble("Percolator q-Value"));
      match.addScore("Percolator PEP", file.getDouble("Percolator PEP"));
      if (file.getDouble("DeltaScore")!=null) match.addScore("DeltaScore", file.getDouble("DeltaScore"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }

    return id_match;
  }
  //allPeptide//Raw file	Type	Charge	m/z	Mass	Uncalibrated m/z	Resolution	Number of data points	Number of scans	Number of isotopic peaks	PIF	Mass fractional part	Mass deficit	Mass precision [ppm]	Max intensity m/z 0	Retention time	Retention length	Retention length (FWHM)	Min scan number	Max scan number	Identified	MS/MS IDs	Sequence	Length	Modifications	Modified sequence	Proteins	Score	Intensity	Intensities	MS/MS Count	MSMS Scan Numbers	MSMS Isotope Indices
  //evidence//Sequence	Length	Modifications	Modified sequence	Oxidation (M) Probabilities	Oxidation (M) Score Diffs	Acetyl (Protein N-term)	Oxidation (M)	Missed cleavages	Proteins	Leading proteins	Leading razor protein	Type	Raw file	Fraction	Experiment	MS/MS m/z	Charge	m/z	Mass	Resolution	Uncalibrated - Calibrated m/z [ppm]	Uncalibrated - Calibrated m/z [Da]	Mass Error [ppm]	Mass Error [Da]	Uncalibrated Mass Error [ppm]	Uncalibrated Mass Error [Da]	Max intensity m/z 0	Retention time	Retention length	Calibrated retention time	Calibrated retention time start	Calibrated retention time finish	Retention time calibration	Match time difference	Match m/z difference	Match q-value	Match score	Number of data points	Number of scans	Number of isotopic peaks	PIF	Fraction of total spectrum	Base peak fraction	PEP	MS/MS Count	MS/MS Scan Number	Score	Delta score	Combinatorics	Intensity	Reporter intensity 0	Reporter intensity 1	Reporter intensity 2	Reporter intensity 3	Reporter intensity 4	Reporter intensity 5	Reporter intensity 6	Reporter intensity 7	Reporter intensity 8	Reporter intensity 9	Reporter intensity not corrected 0	Reporter intensity not corrected 1	Reporter intensity not corrected 2	Reporter intensity not corrected 3	Reporter intensity not corrected 4	Reporter intensity not corrected 5	Reporter intensity not corrected 6	Reporter intensity not corrected 7	Reporter intensity not corrected 8	Reporter intensity not corrected 9	Reporter intensity count 0	Reporter intensity count 1	Reporter intensity count 2	Reporter intensity count 3	Reporter intensity count 4	Reporter intensity count 5	Reporter intensity count 6	Reporter intensity count 7	Reporter intensity count 8	Reporter intensity count 9	Reporter PIF	Reporter fraction	Reverse	Potential contaminant	id	Protein group IDs	Peptide ID	Mod. peptide ID	MS/MS IDs	Best MS/MS	AIF MS/MS IDs	Oxidation (M) site IDs
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMaxquant(String filename) throws IOException
  {
    TabFile file = new TabFile(filename, TabFile.tabb);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = file.get("Raw file"), scan = file.get("MS/MS Scan Number");
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("m/z"));
        id.setPrecursorNeutralMass(file.getDouble("Mass"));
//        id.setPrecursorIntensity(file.getDouble("Intensity"));
        id.setSpectrumFile(run);
        id.addRetentionTime(file.getDouble("Retention time"), TimeUnit.MINUTE);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = new PeptideMatch(file.get("Annotated Sequence").toUpperCase());
      String[] mods = Strs.split(file.get("Modifications"), ';');
      if (Tools.isSet(mods))
        for (String mod : mods)
        {
          String tag = mod.substring(0, mod.indexOf('(')).trim(), m = mod.substring(mod.indexOf('(')+1, mod.indexOf(')'));
          if      (Strs.equals(tag, "N-Term")) match.addModificationMatch(ModAttachment.N_TERM, new Modification(m, new NumericMass(0d)));
          else if (Strs.equals(tag, "C-Term")) match.addModificationMatch(ModAttachment.C_TERM, new Modification(m, new NumericMass(0d)));
          else
          {
            match.addModificationMatch(new Integer(tag.substring(1, tag.length()))-1, new Modification(m, new NumericMass(0d)));
          }
        }
      match.addProteinMatch(new PeptideProteinMatch(file.get("Master Protein Accessions"),
          Optional.of("-"), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("DeltaM [ppm]")*id.getPrecursorMz().get()*id.getAssumedCharge().get()*1E-6);
//      match.setNumMatchedIons(new Integer(matched[0]));
//      match.setTotalNumIons(new Integer(matched[1]));
      match.setNumMissedCleavages(file.getInt("# Missed Cleavages"));
      match.setRank(file.getInt("Rank"));
      Peptides.addScore(match, "XCorr",              file.getDouble("XCorr"));
      Peptides.addScore(match, "Percolator q-Value", file.getDouble("Percolator q-Value"));
      Peptides.addScore(match, "Percolator PEP",     file.getDouble("Percolator PEP"));
      Peptides.addScore(match, "DeltaScore",         file.getDouble("DeltaScore"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMSGFplus(String filename) throws IOException
  {
    TabFile file = new TabFile(filename, TabFile.tabb);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = Strs.stripLastOf(file.get("#SpecFile"), '.'), scan = file.get("ScanNum");
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("Precursor"));
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = Peptides.fromNumModSequence(file.get("Peptide"));
      // parse the protein names
      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
      String[] accs = Strs.split(file.get("Protein"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("XXX_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|');
          int       pre = acc.indexOf("pre=")+4, post = acc.indexOf("post=")+5;
          match.addProteinMatch(new PeptideProteinMatch(acs[1],
              Optional.of(acs[0]), Optional.of(acc.substring(pre,pre+1)), Optional.of(acc.substring(post,post+1)),
              decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
        }

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("PrecursorError(ppm)") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
      Peptides.addScore(match, "DeNovoScore",  file.getDouble("DeNovoScore"));
      Peptides.addScore(match, "IsotopeError", file.getInt(   "IsotopeError"));
      Peptides.addScore(match, "MSGFScore",    file.getDouble("MSGFScore"));
      Peptides.addScore(match, "SpecEValue",   file.getDouble("SpecEValue"));
      Peptides.addScore(match, "EValue",       file.getDouble("EValue"));
      Peptides.addScore(match, "QValue",       file.getDouble("QValue"));
      Peptides.addScore(match, "PepQValue",    file.getDouble("PepQValue"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }

    return id_match;
  }
  // TODO. Need more work in the future
  public static Dataframe fetchMzID(String root, double q, boolean save_decoy)
  {
    PsmReaders readers = new PsmReaders();

    List<String> files = IOs.listFiles(root, new FileFilterExtensions(new String[]{"mzid"}));
    if (Tools.isSet(files))
      for (String f : files)
      {
        System.out.println("Reading " + f);
        readers.parseMSGFplusMzID(f);
        // purge the matches
        if (!save_decoy || q<1)
        {
          System.out.print("Purging the low quality matches: " + readers.getResultMap().keySet().size());
          Iterator<SpectrumIdentifier> idr = readers.getResultMap().keySet().iterator();
          while (idr.hasNext())
          {
            List<PeptideMatch> matches = readers.getResultMap().get(idr.next());
            for (PeptideMatch match : matches)
            {
              if (match.getScore("MS-GF:QValue")>q ||
                  match.getProteinMatches().get(0).getHitType().equals(PeptideProteinMatch.HitType.DECOY)==save_decoy)
                idr.remove();

              break;
            }
          }
          System.out.print(" --> " + readers.getResultMap().keySet().size() + "\n");
        }
        break;
//        System.gc();
      }

    System.out.println("Preparing the dataframe");

    // output the PSMs in PD format
    Dataframe            df = new Dataframe("PSMs from " + root);
    Map<String, Object> row = new HashMap<>();
    long nrow = 0;
    for (SpectrumIdentifier id : readers.getResultMap().keySet())
    {
      // "Peptide","Mods","Accession","Charge","Rank","mz","ppm","Interference","InjectTime","RT","Scan1","Run","Score","q.val","Quan.Usage", channels
      List<PeptideMatch> matches = readers.getResultMap().get(id);
      for (PeptideMatch match : matches)
      {
        row.clear();
        row.put("Peptide", Peptides.toNumModSequence(match));
        Map<String, String> accs = Strs.toStrMap(Strs.split(match.getProteinMatches().get(0).getAccession(), '|'));
        for (String acc : accs.keySet())
          if (acc.indexOf("gi")>=0) row.put("Accession", accs.get(acc));
        if (row.get("Accession")==null) row.put("Accession",match.getProteinMatches().get(0).getAccession());
        row.put("Charge", id.getAssumedCharge().get());
        row.put("Rank", match.getRank());
        if (id.getPrecursorMz().isPresent()) row.put("mz", id.getPrecursorMz());
        row.put("ppm", 1E6 * match.getMassDiff() / id.getPrecursorNeutralMass().get());
        if (Tools.isSet(id.getRetentionTimes())) row.put("RT",    id.getRetentionTimes().getFirst().getTime());
        if (Tools.isSet(id.getScanNumbers()))    row.put("Scan1", id.getScanNumbers().getFirst().getValue());
        row.put("Run", Strs.stripLastOf(id.getSpectrumFile().get(), '.'));
        row.put("Score", match.getScore("MS-GF:RawScore"));
        row.put("q.val", match.getScore("MS-GF:QValue"));
        row.put("Decoy", match.getProteinMatches().get(0).getHitType().equals(PeptideProteinMatch.HitType.DECOY) ? 1:0);
        row.put("EValue", match.getScore("MS-GF:SpecEValue"));

        df.addRow((++nrow) + "", row);
        break; // TODO, keep the top one for now. Will test for homology next
      }
      if (nrow %5000==0) System.out.print(".");
    }
    df.init();

    return df;
  }
}
