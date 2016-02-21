package org.ms2ms.io;

import com.compomics.thermo_msf_parser_API.lowmeminstance.controllers.PeptideLowMemController;
import com.compomics.thermo_msf_parser_API.lowmeminstance.model.MsfFile;
import com.compomics.thermo_msf_parser_API.lowmeminstance.model.PeptideLowMem;
import com.google.common.base.Optional;
import com.google.common.collect.*;
import org.apache.commons.io.FilenameUtils;
import org.expasy.mzjava.core.mol.NumericMass;
import org.expasy.mzjava.core.ms.spectrum.TimeUnit;
import org.expasy.mzjava.proteomics.io.ms.ident.PSMReaderCallback;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodManager;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodMod;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodModificationResolver;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.algo.LCMSMS;
import org.ms2ms.algo.PSMs;
import org.ms2ms.algo.Peptides;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.data.ms.PeptideFeature;
import org.ms2ms.math.Stats;
import org.ms2ms.mzjava.*;
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
  private Engine mEngine     = null;
  private int    mLowestRank = 0;
  private PSMs.DesendScorePeptideMatch mSorter = null;

  private final Map<SpectrumIdentifier, List<PeptideMatch>> searchResultMap = new HashMap<>();
  // TODO should test the approach with a fixed size queue. The problem is with the default comparator of PeptideMatch
  // it's based the rank which is an integer
  // private final Map<SpectrumIdentifier, MinMaxPriorityQueue<PeptideMatch>> searchResultMap = new HashMap<>();
  private final Multimap<SpectrumIdentifier, PeptideMatch>  mIdMatch = HashMultimap.create();

  PSMReaderCallback insertIdResultCB = new PSMReaderCallback() {
    @Override
    public void resultRead(SpectrumIdentifier identifier, PeptideMatch searchResult)
    {
      List<PeptideMatch> results;
      if (searchResultMap.containsKey(identifier)) {
        results = searchResultMap.get(identifier);
      } else {
        results = new ArrayList<>();
        searchResultMap.put(identifier, results);
      }

      if (mLowestRank>0 && mSorter!=null && results.size()>=mLowestRank)
      {
        // let's see if the new one is better than the last one in the list
        if (results.get(results.size()-1).getScore(mEngine.getCanonicalScore()) <
                             searchResult.getScore(mEngine.getCanonicalScore()))
        {
          results.remove(results.size()-1);
          results.add(searchResult);
          // sort the list if it has multiple entries
          if (results.size()>1) Collections.sort(results, mSorter);
        }
      }
      else results.add(searchResult);
    }
  };
  public PSMReaderCallback insertIdMatch = new PSMReaderCallback() {
    @Override
    public void resultRead(SpectrumIdentifier identifier, PeptideMatch searchResult)
    {
      mIdMatch.put(identifier, searchResult);
    }
  };

  public PsmReaders() { super(); }
  public PsmReaders(Engine engine, int tops)
  {
    super();
    mEngine=engine; mLowestRank=tops;
    mSorter = new PSMs.DesendScorePeptideMatch(engine.getCanonicalScore());
  }

  public Map<SpectrumIdentifier, List<PeptideMatch>> getResultMap() { return searchResultMap; }

  public void parseMzID(String filename)
  {
    MzIdentMLReader reader = new MzIdentMLReader(new NumModResolver());
//    MzIdentMlReader reader = new MzIdentMlReader(new NumModResolver());

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
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMascotPD(String filename, String... quant_cols) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, TabFile.comma);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = Strs.stripLastOf(file.get("Spectrum File"), '.'), scan = file.get("First Scan");
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

      String[] mods = Strs.split(file.get("Modifications"), ';');
      PeptideMatch match = Tools.isSet(quant_cols) ?
        new PeptideFeature(file.get("Sequence").toUpperCase()) :
        new PeptideMatch(  file.get("Sequence").toUpperCase());

      if (Tools.isSet(mods))
        for (String mod : mods)
        {
          String tag = mod.substring(0, mod.indexOf('(')).trim(), m = mod.substring(mod.indexOf('(') + 1, mod.indexOf(')'));
          if      (Strs.equals(tag, "N-Term")) match.addModificationMatch(ModAttachment.N_TERM, new Modification(m, new NumericMass(0d)));
          else if (Strs.equals(tag, "C-Term")) match.addModificationMatch(ModAttachment.C_TERM, new Modification(m, new NumericMass(0d)));
          else
          {
            match.addModificationMatch(new Integer(tag.substring(1, tag.length()))-1, new Modification(m, new NumericMass(0d)));
          }
        }
      match.addProteinMatch(new PeptideProteinMatch(Strs.split(file.get("Protein Group Accessions"),';')[0],
          Optional.of(file.get("Protein Descriptions")), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("?M [ppm]")*id.getPrecursorMz().get()*id.getAssumedCharge().get()*1E-6);
      match.setNumMatchedIons(new Integer(matched[0]));
      match.setTotalNumIons(new Integer(matched[1]));
      match.setNumMissedCleavages(file.getInt("# Missed Cleavages"));
      match.setRank(file.getInt("Rank"));
      match.addScore("IonScore", file.getDouble("IonScore"));
      match.addScore("Exp Value", file.getDouble("Exp Value"));

      if (Tools.isSet(quant_cols))
      {
        PeptideFeature feature = (PeptideFeature )match;
        for (String col : quant_cols)
          if (file.getDouble(col)!=null) feature.setAbundance(col, file.getDouble(col));
      }
      // skip the row if quant cols are specified and "Quan Usage" is not "Used"
      if (!Tools.isSet(quant_cols) || Strs.equals("Used", file.get("Quan Usage"))) id_match.put(id, match);
//      if (!id_match.put(id, match))
//      {
////        System.out.println("Duplicated?");
//      }
    }

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readSequestHT(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, TabFile.tabb);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = Strs.stripLastOf(file.get("Spectrum File"), '.'), scan = file.get("First Scan");
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
      String[]      mods = Strs.split(file.get("Modifications"), ';');
      if (Tools.isSet(mods))
        for (String mod : mods)
        {
          if (!Strs.isSet(mod) || mod.indexOf('(')<0 || mod.indexOf('(')<0) continue;
          String tag = mod.substring(0, mod.indexOf('(')).trim(), m = mod.substring(mod.indexOf('(') + 1, mod.indexOf(')'));
          if      (Strs.equals(tag, "N-Term")) match.addModificationMatch(ModAttachment.N_TERM, new Modification(m, new NumericMass(0d)));
          else if (Strs.equals(tag, "C-Term")) match.addModificationMatch(ModAttachment.C_TERM, new Modification(m, new NumericMass(0d)));
          else
          {
            match.addModificationMatch(new Integer(tag.substring(1, tag.length()))-1, new Modification(m, new NumericMass(0d)));
          }
        }
      match.addProteinMatch(new PeptideProteinMatch(file.getNotNull("Master Protein Accessions", "Protein Accessions"),
        Optional.of("-"), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("DeltaM [ppm]") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
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
//  Peptide	-10lgP	Mass	Length	ppm	m/z	RT	Area	Scan	Accession	PTM	AScore
//  S(+229.16)LGPNSC(+57.02)SAN(+.98)GPGLYLIHGPNLYC(+57.02)YSDVEK(+229.16)	149	3740.822	30	2.2	1247.9507	56.24	7.71E+07	16073	gi|11321561:gi|530395269	TMT6plex; Carbamidomethylation; Deamidation (NQ)	S1:TMT6plex:1000.00;C7:Carbamidomethylation:1000.00;N10:Deamidation (NQ):76.17;C24:Carbamidomethylation:1000.00;K30:TMT6plex:1000.00
  public static Multimap<SpectrumIdentifier, PeptideMatch> readPEAKS(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, TabFile.comma);
    String   run = FilenameUtils.getBaseName(filename);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String           scan = file.get("Scan");
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(Math.round(file.getFloat("Mass") / file.getFloat("m/z")));
        id.setPrecursorMz(file.getDouble("m/z"));
        if (file.getDouble("Area")!=null) id.setPrecursorIntensity(file.getDouble("Area"));
        id.setSpectrumFile(run);
        id.addRetentionTime(file.getDouble("RT"), TimeUnit.MINUTE);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      // Peptide: S(+229.16)LGPNSC(+57.02)SAN(+.98)GPGLYLIHGPNLYC(+57.02)YSDVEK(+229.16)
      //     PTM: TMT6plex; Carbamidomethylation; Deamidation (NQ)
      //  AScore: S1:TMT6plex:1000.00;C7:Carbamidomethylation:1000.00;N10:Deamidation (NQ):76.17;C24:Carbamidomethylation:1000.00;K30:TMT6plex:1000.00
      PeptideMatch match = new PeptideMatch(Peptides.keepAAs(file.get("Peptide")));
      String[]      mods = Strs.split(file.get("PTM"), ';', true), ascores = Strs.split(file.get("AScore"), ';', true);

      UnimodModificationResolver unimod = new UnimodModificationResolver();
      unimod.putTranslate("Carbamidomethylation",             "Carbamidomethyl");
      unimod.putTranslate("Deamidation (NQ)",                 "Deamidated");
      unimod.putTranslate("Oxidation (M)",                    "Oxidation");
      unimod.putTranslate("Oxidation (HW)",                   "Oxidation");
      unimod.putTranslate("Dihydroxy",                        "Dioxidation");
      unimod.putTranslate("Phosphorylation (STY)",            "Phospho");
      unimod.putTranslate("Sulphone",                         "Sulfo");
      unimod.putTranslate("Dehydration",                      "Dehydrated");
      unimod.putTranslate("Tryptophan oxidation to kynurenin","Trp->Kynurenin");
      unimod.putTranslate("Methylmalonylation on Serine",     "Methylmalonylation"); // intrim name
      unimod.putTranslate("Replacement of proton with ammonium ion","Ammonium");
      unimod.putTranslate("Methyl ester",                     "Methyl");
      unimod.putTranslate("2-amino-3-oxo-butanoic_acid",      "Didehydro");
      unimod.putTranslate("Hex1HexNAc1NeuAc1",                "Hex(1)HexNAc(1)NeuAc(1)");
      unimod.putTranslate("Deamidation followed by a methylation", "Methyl+Deamidated");
      unimod.putTranslate("Tyrosine oxidation to 2-aminotyrosine", "Amino");
      unimod.putTranslate("Fluorination",                      "Fluoro");
      unimod.putTranslate("Pyro-glu from E",                   "Glu->pyro-Glu");
      unimod.putTranslate("Pyro-glu from Q",                   "Gln->pyro-Glu");
      unimod.putTranslate("Methylation",                       "Methyl");
      unimod.putTranslate("Formylation (Protein N-term)",      "Formyl");
      unimod.putTranslate("Formylation",                       "Formyl");
      unimod.putTranslate("Carbamidomethylation (DHKE  X@N-term)","Carbamidomethyl");
      unimod.putTranslate("Ammonia-loss (N)",                  "Ammonia-loss");
      unimod.putTranslate("2 4-diacetamido-2 4 6-trideoxyglucopyranose","Bacillosamine");
      unimod.putTranslate("Biotinylation",                     "Biotin");
      unimod.putTranslate("Phosphorylation (HCDR)",            "Phospho");
      unimod.putTranslate("Sodium adduct",                     "Cation:Na");
      unimod.putTranslate("Amidation",                         "Amidated");
      unimod.putTranslate("Persulfide",                        "Sulfide");
      unimod.putTranslate("Deamidation (R)",                   "Deamidated");
      unimod.putTranslate("2-OH-ethyl thio-Ser",               "MercaptoEthanol"); // intrim name
      unimod.putTranslate("Aminoethylbenzenesulfonylation",    "AEBS");
      unimod.putTranslate("HexNAcylation (ST)",                "HexNAc");
      unimod.putTranslate("Acetylation (N-term)",              "Acetyl");
      unimod.putTranslate("3-sulfanylpropanoyl",               "Thioacyl");
      unimod.putTranslate("Proline oxidation to pyroglutamic acid","Pro->pyro-Glu");
      unimod.putTranslate("Methylphosphonylation",             "Methylphosphonate"); // intrim name
      unimod.putTranslate("Glycosyl-L-hydroxyproline",         "Glycosyl");
      unimod.putTranslate("Tryptophan oxidation to oxolactone","Trp->Oxolactone");
      unimod.putTranslate("Replacement of proton by lithium",  "Cation:Li");
      unimod.putTranslate("S-Ethylcystine from Serine",        "Delta:H(4)C(2)O(-1)S(1)");
      unimod.putTranslate("Carbamylation",                     "Carbamyl");
      unimod.putTranslate("Hexose (NSY)",                      "Hex");
      unimod.putTranslate("Nucleophilic addition to cytopiloyne+H2O","Cytopiloyne+water");
      unimod.putTranslate("Ethyl amino",                       "ethylamino"); // intrim name
      unimod.putTranslate("Replacement of 2 protons by calcium","Cation:Ca[II]"); // intrim name
      unimod.putTranslate("Flavin mononucleotide",             "FMNH");
      unimod.putTranslate("Glycerylphosphorylethanolamine",    "GlycerylPE");
      unimod.putTranslate("Dithiothreitol (DTT)",              "DTT_C"); // intrim name
      unimod.putTranslate("Hexosamine",                        "HexN");
      unimod.putTranslate("Phosphorylation to pyridyl thiol",  "PET");
      unimod.putTranslate("4-hydroxynonenal (HNE)",            "HNE");
      unimod.putTranslate("Tryptophan oxidation to hydroxykynurenin","Trp->Hydroxykynurenin");
      unimod.putTranslate("Acetaldehyde +26",                  "Delta:H(2)C(2)");
      unimod.putTranslate("Carboxylation (E)","Carboxy");
      unimod.putTranslate("N-glucuronylation","Glucuronyl");
      unimod.putTranslate("Triglutamyl","GluGluGlu");
      unimod.putTranslate("Deamidation followed by esterification with ethanol","Ethyl+Deamidated"); // intrim name
      unimod.putTranslate("Tetraglutamyl","GluGluGluGlu");
      unimod.putTranslate("Beta-methylthiolation (ND)","Methylthio");
      unimod.putTranslate("Ethanolation (KR)","Ethanolyl");
      unimod.putTranslate("Propionaldehyde +40","Delta:H(4)C(3)");
      unimod.putTranslate("Dihydroxy methylglyoxal adduct","Dihydroxyimidazolidine"); // intrim name
      unimod.putTranslate("Aminoethylcysteine","AEC-MAEC"); // intrim name
      unimod.putTranslate("Fucose","dHex");
      unimod.putTranslate("Dichlorination of tyrosine residues","dichlorination"); // intrim name
      unimod.putTranslate("Nitroalkylation by Nitro Oleic Acid","NA-LNO2"); // intrim name
      unimod.putTranslate("Chlorination of tyrosine residues","Chlorination"); // intrim name
      unimod.putTranslate("Acetylation (TSCYH)","Acetyl");
      unimod.putTranslate("Carboxylation (DKW)","Carboxy");
      unimod.putTranslate("Carboxylation (E)","Carboxy");
      unimod.putTranslate("2 3-dihydro-2 2-dimethyl-7-benzofuranol N-methyl carbamate","Carbofuran"); // intrim name
      unimod.putTranslate("Ammonia-loss (C@N-term)","Ammonia-loss");
      unimod.putTranslate("Dihydroxy methylglyoxal adduct","Dihydroxyimidazolidine"); // intrim name
      unimod.putTranslate("Acetylation (Protein N-term)","Acetyl");
      unimod.putTranslate("Phosphorylation to amine thiol","DAET");
      unimod.putTranslate("O-Pinacolylmethylphosphonylation","O-pinacolylmethylphosphonate"); // intrim name
      unimod.putTranslate("HexNAcylation (N)","HexNAc");
      unimod.putTranslate("Pyrrolidone from Proline","Pro->Pyrrolidone");
      unimod.putTranslate("Replacement of 2 protons by nickel","Cation:Ni[II]");  // intrim name
      unimod.putTranslate("Biantennary","Hex(5)HexNAc(4)");
      unimod.putTranslate("Proline oxidation to pyrrolidinone","Pro->Pyrrolidone");
//      unimod.putTranslate("","");

      // parse the mutation
      // Peptide: V(+229.16)INLPQ(sub L)DSMAAPWETGDTFPDVVAIAPDVR
      //     PTM: TMT6plex; Mutation
      // AScores: V1:TMT6plex:1000.00
      if (Tools.isSet(ascores))
        // S1:TMT6plex:1000.00
        for (String mod : ascores)
        {
          if (!Strs.isSet(mod)) continue;
          String[] items = Strs.split(mod, ':');
//          String AA = items[0].substring(0, 1);
//          double ascore = Stats.toDouble(items[2]);
          Optional<Modification> M = unimod.resolve(items[1]);
          if (M.isPresent())
          {
            match.addModificationMatch(Stats.toInt(items[0].substring(1))-1, M.get());
          }
          else
          {
            System.out.println(items[1] + " can not be resolved by UniMod");
          }
//          if      (Strs.equals(tag, "N-Term")) match.addModificationMatch(ModAttachment.N_TERM, new Modification(m, new NumericMass(0d)));
//          else if (Strs.equals(tag, "C-Term")) match.addModificationMatch(ModAttachment.C_TERM, new Modification(m, new NumericMass(0d)));
        }
      // set the protein info
      match.addProteinMatch(new PeptideProteinMatch(file.get("Accession"),
        Optional.of("-"), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

      match.setMassDiff(file.getDouble("ppm") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
      match.addScore("PeakScore", file.getDouble("-10lgP"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }
    return id_match;
  }
  // TODO to be completed
  //allPeptide//Raw file	Type	Charge	m/z	Mass	Uncalibrated m/z	Resolution	Number of data points	Number of scans	Number of isotopic peaks	PIF	Mass fractional part	Mass deficit	Mass precision [ppm]	Max intensity m/z 0	Retention time	Retention length	Retention length (FWHM)	Min scan number	Max scan number	Identified	MS/MS IDs	Sequence	Length	Modifications	Modified sequence	Proteins	Score	Intensity	Intensities	MS/MS Count	MSMS Scan Numbers	MSMS Isotope Indices
  //evidence//Sequence	Length	Modifications	Modified sequence	Oxidation (M) Probabilities	Oxidation (M) Score Diffs	Acetyl (Protein N-term)	Oxidation (M)	Missed cleavages	Proteins	Leading proteins	Leading razor protein	Type	Raw file	Fraction	Experiment	MS/MS m/z	Charge	m/z	Mass	Resolution	Uncalibrated - Calibrated m/z [ppm]	Uncalibrated - Calibrated m/z [Da]	Mass Error [ppm]	Mass Error [Da]	Uncalibrated Mass Error [ppm]	Uncalibrated Mass Error [Da]	Max intensity m/z 0	Retention time	Retention length	Calibrated retention time	Calibrated retention time start	Calibrated retention time finish	Retention time calibration	Match time difference	Match m/z difference	Match q-value	Match score	Number of data points	Number of scans	Number of isotopic peaks	PIF	Fraction of total spectrum	Base peak fraction	PEP	MS/MS Count	MS/MS Scan Number	Score	Delta score	Combinatorics	Intensity	Reporter intensity 0	Reporter intensity 1	Reporter intensity 2	Reporter intensity 3	Reporter intensity 4	Reporter intensity 5	Reporter intensity 6	Reporter intensity 7	Reporter intensity 8	Reporter intensity 9	Reporter intensity not corrected 0	Reporter intensity not corrected 1	Reporter intensity not corrected 2	Reporter intensity not corrected 3	Reporter intensity not corrected 4	Reporter intensity not corrected 5	Reporter intensity not corrected 6	Reporter intensity not corrected 7	Reporter intensity not corrected 8	Reporter intensity not corrected 9	Reporter intensity count 0	Reporter intensity count 1	Reporter intensity count 2	Reporter intensity count 3	Reporter intensity count 4	Reporter intensity count 5	Reporter intensity count 6	Reporter intensity count 7	Reporter intensity count 8	Reporter intensity count 9	Reporter PIF	Reporter fraction	Reverse	Potential contaminant	id	Protein group IDs	Peptide ID	Mod. peptide ID	MS/MS IDs	Best MS/MS	AIF MS/MS IDs	Oxidation (M) site IDs
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMaxquant(String mq) throws IOException
  {
    UnimodModificationResolver modResolver = new UnimodModificationResolver();
    modResolver.putTranslate("de", "Deamidated");
    modResolver.putTranslate("ox", "Oxidation");
    modResolver.putTranslate("ac", "Acetyl");
    modResolver.putTranslate("gl", "Gln->pyro-Glu");

    MaxQuantReader maxQuantReader = new MaxQuantReader(modResolver);
    PsmReaders                psm = new PsmReaders();
    maxQuantReader.parse(new File(mq), psm.insertIdMatch);

    System.out.println(psm.mIdMatch.size());

    return psm.mIdMatch;
  }
/*
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMaxquant(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, TabFile.tabb);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = Strs.stripLastOf(file.get("Raw file"), '.'), scan = file.get("MS/MS Scan Number");
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
      PeptideMatch match = new PeptideMatch(file.get("Sequence").toUpperCase());
      String[]      mods = Strs.split(file.get("Modifications"), ';');

//      Oxidation (M)                         2 Oxidation (M)                       Acetyl (Protein N-term),Oxidation (M)
//      Acetyl (Protein N-term)               Gln->pyro-Glu                         Oxidation (M),Gln->pyro-Glu
//      3 Oxidation (M)
      // TODO to work out the mod parsing code
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
      match.addProteinMatch(new PeptideProteinMatch(file.get("Leading proteins"),
          Optional.of("-"), Optional.of("-"), Optional.of("-"), PeptideProteinMatch.HitType.TARGET));

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("DeltaM [ppm]")*id.getPrecursorMz().get()*id.getAssumedCharge().get()*1E-6);
//      match.setNumMatchedIons(new Integer(matched[0]));
//      match.setTotalNumIons(new Integer(matched[1]));
      match.setNumMissedCleavages(file.getInt("# Missed Cleavages"));
      match.setRank(file.getInt("Rank"));
      PSMs.addScore(match, "XCorr", file.getDouble("XCorr"));
      PSMs.addScore(match, "Percolator q-Value", file.getDouble("Percolator q-Value"));
      PSMs.addScore(match, "Percolator PEP",     file.getDouble("Percolator PEP"));
      PSMs.addScore(match, "DeltaScore",         file.getDouble("DeltaScore"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }

    return id_match;
  }
*/
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMSGFplus(Collection<String> files) throws IOException
  {
    if (Tools.isSet(files))
    {
      // the totals
      Multimap<SpectrumIdentifier, PeptideMatch> matches = HashMultimap.create();
      for (String f : files)
      {
        matches.putAll(readMSGFplus(f, '\t', 0));
      }
      return matches;
    }
    return null;
  }
//  Scan Number     Title   Sequence        Modifications   Protein Accessions      Amanda Score    Weighted Probability    Rank    m/z     Charge  RT      Filename
//  20130510_EXQ1_IgPa_QC_UPS1_01.34.34.4   mFVcSDTDYcRQQSEAKNQ     M1(Oxidation|15.994915|variable);C4(Carbamidomethyl|57.021464|fixed);C10(Carbamidomethyl|57.021464|fixed)       gi|253775272    3.42563717374274        0.454397865483082       1       596.5   4       11.2614798      20130510_EXQ1_IgPa_QC_UPS1_01.mgf
// TODO to be completed
  public static Multimap<SpectrumIdentifier, PeptideMatch> readAmanda(String filename, int lowest_rank) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, TabFile.tabb, 1); int counts=0;
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      // 20130510_EXQ1_IgPa_QC_UPS1_01.34.34.4
      String[] items = Strs.split(file.get("Title"), '.');
      String     run = Strs.stripLastOf(file.get("Filename"), '.'), scan = items[items.length-3];
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("m/z"));
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = new PeptideMatch(file.get("Sequence").toUpperCase());
      // grab the rank up-front so we can skip the rest if needed
      match.setRank(Stats.toInt(file.get("Rank")));
      if (match.getRank()>lowest_rank) continue;

      // Sequence	Modifications
      // mFVcSDTDYcRQQSEAKNQ	M1(Oxidation|15.994915|variable);C4(Carbamidomethyl|57.021464|fixed);C10(Carbamidomethyl|57.021464|fixed)
      String[] mods = Strs.isSet(file.get("Modifications"))?Strs.split(file.get("Modifications"), ';'):null;
      if (Tools.isSet(mods))
        for (String mod : mods)
        {
          int      left=mod.indexOf('('), right=mod.indexOf(')', left); Integer pos=-1;
          String      r=mod.substring(0, 1), tag=mod.substring(0, left);
          try
          {
            String[]  def=Strs.split(mod.substring(left + 1, right), '|');
            double increM=Stats.toDouble(def[1]);

            // TODO not optimal for terminal cases
            if      (Strs.equals(tag,"N-Term")) match.addModificationMatch(ModAttachment.N_TERM, increM);
            else if (Strs.equals(tag,"C-Term")) match.addModificationMatch(ModAttachment.C_TERM, increM);
            else
            {
              pos = Stats.toInt(tag.substring(1));
              //if (pos==null) pos = Stats.toInt(tag);
              if (pos==null)
                System.out.println();
              else
                match.addModificationMatch(pos-1, increM);
            }
          }
          catch (Exception e)
          {
            e.printStackTrace();
          }
        }

      // parse the protein names: gi|253775272
      String[] accs = Strs.split(file.get("Protein Accessions"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("REV_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|');
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
              Optional.of(acs[0]), Optional.of("-"), Optional.of("-"),
              decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
          }
          catch (ArrayIndexOutOfBoundsException e)
          {
            System.out.println(e.toString());
          }
        }

      match.setNeutralPeptideMass(match.toPeptide(PSMs.sNumModResolver).getMolecularMass());
      match.setMassDiff(id.getPrecursorNeutralMass().get() - match.getNeutralPeptideMass());

      PSMs.addScore(match, "AmandaScore", file.getDouble("Amanda Score"));
      PSMs.addScore(match, "WeightedProbability", file.getDouble("Weighted Probability"));

      // add as much information to the match as possible
      PSMs.addScore(match, "^Charge", file.getInt("Charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
      if (++counts%50000==0) System.out.print(".");
    }

    return id_match;
  }

  //  #SpecFile       SpecID  ScanNum FragMethod      Precursor       IsotopeError    PrecursorError(ppm)     Charge  Peptide Protein DeNovoScore     MSGFScore       SpecEValue      EValue  QValue  PepQValue
  //  20130510_EXQ1_IgPa_QC_UPS1_01.mzML      controllerType=0 controllerNumber=1 scan=61474  61474   HCD     1101.2928       0       10.197636       4       QVDPAALTVHYVTALGTDSFSQQMLDAWHGENVDTSLTQR        gi|253771643|ref|YP_003034474.1|(pre=R,post=M)  355     335     1.5431253692232283E-44  4.0382819349887276E-38  0.0     0.0
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMSGFplus(String filename, char delimiter, int lowest_rank) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, delimiter); int counts=0;
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
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = PSMs.fromNumModSequence(file.get("Peptide"));
      // parse the protein names
      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
      String[] accs = Strs.split(file.get("Protein"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("XXX_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|');
          int       pre = acc.indexOf("pre=")+4, post = acc.indexOf("post=")+5;
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
                Optional.of(acs[0]), Optional.of(acc.substring(pre,pre+1)), Optional.of(acc.substring(post,post+1)),
                decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
            // try to min the storage cost, WYU 20160216
            if (match.getProteinMatches() instanceof ArrayList)
              ((ArrayList )match.getProteinMatches()).trimToSize();
          }
          catch (ArrayIndexOutOfBoundsException e)
          {
            System.out.println(e.toString());
          }
        }

      if (file.getDouble("PrecursorError(ppm)")!=null)
           match.setMassDiff(file.getDouble("PrecursorError(ppm)") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
      else match.setMassDiff(file.getDouble("PrecursorError(Da)"));

      PSMs.addScore(match, "DeNovoScore",  file.getDouble("DeNovoScore"));
      PSMs.addScore(match, "IsotopeError", file.getInt("IsotopeError"));
      PSMs.addScore(match, "MSGFScore",    file.getDouble("MSGFScore"));
      PSMs.addScore(match, "SpecEValue",   file.getDouble("SpecEValue"));
      PSMs.addScore(match, "EValue",       file.getDouble("EValue"));
      PSMs.addScore(match, "QValue",       file.getDouble("QValue"));
      PSMs.addScore(match, "PepQValue",    file.getDouble("PepQValue"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
      PSMs.addScore(match, "^Charge",      file.getInt("Charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }

      if (++counts%5000==0) System.out.print(".");
    }
    System.out.println();
    // make sure we clear the temp objects
    file.close(); file=null; Tools.dispose(run_scan_id);

    LCMSMS.rank(id_match, "MSGFScore", true, false);

    if (lowest_rank>0)
      for (SpectrumIdentifier id : id_match.keySet())
      {
        Iterator<PeptideMatch> itr = id_match.get(id).iterator();
        while (itr.hasNext())
          if (itr.next().getRank()>lowest_rank) itr.remove();
      }

    return id_match;
  }
  // CSV file converted from mzID output
  //  Raw data location,Spectrum ID,Spectrum Title,Retention Time (s),PSM_ID,rank,Pass Threshold,Calc m/z,Exp m/z,Charge,Sequence,Modifications,number of matched peaks,number of unmatched peaks,MyriMatch:MVH,MyriMatch:mzFidelity,proteinacc_start_stop_pre_post_;,Is decoy
  // /Biomarker/yuw/data/Joslin/Mouse_Plasma_LIRKO2_01_27Aug12_Lynx_12-06-05.mzML,"controllerType=0 controllerNumber=1 scan=10","",0.00302705555555,"SIR_1_SII_1",1,true,497.735128548155,497.7578337711,4,"KLAQCYQCDELHHIVM",Carbamidomethyl:0,6,82,26.45690389115,20.745448124535,"rev_gi|568977859|ref|XP_006515205.1|Agr3$_89_105_K_-;rev_gi|568977861|ref|XP_006515206.1|Agr3$_89_105_K_-",true
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMyriMatch(String filename, int lowest_rank) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, ','); int counts=0;
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    String[] scores = {"number of matched peaks", "number of unmatched peaks", "MyriMatch:MVH", "MyriMatch:mzFidelity"};

    while (file.hasNext())
    {
      // controllerType=0 controllerNumber=1 scan=3
      String run = new File(file.get("Raw data location")).getName(),
            scan = Tools.back(Strs.split(file.get("Spectrum ID"),'=',true));
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("Exp m/z"));
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }

      PeptideMatch match = new PeptideMatch(file.get("Sequence"));

      // Carbamidomethyl:0;iTRAQ4plex114:14;iTRAQ4plex114:15
      // DVLTLQLEVLMETDSRLHFKIK	iTRAQ4plex114:0;Oxidation:11;iTRAQ4plex114:20;iTRAQ4plex114:22
      if (Strs.isSet(file.get("Modifications")))
      {
        String[] items = Strs.split(file.get("Modifications"), ';', true);
        for (String mod : items)
        {
          String[] tags = Strs.split(mod, ':');
          int       pos = Integer.valueOf(tags[1]);
          Optional<UnimodMod> M = UnimodManager.getModification(tags[0]);

          if (M.isPresent())
            if (pos==0) match.addModificationMatch(ModAttachment.N_TERM, M.get());
            else        match.addModificationMatch(pos-1, M.get());
        }
      }

      // parse the protein names
      // gi|226693367|ref|NP_032090.3|Gaa$_163_185_K_D;gi|226693369|ref|NP_001152796.1|Gaa$_163_185_K_D
      String[] accs = Strs.split(file.get("proteinacc_start_stop_pre_post_;"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("rev_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|'), stops = Strs.split(Tools.back(acs), '_');
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
                Optional.of(acs[0]), Optional.of(stops[stops.length-2]), Optional.of(stops[stops.length-1]),
                decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
            // try to min the storage cost, WYU 20160216
            if (match.getProteinMatches() instanceof ArrayList)
              ((ArrayList )match.getProteinMatches()).trimToSize();
          }
          catch (Exception e)
          {
            System.out.println(e.toString());
          }
        }

      match.setMassDiff((file.getDouble("Exp m/z") - file.getDouble("Calc m/z")) * id.getAssumedCharge().get());

      for (String score : scores)
        if (file.get(score)!=null) PSMs.addScore(match, score, file.getDouble(score));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
      PSMs.addScore(match, "^Charge",      file.getInt("Charge"));

      id_match.put(id, match);

      if (++counts%5000==0) System.out.print(".");
    }
    System.out.println();
    // make sure we clear the temp objects
    file.close(); Tools.dispose(run_scan_id);

    LCMSMS.rank(id_match, Engine.MYRIMATCH.getCanonicalScore(), true, false);

    if (lowest_rank>0)
      for (SpectrumIdentifier id : id_match.keySet())
      {
        Iterator<PeptideMatch> itr = id_match.get(id).iterator();
        while (itr.hasNext())
          if (itr.next().getRank()>lowest_rank) itr.remove();
      }

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readXTandem(String filename, int lowest_rank) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, ','); int counts=0;
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();

    while (file.hasNext())
    {
      // Mouse_Plasma_LIRKO2_01_27Aug12_Lynx_12-06-05.mzXML scan 6484 (charge 3)
      String run = new File(file.get("Raw data location")).getName(),
            scan = Strs.split(Strs.split(file.get("Spectrum ID"),"scan", true).get(1), ' ')[0];
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("Exp m/z"));
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = new PeptideMatch(file.get("Sequence"));

      if (Strs.isSet(file.get("Modifications")))
      {
        // KCCEDGMRDIPMR	:2;:3;:1;:1;Oxidation:7;Oxidation:12
        String[] items = Strs.split(file.get("Modifications"), ';', true);
        for (String mod : items)
        {
          String[] tags = Strs.split(mod, ':');
          int       pos = Integer.valueOf(tags[1]);
          Optional<UnimodMod> M = UnimodManager.getModification(tags[0]);

          if (M.isPresent())
            if (pos==0) match.addModificationMatch(ModAttachment.N_TERM, M.get());
            else        match.addModificationMatch(pos-1, M.get());
        }
      }

      // parse the protein names
      // gi|226693367|ref|NP_032090.3|Gaa$_163_185_K_D;gi|226693369|ref|NP_001152796.1|Gaa$_163_185_K_D
      String[] accs = Strs.split(file.get("proteinacc_start_stop_pre_post_;"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("rev_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|'), stops = Strs.split(Tools.back(acs), '_');
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
                Optional.of(acs[0]), Optional.of(stops[stops.length-2].substring(0,1)), Optional.of(stops[stops.length-1].substring(0,1)),
                decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
            // try to min the storage cost, WYU 20160216
            if (match.getProteinMatches() instanceof ArrayList)
              ((ArrayList )match.getProteinMatches()).trimToSize();
          }
          catch (Exception e)
          {
            System.out.println(e.toString());
          }
        }

      match.setMassDiff((file.getDouble("Exp m/z") - file.getDouble("Calc m/z")) * id.getAssumedCharge().get());

      PSMs.addScore(match, "X\\!Tandem:expect", -10d*Math.log10(file.getDouble("X\\!Tandem:expect")));
      PSMs.addScore(match, "X\\!Tandem:hyperscore", file.getDouble("X\\!Tandem:hyperscore"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
      PSMs.addScore(match, "^Charge",      file.getInt("Charge"));

      id_match.put(id, match);

      if (++counts%5000==0) System.out.print(".");
    }
    System.out.println();
    // make sure we clear the temp objects
    file.close(); Tools.dispose(run_scan_id);

    LCMSMS.rank(id_match, Engine.XTANDEM.getCanonicalScore(), true, false);

    if (lowest_rank>0)
      for (SpectrumIdentifier id : id_match.keySet())
      {
        Iterator<PeptideMatch> itr = id_match.get(id).iterator();
        while (itr.hasNext())
          if (itr.next().getRank()>lowest_rank) itr.remove();
      }

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMSGFplus_mzidlib(String filename, String delimiter) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    TabFile file = new TabFile(filename, delimiter);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String run = Strs.stripLastOf(file.get("Raw data location"), '.'), scan = Strs.split(file.get("Spectrum ID"), "scan=", true).get(1);
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorMz(file.getDouble("Exp m/z"));
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        id.addRetentionTime(file.getDouble("Retention Time (s)"), TimeUnit.SECOND);
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = PSMs.fromNumModSequence(file.get("Sequence"));
      // parse the protein names
      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
      String[] items = Strs.split(file.get("proteinacc_start_stop_pre_post_;"), '_'), accs = Strs.split(items[0], ';');
      Optional<String> pre = Optional.of(items[items.length-2]), post=Optional.of(items[items.length-1]);
      boolean  decoy = Strs.equals("true", file.get("is decoy"));
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|');
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
                Optional.of(acs[0]), pre, post,
                decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
          }
          catch (ArrayIndexOutOfBoundsException e)
          {
            System.out.println(e.toString());
          }
        }

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      match.setMassDiff(file.getDouble("Calc m/z")-file.getDouble("Exp m/z"));

      PSMs.addScore(match, "DeNovoScore",  file.getDouble("MS-GF:DeNovoScore"));
//      PSMs.addScore(match, "IsotopeError", file.getInt("IsotopeError"));
      PSMs.addScore(match, "MSGFScore",    file.getDouble("MS-GF:RawScore"));
      PSMs.addScore(match, "SpecEValue",   file.getDouble("MS-GF:SpecEValue"));
      PSMs.addScore(match, "EValue",       file.getDouble("MS-GF:EValue"));
      PSMs.addScore(match, "QValue",       file.getDouble("MS-GF:QValue"));
      PSMs.addScore(match, "PepQValue",    file.getDouble("MS-GF:PepQValue"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
      PSMs.addScore(match, "^Charge", file.getInt("Charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }
    LCMSMS.rank(id_match, "MSGFScore", true, false);

    return id_match;
  }
//  Filename        Spectrum Number Spectrum ID     Spectrum Title  Retention Time (minutes)        Precursor m/z   Precursor Intensity     Precursor Charge        Precursor Mass (Da)     Experimental Peaks      Total Intensity Peptide Sequence        Base Peptide Sequence   Protein Description     Start Residue Number    Stop Residue Number     Missed Cleavages        Theoretical Mass (Da)   Precursor Mass Error (Da)       Precursor Mass Error (ppm)      Matching Products       Total Products  Ratio of Matching Products      Matching Intensity      Fraction of Intensity Matching  Morpheus Score  Target? Decoy?  Cumulative Target       Cumulative Decoy        Q-Value (%)
//  X:\data\UPS_EColi_1\20130510_EXQ1_IgPa_QC_UPS1_01.mzML  67969   controllerType=0 controllerNumber=1 scan=67969          168.04589       1065.86999511719        7412378.5       3       3194.58815621156        316     17005826.3769531        R.TAGSSGANPFAC[carbamidomethylation of C]IAAGIASLWGPAHGGANEAALK.M       TAGSSGANPFACIAAGIASLWGPAHGGANEAALK      gi|253774310|ref|YP_003037141.1| citrate synthase I [Escherichia coli BL21(DE3)]        241     274     0       3194.5567307142 0.0314254973618517      9.83720121784345        45      66      0.681818181818182       4206046.93017578        0.247329758457134       45.2473297584571        True    False   1       0       0
// TODO to be completed
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMorpheus(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

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
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = PSMs.fromNumModSequence(file.get("Peptide"));
      // parse the protein names
      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
      String[] accs = Strs.split(file.get("Protein"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("XXX_")==0);
          String[]  acs = Strs.split(decoy ? acc.substring(4) : acc, '|');
          int       pre = acc.indexOf("pre=")+4, post = acc.indexOf("post=")+5;
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
              Optional.of(acs[0]), Optional.of(acc.substring(pre,pre+1)), Optional.of(acc.substring(post,post+1)),
              decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
          }
          catch (ArrayIndexOutOfBoundsException e)
          {
            System.out.println(e.toString());
          }
        }

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      if (file.getDouble("PrecursorError(ppm)")!=null)
        match.setMassDiff(file.getDouble("PrecursorError(ppm)") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
      else match.setMassDiff(file.getDouble("PrecursorError(Da)"));

      PSMs.addScore(match, "DeNovoScore",  file.getDouble("DeNovoScore"));
      PSMs.addScore(match, "IsotopeError", file.getInt("IsotopeError"));
      PSMs.addScore(match, "MSGFScore", file.getDouble("MSGFScore"));
      PSMs.addScore(match, "SpecEValue", file.getDouble("SpecEValue"));
      PSMs.addScore(match, "EValue", file.getDouble("EValue"));
      PSMs.addScore(match, "QValue", file.getDouble("QValue"));
      PSMs.addScore(match, "PepQValue", file.getDouble("PepQValue"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
      PSMs.addScore(match, "^Charge", file.getInt("Charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }
    LCMSMS.rank(id_match, "MSGFScore", true, false);

    return id_match;
  }
//  # id, scanNum, RT, mz(data), z, pepMass(denovo), err(data-denovo), ppm(1e6*err/(mz*z)), score, peptide, aaScore,
//  1, 0, 1.1, 419.3200, 2, 836.5848, 0.0407, 48.5, 0.1, PLLLLLR, 1-1-1-1-1-1-32
//  2, 0, 1.3, 596.7500, 4, 2383.0038, -0.0329, -13.8, 5.3, C(Cam)LC(Cam)MSSPDAWVSDRC(Cam)NRNR, 36-1-1-1-1-1-1-7-1-1-1-1-1-1-1-1-1-7-32
//  3, 0, 7.6, 1550.4600, 3, 4648.6858, -0.3276, -70.4, 0.0, HGC(Cam)C(Cam)C(Cam)C(Cam)DYKKLFGENRM(O)HC(Cam)C(Cam)C(Cam)NFGQC(Cam)C(Cam)C(Cam)C(Cam)GVTYM(O)K, 3-1-3-1-1-1-1-2-1-1-1-1-1-1-1-1-1-1-1-2-2-2-1-1-1-1-1-1-1-2-1-1-1-9-4
// TODO to be completed
  public static Multimap<SpectrumIdentifier, PeptideMatch> readNovor(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

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
        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = PSMs.fromNumModSequence(file.get("Peptide"));
      // parse the protein names
      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
      String[] accs = Strs.split(file.get("Protein"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("XXX_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|');
          int       pre = acc.indexOf("pre=")+4, post = acc.indexOf("post=")+5;
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
              Optional.of(acs[0]), Optional.of(acc.substring(pre,pre+1)), Optional.of(acc.substring(post,post+1)),
              decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
          }
          catch (ArrayIndexOutOfBoundsException e)
          {
            System.out.println(e.toString());
          }
        }

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      if (file.getDouble("PrecursorError(ppm)")!=null)
        match.setMassDiff(file.getDouble("PrecursorError(ppm)") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
      else match.setMassDiff(file.getDouble("PrecursorError(Da)"));

      PSMs.addScore(match, "DeNovoScore",  file.getDouble("DeNovoScore"));
      PSMs.addScore(match, "IsotopeError", file.getInt("IsotopeError"));
      PSMs.addScore(match, "MSGFScore",    file.getDouble("MSGFScore"));
      PSMs.addScore(match, "SpecEValue",   file.getDouble("SpecEValue"));
      PSMs.addScore(match, "EValue",       file.getDouble("EValue"));
      PSMs.addScore(match, "QValue", file.getDouble("QValue"));
      PSMs.addScore(match, "PepQValue", file.getDouble("PepQValue"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
      PSMs.addScore(match, "^Charge", file.getInt("Charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }
    LCMSMS.rank(id_match, "MSGFScore", true, false);

    return id_match;
  }
//  Spectrum number, Filename/id, Peptide, E-value, Mass, gi, Accession, Start, Stop, Defline, Mods, Charge, Theo Mass, P-value, NIST score
//  171,20130510_EXQ1_IgPa_QC_UPS1_01.1158.1158.2,HFTAK,0.553698705712872,602.319,0,BL_ORD_ID:56,463,467,gi|253775657|ref|YP_003038488.1| tryptophanase [Escherichia coli BL21(DE3)],,2,602.318,0.00542841868345953,0
//  236,20130510_EXQ1_IgPa_QC_UPS1_01.1602.1602.2,SISASGHK,0.0530103870541474,785.405,0,BL_ORD_ID:2177,269,276,gi|253773536|ref|YP_003036367.1| glutamate decarboxylase [Escherichia coli BL21(DE3)],,2,785.403,0.00151458248726135,0
  public static Multimap<SpectrumIdentifier, PeptideMatch> readOMSSA(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    Map<String, Modification> modmap = new HashMap<>();
    modmap.put("acetylation of protein n-term", UnimodManager.getModification("Acetyl").get());
    modmap.put("carbamidomethyl C", UnimodManager.getModification("Carbamidomethyl").get());
    modmap.put("deamidation of N and Q", UnimodManager.getModification("Deamidation").get());
    modmap.put("dehydro of S and T", UnimodManager.getModification("Dehydrated").get());
    modmap.put("hydroxylation of D", UnimodManager.getModification("Oxidation").get());
    modmap.put("hydroxylation of F", UnimodManager.getModification("Oxidation").get());
    modmap.put("hydroxylation of K", UnimodManager.getModification("Oxidation").get());
    modmap.put("hydroxylation of N", UnimodManager.getModification("Oxidation").get());
    modmap.put("hydroxylation of P", UnimodManager.getModification("Oxidation").get());
    modmap.put("hydroxylation of  Y", UnimodManager.getModification("Oxidation").get());
    modmap.put("iTRAQ114 on nterm", UnimodManager.getModification("iTRAQ4plex114").get());
    modmap.put("iTRAQ114 on K", UnimodManager.getModification("iTRAQ4plex114").get());
    modmap.put("methyl ester of D", UnimodManager.getModification("Methyl").get());
    modmap.put("methyl ester of E (duplicate of 17)", UnimodManager.getModification("Methyl").get());
    modmap.put("methyl ester of peptide c-term (duplicate of 18)", UnimodManager.getModification("Methyl").get());
    modmap.put("oxidation of M",         UnimodManager.getModification("Hydroxylation").get());
    modmap.put("oxidation of H", UnimodManager.getModification("Oxidation").get());
    modmap.put("oxidation of W", UnimodManager.getModification("Oxidation").get());
    modmap.put("oxidation of W to formylkynurenin", UnimodManager.getModification("Trp->Kynurenin").get());
    modmap.put("oxidation of W to hydroxykynurenin", UnimodManager.getModification("Trp->Hydroxykynurenin").get());
    modmap.put("pyro-glu from n-term Q", UnimodManager.getModification("Pyro-glu").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());
//    modmap.put("", UnimodManager.getModification("").get());

    TabFile file = new TabFile(filename, TabFile.comma);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String[] ids = Strs.split(file.get("Filename/id"), '.');
      String run = Strs.toString(Arrays.copyOfRange(ids, 0, ids.length-3), "."), scan = ids[ids.length-3];
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+ scan);
        id.setName(run + "#" + scan);
        id.setAssumedCharge(file.getInt("Charge"));
        id.setPrecursorNeutralMass(file.getDouble("Mass"));
        id.setPrecursorMz((id.getPrecursorNeutralMass().get() + 1.007825035d * id.getAssumedCharge().get()) /id.getAssumedCharge().get());
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      PeptideMatch match = new PeptideMatch(file.get("Peptide").toUpperCase());
      // add the modifications if necessary: YHSETEmmR; oxidation of M:7 ,oxidation of M:8
      if (Strs.isSet(file.get("Mods")))
      {
        String[] mods = Strs.split(file.get("Mods"), ',', true);
        for (String mod : mods)
        {
          String[] mm = Strs.split(mod, ':');
          int pos = Stats.toInt(mm[1]); // 0-based index
          if (modmap.containsKey(mm[0])) match.addModificationMatch(pos-1, modmap.get(mm[0]));
          else
          {
            System.out.println("Unknown mod: " + mod);
          }
        }
      }
      // gi|253775657|ref|YP_003038488.1| tryptophanase [Escherichia coli BL21(DE3)]
      boolean decoy = (file.get("Defline").indexOf("XXX_")==0);
      String accs = Strs.split(file.get("Defline"), ' ')[0];
      if (Strs.isSet(accs))
      {
        String[]  acs = Strs.split(decoy?accs.substring(4):accs, '|');
        try
        {
          match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
            Optional.of(acs[0]), Optional.of("-"), Optional.of("-"),
            decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
        }
        catch (ArrayIndexOutOfBoundsException e)
        {
          System.out.println(e.toString());
        }
      }

      PSMs.addScore(match, "OMSSA:EVal", file.getDouble("E-value"));
      PSMs.addScore(match, "OMSSA:PVal", file.getDouble("P-value"));
      PSMs.addScore(match, "dbEVal", -10d * Math.log10(file.getDouble("E-value")));
      PSMs.addScore(match, "dbPVal", -10d * Math.log10(file.getDouble("P-value")));
      PSMs.addScore(match, "OMSSA:NIST", file.getDouble("NIST score"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(file.getDouble("Theo Mass"));
      match.setMassDiff(id.getPrecursorNeutralMass().get() - match.getNeutralPeptideMass());

      PSMs.addScore(match, "^Charge", file.getInt("Charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }
    LCMSMS.rank(id_match, "OMSSA:EVal", false, false);

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMSF(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    MsfFile msf = null;
    try
    {
      msf = new MsfFile(new File(filename));
      PeptideLowMemController instance = new PeptideLowMemController();
      List<PeptideLowMem> result = instance.getPeptidesWithConfidenceLevel(1, msf);

      // grab the details
      for (PeptideLowMem peptide : result)
      {
        // TODO no additional info comine back. The query resulted in null 'rs' inside the call!
        List info = instance.getInformationForPeptide(peptide.getPeptideId(), msf, true);
        System.out.println(info.size());
      }

      System.out.println(result.size());
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }

    return null;
  }
//  N       Unused  Total   %Cov    %Cov(50)        %Cov(95)        Accessions      Names   Used    Annotation      Contrib Conf    Sequence        Modifications   ProteinModifications    Cleavages       dMass   Obs MW  Obs m/z Theor MW        Theor m/z       Theor z Sc      Spectrum        Acq Time        Intensity (Peptide)     PrecursorIntensityAcquisition   Apex Time (Peptide)     Elution Peak Width (Peptide)    MS2Counts
//  1       187.69  187.69  88.9800012111664        80.0999999046326        76.5500009059906        gi|253775383    DNA-directed RNA polymerase, beta' subunit [Escherichia coli BL21(DE3)]                 2       99.0000009536743        AAAESSIQVK                              -0.00113623996730894    1002.53344726563        502.274 1002.53454589844        502.274566650391        2       11      1.1.1.6275.1    29.61068        3.156882E+08
// TODO to be completed
  public static Multimap<SpectrumIdentifier, PeptideMatch> readProteinPilot(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    ProteinPilotReader ppReader = new ProteinPilotReader(new ProteinPilotModResolver());
    PsmReaders              psm = new PsmReaders();
    ppReader.parse(new File(filename), psm.insertIdMatch);

    System.out.println(psm.mIdMatch.size());

    return psm.mIdMatch;

//    TabFile file = new TabFile(filename, TabFile.tabb);
//    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
//    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
//    while (file.hasNext())
//    {
//      String run = Strs.stripLastOf(file.get("#SpecFile"), '.'), scan = file.get("ScanNum");
//      SpectrumIdentifier id = run_scan_id.get(run, scan);
//      if (id==null)
//      {
//        id = new SpectrumIdentifier(run+"#"+scan);
//        id.setName(run + "#" + scan);
//        id.setAssumedCharge(file.getInt("Charge"));
//        id.setPrecursorMz(file.getDouble("Precursor"));
//        id.setPrecursorNeutralMass((id.getPrecursorMz().get() - 1.007825035d) * id.getAssumedCharge().get());
//        id.setSpectrumFile(run);
//        id.addScanNumber(new Integer(scan));
//        run_scan_id.put(run, scan, id);
//      }
//      PeptideMatch match = PSMs.fromNumModSequence(file.get("Peptide"));
//      // parse the protein names
//      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
//      String[] accs = Strs.split(file.get("Protein"), ';');
//      if (Tools.isSet(accs))
//        for (String acc : accs)
//        {
//          boolean decoy = (acc.indexOf("XXX_")==0);
//          String[]  acs = Strs.split(decoy?acc.substring(4):acc, '|');
//          int       pre = acc.indexOf("pre=")+4, post = acc.indexOf("post=")+5;
//          try
//          {
//            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
//              Optional.of(acs[0]), Optional.of(acc.substring(pre,pre+1)), Optional.of(acc.substring(post,post+1)),
//              decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
//          }
//          catch (ArrayIndexOutOfBoundsException e)
//          {
//            System.out.println(e.toString());
//          }
//        }
//
////      String[] matched = Strs.split(file.get("Ions Matched"), '/');
//      if (file.getDouble("PrecursorError(ppm)")!=null)
//        match.setMassDiff(file.getDouble("PrecursorError(ppm)") * id.getPrecursorMz().get() * id.getAssumedCharge().get() * 1E-6);
//      else match.setMassDiff(file.getDouble("PrecursorError(Da)"));
//
//      PSMs.addScore(match, "DeNovoScore",  file.getDouble("DeNovoScore"));
//      PSMs.addScore(match, "IsotopeError", file.getInt("IsotopeError"));
//      PSMs.addScore(match, "MSGFScore", file.getDouble("MSGFScore"));
//      PSMs.addScore(match, "SpecEValue",   file.getDouble("SpecEValue"));
//      PSMs.addScore(match, "EValue", file.getDouble("EValue"));
//      PSMs.addScore(match, "QValue", file.getDouble("QValue"));
//      PSMs.addScore(match, "PepQValue", file.getDouble("PepQValue"));
//
//      // add as much information to the match as possible
//      match.setNeutralPeptideMass(id.getPrecursorNeutralMass().get() + match.getMassDiff());
//      PSMs.addScore(match, "^Charge", file.getInt("Charge"));
//
//      if (!id_match.put(id, match))
//      {
////        System.out.println("Duplicated?");
//      }
//    }
//    LCMSMS.rank(id_match, "MSGFScore");

//    return id_match;
  }
//  scan    charge  spectrum precursor m/z  spectrum neutral mass   peptide mass    delta_cn        sp score        sp rank xcorr score     xcorr rank      b/y ions matched        b/y ions total  distinct matches/spectrum       sequence        cleavage type   protein id      flanking aa     original target sequence
//  10254   2       300.175 598.335 598.344 0.0402319       15.4584 1       0.188466        1       2       10      1       PNAVAK  trypsin/p-full-digest   gi|253772588|ref|YP_003035419.1|(12)    KN      PNAVAK
//  10254   2       300.175 598.335 598.344 0       15.4584 2       0.148234        2       2       10      1       PANAVK  trypsin/p-full-digest   decoy_gi|253772588|ref|YP_003035419.1|(12)      KN      PNAVAK
// TODO to be completed
  public static Multimap<SpectrumIdentifier, PeptideMatch> readCrux(String filename) throws IOException
  {
    System.out.println("fetching the PSM from " + filename);

    String run = FilenameUtils.getBaseName(filename);
    if (run.indexOf("Crux_")==0) run = run.substring(5);
    // Strs.stripLastOf(filename, '.');

    TabFile file = new TabFile(filename, TabFile.tabb);
    Table<String, String, SpectrumIdentifier> run_scan_id = HashBasedTable.create();
    Multimap<SpectrumIdentifier, PeptideMatch> id_match   = HashMultimap.create();
    while (file.hasNext())
    {
      String scan = file.get("scan");
      SpectrumIdentifier id = run_scan_id.get(run, scan);
      if (id==null)
      {
        id = new SpectrumIdentifier(run+"#"+ scan);
        id.setName(run+"#"+scan);
        id.setAssumedCharge(file.getInt("charge"));
        id.setPrecursorMz(file.getDouble("spectrum precursor m/z"));
        id.setPrecursorNeutralMass(file.getDouble("spectrum neutral mass"));
        id.setSpectrumFile(run);
        id.addScanNumber(new Integer(scan));
        run_scan_id.put(run, scan, id);
      }
      // IETLMRNLM[15.9949]PWRK
      PeptideMatch match = PSMs.fromNumModSequence(file.get("sequence"));
      // add the fixed Cys mods
      for (int i=0; i<match.size(); i++)
      {
        if ("C".equals(match.getSymbol(i).getSymbol())) match.addModificationMatch(i, 57.02);
      }
      // parse the protein names
      // XXX_gi|4826734|ref|NP_004951.1|(pre=K,post=W);XXX_gi|767988385|ref|XP_011544083.1|(pre=K,post=W);XXX_gi|283135173|ref|NP_001164408.1|(pre=K,post=W);XXX_gi|283135201|ref|NP_001164105.1|(pre=K,post=W);XXX_gi|530407875|ref|XP_005255290.1|(pre=K,post=W);XXX_gi|767988388|ref|XP_011544084.1|(pre=K,post=W)
      String flanking = file.get("flanking aa");
      String[]   accs = Strs.split(file.get("protein id"), ';');
      if (Tools.isSet(accs))
        for (String acc : accs)
        {
          boolean decoy = (acc.indexOf("decoy_")==0);
          String[]  acs = Strs.split(decoy?acc.substring(6):acc, '|');
          try
          {
            match.addProteinMatch(new PeptideProteinMatch(acs.length>1?acs[1]:"unknown",
              Optional.of(acs[0]), Optional.of(Strs.isSet(flanking)?flanking.substring(0,1):"-"),
              Optional.of((flanking!=null&&flanking.length()>1)?flanking.substring(1,2):"-"),
              decoy? PeptideProteinMatch.HitType.DECOY:PeptideProteinMatch.HitType.TARGET));
          }
          catch (ArrayIndexOutOfBoundsException e)
          {
            System.out.println(e.toString());
          }
        }

//      String[] matched = Strs.split(file.get("Ions Matched"), '/');
      Peptide p = match.toPeptide(new NumModMatchResolver());

      PSMs.addScore(match, "delta_cn", file.getDouble("delta_cn"));
      PSMs.addScore(match, "sp score", file.getDouble("sp score"));
      PSMs.addScore(match, "xcorr score",    file.getDouble("xcorr score"));
      PSMs.addScore(match, "b/y ions matched", file.getDouble("b/y ions matched"));
      PSMs.addScore(match, "b/y ions total", file.getDouble("b/y ions total"));
      PSMs.addScore(match, "sp rank", file.getInt("sp rank"));

      // add as much information to the match as possible
      match.setNeutralPeptideMass(p.getMolecularMass());
//      match.setMassDiff(match.getNeutralPeptideMass() - p.getMolecularMass());
      match.setMassDiff(file.getDouble("spectrum neutral mass") - file.getDouble("peptide mass"));

      PSMs.addScore(match, "^Charge", file.getInt("charge"));

      if (!id_match.put(id, match))
      {
//        System.out.println("Duplicated?");
      }
    }
    LCMSMS.rank(id_match, "xcorr score", true, false);

    return id_match;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> readMzID(String filename, int lowest_rank, Engine engine)
  {
    Multimap<SpectrumIdentifier, PeptideMatch> id_match = HashMultimap.create();

    System.out.println("Reading " + filename);

    PsmReaders readers = new PsmReaders(engine, lowest_rank);
    readers.parseMzID(filename);

    PSMs.DesendScorePeptideMatch sorter = new PSMs.DesendScorePeptideMatch(engine.getCanonicalScore());
    List<PeptideMatch>                   matches = new ArrayList<>();
    Iterator<SpectrumIdentifier> idr = readers.getResultMap().keySet().iterator();
    while (idr.hasNext())
    {
      SpectrumIdentifier id = idr.next();

      matches.clear();
      for (PeptideMatch m : readers.getResultMap().get(id))
      {
        // fix the absence of the neutral peptide mass
        try
        {
          // provoke the call which will throw the exception if the value is absent
          if (id.getPrecursorNeutralMass().isPresent()) m.getNeutralPeptideMass();
        }
        catch (IllegalStateException e)
        {
          m.setNeutralPeptideMass(id.getPrecursorNeutralMass().get()-m.getMassDiff());
        }
        matches.add(m);
      }
      id_match.putAll(id, PSMs.trimByRank(matches, sorter, lowest_rank, engine.getCanonicalScore(), engine.getDeltaScore()));
      if (id_match.keySet().size()%1000==0) System.out.print(".");
    }
    System.out.print("\n --> " + readers.getResultMap().keySet().size() + "\n");
    return id_match;
  }
//  // TODO to be completed
//  public static Multimap<SpectrumIdentifier, PeptideMatch> readXTandem(String filename)
//  {
//    Multimap<SpectrumIdentifier, PeptideMatch> id_match = HashMultimap.create();
//
//    System.out.println("Reading " + filename);
//
//    PsmReaders readers = new PsmReaders();
//    readers.parseMzID(filename);
//
//    Iterator<SpectrumIdentifier> idr = readers.getResultMap().keySet().iterator();
//    while (idr.hasNext())
//    {
//      SpectrumIdentifier id = idr.next();
//      id_match.putAll(id, readers.getResultMap().get(id));
//      if (id_match.keySet().size()%1000==0) System.out.print(".");
//    }
//    System.out.print("\n --> " + readers.getResultMap().keySet().size() + "\n");
//    return id_match;
//  }
/*

  // TODO. Need more work in the future
  public static Dataframe fetchMzID(String root, double q, boolean save_decoy)
  {
    PsmReaders readers = new PsmReaders();

    List<String> files = IOs.listFiles(root, new FileFilterExtensions(new String[]{"mzid"}));
    if (Tools.isSet(files))
      for (String f : files)
      {
        System.out.println("Reading " + f);
        readers.parseMzID(f);
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
        row.put("Peptide", PSMs.toNumModSequence(match));
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
*/
}
