package org.ms2ms.io;

import com.google.common.base.Optional;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodManager;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodMod;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.pride.jmztab.model.*;

import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;

/** A more flexible way to format a quantitative study into a mzTab output
  */
public class mzTabBuilder extends MZTabFile
{
  public static CVParam HUMAN = new CVParam("NEWT", "9606", "Homo sapiens (Human)", null);

  public mzTabBuilder(Metadata metadata)
  {
    super(metadata);
  }

  public static Metadata TMT10(Metadata mtd)
  {
    mtd.setQuantificationMethod(new CVParam("MS", "MS:1002010", "TMT quantitation analysis", null));
    mtd.setProteinQuantificationUnit(new CVParam("MS", "MS:10001847", "reporter ion intensity", null));
    mtd.setPeptideQuantificationUnit(new CVParam("MS", "MS:10001847", "reporter ion intensity", null));
    mtd.setSmallMoleculeQuantificationUnit(new CVParam("MS", "MS:10001847", "reporter ion intensity", null));

    // http://www.ebi.ac.uk/ontology-lookup/init.do#soft
    // choose "PRIDE Controlled Vocabulary [PRIDE]
    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "PRIDE:0000285", "TMT reagent", "126"));
    mtd.addAssayQuantificationReagent(2, new CVParam("PRIDE", "PRIDE:0000286", "TMT reagent", "127"));
    mtd.addAssayQuantificationReagent(3, new CVParam("PRIDE", "PRIDE:0000287", "TMT reagent", "128"));
    mtd.addAssayQuantificationReagent(4, new CVParam("PRIDE", "PRIDE:0000288", "TMT reagent", "129"));
    mtd.addAssayQuantificationReagent(5, new CVParam("PRIDE", "PRIDE:0000289", "TMT reagent", "130"));
    mtd.addAssayQuantificationReagent(6, new CVParam("PRIDE", "PRIDE:0000290", "TMT reagent", "131"));
// TODO no CV found for TMT 127N, 127C, etc
//    mtd.addAssayQuantificationReagent(7, new CVParam("PRIDE", "PRIDE:0000291", "TMT reagent", "126"));
//    mtd.addAssayQuantificationReagent(8, new CVParam("PRIDE", "PRIDE:0000292", "TMT reagent", "126"));
//    mtd.addAssayQuantificationReagent(9, new CVParam("PRIDE", "PRIDE:0000293", "TMT reagent", "126"));
//    mtd.addAssayQuantificationReagent(10,new CVParam("PRIDE", "PRIDE:0000294", "TMT reagent", "126"));
//    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "MS:1002038", "unlabeled sample", null));

    return mtd;
  }
  public static Metadata addSamples(Metadata mtd, Sample... samples)
  {
    if (Tools.isSet(samples))
      for (Sample sample : samples)
      {
        int assayid = mtd.getAssayMap().size()+1;
        mtd.addAssaySample(assayid, sample);

        mtd.addAssayQuantificationModParam(assayid, 1, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
        mtd.addAssayQuantificationModParam(assayid, 2, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
        mtd.addAssayQuantificationModSite(assayid, 1, "R");
        mtd.addAssayQuantificationModSite(assayid, 2, "K");
        mtd.addAssayQuantificationModPosition(assayid, 1, "Anywhere");
        mtd.addAssayQuantificationModPosition(assayid, 2, "Anywhere");
      }

    return mtd;
  }
  public static Metadata iTraq2(Metadata mtd, Sample... samples)
  {
    mtd.setQuantificationMethod(new CVParam("MS", "MS:1001837", "iTRAQ quantitation analysis", null));
    mtd.setProteinQuantificationUnit(new CVParam("PRIDE", "PRIDE:0000395", "Ratio", null));
    mtd.setPeptideQuantificationUnit(new CVParam("PRIDE", "PRIDE:0000395", "Ratio", null));
    mtd.setSmallMoleculeQuantificationUnit(new CVParam("PRIDE", "PRIDE:0000395", "Ratio", null));

    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "PRIDE:0000114", "iTRAQ reagent", "114"));
    mtd.addAssayQuantificationReagent(2, new CVParam("PRIDE", "PRIDE:0000115", "iTRAQ reagent", "115"));
    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "MS:1002038", "unlabeled sample", null));

    if (Tools.isSet(samples))
      for (Sample sample : samples)
      {
        int id = mtd.getAssayMap().size()+1;
        mtd.addAssaySample(id, sample);

        mtd.addAssayQuantificationModParam(2, 1, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
        mtd.addAssayQuantificationModParam(2, 2, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
        mtd.addAssayQuantificationModSite(2, 1, "R");
        mtd.addAssayQuantificationModSite(2, 2, "K");
        mtd.addAssayQuantificationModPosition(2, 1, "Anywhere");
        mtd.addAssayQuantificationModPosition(2, 2, "Anywhere");
      }

    return mtd;
  }
  public static Metadata addVarMods(Metadata mtd, String... mods)
  {
    for (String mod : mods)
    {
      Optional<UnimodMod> opt = UnimodManager.getModification(mod);
      if (opt.isPresent())
      {
        int      id = mtd.getVariableModMap().size() + 1;
        UnimodMod m = opt.get();
        mtd.addVariableModSite(id, Strs.toString(m.getSites(), ","));
        mtd.addVariableModParam(id, new CVParam("UNIMOD", "UNIMOD:"+m.getRecordId(), m.getFullName(), null));
        // TODO complete the N/C-termus setting
      }
    }
//    mtd.addFixedModSite(1, "M");
//    mtd.addFixedModParam(2, new CVParam("UNIMOD", "UNIMOD:35", "Oxidation", null));
//    mtd.addFixedModSite(2, "N-term");
//    mtd.addFixedModParam(3, new CVParam("UNIMOD", "UNIMOD:1", "Acetyl", null));
//    mtd.addFixedModPosition(3, "Protein C-term");

//    mtd.addVariableModParam(1, new CVParam("UNIMOD", "UNIMOD:21", "Phospho", null));
//    mtd.addVariableModSite(1, "M");
//    mtd.addVariableModParam(2, new CVParam("UNIMOD", "UNIMOD:35", "Oxidation", null));
//    mtd.addVariableModSite(2, "N-term");
//    mtd.addVariableModParam(3, new CVParam("UNIMOD", "UNIMOD:1", "Acetyl", null));
//    mtd.addVariableModPosition(3, "Protein C-term");

    return mtd;
  }
  public static Metadata addMascot(Metadata mtd)
  {
    mtd.addSoftwareParam(1, new CVParam("MS", "MS:1001207", "Mascot", "2.3"));
    mtd.addProteinSearchEngineScoreParam(1, new CVParam("MS", "MS:1001171", "Mascot:score", null));
    mtd.addPeptideSearchEngineScoreParam(1, new CVParam("MS", "MS:1001153", "search engine specific score", null));

    return mtd;
  }
  public static Metadata addXTandem(Metadata mtd)
  {
    mtd.addPsmSearchEngineScoreParam(1, new CVParam("MS", "MS:1001330", "X!Tandem:expect", null));
    mtd.addPsmSearchEngineScoreParam(2, new CVParam("MS", "MS:1001331", "X!Tandem:hyperscore", null));

    return mtd;
  }
  public static Sample addSample(Metadata mtd, CVParam species, CVParam tissue, CVParam cell, CVParam disease, String desc)
  {
    int id = mtd.getSampleMap().size()+1;
    if (species!=null) mtd.addSampleSpecies(    id, species);
    if (tissue !=null) mtd.addSampleTissue(id, tissue);
    if (cell   !=null) mtd.addSampleCellType(id, cell);
    if (disease!=null) mtd.addSampleDisease(id, disease);
    if (desc   !=null) mtd.addSampleDescription(id, desc);

    return mtd.getSampleMap().get(id);
  }
  public static Metadata newQE(String desc)
  {
    Metadata mtd = new Metadata(new MZTabDescription(MZTabDescription.Mode.Summary, MZTabDescription.Type.Quantification));

    // http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MS
    mtd.addInstrumentName(1, new CVParam("MS", "MS:1001911", "Q Exactive", null));
    mtd.addInstrumentSource(1, new CVParam("MS", "MS:1000073", "ESI", null));
    mtd.addInstrumentSource(2, new CVParam("MS", "MS:1000598", "ETD", null));
//    mtd.addInstrumentAnalyzer(1, new CVParam("MS", "MS:1000291", "linear ion trap", null));
    mtd.addInstrumentAnalyzer(2, new CVParam("MS", "MS:1000484", "orbitrap", null));
    mtd.addInstrumentDetector(1, new CVParam("MS", "MS:1000253", "electron multiplier", null));
    mtd.addInstrumentDetector(2, new CVParam("MS", "MS:1000348", "focal plane collector", null));

    mtd.addSampleProcessingParam(2, new CVParam("SEP", "SEP:00142", "enzyme digestion", null));
    mtd.addSampleProcessingParam(2, new CVParam("MS", "MS:1001251", "Trypsin", null));

    mtd.addMsRunFragmentationMethod(2, new CVParam("MS", "MS:1000422", "HCD", null));

    mtd.addFixedModParam(1, new CVParam("UNIMOD", "UNIMOD:4", "Carbamidomethyl", null));

    mtd.addSoftwareSetting(1, "Fragment tolerance = 0.05Da");
    mtd.addSoftwareSetting(1, "Parent tolerance = 20ppm");

    mtd.addFalseDiscoveryRateParam(new CVParam("MS", "MS:1001364", "pep:global FDR", "0.01"));
/*
    mtd.addFalseDiscoveryRateParam(new CVParam("MS", "MS:1001214", "pep:global FDR", "0.08"));

    mtd.addSoftwareParam(2, new CVParam("MS", "MS:1001561", "Scaffold", "1.0"));

    mtd.addSampleProcessingParam(1, new CVParam("SEP", "SEP:00173", "SDS PAGE", null));

    mtd.addSmallMoleculeSearchEngineScoreParam(1, new CVParam("MS", "MS:1001420", "SpectraST:delta", null));

    mtd.addPublicationItem(1, PublicationItem.Type.PUBMED, "21063943");
    mtd.addPublicationItem(1, PublicationItem.Type.DOI, "10.1007/978-1-60761-987-1_6");
    mtd.addPublicationItem(2, PublicationItem.Type.PUBMED, "20615486");
    mtd.addPublicationItem(2, PublicationItem.Type.DOI, "10.1016/j.jprot.2010.06.008");

    mtd.addContactName(1, "James D. Watson");
    mtd.addContactName(2, "Francis Crick");
    mtd.addContactAffiliation(1, "Cambridge University, UK");
    mtd.addContactAffiliation(2, "Cambridge University, UK");
    mtd.addContactEmail(1, "watson@cam.ac.uk");
    mtd.addContactEmail(2, "crick@cam.ac.uk");
*/
    try
    {
      mtd.addUri(new URI("http://www.ebi.ac.uk/pride/url/to/experiment"));
      mtd.addUri(new URI("http://proteomecentral.proteomexchange.org/cgi/GetDataset"));

      mtd.addMsRunLocation(1, new URL("file://ftp.ebi.ac.uk/path/to/file"));
      mtd.addMsRunLocation(2, new URL("ftp://ftp.ebi.ac.uk/path/to/file"));
    }
    catch (URISyntaxException us)
    {
      throw new RuntimeException(us);
    }
    catch (MalformedURLException mf)
    {
      throw new RuntimeException(mf);
    }

//    mtd.addMsRunFormat(1, new CVParam("MS", "MS:1000584", "mzML file", null));
//    mtd.addMsRunFormat(2, new CVParam("MS", "MS:1001062", "Mascot MGF file", null));
//    mtd.addMsRunIdFormat(1, new CVParam("MS", "MS:1001530", "mzML unique identifier", null));
//    mtd.addMsRunFragmentationMethod(1, new CVParam("MS", "MS:1000133", "CID", null));
//    mtd.addMsRunHash(2, "de9f2c7fd25e1b3afad3e85a0bd17d9b100db4b3");
//    mtd.addMsRunHashMethod(2, new CVParam("MS", "MS:1000569", "SHA-1", null));
//        mtd.addMsRunFragmentationMethod(2, new CVParam("MS", "MS:1000422", "HCD", null));

//    mtd.addCustom(new UserParam("MS operator", "Florian"));

/*
    mtd.addSampleSpecies(1, new CVParam("NEWT", "9606", "Homo sapiens (Human)", null));
    mtd.addSampleSpecies(1, new CVParam("NEWT", "573824", "Human rhinovirus 1", null));
    mtd.addSampleSpecies(2, new CVParam("NEWT", "9606", "Homo sapiens (Human)", null));
    mtd.addSampleSpecies(2, new CVParam("NEWT", "12130", "Human rhinovirus 2", null));
    mtd.addSampleTissue(1, new CVParam("BTO", "BTO:0000759", "liver", null));
    mtd.addSampleCellType(1, new CVParam("CL", "CL:0000182", "hepatocyte", null));
    mtd.addSampleDisease(1, new CVParam("DOID", "DOID:684", "hepatocellular carcinoma", null));
    mtd.addSampleDisease(1, new CVParam("DOID", "DOID:9451", "alcoholic fatty liver", null));
    mtd.addSampleDescription(1, "Hepatocellular carcinoma samples.");
    mtd.addSampleDescription(2, "Healthy control samples.");
    mtd.addSampleCustom(1, new UserParam("Extraction date", "2011-12-21"));
    mtd.addSampleCustom(1, new UserParam("Extraction reason", "liver biopsy"));
*/

    Sample sample1 = mtd.getSampleMap().get(1);
    Sample sample2 = mtd.getSampleMap().get(2);
    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "PRIDE:0000114", "iTRAQ reagent", "114"));
    mtd.addAssayQuantificationReagent(2, new CVParam("PRIDE", "PRIDE:0000115", "iTRAQ reagent", "115"));
    mtd.addAssayQuantificationReagent(1, new CVParam("PRIDE", "MS:1002038", "unlabeled sample", null));
    mtd.addAssaySample(1, sample1);
    mtd.addAssaySample(2, sample2);

    mtd.addAssayQuantificationModParam(2, 1, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
    mtd.addAssayQuantificationModParam(2, 2, new CVParam("UNIMOD", "UNIMOD:188", "Label:13C(6)", null));
    mtd.addAssayQuantificationModSite(2, 1, "R");
    mtd.addAssayQuantificationModSite(2, 2, "K");
    mtd.addAssayQuantificationModPosition(2, 1, "Anywhere");
    mtd.addAssayQuantificationModPosition(2, 2, "Anywhere");

    MsRun msRun1 = mtd.getMsRunMap().get(1);
    mtd.addAssayMsRun(1, msRun1);

    Assay assay1 = mtd.getAssayMap().get(1);
    Assay assay2 = mtd.getAssayMap().get(2);

    mtd.addStudyVariableAssay(1, assay1);
    mtd.addStudyVariableAssay(1, assay2);
    mtd.addStudyVariableSample(1, sample1);
    mtd.addStudyVariableDescription(1, "description Group B (spike-in 0.74 fmol/uL)");

    mtd.addStudyVariableAssay(2, assay1);
    mtd.addStudyVariableAssay(2, assay2);
    mtd.addStudyVariableSample(2, sample1);
    mtd.addStudyVariableDescription(2, "description Group B (spike-in 0.74 fmol/uL)");

    mtd.addCVLabel(1, "MS");
    mtd.addCVFullName(1, "MS");
    mtd.addCVVersion(1, "3.54.0");
    mtd.addCVURL(1, "http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");

    mtd.addProteinColUnit(ProteinColumn.RELIABILITY, new CVParam("MS", "MS:00001231", "PeptideProphet:Score", null));

    MZTabColumnFactory peptideFactory = MZTabColumnFactory.getInstance(Section.Peptide);
//    peptideFactory.addDefaultStableColumns();

    PeptideColumn peptideColumn = (PeptideColumn) peptideFactory.findColumnByHeader("retention_time");
    mtd.addPeptideColUnit(peptideColumn, new CVParam("UO", "UO:0000031", "minute", null));

    mtd.addPSMColUnit(PSMColumn.RETENTION_TIME, new CVParam("UO", "UO:0000031", "minute", null));
    mtd.addSmallMoleculeColUnit(SmallMoleculeColumn.RETENTION_TIME, new CVParam("UO", "UO:0000031", "minute", null));


    return mtd;
  }
}
