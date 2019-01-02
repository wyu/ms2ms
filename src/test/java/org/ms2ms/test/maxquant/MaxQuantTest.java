package org.ms2ms.test.maxquant;

import org.junit.Test;
import org.ms2ms.algo.MsStats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;
import uk.ac.liv.pgb.jmzqml.model.mzqml.*;
import uk.ac.liv.pgb.jmzqml.xml.io.MzQuantMLMarshaller;
import uk.ac.liv.mzqlib.maxquant.converter.MaxquantFilesReader;

import java.io.IOException;
import java.math.BigInteger;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: wyu
 * Date: 7/6/14
 * Time: 6:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class MaxQuantTest extends TestAbstract
{
  String root = "/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/";

  @Test
  public void calibration() throws Exception
  {
    Dataframe evidences = new Dataframe(root+"evidence100.txt", '\t');

    // Dataframe pivot(String col, String val, MsStats.Aggregator func, String... rows)
    Dataframe out = evidences.pivot("Raw file", "Retention time calibration", MsStats.Aggregator.MEAN, "Modified sequence", "Retention time");
    //System.out.println("\n" + out.display());
  }

  @Test
  public void parseSummary() throws IOException
  {
    Dataframe summary = new Dataframe(root+"summary.txt", '\t', "Raw file");
    System.out.println(summary.size());
  }
  @Test
  public void parsePeptides() throws IOException
  {
    Dataframe peptides = new Dataframe(root+"allPeptides.txt", '\t');
    System.out.println(peptides.size());

/*
    To facilitate the inter-sample comparison of protein expression ratios or to assemble data into a specific form—for instance into a time series—
    the software needs to know which LC/MS runs belong to which ‘time point.’ This can be specified by uploading an experimental design file in
    Step 30. If no experimental design file is used, SILAC ratios and other protein information are provided as a whole for the entire data set.
    A template file for the experimental design is created automatically in the ‘combined’ folder. Open this template file in Excel, make appropriate
    changes and save it as a tab-delimited text file. If the ‘Experiment’ column is filled in, separate average ratios (and other information) will be
    reported for each different term used in the Experiment column. For instance, suppose one measures a six-point time course with five SILAC
    experiments by using the time zero sample always as the light state and the samples at five later time points in the heavy states. In addition,
    one may have separated proteins into 10 gel slices, resulting in 50 LC/MS runs. In this case one would fill the ‘Experiment’ column with five
    different terms, which would give rise to an expression ratio matrix with five columns. In addition, in the ‘Slice’ column denote the numbers
    1–10 regarding which gel slice an LC/MS run corresponds to. Additional statistics on slice-specific identifications of proteins will then be
    provided, which might, for instance, be useful to detect isoforms differing in molecular weight. In the ‘Invert’ column, one can specify if
    the ratios originating from certain LC/MS runs should be inverted. This is useful when the labels have been swapped.
*/
    String filename = "/media/data/maxquant/20081129_Orbi6_NaNa_SA_FASP_out/combined/txt/allPeptides.txt";
    // Raw file	Type	Charge	m/z	Mass	Resolution	Number of data points	Number of scans	Number of isotopic peaks
    // PIF	Mass fractional part	Mass deficit	Mass precision [ppm]	Max intensity m/z 0	Retention time	Retention Length
    // Retention Length (FWHM)	Min scan number	Max scan number	Identified	MS/MS IDs Sequence	Length	Modifications
    // Modified sequence	Proteins	Score	Intensity	Intensities MS/MS Count	MSMS Scan Numbers	MSMS Precursors	MSMS Isotope Indices
    // TabFile file = new TabFile(filename, TabFile.tabb);
  }
  public static void convert(String in, String out) throws IOException{

    MaxquantFilesReader maxRd = new MaxquantFilesReader(in);

    /**
     * ****************************************************************
     *
     * start marshalling MzQuantML file
     *
     ******************************************************************
     */
    MzQuantML qml = new MzQuantML();

    String version = "1.0.0";
    qml.setVersion(version);

    Calendar rightnow = Calendar.getInstance();
    qml.setCreationDate(rightnow);

    int day = rightnow.get(Calendar.DATE);
    int month = rightnow.get(Calendar.MONTH) + 1;
    int year = rightnow.get(Calendar.YEAR);

        /*
         * set mzQuantML id
         */
    if (maxRd.isLabelFree()) {
      qml.setId("MaxQuant-Label-Free-" + String.valueOf(day) + String.valueOf(month) + String.valueOf(year));
    } else {
      qml.setId("MaxQuant-SILAC-" + String.valueOf(day) + String.valueOf(month) + String.valueOf(year));
    }

    /**
     * *
     * create CvListType
     */
    CvList cvs = new CvList();
    List<Cv> cvList = cvs.getCv();
    // psi-ms
    Cv cv = new Cv();
    cv.setId("PSI-MS");
    cv.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
    cv.setFullName("Proteomics Standards Initiative Mass Spectrometry Vocabularies");
    cv.setVersion("2.25.0");
    cvList.add(cv);

    //unimod
    Cv cv_unimod = new Cv();
    cv_unimod.setId("UNIMOD");
    cv_unimod.setUri("http://www.unimod.org/obo/unimod.obo");
    cv_unimod.setFullName("Unimod");
    cvList.add(cv_unimod);

    //psi-mod
    Cv cv_mod = new Cv();
    cv_mod.setId("PSI-MOD");
    cv_mod.setUri("http://psidev.cvs.sourceforge.net/psidev/psi/mod/data/PSI-MOD.obo");
    cv_mod.setFullName("Proteomics Standards Initiative Protein Modifications Vocabularies");
    cvList.add(cv_mod);

    qml.setCvList(cvs);

    Label label = new Label();
    CvParam labelCvParam = new CvParam();
    labelCvParam.setAccession("MS:1002038");
    labelCvParam.setName("unlabeled sample");
    labelCvParam.setCv(cv);
    List<ModParam> modParams = label.getModification();
    ModParam modParam = new ModParam();
    modParam.setCvParam(labelCvParam);
    modParams.add(modParam);

        /*
         * create AnalysisSummary
         */

    ParamList analysisSummary = new ParamList();

    if (maxRd.isLabelFree()) {

      analysisSummary.getParamGroup().add(createCvParam("LC-MS label-free quantitation analysis", "PSI-MS", "MS:1001834"));

      CvParam analysisSummaryCv = createCvParam("label-free raw feature quantitation", "PSI-MS", "MS:1002019");
      analysisSummaryCv.setValue("true");
      analysisSummary.getParamGroup().add(analysisSummaryCv);

      analysisSummaryCv = createCvParam("label-free peptide level quantitation", "PSI-MS", "MS:1002020");
      analysisSummaryCv.setValue("true");
      analysisSummary.getParamGroup().add(analysisSummaryCv);

      analysisSummaryCv = createCvParam("label-free protein level quantitation", "PSI-MS", "MS:1002021");
      analysisSummaryCv.setValue("true");
      analysisSummary.getParamGroup().add(analysisSummaryCv);

      analysisSummaryCv = createCvParam("label-free proteingroup level quantitation", "PSI-MS", "MS:1002022");
      analysisSummaryCv.setValue("false");
      analysisSummary.getParamGroup().add(analysisSummaryCv);
    } else {
      analysisSummary.getParamGroup().add(createCvParam("MS1 label-based analysis", "PSI-MS", "MS:1002018"));

      CvParam analysisSummaryCv = createCvParam("MS1 label-based raw feature quantitation", "PSI-MS", "MS:1002001");
      analysisSummaryCv.setValue("true");
      analysisSummary.getParamGroup().add(analysisSummaryCv);

      analysisSummaryCv = createCvParam("MS1 label-based peptide level quantitation", "PSI-MS", "MS:1002002");
      analysisSummaryCv.setValue("true");
      analysisSummary.getParamGroup().add(analysisSummaryCv);

      analysisSummaryCv = createCvParam("MS1 label-based protein level quantitation", "PSI-MS", "MS:1002003");
      analysisSummaryCv.setValue("true");
      analysisSummary.getParamGroup().add(analysisSummaryCv);

      analysisSummaryCv = createCvParam("MS1 label-based proteingroup level quantitation", "PSI-MS", "MS:1002004");
      analysisSummaryCv.setValue("false");
      analysisSummary.getParamGroup().add(analysisSummaryCv);
    }
    // TODO need to fix the error below
    //qml.setAnalysisSummary(analysisSummary);

    /**
     * create AuditCollection
     */
    AuditCollection auditCollection = new AuditCollection();

    Organization uol = new Organization();
    uol.setId("ORG_UOL");
    uol.setName("University of Liverpool");

    Person andy = new Person();
    andy.setFirstName("Andy");
    andy.setLastName("Jones");

    Affiliation aff = new Affiliation();
    aff.setOrganization(uol);
    andy.getAffiliation().add(aff);
    andy.setId("PERS_ARJ");
    auditCollection.getPerson().add(andy);

    Person ddq = new Person();
    ddq.setFirstName("Da");
    ddq.setLastName("Qi");
    ddq.setId("PERS_DQ");
    ddq.getAffiliation().add(aff);
    auditCollection.getPerson().add(ddq);

    // the schema require person before organization
    auditCollection.getOrganization().add(uol);

    qml.setAuditCollection(auditCollection);

    /**
     * *
     * create InputFiles
     */
    InputFiles inputFiles = new InputFiles();
    List<RawFilesGroup> rawFilesGroupList = inputFiles.getRawFilesGroup();

    Iterator iRaw = maxRd.getRawFileList().iterator();
    HashMap<String, String> rawFileNameIdMap = new HashMap<String, String>();
    HashMap<String, ArrayList<String>> rgIdrawIdMap = new HashMap<String, ArrayList<String>>();
    int raw_i = 0;
    while (iRaw.hasNext()) {
      String rawFn = (String) iRaw.next();
      String rawId = "raw_" + Integer.toString(raw_i);
      RawFilesGroup rawFilesGroup = new RawFilesGroup();
      List<RawFile> rawFilesList = rawFilesGroup.getRawFile();
      RawFile rawFile = new RawFile();
      rawFile.setName(rawFn);
      rawFile.setId(rawId);
            /*
             * make up some raw file locations as we don't know the real
             * location from progenesis files
             */
      rawFile.setLocation("../msmsdata/" + rawFn);
      rawFileNameIdMap.put(rawFn, rawId);
      rawFilesList.add(rawFile);

      String rgId = "rg_" + Integer.toString(raw_i);
      rawFilesGroup.setId(rgId);
      rawFilesGroupList.add(rawFilesGroup);

            /*
             * build rgIdrawIdMap raw files group id as key raw file ids
             * (ArrayList) as value
             */
      ArrayList<String> rawIds = rgIdrawIdMap.get(rgId);
      if (rawIds == null) {
        rawIds = new ArrayList<String>();
        rgIdrawIdMap.put(rgId, rawIds);
      }
      rawIds.add(rawId);

      raw_i++;
    }

    //add search databases
    List<SearchDatabase> searchDBs = inputFiles.getSearchDatabase();
    SearchDatabase db = new SearchDatabase();
    db.setId("SD1");
    db.setLocation("sgd_orfs_plus_ups_prots.fasta");
    searchDBs.add(db);
    Param dbName = new Param();
    db.setDatabaseName(dbName);
    UserParam dbNameParam = new UserParam();
    dbNameParam.setName("sgd_orfs_plus_ups_prots.fasta");
    dbName.setParam(dbNameParam);

    qml.setInputFiles(inputFiles);


    /**
     * *
     * create AssayList
     */
    AssayList assays = new AssayList();
    assays.setId("AssayList_1");
    List<Assay> assayList = assays.getAssay();
    Iterator iAss = maxRd.getAssayList().iterator();
    HashMap<String, String> assayNameIdMap = new HashMap<String, String>();
    int ass_i = 0;
    while (iAss.hasNext()) {
      Assay assay = new Assay();
      String assName = (String) iAss.next();
      String assId = "ass_" + Integer.toString(ass_i);
      assay.setId(assId);
      assay.setName(assName);
      assayNameIdMap.put(assName, assId);

      //Object rawFilesGroupRef = assay.getRawFilesGroupRef();
            /*
             * find the corresponding rawFilesGroup from inputFiles
             */
      RawFilesGroup rawFilesGroup = new RawFilesGroup();
      String rawFileName = (String) maxRd.getAssayRawFileMap().get(assName);
      rawFilesGroup.setId("rg_" + rawFileNameIdMap.get(rawFileName).substring(4));

      assay.setRawFilesGroup(rawFilesGroup);

      //label free example
      if (maxRd.isLabelFree()) {
        assay.setLabel(label);
        assayList.add(assay);
      } //non label free example
      else {
        //TODO: this is fixed label modification just used for example file
        //TODO: need to find a better way to form this later

        //find out if it is light or heavy label assay
        if (assName.toLowerCase().contains("light")) {
          assay.setLabel(label);
          assayList.add(assay);
        } else if (assName.toLowerCase().contains("heavy")) {
          Label label_heavy = new Label();
          CvParam label_lysine = new CvParam();
          label_lysine.setAccession("MOD:00582");
          label_lysine.setName("6x(13)C,2x(15)N labeled L-lysine");
          label_lysine.setCv(cv_mod);
          label_lysine.setValue("Lys8");

          CvParam label_arginine = new CvParam();
          label_arginine.setAccession("MOD:00587");
          label_arginine.setName("6x(13)C,4x(15)N labeled L-arginine");
          label_arginine.setCv(cv_mod);
          label_arginine.setValue("Arg10");

          ModParam modParam_lysine = new ModParam();
          modParam_lysine.setCvParam(label_lysine);
          label_heavy.getModification().add(modParam_lysine);

          ModParam modParam_arginine = new ModParam();
          modParam_arginine.setCvParam(label_arginine);
          label_heavy.getModification().add(modParam_arginine);

          assay.setLabel(label_heavy);
          assayList.add(assay);
        }
      }

      ass_i++;
    }

    qml.setAssayList(assays);

    /**
     * *
     * create StudyVariableList
     */
    StudyVariableList studyVariables = new StudyVariableList();
    List<StudyVariable> studyVariableList = studyVariables.getStudyVariable();
    Iterator iStudyGroup = maxRd.getStudyGroupMap().entrySet().iterator();
    while (iStudyGroup.hasNext()) {
      StudyVariable studyVariable = new StudyVariable();
      Map.Entry<String, ArrayList<String>> entry =
          (Map.Entry<String, ArrayList<String>>) iStudyGroup.next();
      String key = entry.getKey();

      studyVariable.setName(key);
      studyVariable.setId("SV_" + key);
      ArrayList<String> value = entry.getValue();
      Iterator iV = value.iterator();
      List assayRefList = studyVariable.getAssayRefs();
      while (iV.hasNext()) {
        Assay assay = new Assay();
        String assayName = (String) iV.next();
        String id = assayNameIdMap.get(assayName);
        assay.setId(id);
        assayRefList.add(assay);
      }
      CvParam cvParam = createCvParam("technical replicate", "PSI-MS", "MS:1001808");
      studyVariable.getCvParam().add(cvParam);
      studyVariableList.add(studyVariable);
    }

    qml.setStudyVariableList(studyVariables);

    /**
     * *
     * create RatioList only for non label free example
     */
    if (!maxRd.isLabelFree()) {
      RatioList ratios = new RatioList();
      List<String> ratioTitles = (List<String>) maxRd.getPeptideRatioMap().get("0");
      for (int i = 0; i < ratioTitles.size(); i++) {
        Ratio ratio = new Ratio();
        String ratioTitle = (String) ratioTitles.get(i);
        ratio.setId(ratioTitle.replace(' ', '_').replace('/', '_'));

        ParamList ratioCalculations = new ParamList();
        ratioCalculations.getCvParam().add(createCvParam("simple ratio of two values", "PSI-MS", "MS:1001848"));
        ratio.setRatioCalculation(ratioCalculations);

        CvParamRef denom_ref = new CvParamRef();
        denom_ref.setCvParam(createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
        ratio.setDenominatorDataType(denom_ref);

        CvParamRef numer_ref = new CvParamRef();
        numer_ref.setCvParam(createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
        ratio.setNumeratorDataType(numer_ref);

        // TODO: a not very smart way to map denominator_ref and numerator_ref
        // TODO: to StudyVariable or Assay.
        // TODO: This is only fixed to H/L ratio
        // TODO: need to re-code
        ArrayList<String> primeStudyVars = (ArrayList<String> )maxRd.getPrimeStudyVariableList();
        for (String primeStudyVar : primeStudyVars) {
          if (ratioTitle.contains(primeStudyVar)) {
            for (StudyVariable sv : qml.getStudyVariableList().getStudyVariable()) {
              if (sv.getId().contains(primeStudyVar + "_H")) {
                ratio.setNumerator(sv);
              } else if (sv.getId().contains(primeStudyVar + "_L")) {
                ratio.setDenominator(sv);
              }
            }
          }
        }
        ratios.getRatio().add(ratio);
      }

      qml.setRatioList(ratios);
    }

    /**
     * *
     * create ProteinList
     */
    HashMap<String, ArrayList<String>> proteinPeptidesMap = (HashMap<String, ArrayList<String>> )maxRd.getProteinPeptidesMap();

    ProteinList proteins = new ProteinList();
    List<Protein> proteinList = proteins.getProtein();

        /*
         * cells the protein list from proteinPeptidesMap.keySet()
         */

    ArrayList protList = new ArrayList(proteinPeptidesMap.keySet());
    Iterator iProt = protList.iterator();
    HashMap<String, String> proteinAccessionIdMap = new HashMap<String, String>();

    while (iProt.hasNext()) {
      //Map.Entry<String, ArrayList<String>> entry = (Map.Entry<String, ArrayList<String>>) iPro.next();
      String protAccession = (String) iProt.next();

      ArrayList<String> pepSequences = proteinPeptidesMap.get(protAccession);
      int id = protList.indexOf(protAccession);
      String protId = "prot_" + Integer.toString(id);

      Protein protein = new Protein();
      protein.setId(protId);
      protein.setAccession(protAccession);
      proteinAccessionIdMap.put(protAccession, protId);
      protein.setSearchDatabase(db);

      if (pepSequences != null) {
        Iterator iPep = pepSequences.iterator();
        List peptideConsensusRefList = protein.getPeptideConsensusRefs();
        ArrayList<String> pepIds = new ArrayList<String>();
        while (iPep.hasNext()) {
          PeptideConsensus peptideConsensus = new PeptideConsensus();
          String pepSeq = (String) iPep.next();
          if (!pepSeq.isEmpty()) {
            String pepId = "pep_" + pepSeq;
            peptideConsensus.setId(pepId);
            peptideConsensus.setPeptideSequence(pepSeq);

            // avoid duplicate peptide existing in one protein
            if (!pepIds.contains(pepId)) {
              peptideConsensusRefList.add(peptideConsensus);
              pepIds.add(pepId);
            }
          }
        }
      }
      proteinList.add(protein);

    }

        /*
         * add GlobleQuantLayer to ProteinList row type is protein ref column
         * types are .......
         */

    // AssayQuantLayer for protein intensity if it is label free example
    if (maxRd.isLabelFree()) {
      QuantLayer assayQL_prot_int = new QuantLayer();
      assayQL_prot_int.setId("Prot_Assay_QL1");
      CvParamRef cvParamRef_prot_int = new CvParamRef();
      cvParamRef_prot_int.setCvParam(createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
      assayQL_prot_int.setDataType(cvParamRef_prot_int);

      iAss = maxRd.getAssayList().iterator();
      while (iAss.hasNext()) {
        String assName = (String) iAss.next();
        String assId = assayNameIdMap.get(assName);
        Assay assay = new Assay();
        assay.setId(assId);
        assayQL_prot_int.getColumnIndex().add(assay);
      }

      DataMatrix protIntDM = new DataMatrix();
      Iterator iProtInt = (Iterator )maxRd.getProteinIntensityMap().iterator();
      while (iProtInt.hasNext()) {
        Map.Entry<String, ArrayList<String>> entry = (Map.Entry<String, ArrayList<String>>) iProtInt.next();
        String key = entry.getKey();
        String proteinId = proteinAccessionIdMap.get(key);
        if (proteinId != null) {
          Protein protein = new Protein();
          protein.setId(proteinId);
          Row row = new Row();
          row.setObject(protein);

          ArrayList<String> value = entry.getValue();
          Iterator iV = value.iterator();
          while (iV.hasNext()) {
            String protIntV = (String) iV.next();
            row.getValue().add(protIntV);
          }
          protIntDM.getRow().add(row);
        }
      }
      assayQL_prot_int.setDataMatrix(protIntDM);
      proteins.getAssayQuantLayer().add(assayQL_prot_int);
    }

    // StudyVariableQuantLayer for protein intensity if it is NOT label free example
    if (!maxRd.isLabelFree()) {
      int labelNum = maxRd.getLabelNumber();
      QuantLayer assayQL_prot_int = new QuantLayer();
      assayQL_prot_int.setId("Prot_StudyVariable_QL1");
      CvParamRef cvParamRef_prot_int = new CvParamRef();
      cvParamRef_prot_int.setCvParam(createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
      assayQL_prot_int.setDataType(cvParamRef_prot_int);

      Iterator iPrimeStu = maxRd.getPrimeStudyVariableList().iterator();
      while (iPrimeStu.hasNext()) {
        String studyVariableName = (String) iPrimeStu.next();

        StudyVariable studyVariable_L = new StudyVariable();
        studyVariable_L.setId("SV_" + studyVariableName + "_L");
        assayQL_prot_int.getColumnIndex().add(studyVariable_L);

        StudyVariable studyVariable_H = new StudyVariable();
        studyVariable_H.setId("SV_" + studyVariableName + "_H");
        assayQL_prot_int.getColumnIndex().add(studyVariable_H);
      }

      DataMatrix protIntDM = new DataMatrix();
      Iterator iProtInt = (Iterator )maxRd.getProteinIntensityMap().iterator();
      while (iProtInt.hasNext()) {
        Map.Entry<String, ArrayList<String>> entry = (Map.Entry<String, ArrayList<String>>) iProtInt.next();
        String key = entry.getKey();
        String proteinId = proteinAccessionIdMap.get(key);
        if (proteinId != null) {
          Protein protein = new Protein();
          protein.setId(proteinId);
          Row row = new Row();
          row.setObject(protein);

          ArrayList<String> value = entry.getValue();
          Iterator iV = value.iterator();
          while (iV.hasNext()) {
            String protIntV = (String) iV.next();
            row.getValue().add(protIntV);
          }
          protIntDM.getRow().add(row);
        }
      }
      assayQL_prot_int.setDataMatrix(protIntDM);
      proteins.getStudyVariableQuantLayer().add(assayQL_prot_int);
    }

    //AssayQuantLayer for unique peptides if it is label-free example
    if (maxRd.isLabelFree()) {
      QuantLayer assayQL_prot_uniqpep = new QuantLayer();
      assayQL_prot_uniqpep.setId("Prot_Assay_QL2");
      CvParamRef cvParamRef_prot_uniqpep = new CvParamRef();
      cvParamRef_prot_uniqpep.setCvParam(createCvParam("MaxQuant:peptide counts (unique)", "PSI-MS", "MS:1001897"));
      assayQL_prot_uniqpep.setDataType(cvParamRef_prot_uniqpep);

      iAss = maxRd.getAssayList().iterator();
      while (iAss.hasNext()) {
        String assName = (String) iAss.next();
        String assId = assayNameIdMap.get(assName);
        Assay assay = new Assay();
        assay.setId(assId);
        assayQL_prot_uniqpep.getColumnIndex().add(assay);
      }

      DataMatrix protUniqPepDM = new DataMatrix();
      Iterator iProtUniqPep = (Iterator )maxRd.getProteinUniquePeptiedsMap().iterator();
      while (iProtUniqPep.hasNext()) {
        Map.Entry<String, ArrayList<String>> entry = (Map.Entry<String, ArrayList<String>>) iProtUniqPep.next();
        String key = entry.getKey();
        String proteinId = proteinAccessionIdMap.get(key);
        if (proteinId != null) {
          Protein protein = new Protein();
          protein.setId(proteinId);
          Row row = new Row();
          row.setObject(protein);

          ArrayList<String> value = entry.getValue();
          Iterator iV = value.iterator();
          while (iV.hasNext()) {
            String protUniqPepV = (String) iV.next();
            row.getValue().add(protUniqPepV);
          }
          protUniqPepDM.getRow().add(row);
        }
      }
      assayQL_prot_uniqpep.setDataMatrix(protUniqPepDM);
      proteins.getAssayQuantLayer().add(assayQL_prot_uniqpep);
    }

    proteins.setId("ProtL1");
    qml.setProteinList(proteins);


    /**
     * *
     * create FeatureList
     */
    List<FeatureList> featureLists = qml.getFeatureList();
    HashMap<String, FeatureList> rgIdFeatureListMap = new HashMap<String, FeatureList>();
    HashMap<String, ArrayList<Feature>> peptideFeaturesMap = new HashMap<String, ArrayList<Feature>>();
    HashMap<String, ArrayList<String>> peptideAssaysMap = new HashMap<String, ArrayList<String>>();
    HashMap<String, String> featureAssNameMap = new HashMap<String, String>();

    HashMap<String, ArrayList<String>> evidenceMap = (HashMap<String, ArrayList<String>> )maxRd.getEvidenceMap();
    Iterator iEvd = evidenceMap.entrySet().iterator();
    while (iEvd.hasNext()) {
      Map.Entry<String, ArrayList<String>> entry = (Map.Entry<String, ArrayList<String>>) iEvd.next();
      String key = entry.getKey();
      ArrayList<String> value = entry.getValue();

      //TODO: revisit key_H
      // @key_H is the id number for heavy label created artificial
      String key_H = String.valueOf(Integer.parseInt(key) + evidenceMap.size());

      if (maxRd.isLabelFree()) {

        Feature feature = new Feature();
                /*
                 * The positions for features in label free example are:
                 * (0)Sequence,
                 * (1)Modifications,
                 * (2)Leading Proteins,
                 * (3)Raw File,
                 * (4)Experiment,
                 * (5)Charge,
                 * (6)m/z,
                 * (7)Retention Time,
                 * (8)Intensity.
                 */
        String reten = value.get(7);
        String mz = value.get(6);
        String chr = value.get(5);
        String rfn = value.get(3);
        String intensity = value.get(8);
        if (intensity.isEmpty()) {
          intensity = "0";
        }

        feature.setMz(Double.parseDouble(mz));
        feature.setRt(reten);
        feature.setCharge(chr);

        String ftId = "ft_" + key;
        feature.setId(ftId);

        // create peptide sequence to feture HashMap: peptideFeaturesMap
        String pepSeq = value.get(0);
        if (pepSeq != null) {
          ArrayList<Feature> fList = peptideFeaturesMap.get(pepSeq);
          if (fList == null) {
            fList = new ArrayList<Feature>();
            peptideFeaturesMap.put(pepSeq, fList);
          }
          fList.add(feature);
        }

        String rgId = "rg_" + rawFileNameIdMap.get(rfn + ".raw").substring(4);

        // create peptide sequence to assay id HashMap: peptideAssaysMap
        String assayId = assayNameIdMap.get(rfn);
        if (pepSeq != null) {
          ArrayList<String> aList = peptideAssaysMap.get(pepSeq);
          if (aList == null) {
            aList = new ArrayList<String>();
            peptideAssaysMap.put(pepSeq, aList);
          }
          if (!aList.contains(assayId)) {
            aList.add(assayId);
          }
        }

        FeatureList features = rgIdFeatureListMap.get(rgId);
        if (features == null) {
          features = new FeatureList();
          RawFilesGroup rawFilesGroup = new RawFilesGroup();
          rawFilesGroup.setId(rgId);
          features.setRawFilesGroup(rawFilesGroup);
          String fListId = "Flist_" + rgId.substring(3);
          features.setId(fListId);
          rgIdFeatureListMap.put(rgId, features);
          featureLists.add(features);

          // create feature QuantLayer for label free
          List<GlobalQuantLayer> featureQuantLayerList = features.getFeatureQuantLayer();
          GlobalQuantLayer featureQuantLayer = new GlobalQuantLayer();
          featureQuantLayerList.add(featureQuantLayer);
          featureQuantLayer.setId("FQL_" + fListId.substring(6));
          ColumnDefinition featureColumnIndex = new ColumnDefinition();
          featureQuantLayer.setColumnDefinition(featureColumnIndex);

          // cvParam for intensity
          Column featureColumn_int =
              createColumn(0, createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
          featureColumnIndex.getColumn().add(featureColumn_int);

          Row row = new Row();
          row.setObject(feature);
          row.getValue().add(intensity);

          DataMatrix dataMatrix = new DataMatrix();
          dataMatrix.getRow().add(row);
          featureQuantLayer.setDataMatrix(dataMatrix);
        } else {
          Row row = new Row();
          row.setObject(feature);
          row.getValue().add(intensity);
          features.getFeatureQuantLayer().get(0).getDataMatrix().getRow().add(row);
        }
        features.getFeature().add(feature);
        featureAssNameMap.put(feature.getId(), rfn);
      } /**
       * *
       * not label free
       */
      else {
        Feature feature_L = new Feature();
        Feature feature_H = new Feature();
                /*
                 * The positions for features in non label free example are:
                 * (0)Sequence,
                 * (1)Modifications,
                 * (2)Leading Proteins,
                 * (3)Raw File,
                 * (4)Experiment,
                 * (5)Charge,
                 * (6)m/z,
                 * (7)Retention Time,
                 * (8)Intensity,
                 * (9)Intensity L,
                 * (10)Intensity H.
                 */
        String reten = value.get(7);
        String mz = value.get(6);
        String chr = value.get(5);
        String rfn = value.get(3);
        String intensityL = value.get(9);
        if (intensityL.isEmpty()) {
          intensityL = "0";
        }
        String intensityH = value.get(10);
        if (intensityH.isEmpty()) {
          intensityH = "0";
        }

        // light label feature
        feature_L.setMz(Double.parseDouble(mz));
        feature_L.setRt(reten);
        feature_L.setCharge(chr);

        String ftId_L = "ft_" + key;
        feature_L.setId(ftId_L);

        // heavey label feature

        //TODO: this mass shift is artificial
        double massShift = 18;
        double mz_H = Double.parseDouble(mz) + massShift / Double.parseDouble(chr);
        feature_H.setMz(mz_H);

        //TODO: retention time is als artifical
        feature_H.setRt(reten);

        feature_H.setCharge(chr);

        String ftId_H = "ft_" + key_H;
        feature_H.setId(ftId_H);

        // create peptide sequence to feture HashMap: peptideFeaturesMap
        String pepSeq = value.get(0);
        if (pepSeq != null) {
          ArrayList<Feature> fList = peptideFeaturesMap.get(pepSeq);
          if (fList == null) {
            fList = new ArrayList<Feature>();
            peptideFeaturesMap.put(pepSeq, fList);
          }
          fList.add(feature_L);
          fList.add(feature_H);
        }

        String rgId = "rg_" + rawFileNameIdMap.get(rfn + ".raw").substring(4);

        String assayId = assayNameIdMap.get(rfn + "_Light");
        if (pepSeq != null) {
          ArrayList<String> aList = peptideAssaysMap.get(pepSeq);
          if (aList == null) {
            aList = new ArrayList<String>();
            peptideAssaysMap.put(pepSeq, aList);
          }

          if (!aList.contains(assayId)) {
            aList.add(assayId);
          }
        }

        assayId = assayNameIdMap.get(rfn + "_Heavy");
        if (pepSeq != null) {
          ArrayList<String> aList = peptideAssaysMap.get(pepSeq);
          if (aList == null) {
            aList = new ArrayList<String>();
            peptideAssaysMap.put(pepSeq, aList);
          }

          if (!aList.contains(assayId)) {
            aList.add(assayId);
          }
        }

        FeatureList features = rgIdFeatureListMap.get(rgId);
        if (features == null) {
          features = new FeatureList();
          RawFilesGroup rawFilesGroup = new RawFilesGroup();
          rawFilesGroup.setId(rgId);
          features.setRawFilesGroup(rawFilesGroup);
          String fListId = "Flist_" + rgId.substring(3);
          features.setId(fListId);
          rgIdFeatureListMap.put(rgId, features);
          featureLists.add(features);

          // create feature QuantLayer for non label free
          List<GlobalQuantLayer> featureQuantLayerList = features.getFeatureQuantLayer();
          GlobalQuantLayer featureQuantLayer = new GlobalQuantLayer();
          featureQuantLayerList.add(featureQuantLayer);
          featureQuantLayer.setId("FQL_" + fListId.substring(6));
          ColumnDefinition featureColumnIndex = new ColumnDefinition();
          featureQuantLayer.setColumnDefinition(featureColumnIndex);

          // cvParam for intensity
          Column featureColumn_int =
              createColumn(0, createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
          featureColumnIndex.getColumn().add(featureColumn_int);

          Row row_L = new Row();
          row_L.setObject(feature_L);
          row_L.getValue().add(intensityL);

          //TODO: this is artificial row for heavy label
          Row row_H = new Row();
          row_H.setObject(feature_H);
          row_H.getValue().add(intensityH);

          DataMatrix dataMatrix = new DataMatrix();
          dataMatrix.getRow().add(row_L);
          dataMatrix.getRow().add(row_H);
          featureQuantLayer.setDataMatrix(dataMatrix);
        } else {
          Row row_L = new Row();
          row_L.setObject(feature_L);
          row_L.getValue().add(intensityL);
          features.getFeatureQuantLayer().get(0).getDataMatrix().getRow().add(row_L);

          //TODO: this is artificial row for heavy label
          Row row_H = new Row();
          row_H.setObject(feature_H);
          row_H.getValue().add(intensityH);
          features.getFeatureQuantLayer().get(0).getDataMatrix().getRow().add(row_H);
        }
        features.getFeature().add(feature_L);
        features.getFeature().add(feature_H);
        featureAssNameMap.put(feature_L.getId(), rfn + "_Light");
        featureAssNameMap.put(feature_H.getId(), rfn + "_Heavy");
      } // end for non label free
    }

    /**
     * *
     * create PeptideConsensusList
     */
    List<PeptideConsensusList> peptideConsensusListList = qml.getPeptideConsensusList();
    PeptideConsensusList peptideConsensuses = new PeptideConsensusList();

    List<PeptideConsensus> peptideList = peptideConsensuses.getPeptideConsensus();

    Iterator iPep = maxRd.getPeptideList().iterator();
    // peptide feature id map
    HashMap<String, ArrayList<String>> peptideFeatureIdsMap = (HashMap<String, ArrayList<String>> )maxRd.getPeptideEvidenceIdsMap();

    DataMatrix pep_IntDM = new DataMatrix();
    DataMatrix pep_RatioDM = new DataMatrix();

    while (iPep.hasNext()) {
      PeptideConsensus peptideConsensus = new PeptideConsensus();
      String key = (String) iPep.next();

      if (!key.isEmpty()) {
        ArrayList chrKeys = (ArrayList) peptideFeatureIdsMap.get(key);
        Iterator iChr = chrKeys.iterator();
        while (iChr.hasNext()) {
          String chrKey = (String) iChr.next();
          String charge = evidenceMap.get(chrKey).get(5);
          if (!peptideConsensus.getCharge().contains(charge)) {
            peptideConsensus.getCharge().add(charge);
          }

          peptideConsensus.setId("pep_" + key);
          peptideConsensus.setPeptideSequence(key);
          peptideConsensus.setSearchDatabase(db);
        }
        // add Feature_refs to individual peptideConsensus
        ArrayList<Feature> fList = peptideFeaturesMap.get(key);
        Iterator iFList = fList.iterator();
        while (iFList.hasNext()) {
          Feature f = (Feature) iFList.next();

          EvidenceRef evRef = new EvidenceRef();
          evRef.setFeature(f);
          peptideConsensus.getEvidenceRef().add(evRef);

          //add assay_refs
          String assN = featureAssNameMap.get(f.getId());
          String ass_id = assayNameIdMap.get(assN);

          Assay tempAssay = new Assay();
          tempAssay.setId(ass_id);
          tempAssay.setName(assN);
          evRef.getAssays().add(tempAssay);
        }

        // create DataMatrix of AssayQuantLayer for intensity
        ArrayList intKeys = (ArrayList) maxRd.getPeptideIntensityMap().get(key);
        Iterator iInt = intKeys.iterator();
        Row row = new Row();
        row.setObject(peptideConsensus);
        while (iInt.hasNext()) {
          String intKey = (String) iInt.next();
          row.getValue().add(intKey);

          //TODO: this is artificial (duplicate) intensity value for heavy label assay
          if (!maxRd.isLabelFree()) {
            row.getValue().add(intKey);
          }
        }
        pep_IntDM.getRow().add(row);

        // create DataMatrix of RatioQuantLayer for non label free example
        if (!maxRd.isLabelFree()) {
          HashMap<String, List<String>> peptideToRatioMap = (HashMap<String, List<String>> )maxRd.getPeptideRatioMap();
          List<String> ratioValues = peptideToRatioMap.get(key);
          Row row_ratio = new Row();
          row_ratio.setObject(peptideConsensus);
          row_ratio.getValue().addAll(ratioValues);
          pep_RatioDM.getRow().add(row_ratio);
        }

        peptideList.add(peptideConsensus);
      }
    }

        /*
         * add AssayQuantLayer to peptides
         */
    QuantLayer pepAQL_int = new QuantLayer();
    pepAQL_int.setId("Pep_AQL_1");
    CvParamRef cvParamRef_pep_aql_int = new CvParamRef();
    cvParamRef_pep_aql_int.setCvParam(createCvParam("MaxQuant:feature intensity", "PSI-MS", "MS:1001903"));
    pepAQL_int.setDataType(cvParamRef_pep_aql_int);

    iAss = maxRd.getAssayList().iterator();
    while (iAss.hasNext()) {
      String assName = (String) iAss.next();
      String assId = assayNameIdMap.get(assName);
      Assay assay = new Assay();
      assay.setId(assId);
      pepAQL_int.getColumnIndex().add(assay);
    }
    pepAQL_int.setDataMatrix(pep_IntDM);
    peptideConsensuses.getAssayQuantLayer().add(pepAQL_int);

        /*
         * add RatioQuantLayer to peptides if not lable free
         */
    if (!maxRd.isLabelFree()) {
      QuantLayer pepRQL = new QuantLayer();
      pepRQL.setId("Pep_RQL_1");
      pepRQL.setDataType(cvParamRef_pep_aql_int);

      pepRQL.setDataMatrix(pep_RatioDM);
      for (Ratio ratio : qml.getRatioList().getRatio()) {
        pepRQL.getColumnIndex().add(ratio);
      }
      // TODO need to fix this following error
      //peptideConsensuses.setRatioQuantLayer(pepRQL);
    }

    peptideConsensuses.setId("PepList1");
    peptideConsensuses.setFinalResult(true);

    peptideConsensusListList.add(peptideConsensuses);


    /**
     * *
     * create SoftwareList
     */
    SoftwareList softwareList = new SoftwareList();
    Software software = new Software();
    softwareList.getSoftware().add(software);
    software.setId("MaxQuant");
    software.setVersion("1.2.0.18");
    software.getCvParam().add(createCvParam("MaxQuant", "PSI-MS", "MS:1001583"));
    qml.setSoftwareList(softwareList);

    /**
     * *
     * create DataProcessingList
     */
    DataProcessingList dataProcessingList = new DataProcessingList();
    DataProcessing dataProcessing = new DataProcessing();
    dataProcessing.setId("feature_quantification");
    dataProcessing.setSoftware(software);
    dataProcessing.setOrder(BigInteger.ONE);
    ProcessingMethod processingMethod = new ProcessingMethod();
    processingMethod.setOrder(BigInteger.ONE);
    processingMethod.getParamGroup().add(createCvParam("quantification data processing", "PSI-MS", "MS:1001861"));
    dataProcessing.getProcessingMethod().add(processingMethod);

    dataProcessingList.getDataProcessing().add(dataProcessing);
    qml.setDataProcessingList(dataProcessingList);

    /**
     * *
     * create a Marshaller and marshal to File
     */
    //String mzqFn = in + "//maxquant.mzq";

    MzQuantMLMarshaller marshaller = new MzQuantMLMarshaller(out);
    marshaller.marshall(qml);

  }
// create a CvParamType instance

  private static CvParam createCvParam(String name, String cvRef,
                                       String accession) {
    CvParam cp = new CvParam();
    cp.setName(name);
    Cv cv = new Cv();
    cv.setId(cvRef);
    cp.setCv(cv);
    cp.setAccession(accession);
    return cp;
  }

  private static Column createColumn(long index, CvParam cvParam) {
    Column column = new Column();
    column.setIndex(BigInteger.valueOf(index));
    CvParamRef cvParamRef = new CvParamRef();
    cvParamRef.setCvParam(cvParam);
    column.setDataType(cvParamRef);
    return column;
  }
}
