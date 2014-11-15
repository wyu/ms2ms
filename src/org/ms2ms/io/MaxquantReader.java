package org.ms2ms.io;

/**Adapted from MaxQuantFilesReader by Da Qi for maxquantmzquantmlconvertor project
 *
 * Change: set isLebelFree=true as the default
 * Change: allow the absence of experimental design file.
 *         quoting from the protocol from NBT: If no experimental design file is used, SILAC ratios and other protein
 *         information are provided as a whole for the entire data set.
 */

import au.com.bytecode.opencsv.CSVReader;
import org.ms2ms.utils.Tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author Da Qi
 */
public class MaxquantReader
{
  private static final String EVIDENCE = "evidence.txt";
  private static final String PEPTIDES = "peptides.txt";
  private static final String PROTEINS = "proteinGroups.txt";
  private static final String DESIGN = "experimentalDesignTemplate.txt";
  private static final String PARAM = "parameters.txt";
  private BufferedReader br;
  private FileReader rd;
  private CSVReader csvReader;
  private ArrayList assays;
  private ArrayList experiments;
  private ArrayList studyVars;
  private ArrayList primeStudyVars;
  private ArrayList peptides;
  private ArrayList rawFiles;
  private HashMap<String, String> assayExpMap;
  private HashMap<String, String> assayRawFileMap;
  private HashMap<String, ArrayList<String>> evidenceMap;
  private HashMap<String, ArrayList<String>> primeStudyGroupMap;
  private HashMap<String, ArrayList<String>> studyGroupMap;
  private HashMap<String, ArrayList<String>> peptideToIntMap;
  private HashMap<String, ArrayList<String>> peptideToEvdMap;
  private HashMap<String, List<String>> peptideToRatioMap;
  private HashMap<String, ArrayList<String>> peptideToProtMap;
  private HashMap<String, ArrayList<String>> proteinToPepMap;
  private HashMap<String, ArrayList<String>> proteinToIntMap;
  private HashMap<String, ArrayList<String>> proteinToUniqPepMap;
  private Boolean hasProteinGroupsFile = false;
  private Boolean isLabelFree = true; // CHANGE: was false
  private int multiplicity;


  public void readPeptides() throws IOException
  {

  }
  public MaxquantReader(String fn) throws IOException
  {
    File path = new File(fn);
    String[] fileList = path.list();
    String[] requireList = {EVIDENCE, PEPTIDES, PROTEINS, PARAM};

    String parentPath = path.getPath();

    /**
     * check if the file folder contains all the required files
     */
    if (!containsAll(fileList, requireList)) {
      System.out.println("The selected folder doesn't contain one of the following files requried:");
      System.out.println("(" + EVIDENCE + ", " + PEPTIDES + ", " + PROTEINS + ", " + PARAM + ").");
    } else {

      /**
       * check if proteinGroups.txt exists in the path
       */
      List<String> fList = Arrays.asList(fileList);
      if (fList.contains(PROTEINS)) {
        hasProteinGroupsFile = true;
      }

      String nextLine[];

      /**
       * *
       * read from parameters.txt and find out the type of the experiment
       * (e.g. label free or SILAC) There are two parameters indicating
       * the difference. One is 'Label-free protein quantification' with
       * value TRUE or FALSE associated. The other is 'Multiplicity' with
       * positive integer associated. One means label-free, two means
       * SILAC with two different labelling.
       */
      multiplicity = 0;

      rd = new FileReader(parentPath + "//" + PARAM);
      br = new BufferedReader(rd);
      csvReader = new CSVReader(br, '\t');

      while ((nextLine = csvReader.readNext()) != null) {
        // TODO need a different mechanism to detect the labeling
        if (nextLine[0].toLowerCase().contains("label-free protein")) {
          isLabelFree = Boolean.valueOf(nextLine[1]);
        } else if (nextLine[0].toLowerCase().contains("multiplicity")) {
          multiplicity = Integer.valueOf(nextLine[1]);
        }
      }

            /*
             * Adjust multiplicity for label-free
             */
      if (isLabelFree) {
        multiplicity = 0;
      }

      csvReader.close();

      /**
       * *
       * read from experimentalDesignTemplate.txt
       */
      if (Tools.contains(fileList, "", 0))
      rd = new FileReader(parentPath + "//" + DESIGN);
      br = new BufferedReader(rd);
      csvReader = new CSVReader(br, '\t');

      /**
       * build the assay list and raw file list
       */
      csvReader.readNext();
      assays = new ArrayList();
      experiments = new ArrayList();
      rawFiles = new ArrayList();
      assayExpMap = new HashMap<String, String>();
      assayRawFileMap = new HashMap<String, String>();
      String assayName;
      String rawFileName;
      int expId = 0;
      while ((nextLine = csvReader.readNext()) != null) {
        // build different assay lists for label-free or SILAC method
        if (isLabelFree) {
          assayName = nextLine[0];
          rawFileName = nextLine[0] + ".raw";
          assays.add(assayName);
          assayRawFileMap.put(assayName, rawFileName);
        } else if (multiplicity == 2) {
          rawFileName = nextLine[0] + ".raw";
          assayName = nextLine[0] + "_Light";
          assays.add(assayName);
          assayRawFileMap.put(assayName, rawFileName);
          assayName = nextLine[0] + "_Heavy";
          assays.add(assayName);
          assayRawFileMap.put(assayName, rawFileName);
        }
        // end of build different assay lists for .....
        if (nextLine[2].isEmpty()) {
          experiments.add(expId++);
        } else {
          experiments.add(nextLine[2]);
        }
        rawFiles.add(nextLine[0] + ".raw");
        assayExpMap.put(nextLine[0], nextLine[2]);
      }
      csvReader.close();

      /**
       * build a study variable group list and map if available
       */
      primeStudyGroupMap = new HashMap<String, ArrayList<String>>();
      ArrayList tempAssays;
      primeStudyVars = new ArrayList<String>();
      for (int i = 0; i < experiments.size(); i++) {
        String studyV;
        int prefixLen = 0;
        tempAssays = new ArrayList();
        studyV = new String();
        for (int j = 0; j < experiments.size(); j++) {
          if (!(i == j)) {
            int tempLen = lengthOfPrefix((String) experiments.get(i), (String) experiments.get(j));
            if ((tempLen > 0) && (prefixLen < tempLen)) {
              prefixLen = tempLen;
              studyV = ((String) experiments.get(i)).substring(0, prefixLen);
              tempAssays.clear();
              if (isLabelFree) {
                tempAssays.add(assays.get(j));
              } else if (multiplicity == 2) {
                tempAssays.add(assays.get(j * multiplicity));
                tempAssays.add(assays.get(j * multiplicity + 1));
              }
            } else if ((prefixLen > 0) && (prefixLen == tempLen)) {
              if (isLabelFree) {
                tempAssays.add(assays.get(j));
              } else if (multiplicity == 2) {
                tempAssays.add(assays.get(j * multiplicity));
                tempAssays.add(assays.get(j * multiplicity + 1));
              }
            }
          }
        }
        if (prefixLen == 0) {
          studyV = (String) experiments.get(i);
        }
        if (isLabelFree) {
          tempAssays.add(assays.get(i));
        } else if (multiplicity == 2) {
          tempAssays.add(assays.get(i * multiplicity));
          tempAssays.add(assays.get(i * multiplicity + 1));
        }
        if (primeStudyGroupMap.get(studyV) == null) {
          primeStudyGroupMap.put(studyV, new ArrayList(tempAssays));
        }
        if (!primeStudyVars.contains(studyV)) {
          primeStudyVars.add(studyV);
        }
      }

      studyGroupMap = (HashMap) primeStudyGroupMap.clone();
      HashMap cloneStudyGroupMap = new HashMap();
      cloneStudyGroupMap = (HashMap) primeStudyGroupMap.clone();
      if (!isLabelFree && multiplicity == 2) {
        Iterator iStudy = cloneStudyGroupMap.keySet().iterator();
        while (iStudy.hasNext()) {
          String studyV = (String) iStudy.next();
          String studyV_L = studyV + "_L";
          String studyV_H = studyV + "_H";
          ArrayList assListLight = new ArrayList();
          ArrayList assListHeavy = new ArrayList();
          ArrayList value = studyGroupMap.get(studyV);
          Iterator iV = value.iterator();
          while (iV.hasNext()) {
            String assName = (String) iV.next();
            if (assName.toLowerCase().contains("light")) {
              assListLight.add(assName);
            } else if (assName.toLowerCase().contains("heavy")) {
              assListHeavy.add(assName);
            }
          }
          studyGroupMap.put(studyV_L, new ArrayList(assListLight));
          studyGroupMap.put(studyV_H, new ArrayList(assListHeavy));
        }
      }
      studyVars = new ArrayList(studyGroupMap.keySet());

            /*
             * multiplicity indicates the number of labels (e.g. heavy vs.
             * light)
             * number of condition groups is number of assay divided by
             * multiplicity
             * the intensity value of each assay location is show below
             * (e.g. multiplicity=2, assay number=4):
             * intensity of all assay + intensity of each label + intensity of
             * group 1
             * + intensity of assay 1 and assay 2 + intensity of group 2 +
             * intensity of assay 3 and assay 4
             */

      int groupNum = 0;
      if (multiplicity != 0) {
        groupNum = primeStudyGroupMap.size();
      } else {
        groupNum = 1;
      }

      /**
       * *
       * read from evidence.txt
       */
      rd = new FileReader(parentPath + "//" + EVIDENCE);
      br = new BufferedReader(rd);
      csvReader = new CSVReader(br, '\t');

            /*
             * read the column title and acquire column number
             */
      int posSeq = 0, posMod = 0, posProt = 0, posRaw = 0, posExp = 0, posChr = 0, posMz = 0, posRet = 0, posInt = 0;
      nextLine = csvReader.readNext();
      for (int i = 0; i < nextLine.length; i++) {
        if (nextLine[i].toLowerCase().equals("sequence")) {
          posSeq = i;
        } else if (nextLine[i].toLowerCase().equals("modifications")) {
          posMod = i;
        } else if (nextLine[i].toLowerCase().equals("leading proteins")) {
          posProt = i;
        } else if (nextLine[i].toLowerCase().equals("raw file")) {
          posRaw = i;
        } else if (nextLine[i].toLowerCase().equals("experiment")) {
          posExp = i;
        } else if (nextLine[i].toLowerCase().equals("charge")) {
          posChr = i;
        } else if (nextLine[i].toLowerCase().equals("m/z")) {
          posMz = i;
        } else if (nextLine[i].toLowerCase().equals("retention time")) {
          posRet = i;
        } else if (nextLine[i].toLowerCase().equals("intensity")) {
          posInt = i;
        }
      }
      // For Java 7 version
//        for (int i=0; i<nextLine.length; i++){
//            switch (nextLine[i].toLowerCase()){
//                case "sequence": posSeq = i; break;
//                case "modifications":  posMod = i; break;
//                case "raw file": posRaw = i; break;
//                case "experiment": posExp = i; break;
//                case "charge": posChr = i; break;
//                case "m/z": posMz = i; break;
//                case "retention time": posRet = i; break;
//                case "intensity": posInt = i; break;
//                default: break;
//            }
//        }

            /*
             * read specific column in to an intermediate HashMap with id as
             * key, an array list as value
             */
      evidenceMap = new HashMap<String, ArrayList<String>>();
      while ((nextLine = csvReader.readNext()) != null) {
        String key = nextLine[0];
                /*
                 * The positions for features are: (0)Sequence,
                 * (1)Modifications, (2)Leading Proteins,
                 * (3)Raw File, (4)Experiment, (5)Charge, (6)m/z, (7)Retention
                 * Time, (8)Intensity
                 * (9)Intensity L, (10)Intensity H
                 */
        List<String> valueList = new ArrayList<String>();
        if (isLabelFree) {
          String[] values = {nextLine[posSeq], nextLine[posMod], nextLine[posProt], nextLine[posRaw], nextLine[posExp],
              nextLine[posChr], nextLine[posMz], nextLine[posRet], nextLine[posInt]};
          valueList = Arrays.asList(values);
        } else {
          String[] values = {nextLine[posSeq], nextLine[posMod], nextLine[posProt], nextLine[posRaw], nextLine[posExp],
              nextLine[posChr], nextLine[posMz], nextLine[posRet], nextLine[posInt], nextLine[posInt + 1], nextLine[posInt + 2]};
          valueList = Arrays.asList(values);
        }
        evidenceMap.put(key, new ArrayList(valueList));
      }
      csvReader.close();


      /**
       * *
       * read from proteinGroups.txt if hasProteinGroupsFile = true
       */
      if (hasProteinGroupsFile) {

        rd = new FileReader(parentPath + "//" + PROTEINS);
        br = new BufferedReader(rd);
        csvReader = new CSVReader(br, '\t');

                /*
                 * cells the column number for required fields
                 */
        int posMajProtIds = 0, posProtInt = 0, posUniqPep = 0;
        nextLine = csvReader.readNext();
        for (int i = 0; i < nextLine.length; i++) {
          if (nextLine[i].toLowerCase().equals("majority protein ids")) {
            posMajProtIds = i;
          } else if (nextLine[i].toLowerCase().equals("intensity")) {
            posProtInt = i;
          } else if (nextLine[i].toLowerCase().contains("unique peptides")) {  // this is the position of the LAST column title containing "unique peptides";
            posUniqPep = i;
          }
        }
        posUniqPep = posUniqPep - assays.size();
        // adjust posUniqPep to point one position AHEAD OF the first title containing "unique peptides"; this trick only works for Label free example as SILAC doesn't need to report "unique peptides"

        proteinToIntMap = new HashMap<String, ArrayList<String>>();
        proteinToUniqPepMap = new HashMap<String, ArrayList<String>>();

        while ((nextLine = csvReader.readNext()) != null) {
          String key = nextLine[posMajProtIds];
          //only take the first ID in Ids as prot id
          if (key.contains(";")) {
            int endIndex = key.indexOf(";");
            key = key.substring(0, endIndex);
          }

                    /*
                     * build a intensity list of each assay on each protein
                     * assay number comes from assays.size()
                     * first id of major protein ids as the key, intensity list
                     * as value
                     */
          ArrayList intList = new ArrayList();

          if (isLabelFree) {
            for (int i = 0; i < assays.size(); i++) {
              intList.add(nextLine[posProtInt + i + 1]);
            }
          } else {
            int i = 0; // i is the count of assays
            for (int j = 0; j < groupNum; j++) {  // j is the count of groups
              int k = 0; // k is the count of labels in each group
              while (k < multiplicity) {
                if (groupNum == 1) {
                  intList.add(nextLine[posProtInt + multiplicity + (i + 1) + j]);
                } else {
                  intList.add(nextLine[posProtInt + multiplicity + (i + 1) + (j + 1)]);
                }
                i++;
                k++;
              }
            }
          }
          proteinToIntMap.put(key, new ArrayList(intList));

                    /*
                     * build a uniqe peptides list of each assay on each protein
                     * this is only for label-free examples, not included in
                     * SILAC examples
                     * first id of major protein ids as the key, unique peptide
                     * list as value
                     */
          if (isLabelFree) {
            ArrayList uniqPepList = new ArrayList();
            for (int j = 0; j < assays.size(); j++) {
              uniqPepList.add(nextLine[posUniqPep + j + 1]);
            }
            proteinToUniqPepMap.put(key, new ArrayList(uniqPepList));
          }
        }
      }


      /**
       * *
       * read from peptides.txt
       */
      rd = new FileReader(parentPath + "//" + PEPTIDES);
      br = new BufferedReader(rd);
      csvReader = new CSVReader(br, '\t');

            /*
             * cells the column number for required fields
             */
      nextLine = csvReader.readNext();
      int posEvd = 0;
      int posPepSeq = 0;
      int posPepInt = 0;

      //Maxquant only reports H/L ratio
      int[] posPepRatio = new int[primeStudyVars.size()];

      //titles of ratio are stored in String[] titlePepRatio
      String[] titlePepRatio = new String[primeStudyVars.size()];

            /*
             * Initialise peptideToRatioMap because
             * the first entry is used to record ratio title,
             * such as "Ratio H/L config_1", set key as "0"
             */
      peptideToRatioMap = new HashMap<String, List<String>>();

      for (int i = 0; i < nextLine.length; i++) {
        if (nextLine[i].toLowerCase().equals("evidence ids")) {
          posEvd = i;
        } else if (nextLine[i].toLowerCase().equals("sequence")) {
          posPepSeq = i;
        } else if (nextLine[i].toLowerCase().equals("proteins")) {
          posProt = i;
        } else if (nextLine[i].toLowerCase().contains("h/l")) {
          for (int k = 0; k < primeStudyVars.size(); k++) {
            if (nextLine[i].toLowerCase().equals("ratio h/l " + (String) primeStudyVars.get(k))) {
              posPepRatio[k] = i;
              titlePepRatio[k] = nextLine[i];
            }
          }
        } else if (nextLine[i].toLowerCase().equals("intensity")) {
          posPepInt = i;
          break;
        }
      }

      peptideToRatioMap.put("0", Arrays.asList(titlePepRatio));

            /*
             * A peptide sequence to feature ids HashMap<String, ArrayList>
             */
      peptideToEvdMap = new HashMap<String, ArrayList<String>>();
      peptideToIntMap = new HashMap<String, ArrayList<String>>();
      peptideToProtMap = new HashMap<String, ArrayList<String>>();
      peptides = new ArrayList();
      while ((nextLine = csvReader.readNext()) != null) {

                /*
                 * build a peptide to evidence ids map, with peptide sequence as
                 * key,
                 * evidence id list as value HashMap<String, ArrayList>
                 */
        String key = nextLine[posPepSeq];
        String[] evidenceIds = nextLine[posEvd].split(";");
        List<String> evidenceIdList = Arrays.asList(evidenceIds);
        peptideToEvdMap.put(key, new ArrayList(evidenceIdList));

                /*
                 * build a intensity list of each assay on each peptide
                 * assay number comes from assays.size()
                 * peptide sequence as key, intensity list as value
                 */
        ArrayList intList = new ArrayList();

        if (isLabelFree) {
          for (int i = 0; i < assays.size(); i++) {
            intList.add(nextLine[posPepInt + i + 1]);
          }
        } else {
          int i = 0; // i is the count of assays
          for (int j = 0; j < groupNum; j++) {  // j is the count of groups
            int k = 0; // k is the count of labels in each group
            while (k < multiplicity) {
              if (groupNum == 1) {
                intList.add(nextLine[posPepInt + multiplicity + (i + 1) + j]);
              } else {
                intList.add(nextLine[posPepInt + multiplicity + (i + 1) + (j + 1)]);
              }
              i++;
              k++;
            }
          }
        }

        peptideToIntMap.put(key, new ArrayList(intList));

                /*
                 * build a peptide sequence to ratio list HashMap
                 * peptide sequence as key, ratio list as value
                 * The first entry used to record ratio title,
                 * such as "Ratio H/L config_1", set key as "0"
                 */
        if (!isLabelFree) {
          List<String> ratioList = peptideToRatioMap.get(key);
          if (ratioList == null) {
            ratioList = new ArrayList<String>();
            for (int i = 0; i < posPepRatio.length; i++) {
              // TODO: treat empty cell as ZERO, correct?
              String ratioValue = (nextLine[posPepRatio[i]].isEmpty() ? "0" : nextLine[posPepRatio[i]]);
              ratioList.add(ratioValue);
            }
          }
          peptideToRatioMap.put(key, ratioList);
        }

                /*
                 * build a peptide sequence to protein accessions HashMap
                 * peptide sequence as key, protein accession list as value
                 */
        String[] proteins = nextLine[posProt].split(";");
        List<String> proteinList = Arrays.asList(proteins);
        peptideToProtMap.put(key, new ArrayList(proteinList));

                /*
                 * build a peptide list
                 */
        peptides.add(key);
      }
      csvReader.close();

            /*
             * build a protein to peptide HashMap,
             * with protein accession as key,
             * and peptide sequence list as value
             */

      proteinToPepMap = mapTranspose(peptideToProtMap);
    }

  }

  public ArrayList<String> getAssayList() {
    return this.assays;
  }

  public HashMap getAssayExpMap() {
    return this.assayExpMap;
  }

  public ArrayList<String> getPrimeStudyVariableList() {
    return this.primeStudyVars;
  }

  public ArrayList<String> getStudyVariableList() {
    return this.studyVars;
  }

  public HashMap getPrimeStudyGroupMap() {
    return this.primeStudyGroupMap;
  }

  public HashMap getStudyGroupMap() {
    return this.studyGroupMap;
  }

  public ArrayList getPeptideList() {
    return this.peptides;
  }

  public ArrayList getRawFileList() {
    return this.rawFiles;
  }

  public HashMap getAssayRawFileMap() {
    return this.assayRawFileMap;
  }

  public HashMap getEvidenceMap() {
    return this.evidenceMap;
  }

  public HashMap getPeptideFeatureIdsMap() {
    return this.peptideToEvdMap;
  }

  public HashMap getPeptideIntensityMap() {
    return this.peptideToIntMap;
  }

  public HashMap getPeptideRatioMap() {
    return this.peptideToRatioMap;
  }

  public HashMap getProteinIntensityMap() {
    return this.proteinToIntMap;
  }

  public HashMap getProteinUniquePeptiedsMap() {
    return this.proteinToUniqPepMap;
  }

  public HashMap getProteinPeptidesMap() {
    return this.proteinToPepMap;
  }

  public Boolean hasProteinGroupsFile() {
    return this.hasProteinGroupsFile;
  }

  public Boolean isLabelFree() {
    return this.isLabelFree;
  }

  public int getLabelNumber() {
    return this.multiplicity;
  }

  private boolean containsAll(String[] list, String[] require) {
    boolean contain = true;
    List<String> aList = Arrays.asList(list);
    for (int i = 0; i < require.length; i++) {
      if (!aList.contains(require[i])) {
        contain = false;
        break;
      }
    }
    return contain;
  }

  private int lengthOfPrefix(String a, String b) {
    int length = (a.length() < b.length()) ? a.length() : b.length();
    int len = 1;
    if ((length == 1) && (a.charAt(0) == b.charAt(0))) {
      len = 1;
    } else if ((a.charAt(0) != b.charAt(0))) {
      len = 0;
    } else {
      len = len + lengthOfPrefix(a.substring(1), b.substring(1));
    }
    return len;
  }

  /*
   * transpose a HashMap<String, String> to a HashMap<String,
   * ArrayList<String>>
   */
  private HashMap<String, ArrayList<String>> mapTranspose(
      HashMap<String, ArrayList<String>> m) {
    HashMap<String, ArrayList<String>> mTr = new HashMap<String, ArrayList<String>>();
    Iterator iM = m.entrySet().iterator();
    while (iM.hasNext()) {
      Map.Entry<String, ArrayList<String>> entry = (Map.Entry<String, ArrayList<String>>) iM.next();
      String key = entry.getKey();
      ArrayList val = entry.getValue();
      Iterator iVal = val.iterator();
      while (iVal.hasNext()) {
        String s = (String) iVal.next();
        ArrayList<String> aL = mTr.get(s);
        if (aL == null) {
          aL = new ArrayList<String>();
          mTr.put(s, aL);
        }
        aL.add(key);
      }

    }
    return mTr;
  }
}
