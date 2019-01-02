package org.ms2ms.test.io;

import com.google.common.collect.Multimap;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.junit.Test;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.io.MsIO;
import org.ms2ms.io.MsReaders;
import org.ms2ms.io.MzIDs;
import org.ms2ms.test.TestAbstract;
import uk.ac.ebi.jmzidml.MzIdentMLElement;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLUnmarshaller;
import uk.ac.ebi.jmzml.model.mzml.BinaryDataArray;
import uk.ac.ebi.jmzml.model.mzml.CVParam;
import uk.ac.ebi.jmzml.model.mzml.Spectrum;
import uk.ac.ebi.jmzml.xml.io.MzMLObjectIterator;
import uk.ac.ebi.jmzml.xml.io.MzMLUnmarshaller;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

/**
 * Created by yuw on 10/2/2015.
 *
 * adapted from https://code.google.com/p/psi-standard-formats-tutorial/source/browse/trunk/src/main/java/standards/TestGUI.java
 */
public class PsiStdTest extends TestAbstract
{
  @Test
  public void writeMzIDviaPeptideMatches() throws Exception
  {
    String data = "/Users/yuw/Documents/Data/Joslin/PSMs/default/Mouse_Plasma_LIRKO2_24_27Aug12_Lynx_12-07-02.omssa.csv";

    // produce a list of peptide matches
    Multimap<SpectrumIdentifier, PeptideMatch> spec_match = Engine.OMSSA.readPSMs(data, 1);

    // a basic setup to write out an mzID file from a set of PeptideMatches
    MzIDs  mzid = new MzIDs("data/TMT10QE.param").addSourceFile(data, "CSV");

    // build the data set
    mzid.buildPeptideMatches(spec_match);
    // write the output
    mzid.write(data+".mzid");
  }
  @Test
  public void testSplitMS2MS3() throws Exception
  {
    //String   root = "C:/local/data/TMT_MS3/TMT_MS3120_1800_40CE";
//    String   root = "C:/local/data/TMT_MS3/TMT_MS3120_1800_55CE";

//    File[] files = new File("/Users/yuw/Documents/Data/MS2MS3/mzML").listFiles();
//    File[] files = new File("/Users/yuw/Documents/Data/MS2MS3/mzML150").listFiles();
//    File[] files = new File("/Users/yuw/Documents/Data/MS2MS3/mzML.refs").listFiles();
    File[] files = new File("/Users/yuw/Documents/Data/MS2MS3/mzMLCHO").listFiles();

//    File[] files = new File("/Volumes/PPP_Data/Fusion-03 Data/151107_Ecoli_12protein_Spike_MS2_MS3/working/").listFiles();
    for (File mzML : files)
    {
      if (mzML.getName().indexOf(".mzML")<=0) continue;
      String mzml_root = mzML.getCanonicalPath().substring(0, mzML.getCanonicalPath().lastIndexOf("."));
      System.out.println("Splitting "+mzml_root);
      MsReaders.splitMs2MS3(mzml_root);
    }
  }

  @Test
  public void readSpectra() throws Exception
  {
    //create a new unmarshaller object, you can use a URL or a File to initialize the unmarshaller
    File xmlFile = new File("C:\\local\\data\\TMT_MS3\\TMT_MS3120_1800_40CE.mzML");
    MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(xmlFile);

    // looping through the scans
    MzMLObjectIterator<Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);
    while (spectrumIterator.hasNext())
    {
      Spectrum spectrum = spectrumIterator.next();
      String mslevel = MsIO.get(spectrum.getCvParam(), "MS:1000511");
      String rt = MsIO.getIf(spectrum.getScanList().getScan().get(0).getCvParam(), "MS:1000016", "UO:0000031");
    }
  }
  public void mzIdentMLExample()
  {
    File file = new File("");
    MzIdentMLUnmarshaller unmarshaller = new MzIdentMLUnmarshaller(file);
    HashMap<String, Peptide> peptideIdHashMap = new HashMap();
    HashMap<String, PeptideEvidence> peptideEvidHashMap = new HashMap();

    Iterator<Peptide> iterPeptide = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.Peptide);
    while (iterPeptide.hasNext()) {
      Peptide peptide = iterPeptide.next();
      peptideIdHashMap.put(peptide.getId(), peptide);
    }

    Iterator<PeptideEvidence> iterPeptideEvidence = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.PeptideEvidence);
    while (iterPeptideEvidence.hasNext()) {
      PeptideEvidence peptideEvidence = iterPeptideEvidence.next();
      peptideEvidHashMap.put(peptideEvidence.getId(), peptideEvidence);
    }

    HashMap<String, DBSequence> dbseqMap = new HashMap();
    Iterator<DBSequence> iterDbSeq = unmarshaller.unmarshalCollectionFromXpath(MzIdentMLElement.DBSequence);
    while (iterDbSeq.hasNext())
    {
      DBSequence dbseq = iterDbSeq.next();
      dbseqMap.put(dbseq.getId(), dbseq);
    }

    AnalysisData ad = unmarshaller.unmarshal(DataCollection.class).getAnalysisData();
    List<SpectrumIdentificationList> sil = ad.getSpectrumIdentificationList();
    SpectrumIdentificationList sIdentList = sil.get(0);

    for (SpectrumIdentificationResult spectrumIdentResult : sIdentList.getSpectrumIdentificationResult())
    {
      String spectrumID = spectrumIdentResult.getSpectrumID();
      for (SpectrumIdentificationItem sii : spectrumIdentResult.getSpectrumIdentificationItem())
      {
        int rank = sii.getRank();
        int charge = sii.getChargeState();
        Peptide pep = peptideIdHashMap.get(sii.getPeptideRef());
//        model.insertRow(model.getRowCount(), new Object[]{spectrumID,
//            pep.getPeptideSequence(),
//            rank,
//            charge, ""});

        for (PeptideEvidenceRef pepEvidRef : sii.getPeptideEvidenceRef())
        {
          PeptideEvidence pepEvid = peptideEvidHashMap.get(pepEvidRef.getPeptideEvidenceRef());
          DBSequence dbSeq = dbseqMap.get(pepEvid.getDBSequenceRef());
//          model.insertRow(model.getRowCount(), new Object[]{"",
//              "",
//              "",
//              "",
//              dbSeq.getAccession()});
        }
      }
    }
  }
  public void XIC() throws Exception
  {
    double frtStarts = 0;
    double frtEnds = 1000;
    double fmzStarts = 0;
    double fmzEnds = 10000;

    File file = new File("");
    MzMLUnmarshaller unmarshaller = new MzMLUnmarshaller(file);

    //... Iterate in the rt dimension ...//
    MzMLObjectIterator<Spectrum> spectrumIterator = unmarshaller.unmarshalCollectionFromXpath("/run/spectrumList/spectrum", Spectrum.class);
    while (spectrumIterator.hasNext())
    {
      Spectrum spectrum = spectrumIterator.next();

      //... Identify MS1 data ...//
      String mslevel = "";
      List<CVParam> specParam = spectrum.getCvParam();
      for (Iterator lCVParamIterator = specParam.iterator(); lCVParamIterator.hasNext(); )
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals("MS:1000511"))
        {
          mslevel = lCVParam.getValue().trim();
        }
      }
      if (mslevel.equals("1"))
      {
        //... Get rt from spectrum ...//
        double rt = 0.0;
        String unitRT = "";
        List<CVParam> scanParam = spectrum.getScanList().getScan().get(0).getCvParam();
        for (Iterator lCVParamIterator = scanParam.iterator(); lCVParamIterator.hasNext(); )
        {
          CVParam lCVParam = (CVParam) lCVParamIterator.next();
          if (lCVParam.getAccession().equals("MS:1000016"))
          {
            unitRT = lCVParam.getUnitAccession().trim();
            if (unitRT.equals("UO:0000031")) //... Validating rt unit (mins or secs) ...//
            {
              rt = Float.parseFloat(lCVParam.getValue().trim()) * 60;
            } else {
              rt = Float.parseFloat(lCVParam.getValue().trim());
            }
          }
        }
        //... Get XIC across intervals ...//
        if (rt >= frtStarts && rt <= frtEnds)
        {
          Number[] mzNumbers = null;
          Number[] intenNumbers = null;
          //... Reading mz Values ...//
          List<BinaryDataArray> bdal = spectrum.getBinaryDataArrayList().getBinaryDataArray();
          for (BinaryDataArray bda : bdal)
          {
            List<CVParam> cvpList = bda.getCvParam();
            for (CVParam cvp : cvpList)
            {
              if (cvp.getAccession().equals("MS:1000514"))
              {
                mzNumbers = bda.getBinaryDataAsNumberArray();
              }
              if (cvp.getAccession().equals("MS:1000515"))
              {
                intenNumbers = bda.getBinaryDataAsNumberArray();
              }
            }
          }
          //... Generating XIC ...//
          double intensXIC = 0.0;
          if (mzNumbers != null)
          {
            for (int iI = 0; iI < mzNumbers.length; iI++)
            {
              if (mzNumbers[iI].doubleValue() >= fmzStarts && mzNumbers[iI].doubleValue() <= fmzEnds)
              {
                intensXIC += intenNumbers[iI].doubleValue();
              }
            }
//            model.insertRow(model.getRowCount(), new Object[]{
//                Float.parseFloat(String.format("%.2f", rt)),
//                Float.parseFloat(String.format("%.4f", intensXIC))});
          }
        }
      }
    }
  }
}
