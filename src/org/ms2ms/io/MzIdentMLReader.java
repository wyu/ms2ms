package org.ms2ms.io;

import org.expasy.mzjava.core.ms.spectrum.RetentionTimeList;
import org.expasy.mzjava.core.ms.spectrum.ScanNumberList;
import org.expasy.mzjava.proteomics.io.ms.ident.MzIdentMlReader;
import org.expasy.mzjava.proteomics.io.ms.ident.NamespaceFilterXMLReader;
import org.expasy.mzjava.proteomics.io.ms.ident.PSMReaderCallback;
import org.expasy.mzjava.proteomics.io.ms.ident.mzidentml.v110.*;
import org.expasy.mzjava.proteomics.mol.AAMassCalculator;
import org.expasy.mzjava.proteomics.mol.modification.ModificationResolver;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.math.Stats;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.sax.SAXSource;
import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by yuw on 10/10/2015.
 */
public class MzIdentMLReader extends MzIdentMlReader
{
  public static Map<String, String> mScoreMap = new HashMap<>();

  static
  {
    mScoreMap.put("MS:1001330", "X\\!Tandem:expect");
    mScoreMap.put("MS:1001331", "X\\!Tandem:hyperscore");
    mScoreMap.put("MS:1001328", "OMSSA:evalue");
    mScoreMap.put("MS:1001329", "OMSSA:pvalue");
    mScoreMap.put("MS:1002049", "MS-GF:RawScore");
    mScoreMap.put("MS:1002050", "MS-GF:DeNovoScore");
    mScoreMap.put("MS:1002052", "MS-GF:SpecEValue");
    mScoreMap.put("MS:1002053", "MS-GF:EValue");
    mScoreMap.put("MS:1002054", "MS-GF:QValue");
    mScoreMap.put("MS:1002055", "MS-GF:PepQValue");
    mScoreMap.put("MS:1001171", "Mascot:score");
    mScoreMap.put("MS:1001172", "Mascot:expectation value");
    mScoreMap.put("MS:1002012", "Mascot:delta score");
    mScoreMap.put("MS:1001589", "MyriMatch:MVH");
    mScoreMap.put("MS:1001590", "MyriMatch:mzFidelity");
    mScoreMap.put("MS:1001121", "number of matched peaks");
    mScoreMap.put("MS:1001362", "number of unmatched peaks");
    mScoreMap.put("user:xcorr", "xcorr");
//    mScoreMap.put("", "");
  }
  public MzIdentMLReader(ModificationResolver modResolver) { super(modResolver); }

  public void parse(InputSource is, PSMReaderCallback callback) throws SAXException, JAXBException
  {
    //peptideRefMap = null;

    Class<MzIdentMLType> docClass = MzIdentMLType.class;

    XMLReader reader = XMLReaderFactory.createXMLReader();
    reader.setFeature("http://apache.org/xml/features/allow-java-encodings", true);

    NamespaceFilterXMLReader filter = new NamespaceFilterXMLReader("http://psidev.info/psi/pi/mzIdentML/1.1");
    filter.setParent(reader);

    SAXSource ss = new SAXSource(filter, is);

    String packageName = docClass.getPackage().getName();
    JAXBContext jc     = JAXBContext.newInstance(packageName);
    Unmarshaller u     = jc.createUnmarshaller();
    MzIdentMLType doc  = (MzIdentMLType) ((JAXBElement) u.unmarshal(ss)).getValue();

    readPeptideIDProteinMap(doc);
    readPeptideIDProteinPositionMap(doc);

    SpectraDataType spectraDataType = doc.getDataCollection().getInputs().getSpectraData().get(0);

    String msFileLocation = spectraDataType.getLocation();

    String fileName = null;
    if (msFileLocation != null)
    {
      File tmp = new File(msFileLocation);
      fileName = tmp.getName();
    }

    for (SpectrumIdentificationListType identList : doc.getDataCollection().getAnalysisData().getSpectrumIdentificationList())
    {
      for (SpectrumIdentificationResultType idResult : identList.getSpectrumIdentificationResult())
      {
        Map<String, AbstractParamType> idResultCVParamMap = MsReaders.toCVMap(idResult.getParamGroup());

        ScanNumberList       scanNumbers = MsReaders.newScanNumberList(idResultCVParamMap);
        RetentionTimeList retentionTimes = MsReaders.newRetentionTimeList(idResultCVParamMap);
//        int                        index = Stats.toInt(idResult.getSpectrumID());
        String                     title = idResultCVParamMap.containsKey(TITLE_CV)?idResultCVParamMap.get(TITLE_CV).getValue():"unTitled";

        for (SpectrumIdentificationItemType idItem : idResult.getSpectrumIdentificationItem())
        {
          SpectrumIdentifier identifier = new SpectrumIdentifier(title);

          identifier.setSpectrumFile(fileName);
          identifier.addRetentionTime(retentionTimes);
          identifier.addScanNumber(scanNumbers);
//          identifier.setIndex(index);
          identifier.setAssumedCharge(idItem.getChargeState());
          identifier.setPrecursorNeutralMass(AAMassCalculator.getInstance().calculateNeutralMolecularMass(idItem.getExperimentalMassToCharge(), idItem.getChargeState()));

          PeptideMatch peptideMatch = readPeptide(idItem.getPeptideRef(), doc.getSequenceCollection());

          Map<String, AbstractParamType> idItemCVParamMap = MsReaders.toCVMap(idItem.getParamGroup());

          for (String cv : mScoreMap.keySet())
            if (idItemCVParamMap.containsKey(cv))
              peptideMatch.addScore(mScoreMap.get(cv), Stats.toDouble(idItemCVParamMap.get(cv).getValue()));

          peptideMatch.setRank(idItem.getRank());
          peptideMatch.setMassDiff(idItem.getExperimentalMassToCharge() - idItem.getCalculatedMassToCharge());
          peptideMatch.setRejected(false);

          callback.resultRead(identifier, peptideMatch);
        }
      }
    }
  }
  @Override
  protected void readPeptideIDProteinPositionMap(MzIdentMLType doc)
  {
    // fix the syntax problem before calling the super()
    SequenceCollectionType seqCollection = doc.getSequenceCollection();
    if (seqCollection == null) return;

    List<PeptideEvidenceType> peptideEvidences = seqCollection.getPeptideEvidence();
    if (peptideEvidences == null || peptideEvidences.isEmpty()) return;

    for (PeptideEvidenceType entry : peptideEvidences)
    {
      if (entry.getPost().length()>1) entry.setPost(entry.getPost().substring(0, 1));
      if (entry.getPre( ).length()>1) entry.setPre( entry.getPre( ).substring(0, 1));
    }

    super.readPeptideIDProteinPositionMap(doc);
  }
}
