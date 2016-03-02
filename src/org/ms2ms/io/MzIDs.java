package org.ms2ms.io;

import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;

/** Utils to deal with the mzID file
 *
 */
public class MzIDs
{
  //Some IDs to be used throughout;
  static String siiListID         = "SII_LIST_1";
  static String spectraDataID     = "SID_1";
  static String psiCvID           = "PSI-MS";
  static String siProtocolID      = "SearchProtocol_1";
  static String searchDBID        = "SearchDB_1";
  static String pepEvidenceListID = "PepEvidList_1";
  static String analysisSoftID    = "ID_software";
  static String specIdentID       = "SpecIdent_1";
  static String unimodID          = "UNIMOD";
  static String unitCvID          = "UO";
  static String measureMzID       = "Measure_MZ";
  static String measureIntID      = "Measure_Int";
  static String measureErrorID    = "Measure_Error";
  static String sourceFileID      = "SourceFile_1";

  //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
  //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
  //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>

  static Cv sPsiCV=new Cv(), sUnimodCV=new Cv(), sUnitCV=new Cv();

  static
  {
    sPsiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
    sPsiCV.setId(psiCvID);
    sPsiCV.setVersion("2.25.0");
    sPsiCV.setFullName("PSI-MS");

    sUnimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
    sUnimodCV.setId(unimodID);
    sUnimodCV.setFullName("UNIMOD");

    sUnitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
    sUnitCV.setId(unitCvID);
    sUnitCV.setFullName("UNIT-ONTOLOGY");
  }

  // mzIdentML
  //     cvList
  //     AnalysisSoftwareList
  //     Provider
  //     AuditCollection
  //     AnalysisSampleCollection
  //     SequenceCollection
  //     AnalysisCollection
  //     AnalysisProtocolCollection
  //     DataCollection
  //         Inputs
  //         AnalysisData
  //             SpectrumIdentificationList
  //             ProteinDetectionList
  //         /AnalysisData
  //     /DataCollection
  //     BibliographicReference
  // /mzIdentML

  public static void write(Writer writer, Person owner,
                           AnalysisSampleCollection samples,
                           SequenceCollection sequences,
                           AnalysisCollection analysis,
                           AnalysisProtocolCollection protocols,
                           Inputs inputs,
                           SpectrumIdentificationList psms) throws Exception
  {
    MzIdentMLMarshaller m = new MzIdentMLMarshaller();

    // Note: writing of '\n' characters is optional and only for readability of the produced XML document
    // Also note: since the XML is produced in individual parts, the overall formatting of the document
    //            is not as nice as it would be when marshalling the whole structure at once.

    // XML header
    writer.write(m.createXmlHeader() + "\n");

    // mzIdentML start tag
    writer.write(m.createMzIdentMLStartTag("12345") + "\n");

    // create the software list
    m.marshall(newCvList(sPsiCV, sUnimodCV, sUnitCV),             writer); writer.write("\n");
    m.marshal(addAnalysisSoftware(null, newConvertor()),          writer); writer.write("\n");
    m.marshall(newProvider(owner),                                writer); writer.write("\n");
    m.marshall(newAuditCollection(owner,"address","affiliation"), writer); writer.write("\n");
    m.marshall(samples,                                           writer); writer.write("\n");
    m.marshall(sequences,                                         writer); writer.write("\n");
    m.marshall(analysis,                                          writer); writer.write("\n");
    m.marshall(protocols,                                         writer); writer.write("\n");

    writer.write(m.createDataCollectionStartTag() + "\n");

    m.marshall(inputs,                                            writer); writer.write("\n");

    writer.write(m.createAnalysisDataStartTag()   + "\n");

    m.marshall(psms, writer); writer.write("\n");

    writer.write(m.createProteinDetectionListStartTag("PDL_1", null) + "\n");
    writer.write(m.createProteinDetectionListClosingTag()            + "\n");
    writer.write(m.createAnalysisDataClosingTag()                    + "\n");
    writer.write(m.createDataCollectionClosingTag()                  + "\n");
    writer.write(m.createMzIdentMLClosingTag());

    writer.close();
  }
  public static CvList newCvList(Cv... cvs)
  {
    CvList cvList = new CvList();
    List<Cv> localCvList = cvList.getCv();

    if (Tools.isSet(cvs))
      for (Cv cv : cvs) localCvList.add(cv);

    return cvList;
  }
  /**
   *
   * Aim is to write out set up the analysisSoftwareList following this
   * structure: <AnalysisSoftware id="ID_software" name="xtandem"
   * version="2008.12.1.1" > <SoftwareName> <cvParam accession="MS:1001476"
   * name="xtandem" cvRef="PSI-MS" /> </SoftwareName>
   *
   */
  public static AnalysisSoftwareList addAnalysisSoftware(AnalysisSoftwareList analysisSoftwareList, AnalysisSoftware... softwares)
  {
    if (analysisSoftwareList==null) analysisSoftwareList = new AnalysisSoftwareList();

    List<AnalysisSoftware>  analysisSoftwares = analysisSoftwareList.getAnalysisSoftware();
        /*
         * TO DO - need to work out how to use Param CvParam cvParam = new
         * CvParam(); cvParam.setName("xtandem"); cvParam.setCvRef(psiCvID);
         * cvParam.setAccession("MS:1001476"); ParamAlternative paramAlt = new
         * ParamAlternative(); paramAlt.setCvParam(cvParam);
         *
         * analysisSoftware.setSoftwareName(makeCvParam("MS:1001476","xtandem",psiCV));
         * analysisSoftware.setSoftwareName(paramAlt);
         */
    if (Tools.isSet(softwares))
      for (AnalysisSoftware software : softwares) analysisSoftwares.add(software);

    return analysisSoftwareList;
  }

//  public void handleCVs()
//  {
//    //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
//    //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
//    //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>
//
//    CvList cvList = new CvList();
//    List<Cv> localCvList = cvList.getCv();
//    Cv psiCV = new Cv();
//
//    psiCV.setUri("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
//    psiCV.setId(psiCvID);
//    psiCV.setVersion("2.25.0");
//    psiCV.setFullName("PSI-MS");
//
//    Cv unimodCV = new Cv();
//    unimodCV.setUri("http://www.unimod.org/obo/unimod.obo");
//    unimodCV.setId(unimodID);
//    unimodCV.setFullName("UNIMOD");
//
//    Cv unitCV = new Cv();
//    unitCV.setUri("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
//    unitCV.setId(unitCvID);
//    unitCV.setFullName("UNIT-ONTOLOGY");
//
//    localCvList.add(psiCV);
//    localCvList.add(unimodCV);
//    localCvList.add(unitCV);
//  }
  /**
   * Helper method to create and return a CvParam from accession, name and CV
   *
   * @return CvParam
   */
  public static CvParam makeCvParam(String accession, String name, Cv cv)
  {
    CvParam cvParam = new CvParam();
    cvParam.setAccession(accession);
    cvParam.setName(name);
    cvParam.setCv(cv);
    return cvParam;
  }

  /**
   * Helper method to create and return a CvParam from accession, name and CV
   *
   * @return CvParam
   */
  public static CvParam makeCvParam(String accession, String name, Cv cv, String value)
  {
    CvParam cvParam = new CvParam();
    cvParam.setAccession(accession);
    cvParam.setName(name);
    cvParam.setCv(cv);
    cvParam.setValue(value);
    return cvParam;
  }

  /**
   * Helper method to create and return a CvParam from accession, name, CV,
   * unitAccession and unitName (unitCV is automatically provided)
   *
   * @return CvParam
   */
//  public static CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName)
//  {
//    CvParam cvParam = new CvParam();
//    cvParam.setAccession(accession);
//    cvParam.setName(name);
//    cvParam.setCv(cv);
//    cvParam.setUnitAccession(unitAccession);
//    cvParam.setUnitCv(unitCV);
//    cvParam.setUnitName(unitName);
//    return cvParam;
//  }

  /**
   * Helper method to create and return a CvParam from accession, name, CV,
   * unitAccession, unitName and unitCV
   *
   * @return CvParam
   */
  public static CvParam makeCvParam(String accession, String name, Cv cv, String unitAccession, String unitName, Cv alternateUnitCV)
  {
    CvParam cvParam = new CvParam();
    cvParam.setAccession(accession);
    cvParam.setName(name);
    cvParam.setCv(cv);
    cvParam.setUnitAccession(unitAccession);
    cvParam.setUnitCv(alternateUnitCV);
    cvParam.setUnitName(unitName);
    return cvParam;
  }
  public static AnalysisSoftware newXTandem(String version)
  {
    AnalysisSoftware analysisSoftware = new AnalysisSoftware();
    analysisSoftware.setName("xtandem");

    Param tempParam = new Param();
    tempParam.setParam(makeCvParam("MS:1001476", "X\\!Tandem", sPsiCV));
    analysisSoftware.setSoftwareName(tempParam);

    analysisSoftware.setId(analysisSoftID);
    analysisSoftware.setVersion(version);

    return analysisSoftware;
  }
  public static AnalysisSoftware newConvertor()
  {
    AnalysisSoftware analysisSoftware = new AnalysisSoftware();
    Date date = new Date() ;
    SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH-mm-ss") ;
    //TODO: this seems strange...a date in the id and name fields? Check if there
    //is perhaps a date field where the date information could be stored.
    analysisSoftware.setName("MzIDs_"+dateFormat.format(date)); analysisSoftware.setId("MzIDs_"+dateFormat.format(date));
    Param param = new Param();
    param.setParam(makeCvParam("MS:1002237", "mzidLib", sPsiCV));
    analysisSoftware.setSoftwareName(param);

    return analysisSoftware;
  }
  /**
   * Setup Provider element as follows <Provider id="PROVIDER"> <ContactRole
   * Contact_ref="PERSON_DOC_OWNER"> <role> <cvParam accession="MS:1001271"
   * name="researcher" cvRef="PSI-MS"/> </role> </ContactRole> </Provider>
   *
   */
  public static Provider newProvider(Person docOwner)
  {
    Provider provider = new Provider();
    provider.setId("PROVIDER");

    ContactRole contactRole = new ContactRole();
    contactRole.setContact(docOwner);

    Role role = new Role();
    role.setCvParam(makeCvParam("MS:1001271", "researcher", sPsiCV));
    contactRole.setRole(role);

    provider.setContactRole(contactRole);

    return provider;
  }

  public static Person newPerson(String first, String last, String address)
  {
    Person docOwner = new Person();
    docOwner.setId("PERSON_DOC_OWNER");
    docOwner.setFirstName(first);
    docOwner.setLastName(last);
    docOwner.getCvParam().add(makeCvParam("MS:1000587", "contact address", sPsiCV, address));

    return docOwner;
  }
  /**
   * TO DO Capture name and email of the user <AuditCollection> <Person
   * id="PERSON_DOC_OWNER" firstName="Andy" lastName="Jones"
   * email="someone@someuniversity.com"> <affiliations
   * Organization_ref="ORG_DOC_OWNER"/> </Person> <Organization
   * id="ORG_DOC_OWNER" address="Some address" name="Some place" />
   * </AuditCollection>
   *
   *
   */
  public static AuditCollection newAuditCollection(Person docOwner, String address, String affiliationName)
  {
    AuditCollection auditCollection = new AuditCollection();
    //List<Contact> contactList = auditCollection.getContactGroup();
    List<AbstractContact> contactList = auditCollection.getPersonOrOrganization();

    Organization org = new Organization();
    org.setId("ORG_DOC_OWNER");
    org.setName(affiliationName);
    org.getCvParam().add(makeCvParam("MS:1000586", "contact name", sPsiCV, address));
    //org.setAddress(address);

    List<Affiliation> affList = docOwner.getAffiliation();
    Affiliation aff = new Affiliation();
    aff.setOrganization(org);
    affList.add(aff);
    contactList.add(docOwner);
    contactList.add(org);

    return auditCollection;
  }
}
