package org.ms2ms.io;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import org.expasy.mzjava.proteomics.mol.modification.*;
import org.expasy.mzjava.proteomics.mol.modification.unimod.UnimodMod;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.algo.PSMs;
import org.ms2ms.data.ms.Engine;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzidml.model.mzidml.*;
import uk.ac.ebi.jmzidml.model.mzidml.Modification;
import uk.ac.ebi.jmzidml.xml.io.MzIdentMLMarshaller;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.*;

/** The class to create an instance of mzID file. To get started, create the object and populate the required parameters.
 *  Then call the "write" method with key IDs to complete the output.
 *
 */
public class MzIDs
{
  //Some IDs to be used throughout;
  static final String siiListID         = "SII_LIST_1";
  static final String psiCvID           = "PSI-MS";
  static final String siProtocolID      = "SearchProtocol_1";
  static final String searchDBID        = "SearchDB_1";
  static final String pepEvidenceListID = "PepEvidList_1";
  static final String analysisSoftID    = "ID_software";
  static final String specIdentID       = "SpecIdent_1";
  static final String unimodID          = "UNIMOD";
  static final String unitCvID          = "UO";
  static final String measureMzID       = "Measure_MZ";
  static final String measureIntID      = "Measure_Int";
  static final String measureErrorID    = "Measure_Error";
  static final String sourceFileID      = "SourceFile_1";

  static final String TYPE_PARENT_MASS  = "Parent mass type";
  static final String TYPE_FRAG_MASS    = "Fragment mass type";
  static final String TOL_FRAG_PLUS     = "Fragment search tolerance plus";
  static final String TOL_FRAG_MINUS    = "Fragment search tolerance minus";
  static final String TOL_PARENT_PLUS   = "Parent search tolerance plus";
  static final String TOL_PARENT_MINUS  = "Parent search tolerance minus";
  static final String THLD_PSM          = "PSM threshold";
  static final String FMT_FILE_INPUT    = "Input file format";
  static final String FMT_FILE_DB       = "Database file format";
  static final String FMT_DATA_SPEC     = "Spectra data file format";
  static final String FMT_ID_SPEC       = "Spectrum ID format";

  static final String ENZYME            = "Enzyme";

  static final String PATH_DB_LOCAL     = "Local database path";
  static final String SPEC_SEARCHED     = "Searched spectrum";
  static final String VER_SOFTWARE      = "Software version";
  static final String MISSED_CLEAVAGE   = "Missed cleavages";

  static final String NAME_DB           = "Database name";
  static final String NAME_SOFTWARE     = "Software name";
  static final String NAME_FIRST        = "File contact first name";
  static final String NAME_LAST         = "File contact last name";
  static final String NAME_ORG          = "File contact organization name";
  static final String ADDRESS           = "File contact address";
  static final String DECOY_COMP        = "Decoy database composition";
  static final String DECOY_REGEX       = "Decoy database regex";
  static final String DECOY_TYPE        = "Decoy database type";

  static final String CAT_SEARCH        = "SEARCH";
  static final String CAT_MOD           = "MOD";
  static final String CAT_SCORE         = "SCORE";

  //<cv id="PSI-MS" fullName="PSI-MS" URI="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo" version="2.25.0"/>
  //<cv id="UNIMOD" fullName="UNIMOD" URI="http://www.unimod.org/obo/unimod.obo" />
  //<cv id="UO" fullName="UNIT-ONTOLOGY" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"></cv>

  // system wide parameters
  static Cv sPsiCV=new Cv(), sUnimodCV=new Cv(), sUnitCV=new Cv();

  static
  {
    sPsiCV    = newCV("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo",psiCvID,"PSI-MS","2.25.0");
    sUnimodCV = newCV("http://www.unimod.org/obo/unimod.obo", unimodID,"UNIMOD");
    sUnitCV   = newCV("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo",unitCvID,"UNIT-ONTOLOGY");
  }

  // score mapping
  private Table<String, String, String>  mParams=null;
  private Map<String, CvParam>           mScoreCV=null;
  private Person                         mOwner=newPerson("Dummy","Dumb","One Main St, Somewhereville, OO");
  private AnalysisSoftwareList           mAnalysisSoftwareList=null;
  private AuditCollection                mAuditCollection=null;
  private Provider                       mProvider=null;
  private SpectrumIdentificationProtocol mSiProtocol=null;
  private SearchDatabase                 mSearchDB=null;
  private AnalysisProtocolCollection     mAnalysisProtocolCollection=null;
  private Inputs                         mInputs=null;
  private AnalysisSampleCollection       mSamples=null;
  private SpectraData                    mSpectraData=null;
  private SequenceCollection             mSequences=null;
  private AnalysisCollection             mAnalysis=null;
  private DataCollection                 mData=null;
  private AnalysisData                   mAnalysisData=null;
  private SpectrumIdentificationList     mSiList=null;
  private Engine                         mEngine=null;
  private Map<String, SpectraData>       mRun2SpectraData=new HashMap<>();
  private int                            mSpecDataCounter=1;
  private String                         mSpectraDataID="SID_1";
  private SpectrumIdentification         mSpecIdent=null;

  public MzIDs(String param_File) { init(param_File); }

  public void init(String param_file)
  {
    init();

    mParams = readParams(param_file);
    // make an Inputs object first
    mInputs = newInputs();

    mSiProtocol = newSpectrumIdentificationProtocol();
    addAnalysisProtocol(mSiProtocol);

    // add the convertor as the first
    addAnalysisSoftware(newConvertor());
    // parse the software names
    String[] softwares = Strs.split(getParam(CAT_SEARCH, NAME_SOFTWARE), ';', true);
    if (Tools.isSet(softwares))
      for (String software : softwares)
      {
        mEngine = Engine.valueOf(Tools.back(Strs.split(software, '#', true)));
        addAnalysisSoftware(mEngine!=null?newEngine(mEngine):newSoftware(software));
      }

    newAnalysisCollection();
  }
  public void init()
  {
    addAudit(   mOwner, "self");
    addProvider(mOwner);
  }
  public String getParam(String cat, String key) { return mParams!=null?mParams.get(cat,key):null; }
  public Table<String,String,String> getParams() { return mParams; }
  public CvParam getScoreCV(String s) { return mScoreCV!=null?mScoreCV.get(s):null; }
  public Engine getEngine() { return mEngine; }
  public SpectraData getSpectraData() { return mSpectraData; }
  public SearchDatabase getSearchDB() { return mSearchDB; }
  public SpectrumIdentificationProtocol getSiProtocol() { return mSiProtocol; }
  public SpectrumIdentificationList getSpecIDs() { return mSiList; }
  public AnalysisProtocolCollection getAnalysisProtocols() { return mAnalysisProtocolCollection; }
  public AnalysisCollection getAnalysisCollection() { return mAnalysis; }
  public AnalysisSampleCollection getSamples() { return mSamples; }
  public SequenceCollection getSequenceCollection() { return mSequences; }
  public Inputs getInputs() { return mInputs; }
  public Person getOwner() { return mOwner; }
  public MzIDs  setOwner(String first, String last, String address) { mOwner=newPerson(first, last, address); return this; }

  public void write(String out)
  {
    Writer w=null;
    try
    {
      try
      {
        w = new FileWriter(out);
        write(w);
      }
      finally { if (w!=null) w.close(); }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  public void write(Writer writer) throws IOException
  {
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

    MzIdentMLMarshaller m = new MzIdentMLMarshaller();

    // Note: writing of '\n' characters is optional and only for readability of the produced XML document
    // Also note: since the XML is produced in individual parts, the overall formatting of the document
    //            is not as nice as it would be when marshalling the whole structure at once.

    // XML header
    writer.write(m.createXmlHeader() + "\n");

    // mzIdentML start tag
    writer.write(m.createMzIdentMLStartTag("12345") + "\n");

    // create the software list
    m.marshall(newCvList(sPsiCV, sUnimodCV, sUnitCV), writer); writer.write("\n");
    m.marshal(mAnalysisSoftwareList,                  writer); writer.write("\n");
    m.marshall(mProvider,                             writer); writer.write("\n");
    m.marshall(mAuditCollection,                      writer); writer.write("\n");
    // This part is optional in the file - not yet completed
//    m.marshall(getSamples(),                               writer); writer.write("\n");
    m.marshall(getSequenceCollection(),               writer); writer.write("\n");
    m.marshall(getAnalysisCollection(),               writer); writer.write("\n");
    m.marshall(getAnalysisProtocols(),                writer); writer.write("\n");

    writer.write(m.createDataCollectionStartTag()        +"\n"); m.marshall(getInputs(),  writer); writer.write("\n");
    writer.write(m.createAnalysisDataStartTag()          +"\n");

    m.marshall(newFragmentationTable(),               writer); writer.write("\n");
    m.marshall(getSpecIDs(),                          writer); writer.write("\n");
//    writer.write(m.createProteinDetectionListStartTag("PDL_1",null)+"\n");
//    writer.write(m.createProteinDetectionListClosingTag()+"\n");
    writer.write(m.createAnalysisDataClosingTag()        +"\n");
    writer.write(m.createDataCollectionClosingTag()      +"\n");
    writer.write(m.createMzIdentMLClosingTag());

    writer.close();
  }
  //    <Measure id="Measure_MZ">
  //    <cvParam accession="MS:1001225" cvRef="PSI-MS" unitCvRef="PSI-MS" unitName="m/z"
  //    unitAccession="MS:1000040" name="product ion m/z"/>
  //    </Measure>
  public static FragmentationTable newFragmentationTable()
  {
    FragmentationTable tbl = new FragmentationTable();

    Measure M = new Measure();
    CvParam cv = newCvParam("MS:1001225", "product ion m/z");
    cv.setUnitAccession("MS:1000040");
    cv.setUnitName("m/z");
    M.getCvParam().add(cv);
    tbl.getMeasure().add(M);
    return tbl;
  }
//  public void write(Writer writer,
//                           AnalysisSampleCollection samples,
//                           SequenceCollection sequences,
//                           AnalysisCollection analysis,
//                           AnalysisProtocolCollection protocols,
//                           Inputs inputs,
//                           SpectrumIdentificationList psms) throws Exception
//  {
//    MzIdentMLMarshaller m = new MzIdentMLMarshaller();
//
//    // Note: writing of '\n' characters is optional and only for readability of the produced XML document
//    // Also note: since the XML is produced in individual parts, the overall formatting of the document
//    //            is not as nice as it would be when marshalling the whole structure at once.
//
//    // XML header
//    writer.write(m.createXmlHeader() + "\n");
//
//    // mzIdentML start tag
//    writer.write(m.createMzIdentMLStartTag("12345") + "\n");
//
//    // create the software list
//    m.marshall(newCvList(sPsiCV, sUnimodCV, sUnitCV), writer); writer.write("\n");
//    m.marshal(mAnalysisSoftwareList,                  writer); writer.write("\n");
//    m.marshall(mProvider,                             writer); writer.write("\n");
//    m.marshall(mAuditCollection,                      writer); writer.write("\n");
//    // This part is optional in the file - not yet completed
//    m.marshall(samples,                               writer); writer.write("\n");
//    m.marshall(sequences,                             writer); writer.write("\n");
//    m.marshall(analysis,                              writer); writer.write("\n");
//    m.marshall(protocols,                             writer); writer.write("\n");
//
//    writer.write(m.createDataCollectionStartTag() + "\n");
//
//    m.marshall(inputs,                                       writer); writer.write("\n");
//
//    writer.write(m.createAnalysisDataStartTag()   + "\n");
//
//    m.marshall(psms, writer); writer.write("\n");
//
//    writer.write(m.createProteinDetectionListStartTag("PDL_1", null) + "\n");
//    writer.write(m.createProteinDetectionListClosingTag()            + "\n");
//    writer.write(m.createAnalysisDataClosingTag()                    + "\n");
//    writer.write(m.createDataCollectionClosingTag()                  + "\n");
//    writer.write(m.createMzIdentMLClosingTag());
//
//    writer.close();
//  }
  public static CvList newCvList(Cv... cvs)
  {
    CvList cvList = new CvList();
    List<Cv> localCvList = cvList.getCv();

    if (Tools.isSet(cvs))
      for (Cv cv : cvs) localCvList.add(cv);

    return cvList;
  }
  public MzIDs addAnalysisSample(Sample s)
  {
    if (mSamples==null) mSamples = new AnalysisSampleCollection();

    if (s!=null)
    {
      List<Sample> samples = mSamples.getSample();
      samples.add(s);
    }
    return this;
  }

  /* Aim is to write out set up the analysisSoftwareList following this
   * structure: <AnalysisSoftware id="ID_software" name="xtandem"
   * version="2008.12.1.1" > <SoftwareName> <cvParam accession="MS:1001476"
   * name="xtandem" cvRef="PSI-MS" /> </SoftwareName>
   */
  public MzIDs addAnalysisSoftware(AnalysisSoftware... softwares)
  {
    if (mAnalysisSoftwareList==null) mAnalysisSoftwareList = new AnalysisSoftwareList();

    List<AnalysisSoftware>  analysisSoftwares = mAnalysisSoftwareList.getAnalysisSoftware();
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

    return this;
  }
  // abbriviated software entry without the CV!
  public static AnalysisSoftware newSoftware(String name)
  {
    AnalysisSoftware analysisSoftware = new AnalysisSoftware();
    analysisSoftware.setName(name);

//    Param tempParam = new Param();
//    tempParam.setParam(makeCvParam(engine.getPsiCV(), engine.getPsiName(), sPsiCV));
//    analysisSoftware.setSoftwareName(tempParam);

    analysisSoftware.setId(analysisSoftID);
//    analysisSoftware.setVersion(engine.getVersion());

    return analysisSoftware;
  }

  public static AnalysisSoftware newEngine(Engine engine)
  {
    if (engine==null) return null;

    AnalysisSoftware analysisSoftware = new AnalysisSoftware();
    analysisSoftware.setName(engine.getName());

    Param tempParam = new Param();
    tempParam.setParam(makeCvParam(engine.getPsiCV(), engine.getPsiName(), sPsiCV));
    analysisSoftware.setSoftwareName(tempParam);

    analysisSoftware.setId(analysisSoftID);
    analysisSoftware.setVersion(engine.getVersion());

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
  /* Setup Provider element as follows <Provider id="PROVIDER"> <ContactRole
   * Contact_ref="PERSON_DOC_OWNER"> <role> <cvParam accession="MS:1001271"
   * name="researcher" cvRef="PSI-MS"/> </role> </ContactRole> </Provider>
   */
  public MzIDs addProvider(Person docOwner)
  {
    if (mProvider==null) mProvider = new Provider();
    mProvider.setId("PROVIDER");

    ContactRole contactRole = new ContactRole();
    contactRole.setContact(docOwner);

    Role role = new Role();
    role.setCvParam(makeCvParam("MS:1001271", "researcher", sPsiCV));
    contactRole.setRole(role);

    mProvider.setContactRole(contactRole);

    return this;
  }

  public static Person newPerson(String first, String last, String address)
  {
    Person docOwner = new Person();
    docOwner.setId("PERSON_DOC_OWNER");
    if (first  !=null) docOwner.setFirstName(first);
    if (last   !=null) docOwner.setLastName(last);
    if (address!=null) docOwner.getCvParam().add(makeCvParam("MS:1000587", "contact address", sPsiCV, address));

    return docOwner;
  }
  /* TO DO Capture name and email of the user <AuditCollection> <Person
   * id="PERSON_DOC_OWNER" firstName="Andy" lastName="Jones"
   * email="someone@someuniversity.com"> <affiliations
   * Organization_ref="ORG_DOC_OWNER"/> </Person> <Organization
   * id="ORG_DOC_OWNER" address="Some address" name="Some place" />
   * </AuditCollection>
   */
  public MzIDs addAudit(Person docOwner, String affiliationName)
  {
    if (mAuditCollection==null) mAuditCollection=new AuditCollection();
    //List<Contact> contactList = auditCollection.getContactGroup();
    List<AbstractContact> contactList = mAuditCollection.getPersonOrOrganization();

    Organization org = new Organization();
    org.setId("ORG_DOC_OWNER");
    org.setName(affiliationName);
    org.getCvParam().add(makeCvParam("MS:1000586", "contact name", sPsiCV, get(docOwner.getCvParam(), "MS:1000587").getValue()));
    //org.setAddress(address);

    List<Affiliation> affList = docOwner.getAffiliation();
    Affiliation aff = new Affiliation();
    aff.setOrganization(org);
    affList.add(aff);
    contactList.add(docOwner);
    contactList.add(org);

    return this;
  }
  /**
   * <AnalysisProtocolCollection> <SpectrumIdentificationProtocol
   * id="SearchProtocol" AnalysisSoftware_ref="ID_software"> <SearchType>
   * <cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
   * </SearchType> <AdditionalSearchParams> <cvParam accession="MS:1001211"
   * name="parent mass type mono" cvRef="PSI-MS"/> <cvParam
   * accession="MS:1001256" name="fragment mass type mono" cvRef="PSI-MS"/>
   * </AdditionalSearchParams> <ModificationParams> <SearchModification
   * fixedMod="true"> <ModParam massDelta="57.021464" residues="C"> <cvParam
   * accession="UNIMOD:4" name="Carbamidomethyl" cvRef="UNIMOD" /> </ModParam>
   * </SearchModification> <SearchModification fixedMod="false"> <ModParam
   * massDelta="15.994919" residues="M"> <cvParam accession="UNIMOD:35"
   * name="Oxidation" cvRef="UNIMOD" /> </ModParam> </SearchModification>
   * </ModificationParams> <Enzymes independent="0"> <Enzyme id="ENZ_1"
   * CTermGain="OH" NTermGain="H" missedCleavages="1" semiSpecific="0">
   * <EnzymeName> <cvParam accession="MS:1001251" name="Trypsin"
   * cvRef="PSI-MS" /> </EnzymeName> </Enzyme> </Enzymes> <MassTable id="0"
   * msLevel="2"> </MassTable> <FragmentTolerance> <cvParam
   * accession="MS:1001412" name="search tolerance plus value" value="0.5"
   * cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton"
   * unitCvRef="UO" /> <cvParam accession="MS:1001413" name="search tolerance
   * minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221"
   * unitName="dalton" unitCvRef="UO" /> </FragmentTolerance>
   * <ParentTolerance> <cvParam accession="MS:1001412" name="search tolerance
   * plus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221"
   * unitName="dalton" unitCvRef="UO" /> <cvParam accession="MS:1001413"
   * name="search tolerance minus value" value="0.5" cvRef="PSI-MS"
   * unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
   * </ParentTolerance> <Threshold> <cvParam accession="MS:1001494" name="no
   * threshold" cvRef="PSI-MS" /> </Threshold>
   * </SpectrumIdentificationProtocol> </AnalysisProtocolCollection>
   *
   *
   */
  public MzIDs addAnalysisProtocol(SpectrumIdentificationProtocol p)
  {
    if (mAnalysisProtocolCollection==null) mAnalysisProtocolCollection = new AnalysisProtocolCollection();

    List<SpectrumIdentificationProtocol> sipList = mAnalysisProtocolCollection.getSpectrumIdentificationProtocol();

    sipList.add(p);
    return this;
  }
  /**
   * <AnalysisCollection> <SpectrumIdentification id="SI_1"
   * SpectrumIdentificationProtocol_ref="SearchProtocol"
   * SpectrumIdentificationList_ref="siiListID"
   * activityDate="2008-02-27T08:22:12"> <InputSpectra
   * SpectraData_ref="SD_1"/> <SearchDatabase
   * SearchDatabase_ref="search_database"/> </SpectrumIdentification>
   * </AnalysisCollection>
   *
   */
  public void newAnalysisCollection()
  {
    mAnalysis = new AnalysisCollection();
    List<SpectrumIdentification> specIdentList = getAnalysisCollection().getSpectrumIdentification();
    if (mSpecIdent==null) mSpecIdent = new SpectrumIdentification();
    mSpecIdent.setId(specIdentID);
    mSpecIdent.setSpectrumIdentificationProtocol(getSiProtocol());
    specIdentList.add(mSpecIdent);
    List<SearchDatabaseRef> searchDBRefList = mSpecIdent.getSearchDatabaseRef();
    SearchDatabaseRef searchDBRef = new SearchDatabaseRef();
    searchDBRef.setSearchDatabase(getSearchDB());
    searchDBRefList.add(searchDBRef);

    if (getParam(CAT_SEARCH, SPEC_SEARCHED)!= null)
    {
      List<InputSpectra> inputSpecList = mSpecIdent.getInputSpectra();
      InputSpectra inputSpec = new InputSpectra();
      inputSpec.setSpectraData(getSpectraData());
      inputSpecList.add(inputSpec);
      //specIdentList.add(specIdent);     //bug caused by this line - removed by ARJ 27/03/2013
    }
  }

  public SearchModification newSearchModification(UnimodMod m, boolean fixed, ModAttachment attachment)
  {
    SearchModification searchMod = new SearchModification();
    searchMod.setFixedMod(fixed);
    List<CvParam> modCvParamList = searchMod.getCvParam();
    modCvParamList.add(newCvParam(m.getLabel(), m.getFullName()));
    searchMod.setMassDelta((float )m.getMass().getMolecularMass());
    List<String> residueList = searchMod.getResidues();
    residueList.addAll(m.getSites());

    if (ModAttachment.N_TERM.equals(attachment))
    {
      modCvParamList.add(makeCvParam("MS:1001189", "modification specificity key N-term", sPsiCV));
    }
    else if (ModAttachment.C_TERM.equals(attachment))
    {
      modCvParamList.add(makeCvParam("MS:1001190", "modification specificity key C-term", sPsiCV));
    }
    // no provision for protein N/C-terminus?

    return searchMod;
  }
  public static Table<String, String, String> readParams(String params)
  {
    TabFile param = null;
    Table<String, String, String> cat_key_val = HashBasedTable.create();
    try
    {
      try
      {
        param = new TabFile(params, TabFile.comma, "//");
        while (param.hasNext())
        {
          if  (!CAT_MOD.equals(param.get("Category"))) cat_key_val.put(param.get("Category"), param.get("Key"), param.get("Value"));
          else cat_key_val.put(param.get("Category"), param.get("Key")+"#"+param.get("Required"), param.get("Value"));
        }
      }
      finally
      {
        if (param!=null) param.close();
      }
    }
    catch (IOException e) {}

    return cat_key_val;
  }
//  MOD, N, oxidation of M,                     UNIMOD:35#Oxidation#Var#Any#M#15.9949
//  MOD, Y, itraq n-term,                       UNIMOD:529#iTRAQ4plex114#Fixed#Nt#-#144.106
//  MOD, Y, itraq k,                            UNIMOD:529#iTRAQ4plex114#Fixed#Any#K#144.106
//  MOD, N, pyro-glu from n-term Q,             UNIMOD:28#Gln-&gt;pyro-Glu#Var#Nt#Q#-17.027
//  MOD, N, acetyl protein n term,              UNIMOD:1#Acetyl#Var#NT#-#42.01
//  MOD, Y, carbamidomethyl C,                  UNIMOD:4#Carbamidomethyl#Fixed#Any#C#57.021464
  public SearchModification newSearchModification(String name, String mod)
  {
    String[] strs = Strs.split(mod, '#', true), names=Strs.split(name, '#',true);

    SearchModification searchMod = new SearchModification();
    searchMod.setFixedMod("Fixed".equals(strs[2]));

    List<CvParam> modCvParamList = searchMod.getCvParam();
    modCvParamList.add(newCvParam(strs[0], strs[1], names[0]));
    searchMod.setMassDelta(Stats.toFloat(strs[5]));

    List<String> residueList = searchMod.getResidues();
    if (Strs.isSet(strs[4]))
      for (int i=0; i<strs[4].length(); i++) residueList.add(strs[4].substring(i, i + 1));

    if      ("Nt".equals(strs[3])) modCvParamList.add(newCvParam("MS:1001189", "modification specificity key N-term"));
    else if ("Ct".equals(strs[3])) modCvParamList.add(newCvParam("MS:1001190", "modification specificity key C-term"));
    else if ("NT".equals(strs[3])) modCvParamList.add(newCvParam("MS:1002057", "modification specificity protein N-term"));
    else if ("CT".equals(strs[3])) modCvParamList.add(newCvParam("MS:1002058", "modification specificity protein C-term"));

    return searchMod;
  }
  public SpectrumIdentificationProtocol newSpectrumIdentificationProtocol()
  {
    SpectrumIdentificationProtocol siProtocol = new SpectrumIdentificationProtocol();

    siProtocol.setId(siProtocolID);
    siProtocol.setAnalysisSoftware(newEngine(Engine.valueOf(getParam(CAT_SEARCH, NAME_SOFTWARE))));

    //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
    siProtocol.setSearchType(newParam("MS:1001083", "ms-ms search"));

    //List<CvParam> cvParamList = siProtocol.getAdditionalSearchCvParams();
    ParamList paramList = siProtocol.getAdditionalSearchParams();
    if (paramList==null)
    {
      paramList = new ParamList();
      siProtocol.setAdditionalSearchParams(paramList);
    }
    List<CvParam> cvParamList = paramList.getCvParam();
    cvParamList.add(newCvParam(TYPE_PARENT_MASS, getParam(CAT_SEARCH, TYPE_PARENT_MASS), '#'));
    cvParamList.add(newCvParam(TYPE_FRAG_MASS,   getParam(CAT_SEARCH, TYPE_FRAG_MASS),   '#'));

    ModificationParams           modParams = new ModificationParams();
    List<SearchModification> searchModList = modParams.getSearchModification();

    for (String mod : getParams().row(CAT_MOD).keySet())
    {
      searchModList.add(newSearchModification(mod, getParam(CAT_MOD, mod)));
    }

    siProtocol.setModificationParams(modParams);

    Enzymes enzymes = siProtocol.getEnzymes();

    if (enzymes == null) {
      enzymes = new Enzymes();
      siProtocol.setEnzymes(enzymes);
    }
    enzymes.setIndependent(false);

    List<Enzyme> enzymeList = enzymes.getEnzyme();

    Enzyme enzyme = new Enzyme();
    //[KR]|{P}

    enzyme.setId("Enz1");
    enzyme.setCTermGain("OH");
    enzyme.setNTermGain("H");
    enzyme.setMissedCleavages(Integer.parseInt(getParam(CAT_SEARCH, MISSED_CLEAVAGE)));
    enzyme.setSemiSpecific(false);
    ParamList eParamList = enzyme.getEnzymeName();
    if (eParamList == null) {
      eParamList = new ParamList();
      enzyme.setEnzymeName(eParamList);
    }
    List<CvParam> eCvParamList = eParamList.getCvParam();
    eCvParamList.add(newCvParam(ENZYME, getParam(CAT_SEARCH, ENZYME), '#'));
    enzymeList.add(enzyme);

    Tolerance fragTol = new Tolerance();
    Tolerance parTol  = new Tolerance();

    List<CvParam> fragCvList = fragTol.getCvParam();
    fragCvList.add(newCvParamMass(TOL_FRAG_MINUS,  getParam(CAT_SEARCH, TOL_FRAG_MINUS),   '#'));
    fragCvList.add(newCvParamMass(TOL_FRAG_PLUS,   getParam(CAT_SEARCH, TOL_FRAG_PLUS),    '#'));

    List<CvParam> parCvList = parTol.getCvParam();
    parCvList.add(newCvParamMass(TOL_PARENT_MINUS, getParam(CAT_SEARCH, TOL_PARENT_MINUS), '#'));
    parCvList.add(newCvParamMass(TOL_PARENT_PLUS,  getParam(CAT_SEARCH, TOL_PARENT_PLUS),  '#'));

    siProtocol.setFragmentTolerance(fragTol);
    siProtocol.setParentTolerance(parTol);

    // siProtocol.getThresholdCvParams();
    ParamList sip_paramList = siProtocol.getThreshold();
    if (sip_paramList == null) {
      sip_paramList = new ParamList();
      siProtocol.setThreshold(sip_paramList);
    }
    cvParamList = sip_paramList.getCvParam();
    cvParamList.add(newCvParam(THLD_PSM, getParam(CAT_SEARCH, THLD_PSM)));

    return siProtocol;
  }
  // SEARCH, Y, Fragment search tolerance plus,  MS:1001412#search tolerance plus value#Da
  // SEARCH, Y, Fragment search tolerance minus, MS:1001413#search tolerance minus value#DA
  public static CvParam newCvParamMass(String name, String val, char t)
  {
    CvParam cvParam = new CvParam();

    String[] strs = Strs.split(val, t, true);

    //<cvParam accession="MS:1001413" name="search tolerance minus value" value="0.5" cvRef="PSI-MS" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO" />
    cvParam.setCv(sPsiCV);
    cvParam.setUnitCv(sUnitCV);

    cvParam.setAccession(strs[0]);
    cvParam.setName(strs[1]);
    cvParam.setValue(name);

    if      (strs.length<3 || "Da".equals(strs[2]))
    {
      cvParam.setUnitAccession("UO:0000221");
      cvParam.setUnitName("dalton");
    }
    else if (strs.length>2 && "ppm".equals(strs[2]))
    {
      cvParam.setUnitAccession("UO:0000169");
      cvParam.setUnitName("parts per million");
    }
    return cvParam;
  }
  public MzIDs addSourceFile(String filename, String format)
  {
    List<SourceFile> sourceFileList = mInputs.getSourceFile();
    SourceFile sourceFile = new SourceFile();
    sourceFile.setLocation(filename);
    sourceFile.setId(sourceFileID);
    FileFormat ff = new FileFormat();
    ff.setCvParam(newCvParam(format, format));

    sourceFile.setFileFormat(ff);
    sourceFileList.add(sourceFile);

    return this;
  }
  public Inputs newInputs() {

    mInputs = new Inputs();
    List<SearchDatabase> searchDBList = getInputs().getSearchDatabase();

    if (getParam(CAT_SEARCH, NAME_FIRST)!=null)
      setOwner(getParam(CAT_SEARCH, NAME_FIRST), getParam(CAT_SEARCH, NAME_LAST), getParam(CAT_SEARCH, ADDRESS));

    mSearchDB = new SearchDatabase();
    getSearchDB().setId(searchDBID);

    UserParam param = new UserParam();
    //param.setName(settings.MSSearchSettings_db);
    Param tempParam = new Param();

    param.setName(getParam(CAT_SEARCH, NAME_DB));
    tempParam.setParam(param);
    getSearchDB().setDatabaseName(tempParam);
    getSearchDB().setLocation(getParam(CAT_SEARCH, PATH_DB_LOCAL));

    List<CvParam> searchDBCvParamList = getSearchDB().getCvParam();

    addCvParam(searchDBCvParamList, DECOY_COMP, getParam(CAT_SEARCH, DECOY_COMP),  '#');
    addCvParam(searchDBCvParamList, DECOY_REGEX,getParam(CAT_SEARCH, DECOY_REGEX), '#');
    addCvParam(searchDBCvParamList, DECOY_TYPE, getParam(CAT_SEARCH, DECOY_TYPE),  '#');

    FileFormat ff = new FileFormat();
    ff.setCvParam(newCvParam(FMT_FILE_DB, getParam(CAT_SEARCH, FMT_FILE_DB), '#'));
    getSearchDB().setFileFormat(ff);
    searchDBList.add(getSearchDB());

    if (getParam(CAT_SEARCH, SPEC_SEARCHED) != null)
    {
//      List<SpectraData> spectraDataList = inputs.getSpectraData();
      mSpectraData = newSpectraData(getParam(CAT_SEARCH, SPEC_SEARCHED));
//      SpectrumIDFormat sif = new SpectrumIDFormat();
//      sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", sPsiCV));
//      getSpectraData().setSpectrumIDFormat(sif);
//
//      ff = new FileFormat();
//      ff.setCvParam(newCvParam(FMT_DATA_SPEC, params.get(CAT_SEARCH, FMT_DATA_SPEC), '#'));
//      getSpectraData().setFileFormat(ff);
//
//      getSpectraData().setId(mSpectraDataID);
//      getSpectraData().setLocation(params.get(CAT_SEARCH, SPEC_SEARCHED));
//      spectraDataList.add(getSpectraData());
    }
    if (mScoreCV==null) mScoreCV = new HashMap<>();
    for (String score : getParams().row(CAT_SCORE).keySet())
    {
      mScoreCV.put(score, newCvParam(score, getParam(CAT_SCORE, score), '#'));
    }

    return getInputs();
  }
  public SpectraData newSpectraData(String run)
  {
    List<SpectraData> spectraDataList = getInputs().getSpectraData();
    SpectraData spectraData = new SpectraData();
    SpectrumIDFormat sif = new SpectrumIDFormat();
    sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", sPsiCV));
    spectraData.setSpectrumIDFormat(sif);

    FileFormat ff = new FileFormat();
    ff.setCvParam(newCvParam(FMT_DATA_SPEC, getParam(CAT_SEARCH, FMT_DATA_SPEC), '#'));
    spectraData.setFileFormat(ff);

    mSpectraDataID = "SD_" + mSpecDataCounter;
    mSpecDataCounter++;
    spectraData.setId(mSpectraDataID);
    spectraData.setLocation(run);
    spectraDataList.add(spectraData);

    if (mSpecIdent==null) mSpecIdent = new SpectrumIdentification();

    List<InputSpectra> inputSpecList = mSpecIdent.getInputSpectra();
    InputSpectra inputSpec = new InputSpectra();
    inputSpec.setSpectraData(spectraData);
    inputSpecList.add(inputSpec);
    //specIdentList.add(specIdent);
    mRun2SpectraData.put(run, spectraData);

    return spectraData;
  }
  // convert a PeptideMatch to PeptideEvdence for MzID output
  public static PeptideEvidence newPeptideEvidence(PeptideMatch pm, SearchDatabase db)
  {
    PeptideProteinMatch ppm = Tools.front(pm.getProteinMatches());

    DBSequence dbSequence = new DBSequence();
    dbSequence.setId("dbseq_" + ppm.getAccession());
    dbSequence.setAccession(ppm.getAccession());
//    dbSequenceList.add(dbSequence);
//    dbSequence.getCvParam().add(makeCvParam("MS:1001088", "protein description", sPsiCV, ppm.getAccession()));
    dbSequence.setSearchDatabase(db);

    String pepID = pm.toSymbolString() + "_" + PSMs.toNumModSequence(pm);

    Peptide pep = new Peptide();
    pep.setPeptideSequence(pm.toSymbolString());
    pep.setId(pepID);

    PeptideEvidence peptideEvidence = new PeptideEvidence();
    peptideEvidence.setId(pm.toSymbolString()+"_"+ppm.getAccession()+"_"+ppm.getStart()+"_"+ppm.getEnd());
    peptideEvidence.setDBSequence(dbSequence);
    peptideEvidence.setPeptide(pep);
    peptideEvidence.setStart(ppm.getStart());
    peptideEvidence.setEnd(ppm.getEnd());
    peptideEvidence.setIsDecoy(PeptideProteinMatch.HitType.DECOY.equals(ppm.getHitType()));

    if (ppm.getPreviousAA().isPresent()) peptideEvidence.setPre(ppm.getPreviousAA().get());
    if (ppm.getNextAA(    ).isPresent()) peptideEvidence.setPre(ppm.getNextAA(    ).get());

    List<Modification> modList = pep.getModification();

    for (ModificationMatch mm : pm.getModifications(ModAttachment.all))
    {
      Modification mzidMod = new Modification();
      mzidMod.setMonoisotopicMassDelta(mm.getMassShift());
      mzidMod.setLocation(mm.getPosition());

      if (ModAttachment.N_TERM.equals(mm.getModAttachment()))
      {
        mzidMod.getCvParam().add(makeCvParam("MS:1001189", "modification specificity key N-term", sPsiCV));
      }
      else if (ModAttachment.C_TERM.equals(mm.getModAttachment()))
      {
        mzidMod.getCvParam().add(makeCvParam("MS:1001190", "modification specificity key C-term", sPsiCV));
      }
      else
      {
        mzidMod.getResidues().add(mm.getResidue().getSymbol());
      }
      for (int i=0; i<mm.getCandidateCount(); i++)
      {
        org.expasy.mzjava.proteomics.mol.modification.Modification m = mm.getModificationCandidate(i);

        // TODO we need to get the UNIMOD accession instead of using the label
        mzidMod.getCvParam().add(newCvParam(m.getLabel(), m.getLabel()));
        modList.add(mzidMod);
      }
    }
    return peptideEvidence;
  }
  // Spectrum number, Filename/id,                                          Peptide,      E-value,          Mass,     gi, Accession,        Start,  Stop, Defline,                                                                                                        Mods,                                         Charge, Theo Mass,  P-value,            NIST score
  // 11,              Mouse_Plasma_LIRKO2_01_27Aug12_Lynx_12-06-05.14.14.,  EAAILYKPAVSTK,783.759815794219, 1822.118, 0,  BL_ORD_ID:8530,   25444,  25456,gi|568916279|ref|XP_006499221.1|Ttn$ PREDICTED: titin isoform X2 [Mus musculus],                             ,                                                3,      1822.105,   0.773701693775142,  0
  //160,              Mouse_Plasma_LIRKO2_01_27Aug12_Lynx_12-06-05.185.185.,NDTSQTmHSnM,  22.4967193416469, 1425.578, 0,  BL_ORD_ID:102317, 308,    318,  XXX_gi|568954962|ref|XP_006509551.1|Cdkn2aip$ PREDICTED: CDKN2A-interacting protein isoform X1 [Mus musculus]," oxidation of M:7 ,deamidation of N and Q:10", 3,      1425.567,   1.24981774120261,   0
  public void buildPeptideMatches(Multimap<SpectrumIdentifier, PeptideMatch> spec_matches)
  {
    mSequences = new SequenceCollection();
    List<Peptide> peptideList = getSequenceCollection().getPeptide();

    mSiList = new SpectrumIdentificationList();
    getSpecIDs().setId(siiListID);

    SpectrumIdentification specIdent = new SpectrumIdentification();
    specIdent.setSpectrumIdentificationList(getSpecIDs());
    List<DBSequence>           dbSequenceList = getSequenceCollection().getDBSequence();
    List<PeptideEvidence> peptideEvidenceList = getSequenceCollection().getPeptideEvidence();

    AnalysisData analysisData = new AnalysisData();

    List<SpectrumIdentificationResult> sirList = getSpecIDs().getSpectrumIdentificationResult();

    int sirCounter = 1;

    // go through the key matches
    for (SpectrumIdentifier ms2 : spec_matches.keySet())
    {
      SpectrumIdentificationResult sir = new SpectrumIdentificationResult();
      sirList.add(sir);
      sir.setId("SIR_" + sirCounter);
      sir.setSpectrumID(ms2.getSpectrum());
      sirCounter++;

      // get a run specific spectraData object
      mSpectraData = mRun2SpectraData.get(ms2.getSpectrumFile().get());

      //Create a new spectra data object and add input spectra to the SpectrumIdentification object
      if (mSpectraData == null) {
        mSpectraData = newSpectraData(ms2.getSpectrumFile().get());
//        List<SpectraData> spectraDataList = inputs.getSpectraData();
//        spectraData = new SpectraData();
//        SpectrumIDFormat sif = new SpectrumIDFormat();
//        sif.setCvParam(makeCvParam("MS:1000774", "multiple peak list nativeID format", psiCV));
//        spectraData.setSpectrumIDFormat(sif);
//
//        FileFormat ff = new FileFormat();
//        ff.setCvParam(makeCvParam(paramToCvParamAcc.get("Spectra data file format"), paramToCvParamName.get("Spectra data file format"), psiCV));
//        spectraData.setFileFormat(ff);
//
//        spectraDataID = "SD_" + specDataCounter;
//        specDataCounter++;
//        spectraData.setId(spectraDataID);
//        spectraData.setLocation(spectraLocation);
//        spectraDataList.add(spectraData);
//        List<InputSpectra> inputSpecList = specIdent.getInputSpectra();
//        InputSpectra inputSpec = new InputSpectra();
//        inputSpec.setSpectraData(spectraData);
//        inputSpecList.add(inputSpec);
//        //specIdentList.add(specIdent);
//        locationToSpectraDataMap.put(spectraLocation, spectraData);
      }

      sir.setSpectraData(getSpectraData());
      // setup the rest of the spectrum header
      List<CvParam> sir_cvParamList = sir.getCvParam();
      if (ms2.getSpectrumFile().isPresent())
        sir_cvParamList.add(newCvParam("MS:1000796", "spectrum title", ms2.getSpectrumFile().get()));
      if (ms2.getRetentionTimes()!=null && ms2.getRetentionTimes().size()>0)
        sir_cvParamList.add(newCvParam("MS:1000016", "scan start time", ms2.getRetentionTimes().getFirst().getTime()+""));
      if (ms2.getScanNumbers()!=null && ms2.getScanNumbers().size()>0)
        sir_cvParamList.add(newCvParam("MS:1001115", "scan number(s)", ms2.getScanNumbers().getFirst().getValue()+""));

      List<SpectrumIdentificationItem> siiList = sir.getSpectrumIdentificationItem();

      for (PeptideMatch pm : spec_matches.get(ms2))
      {
        PeptideEvidence pv = newPeptideEvidence(pm, Tools.front(getInputs().getSearchDatabase()));
        dbSequenceList.add(pv.getDBSequence());
        peptideEvidenceList.add(pv);

        peptideList.add(pv.getPeptide());

        SpectrumIdentificationItem sii = new SpectrumIdentificationItem();
        if (ms2.getAssumedCharge(       ).isPresent()) sii.setChargeState(             ms2.getAssumedCharge().get());
        if (ms2.getPrecursorMz(         ).isPresent()) sii.setExperimentalMassToCharge(ms2.getPrecursorMz().get());
        if (ms2.getPrecursorNeutralMass().isPresent()) sii.setCalculatedMassToCharge(  ms2.getPrecursorNeutralMass().get());

        // deposit the scores
        List<CvParam> cvParamList = sii.getCvParam();
        for (String score : pm.getScoreMap().keySet())
        {
          CvParam scoreCV = getScoreCV(score);
          // MyriMatch:MVH
          if (scoreCV!=null)
          {
            scoreCV.setValue(pm.getScore(score)+"");
            cvParamList.add(scoreCV);
          }
          else cvParamList.add(newCvParam(score, Tools.back(Strs.split(score, ':')), pm.getScore(score)+""));
        }
        sii.setPeptide(pv.getPeptide());
        // set to PASS by default
        sii.setPassThreshold(true);

        PeptideEvidenceRef pepEvidRef = new PeptideEvidenceRef();
        pepEvidRef.setPeptideEvidence(pv);
        sii.getPeptideEvidenceRef().add(pepEvidRef);

        addSIIToListAndSetRank(siiList, sii, getEngine().getCanonicalScore(), true, sir.getId());
      }
    }
  }

  /* Helper method to retrive a particular score type from an SII, based on the CV accession. Will return 0 if the accession is not found. */
  public static double getScoreFromSII(SpectrumIdentificationItem sii, String cvParamAccForScore)
  {
    double score = 0.0;
    for (CvParam cvParam : sii.getCvParam())
      if (cvParam.getAccession().equals(cvParamAccForScore)) score = Double.parseDouble(cvParam.getValue());

    return score;
  }

  /* Accepts the SIIs associated with any SIR and inserts the correct rank values */
  private void addSIIToListAndSetRank(List<SpectrumIdentificationItem> siiList, SpectrumIdentificationItem sii,
                                      String scoreCvParamNameToOrderBy, boolean orderLowToHigh, String sirID)
  {
    double scoreOfNewSII = getScoreFromSII(sii, scoreCvParamNameToOrderBy);
    final String cvParamScore = scoreCvParamNameToOrderBy;
    final boolean lowToHigh = orderLowToHigh;

    siiList.add(sii);
    Collections.sort(siiList, new Comparator<SpectrumIdentificationItem>() {

      @Override
      public int compare(SpectrumIdentificationItem sii1, SpectrumIdentificationItem sii2)
      {
        double sii1Score = getScoreFromSII(sii1, cvParamScore);
        double sii2Score = getScoreFromSII(sii2, cvParamScore);
        int i = 0;

        if (lowToHigh)
        {
          if (sii1Score<sii2Score) { i=-1; } else if (sii1Score>sii2Score) { i=+1; } else {i = 0;}
        }
        else
        {
          if (sii2Score<sii1Score) { i=-1; } else if (sii2Score>sii1Score) { i=+1; } else { i = 0; }
        }
        return i;
      }
    });

    int rank = 0;
    int innerRankCounter = 2;
    double lastScore = -999.0;

    for (SpectrumIdentificationItem spii : siiList)
    {
      double score = getScoreFromSII(spii, scoreCvParamNameToOrderBy);
      boolean sameScore = false;
      if (score != lastScore)
      { //Otherwise set the same rank as previous
        rank++;
        sameScore = true;
        spii.setId(sirID + "_SII_" + rank);
        spii.setRank(rank);
        innerRankCounter = 2;
      }
      else
      {
        spii.setId(sirID + "SII_" + rank + "_" + innerRankCounter);
        spii.setRank(rank);
        innerRankCounter++;
      }
    }
  }
  public static List<CvParam> addCvParam(List<CvParam> cvs, String name, String val, char t)
  {
    if (cvs!=null && val!=null) cvs.add(newCvParam(name, val, t));
    return cvs;
  }
  public static Cv newCV(String url, String id, String name, String version)
  {
    Cv unitCV = newCV(url,id,name);
    unitCV.setVersion(version);
    return unitCV;
  }
  public static Cv newCV(String url, String id, String name)
  {
    Cv unitCV = new Cv();
    unitCV.setUri(url);
    unitCV.setId(id);
    unitCV.setFullName(name);

    return unitCV;
  }
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
  public static CvParam newCvParam(String name, String vals, char t)
  {
    String[] strs = Strs.split(vals, t, true);
    return newCvParam(strs[0], strs[1], name);
  }

  public static CvParam newCvParam(String accession, String name)
  {
    return makeCvParam(accession, name, sPsiCV);
  }
  public static CvParam newCvParam(String accession, String name, String value)
  {
    return makeCvParam(accession, name, sPsiCV, value);
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
  public static CvParam get(List<CvParam> params, String acc)
  {
    if (Tools.isSet(params))
      for (CvParam param : params)
        if (Strs.equals(acc, param.getAccession())) return param;

    return null;
  }
  public static Param newParam(String acc, String name)
  {
    //<cvParam accession="MS:1001083" name="ms-ms search" cvRef="PSI-MS"/>
    Param tempParam = new Param();
    tempParam.setParam(makeCvParam(acc, name, sPsiCV));

    return tempParam;
  }
}
