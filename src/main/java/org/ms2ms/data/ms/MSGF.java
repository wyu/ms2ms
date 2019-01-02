package org.ms2ms.data.ms;

import com.google.common.collect.*;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by yuw on 5/18/17.
 */
public class MSGF
{
  public static Table<String, String, Map<String, String>> parseMSGF(String search) throws IOException
  {
    Table<String, String, Map<String, String>> out = HashBasedTable.create();

    TabFile  msgf = new TabFile(search, TabFile.tabb);
    while(msgf.hasNext())
    {
      String  row = "#"+msgf.get("ScanNum"),
          peptide = msgf.get("Peptide").replaceAll("[0|1|2|3|4|5|6|7|8|9|\\\\+|\\\\.|\\\\-]", "");

      Map<String, String> props = out.get(row, peptide);
      if (props==null || Tools.getDouble(props, "MSGFScore")<msgf.getDouble("MSGFScore"))
        out.put(row, peptide, new HashMap<>(msgf.getMappedRow()));
    }
    msgf.close();

    return out;
  }
  public static void parseMSGFSearches(String search, String assignments) throws IOException
  {
    Map<String, String> runscan_seq = readAssignments(assignments);

//    Map<String, String> runscan_seq = new HashMap<>();
//
//    if (Strs.isSet(assignments))
//    {
//      TabFile refs = new TabFile(assignments, TabFile.comma);
//      while (refs.hasNext())
//      {
//        if (refs.get("Protein")==null || refs.get("Protein").indexOf("XXX")<0)
//          runscan_seq.put(refs.get("Run")+"#"+refs.getInt("Scan"), refs.get("Sequence")+" m/z"+refs.get("mz") + "@z"+refs.get("z"));
//      }
//      refs.close();
//    }

    Dataframe out = readMSGF(null,search, 0.01, runscan_seq, null);

    out.init(false).write(new FileWriter("/Users/yuw/Apps/pipeline/Chorus/data/cKit/NR/NR_20170513_qual.tsv"), "\t", true, "NA");

//    TreeMultimap<Integer, Double> all = TreeMultimap.create(), annot = TreeMultimap.create();
//    Multimap<String, String> distincts = HashMultimap.create();
//
//    Table<String, String, Map<String, String>> psm = parseMSGF(search);
//
//    int passed=0, totals=0, ids=0, specs=0, good_miss=0;
//    Dataframe out = new Dataframe("");
//    for (Map<String, String> props : psm.values())
//    {
//      totals++;
//
//      double     fdr = Tools.getDouble(props, "QValue");
//      int        row = Tools.getInt(props, "ScanNum");
//      String peptide = props.get("Peptide");
//      String[] scans = Strs.split(props.get("Title"), '+');
//
//      specs += scans.length;
//
//      Set<String> scns = new HashSet<>();
//      for (String scan : scans)
//      {
//        scns.add(Strs.split(scan, "m", true).get(0));
//      }
//      all.put(scans.length, fdr);
//      if (fdr<=0.01)
//      {
//        passed++; ids += scans.length;
//
//        annot.put(scans.length, fdr);
//        distincts.put(peptide.replaceAll("[0|1|2|3|4|5|6|7|8|9|\\\\+|\\\\.|\\\\-]", ""), row+"*"+scans.length+"*"+scans[0]);
//
//        out.put(row, "Charge",       Tools.getInt(props, "Charge"));
//        out.put(row, "DeNovoScore",  Tools.getDouble(props,"DeNovoScore"));
//        out.put(row, "IsotopeError", Tools.getInt(props,"IsotopeError"));
//        out.put(row, "MSGFScore",    Tools.getDouble(props,"MSGFScore"));
//        out.put(row, "PepQValue",    Tools.getDouble(props,"PepQValue"));
//        out.put(row, "m.z",          Tools.getDouble(props,"Precursor"));
//        out.put(row, "ppm",          Tools.getDouble(props,"PrecursorError(ppm)"));
//        out.put(row, "QValue",       Tools.getDouble(props,"QValue"));
//        out.put(row, "SpecEValue",   Tools.getDouble(props,"SpecEValue"));
//        out.put(row, "Peptide",      props.get("Peptide"));
//        out.put(row, "Protein",      props.get("Protein"));
//      }
//      else if (Sets.intersection(scns, runscan_seq.keySet()).size()>0)
//      {
//        good_miss+=scans.length;
//        if (scans.length>2)
//        {
//          System.out.println(row+", "+peptide+"@"+props.get("QValue")+" m/z"+Tools.d2s(Tools.getDouble(props, "Precursor"), 4)+" @z"+Tools.getInt(props, "Charge"));
//          for (String scn : scns)
//            System.out.println("    "+scn+"\t"+runscan_seq.get(scn));
//        }
//      }
//    }
//    out.init(false);
//    out.write(new FileWriter("/Users/yuw/Apps/pipeline/Chorus/data/cKit/NR/NR_20170513_qual.tsv"), "\t");
//
//    int annotated=0, alls=0;
//    for (Integer cnt : annot.keySet()) annotated+=cnt*annot.get(cnt).size();
//    for (Integer cnt :   all.keySet())      alls+=cnt*  all.get(cnt).size();
//
//    System.out.println("Annot/All: " + annotated+"/"+alls + "--> " + Tools.d2s(100d*annotated/alls, 2));
  }
  public static Map<String, String> readAssignments(String assignments) throws IOException
  {
    System.out.println("Reading the assignments from " + assignments);

    Map<String, String> runscan_seq=new HashMap<>();

    if (Strs.isSet(assignments)) {
      TabFile refs=new TabFile(assignments, TabFile.comma);
      while (refs.hasNext()) {
        if (refs.get("Protein")==null||refs.get("Protein").indexOf("XXX")<0)
          runscan_seq.put(refs.get("Run")+"#"+refs.getInt("Scan"), refs.get("Sequence")+" m/z"+refs.get("mz")+"@z"+refs.get("z"));
      }
      refs.close();
    }
    return runscan_seq;
  }
  public static StringBuffer countMSGFSearches(String search, Map<String, String> runscan_seq, String name) throws IOException
  {
    System.out.print("Counting the matches: "+search);

    TreeMultimap<Integer, Double> all = TreeMultimap.create(), annot = TreeMultimap.create();
    Multimap<String, String> distincts = HashMultimap.create();

    int passed=0, ids=0, specs=0, good_miss=0, counts=0;

    TabFile  msgf = new TabFile(search, TabFile.tabb);
    Map<Integer, Double> scores = new HashMap<>();
    while(msgf.hasNext())
    {
      int        row = msgf.getInt("ScanNum");
      String peptide = msgf.get("Peptide").replaceAll("[0|1|2|3|4|5|6|7|8|9|\\\\+|\\\\.|\\\\-]", "");

      // total number of MSGF hits
      if (++counts%10000==0) System.out.print(".");

      if (scores.get(row)==null || scores.get(row)<msgf.getDouble("MSGFScore"))
      {
        double     fdr = msgf.getDouble("QValue");
        String[] scans = Strs.split(msgf.get("Title"), '+');

        // only count it once
        if (!scores.containsKey(row)) specs += scans.length;

        Set<String> scns = new HashSet<>();
        for (String scan : scans)
        {
          scns.add(Strs.split(scan, "m", true).get(0));
        }
        all.put(scans.length, fdr);
        if (fdr<=0.01)
        {
          passed++; ids += scans.length;

          annot.put(scans.length, fdr);
          distincts.put(peptide, row+"*"+scans.length+"*"+scans[0]);
        }
        else if (Sets.intersection(scns, runscan_seq.keySet()).size()>0)
        {
          good_miss+=scans.length;
        }
        scores.put(row, msgf.getDouble("MSGFScore"));
      }
    }
    msgf.close();
    System.out.println();

    int annotated=0, alls=0;
    for (Integer cnt : annot.keySet()) annotated+=cnt*annot.get(cnt).size();
    for (Integer cnt :   all.keySet())      alls+=cnt*  all.get(cnt).size();

    return new StringBuffer("|"+name+"|"+passed+"|"+ids+"|"+good_miss+"|"+Tools.d2s(100d*ids/specs, 2)+"|"+specs+"|"+counts+"|\n");
  }

  public static Dataframe readMSGFs(Dataframe out, String root, String search, double min_fdr, Map<String, String> refs, String tag) throws IOException
  {
    List<String> searches = IOs.listFiles(root, new WildcardFileFilter(search),1);

    if (Tools.isSet(searches))
      for (String s : searches)
        out = readMSGF(out, s, min_fdr, refs, tag);

    return out;
  }
  public static Dataframe readMSGF(Dataframe out, String search, double min_fdr, Map<String, String> refs, String tag) throws IOException
  {
    Table<String, String, Map<String, String>> psm = parseMSGF(search);

    TreeMultimap<Integer, Double> all = TreeMultimap.create(), annot = TreeMultimap.create();
    Multimap<String, String> distincts = HashMultimap.create();

    if (out==null) out = new Dataframe("");

    int passed=0, totals=0, ids=0, specs=0, good_miss=0;
    for (Map<String, String> props : psm.values())
    {
      totals++;

      double     fdr = Tools.getDouble(props, "QValue");
      String     row = search+"#"+props.get("ScanNum")+"$"+props.hashCode();
      String peptide = props.get("Peptide");
      String[] scans = Strs.split(props.get("Title"), '+');

      Set<String> scns = new HashSet<>();
      if (Tools.isSet(scans))
      {
        specs += scans.length;

        for (String scan : scans)
        {
          scns.add(Strs.split(scan, "m", true).get(0));
        }
        all.put(scans.length, fdr);
      }
      if (fdr<=min_fdr)
      {
        passed++;

        if (Tools.isSet(scans))
        {
          ids += scans.length;
          annot.put(scans.length, fdr);
          distincts.put(peptide.replaceAll("[0|1|2|3|4|5|6|7|8|9|\\\\+|\\\\.|\\\\-]", ""), row+"*"+scans.length+"*"+scans[0]);
        }

        out.put(row, "Charge",       Tools.getInt(props, "Charge"));
        out.put(row, "DeNovoScore",  Tools.getDouble(props,"DeNovoScore"));
        out.put(row, "IsotopeError", Tools.getInt(props,"IsotopeError"));
        out.put(row, "MSGFScore",    Tools.getDouble(props,"MSGFScore"));
        out.put(row, "PepQValue",    Tools.getDouble(props,"PepQValue"));
        out.put(row, "m.z",          Tools.getDouble(props,"Precursor"));
        out.put(row, "ppm",          Tools.getDouble(props,"PrecursorError(ppm)"));
        out.put(row, "QValue",       Tools.getDouble(props,"QValue"));
        out.put(row, "SpecEValue",   Tools.getDouble(props,"SpecEValue"));
        out.put(row, "Peptide",      props.get("Peptide"));
        out.put(row, "Sequence",     props.get("Peptide").replaceAll("[0|1|2|3|4|5|6|7|8|9|\\\\+|\\\\.|\\\\-]", ""));
        out.put(row, "Protein",      props.get("Protein"));
        out.put(row, "isDecoy",      props.get("Protein").indexOf("XXX_")==0);
        out.put(row, "Search", search);
        out.put(row, "Tag", tag);
        out.put(row, "Run", props.get("#SpecFile").substring(0,props.get("#SpecFile").indexOf('.')));

        String title = props.get("Title");
        if (props.get("ScanNum").equals("-1") && title!=null && title.indexOf("NativeID")>0)
        {
          // GR_LCL_1D_ETHCD_immunopeptidome.6327.6327.2 File:"GR_LCL_1D_ETHCD_immunopeptidome.raw", NativeID:"controllerType=0 controllerNumber=1 scan=6327"
          String ID = title.substring(title.indexOf("NativeID:")+9);
          // remove the quotes
          if (ID.indexOf('"')>=0) ID = ID.substring(ID.indexOf('"')+1, ID.lastIndexOf('"'));
          out.put(row, "Scan",   ID.substring(ID.indexOf("scan=")+5));
          out.put(row, "SpecID", ID);
        }
        else
        {
          out.put(row, "Scan",   Tools.getInt(props, "ScanNum"));
          out.put(row, "SpecID", props.get("SpecID"));
        }
      }
      else if (Tools.isSet(refs) && Tools.isSet(scans) && Sets.intersection(scns, refs.keySet()).size()>0)
      {
        good_miss+=scans.length;
        if (scans.length>2)
        {
          System.out.println(row+", "+peptide+"@"+props.get("QValue")+" m/z"+Tools.d2s(Tools.getDouble(props, "Precursor"), 4)+" @z"+Tools.getInt(props, "Charge"));
          for (String scn : scns)
            System.out.println("    "+scn+"\t"+refs.get(scn));
        }
      }
    }

    int annotated=0, alls=0;
    for (Integer cnt : annot.keySet()) annotated+=cnt*annot.get(cnt).size();
    for (Integer cnt :   all.keySet())      alls+=cnt*  all.get(cnt).size();

    System.out.println("Annot/All: " + annotated+"/"+alls + "--> " + Tools.d2s(100d*annotated/alls, 2));

    return out;
  }
}
