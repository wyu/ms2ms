package org.ms2ms.io;

import org.expasy.mzjava.proteomics.io.mol.FastaProteinReader;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.ms2ms.utils.Strs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.Reader;
import java.util.EnumMap;

import static org.ms2ms.io.CustomFastaProteinReader.SeqInfo.*;

/**
 * Created by yuw on 1/14/17.
 */
public class CustomFastaProteinReader extends FastaProteinReader
{
  public enum SeqInfo { AC, NAME, DESC, GN, OS };

  public CustomFastaProteinReader(Reader reader)
  {
    super(reader);
  }

  public CustomFastaProteinReader(File file) throws FileNotFoundException
  {
    super(file);
  }

  @Override
  protected void parseHeader(String header, Protein.Builder builder)
  {
    int left = header.indexOf('>')+1, right = header.indexOf(' ');
//    builder.setAccessionId((right>0?header.substring(left, right):header.substring(left)).replaceAll("\\|", "_"));
    builder.setAccessionId(header.substring(left).replaceAll("\\|", "_"));
  }

  public static EnumMap<SeqInfo, String> parseUniprotHeader(String header, boolean decoy)
  {
    if (!Strs.isSet(header)) return null;

    EnumMap<SeqInfo, String> parcels = new EnumMap<>(SeqInfo.class);

    // sp_O75396_SC22B_HUMAN Vesicle-trafficking protein SEC22b OS=Homo sapiens GN=SEC22B PE=1 SV=4
    int r0 = header.indexOf(" "), r1 = header.indexOf("OS=");
    parcels.put(AC, r0>0?header.substring(0, r0):header);
    // tag the decoy if not already
    if (decoy && parcels.get(AC).indexOf("DECOY")<0) parcels.put(AC, "DECOY_"+parcels.get(AC));

    if (r0>0 && r1>r0) parcels.put(DESC, header.substring(r0, r1).trim());

    String[] items = r1>0? Strs.split(header.substring(r1), ' '):null;
    if (items!=null)
      for (String item : items)
      {
        String[] nv=Strs.split(item, '=', true);
        if (nv!=null && nv.length>1 && "GN".equals(nv[0])) parcels.put(GN, decoy?"_"+nv[1]+"_":nv[1]);
        if (nv!=null && nv.length>1 && "OS".equals(nv[0])) parcels.put(OS, nv[1].trim());
      }

    return (parcels);
  }
}
