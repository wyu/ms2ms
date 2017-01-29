package org.ms2ms.io;

import org.expasy.mzjava.proteomics.io.mol.FastaProteinReader;
import org.expasy.mzjava.proteomics.mol.Protein;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.Reader;

/**
 * Created by yuw on 1/14/17.
 */
public class CustomFastaProteinReader extends FastaProteinReader
{
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
}
