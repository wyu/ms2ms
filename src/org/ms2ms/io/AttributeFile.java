package org.ms2ms.io;

import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.ms2ms.alg.Parser;
import org.ms2ms.data.NameValue;
import org.ms2ms.utils.TabFile;
import org.ms2ms.utils.Tools;

import java.io.*;
import java.util.Collection;

/**
 * Authos: wyu
 * <p/>
 * Synposis: a parameter file consists of name-value pair
 *
 * <p/>
 * <p/>
 * Comments:
 */
public class AttributeFile
{
  protected Multimap<String, String> mParamTemplate = TreeMultimap.create();
  protected String delimiter = TabFile.semicolumn;

  public AttributeFile()                { super(); }
  public AttributeFile(String cfg_file) { super(); setTemplate(cfg_file);}
  public AttributeFile(File cfg)        { super(); setTemplate(cfg.getPath());}
  public AttributeFile(String cfg_file, String delim)
  {
    super();
    delimiter = delim;
    setTemplate(cfg_file);
  }

  public void set(String name, String  val) { mParamTemplate.removeAll(name); mParamTemplate.put(name, val); }
  public void add(String name, String  val) { mParamTemplate.put(name, val); }
  public void set(String name, Float   val) { if (val != null) set(name, val.toString()); }
  public void add(String name, Float   val) { if (val != null) add(name, val.toString()); }
  public void set(String name, Double  val) { if (val != null) set(name, val.toString()); }
  public void add(String name, Double  val) { if (val != null) add(name, val.toString()); }
  public void set(String name, Integer val) { if (val != null) set(name, val.toString()); }
  public void add(String name, Integer val) { if (val != null) add(name, val.toString()); }

  public Multimap<String, String> getParams() { return mParamTemplate; }
  /**
   * initiate the class with the default configuration
   *
   * @param cfg_file = location of the default configuration
   */
  public void setTemplate(String cfg_file)
  {
    if (!Tools.isSet(cfg_file))
    {
      throw new RuntimeException("Couldn't figure out which chaperon config file to use!");
    }

    File defaultParamFile = new File(cfg_file);
    if (!defaultParamFile.exists())
    {
      throw new RuntimeException("The default template file '" + cfg_file + "' doesn't exist!");
    }

    BufferedReader is = null;
    String line = null;
    try
    {
      try
      {
        is = new BufferedReader(new FileReader(cfg_file));
        while ((line = is.readLine()) != null)
        {
          // skip the line with "//" at the front
          if (line.indexOf("//") == 0) continue;

          // parse the fields
          NameValue nval = Parser.newNameValue(line, ':');
          if (nval != null && Tools.isSet(nval.val))
          {
            set(nval.name.toLowerCase(), nval.val);
          }
          /*String[] items = line.split(delimiter);
          trim(items);
          if (items.length > 1)
          {
            mParamTemplate.set(items[0].toLowerCase(), items[1]);
          } */
        }
      }
      finally { if (is != null) is.close(); }
    }
    catch (Exception e) { e.printStackTrace(); }
  }
  public void write(String filename) throws IOException
  {
    write(new File(filename));
  }
  public void write(File file) throws IOException
  {
    // prepare for the output
    BufferedWriter writer = null;
    try
    {
       writer = new BufferedWriter(new FileWriter(file));
       // specification about the location of the msms input
       for (String key : mParamTemplate.keySet())
       {
         Collection<String> vals = mParamTemplate.get(key);
         for (String val : vals)
         {
            writer.write(key);
            writer.write(delimiter + " ");
            writer.write(val);
            writer.newLine();
         }
       }
    }
    finally
    {
       if (writer != null) writer.close();
    }
  }
}
