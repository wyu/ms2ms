package org.ms2ms.utils;

import java.io.*;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * User: wyu
 * Date: 6/12/14
 * Time: 11:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class TabFile
{
  private int                mLineSkip   = 0;
  private String[]           mHds        = null;
  private Map<String,String> mCol = null, mTransforms = null;
  private BufferedReader     mFileReader = null;
  private String             mDelimiter  = tabb;
  private char               mTokenChar  = 0;
  private Pattern mToken      = Pattern.compile(tabb);
  private String             mFilename   = null;
  private String             mCurrentLine = null;
  private String             mHeaderLine = null;
  private String             mSkip       = null;

  public static String       comma       = ",(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))";
  public static String       tabb        = "\\t";
  public static String       dot         = "\\.";
  public static String       semicolumn  = ":";
  public static String       column      = ";";
  public static String       space       = "\\s";
  public static String       slash       = "/";
  public static String       at          = "@";
  public static String       pipe        = "\\|";
  public static String       pminus      = "?";

  public TabFile(String              infile,
                 String              delimiter) throws IOException
  {
    // no need to trap the exceptions since they are
    // delegated to the calling party
    mFilename   = infile;
    setDelimiter(delimiter);

    init();
  }
  public TabFile(String              infile,
                 char                delimiter) throws IOException
  {
    // no need to trap the exceptions since they are
    // delegated to the calling party
    mFilename   = infile;
    setDelimiter(null);
    mTokenChar = delimiter;

    init();
  }
  public TabFile(String              infile,
                 String              delimiter,
                 int                 line_skip) throws IOException
  {
    // no need to trap the exceptions since they are
    // delegated to the calling party
    mFilename   = infile;
    setDelimiter(delimiter);
    mLineSkip = line_skip;

    init();
  }
  public TabFile(String              infile,
                 String              delimiter,
                 String              skip) throws IOException
  {
    // no need to trap the exceptions since they are
    // delegated to the calling party
    mFilename   = infile;
    setDelimiter(delimiter);
    mSkip       = skip;

    init();
  }
  public TabFile(String              infile,
                   String              delimiter,
                   String[]            cols)
      throws IOException
  {
    // no need to trap the exceptions since they are
    // delegated to the calling party
    mFilename   = infile;
    setDelimiter(delimiter);
    mHds        = cols;

    init();
  }
  protected void init() throws IOException
  {
    mCol           = new TreeMap<String, String>();
    InputStream is = new FileInputStream(mFilename);
    // Gracefully handle gzipped files.
    if (mFilename.endsWith(".gz"))
    {
      is = new GZIPInputStream(is);
    }
    mFileReader = new BufferedReader(new InputStreamReader(is));

    // skip as many lines as indicated
    if (mLineSkip > 0)
      for (int i = 0; i < mLineSkip; i++)  mFileReader.readLine();

    if (mHds == null)
    {
      // skip line starting with 'skip' if asked
      if (Tools.isSet(mSkip))
      {
        while (mFileReader.ready())
        {
          mHeaderLine = mFileReader.readLine();
          if (mHeaderLine.trim().indexOf(mSkip) != 0) break;
        }
      }
      else mHeaderLine = mFileReader.readLine();

      if (Tools.isSet(mHeaderLine))
      {
        mHds = split(mHeaderLine);
        Strs.trim(     mHds);
        Strs.dequotes( mHds, '"');
      }
    }
  }
  public String              getFileName()  { return mFilename; }
  public String              getDelimiter() { return mDelimiter; }
  public String[]            getHeaders()   { return mHds; }
  public String              getHeaderLine() { return mHeaderLine; }

  //--------------------------------------------------------------------------
  public Map<String, String> getMappedRow()
  {
    return mCol;
  }

//  public Double number( String key) { return Double.parseDouble(cells(key)); }
//  public Long   integer(String key) { return Long.parseLong(cells(key)); }

  public String getCurrentLine() { return mCurrentLine; }

  public void setDelimiter(String s)
  {
    mDelimiter =  s;
    mToken     = (s != null ? Pattern.compile(s) : null);
  }
  //--------------------------------------------------------------------------
  public boolean hasNext() throws IOException
  {
    // read-in the line in a separated step
    // skip line starting with 'skip' if asked
    if (Tools.isSet(mSkip))
    {
      while (mFileReader.ready())
      {
        mCurrentLine = mFileReader.readLine();
        if (mCurrentLine.trim().indexOf(mSkip) != 0) break;
      }
    }
    else mCurrentLine = mFileReader.readLine();

    // split the fields according to the delimiter
    String[] cols = split(mCurrentLine);
    if (mCurrentLine != null && Tools.isSet(cols))
    {
      // forget the duplicated header
      while (Tools.isSet(cols[0]) && Tools.isSet(mHds[0]) && cols[0].equals(mHds[0]))
      {
        String hd = mFileReader.readLine();
        if (hd != null) cols = split(hd);
      }
      Strs.dequotes(cols, '"');
      Strs.trim(    cols);
      // populate the col map
      mCol.clear();
      for (int i = 0; i != cols.length; i++)
        // do we need to do the transformation?
        if (i < getHeaders().length) {
          String key = getHeaders()[i];
          mCol.put(key, cols[i]);
        }
    }

    return (mCurrentLine != null);
  }

  public Map<String, String> nextRow() { return mCol; }

  protected String[] split(String s)
  {
    return s != null ? (mTokenChar != 0 ? Strs.split(s, mTokenChar) : mToken.split(s)) : null;
  }

  //--------------------------------------------------------------------------
  public boolean ready() throws IOException { return mFileReader.ready(); }

  //--------------------------------------------------------------------------
  public void close() throws IOException
  {
    if (mFileReader != null) mFileReader.close();
  }
  protected void finalize() throws IOException { close(); }
  public String get(String key)
  {
    String value = null;
    if (mCol != null && key != null)
    {
      value = mCol.get(key);
    }

    return value;
  }
}
