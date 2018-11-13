package org.ms2ms.io;

import com.google.common.base.Optional;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.io.filefilter.WildcardFileFilter;
import org.expasy.mzjava.core.io.ms.spectrum.MgfReader;
import org.expasy.mzjava.core.io.ms.spectrum.MgfWriter;
import org.expasy.mzjava.core.mol.NumericMass;
import org.expasy.mzjava.core.ms.peaklist.Peak;
import org.expasy.mzjava.core.ms.peaklist.PeakList;
import org.expasy.mzjava.core.ms.spectrum.*;
import org.expasy.mzjava.proteomics.io.mol.FastaProteinReader;
import org.expasy.mzjava.proteomics.mol.AminoAcid;
import org.expasy.mzjava.proteomics.mol.Protein;
import org.expasy.mzjava.proteomics.mol.modification.ModAttachment;
import org.expasy.mzjava.proteomics.mol.modification.Modification;
import org.expasy.mzjava.proteomics.ms.ident.ModificationMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.PeptideProteinMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.data.Binary;
import org.ms2ms.data.Features;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.data.ms.MsSpectrum;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.utils.BufferedRandomAccessFile;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.Tools;
import uk.ac.ebi.jmzml.model.mzml.*;

import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.net.URI;
import java.util.*;

/** Responsible for writing out the MS objects in binary format. Not to be confused with import duty
 *
 * User: wyu
 * Date: 9/29/14
 */
public class MsIO extends IOs
{
  public static void writeIntStrFeatures(DataOutput ds, MultiTreeTable<Integer, String, Features> data) throws IOException
  {
    write(ds, Tools.isSet(data) ? data.keySet().size() : 0);
    if (Tools.isSet(data))
      for (Integer row : data.keySet())
        for (String col : data.row(row).keySet())
        {
          write(ds, row);
          write(ds, col);
          write(ds, data.get(row, col));
        }
  }

  public static MsSpectrum readSpectrumIdentifier(RandomAccessFile w, long offset)
  {
    try
    {
      w.seek(offset);
      return readSpectrumIdentifier(w);
    }
    catch(IOException i)
    {
      i.printStackTrace();
    }
    return null;
  }

  public static MsSpectrum readSpectrumIdentifier(RandomAccessFile w)
  {
    try
    {
      byte[] bs = new byte[w.readInt()]; w.read(bs);
      MsSpectrum spec = MsSpectrum.fromBytes(bs); bs=null;
      return spec;
    }
    catch(IOException i)
    {
//      i.printStackTrace();
    }
    return null;
//    catch(ClassNotFoundException c)
//    {
//      System.out.println("Employee class not found");
//      c.printStackTrace();
//    }
  }
  public static long write(RandomAccessFile w, MsSpectrum ms)
  {
    try
    {
      long p1 = w.getFilePointer();
      byte[] bs = MsSpectrum.toBytes(ms);
      w.writeInt(bs.length); w.write(bs); bs=null;
      return p1;
    }
    catch (IOException e)
    { throw new RuntimeException("Error during persistence", e); }
  }

  // BufferedWriter
  public static void write(DataOutput w, PeakList ms) throws IOException
  {
    w.writeDouble(ms.getPrecursor()!=null?ms.getPrecursor().getMz():0d);
    w.writeDouble(ms.getPrecursor()!=null?ms.getPrecursor().getIntensity():0d);
    w.writeInt(   ms.getPrecursor()!=null?ms.getPrecursor().getCharge():0);

    w.writeInt(   ms.size());
    if (ms.size()>0)
    {
      for (int i=0; i<ms.size(); i++)
      {
        w.writeDouble(ms.getMz(i));
        w.writeDouble(ms.getIntensity(i));
      }
    }
  }
  public static void write(DataOutput w, UUID s) throws IOException
  {
    write(w, s!=null?1:0);
    if (s!=null) write(w, s.toString());
  }
  public static UUID readUUID(DataInput w) throws IOException
  {
    int n = read(w, 0);
    return (n==1?UUID.fromString(read(w,"")):null);
  }
  public static void write(DataOutput w, URI s) throws IOException
  {
    write(w, s!=null?1:0);
    if (s!=null) write(w, s.toString());
  }
  public static URI readURI(DataInput w) throws IOException
  {
    int n = read(w, 0);
    return (n==1?URI.create(read(w,"")):null);
  }
  public static void write(DataOutput w, PeakList.Precision s) throws IOException
  {
    write(w, s!=null?1:0);
    if (s!=null) write(w, s.toString());
  }
  public static PeakList.Precision readPrecision(DataInput w) throws IOException
  {
    int n = read(w, 0);
    return (n==1? PeakList.Precision.valueOf(read(w,"")):null);
  }

  public static void writeMsnMS(DataOutput w, MsnSpectrum ms) throws IOException
  {
    // peaks and RT only
    write(w, ms);
    // write the rest of the info. use a diff call signature to avoid messing up the old codes
    write(w, ms.getComment());
    write(w, ms.getScanNumbers().getFirst().getValue());
    write(w, ms.getMsLevel());
    write(w, ms.getSpectrumIndex());
    write(w, ms.getFragMethod());
    write(w, ms.getId());
    write(w, ms.getSpectrumSource());
  }
  public static MsnSpectrum readMsnMS(DataInput w) throws IOException
  {
    // peaks and RT only
    MsnSpectrum ms = readSpectrumIdentifier(w, new MsnSpectrum());

    // write the rest of the info. use a diff call signature to avoid messing up the old codes
    ms.setComment(read(w,""));
    ms.addScanNumber(read(w,0));
    ms.setMsLevel(read(w,0));
    ms.setSpectrumIndex(read(w,0));
    ms.setFragMethod(read(w,""));
    ms.setId(readUUID(w));
    ms.setSpectrumSource(readURI(w));

    return ms;
  }
  public static void write(DataOutput w, MsnSpectrum ms) throws IOException
  {
    write(w, (PeakList) ms);

    int counts = ms.getRetentionTimes()!=null&&ms.getRetentionTimes().size()>0?ms.getRetentionTimes().size():0;
    w.writeInt(counts);
    if (counts>0)
      for (int i=0; i<counts; i++) w.writeDouble(ms.getRetentionTimes().get(i).getTime());
  }
  public static void
  writeStrSpectrumMap(DataOutput ds, Map<String, MsnSpectrum> data) throws IOException
  {
    write(ds, Tools.isSet(data) ? data.keySet().size() : 0);

    if (Tools.isSet(data))
      for (String key : data.keySet())
      {
        write(ds, key);
        writeMsnMS(ds, data.get(key));
      }
  }
  public static Map<String, MsnSpectrum> readStrSpectrumMap(DataInput ds) throws IOException
  {
    int n = read(ds, 0);
    if (n > 0)
    {
      Map<String, MsnSpectrum> data = new TreeMap<>();
      for (int i = 0; i < n; i++)
      {
        String   K =read(ds, "");
        data.put(K, readMsnMS(ds));
      }
      return data;
    }
    return null;
  }

  public static PeakList readSpectrumIdentifier(DataInput w, PeakList ms) throws IOException
  {
    if (ms!=null)
    {
      ms.setPrecursor(new Peak(w.readDouble(), w.readDouble(), w.readInt()));

      int counts = w.readInt(); ms.clear();
      if (counts>0)
        for (int i=0; i<counts; i++)
          ms.add(w.readDouble(), w.readDouble());
    }
    return ms;
  }
  public static Map<Float,Float> readSpectrumIdentifierAsIonMap(DataInput w, Map<Float,Float> ms) throws IOException
  {
    if (ms!=null)
    {
      w.readDouble(); w.readDouble(); w.readInt();

      int counts = w.readInt(); ms.clear();
      if (counts>0)
        for (int i=0; i<counts; i++)
          ms.put((float )w.readDouble(), (float )w.readDouble());
    }
    return ms;
  }
  public static MsnSpectrum readSpectrumIdentifier(BufferedRandomAccessFile w, MsnSpectrum ms, long offset) throws IOException
  {
    w.seek(offset);
    return readSpectrumIdentifier(w, ms);
  }
  public static Map<Float,Float> readSpectrumIdentifierAsIonMap(BufferedRandomAccessFile w, long offset) throws IOException
  {
    w.seek(offset);
    return readSpectrumIdentifierAsIonMap(w);
  }
  public static MsnSpectrum readSpectrumIdentifier(DataInput w, MsnSpectrum ms) throws IOException
  {
    if (ms!=null)
    {
      ms = (MsnSpectrum ) readSpectrumIdentifier(w, (PeakList) ms);

      int counts = w.readInt();
      // clear the content of the retention times if necessary
      if (ms.getRetentionTimes()!=null) ms.getRetentionTimes().clear();
      if (counts>0)
        for (int i=0; i<counts; i++)
          try
          {
            ms.addRetentionTime(new RetentionTimeDiscrete(w.readDouble(), TimeUnit.SECOND));
          }
          catch (Exception e)
          {
            e.printStackTrace();
          }
    }
    return ms;
  }
  public static Map<Float,Float> readSpectrumIdentifierAsIonMap(DataInput w) throws IOException
  {
    Map<Float,Float> ms = readSpectrumIdentifierAsIonMap(w, new HashMap<Float,Float>());

    int counts = w.readInt();
    // clear the content of the retention times if necessary
    if (counts>0)
      for (int i=0; i<counts; i++)
        try
        {
          w.readDouble();
        }
        catch (Exception e)
        {
          e.printStackTrace();
        }
    return ms;
  }
  public static List<MsnSpectrum> readSpectra(String s) { return readSpectra(s, null); }
  public static List<MsnSpectrum> readSpectra(String s, long[] offsets)
  {
    RandomAccessFile F = null;
    try
    {
      try
      {
        F = new RandomAccessFile(s, "r");
        List<MsnSpectrum> spec = Tools.isSet(offsets) ? readSpectra(F, offsets) : readSpectra(F);
        F.close();
        return spec;
      }
      catch (FileNotFoundException fne)
      {
        throw new RuntimeException("Not able to locate the file: " + s, fne);
      }
      finally
      {
        if (F!=null) F.close();
      }
    }
    catch (IOException ie) { throw new RuntimeException("Error while reading the spectra", ie); }
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s) throws IOException
  {
    if (s==null) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>();
    try
    {
      while (1==1)
      {
        spectra.add(readSpectrumIdentifier(s).toMsnSpectrum());
      }
    }
    catch (Exception ie) {}
    return spectra;
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s, long[] offsets) throws IOException
  {
    if (s==null || !Tools.isSet(offsets)) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>(offsets.length);
    for (long offset : offsets)
    {
      s.seek(offset);
      spectra.add(readSpectrumIdentifier(s).toMsnSpectrum());
    }
    return spectra;
  }
  public static List<MsnSpectrum> readSpectra(RandomAccessFile s, Map<Long, String> offset_row) throws IOException
  {
    if (s==null || !Tools.isSet(offset_row)) return null;
    // the output
    List<MsnSpectrum> spectra = new ArrayList<>(offset_row.size());
    for (Long offset : offset_row.keySet())
    {
      s.seek(offset);
      MsnSpectrum spec = readSpectrumIdentifier(s).toMsnSpectrum();
      spec.setComment(offset_row.get(offset)); spectra.add(spec);
    }
    return spectra;
  }
  public static void writeSpectra(String s, Collection<MsnSpectrum> spectra) throws IOException
  {
    if (s==null || !Tools.isSet(spectra)) return;
    // the output
    RandomAccessFile F = null;
    try
    {
      try
      {
        F = new RandomAccessFile(s, "rw");
        for (MsnSpectrum m : spectra)
        {
          write(F, MsSpectrum.adopt(m));
        }
        F.close();
        return;
      }
      catch (FileNotFoundException fne)
      {
        throw new RuntimeException("Not able to locate the file: " + s, fne);
      }
      finally
      {
        if (F!=null) F.close();
      }
    }
    catch (IOException ie) { throw new RuntimeException("Error while writing the spectra", ie); }
  }
  public static void writeAAs(DataOutput w, List<AminoAcid> seqs) throws IOException
  {
    write(w, seqs.size());
    for (int i=0; i<seqs.size(); i++) write(w, seqs.get(i).getSymbol());
  }
  public static void write(DataOutput w, Modification m) throws IOException
  {
    write(w, m.getLabel());
    write(w, m.getMolecularMass());
  }
  public static Modification readModification(DataInput w) throws IOException
  {
    String label = read(w, "");
    double  mass = read(w, 0d);

    return new Modification(label, new NumericMass(mass));
  }
  public static void writeMods(DataOutput w, List<Modification> mods) throws IOException
  {
    write(w, mods.size());
    for (int i=0; i<mods.size(); i++) write(w, mods.get(i));
  }
  public static List<AminoAcid> readAAs(DataInput w) throws IOException
  {
    int size = read(w, 0);
    List<AminoAcid> seqs = new ArrayList<AminoAcid>(size);
    for (int i=0; i<size; i++) seqs.add(AminoAcid.valueOf(read(w, "")));
    return seqs;
  }
  public static void write(DataOutput w, ModificationMatch m) throws IOException
  {
    write(w, m.getMassShift());
    write(w, m.getResidue().getSymbol());
    write(w, m.getPosition());

    write(w, m.getCandidateCount());
    if (m.getCandidateCount()>0)
      for (int i=0; i<m.getCandidateCount(); i++) write(w, m.getModificationCandidate(i));

    write(w, m.getModAttachment().toString());
  }
  public static ModificationMatch readModificationWatch(DataInput w) throws IOException
  {
    double mass_shift = read(w, 0d);
    String     symbol = read(w, "");
    int      position = read(w, 0), counts = read(w, 0);

//    ModificationMatch m = new ModificationMatch(mass_shift, AminoAcid.valueOf(symbol), position, ModAttachment.SIDE_CHAIN);
//
    List<Modification> mods = new ArrayList<>();
    if (counts>0)
      for (int i=0; i<counts; i++)
        mods.add(readModification(w));

    // has to make another copy because the ModAttachment was written in the wrong order
    ModificationMatch m = new ModificationMatch(mass_shift, AminoAcid.valueOf(symbol), position, ModAttachment.valueOf(read(w, "")));

    if (counts>0)
      for (int i=0; i<counts; i++)
        m.addPotentialModification(mods.get(i));

    return m;
  }
  public static void writeModMatches(DataOutput w, Collection<ModificationMatch> mods) throws IOException
  {
    write(w, mods.size());
    if (mods.size()>0)
      for (ModificationMatch m : mods) write(w, m);
  }
  public static Collection<ModificationMatch> readModMatches(DataInput w) throws IOException
  {
    Collection<ModificationMatch> mods = new ArrayList<>();
    int size = read(w, 0);
    if (size>0)
      for (int i=0; i<size; i++) mods.add(readModificationWatch(w));

    return mods;
  }
/*
  public static void writeSidechainMatches(DataOutput w, ListMultimap<Integer, ModificationMatch> map) throws IOException
  {
    write(w, map.keySet().size());

    if (map.size()>0)
      for (Integer i : map.keySet())
      {
        write(w, i); write(w, map.get(i).size());
        for (ModificationMatch m : map.get(i)) write(w, m);
      }
  }
  public static ListMultimap<Integer, ModificationMatch> readSidechainMatches(DataInput w) throws IOException
  {
    int size = readSpectrumIdentifier(w, 0);
    ListMultimap<Integer, ModificationMatch> map = ArrayListMultimap.create();

    if (size>0)
      for (int i=0; i<size; i++)
      {
        Integer k=readSpectrumIdentifier(w, 0), s=readSpectrumIdentifier(w, 0);
        for (int j=0; j<s; j++)
          map.put(k, readModificationWatch(w));
      }

    return map;
  }
  public static void writeTermMatches(DataOutput w, ListMultimap<ModAttachment, ModificationMatch> map) throws IOException
  {
    write(w, map.keySet().size());

    if (map.size()>0)
      for (ModAttachment i : map.keySet())
      {
        write(w, i.toString()); write(w, map.get(i).size());
        for (ModificationMatch m : map.get(i)) write(w, m);
      }
  }
  public static ListMultimap<ModAttachment, ModificationMatch> readTermMatches(DataInput w) throws IOException
  {
    int size = readSpectrumIdentifier(w, 0);
    ListMultimap<ModAttachment, ModificationMatch> map = ArrayListMultimap.create();

    if (size>0)
      for (int i=0; i<size; i++)
      {
        ModAttachment mt = ModAttachment.valueOf(readSpectrumIdentifier(w, "")); Integer s=readSpectrumIdentifier(w, 0);
        for (int j=0; j<s; j++)
          map.put(mt, readModificationWatch(w));
      }

    return map;
  }
*/
  public static void write(DataOutput w, PeptideProteinMatch pm) throws IOException
  {
    write(w, pm.getHitType().toString());
    write(w, pm.getAccession());
    write(w, pm.getStart());
    write(w, pm.getEnd());

    writeOpStr(w, pm.getSearchDatabase());
    writeOpStr(w, pm.getPreviousAA());
    writeOpStr(w, pm.getNextAA());
  }
  public static PeptideProteinMatch readPeptideProteinMatch(DataInput w) throws IOException
  {
    PeptideProteinMatch.HitType hitType = PeptideProteinMatch.HitType.valueOf(read(w, ""));
    String accession = read(w, "");
    int start=read(w, 0), end=read(w, 0);

    Optional<String> db=readOpStr(w), prev=readOpStr(w), next=readOpStr(w);

    // make sure we have a valid range
    if (start==-1 && end<=start) { start=1; end=Integer.MAX_VALUE; }

    return new PeptideProteinMatch(accession, db, prev, next, start, end, hitType);
  }
  public static void write(DataOutput w, PeptideMatch m) throws IOException
  {
    write(w, m.toSymbolString());
    write(w, m.getRank());
    write(w, m.getNumMatchedIons());
    write(w, m.getTotalNumIons());
    write(w, m.getMassDiff());
    write(w, m.getNumMissedCleavages());
    write(w, m.isRejected());
    try
    {
      // there is no way to get the Optional(Double) directly
      write(w, m.getNeutralPeptideMass());
    }
    catch (IllegalStateException e)
    {
      write(w, 0d);
    }
    writeModMatches(w, m.getModifications(ModAttachment.all));

    List<PeptideProteinMatch> matches = m.getProteinMatches();
    write(w, matches!=null?matches.size():0);
    for (PeptideProteinMatch pm : matches) write(w, pm);

    write(w, m.getScoreMap()!=null?m.getScoreMap().keySet().size():0);
    for (String tag : m.getScoreMap().keySet())
    {
      write(w, tag); write(w, m.getScore(tag));
    }
  }
  public static PeptideMatch readPeptideMatch(DataInput w) throws IOException
  {
    PeptideMatch pm = new PeptideMatch(read(w, ""));

    pm.setRank(read(w, 0));
    pm.setNumMatchedIons(read(w, 0));
    pm.setTotalNumIons(read(w, 0));
    pm.setMassDiff(read(w, 0d));
    pm.setNumMissedCleavages(read(w, 0));
    pm.setRejected(read(w, false));

    double nm = read(w, 0d);
    if (nm>0) pm.setNeutralPeptideMass(nm);

    Collection<ModificationMatch> matches = readModMatches(w);
    if (Tools.isSet(matches))
      for (ModificationMatch m : matches)
        // TODO, not sure whether position is the same as index
        pm.addModificationMatch(m.getPosition(), m);

    int n_pm = read(w, 0);
    if (n_pm>0)
    {
      for (int i=0; i<n_pm; i++) pm.addProteinMatch(readPeptideProteinMatch(w));
    }

    int n_score = read(w, 0);
    if (n_score>0)
    {
      for (int i=0; i<n_score; i++) pm.addScore(read(w, ""), read(w, 0d));
    }
    return pm;
  }
  public static <R extends RetentionTime> void write(DataOutput w, R rt) throws IOException
  {
    write(w, rt.getMinRetentionTime());
    write(w, rt.getMaxRetentionTime());
    write(w, rt.getTime());
  }
  public static void writeRetentionTimeList(DataOutput w, RetentionTimeList RL) throws IOException
  {
    write(w, RL!=null?RL.size():0);
    for (int i=0; i<RL.size(); i++) write(w, RL.get(i));
  }
  public static <S extends ScanNumber> void write(DataOutput w, S scan) throws IOException
  {
    write(w, scan.getMinScanNumber());
    write(w, scan.getMaxScanNumber());
    write(w, scan.getValue());
  }
  public static void writeScanNumberList(DataOutput w, ScanNumberList scans) throws IOException
  {
    write(w, scans != null ? scans.size() : 0);
    for (int i=0; i<scans.size(); i++) write(w, scans.get(i));
  }
  public static ScanNumber readScanNumber(DataInput w) throws IOException
  {
    // still need to readSpectrumIdentifier the int
    int min=read(w, 0), max=read(w, 0), mid=read(w, 0);
    return min==max?new ScanNumberDiscrete(mid):new ScanNumberInterval(min, max);
  }
  public static RetentionTime readRetentionTime(DataInput w) throws IOException
  {
    // still need to readSpectrumIdentifier the int
    double min=read(w, 0d), max=read(w, 0d), mid=read(w, 0d);
    return min==max?new RetentionTimeDiscrete(mid,      TimeUnit.SECOND):
                    new RetentionTimeInterval(min, max, TimeUnit.SECOND);
  }
  public static RetentionTimeList readRetentionTimeList(DataInput w) throws IOException
  {
    int size = read(w, 0);
    RetentionTimeList RT = new RetentionTimeList();
    for (int i=0; i<size; i++) RT.add(readRetentionTime(w));

    return RT;
  }
  public static ScanNumberList readScanNumberList(DataInput w) throws IOException
  {
    int size = read(w, 0);
    ScanNumberList scans = new ScanNumberList();
    for (int i=0; i<size; i++) scans.add(readScanNumber(w));

    return scans;
  }
  public static void write(DataOutput w, SpectrumIdentifier id) throws IOException
  {
    write(        w, id.getSpectrum());
    writeOpStr(   w, id.getName());
    writeOpDouble(w, id.getPrecursorNeutralMass());
    writeOpDouble(w, id.getPrecursorIntensity());
    writeOpDouble(w, id.getPrecursorMz());
    writeOpInt(   w, id.getAssumedCharge());
    writeOpInt(   w, id.getIndex());
    writeOpStr(   w, id.getSpectrumFile());
    writeScanNumberList(   w, id.getScanNumbers());
    writeRetentionTimeList(w, id.getRetentionTimes());
  }
  public static SpectrumIdentifier readSpectrumIdentifier(DataInput w) throws IOException
  {
    SpectrumIdentifier id = new SpectrumIdentifier(read(w, ""));

    if (read(w, true)) id.setName(                read(w, "unTitled"));
    if (read(w, true)) id.setPrecursorNeutralMass(read(w, 0d));
    if (read(w, true)) id.setPrecursorIntensity(  read(w, 0d));
    if (read(w, true)) id.setPrecursorMz(         read(w, 0d));
    if (read(w, true)) id.setAssumedCharge(       read(w, 0));
    if (read(w, true)) id.setIndex(               read(w, 0));
    if (read(w, true)) id.setSpectrumFile(        read(w, ""));
    id.addScanNumber(   readScanNumberList(w));
    id.addRetentionTime(readRetentionTimeList(w));

    return id;
  }
  public static String get(Collection<CVParam> params, String key)
  {
    if (Tools.isSet(params) && Strs.isSet(key))
      for (Iterator lCVParamIterator = params.iterator(); lCVParamIterator.hasNext();)
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals(key))
        {
          return lCVParam.getValue().trim();
        }
      }
    return null;
  }
  public static Integer getInt(   Collection<CVParam> params, String key) { return Stats.toInt(   get(params, key)); }
  public static Double  getDouble(Collection<CVParam> params, String key) { return Stats.toDouble(get(params, key)); }

  public static String getIf(Collection<CVParam> params, String tag_required, String val_required)
  {
    if (Tools.isSet(params) && Strs.isSet(tag_required) && Strs.isSet(val_required))
      for (Iterator lCVParamIterator = params.iterator(); lCVParamIterator.hasNext();)
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals(tag_required))
        {
          String unitRT = lCVParam.getUnitAccession().trim();
          if (unitRT.equals(val_required)) //... Validating rt unit (mins or secs) ...//
          {
            return lCVParam.getValue().trim();
          }
        }
      }
    return null;
  }
  public static boolean hasAccession(Collection<CVParam> params, String tag_required)
  {
    if (Tools.isSet(params) && Strs.isSet(tag_required))
      for (Iterator lCVParamIterator = params.iterator(); lCVParamIterator.hasNext();)
      {
        CVParam lCVParam = (CVParam) lCVParamIterator.next();
        if (lCVParam.getAccession().equals(tag_required)) return true;
      }
    return false;
  }
  public static Integer ScanNumberFromSpectrumRef(String ref)
  {
    String[] ids = Strs.split(ref, ' ', true);
    for (String id : ids)
      if (id.indexOf("scan=")==0) return Stats.toInt(id.substring(5));

    return null;
  }
  public static Map<String, MsnSpectrum> readMGF(String mgffile) throws IOException
  {
    MgfReader                    mgf = new MgfReader(new File(mgffile), PeakList.Precision.DOUBLE);
    Map<String, MsnSpectrum> spectra = new HashMap<>();

    while (mgf.hasNext())
    {
      MsnSpectrum ms = mgf.next();
      spectra.put(ms.getScanNumbers().getFirst().getValue()+"",ms);
    }
    mgf.close();
    return spectra;
  }
  public static void writeMGF(String mgf, Collection<MsnSpectrum> spectra, int min_dup) throws IOException
  {
    System.out.println("Writing "+spectra.size()+" MS/MS spectra with "+min_dup+" or more members to "+mgf);

    MgfWriter MGF = new MgfWriter(new File(mgf), PeakList.Precision.DOUBLE);
    for (MsnSpectrum ms : spectra)
      if (ms.getSpectrumIndex()>=min_dup) MGF.write(ms);

    MGF.close();
  }

  //  // write the content of the spectrum from an mzML file to another MGF file
//  public static void writeMGF(MgfWriter w, uk.ac.ebi.jmzml.model.mzml.Spectrum ms) throws IOException
//  {
//    if (w==null || ms==null) return;
//
//
//
//
//  }
  public static Multimap<String, Protein> readFASTA(String filename)
  {
    System.out.println("Reading the protein sequences from " + filename);

    FastaProteinReader fasta = null;
    try
    {
      Multimap<String, Protein> seq_proteins = HashMultimap.create();
      try
      {
        fasta = new FastaProteinReader(new File(filename));
        while (fasta.hasNext())
        {
          Protein p = fasta.next();
          seq_proteins.put(p.toSymbolString(), p);
        }
      }
      finally
      {
        if (fasta!=null) fasta.close();
      }
      return seq_proteins;
    }
    catch (IOException e) {}

    return null;
  }
  public static Dataframe writeScanPeptideMatchesForRun(
       RandomAccessFile w, Multimap<SpectrumIdentifier, PeptideMatch> scan_matches) throws IOException
  {
    // create the binary dump and index
    Dataframe index = new Dataframe("ScanIndex");
    // order the scans first
    SortedMap<Integer, SpectrumIdentifier> scan_id = new TreeMap<>();
    for (SpectrumIdentifier id : scan_matches.keySet())
      scan_id.put(id.getScanNumbers().getFirst().getValue(), id);

    // write the number of spectrum ID first
    write(w, scan_matches.keySet().size());
    // write the scans sequentially
    for (Integer scan : scan_id.keySet())
    {
      SpectrumIdentifier id = scan_id.get(scan);
      // add the row to the index
      index.addRow(
          "Scan", scan,
          "RT", Tools.isSet(id.getRetentionTimes())?id.getRetentionTimes().getFirst().getTime():-1,
          "mz", id.getPrecursorMz(),
          "z", id.getAssumedCharge().isPresent() ? id.getAssumedCharge().get() : 0,
          "FileOffset", w.getFilePointer(),
          "N.matches", scan_matches.get(id).size());
      // dump the contents to the binary file
      write(w, scan);
      write(w, scan_matches.get(id).size());
      write(w, id);
      for (PeptideMatch m : scan_matches.get(id)) MsIO.write(w, m);
    }
    return index;
  }
  public static long writeScanPeptideMatches(
      RandomAccessFile w, Multimap<SpectrumIdentifier, PeptideMatch> scan_matches) throws IOException
  {
    // write the number of spectrum ID first
    write(w, scan_matches.keySet().size());

    long counts=0;
    // write the scans sequentially
    for (SpectrumIdentifier id : scan_matches.keySet())
    {
      // dump the contents to the binary file
      write(w, id.getScanNumbers().getFirst().getValue());
      write(w, scan_matches.get(id).size());
      write(w, id);
      for (PeptideMatch m : scan_matches.get(id)) { counts++; MsIO.write(w, m); }
    }
    return counts;
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch> loadScanPeptideMatches(String filename)
  {
    System.out.print("Loading from " + filename);
    try (RandomAccessFile w = new RandomAccessFile(filename, "r"))
    {
      Multimap<SpectrumIdentifier, PeptideMatch> matches = MsIO.readScanPeptideMatches(w, Integer.MAX_VALUE);
      System.out.print(" (" + matches.keySet().size() + "/" + matches.size() + ") ");

      return matches;
    }
    catch (Exception e)
    {
      throw new RuntimeException(e);
    }
  }
  public static Multimap<SpectrumIdentifier, PeptideMatch>
      readScanPeptideMatches(DataInput w, int n) throws IOException
  {
    Multimap<SpectrumIdentifier, PeptideMatch> scan_matches = HashMultimap.create();
    // fetch as many batches as asked till the end of the file
    try
    {
      int spectra=0;
      spectra = read(w, spectra);
      for (int i=0; i<spectra; i++)
      {
        // dump the contents to the binary file
        Integer scan = read(w, 0), nmatches = read(w, 0);
        SpectrumIdentifier id = readSpectrumIdentifier(w);
        for (int k=0; k<nmatches; k++)
          scan_matches.put(id, readPeptideMatch(w));
      }
    }
    catch (EOFException eof) {}
    catch (IOException e)
    {
      e.printStackTrace();
    }
    return scan_matches;
  }
  // expand the file with wildcard to individual file names with full path
  public static List<String> listFiles(String root)
  {
    File template = new File(root);
    return listFiles(template.getParent(), new WildcardFileFilter(template.getName()), 1);
  }
  public static void writeTsv(Writer w, MsnSpectrum ms) throws IOException
  {
    w.write("mz\tai\n");
    for (int i=0; i<ms.size(); i++)
      w.write(ms.getMz(i)+"\t"+ms.getIntensity(i)+"\n");
  }
  public static void writeTsv(String out, MsnSpectrum ms)
  {
    FileWriter w = null;
    try
    {
      w = new FileWriter(out);
      writeTsv(w, ms);
      w.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
}
