package org.ms2ms.data.ms;

import com.google.common.collect.TreeBasedTable;
import com.google.common.collect.TreeMultimap;
import org.ms2ms.data.Binary;
import org.ms2ms.data.Features;
import org.ms2ms.data.collect.MultiTreeTable;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

public class LcMsMsNeighborhood implements Binary
{
  private MultiTreeTable<String, String, Features> mRunPeerLibrary;
  private TreeMultimap<Double, Features> mMassLibrary;

  public LcMsMsNeighborhood() {super(); }
  public LcMsMsNeighborhood(LcMsMsFeatures s) {super(); indice(s); }


  public SortedMap<Double, Collection<Features>> slice(Double m0, Double m1)
  {
     return mMassLibrary.asMap().subMap(m0,m1);
  }
  public Collection<Features> get(String run, String peer)
  {
    return mRunPeerLibrary.get(run, peer);
  }
  public LcMsMsNeighborhood indice(LcMsMsFeatures lcmsms)
  {
    mRunPeerLibrary = MultiTreeTable.create();
    mMassLibrary    = TreeMultimap.create();

    Map<String, TreeBasedTable<Double, Double, Features>> RunRtAiIons = lcmsms.indexing();
    int rows=0;
    for (String run : RunRtAiIons.keySet())
      for (Features F : RunRtAiIons.get(run).values())
      {
        if (F.get(LcMsMsFeatures.COL_PEPTIDE)!=null)
        {
          if (++rows%100==0) System.out.print(  ".");
          if (rows%10000==0) System.out.println("");

          mRunPeerLibrary.put(run, F.getStr(LcMsMsFeatures.COL_PEPTIDE), F.clone());
          mMassLibrary.put(        F.getDouble(LcMsMsFeatures.COL_MASS), F.clone());
        }
      }

    System.out.println();

    RunRtAiIons=(TreeMap) Tools.dispose(RunRtAiIons);
    return this;
  }

  @Override
  public void write(DataOutput ds) throws IOException
  {
    IOs.writeStr2StrMultiTable(ds, mRunPeerLibrary);
    IOs.writeDoubleMultimap(   ds, mMassLibrary);
  }

  @Override
  public void read(DataInput ds) throws IOException
  {
    mRunPeerLibrary = IOs.readStr2StrMultiTable(ds, Features.class);
    mMassLibrary    = (TreeMultimap)IOs.readDoubleMaps(ds, mMassLibrary, Features.class);
  }
}
