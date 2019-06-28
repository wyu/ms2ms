package org.ms2ms.io;

import com.google.common.collect.Range;
import org.ms2ms.utils.Strs;

public class PeptideMutant
{
  private Boolean mIsMutated=null;
  private Range<Integer> mMutated=null;
  private String mSequenceWT, mSequenceMT, mHeader, mGene, mTranscriptAccession, mMutationType, mMutation;

  public PeptideMutant() { super(); }
  public PeptideMutant(String header, String sequence)
  {
    // >MT.A1BG.ENST00000263100.missense.52H/R
    // IFYETQPSLWAESESLLKPLANVTLTCQARLETPDFQLFKNGVAQEPVHLDSPAIKHQF
    mHeader=header; init(header, sequence);
  }

  public PeptideMutant init(String s, String seq)
  {
    // MT.A1BG.ENST00000263100.missense.52H/R
    // WT.FOXD2.ENST00000334793.inframe_ins.287-288-/PH
    // MT.A3GALT2.ENST00000330379.FS.50-51TAGTC/T
    // MT.AC006116.20.ENST00000586573.inframe_del.54L/-
    // MT.AGAP4.ENST00000448048.inframe_del.36-43AGDRMAGA/A

    String[] items = Strs.split(s, '.', true);
    if (items.length>0)
    {
      if      ("MT".equals(items[0])) { mSequenceMT=seq; mIsMutated=true; }
      else if ("WT".equals(items[0])) { mSequenceWT=seq; mIsMutated=false; }
    }
    items = s.substring(s.indexOf('.')+1).split(".ENST");
    if (items.length>1)
    {
      mGene=items[0];
      items = Strs.split("ENST"+items[1], '.', true);
      if (items.length>0) mTranscriptAccession=items[0];
      if (items.length>1) mMutationType       =items[1];
      if (items.length>2) mMutation           =items[2];
    }

    return this;
  }
  public Range<Integer> getMutationRange()
  {
    if (mMutated==null && mSequenceMT!=null || mSequenceWT!=null)
    {
      // SampleDerived$ grep ">MT" SCVhIG_84LN.pep.dedup.fasta | grep -v missense | grep -v inframe_ins | grep -v .FS. | grep -v inframe_del | more
      int start=-1, stop=-1;
      if (Strs.isA(mMutationType, "missense"))
      {
        for (int i=0; i<mSequenceMT.length(); i++)
        {
          if (i>=mSequenceWT.length() || mSequenceWT.charAt(i)!=mSequenceMT.charAt(i))
          {
            if (start==-1) start=i;
            stop=i;
          }
          else if (start>=0) { break; }
        }
      }
      else if (Strs.isA(mMutationType, "inframe_ins"))
      {
        //>MT.AC004076.9.ENST00000596831.inframe_ins.4A/AA
        //MAAAAEPMGPAQVPMNSEVIVDPIQGQVNFEDVFVYFSQEEWVLLDEAQRLLYRDVMLEN
        //MAAA EPMGPAQVPMNSEVIVDPIQGQVNFEDVFVYFSQEEWVLLDEAQRLLYRDVMLEN
        int offset=0;
        for (int i=0; i<mSequenceMT.length(); i++)
        {
          if (i+offset<mSequenceWT.length() && mSequenceWT.charAt(i+offset)!=mSequenceMT.charAt(i))
          {
            if (start==-1) start=i;
            stop=i; offset=-1;
          }
          else if (start>=0) { break; }
        }
      }
      else if (Strs.isA(mMutationType, "inframe_del"))
      {
        //>MT.AC006116.20.ENST00000586573.inframe_del.54L/-
        //LPSGIIGLMSRLSPDDGLNPNRCSCCIYAPPEPPHLPFFKWTYSFL KSINLKKLLYTA
        //LPSGIIGLMSRLSPDDGLNPNRCSCCIYAPPEPPHLPFFKWTYSFLLKSINLKKLLYTA
        int offset=0;
        for (int i=0; i<mSequenceMT.length(); i++)
        {
          if (i+offset<mSequenceWT.length() && mSequenceWT.charAt(i+offset)!=mSequenceMT.charAt(i))
          {
            if (start==-1) start=i;
            stop=i; offset=1;
          }
          else if (start>=0) { break; }
        }

      }
      else if (Strs.isA(mMutationType, "FS"))
      {
        //>MT.AC010536.1.ENST00000538868.FS.138C/CA
        //LTTEPFPGAPFCSDSCSPALCCVLLTPGRGERERGPSGKLPGRQGWVLPYRYSTCCTGNG
        //LTTEPFPGAPFCSDSCSPALCCVLLTPGRG
        //                              RGSVAPLGSSQDGRGGCCPTDTVPAAQG

        //>MT.AKR7L.ENST00000429712.FS.23-30CTGCGGCGCTGGTGGGCGCGTCCA/C
        //MSRQLSRARPATVLGAMEMGRR
        //                      SHARLPGARPHRDRHGLPVQRRPVRDHPWRPGAPNGQQRLQSENCYQGQSMDWELPEA
        //MSRQLSRARPATVLGAMEMGRRMDAPTSAAVTRAFLERGHTEIDTAFLYSD
        for (int i=0; i<mSequenceMT.length(); i++)
        {
          if (mSequenceWT.charAt(i)!=mSequenceMT.charAt(i))
          {
            if (start==-1) start=i;
            stop=mSequenceMT.length()-1;
            break;
          }
        }
      }
      else
      {
        System.out.println();
      }
      mMutated = Range.closed(start, stop);
    }
    return mMutated;
  }

  public Boolean isMutated()              { return mIsMutated; }
  public String  getGene()                { return mGene; }
  public String  getHeading()             { return mTranscriptAccession+"."+mMutationType+"."+mMutation; }
  public String  getTranscriptAccession() { return mTranscriptAccession; }
  public String  getMutationType()        { return mMutationType; }
  public String  getMutation()            { return mMutation; }
  public String  getMutatedSequence()     { return mSequenceMT; }
  public String  getWildTypeSequence()    { return mSequenceWT; }

  public PeptideMutant setMutatedSequence( String s) { mSequenceMT=s; return this; }
  public PeptideMutant setWildTypeSequence(String s) { mSequenceWT=s; return this; }
}
