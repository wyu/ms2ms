package org.ms2ms.test.io;

import org.junit.Test;
import org.ms2ms.math.Stats;
import org.ms2ms.test.TestAbstract;
import org.ms2ms.utils.Strs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;

/**
 * Created by yuw on 1/23/16.
 */
public class CruxIO extends TestAbstract
{
  @Test
  public void writeBulleyes() throws Exception
  {
    parseBulleyesHeaders("/Users/yuw/Documents/Apps/QC/data/160113_Clare_CHO_direct_10.bull/");
    parseBulleyesHeaders("/Users/yuw/Documents/Apps/QC/data/160113_Clare_OPC_direct_10.bull/");
    parseBulleyesHeaders("/Users/yuw/Documents/Apps/QC/data/160113_Clare_microglia_direct_10.bull/");
  }

  public static void parseBulleyesHeaders(String root) throws Exception
  {
    System.out.println("Reading from " + root);
//    String root = "/Users/yuw/Documents/Apps/QC/data/160113_Clare_CHO_direct_10.bull/";
    BufferedReader nopid = new BufferedReader(new FileReader(root+"bullseye.no-pid.ms2")),
                     pid = new BufferedReader(new FileReader(root+"bullseye.pid.ms2"));

    FileWriter tsv = new FileWriter(root+"bulleyes_summary.tsv");
    // write the header
    tsv.write("pid\tscan1\tscan2\tmz\tRT\tz\tintensity\tTIC\tions\n");

    String scan1=null, scan2=null, mz=null, RT=null;
    double tic=0; int ions=0;
    //  bullseye.no-pid. – a file containing the fragmentation spectra for which accurate masses were not inferred.
    //
    //  S       30      30      1188.9492
    //  I       RTime   25.2458
    //  Z       2       2376.8912
    System.out.println("Parsing no-pid.ms2 ...");
    while (nopid.ready())
    {
      String[] line = Strs.split(nopid.readLine(), '\t'), nums=Strs.split(line[0], ' ');
      if (Strs.equals("S", line[0]))
      {
        if (tic>0) { tsv.write(tic+"\t"+ions+"\n"); tic=0; ions=0; }
        scan1=line[1]; scan2=line[2]; mz=line[3];
      }
      else if (Strs.equals("I", line[0]))
      {
        RT=line[2];
      }
      else if (Strs.equals("Z", line[0]))
      {
        tsv.write("no\t"+scan1+"\t"+scan2+"\t"+mz+"\t"+RT+"\t"+line[1]+"\t\t");
      }
      else if (line.length==1 && Stats.toDouble(nums[1])!=null) { tic+=Stats.toDouble(nums[1]); ions++; }
    }
    nopid.close();

    //  bullseye.pid. – a file containing the fragmentation spectra for which accurate masses were successfully inferred.
    //  Unless otherwise specified (with the --spectrum-format option), the output file format is ".ms2".
    //  Note that if the output format is ".ms2," then a single spectrum may have multiple "Z" lines, each indicating
    //  a charge state and accurate mass. In addition, Bullseye inserts an "I" line (for charge-dependent analysis)
    //  corresponding to each "Z" line. The "I" line contains "EZ" in the second column, the charge and mass from the
    //  associated "Z" line in the third and fourth colummns, followed by the chromatographic apex and the intensity at the chromatographic apex.
//    S       1800    1800    623.9537
//    I       RTime   38.8216
//    I       EZ      1       620.9632        38.4822 1406121.6
//    I       EZ      1       621.0686        38.4666 1150576.4
//    Z       1       620.9632
//    Z       1       621.0686
    System.out.println("Parsing pid.ms2 ...");
    tic=0; ions=0;
    while (pid.ready())
    {
      String[] line = Strs.split(pid.readLine(), '\t'), nums=Strs.split(line[0], ' ');
      if (Strs.equals("S", line[0]))
      {
        if (tic>0) { tsv.write(tic+"\t"+ions+"\n"); tic=0; ions=0; }
        scan1=line[1]; scan2=line[2]; mz=line[3];
      }
      else if (Strs.equals("I", line[0]) && Strs.equals("EZ",line[1]))
      {
        tsv.write("yes\t"+scan1+"\t"+scan2+"\t"+mz+"\t"+line[4]+"\t"+line[2]+"\t"+line[5]+"\n");
      }
      else if (line.length==1 && Stats.toDouble(nums[1])!=null) { tic+=Stats.toDouble(nums[1]); ions++; }
    }
    pid.close(); tsv.close();
  }
}
