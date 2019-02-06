package org.ms2ms.test.ms;

import com.google.common.collect.Multimap;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.junit.Test;
import org.ms2ms.algo.LCMSMS;
import org.ms2ms.io.PsmReaders;
import org.ms2ms.r.Dataframe;
import org.ms2ms.test.TestAbstract;

/**
 * Created by yuw on 10/3/2015.
 */
public class SearchOutputTest extends TestAbstract
{
  @Test
  public void RocTest() throws Exception
  {
    String root = "c:/local/data/TMT_MS3/";
    double[] scale = new double[] {0.95,0.8,0.64, 0.32,0.16,0.08,0.04,0.02,0.01,0.005,0.0025,0.00125,0.00075,0.00025};
    Multimap<SpectrumIdentifier, PeptideMatch>
        ms2ce40 = PsmReaders.readMSGFplus(root+"TMT_MS3120_1800_40CE.ms2.tsv", '\t', 0),
        ms3ce40 = PsmReaders.readMSGFplus(root+"TMT_MS3120_1800_40CE.ms3.tsv", '\t', 0);
//    Multimap<SpectrumIdentifier, PeptideMatch> psm = PsmReaders.readMSGFplus(root+"ms2_only/CATSIgGTMTMS3.tsv");

    Dataframe cuts = LCMSMS.cut(ms2ce40, "QValue", scale);

    cuts.init(true);
  }
}
