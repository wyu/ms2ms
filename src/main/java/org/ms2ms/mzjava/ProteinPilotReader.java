package org.ms2ms.mzjava;

import com.google.common.collect.Lists;
import org.expasy.mzjava.proteomics.io.ms.ident.ProteinPilotPsmReader;
import org.expasy.mzjava.proteomics.mol.modification.ModificationResolver;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.math.Stats;
import org.ms2ms.utils.Strs;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collection;
import java.util.List;

/**
 * Created by yuw on 10/28/2015.
 */
public class ProteinPilotReader extends ProteinPilotPsmReader
{
  public ProteinPilotReader()                                 { super(); }
  public ProteinPilotReader(ModificationResolver modResolver) { super(modResolver); }

  @Override
  protected SpectrumIdentifier getSpectrumId(ResultSet results) throws SQLException
  {
    int charge = results.getInt("Theor z");
    double precMw = results.getDouble("Theor MW");
    String spectrum = results.getString("Spectrum");
//    double rt = results.getDouble("Time");

    // let's parse out the scan number: "1.1.1.6275.1"
    SpectrumIdentifier identifier = new SpectrumIdentifier(spectrum);
    identifier.setAssumedCharge(charge);
    identifier.setPrecursorNeutralMass(precMw);
    identifier.addScanNumber(Stats.toInt(Strs.split(spectrum, '.')[3]));
//    identifier.addRetentionTime(rt, TimeUnit.MINUTE);

    return identifier;
  }
  @Override
  protected Collection<String> getAccessionNumbers(ResultSet results) throws SQLException
  {
    return parseAccessionNumbers(results.getString("Accessions"));
  }
  // correct a mistake in the base class
  private List<String> parseAccessionNumbers(String headers)
  {
    List<String> accessions = Lists.newArrayList();
    if (headers != null && headers.length() > 0)
      for (String header : headers.split(";"))
        accessions.add(header.contains("\\|")?header.split("\\|")[1]:header);

    return accessions;
  }

}
