package org.ms2ms.mzjava;

import com.google.common.collect.Multimap;
import org.expasy.mzjava.proteomics.io.ms.ident.MaxQuantPsmReader;
import org.expasy.mzjava.proteomics.mol.modification.ModificationResolver;
import org.expasy.mzjava.proteomics.ms.ident.PeptideMatch;
import org.expasy.mzjava.proteomics.ms.ident.SpectrumIdentifier;
import org.ms2ms.data.ms.MaxQuant;
import org.ms2ms.data.ms.PSM;
import org.ms2ms.math.Stats;
import org.ms2ms.r.Dataframe;
import org.ms2ms.r.Var;
import org.ms2ms.utils.IOs;
import org.ms2ms.utils.Strs;
import org.ms2ms.utils.TabFile;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;

/**
 * Created by yuw on 10/28/2015.
 */
public class MaxQuantReader extends MaxQuantPsmReader
{
  public MaxQuantReader()                                 { super(); }
  public MaxQuantReader(ModificationResolver modResolver) { super(modResolver); }

  @Override
  protected SpectrumIdentifier getSpectrumId(ResultSet results) throws SQLException
  {
    SpectrumIdentifier identifier = super.getSpectrumId(results);

    try
    {
      identifier.setPrecursorIntensity(Stats.toDouble(Strs.split(results.getString("Intensities"), ';')[0]));
      String[] scans = Strs.split(identifier.getSpectrum(), '.');
      identifier.setName(scans[0]+"#"+scans[1]);
      identifier.setSpectrumFile(scans[0]);
    }
    catch (Exception e) {}

    return identifier;
  }
  @Override
  protected void setValues(PeptideMatch peptideMatch, ResultSet results) throws SQLException
  {
    // set the attributes
    try { peptideMatch.setMassDiff(results.getDouble("Mass Error [ppm]") * results.getDouble("Mass")*1E-6); } catch (Exception e) {}
    try { peptideMatch.setTotalNumIons(         results.getInt("Number of Matches"));     } catch (Exception e) {}
    try { peptideMatch.addScore("Score",        results.getDouble("Score"));              } catch (Exception e) {}
    try { peptideMatch.addScore("dScore",       results.getDouble("Delta score"));        } catch (Exception e) {}
    try { peptideMatch.addScore("PEP",          results.getDouble("PEP"));                } catch (Exception e) {}
    try { peptideMatch.addScore("IntCoverage",  results.getDouble("Intensity coverage")); } catch (Exception e) {}
    try { peptideMatch.addScore("PeakCoverage", results.getDouble("Peak coverage"));      } catch (Exception e) {}
  }
}
