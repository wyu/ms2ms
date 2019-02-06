package org.ms2ms.algo.SuperQuant;

/** Adopted from SuperQuantNode on GitHub
 *
 * Created by yuw on 11/5/2015.
 */
public class ComplementaryFinder
{
/*
  private static int PackageSize = 1000;
  //private MassSpectrumCollection spectrumBuffer = new MassSpectrumCollection(PackageSize); //send buffer
  private Collection<MsnSpectrum> spectrumBuffer = new ArrayList<>(); //send buffer

  static double Proton = 1.00727646677;

//  private delegate MassCentroid PeakIntensifier(MassCentroid peak, double rank);
//  private PeakIntensifier intensify = null;
//  private MultiRange exclusionList = null;

  private int spectra_in = 0;
  private int spectra_out = 0;
  private int magic = 0;

*/
/*
  [IntSelectionParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Secondary peptide charge",
  Description = "Possible charge states of secondary peptides",
  SelectionValues = new int[] { 1, 2, 3, 4, 5, 6, 7, 8 },
  DefaultValue = "2, 3, 4",
  IsMultiSelect = true,
  ValueRequired = true,
  Position = 1)]

  public SimpleSelectionParameter<int> chargeStates;
*//*

  public int[] chargeStates;

*/
/*
  [IntegerParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Min number of peaks",
  Description = "Minimum number of peaks per secondary spectrum to extract it",
  MinimumValue = "2",
  DefaultValue = "6",
  ValueRequired = true,
  Position = 2)]

  public IntegerParameter minPeaks;
*//*

  public int minPeaks;

*/
/*
  [MassToleranceParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Parent mass tolerance",
  Description = "Mass tolerance to consider parent masses the same",
  MinimumValue = "0.0 ppm",
  MaximumValue = "100.0 ppm",
  DefaultValue = "5.0 ppm",
  Subset = "ppm",
  ValueRequired = true,
  Position = 3)]

  public MassToleranceParameter mergePPM;
*//*

  public PpmTolerance mergePPM;

*/
/*
  [BooleanParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Relative borders",
  Description = "If True the lower and upper borders of allowed coisolation window (below) are set relative to the instrument isolation window, absolute value is used otherwise",
  ValueRequired = true,
  DefaultValue = "True",
  IsAdvanced = false,
  Position = 4)]

  public BooleanParameter relativeBorders;
*//*

  public boolean relativeBorders;

*/
/*
  [DoubleParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Lower border of coisolation window (u)",
  Description = @"Lower border of the mass window to consider coisolated peptides.
    Positive values move the border to lower mass values (increase isolation window). Negative values should be used only with relative borders.
  If not set the isolation window saved in MS spectrum is used",
  MinimumValue = "-10.0",
  MaximumValue = "10.0",
  DefaultValue = "0.6",
  ValueRequired = false,
  IsAdvanced = false,
  Position = 5)]

  public DoubleParameter lowerBorder;
*//*

  public double lowerBorder;

*/
/*
  [DoubleParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Upper border of coisolation window (u)",
  Description = @"Upper border of the mass window to consider coisolated peptides.
    Positive values move the border to higher mass values (increase isolation window). Negative values should be used only with relative borders.
  If not set the isolation window saved in MS spectrum is used",
  MinimumValue = "-10.0",
  MaximumValue = "10.0",
  ValueRequired = false,
  IsAdvanced = false,
  Position = 6)]

  public DoubleParameter upperBorder;
*//*

  public double upperBorder;

*/
/*
  [StringSelectionParameter(
  Category = "1. Spectrum extraction",
  DisplayName = "Spectra to return",
  Description = "Which spectra should be returned for further processing",
  SelectionValues = new string[] { "Primary and Secondary", "Primary Only", "Secondary Only", "Intensify Only" },
  DefaultValue = "Primary and Secondary",
  IsMultiSelect = false,
  ValueRequired = true,
  IsAdvanced = true,
  Position = 7)]

  public SimpleSelectionParameter<string> returnSelection;
*//*

  public String returnSelection;

*/
/*
  [BooleanParameter(
  Category = "4. Secondary spectra mass verification",
  DisplayName = "Enable MS1 verification",
  Description = "If enabled checks if the predicted mass peak for secondary spectra is present in parent MS1 spectrum and removes spectra that don't meet this criteria",
  DefaultValue = "False",
  ValueRequired = true,
  IsAdvanced = true,
  Position = 1)]

  public BooleanParameter verifyMS1;
*//*

  public boolean verifyMS1;

*/
/*
  [MassToleranceParameter(
  Category = "4. Secondary spectra mass verification",
  DisplayName = "Mass tolerance",
  Description = "Mass tolerance used for MS1 verification",
  MinimumValue = "0.0 mmu| 0.0 ppm",
  MaximumValue = "100.0 mmu| 100.0 ppm",
  DefaultValue = "5.0 ppm",
  Subset = "mmu|ppm",
  ValueRequired = true,
  IsAdvanced = true,
  Position = 2)]

  public MassToleranceParameter VerifyPPM;
*//*

  public PpmTolerance VerifyPPM;

*/
/*
  [BooleanParameter(
  Category = "2. Peak intensification",
  DisplayName = "Enable intensification",
  Description = "Set true if peak intensification is necessary",
  DefaultValue = "True",
  ValueRequired = true,
  Position = 1)]

  public BooleanParameter doIntensification;
*//*

  public boolean doIntensification;

*/
/*
  [IntegerParameter(
  Category = "2. Peak intensification",
  DisplayName = "Intensification factor: b - y",
  Description = "Intensifiaction factor for b - y ion pairs",
  DefaultValue = "1",
  ValueRequired = true,
  Position = 2)]

  public IntegerParameter byRank;
*//*

  public int byRank;

*/
/*
  [StringParameter(
  Category = "3. Peak exclusion",
  DisplayName = "Exclusion list",
  Description = "Path to the exclusion list (the default one contains all immonium ions > 50 Da)",
  DefaultValue = ".\\cf.exclude",
  ValueRequired = true,
  Position = 1)]

  public StringParameter ignoreListPath;
*//*

  public String ignoreListPath;

*/
/*
  [MassToleranceParameter(
  Category = "3. Peak exclusion",
  DisplayName = "Exclusion tolerance",
  Description = "Exclusion window is the doubled tolerance, e.g. (-5 .. +5) ppm",
  MinimumValue = "0.0 mmu| 0.0 ppm",
  MaximumValue = "100.0 mmu| 100.0 ppm",
  DefaultValue = "10.0 ppm",
  Subset = "mmu|ppm",
  ValueRequired = true,
  Position = 2)]

  public MassToleranceParameter ignoreListPPM;
*//*

  public PpmTolerance ignoreListPPM;

  */
/* Change intensity of peak according to rank *//*

  private Peak validIntensify(Peak peak, double maxIntensity)
  {
    Peak newPeak = peak.copy();
    newPeak.setIntensity(peak.getIntensity() + this.byRank.Value * maxIntensity);

    return newPeak;
  }

  */
/* Void intensification *//*

  private Peak voidIntensify(Peak peak, double maxIntensity) { return peak; }
  */
/* Translate index of MultiRange to charge state *//*

  private short getCharge (int index) { return (short) this.chargeStates.Values[index / 2]; }

*/
/* contents of the cf.exclude file
    60.044938
    70.065674
    72.081324
    74.060588
    76.022094
    86.096974
    87.055837
    88.039853
    101.071487
    101.107873
    102.055503
    104.053395
    110.071822
    120.081324
    129.114021
    133.043558
    136.076238
    159.092223
*//*

  private MultiRange loadExcludes ()
  {
    MultiRange exclusionList = new MultiRange();

    BufferedReader reader=null;
    try
    {
      try
      {
        reader = new BufferedReader(new FileReader(ignoreListPath));
        while (reader.ready())
        {
          double mass = Stats.toDouble(reader.readLine());
          exclusionList.add(ignoreListPPM.getMin(mass), ignoreListPPM.getMax(mass));
        }
      }
      finally
      {
        if (reader!=null) reader.close();
      }
    }
    catch (IOException ie) {}

    return exclusionList;
  }

  */
/* Search for peak pairs that reside inside MultiRange and annotate b-y pairs in the spectrum with unique parent mass ID
   * Return: MassCentroidCollection of all parent masses *//*

  private List<AnnotatedPeak> annotateSpectrum(MsnSpectrum spectrum, MultiRange range)
  {
    int start = 0;                   //lower seeker
    int end   = spectrum.size() - 1; //upper seeker
    int k;                           //intermediate seeker
    int lowIndex;                    //position in multirange
    double mass;
//    int spectrumID;                  //unique ID for each subspectrum (set of ion pairs) in spectrum

    Table<Integer, Integer, Boolean> index_prop = HashBasedTable.create();

//    MassAggregator masses = new MassAggregator(this.mergePPM.Value.Tolerance); //collection of possible masses
    List<AnnotatedPeak> masses = new LinkedList<>();

    AnnotatedPeak instrumentMass = new AnnotatedPeak(
      (spectrum.getPrecursor().getMz()-Proton)*spectrum.getPrecursor().getCharge(),
       spectrum.getPrecursor().getIntensity(),
       spectrum.getPrecursor().getCharge(), 0);

    //always add initial mass
    masses.add(instrumentMass);
    //continue untill two seekers are met
    while (start < end)
    {
      //matrix_sum of two peaks converted to neutral M
      mass = spectrum.getMz(start)+spectrum.getMz(end)-2*Proton;

      //mass is too big, make smaller = move end backwards
      if (mass>range.getHighest())        end--;
      //mass is too small, make bigger = move start forward
      else if (mass < range.getLowest())  start++;
      else
      {
        //check all possible pairs with masses higer or equal than current
        k=0; lowIndex=0;
        do
        {
          lowIndex = range.lastIndexBelow(lowIndex, mass); //find the highest index in range smaller than mass
          if (lowIndex % 2 == 0) //if *lowIndex* is even (it is the begginning of subrange) add peak
          {
            masses.add(new AnnotatedPeak(mass, spectrum.getIntensity(start + k) + spectrum.getIntensity(end), getCharge(lowIndex), 2));
            // create the tags
            index_prop.put(start+k, masses.size()-1, true);
            index_prop.put(end, masses.size() - 1, true);
//            spectrumID = masses.AddCentroid(
//              new MassCentroid(mass, spectrum.PeakCentroids[start + k].Intensity + spectrum.PeakCentroids[end].Intensity,
//                getCharge(lowIndex), spectrum.Precursor.Resolution, 2));
            //signal-to-noise is the number of peaks supporting this mass
//            spectrum.PeakCentroids[start + k].SetAnnotation(spectrumID.ToString(), true);
//            spectrum.PeakCentroids[end].SetAnnotation(spectrumID.ToString(), true);
            k++; //try next peak
          }
          //find the index of the first peak above the beginning of next subrange
          else
          {
            for (int i=start+k; i<spectrum.size(); i++)
            {
              if (spectrum.getMz(i)>spectrum.getMz(end) || spectrum.getMz(i)>range.get(lowIndex+1) - spectrum.getMz(end)+2*Proton)
              {
                k=i-start; break;
              }
            }
//            k = spectrum.PeakCentroids.FindIndex(start + k, //start from current peak
//              x => x.Position > spectrum.PeakCentroids[end].Position //peaks above [end] were checked earlier
//              || x.Position > range[lowIndex + 1] - spectrum.PeakCentroids[end].Position + 2 * Proton) - start;
          }

          if (start + k >= end || k < 0) break; // no such peak

          mass = spectrum.getMz(start+k)+spectrum.getMz(end) - 2 * Proton;

        } while (mass < range.getHighest() && k < spectrum.size() - 1);
        //do so untill mass exeeds the upper border or there is no peaks in the spectrum

        end--; // make initial mass smaller, hiher masses were checked earlier
      }
    }

    return masses;
  }

  */
/* Extract spectra according to annotation. Only parent masses supported by at least minPeaks peaks are used,
   * parentMasses is a collection of masses
   * Return: Collection of spectra, position 0 - primary, other - secondary
   *//*

  private Collection<MsnSpectrum> extractSpectra(MsnSpectrum spectrum, Collection<Peak> parentMasses, int minPeaks)
  {
    Collection<MsnSpectrum> outSpectra = new ArrayList<>();

    // always put the original one at first
    outSpectra.add(spectrum.copy(new IdentityPeakProcessor<>())); //clone initial spectrum to position 0; all unassigned peaks will be here

    // assume each parent masses are already validated (diff from C# implementation)
    for (Peak parent : parentMasses)
    {
      MsnSpectrum additionalSpectrum = spectrum.copy(new IdentityPeakProcessor<PeakAnnotation>());

      additionalSpectrum.clear(); //remove all peaks
      additionalSpectrum.clearAnnotations();;
      additionalSpectrum.setPrecursor(parent); //correct mass
      outSpectra.add(additionalSpectrum);
    }
    for (MassCentroid peak in spectrum.PeakCentroids) //try each Centroid in initial spectrum
    {
      if (peak.HasAnnotations) //if it was annotated at all (most of peaks aren't annotated)
      {
        for (int i = 0; i < validIDs.Count; i++) //check if it was annotated with valid ID
        {
          if (peak.GetAnnotation(validIDs[i].ToString()) != null)
          {
            outSpectra[i + 1].PeakCentroids.Add(intensify(peak, spectrum.Header.BasePeakIntensity)); //add it to the corresponding outSpectra element
            outSpectra[0].PeakCentroids.Remove(peak); //remove from initial
          }
        }
      }
    }

    for (int i = 1; i < outSpectra.Count; i++) //add unassigned peaks to each spectrum
    {
      //Log.DebugFormat("Extraction: Subspectrum {0} - {1} peaks", i, outSpectra[i].PeakCentroids.Count);
      outSpectra[i].PeakCentroids.AddRange(outSpectra[0].PeakCentroids.AsEnumerable<MassCentroid>());
      outSpectra[i].PeakCentroids.Sort();
    }

    outSpectra.RemoveAt(0); //remove the initial spectrum

    List<Integer>             validIDs = new ArrayList<>();

    //always include the initial mass (index 0)
    validIDs.add(0);

    if (verifyMS1.Value) //verify parent masses according to MS1 if necessary
    {
      foreach (int i in areMassesInSpectrum(parentMasses, spectrum)) //check IDs, supported by MS1
      {
        if (parentMasses[i].SignalToNoise >= minPeaks) validIDs.Add(i); //signal-to-noise is used to store number of ions in a group
      }
    }
    else
    {
      for (int i = 1; i < parentMasses.Count; i++) //check all IDs
      {
        if (parentMasses[i].SignalToNoise >= minPeaks) validIDs.Add(i);
      }
    }

    outSpectra.add(spectrum.copy(new IdentityPeakProcessor<>())); //clone initial spectrum to position 0; all unassigned peaks will be here

    if (validIDs.Count > 0) //if at least one valid spectrum was found
    {
      for (int i = 0; i < validIDs.Count; i++)//add one clone for each additional spectrum
      {
        MassSpectrum additionalSpectrum = spectrum.Clone();

        additionalSpectrum.ProfilePoints.Clear(); //remove all peaks
        additionalSpectrum.PeakCentroids.Clear();
        additionalSpectrum.Precursor.Charge = parentMasses[validIDs[i]].Charge; //correct mass
        additionalSpectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0] = parentMasses[validIDs[i]];
        additionalSpectrum.Precursor.SinglyChargedMass = new MassValue(parentMasses[validIDs[i]].Position + Proton);
        additionalSpectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0].Position = (parentMasses[validIDs[i]].Position //calculate m/z
          + parentMasses[validIDs[i]].Charge * Proton) / parentMasses[validIDs[i]].Charge;
        outSpectra.Add(additionalSpectrum);
      }

      foreach (MassCentroid peak in spectrum.PeakCentroids) //try each Centroid in initial spectrum
      {
        if (peak.HasAnnotations) //if it was annotated at all (most of peaks aren't annotated)
        {
          for (int i = 0; i < validIDs.Count; i++) //check if it was annotated with valid ID
          {
            if (peak.GetAnnotation(validIDs[i].ToString()) != null)
            {
              outSpectra[i + 1].PeakCentroids.Add(intensify(peak, spectrum.Header.BasePeakIntensity)); //add it to the corresponding outSpectra element
              outSpectra[0].PeakCentroids.Remove(peak); //remove from initial
            }
          }
        }
      }

      for (int i = 1; i < outSpectra.Count; i++) //add unassigned peaks to each spectrum
      {
        //Log.DebugFormat("Extraction: Subspectrum {0} - {1} peaks", i, outSpectra[i].PeakCentroids.Count);
        outSpectra[i].PeakCentroids.AddRange(outSpectra[0].PeakCentroids.AsEnumerable<MassCentroid>());
        outSpectra[i].PeakCentroids.Sort();
      }

      outSpectra.RemoveAt(0); //remove the initial spectrum
    }

    return outSpectra;
  }

  private MassSpectrumCollection intensifyOnly(MassSpectrum spectrum)
            */
/* Do intensification of peaks corresponding to targeted mass only, all secondary masses are ignored
             * spectrum - annotated spectrum
             * Return: MassSpectrumCollection with one spectrum
             *//*

  {
    MassSpectrumCollection outSpectra = new MassSpectrumCollection(1);
    MassSpectrum unintSpectrum = spectrum.Clone(); //position 0 is for unintensified peaks

    MassSpectrum intSpectrum = spectrum.Clone(); //create new entity of spectrum
    intSpectrum.ProfilePoints.Clear(); //remove all peaks
    intSpectrum.PeakCentroids.Clear();

    foreach (MassCentroid peak in spectrum.PeakCentroids) //try each Centroid in initial spectrum
    {
      if (peak.HasAnnotations) //if it was annotated at all (most of peaks aren't annotated)
      {
        if (peak.GetAnnotation("0") != null) //if the peak is annotated to the targeted mass
        {
          intSpectrum.PeakCentroids.Add(intensify(peak, spectrum.Header.BasePeakIntensity)); //add intensified peak
          unintSpectrum.PeakCentroids.Remove(peak); //remove from unintensified peaks
        }
      }
    }

    intSpectrum.PeakCentroids.AddRange(unintSpectrum.PeakCentroids.AsEnumerable<MassCentroid>()); //add all unintensified peaks
    intSpectrum.PeakCentroids.Sort();

    outSpectra.Add(intSpectrum);//add spectrum to the output
    return outSpectra;
  }

  private MassCentroidCollection applyExclusionList (MassSpectrum spectrum)
            */
/* Process the spectrum (IT WILL BE MODIFIED!) by removing all peaks that are inside exclusion list
             * collect all filtered peaks in a separate collection and return it
             *//*

  {
    MassCentroidCollection filtered = new MassCentroidCollection(); //filtered peaks collection

    for (int i = 0; i < exclusionList.Count; i += 2) //check all ranges in exclusion list
    {
      foreach (MassCentroid peak in spectrum.PeakCentroids.FindAllPeaksWithinMassRange(exclusionList[i], exclusionList[i + 1]))
      //find all peaks inside the range; remove them from the spectrum and save in collection
      {
        spectrum.PeakCentroids.Remove(peak);
        filtered.Add(peak);
      }
    }

    return filtered;
  }

  private MassSpectrumCollection rejoinSpectra (MassSpectrumCollection inSpectra, MassCentroidCollection peaksToJoin)
            */
/* Add peaks from MassCentroidCollection to every spectrum in MassSpectrumCollection
             * return modified MassSpectrumCollection
             *//*

  {
    foreach (MassSpectrum spectrum in inSpectra)
    {
      spectrum.PeakCentroids.AddRange(peaksToJoin.AsEnumerable());
      spectrum.PeakCentroids.Sort();
    }

    return inSpectra;
  }

  private List<int> areMassesInSpectrum(MassCentroidCollection masses, MassSpectrum spectrum)
            */
/* For each mass in mass list search if the peak is present in the spectrum
             * Collect and return indices of valid peaks *//*

  {
    double tol; //tolerance in amu
    double mz, cmz; //predicted m/z and the closest match from MS1
    List<int> valid = new List<int>();

    for (int c = 1; c < masses.Count; c++) //check each mass except the first one (first mass - primary spectrum)
    {
      mz = (masses[c].Position + Proton * masses[c].Charge) / masses[c].Charge; //calculate m/z value to search and mass tolerance
      tol = this.VerifyPPM.Value.GetToleranceInU(mz);

      //check if there is at least one peak in precursor spectrum and find the closest peak to m/z value
      if (spectrum.Precursor.IsotopeClusterPeakCentroids.Count > 0) cmz = spectrum.Precursor.IsotopeClusterPeakCentroids.FindClosestPeak(mz).Position;
      else cmz = 0;

      if (Math.Abs(cmz - mz) <= tol) valid.Add(c); //add valid index
    }

    return valid;
  }

  private void ProcessBatch(MassSpectrumCollection results)
        */
/* Process a set of mass spectra using complementary finder *//*

  {
    foreach(MassSpectrum iSpec in results)
    {
      double low, high;
      MultiRange range = new MultiRange();
      MassSpectrum spectrum = iSpec.Clone(); //make local copy of the input spectrum

      if (!this.lowerBorder.IsValueSet) low = spectrum.ScanEvent.IsolationMass - spectrum.ScanEvent.IsolationWindow.LowerLimit;//get right value of lower border
      else if(this.relativeBorders.Value) low = spectrum.ScanEvent.IsolationMass - spectrum.ScanEvent.IsolationWindow.LowerLimit + this.lowerBorder.Value;
      else low =  this.lowerBorder.Value;

      if (!this.upperBorder.IsValueSet) high = spectrum.ScanEvent.IsolationWindow.UpperLimit - spectrum.ScanEvent.IsolationMass;//get right value of upper border
      else if (this.relativeBorders.Value) high = spectrum.ScanEvent.IsolationWindow.UpperLimit - spectrum.ScanEvent.IsolationMass + this.upperBorder.Value;
      else high =  this.upperBorder.Value;

      //mgf-file input compatability
      if (spectrum.ScanEvent.IsolationMass == 0.0) spectrum.ScanEvent.IsolationMass = spectrum.Precursor.MeasuredMonoisotopicPeakCentroids[0].Position;
      spectrum.PeakCentroids.Sort();

      foreach (int z in this.chargeStates.Values)
      {
        range.Add((spectrum.ScanEvent.IsolationMass - low - Proton) * z, (spectrum.ScanEvent.IsolationMass + high - Proton) * z); //neutral MolMass
      }

      if (range.Count == 0)//check if parameters are valid
      {
        SendAndLogErrorMessage("The window for possible parent masses is empty! Please, check the settings of allowed coisolation window", true);
        throw new Exception("Allowed coisolation window is empty");
      }

      //for debugging purposes
      if (spectrum.Header.ScanNumbers[0] == magic)
      {
        Log.DebugFormat("Magic!");
      }

      MassCentroidCollection excludedPeaks = applyExclusionList(spectrum); //apply exclusion list and put all excluded peaks in separate collection

      MassCentroidCollection parentMasses = annotateSpectrum(spectrum, range); //annotate complementary pairs

      MassSpectrumCollection extractedSpectra = new MassSpectrumCollection();

      switch (this.returnSelection.Value) //filter spectra acording to return selection
      {
        case "Primary and Secondary":
          extractedSpectra = rejoinSpectra(extractSpectra(spectrum, parentMasses, this.minPeaks.Value), excludedPeaks); //extract spectra
          break;
        case "Primary Only":
          extractedSpectra = rejoinSpectra(extractSpectra(spectrum, parentMasses, this.minPeaks.Value), excludedPeaks); //extract spectra
          extractedSpectra.RemoveRange(1, extractedSpectra.Count - 1);
          break; //remove all secondary
        case "Secondary Only":
          extractedSpectra = rejoinSpectra(extractSpectra(spectrum, parentMasses, this.minPeaks.Value), excludedPeaks); //extract spectra
          extractedSpectra.RemoveAt(0);
          break; //remove primary
        case "Intensify Only":
          extractedSpectra = rejoinSpectra(intensifyOnly(spectrum), excludedPeaks);
          break;//just intensify ions related to targeted mass
        default:
          throw new Exception("Unrecognized return selection parameter"); //control
      }

      toSend(extractedSpectra); //transfer spectra to sending buffer
    }
  }

  public void OnResultsSent(IProcessingNode sender, MassSpectrumCollection result)
  {
    if (intensify == null) //initialize intensification method
    {
      if (this.doIntensification.Value) intensify = validIntensify;
      else intensify = voidIntensify;
    }

    if (exclusionList == null) //read exclusion list
    {
      try
      {
        exclusionList = loadExcludes();
      }
      catch (Exception ex)
      {
        SendAndLogErrorMessage(String.Format("Can't load exlusion list from '{0}': {1}", ignoreListPath.Value, ex.Message));
        throw;
      }
    }

    Log.Debug(String.Format("{0} spectra received", result.Count));
    spectra_in += result.Count;

    try
    {
      ProcessBatch(result);
    }
    catch (Exception ex)
    {
      SendAndLogErrorMessage("Error: " + ex.Message);
      throw;
    }
  }

  public override void OnParentNodeFinished(IProcessingNode sender, ResultsArguments eventArgs)
{
  flushSend();
  FireProcessingFinishedEvent(new ResultsArguments(MassSpecDataTypes.MSnSpectra));
}

  private void toSend (MassSpectrumCollection inspectra)
        */
/* Check if results can still fit in buffer and add them, otherwise send buffer and flush it *//*

  {
    ProcessingServices.SpectrumProcessingService.InitializeSpectra(this, inspectra);
    if (inspectra.Count + spectrumBuffer.Count > PackageSize)
    {
      SendResults(spectrumBuffer);

      Log.DebugFormat(String.Format("{0} spectra sent", spectrumBuffer.Count));
      SendAndLogTemporaryMessage(String.Format("{0} spectra received, {1} spectra sent", spectra_in, spectra_out), true);
      spectra_out += spectrumBuffer.Count;

      spectrumBuffer = new MassSpectrumCollection(PackageSize);
    }

    spectrumBuffer.AddRange(inspectra.AsEnumerable());
  }

  private void flushSend()
        */
/* Flush the all spectra in buffer *//*

  {
    SendResults(spectrumBuffer);
    Log.DebugFormat(String.Format("{0} spectra sent", spectrumBuffer.Count));
    spectra_out += spectrumBuffer.Count;
    SendAndLogMessage(String.Format("Processing finished: {0} spectra received, {1} spectra sent", spectra_in, spectra_out), true);
  }
}

*/
/* Represents an range consisting of several subranges
 * eg. (1.0, 2.1) and (5.02, 50.04) and ...
 * NOTE: borders aren't included
 *//*

class MultiRange extends Object
{
  private LinkedList<Double> points = new LinkedList<Double>();

  public double getLowest()
  {
    //lowest border of all ranges
    if (points != null && points.size() > 0) return points.get(0);
    else throw new RuntimeException("Empty MultiRange");
  }

  //highest border of all ranges
  public double getHighest()
  {
    if (points != null && points.size() > 0) return points.get(points.size() - 1);
    else throw new RuntimeException("Empty MultiRange");
  }

  public MultiRange() { super(); }

  //add a range to multirange
  public void add(double lower, double upper)
  {
    if (lower >= upper) return; //length of the range should be positive
    //position to insert values in points list
    int first = 0, last = 0;
    //find position of lower
    while (first < points.size() && lower > points.get(first)) first++;

    last = first;
    //find position of upper
    while (last < points.size() && upper > points.get(last)) last++;

    // complete pair to the left and to the right
    if (first % 2 == 0 && last % 2 == 0)
    {
      points.subList(first, last - first).clear();
      points.add(first, upper);
      points.add(first, lower);
    } else if (first % 2 == 0 && last % 2 == 1) //complete pair to the left, incomplete to the right
    {
      points.subList(first, last - first).clear();
      points.add(first, lower);
    } else if (first % 2 == 1 && last % 2 == 0) //incomplete pair to the left, complete to the right
    {
      points.subList(first, last - first).clear();
      points.add(first, upper);
    } else //incomplete pair to the left and to the right
    {
      points.subList(first, last - first).clear();
    }
  }

  //check if value is in multirnage
  public boolean isIn(double value)
  {
    //internal check
    if (points.size() % 2 == 1) throw new RuntimeException("Invalid range");

    for (int i = 0; i < points.size(); i += 2)
      if (value > points.get(i) && value < points.get(i + 1)) return true;

    return false;
  }

  public int lastIndexBelow(int start, double value)
  {
    return Tools.lastIndexBelow(points, start, value);
//    return points.FindIndex(start, x => x > value) - 1;
  }

  public int firstIndexAbove(int start, double value)
  {
    return Tools.firstIndexAbove(points, start, value);
//    return points.FindIndex(start, x => x > value);
  }

  public double get(int index) { return points.get(index); }

  public int size() { return points.size(); }
}
*/
/* Encloses list of MassCentroids, that are kept unique and each has a unique ID
 * Centroids withing certian precision are merged on addition and total intensity is stored
 *//*

class MassAggregator extends Object
{
  private MsnSpectrum massList = new MsnSpectrum();
  private double ppm;

  */
/* Initialize with certain ppm precision *//*

  public MassAggregator(double precision)
  {
    ppm = precision * 1e-6;
  }

  */
/* Check if mass can be merged with some of existing and if so the ID is returned;
   * otherwise mass is added to the list and new ID is generated *//*

  public int addPeak(Peak newMass)
  {
    for (int i = 0; i < massList.size(); i++) //check all masses, that are already in the list
    {
      if (Math.abs(massList.getMz(i)-newMass.getMz()) < massList.getMz(i)*ppm)
      {
        // remove the existing peak

        //merge peaks
//        massList.[i].Position = (massList[i].Position * massList[i].Intensity + newMass.Position * newMass.Intensity) /
//                               (massList[i].Intensity + newMass.Intensity);
//        massList[i].Intensity = massList[i].Intensity + newMass.Intensity;
//        massList[i].SignalToNoise = massList[i].SignalToNoise + newMass.SignalToNoise;
        return i;
      }
    }
    massList.Add(newMass);
    return massList.Count - 1;
  }

  public Peak get(int index) { return new Peak(massList.getMz(index), massList.getIntensity(index)); }
  public int size() { return massList.size(); }
  public PeakList getAll() { return massList;  }

*/
}
