package org.ms2ms.test;

import org.expasy.mzjava.core.ms.peaklist.*;
import org.expasy.mzjava.core.ms.spectrum.IonType;
import org.expasy.mzjava.proteomics.mol.Peptide;
import org.expasy.mzjava.proteomics.ms.spectrum.PepFragAnnotation;
import org.junit.Test;

/**
 * Created by wyu on 4/13/14.
 */
public class SpectrumTest extends TestAbstract
{
    @Test
    public void recipe1_1()
    {
    /*  A PeakList is an interface that handles a list of peaks defined by mz and intensity values with optional
        PeakAnnotations.

        Spectrum is the base class that is used to represent the peak list and meta data associated with an
        msn spectrum. The peak list data is held by the classes that implement PeakList. For each peak an mz and
        intensity value as well as a list of peak annotations is stored. MzJava provides 5 flavors of PeakList
        which differ in the precision that is used to store the m/z and intensity values. M/z values are either
        stored in a float or double array and intensity values can be constant, an array of floats or an array
        of doubles. All implementations of PeakList guarantee that the peaks are sorted. Peak annotations are
        represented by classes that implement PeakAnnotation. The type of PeakAnnotation that a PeakList holds is
        specified using a generic. Each PeakList contains a randomly generated UUID to simplify identification. We
        provide a bunch of PeakList implementations given the precision of m/zs and intensities.
     */

        PeakList<PeakAnnotation> peaklist;

        // In the following, peaks are stored in a single-precision floating-point format:
        peaklist = new FloatPeakList<PeakAnnotation>(); // single-precision floating-point m/zs and intensities
        peaklist = new DoublePeakList<PeakAnnotation>(); // double-precision floating-point m/zs and intensities

        // Next are some implementations where peaks are stored in a mixed mz/intensity precision:
        peaklist = new DoubleConstantPeakList<PeakAnnotation>(1); // double-precision floating-point m/zs and constant intensities of 1

        // We also provide factory to ease PeakList instanciation:
        peaklist = PeakListFactory.newPeakList(PeakList.Precision.DOUBLE_FLOAT);

        // Adding peaks:
        peaklist.add(1, 2);
        peaklist.add(2, 4);
        peaklist.add(3, 6);

        // Adding peak annotations:
        peaklist.add(1, 2);
        peaklist.addAnnotation(0, new PepFragAnnotation(IonType.b, 2, Peptide.parse("CE")));

        // Traversing a Peaklist:
        for (int i=0 ; i<peaklist.size() ; i++)
        {
            peaklist.getMz(i);
            peaklist.getIntensity(i);
            peaklist.getAnnotations(i);
        }

        // Manipulating a Peaklist using a cursor:
        peaklist = new DoublePeakList<PeakAnnotation>();
        peaklist.addSorted(new double[] {101.0909, 202.1079, 203.1045, 244.128, 254.1056, 255.191, 270.2388, 272.2092},
                new double[]{15.4762, 21.4762, 5.3333, 14.381, 3.4286, 7.2381, 2.1905, 27.8095});

        PeakCursor<PeakAnnotation> cursor = peaklist.cursor();

        // move to next peak (peak.mz=101.0909)
        cursor.next();

        // move to next peak with constraint. Here the cursor moves as next peak.m/z <= 250 (current.mz=202.1079).
        cursor.next(250);

        // accessors
        cursor.currMz();
        cursor.currIntensity();
        cursor.currAnnotations();

        // move to next peak with constraint. Here the cursor does not move as next peak.m/z > 200 (current.mz=202.1079).
        cursor.next(200);
        // move to peak where m/z < 240 (current.mz=203.1045)
        cursor.moveBefore(240);
        // move to peak where m/z > 240 (current.mz=244.128)
        cursor.movePast(240);
        // move to the closest peak (current.mz=202.1079)
        cursor.moveToClosest(200);

        // methods to retrieve a peak that is at n peaks distant from the current position without moving the cursor.
        if (cursor.canPeek(4))
        {
            // here peeked peak contains mz 255.191 and intensity 7.2381
            double mz = cursor.peekMz(4);
            double intensity = cursor.peekIntensity(4);
        }
    }
    @Test
    public void recipe1_2()
    {

    }
}
