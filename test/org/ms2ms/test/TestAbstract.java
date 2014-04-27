package org.ms2ms.test;

import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** A skeleton class for JUnit4 test
 *
 * Created by wyu on 4/13/14.
 */
@RunWith(JUnit4.class)
abstract public class TestAbstract
{
    @Ignore("Test is ignored as a demonstration")
    @Test
    public void testSane()
    {
        assert(1==1);
    }
}
