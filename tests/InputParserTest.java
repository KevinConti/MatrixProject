import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class InputParserTest {

    @Test
    public void whenInputParserIsConstructedArrayShouldHaveProperSize(){
        InputParser ip = new InputParser();
        Vector[][] parserResult = ip.convertFile("data/data.txt");

        final int LENGTH_OF_CLASS = 110;
        final int NUM_CLASSES = 2;

        assertEquals(parserResult.length, NUM_CLASSES, 0);
        assertEquals(parserResult[0].length, LENGTH_OF_CLASS, 0);
        assertEquals(parserResult[1].length, LENGTH_OF_CLASS, 0);
    }
}