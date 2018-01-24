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

    @Test
    public void whenConvertFileIsCalledTwoCorrectClassesAreCreated(){
        InputParser ip = new InputParser();
        Vector[][] parserResult = ip.convertFile("data/data.txt");

        Vector classOneVector = parserResult[0][0];
        assertEquals(classOneVector.getX(), 2.566613444, 0);
        assertEquals(classOneVector.getY(), 3.244925222, 0);

        Vector classTwoVector = parserResult[1][0];
        assertEquals(classTwoVector.getX(), 0.587942192, 0);
        assertEquals(classTwoVector.getY(), -1.256260722, 0);

        classOneVector = parserResult[0][49];
        classTwoVector = parserResult[1][49];

        assertEquals(classOneVector.getX(), 2.365301895, 0);
        assertEquals(classOneVector.getY(), 3.16735623, 0);

        assertEquals(classTwoVector.getX(), 0.579117114, 0);
        assertEquals(classTwoVector.getY(), 1.508560248, 0);
    }
}