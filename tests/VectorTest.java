import org.junit.Test;

import java.util.ArrayList;

import static org.junit.Assert.*;

public class VectorTest {
    @Test
    public void whenInitializedProperlyAVectorIsReturned(){
        double x = 3.25;
        double y = -12.43;

        Vector testVector = new Vector(x,y);
        System.out.println(testVector);
        assertEquals(testVector.getX(), 3.25, 0);
        assertEquals(testVector.getY(), -12.43, 0);
    }

    @Test
    public void whenMeanMethodIsCalledTheMeanIsReturned() {
        double[] xValues = {1, 2.5, 3, 4, 5.0, 6};
        double[] yValues = {2, 4, 6, 8.8, 10.5, 12};

        int numOfTestValues = 6;

        double actualXMean;
        double actualYMean;

        //Temp variables for calculating actualX and actualY
        double sumX = 0;
        double sumY = 0;

        for(int i = 0; i < numOfTestValues ; i++){
            sumX = sumX + xValues[i];
            sumY = sumY + yValues[i];
        }
        actualXMean = sumX/numOfTestValues;
        actualYMean = sumY/numOfTestValues;

        System.out.println("X Mean: "+ actualXMean);
        System.out.println("Y Mean: "+ actualYMean);

        Vector[] vectors = new Vector[numOfTestValues];

        for(int i = 0; i < numOfTestValues; i++){
            Vector myVector = new Vector(xValues[i], yValues[i]);
            vectors[i] = myVector;
        }

        Vector meanVector = Vector.mean(vectors);
        assertEquals(meanVector.getX(), actualXMean, 0);
        assertEquals(meanVector.getY(), actualYMean, 0);
    }
}