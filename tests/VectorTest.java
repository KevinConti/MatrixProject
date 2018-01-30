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

    @Test
    public void whenSubtractIsCalledWithTwoVectorsAsParamsThenTheTwoVectorsAreSubtractedFromEachOther(){
        Vector testOne = new Vector(2.0, 6.0);
        Vector testTwo = new Vector(3.0, 6.0);

        Vector result = Vector.subtract(testOne, testTwo);
        assertEquals(result.getX(), -1.0, 0);
        assertEquals(result.getY(), 0, 0);
    }

    @Test
    public void whenSubtractIsCalledWithAnArrayOfVectorsAndASingleVectorAsParamsThenEachVectorInArrayIsSubtractedByTheSingleVector(){
        Vector[] vectors = initializeTestVectors();
        Vector subtrahend = new Vector(3.0, 6.0);
        Vector[] subtractedVectors = Vector.subtract(vectors, subtrahend);

        Vector[] correctVectors = new Vector[4];
        correctVectors[0] = new Vector(-1.0, 0.0);
        correctVectors[1] = new Vector(0.0, -2.0);
        correctVectors[2] = new Vector(0.0, 2.0);
        correctVectors[3] = new Vector(1.0, 0.0);

        for (int i = 0; i < subtractedVectors.length; i++){
            assertEquals(subtractedVectors[i].getX(), correctVectors[i].getX(), 0);
            assertEquals(subtractedVectors[i].getY(), correctVectors[i].getY(), 0);
        }
    }

    private Vector[] initializeTestVectors(){
        Vector[] vectors = new Vector[4];
        vectors[0] = new Vector(2.0, 6.0);
        vectors[1] = new Vector(3.0, 4.0);
        vectors[2] = new Vector(3.0, 8.0);
        vectors[3] = new Vector(4.0, 6.0);

        return vectors;

    }
}