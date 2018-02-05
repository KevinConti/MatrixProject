import org.junit.Test;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

public class MatrixTest {

    @Test
    public void givenAnAppropriateArrayAMatrixCanBeFormed(){
        double[][] table = new double[3][3];
        for (int row = 0; row < 3; row ++) {
            for (int col = 0; col < 3; col++) {
                table[row][col] = 1 + row * 3 + col;
            }
        }
        Matrix myMatrix = new Matrix(table);
        assertNotNull(myMatrix);
    }

    @Test(expected = ArrayIndexOutOfBoundsException.class)
    public void whenAddReceivesAWronglyFormattedMatrixAnErrorIsThrown(){
        double[][] tableOne = new double[][]{
                { 1, 2, 3 },
                { 4, 5, 6 }
        };
        double[][] tableTwo = new double[][]{
                { 3 },
                { 10 }
        };
        Matrix matrixOne = new Matrix(tableOne);
        Matrix matrixTwo = new Matrix(tableTwo);
        Matrix.add(matrixOne, matrixTwo);
    }

    @Test
    public void whenAddReceivesAppropriateMatricesReturnsTheCorrectResult(){
        double[][] tableOne = new double[][]{
                { 1, 2, 3 },
                { 4, 5, 6 }
        };
        double[][] tableTwo = new double[][]{
                { 3, 4, 5 },
                { 10, 12, 14 }
        };
        Matrix matrixOne = new Matrix(tableOne);
        Matrix matrixTwo = new Matrix(tableTwo);
        Matrix addedMatrix = Matrix.add(matrixOne, matrixTwo);

        double[][] correctResultTable = new double[][]{
                { 4, 6, 8 },
                { 14, 17, 20 }
        };
        Matrix correctResult = new Matrix(correctResultTable);
        assertArrayEquals(addedMatrix.getMatrix(), correctResult.getMatrix());
    }

    @Test
    public void whenParseVectorIsCalledOnAOneByTwoVectorTheAppropriateMatrixIsCreated(){
        Vector testVectorOne = new Vector(2.0, 6.0);
        Vector testVectorTwo = new Vector(3.0, 4.0);
        Vector testVectorThree = new Vector(3.0, 8.0);
        Vector testVectorFour = new Vector(4.0, 6.0);

        Matrix testMatrixOne = Matrix.parseVector(testVectorOne);
        Matrix testMatrixTwo = Matrix.parseVector(testVectorTwo);
        Matrix testMatrixThree = Matrix.parseVector(testVectorThree);
        Matrix testMatrixFour = Matrix.parseVector(testVectorFour);

        Matrix correctMatrixOne = new Matrix(new double[][]{
                {2.0, 6.0}
        });
        Matrix correctMatrixTwo = new Matrix(new double[][]{
                {3.0, 4.0}
        });
        Matrix correctMatrixThree = new Matrix(new double[][]{
                {3.0, 8.0}
        });
        Matrix correctMatrixFour = new Matrix(new double[][]{
                {4.0, 6.0}
        });

        assertArrayEquals(testMatrixOne.getMatrix(), correctMatrixOne.getMatrix());
        assertArrayEquals(testMatrixTwo.getMatrix(), correctMatrixTwo.getMatrix());
        assertArrayEquals(testMatrixThree.getMatrix(), correctMatrixThree.getMatrix());
        assertArrayEquals(testMatrixFour.getMatrix(), correctMatrixFour.getMatrix());
    }

    @Test
    public void whenParseVectorsIsCalledTheAppropriateMatricesareReturned(){
        Vector testVectorOne = new Vector(2.0, 6.0);
        Vector testVectorTwo = new Vector(3.0, 4.0);
        Vector testVectorThree = new Vector(3.0, 8.0);
        Vector testVectorFour = new Vector(4.0, 6.0);

        Vector[] vectors = new Vector[4];
        vectors[0] = testVectorOne;
        vectors[1] = testVectorTwo;
        vectors[2] = testVectorThree;
        vectors[3] = testVectorFour;

        Matrix correctMatrixOne = new Matrix(new double[][]{
                {2.0, 6.0}
        });
        Matrix correctMatrixTwo = new Matrix(new double[][]{
                {3.0, 4.0}
        });
        Matrix correctMatrixThree = new Matrix(new double[][]{
                {3.0, 8.0}
        });
        Matrix correctMatrixFour = new Matrix(new double[][]{
                {4.0, 6.0}
        });

        Matrix[] testMatrices = Matrix.parseVectors(vectors);

        assertArrayEquals(testMatrices[0].getMatrix(), correctMatrixOne.getMatrix());
        assertArrayEquals(testMatrices[1].getMatrix(), correctMatrixTwo.getMatrix());
        assertArrayEquals(testMatrices[2].getMatrix(), correctMatrixThree.getMatrix());
        assertArrayEquals(testMatrices[3].getMatrix(), correctMatrixFour.getMatrix());
    }

    @Test
    public void whenToCovarianceIsCalledThenASingleMatrixIsReturned(){
        Vector[] vectors = new Vector[2];
        vectors[0] = new Vector(2.0, 6.0);
        vectors[1] = new Vector(3.0, 4.0);

        Matrix covarianceMatrix = Matrix.toCovariance(vectors);
        assertNotNull(covarianceMatrix);
    }

    @Test
    public void whenInvertMatrixIsCalledTheReturnedMatrixIsInverted(){
        double[][] table = new double[][]{
                {1.0,2.0}
        };
        Matrix myMatrix = new Matrix(table);
        Matrix transposedMatrix = Matrix.transpose(myMatrix);
        assertEquals(transposedMatrix.getMatrix()[0][0], 1.0, 0);
        assertEquals(transposedMatrix.getMatrix()[1][0], 2.0, 0);
    }

    @Test
    public void whenMultiplicationIsCalledTheMatrixIsCorrectlyMultiplied(){
        Matrix[] matrices = initializeTestMatrices();
        Matrix toMultiply = new Matrix(new double[][]{{-1.0},{0.0}});
        Matrix multipliedMatrix;
        try {
            multipliedMatrix = Matrix.multiply(toMultiply, matrices[0]);
        } catch (Exception e){
            e.printStackTrace();
        }

    }

    @Test
    public void testDivideByScalarDestructive(){
        Matrix[] matrices = initializeTestMatrices();

        //Initialize correct tables for assertions. Assumes division by 2.0
        double[][] tableOne = new double[][]{
                {-0.5, 0.0}
        };
        double[][] tableTwo = new double[][]{
                {0.0, -1.0}
        };
        double[][] tableThree = new double[][]{
                {0.0, 1.0}
        };
        double[][] tableFour = new double[][]{
                {0.5, 0.0}
        };

        for(Matrix matrix: matrices){
            Matrix.divideByScalarDestructive(matrix, 2.0);
        }
        assertEquals(tableOne[0][0], matrices[0].getMatrix()[0][0],0);
        assertEquals(tableOne[0][1], matrices[0].getMatrix()[0][1],0);
        assertEquals(tableTwo[0][0], matrices[1].getMatrix()[0][0],0);
        assertEquals(tableTwo[0][1], matrices[1].getMatrix()[0][1],0);
        assertEquals(tableThree[0][1], matrices[2].getMatrix()[0][1],0);
        assertEquals(tableThree[0][1], matrices[2].getMatrix()[0][1],0);
        assertEquals(tableFour[0][0], matrices[3].getMatrix()[0][0],0);
        assertEquals(tableFour[0][1], matrices[3].getMatrix()[0][1],0);
    }

    @Test
    public void testDivideByScalarDestructiveRows(){
        double[][] table = new double[][]{
                {2.0, 4.0, 6.0},
                {3.0, 5.0, 7.0}
        };
        Matrix matrix = new Matrix(table);
        Matrix.divideByScalarDestructive(matrix, 2.0, 0);

        assertEquals(1.0, matrix.getValueAt(0,0), 0);
        assertEquals(2.0, matrix.getValueAt(0,1), 0);
        assertEquals(3.0, matrix.getValueAt(0,2), 0);
        assertEquals(3.0, matrix.getValueAt(1,0), 0);
        assertEquals(5.0, matrix.getValueAt(1,1), 0);
        assertEquals(7.0, matrix.getValueAt(1,2), 0);

        Matrix.divideByScalarDestructive(matrix, -2.0, 1);

        assertEquals(-1.5, matrix.getValueAt(1,0), 0);
        assertEquals(-2.5, matrix.getValueAt(1,1), 0);
        assertEquals(-3.5, matrix.getValueAt(1,2), 0);
    }

    @Test
    public void testInverse(){
        //Initialize test matrices
        double[][] squareTable = new double[][]{
                {1, -1, 0},
                {-2, 2, -1},
                {0, 1, -2}
        };
        double[][] coefficientTable = new double[][]{
                {2},
                {-1},
                {6}
        };
        Matrix squareMatrix = new Matrix(squareTable);
        Matrix coefficientMatrix = new Matrix(coefficientTable);
        Matrix augmentedMatrix = null;
        try {
            augmentedMatrix = Matrix.createAugmentedMatrix(squareMatrix, coefficientMatrix);
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        //test matrix params
        assertEquals(4, augmentedMatrix.numColumns(), 0);
        assertEquals(3, augmentedMatrix.numRows(), 0);

        //Create answerMatrix
        double[][]answerTable = new double[][]{
                {1, 0, 0, 2},
                {0, 1, 0, 0},
                {0, 0, 1, -3}
        };
        Matrix answerMatrix = new Matrix(answerTable);

        //Create inverseMatrices with both versions of the method
        try {
            augmentedMatrix.inverse();
        } catch (InversionException e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
        Matrix secondInverseMatrix = null;
        try {
            secondInverseMatrix = squareMatrix.inverse(coefficientMatrix);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
        }
        //test matrix params
        assertEquals(4, augmentedMatrix.numColumns(), 0);
        assertEquals(3, augmentedMatrix.numRows(), 0);
        //test equalities
        for (int i = 0; i < answerMatrix.numRows(); i++){
            for(int j = 0; j < answerMatrix.numColumns(); j++){
                assertEquals(answerMatrix.getValueAt(i, j), augmentedMatrix.getValueAt(i, j), 0);
                assertEquals(answerMatrix.getValueAt(i, j), secondInverseMatrix.getValueAt(i, j), 0);
            }
        }

    }

    @Test
    public void testPivot(){
        Matrix augmentedMatrix = createTestAugmentedMatrix();
        augmentedMatrix = augmentedMatrix.pivot(1,0);

        assertEquals(-2.0, augmentedMatrix.getValueAt(0,0), 0);
        assertEquals(2.0, augmentedMatrix.getValueAt(0,1), 0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(0,2), 0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(0,3), 0);

        assertEquals(1.0, augmentedMatrix.getValueAt(1,0), 0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(1,1), 0);
        assertEquals(0.0, augmentedMatrix.getValueAt(1,2), 0);
        assertEquals(2.0, augmentedMatrix.getValueAt(1,3), 0);

        assertEquals(0.0, augmentedMatrix.getValueAt(2,0), 0);
        assertEquals(1.0, augmentedMatrix.getValueAt(2,1), 0);
        assertEquals(-2.0, augmentedMatrix.getValueAt(2,2), 0);
        assertEquals(6.0, augmentedMatrix.getValueAt(2,3), 0);

        augmentedMatrix = augmentedMatrix.pivot(2,0);

        assertEquals(0.0, augmentedMatrix.getValueAt(0,0), 0);
        assertEquals(1.0, augmentedMatrix.getValueAt(0,1), 0);
        assertEquals(-2.0, augmentedMatrix.getValueAt(0,2), 0);
        assertEquals(6.0, augmentedMatrix.getValueAt(0,3), 0);

        assertEquals(1.0, augmentedMatrix.getValueAt(1,0), 0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(1,1), 0);
        assertEquals(0.0, augmentedMatrix.getValueAt(1,2), 0);
        assertEquals(2.0, augmentedMatrix.getValueAt(1,3), 0);

        assertEquals(-2.0, augmentedMatrix.getValueAt(2,0), 0);
        assertEquals(2.0, augmentedMatrix.getValueAt(2,1), 0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(2,2), 0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(2,3), 0);
    }

    @Test
    public void testMatrixMean(){
        Matrix[] matrices = initializeTestMatrices();
        Matrix matrixMean = Matrix.matrixMean(matrices);
        double[] arrayValue = {0.0,0.0};
        assertEquals(arrayValue[0], matrixMean.getMatrix()[0][0], 0);
        assertEquals(arrayValue[1], matrixMean.getMatrix()[0][1], 0);
    }

    @Test
    public void testCreateAugmentedMatrix(){

        double[][] coefficientMatrixTable = new double[][]{
                {2},
                {-1},
                {6}
        };

        Matrix augmentedMatrix = createTestAugmentedMatrix();

        assertEquals(coefficientMatrixTable[0][0], augmentedMatrix.getMatrix()[0][3], 0);
        assertEquals(coefficientMatrixTable[1][0], augmentedMatrix.getMatrix()[1][3], 0);
        assertEquals(coefficientMatrixTable[2][0], augmentedMatrix.getMatrix()[2][3], 0);
    }

    @Test
    public void testLargestAbsoluteValue(){
        Matrix[] matrices = initializeTestMatrices();
        double[][] tableOne = new double[][]{
                {-1.0, 0.0},
                {5.0, -2.0},
                {1.0, -7.0}
        };
        Matrix testMatrix = new Matrix(tableOne);
        assertEquals(1, testMatrix.largestAbsoluteValueIndex(0), 0);
        assertEquals(2, testMatrix.largestAbsoluteValueIndex(1), 0);

    }

    private Matrix[] initializeTestMatrices(){
        Matrix[] matrices = new Matrix[4];
        double[][] tableOne = new double[][]{
                {-1.0, 0.0}
        };
        double[][] tableTwo = new double[][]{
                {0.0, -2.0}
        };
        double[][] tableThree = new double[][]{
                {0.0, 2.0}
        };
        double[][] tableFour = new double[][]{
                {1.0, 0.0}
        };
        matrices[0] = new Matrix(tableOne);
        matrices[1] = new Matrix(tableTwo);
        matrices[2] = new Matrix(tableThree);
        matrices[3] = new Matrix(tableFour);

        return matrices;
    }

    private Matrix createTestAugmentedMatrix(){
        //Initialize test matrices
        double[][] squareTable = new double[][]{
                {1, -1, 0},
                {-2, 2, -1},
                {0, 1, -2}
        };
        double[][] coefficientTable = new double[][]{
                {2},
                {-1},
                {6}
        };
        Matrix squareMatrix = new Matrix(squareTable);
        Matrix coefficientMatrix = new Matrix(coefficientTable);
        Matrix augmentedMatrix = null;
        try {
            augmentedMatrix = Matrix.createAugmentedMatrix(squareMatrix, coefficientMatrix);
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        return  augmentedMatrix;
    }
}
