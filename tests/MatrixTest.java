import org.junit.Test;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.*;

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
    public void testSubtract(){
        Matrix one = new Matrix(new double[][]{
                {1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}
        });

        Matrix two = new Matrix(new double[][]{
                {1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}
        });

        Matrix subtracted = Matrix.subtract(one, two);
        for(int i = 0; i < subtracted.numRows(); i++){
            for(int j = 0; j < subtracted.numColumns(); j++){
                assertEquals(0, subtracted.getValueAt(i, j), 0);
            }
        }
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
    public void testTranspose(){
        double[][] table = new double[][]{
                {1.0,2.0}
        };
        double[][] tableTwo = new double[][]{
                {-2, 5, -1},
                {5, -7, 8},
                {6, -4, 3}
        };
        Matrix myMatrix = new Matrix(table);
        Matrix myMatrixTwo = new Matrix(tableTwo);
        Matrix transposedMatrixOne = Matrix.transpose(myMatrix);
        assertEquals(transposedMatrixOne.getMatrix()[0][0], 1.0, 0);
        assertEquals(transposedMatrixOne.getMatrix()[1][0], 2.0, 0);

        myMatrixTwo = Matrix.transpose(myMatrixTwo);
        assertEquals(-2.0, myMatrixTwo.getValueAt(0,0), 0);
        assertEquals(5.0, myMatrixTwo.getValueAt(0,1), 0);
        assertEquals(6.0, myMatrixTwo.getValueAt(0,2), 0);
        assertEquals(5.0, myMatrixTwo.getValueAt(1,0), 0);
        assertEquals(-7.0, myMatrixTwo.getValueAt(1,1), 0);
        assertEquals(-4.0, myMatrixTwo.getValueAt(1,2), 0);
        assertEquals(-1.0, myMatrixTwo.getValueAt(2,0), 0);
        assertEquals(8.0, myMatrixTwo.getValueAt(2,1), 0);
        assertEquals(3.0, myMatrixTwo.getValueAt(2,2), 0);
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
            fail();
        }

        //Multiply a 1x2 * 2x1 matrix and return a 2x2 matrix
        Matrix oneByTwo = new Matrix(1,2);
        oneByTwo.setMatrix(0,0, 5);
        oneByTwo.setMatrix(0,1,10);
        Matrix twoByTwo = new Matrix(2,2);
        for(int i = 0; i < twoByTwo.numRows(); i++){
            twoByTwo.setMatrix(i, 0, i);
            twoByTwo.setMatrix(i, 1, i+1);
        }
        try {
            Matrix.multiply(oneByTwo, twoByTwo);
        } catch (Exception e) {
            e.printStackTrace();
            fail();
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
    public void testMultiplyByScalar(){
        //Create test matrix
        Matrix matrix = new Matrix(new double[][]{
                {1,2,3},
                {4,5,6},
                {7,8,9}
        });
        Matrix copyMatrix = matrix.copy();
        //Create test scalar
        double scalar = 2.5;
        //Do multiplication
        Matrix resultMatrix = matrix.multiplyByScalar(scalar);
        //Verify multiplication was successful
        for(int i = 0; i < matrix.numColumns(); i++){
            for(int j = 0; j < matrix.numRows(); j++){
                assertEquals(matrix.getValueAt(i, j) * 2.5, resultMatrix.getValueAt(i, j), 0);
            }
        }
        //Verify original matrix was not modified
        for(int i = 0; i < matrix.numColumns(); i++){
            for(int j = 0; j < matrix.numRows(); j++){
                assertEquals(matrix.getValueAt(i, j), copyMatrix.getValueAt(i, j), 0);
            }
        }
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
        Matrix secondInverseMatrix = new Matrix(squareTable);
        try {
            secondInverseMatrix.inverse(coefficientMatrix);
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
    public void testRemoveFirstColumn(){
        Matrix[] matrices = initializeTestMatrices();
        double[][] answerTable = new double[][]{
                {0.0},
                {-2.0},
                {2.0},
                {0.0}
        };
        for (int i = 0; i < answerTable.length; i++){
            matrices[i].removeFirstColumn();

            assertEquals(answerTable[i][0], matrices[i].getValueAt(0,0), 0);
        }

        //Test a larger matrix
        Matrix augmented = createTestAugmentedMatrix();
        //Remove three columns so just the coefficients remain
        augmented.removeFirstColumn();
        augmented.removeFirstColumn();
        augmented.removeFirstColumn();
        //Answers should be 2, -1, 6
        answerTable = new double[][]{
                {2.0},
                {-1.0},
                {6.0}
        };
        for(int i = 0; i < answerTable.length; i++){
            assertEquals(answerTable[i][0], augmented.getValueAt(i,0), 0);
        }
    }

    @Test
    public void testInverseWithIdentity(){
        double[][] squareTable = new double[][]{
                {1, -1, 0},
                {-2, 2, -1},
                {0, 1, -2}
        };
        double[][] inverseTable = new double[][]{
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}
        };

        Matrix squareMatrix = new Matrix(squareTable);
        Matrix inverseMatrix = new Matrix(inverseTable);
        Matrix augmentedMatrix = null;
        try {
            augmentedMatrix = Matrix.createAugmentedMatrix(squareMatrix, inverseMatrix);
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        //test matrix params
        assertEquals(6, augmentedMatrix.numColumns(), 0);
        assertEquals(3, augmentedMatrix.numRows(), 0);

        try {
            augmentedMatrix.inverse();
        } catch (InversionException e) {
            e.printStackTrace();
        }

        assertEquals(-3.0, augmentedMatrix.getValueAt(0,3),0);
        assertEquals(-2.0, augmentedMatrix.getValueAt(0,4),0);
        assertEquals(1.0, augmentedMatrix.getValueAt(0,5),0);

        assertEquals(-4.0, augmentedMatrix.getValueAt(1,3),0);
        assertEquals(-2.0, augmentedMatrix.getValueAt(1,4),0);
        assertEquals(1.0, augmentedMatrix.getValueAt(1,5),0);

        assertEquals(-2.0, augmentedMatrix.getValueAt(2,3),0);
        assertEquals(-1.0, augmentedMatrix.getValueAt(2,4),0);
        assertEquals(0.0, augmentedMatrix.getValueAt(2,5),0);
    }

    @Test
    public void testDeterminate(){
        double[][] table = new double[][]{
                {5, -2, 3},
                {0, 7, -8},
                {4, 3, 1}
        };
        Matrix testMatrix = new Matrix(table);
        double determinate = 0;
        try {
            determinate = testMatrix.determinate();
        } catch (InversionException e) {
            e.printStackTrace();
        }
        assertEquals(-135.0, determinate, 0);
    }

    @Test
    public void testCopy(){
        Matrix original = initializeTestMatrices()[0];
        Matrix copy = original.copy();

        //Test copy is correct values
        assertEquals(original.getValueAt(0,0), copy.getValueAt(0,0),0);
        assertEquals(original.getValueAt(0,1), copy.getValueAt(0,1),0);
        //Test copy isn't a reference to original
        copy.setMatrix(0,0, 10.5);
        assertEquals(-1.0, original.getValueAt(0,0),0);
        assertEquals(10.5, copy.getValueAt(0,0), 0);
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
    public void testDivide(){
        Matrix testMatrix = new Matrix(new double[][]{
                {2, 3, 4},
                {-2, 3, 4},
                {5, 6, 7}
        });

        Matrix multMatrix = new Matrix(new double[][]{
                {1.0/2, 1.0/3, 1.0/4},
                {1.0/5, 1.0/6, 1.0/7},
                {1.0/8, 1.0/9, 1.0/10}
        });

        Matrix divMatrix = new Matrix(new double[][]{
                {2, 3, 4},
                {5, 6, 7},
                {8, 9, 10}
        });
        try {
            multMatrix = Matrix.multiply(testMatrix, multMatrix);
            divMatrix = Matrix.divide(testMatrix, divMatrix);


            for(int i = 0; i < testMatrix.numRows(); i++){
                for(int j = 0; j < testMatrix.numColumns(); j++){
                    try {
                        assertEquals(multMatrix.getValueAt(i,j),divMatrix.getValueAt(i,j), 0);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }
        }
        catch (Exception e){
            e.printStackTrace();
        }
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

    @Test
    public void testMaximum(){
        double[] values = new double[]{-3, 4, 10, -999, 400, 32, 44};
        double largestValue = Matrix.maximum(values);
        assertEquals(400.0, largestValue, 0);
    }

    @Test
    public void testRowSums(){
        Matrix coefficientMatrix = createTestAugmentedMatrix();
        coefficientMatrix.removeLastColumn();
        double[] rowSums = coefficientMatrix.rowSums();
        assertEquals(2.0, rowSums[0], 0);
        assertEquals(5.0, rowSums[1], 0);
        assertEquals(3.0, rowSums[2], 0);
    }

    @Test
    public void testConditionNumber(){
        Matrix coefficientMatrix = new Matrix(new double[][]{
                {100, -200},
                {-200, 401}
        });
        double conditionNumber = coefficientMatrix.conditionNumber();
        assertEquals(3612.0, conditionNumber, .01);

    }

    @Test
    public void testLeverriersMethod(){
        //Test 3x3 Matrix
        Matrix testMatrix = new Matrix(new double[][]{
                {1, -1, 0},
                {0, 2, -1},
                {-1, 0, 1}
        });

        Matrix correctAnswer = new Matrix(new double[][]{
                {-4, 0, 0},
                {0, 5, 0},
                {0, 0, -1}
        });

        Matrix resultMatrix = new Matrix(3,3);
        try {
            resultMatrix = testMatrix.leverriersMethod();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        for(int i = 0; i < correctAnswer.numRows(); i++){
            for(int j = 0; j < correctAnswer.numColumns(); j++) {
                assertEquals(correctAnswer.getValueAt(i, j), resultMatrix.getValueAt(i, j), 0);
            }
        }
        //test 1x1 matrix (assert error is thrown)
        testMatrix = new Matrix(1,1);
        boolean flag = true;
        try {
            resultMatrix = testMatrix.leverriersMethod();
            flag = false; //Throws if no exception is thrown
        } catch (Exception e) {
            assertTrue(flag);
        }
        //test non-square matrix (assertFail)
        testMatrix = new Matrix(1,3);
        flag = true;
        try {
            resultMatrix = testMatrix.leverriersMethod();
            flag = false; //Throws if no exception is thrown
        } catch (Exception e) {
            assertTrue(flag);
        }

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
