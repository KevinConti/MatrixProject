import org.junit.Test;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertEquals;

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
}
