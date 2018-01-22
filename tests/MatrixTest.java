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
        System.out.println(myMatrix);
        assertNotNull(myMatrix);
    }
}
