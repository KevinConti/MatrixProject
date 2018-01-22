import java.util.Arrays;

public class Matrix {

    double[][] matrix;

    public Matrix(double[][] matrix) {
        this.matrix = matrix;
    }

    public double[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(double[][] matrix) {
        this.matrix = matrix;
    }

    //Class Methods
    public static Matrix add(Matrix matrixOne, Matrix matrixTwo){
        //Check if matrices are appropriately sized (both must be MxN, returns a MxN matrix
        //Add the two together, store the result in the a new matrix that is returned
        double[][] addedMatrixArray = new double[matrixOne.getMatrix()[0].length][matrixOne.getMatrix()[1].length];
        for(int row = 0; row < addedMatrixArray[0].length; row++){
            for(int column = 0; column < addedMatrixArray[1].length; column++){
                addedMatrixArray[row][column] = matrixOne.getMatrix()[row][column] + matrixTwo.getMatrix()[row][column];
            }
        }
        return new Matrix(addedMatrixArray);
    }

    @Override
    public String toString() {
        StringBuilder fullArray = new StringBuilder();
        for(int row = 0; row < matrix.length; row++){
            fullArray.append(Arrays.toString(matrix[row]));
        }
        return fullArray.toString();
    }
}
